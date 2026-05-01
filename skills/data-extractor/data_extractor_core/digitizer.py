"""Digitize plot data from figure images using Claude vision."""

from __future__ import annotations

import asyncio
import base64
import io
import json
import logging
import re

from PIL import Image

from .models import Confidence, DataSeries, ExtractedData, Figure, PlotType

logger = logging.getLogger(__name__)

# ── Plot-type-specific prompts ──

BASE_PROMPT = """\
You are a scientific data extraction specialist performing digitization of \
figures from published papers for meta-analysis. Your goal is to extract \
ALL numerical data as accurately as possible.

{legend_context}

{panel_focus}

STEP 1 — IDENTIFY AXES AND SCALE:
- Read ALL tick mark values on both axes. Write them down mentally.
- Determine the scale: what is the range? What does each gridline represent?
- CRITICAL — detect LOG SCALE: If tick marks are spaced like 0.001, 0.01, 0.1, 1, 10 \
  (powers of 10) or 1, 2, 5, 10, 20, 50, 100 (multiplicative), the axis is LOG-SCALED. \
  Set x_scale or y_scale to "log" accordingly. Common in dose-response, survival, \
  forest plots, and any plot with orders-of-magnitude range.
- Read axis labels and units.
- Identify all series/groups using the color legend if visible in the image.

STEP 2 — NAMING: Use the REAL names for each series:
- If a color legend is visible in the image, match each bar/line/point color \
to its label from the legend. NEVER use generic names like "Bar 1", \
"red bar", or "Group 1" if a legend is present.
- If the figure legend text names the series, use those names.
- If no legend is visible and no names are given, use descriptive colors as a last resort.

STEP 3 — READ VALUES PRECISELY:
- For EACH data point, anchor your reading to the nearest axis tick marks.
- Example: if ticks are at 0, 100, 200, 300 and a bar top is halfway between \
200 and 300, the value is 250 (not "approximately 200" or "around 300").
- For bars: find the top of each bar, project horizontally to the y-axis, \
and read the value relative to the closest tick marks above and below.
- If COMPUTER VISION CALIBRATION DATA is provided below, USE IT:
  The calibration gives pixel coordinates of detected markers/bars.
  Map each pixel position to an axis value by interpolating between ticks.
  This is MORE ACCURATE than visual estimation — trust the pixel coordinates.
- For line plots: read the y-value at each x-axis tick mark by projecting \
vertically up to the curve, then horizontally to the y-axis.
- For scatter points: read both x and y by projecting to the nearest ticks.

ERROR BARS — THIS IS CRITICAL, READ CAREFULLY:
- error_bars_lower and error_bars_upper are the ± EXTENT (delta), NOT absolute positions.
- First, read the ABSOLUTE y-position of the top and bottom of the error bar from the axis.
- Then SUBTRACT the mean to get the extent:
  error_bars_upper = (top of error bar) − (mean value)
  error_bars_lower = (mean value) − (bottom of error bar)
- Both values should be POSITIVE numbers (they are distances from the mean).
- Example: if a bar height (mean) is 250 and error bar goes from 220 to 280:
  error_bars_lower = 250 − 220 = 30
  error_bars_upper = 280 − 250 = 30
- Example: if mean = 5.2, error bar bottom = 4.1, error bar top = 6.8:
  error_bars_lower = 5.2 − 4.1 = 1.1
  error_bars_upper = 6.8 − 5.2 = 1.6
- For SYMMETRIC error bars (same distance above and below), both values are equal.
- For ASYMMETRIC error bars (e.g., confidence intervals), read EACH end
  INDEPENDENTLY from the axis — do NOT assume they are equal.
- NEVER report the absolute axis position of the error bar tip as the error value.
- If NO error bars are visible ANYWHERE in the figure, OMIT the \
"error_bars_lower" and "error_bars_upper" keys entirely from the series. \
Do NOT include them as empty arrays or arrays of nulls.
- Only include error bar arrays when error bars are actually visible.
- If some points have error bars and others don't, use null for the missing ones.

{plot_specific_guidance}

STEP 4 — OUTPUT: Return a JSON object with this exact structure:
{{
  "plot_type": "<type>",
  "title": "<descriptive title from the legend/caption — NEVER just 'Table 1' or 'Figure 3', always describe what the data shows>",
  "x_label": "<x-axis label>",
  "y_label": "<y-axis label>",
  "x_unit": "<x-axis unit if shown, or null>",
  "y_unit": "<y-axis unit if shown, or null>",
  "x_min": <first (lowest) tick value on x-axis, or null if categorical>,
  "x_max": <last (highest) tick value on x-axis, or null if categorical>,
  "y_min": <first (lowest) tick value on y-axis>,
  "y_max": <last (highest) tick value on y-axis>,
  "x_scale": "<MUST be exactly 'linear' or 'log'>",
  "y_scale": "<MUST be exactly 'linear' or 'log'>",
  "series": [
    {{
      "name": "<series/group name from legend — use real names, not colors>",
      "x_values": [...],
      "y_values": [...],
      "error_bars_lower": [<positive ± extent BELOW mean>],
      "error_bars_upper": [<positive ± extent ABOVE mean>]
    }}
  ],
  "confidence": "<high|medium|low>",
  "notes": "<error bar type (SD/SEM/CI), sample sizes if shown, any ambiguity>"
}}

IMPORTANT — AXIS RANGES: x_min/x_max/y_min/y_max must match the ORIGINAL figure's \
axis range (the first and last tick marks), NOT the range of your extracted data. \
This ensures the reproduced plot uses the same scale as the original.

Confidence guide:
- "high": Clear axes, well-separated points, no ambiguity
- "medium": Some interpolation needed, minor overlaps
- "low": Dense overlapping points, unclear axes, or very small figure

Return ONLY the JSON object, no other text."""

PLOT_GUIDANCE = {
    PlotType.SCATTER: """Scatter plot guidance:
- Extract EVERY visible data point, even if there are many (hundreds is fine).
- For each point, report precise (x, y) coordinates.
- If there are distinct groups/colors, create separate series for each.
- If a regression/trend line is shown, note its equation and R² in notes but
  extract the raw data points — do NOT extract points along the trend line.
- LOG SCALE: If either axis uses log scale (ticks at 1, 10, 100, 1000 or
  0.01, 0.1, 1, 10), report the ACTUAL values (not log-transformed).
  E.g., a point halfway between 10 and 100 on a log axis ≈ 30, not 55.
- For dense clusters where individual points overlap, estimate the count of
  overlapping points and note it, but still extract representative coordinates.
- If point sizes vary (bubble-like), note this and see bubble guidance.
- Significance annotations (* p<0.05 etc.) between groups → note in notes.""",

    PlotType.BAR: """Bar chart guidance:
- x_values should be category labels (strings), not numbers.
- HORIZONTAL BARS: If bars are horizontal, swap your reading — the bar length
  is the value (read from x-axis), categories are on the y-axis. Still report
  categories in x_values and bar lengths in y_values.

INDIVIDUAL DATA POINTS ON BARS — THIS IS THE MOST IMPORTANT RULE:
If dots/circles/jitter points are overlaid on the bars, the individual data
points are the REAL DATA — bars are just summaries (means). Extract the points.

Step 1: IDENTIFY THE GROUPING STRUCTURE.
  Look at the legend carefully. There are typically TWO layers of grouping:
  (a) X-axis CATEGORIES: labels under bars (e.g., "TH+M-", "TH+M+", "TH-M+")
  (b) LEGEND GROUPS: different colors/symbols per bar cluster (e.g., "PD" open
      circles, "PDD" filled circles). These are shown in the legend box.
  READ THE LEGEND AND AXIS LABELS VERY CAREFULLY — scientific labels often
  contain special characters like +, -, superscripts. Copy them exactly.

Step 2: CREATE ONE SERIES PER LEGEND GROUP (not per category).
  E.g., if legend shows PD and PDD, create series "PD" and series "PDD".

Step 3: FOR EACH SERIES, extract every individual data point.
  - x_values = repeat the category label for each dot in that category.
  - y_values = the y-coordinate of each individual dot (read from y-axis).
  Example for grouped bars with 3 categories and individual dots:
    Series "PD": x_values=["TH+M-","TH+M-","TH+M-","TH+M-","TH+M-",
                            "TH+M+","TH+M+","TH+M+","TH+M+","TH+M+",
                            "TH-M+","TH-M+","TH-M+","TH-M+","TH-M+"],
                 y_values=[128,48,25,15,12, 95,30,25,20,15, 120,115,25,18,15]
    Series "PDD": x_values=[...same categories...], y_values=[...dots...]

Step 4: PUT SUMMARY STATS IN NOTES (not as separate series).
  Format: "Bar heights (means): PD TH+M-=48, PD TH+M+=30, ...; Error bars: SEM"
  Also note n per group and significance annotations (* ** ns).

DISTINGUISHING OPEN vs FILLED CIRCLES:
  - Open circles (○) and filled circles (●) usually represent different groups
    (check the legend). Assign each dot to the correct series.
  - Dots that are horizontally jittered (spread left-right) within a bar are
    still the SAME x-category — the jitter is just for visibility.

- If NO individual dots are visible (just bars + error bars), then:
  y_values = [bar height], error_bars = [error bar extent], one entry per category.
- If grouped bars (multiple bars per category), create one series per legend group.
- NEGATIVE BARS: Report negative values as negative numbers.
- Report error bars if present (usually SEM or SD — note which in notes).
- LOG SCALE y-axis: read actual values from axis, not visual bar height.
- Significance annotations (* ** *** ns) → note in notes with group comparisons.""",

    PlotType.LINE: """Line plot guidance:
- Extract the value at each tick mark AND at each visible data point marker.
- If lines have no markers, sample at each x-axis tick mark AND at any
  inflection points where the curve changes direction or slope significantly.
- For curves that are steep in one region and flat in another, sample more
  densely in the steep region (e.g., every half tick interval).
- Create separate series for each line (use legend labels from the color key).
- DUAL Y-AXES: If the plot has two y-axes (left and right), each line uses
  its own axis. Check which axis each series belongs to by matching colors
  or by reading the legend. Report y_values using the correct axis scale.
  Note "dual y-axis: [left series] use left axis, [right series] use right axis".
- Report error bands/ribbons if present — use error_bars_lower and
  error_bars_upper for the band boundaries at each x-point.
- Report error bars (vertical bars at each point) the same way.
- LOG SCALE: If an axis is log-scaled, read actual values (see scatter guidance).
- Include at least 10 points per series for smooth curves.
- TIME SERIES: If x-axis is dates/time, use consistent format (e.g., "2020-01").""",

    PlotType.BOX: """Box and whisker plot guidance:
YOU MUST report plot_type as "box" in your output — not "bar".

For EACH group/category, create ONE series with:
  - name = the group label
  - x_values = [group_label]  (single entry)
  - y_values = [median]       (the middle line of the box)
  - error_bars_lower = [median − whisker_min]  (distance from median to bottom whisker)
  - error_bars_upper = [whisker_max − median]   (distance from median to top whisker)

In notes, report the FULL box statistics for each group:
  "GroupName: min=X, Q1=X, median=X, Q3=X, max=X"
If Q1 and Q3 are visible (the box edges), always include them.

OUTLIERS/FLIERS: If individual outlier points are shown beyond the whiskers,
  list them in notes: "GroupName outliers: [val1, val2, ...]"

- If individual points are overlaid (jittered/beeswarm), ALSO extract
  the raw data points as additional series named "GroupName (points)" with
  x_values = repeat the group label for each dot,
  y_values = individual y-axis values. Keep the box summary series too.
- x_values should be group/category labels.
- If boxes are horizontal, the same logic applies but read from x-axis.
- Significance annotations → note in notes with comparisons.""",

    PlotType.VIOLIN: """Violin plot guidance:
- If individual points are overlaid (jittered, beeswarm, or strip),
  PRIORITIZE the raw data points. Create ONE series per group with:
  x_values = repeat the group label for each dot,
  y_values = the individual y-axis values of each dot.
  Report violin summary stats in notes: "Control: median=12.5, Q1=10.2, Q3=14.8"
- If NO individual points are overlaid, extract summary statistics:
  - One series per group: y_values = [median]
  - error_bars_lower = [Q1], error_bars_upper = [Q3]
- If a box plot is drawn inside the violin, read values from that box.
- If only the median line is shown (no inner box), estimate Q1/Q3 from
  the violin shape where the width starts to narrow significantly.
- Note the approximate distribution shape (unimodal, bimodal, skewed) in notes.
- HALF-VIOLIN / RAINCLOUD plots: treat the dot strip and violin halves
  as a combined figure — extract individual points AND summary stats.
- x_values should be group/category labels.""",

    PlotType.FOREST: """Forest plot guidance:
- Each row is a study/subgroup. Use study names as x_values (strings).
- y_values = effect size (OR, HR, RR, SMD, ROM, or MD — note which in notes).
- ASYMMETRIC CI BOUNDS — read both ends of EACH horizontal line:
  Step 1: Read the position of the point/square/diamond on the x-axis (= estimate).
  Step 2: Read where the LEFT end of the CI line falls on the x-axis (= lower bound).
  Step 3: Read where the RIGHT end of the CI line falls on the x-axis (= upper bound).
  Step 4: error_bars_lower = estimate − lower_bound (always positive)
  Step 5: error_bars_upper = upper_bound − estimate (always positive)
  Example: estimate at 0.5, CI line goes from 0.2 to 0.85 →
    error_bars_lower = 0.5 − 0.2 = 0.3
    error_bars_upper = 0.85 − 0.5 = 0.35
- Include the overall/pooled estimate as the LAST entry (usually a diamond shape).
  For diamonds: y_value = center, error_bars = diamond left/right edges as deltas.
- SUBGROUP ANALYSES: If the forest plot has subgroup headings with their own
  subtotals, create separate series per subgroup. Name each series by the
  subgroup heading. Include subtotal rows within each series.
- Note heterogeneity statistics (I², τ², Q, p-value) if shown — per subgroup
  and overall.
- Note the number of events and total N per study if shown in columns.
- WEIGHT: If study weights (%) are shown, note them in notes.
- The x-axis label tells you the effect measure — report it in notes.
- If the scale is logarithmic (common for OR/HR/RR), read ACTUAL values
  from the axis (e.g., 0.5, 1, 2), not distances.
- Reference line: note the position of the dashed vertical line (usually at
  1.0 for ratios or 0 for differences).""",

    PlotType.KAPLAN_MEIER: """Kaplan-Meier survival curve guidance:
- The curve is a STEP FUNCTION — it drops vertically at event times and stays
  flat between events. Sample at EVERY visible step (drop point), not just ticks.
- x_values = time points, y_values = survival probability (0-1 or 0-100%).
- At minimum, extract: start point, every visible step, and the last point.
- Create separate series for each group/arm.
- CENSORING MARKS: Small vertical ticks on the curve indicate censored
  observations. Note approximate censoring times in notes if visible.
- NUMBER-AT-RISK TABLE: If shown below the plot, extract it in notes as
  a formatted table (group: time1=N1, time2=N2, ...).
- Note median survival time per group if marked (where curve crosses 0.5).
- If confidence bands (shading) are shown, create error_bars_lower and
  error_bars_upper series with the band boundaries at each time point.
- Note the log-rank p-value or HR with CI if displayed on the figure.
- If x-axis goes beyond the last event, don't extrapolate.""",

    PlotType.HISTOGRAM: """Histogram guidance:
- x_values = bin centers (or left edges if centers unclear), y_values = counts/frequencies.
- DENSITY PLOTS: If the y-axis shows density (area = 1) rather than counts,
  note "density" in notes. Sample the curve at each tick and inflection point.
- If a kernel density estimate (smooth curve) is overlaid, extract it as a
  separate series with dense sampling (15+ points).
- If multiple histograms overlap (semi-transparent), create separate series.
  Read the height of each group's bar carefully at each bin.
- Note bin width and total N in notes.
- CUMULATIVE: If the histogram is cumulative, note this in notes.
- FREQUENCY POLYGON: If shown as connected points instead of bars, treat
  like a line plot with x = bin centers, y = frequency.""",

    PlotType.HEATMAP: """Heatmap guidance:
- Extract the matrix of values if numbers are shown in cells.
- x_values = column labels, y_values = cell values (row by row).
- Create one series per row, with the row label as series name.
- If no numbers shown, estimate values from the color scale bar.
  Map each cell's color to the nearest value on the scale and note
  "values estimated from color scale" in notes.
- CLUSTERED HEATMAPS: If rows/columns are reordered by a dendrogram,
  use the reordered labels (as displayed), not the original order.
- CORRELATION HEATMAPS: If the matrix is symmetric (correlation matrix),
  you can extract just the lower triangle — note this in notes.
- Note the color scale range (e.g., -2 to 2, 0 to 100) in notes.""",

    PlotType.DOT_STRIP: """Dot/strip plot guidance:
- This shows individual data points for each group/condition (no bars).
- Create ONE series per group (use the group name from the legend or axis).
- For each series: x_values = repeat the group label for each dot,
  y_values = the individual y-axis values of each dot in that group.
  E.g., "Control" with 6 dots: x_values = ["Control"]*6, y_values = [12.3, 14.1, ...]
- If points are color-coded by a second variable (not the x-axis groups),
  create separate series per color instead, named by the color legend.
- BEESWARM: If points are spread horizontally to avoid overlap (beeswarm),
  still read the y-value — the x-jitter is just for visibility.
- If a horizontal line marks the mean or median for each group, report it
  in notes: "Group means: Control=12.9, Treatment=8.5" (NOT as a separate series).
- If error bars extend from the mean line, note those too in notes.
- Note the number of data points (n) per group in notes.
- Read each point carefully relative to y-axis ticks — do not cluster points
  at round numbers if they are actually spread between ticks.
- Significance annotations → note in notes.""",

    PlotType.STACKED_BAR: """Stacked bar chart guidance:
- x_values = category labels (strings).
- Create one series per stack segment (use legend labels for names).
- y_values = the ABSOLUTE value of each segment (not cumulative top).
  E.g., if a bar has segments 0-30, 30-55, 55-100, the values are 30, 25, 45.
- If the y-axis shows percentages (0-100%), note "proportions" in notes.
- Read segment boundaries carefully by projecting to y-axis ticks.
- If the chart shows counts, report raw counts; if percentages, report percentages.
- 100% STACKED: If all bars reach 100%, values are proportions. Check that
  segments within each bar sum to ~100%.
- DIVERGING STACKED: If bars extend in both directions from a center line
  (e.g., Likert scale), report negative segments as negative values.
- Note total bar height per category if visible.""",

    PlotType.FUNNEL: """Funnel plot guidance:
- This is a meta-analysis funnel plot (effect size vs. precision/SE).
- x_values = effect size (OR, RR, HR, SMD, or MD — note which in notes).
- y_values = standard error, precision, or sample size (note which on y-axis).
- Extract EVERY visible data point as (x, y).
- If the y-axis is inverted (0 at top, larger SE at bottom), preserve the
  actual values — do NOT flip them.
- Note any visible asymmetry or outliers.
- If trim-and-fill imputed studies are shown (often open/hollow circles),
  create a separate series "Imputed studies".
- CONTOUR FUNNEL: If shaded regions show significance contours (e.g.,
  p < 0.01, p < 0.05, p < 0.1), note the contour boundaries in notes.
- If a vertical line marks the pooled estimate, note its value.
- If Egger's regression line is shown, note its slope/intercept.""",

    PlotType.ROC: """ROC curve guidance:
- x_values = false positive rate (1 - specificity), typically 0 to 1.
- y_values = true positive rate (sensitivity), typically 0 to 1.
- Sample the curve densely: at EVERY visible data point marker, plus at
  each axis tick mark, plus at any inflection points. Aim for 15-20+ points.
- Create separate series for each classifier/model if multiple curves shown.
- Note the AUC (area under curve) value if displayed on the figure.
- If a diagonal reference line (chance line) is shown, do NOT extract it.
- If confidence bands are shown, use error_bars for the bounds.
- OPTIMAL CUTOFF: If a specific operating point is marked (e.g., Youden's J),
  note its coordinates and the corresponding threshold value.
- PRECISION-RECALL: If this is actually a precision-recall curve (y = precision,
  x = recall), note this and adjust axis labels accordingly.""",

    PlotType.VOLCANO: """Volcano plot guidance:
- x_values = log2 fold change (or similar effect size).
- y_values = -log10(p-value) or -log10(adjusted p-value). Check the axis label.
- There may be hundreds or thousands of points. Prioritize extraction:
  1. ALL labeled/annotated points (genes, proteins, etc.) — these are most
     important. Create a series per label or a single "Labeled" series.
  2. All points in the significant regions (colored differently).
  3. Representative non-significant points if feasible.
- If points are colored by significance threshold, create series:
  "Significant up", "Significant down", "Non-significant".
- Note the significance thresholds (fold-change cutoff lines, p-value
  cutoff line) and their values.
- If specific genes/proteins are labeled with text annotations, include
  each labeled point's name, x, and y in the notes or as a named series.
- MA PLOTS: If the x-axis is average expression (not fold change), this
  is an MA plot — note this and label axes accordingly.""",

    PlotType.WATERFALL: """Waterfall plot guidance:
- This shows ordered response values (e.g., % tumor change per patient).
- x_values = patient/sample identifiers or sequential indices (1, 2, 3...).
- y_values = response value (e.g., % change from baseline).
- Extract EVERY bar from left to right in order.
- If bars are colored by response category (e.g., responder vs non-responder,
  or by mutation status), create separate series per category.
- Note the response threshold lines if shown (e.g., -30% for partial response,
  +20% for progressive disease in RECIST criteria).
- Values above 0 = growth/increase; below 0 = shrinkage/decrease.
- If individual patient IDs are labeled, use those as x_values.
- Note total number of patients and % in each response category.""",

    PlotType.BLAND_ALTMAN: """Bland-Altman (difference) plot guidance:
- x_values = mean of the two measurements ((Method A + Method B) / 2).
- y_values = difference between methods (Method A - Method B).
- Extract EVERY visible data point as (x, y).
- Note the mean difference (bias) line value if shown (usually a solid
  horizontal line).
- Note the limits of agreement (mean ± 1.96 SD) if shown as dashed
  horizontal lines — report both upper and lower limit values.
- If confidence intervals around the bias or limits are shown (shaded
  bands), note those boundary values.
- Create separate series if points are grouped by a variable (e.g., color).
- PROPORTIONAL BIAS: If the scatter shows a trend (difference increasing
  with mean), note this.
- If a regression line through the points is shown, note its slope.""",

    PlotType.PAIRED: """Paired/spaghetti plot guidance:
- This shows before-after or repeated measurements connected by lines.
- x_values = condition/time-point labels (e.g., "Pre", "Post" or time points).
- Create one series per subject/sample if lines are individually identifiable
  and there are ≤15 subjects. Name them "Subject 1", "Subject 2", etc.
  (or use labels if shown). Each series has 2+ values (one per time point).
- If there are too many individual lines to distinguish (>15), instead:
  1. Create a "Mean" series with the mean at each time point.
  2. Note the approximate range (min, max at each time point) in notes.
  3. Note n (number of lines).
- If a bold/thick line shows the group mean, extract it as a "Mean" series.
- If individual points are also shown at each time point (dots on lines),
  read the values from the dots, not the connecting lines.
- Note whether most lines go up, down, or mixed.
- Significance annotations between time points → note in notes.""",

    PlotType.BUBBLE: """Bubble/balloon plot guidance:
- This is a scatter plot where point SIZE encodes a third variable.
- x_values = x-axis coordinates, y_values = y-axis coordinates.
- Create one series per group/color if groups are present.
- SIZE ENCODING: If a size legend is shown, estimate each bubble's value
  from the legend. Add a note "bubble_sizes: [v1, v2, ...]" in notes,
  matching the order of points in the series.
- If bubble sizes represent sample size, p-value, or another metric, note which.
- If no size legend is visible, estimate relative sizes (small/medium/large)
  and note this in notes.
- COLOR ENCODING: If bubble color represents a continuous variable (not just
  groups), note the color scale range and estimate values.
- Read bubble positions from their CENTER, not their edges.""",

    PlotType.AREA: """Area chart / filled line plot guidance:
- This is like a line plot but the area under the curve is filled.
- Extract the same way as line plots: x = tick marks/time points,
  y = values at each point.
- STACKED AREA: If multiple filled areas are stacked on top of each other,
  extract the ABSOLUTE value for each series (subtract the cumulative
  baseline). E.g., if bottom series goes 0-30 and top series 30-50,
  the top series value is 20.
- Create one series per colored area (use legend labels).
- 100% STACKED AREA: If y-axis is percentage (0-100%), the values are
  proportions. Check that all series sum to ~100% at each x-point.
- If confidence bands are shown as shaded areas around a line, extract
  the center line as y_values and band edges as error_bars.
- Sample at each x-axis tick mark plus any visible inflection points.""",

    PlotType.DOSE_RESPONSE: """Dose-response / sigmoidal curve guidance:
- This shows a sigmoidal (S-shaped) or log-linear relationship between
  dose/concentration (x-axis) and response (y-axis).
- x_values = dose or concentration values. The x-axis is often LOG-SCALED —
  read the actual values (e.g., 0.01, 0.1, 1, 10 µM), not log-distances.
- y_values = response (% inhibition, % viability, fold-change, etc.).
- Sample densely, especially in the transition region (steep part of the
  sigmoid). Aim for 15+ points per curve.
- If data points with error bars are shown (not just the fitted curve),
  extract THOSE points — they are the actual data. Report the fitted curve
  parameters in notes instead.
- Create separate series for each drug/compound/condition.
- EC50/IC50: Note the EC50/IC50 value if marked on the figure or annotated.
- Note Hill slope, top/bottom plateaus if shown in an inset or annotation.
- If the y-axis is normalized (0-100% or 0-1), note what 100% represents.""",

    PlotType.MANHATTAN: """Manhattan plot guidance:
- This shows -log10(p-value) across genomic positions, typically by chromosome.
- x_values = chromosomal position or SNP/gene identifiers.
- y_values = -log10(p-value).
- There are usually thousands of points — prioritize:
  1. ALL points ABOVE the genome-wide significance threshold line
     (typically -log10(5e-8) ≈ 7.3). These are the key findings.
  2. Points above the suggestive significance line if shown.
  3. Any specifically labeled/annotated SNPs or genes.
- Create one series for "Genome-wide significant" with labeled SNPs.
- Create a second series "Suggestive" for points between thresholds.
- For labeled peaks, x_value = gene/SNP name, y_value = -log10(p).
- Note the significance thresholds used (p-value and -log10 equivalent).
- Note chromosome numbers for the significant hits in notes.
- Do NOT try to extract every point — focus on significant findings.""",

    PlotType.CORRELATION_MATRIX: """Correlation matrix / correlogram guidance:
- This is a matrix showing pairwise correlations (typically -1 to +1).
- Create one series per ROW of the matrix, named by the row label.
- x_values = column labels, y_values = correlation coefficients.
- If numbers are shown in cells, extract those exact values.
- If only colors are shown, estimate from the color scale bar.
- SYMMETRIC MATRIX: If the matrix is symmetric (same variables on both axes),
  you can extract just the lower triangle — note this in notes.
- SIGNIFICANCE MARKERS: If cells are marked with * or × for significance,
  note which correlations are significant in notes.
- If cells are blank or crossed out (non-significant after correction),
  report null for those and note the correction method.
- Note the correlation method (Pearson, Spearman) if stated.""",

    PlotType.ERROR_BAR: """Error bar plot (means with error bars only) guidance:
- This shows group means as points/markers with error bars (confidence intervals),
  but WITHOUT bars underneath (unlike a bar chart).
- x_values = group/category labels or time points.
- y_values = mean/estimate values (the center point/marker).
- ASYMMETRIC ERROR BARS — read BOTH ends of each error bar INDEPENDENTLY:
  Step 1: Read the y-position of the point/marker (= mean).
  Step 2: Read where the BOTTOM of the error bar ends on the y-axis.
  Step 3: Read where the TOP of the error bar ends on the y-axis.
  Step 4: error_bars_lower = mean − bottom_end  (always positive)
  Step 5: error_bars_upper = top_end − mean  (always positive)
  These are often NOT equal! CIs from meta-analyses, ratios, and log-scale
  data are frequently asymmetric.
  Example: point at 0.5, CI bottom at 0.2, CI top at 0.85 →
    error_bars_lower = 0.5 − 0.2 = 0.3
    error_bars_upper = 0.85 − 0.5 = 0.35
- IMPORTANT: Read whether the error bars represent SD, SEM, 95% CI, or
  range — this is usually stated in the legend or y-axis label. Note it.
- If multiple groups are shown (different colors/markers), create one
  series per group.
- If this is a meta-analysis summary (ratio of means, odds ratio, etc.),
  note the dashed reference line value (often at 1.0 or 0) in notes.
- If lines connect the points (making it look like a line plot), it may
  overlap with the line plot type — treat as error_bar if the emphasis
  is on discrete group comparisons rather than continuous trends.
- If individual data points are overlaid, extract those as a separate series.
- Significance annotations (p-values, brackets) → note in notes.""",

    PlotType.TABLE: """Table guidance:
- TITLE: The "title" field must be a DESCRIPTIVE title derived from the figure
  legend or table caption — e.g., "Proteasomal subunit levels in substantia
  nigra". NEVER use generic identifiers like "Table 1", "Table 17", or
  "Data table". If the legend says "Table 3. Demographic and clinical
  characteristics", the title should be "Demographic and clinical characteristics".
- This is a data table from a scientific paper. Extract ALL numerical data.
- Create one series per data column that contains numerical values.
- x_values = row labels (first column, or row identifiers like group names,
  brain regions, etc.).
- y_values = the numerical values in that column.
- If a column contains mean ± SD/SEM, split into y_values (mean) and
  error_bars (the ± value). Note whether it's SD, SEM, or CI.
- Use the column header as the series name.
- If the table has sub-headers or grouped rows, note the grouping in notes.
- For p-values, sample sizes (n), or statistical columns, include them in
  notes rather than as series (e.g., "p-values: Group A vs B = 0.003").
- MERGED CELLS: If cells span multiple rows/columns, repeat the value for
  each row/column it covers.
- Extract EVERY row and column — do not skip any data.""",
}

DEFAULT_GUIDANCE = """General guidance:
- First IDENTIFY the plot type from the image. Set plot_type accurately.
- Extract all visible numerical data systematically.
- Use your best judgment for the data structure.
- LOG SCALE: If any axis is log-scaled, report actual values (not log-transformed).
  A point halfway between 10 and 100 on a log axis ≈ 30, not 55.
- DUAL Y-AXES: If two y-axes exist, match each series to its correct axis.
- Note any unusual features, significance annotations, or ambiguities.

BOX AND WHISKER PLOTS — if you see a box plot, you MUST:
  - Set plot_type to "box" (NOT "bar").
  - Create ONE series per group/category with:
    name = group label, x_values = [group_label], y_values = [median].
    error_bars_lower = [median − whisker_min], error_bars_upper = [whisker_max − median].
  - In notes, report: "GroupName: min=X, Q1=X, median=X, Q3=X, max=X" for each group.
  - Outliers beyond whiskers: "GroupName outliers: [val1, val2, ...]" in notes.

VIOLIN PLOTS — if you see a violin plot with overlaid points:
  - Set plot_type to "violin". Prioritize extracting the individual data points.

BAR CHARTS WITH INDIVIDUAL DATA POINTS — if you see bars with dots/circles overlaid:
  - The individual data points (dots) are the REAL DATA — extract those, not bar heights.
  - Check the legend: open circles (○) vs filled circles (●) usually = different groups.
  - Create ONE series per LEGEND GROUP (e.g., "PD", "PDD"), not per category.
  - x_values = repeat category label for each dot; y_values = each dot's y-coordinate.
  - Put bar heights (means) and error bar type in notes.
  - only use plot_type "bar" for actual rectangular bars, never for box plots."""


def _extract_panel_legend(full_legend: str, panel_label: str) -> str | None:
    """Try to extract a panel-specific description from the full figure legend.

    E.g., from "Fig 3. (a) Precision across replicates. (b) Significant associations."
    extract "Precision across replicates" for panel "a".
    """
    if not full_legend or not panel_label or panel_label == "main":
        return None

    # Match patterns like "(a) ...", "(a,b) ...", "a, ...", "a: ..."
    label = re.escape(panel_label.lower())
    pattern = rf"\({label}[,\)][^)]*\)\s*(.+?)(?=\([a-z](?:[,\)])|$)"
    match = re.search(pattern, full_legend, re.IGNORECASE | re.DOTALL)
    if match:
        return match.group(1).strip().rstrip(".")

    # Simpler pattern: "a, Description" or "a: Description"
    pattern2 = rf"(?:^|\s){label}[,:]\s*(.+?)(?=\s+[b-z][,:]|\Z)"
    match2 = re.search(pattern2, full_legend, re.IGNORECASE | re.DOTALL)
    if match2:
        return match2.group(1).strip().rstrip(".")

    return None


def _get_prompt(figure: Figure, panel_label: str | None = None, pre_analysis: dict | None = None) -> str:
    legend_context = ""
    if figure.legend:
        legend_context = f"Figure legend: {figure.legend}"
        if panel_label:
            panel_desc = _extract_panel_legend(figure.legend, panel_label)
            if panel_desc:
                legend_context += f"\n\nThis is panel ({panel_label}): {panel_desc}"
            else:
                legend_context += f"\n\nThis is panel ({panel_label}) of the figure."

    panel_focus = ""
    if panel_label and panel_label != "main":
        panel_focus = (
            "CONTEXT: This is a CROPPED sub-panel (panel " + panel_label + ") from a larger "
            "multi-panel figure. A full-figure context image was provided above.\n"
            "- Use the FULL FIGURE to read: color legend labels, axis scales from adjacent "
            "panels that share the same y-axis, and any shared legends.\n"
            "- Extract data ONLY from panel " + panel_label + " in the cropped image.\n"
            "- If the y-axis tick marks are not visible in the crop, look at the full figure "
            "to find the axis scale for this panel's row.\n"
            "- The crop may show edges of adjacent panels — ignore those completely."
        )

    guidance = PLOT_GUIDANCE.get(figure.plot_type, DEFAULT_GUIDANCE)

    # If pre-analysis was done, inject its findings as context
    pre_analysis_context = ""
    if pre_analysis:
        parts = ["PRE-ANALYSIS RESULTS (use these to guide your extraction):"]
        pt = pre_analysis.get("plot_type")
        if pt:
            parts.append(f"- Plot type identified as: {pt}")
        x_ax = pre_analysis.get("x_axis", {})
        y_ax = pre_analysis.get("y_axis", {})
        if x_ax.get("scale") == "log":
            parts.append(f"- X-axis is LOG SCALED (ticks: {x_ax.get('tick_values', [])})")
        if y_ax.get("scale") == "log":
            parts.append(f"- Y-axis is LOG SCALED (ticks: {y_ax.get('tick_values', [])})")
        if x_ax.get("label"):
            parts.append(f"- X-axis label: {x_ax['label']}" + (f" ({x_ax['unit']})" if x_ax.get("unit") else ""))
        if y_ax.get("label"):
            parts.append(f"- Y-axis label: {y_ax['label']}" + (f" ({y_ax['unit']})" if y_ax.get("unit") else ""))
        legend_entries = pre_analysis.get("legend_entries", [])
        if legend_entries:
            parts.append(f"- Legend entries: {', '.join(legend_entries)}")
        desc = pre_analysis.get("description")
        if desc:
            parts.append(f"- Description: {desc}")
        if pre_analysis.get("has_error_bars"):
            parts.append("- Error bars are PRESENT — extract them carefully")
        else:
            parts.append("- No error bars detected — OMIT error_bars_lower/upper keys")
        pre_analysis_context = "\n".join(parts)

    return BASE_PROMPT.format(
        legend_context=legend_context,
        plot_specific_guidance=guidance,
        panel_focus=panel_focus,
    ) + ("\n\n" + pre_analysis_context if pre_analysis_context else "")


def _parse_response(raw_text: str, figure_id: str) -> ExtractedData:
    """Parse Claude's JSON response into ExtractedData."""
    # Strip markdown code fences
    text = raw_text.strip()
    if text.startswith("```"):
        text = text.split("\n", 1)[1]
    if text.endswith("```"):
        text = text[:-3]
    text = text.strip()

    data = json.loads(text)

    series = []
    for s in data.get("series", []):
        y_vals = s.get("y_values") or []
        x_vals = s.get("x_values") or []
        err_lo = s.get("error_bars_lower") or []
        err_hi = s.get("error_bars_upper") or []

        # Pad error arrays to match y length
        while len(err_lo) < len(y_vals):
            err_lo.append(None)
        while len(err_hi) < len(y_vals):
            err_hi.append(None)

        series.append(DataSeries(
            name=s.get("name", "Series"),
            x_values=x_vals,
            y_values=y_vals,
            error_bars_lower=err_lo,
            error_bars_upper=err_hi,
        ))

    try:
        plot_type = PlotType(data.get("plot_type", "other"))
    except ValueError:
        plot_type = PlotType.OTHER

    try:
        confidence = Confidence(data.get("confidence", "medium"))
    except ValueError:
        confidence = Confidence.MEDIUM

    return ExtractedData(
        figure_id=figure_id,
        plot_type=plot_type,
        title=data.get("title"),
        x_label=data.get("x_label"),
        y_label=data.get("y_label"),
        x_unit=data.get("x_unit"),
        y_unit=data.get("y_unit"),
        x_min=data.get("x_min"),
        x_max=data.get("x_max"),
        y_min=data.get("y_min"),
        y_max=data.get("y_max"),
        x_scale=data.get("x_scale"),
        y_scale=data.get("y_scale"),
        series=series,
        confidence=confidence,
        notes=data.get("notes"),
    )


PANEL_DETECT_PROMPT = """\
This scientific figure may contain multiple sub-panels (a, b, c, d...) arranged in a grid.

Your task: describe the GRID LAYOUT, then assign each panel to a grid cell.

STEP 1 — Determine the grid:
- How many rows and columns of panels are there?
- A 2×2 figure has 2 rows and 2 columns. A figure with 3 panels on top and 1 wide \
panel on the bottom is 2 rows × 3 columns (bottom panel spans all 3 columns).
- Count rows by distinct vertical positions of panels. Count columns by distinct \
horizontal positions.

STEP 2 — Assign each panel to its grid cell(s):
- row and col are 1-indexed (top-left = row 1, col 1).
- If a panel spans multiple cells, use rowspan/colspan.

Return a JSON object:
{{
  "rows": 2,
  "cols": 2,
  "panels": [
    {{"label": "a", "plot_type": "bar", "row": 1, "col": 1, "rowspan": 1, "colspan": 1}},
    {{"label": "b", "plot_type": "scatter", "row": 1, "col": 2, "rowspan": 1, "colspan": 1}},
    {{"label": "c", "plot_type": "line", "row": 2, "col": 1, "rowspan": 1, "colspan": 2}}
  ]
}}

If the figure is a single panel (not multi-panel), return:
{{"rows": 1, "cols": 1, "panels": [{{"label": "main", "plot_type": "<type>", "row": 1, "col": 1, "rowspan": 1, "colspan": 1}}]}}

Plot type must be one of: scatter, bar, line, box, violin, histogram, heatmap, \
forest, kaplan_meier, dot_strip, stacked_bar, funnel, roc, volcano, waterfall, \
bland_altman, paired, bubble, area, dose_response, manhattan, correlation_matrix, \
error_bar, table, other.
Return ONLY the JSON object, no other text."""


_PLOT_TYPE_ENUM = [
    "scatter", "bar", "line", "box", "violin", "histogram", "heatmap",
    "forest", "kaplan_meier", "dot_strip", "stacked_bar", "funnel", "roc",
    "volcano", "waterfall", "bland_altman", "paired", "bubble", "area",
    "dose_response", "manhattan", "correlation_matrix", "error_bar", "table", "other",
]

# ── Tool definitions for structured Claude API calls ──

PRE_ANALYSIS_TOOL = {
    "name": "report_figure_structure",
    "description": (
        "Report the structure of a scientific figure: plot type, axis info, "
        "scale (linear/log), legend entries, and whether error bars are present."
    ),
    "input_schema": {
        "type": "object",
        "required": ["plot_type", "x_axis", "y_axis", "legend_entries", "num_series",
                      "has_error_bars", "description"],
        "properties": {
            "plot_type": {
                "type": "string",
                "enum": _PLOT_TYPE_ENUM,
                "description": "The type of plot shown in the figure.",
            },
            "x_axis": {
                "type": "object",
                "required": ["label", "scale", "is_categorical"],
                "properties": {
                    "label": {"type": "string", "description": "X-axis label text"},
                    "unit": {"type": ["string", "null"], "description": "X-axis unit if shown"},
                    "scale": {"type": "string", "enum": ["linear", "log"],
                              "description": "linear or log. Ticks like 0.01,0.1,1,10,100 = log."},
                    "tick_values": {
                        "type": "array",
                        "items": {"type": ["number", "string"]},
                        "description": "Visible tick mark values on the x-axis",
                    },
                    "is_categorical": {"type": "boolean", "description": "True if x-axis shows categories/labels"},
                },
            },
            "y_axis": {
                "type": "object",
                "required": ["label", "scale"],
                "properties": {
                    "label": {"type": "string", "description": "Y-axis label text"},
                    "unit": {"type": ["string", "null"], "description": "Y-axis unit if shown"},
                    "scale": {"type": "string", "enum": ["linear", "log"],
                              "description": "linear or log. Ticks like 0.01,0.1,1,10,100 = log."},
                    "tick_values": {
                        "type": "array",
                        "items": {"type": ["number", "string"]},
                        "description": "Visible tick mark values on the y-axis",
                    },
                },
            },
            "legend_entries": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Series/group names from the color/shape legend. Empty if no legend visible.",
            },
            "num_series": {
                "type": "integer",
                "description": "Number of distinct data series/groups in the figure",
            },
            "has_error_bars": {
                "type": "boolean",
                "description": "True if any error bars, confidence intervals, or uncertainty indicators are visible",
            },
            "description": {
                "type": "string",
                "description": "One sentence describing what this plot shows",
            },
        },
    },
}

EXTRACT_DATA_TOOL = {
    "name": "report_extracted_data",
    "description": "Report all numerical data extracted from a scientific figure.",
    "input_schema": {
        "type": "object",
        "required": ["plot_type", "title", "x_label", "y_label", "x_scale", "y_scale",
                      "series", "confidence"],
        "properties": {
            "plot_type": {"type": "string", "enum": _PLOT_TYPE_ENUM},
            "title": {
                "type": "string",
                "description": "Descriptive title from the legend/caption — never just 'Figure 1' or 'Table 2'.",
            },
            "x_label": {"type": ["string", "null"]},
            "y_label": {"type": ["string", "null"]},
            "x_unit": {"type": ["string", "null"]},
            "y_unit": {"type": ["string", "null"]},
            "x_min": {"type": ["number", "null"], "description": "First tick value on x-axis, null if categorical"},
            "x_max": {"type": ["number", "null"], "description": "Last tick value on x-axis, null if categorical"},
            "y_min": {"type": ["number", "null"], "description": "First tick value on y-axis"},
            "y_max": {"type": ["number", "null"], "description": "Last tick value on y-axis"},
            "x_scale": {"type": "string", "enum": ["linear", "log"]},
            "y_scale": {"type": "string", "enum": ["linear", "log"]},
            "series": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["name", "x_values", "y_values"],
                    "properties": {
                        "name": {"type": "string", "description": "Series/group name from legend — use real names"},
                        "x_values": {"type": "array", "items": {"type": ["number", "string"]}},
                        "y_values": {"type": "array", "items": {"type": ["number", "null"]}},
                        "error_bars_lower": {
                            "type": "array",
                            "items": {"type": ["number", "null"]},
                            "description": "Positive distance BELOW mean. Omit entirely if no error bars.",
                        },
                        "error_bars_upper": {
                            "type": "array",
                            "items": {"type": ["number", "null"]},
                            "description": "Positive distance ABOVE mean. Omit entirely if no error bars.",
                        },
                    },
                },
            },
            "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
            "notes": {
                "type": ["string", "null"],
                "description": (
                    "Error bar type (SD/SEM/CI), sample sizes, box plot stats "
                    "(GroupName: min=X, Q1=X, median=X, Q3=X, max=X), outliers, ambiguities."
                ),
            },
        },
    },
}

PRE_ANALYSIS_SYSTEM = """\
You are analyzing a scientific figure BEFORE extracting numerical data.
Your job is to identify the structure so the extraction step knows exactly what to do.

You MUST call the report_figure_structure tool with your analysis.

CRITICAL for axis scale detection:
- If tick values increase multiplicatively (0.01, 0.1, 1, 10, 100) → "log"
- If tick values increase additively (0, 20, 40, 60, 80) → "linear"
- Forest plots with OR/HR often have log x-axis
- Dose-response curves often have log x-axis
- Box-and-whisker plots have plot_type "box" — NOT "bar"

Read every tick mark, every legend entry, and every axis label carefully."""


def _image_from_base64(b64: str) -> Image.Image:
    return Image.open(io.BytesIO(base64.b64decode(b64)))


def _image_to_base64(img: Image.Image) -> str:
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()


def _extract_tool_input(response, tool_name: str) -> dict | None:
    """Extract the input dict from a tool_use content block in a Claude response."""
    for block in response.content:
        if block.type == "tool_use" and block.name == tool_name:
            return block.input
    return None


async def pre_analyze(figure: Figure) -> dict:
    """Phase 2: Analyze figure structure before data extraction.

    Uses Claude tool calling for structured, reliable output.
    Returns a dict with plot_type, axis info, scale, legend entries, etc.
    """
    import anthropic
    client = anthropic.AsyncAnthropic()

    response = await client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=2048,
        system=PRE_ANALYSIS_SYSTEM,
        tools=[PRE_ANALYSIS_TOOL],
        tool_choice={"type": "tool", "name": "report_figure_structure"},
        messages=[{
            "role": "user",
            "content": [
                {
                    "type": "image",
                    "source": {
                        "type": "base64",
                        "media_type": "image/png",
                        "data": figure.image_base64,
                    },
                },
                {"type": "text", "text": (
                    "Analyze this scientific figure and report its structure. "
                    "Pay close attention to:\n"
                    "- The EXACT plot type (box vs bar, scatter vs dot_strip, etc.)\n"
                    "- Whether axes are LOG-scaled (look for exponential tick spacing: 1, 10, 100, 1000)\n"
                    "- All legend entries and their visual encoding (color, shape, fill)\n"
                    "- Whether individual data points are overlaid on bars/boxes\n"
                    "- Whether error bars are present and what kind they appear to be"
                )},
            ],
        }],
    )

    analysis = _extract_tool_input(response, "report_figure_structure") or {}

    logger.info(
        f"Pre-analysis for {figure.figure_id}: "
        f"type={analysis.get('plot_type', '?')}, "
        f"x_scale={analysis.get('x_axis', {}).get('scale', '?')}, "
        f"y_scale={analysis.get('y_axis', {}).get('scale', '?')}, "
        f"series={analysis.get('num_series', '?')}, "
        f"legend={analysis.get('legend_entries', [])}"
    )

    return analysis


async def detect_panels(figure: Figure) -> list[dict]:
    """Ask Claude to identify sub-panels via grid layout, then compute bboxes.

    Uses a grid-based approach: Claude reports rows/cols and each panel's cell
    position, then bboxes are computed mathematically from the grid. This is
    far more reliable than asking Claude to estimate pixel coordinates or
    percentages, which it consistently gets wrong.
    """
    import anthropic
    client = anthropic.AsyncAnthropic()

    # Use actual image dimensions (may differ from figure.width/height if resized)
    img = _image_from_base64(figure.image_base64)
    actual_w, actual_h = img.width, img.height

    response = await client.messages.create(
        model="claude-opus-4-6",
        max_tokens=2048,
        messages=[{
            "role": "user",
            "content": [
                {
                    "type": "image",
                    "source": {
                        "type": "base64",
                        "media_type": "image/png",
                        "data": figure.image_base64,
                    },
                },
                {"type": "text", "text": PANEL_DETECT_PROMPT},
            ],
        }],
    )

    raw = response.content[0].text.strip()
    if raw.startswith("```"):
        raw = raw.split("\n", 1)[1]
    if raw.endswith("```"):
        raw = raw[:-3]
    raw = raw.strip()

    data = json.loads(raw)
    panels_info = data["panels"]
    nrows = data.get("rows", 1)
    ncols = data.get("cols", 1)

    logger.info(
        f"Detected {len(panels_info)} panels in {figure.figure_id} "
        f"(grid: {nrows}×{ncols})"
    )

    # Compute bboxes mathematically from grid positions
    cell_w = actual_w / ncols
    cell_h = actual_h / nrows

    panels = []
    for p in panels_info:
        row = p.get("row", 1)
        col = p.get("col", 1)
        rowspan = p.get("rowspan", 1)
        colspan = p.get("colspan", 1)

        # Grid cell to pixel bbox (0-indexed internally)
        x_min = int((col - 1) * cell_w)
        y_min = int((row - 1) * cell_h)
        x_max = int((col - 1 + colspan) * cell_w)
        y_max = int((row - 1 + rowspan) * cell_h)

        bbox = [x_min, y_min, x_max, y_max]
        panels.append({
            "label": p["label"],
            "plot_type": p["plot_type"],
            "bbox": bbox,
        })

    return panels


def _crop_panel(figure: Figure, bbox: list[int], label: str, plot_type_str: str, padding: int = 40) -> Figure:
    """Crop a sub-panel from the full figure image and return as a new Figure."""
    img = _image_from_base64(figure.image_base64)
    x_min, y_min, x_max, y_max = bbox

    # Scale bbox if figure dimensions don't match actual image
    # (can happen if width/height fields weren't updated after resizing)
    if figure.width != img.width or figure.height != img.height:
        sx = img.width / max(figure.width, 1)
        sy = img.height / max(figure.height, 1)
        x_min = int(x_min * sx)
        y_min = int(y_min * sy)
        x_max = int(x_max * sx)
        y_max = int(y_max * sy)
        logger.info(f"Scaled bbox for {label}: figure={figure.width}x{figure.height} img={img.width}x{img.height}")

    # Add padding and clamp to image bounds
    x_min = max(0, x_min - padding)
    y_min = max(0, y_min - padding)
    x_max = min(img.width, x_max + padding)
    y_max = min(img.height, y_max + padding)

    cropped = img.crop((x_min, y_min, x_max, y_max))
    b64 = _image_to_base64(cropped)

    try:
        pt = PlotType(plot_type_str)
    except ValueError:
        pt = PlotType.OTHER

    panel_id = f"{figure.figure_id}_{label}"
    return Figure(
        figure_id=panel_id,
        paper_id=figure.paper_id,
        page_number=figure.page_number,
        image_index=figure.image_index,
        width=cropped.width,
        height=cropped.height,
        image_base64=b64,
        legend=figure.legend,
        plot_type=pt,
        plot_type_confidence=Confidence.MEDIUM,
    )


async def _digitize_single(
    figure: Figure,
    panel_label: str | None = None,
    full_figure_b64: str | None = None,
    pre_analysis_result: dict | None = None,
) -> tuple[ExtractedData, dict]:
    """Phase 2+3: Pre-analyze then extract data from a single-panel figure image.

    Returns (result, pre_analysis_dict) so validation can cross-check.

    If full_figure_b64 is provided (for multi-panel splits), it is sent as
    a context image so the model can reference the full figure's axes, legends,
    and labels even when the cropped panel is missing them.
    """
    import anthropic
    client = anthropic.AsyncAnthropic()

    # Phase 2: Pre-analyze if not already done
    if pre_analysis_result is None:
        try:
            pre_analysis_result = await pre_analyze(figure)
        except Exception as e:
            logger.warning(f"Pre-analysis failed for {figure.figure_id}: {e}")
            pre_analysis_result = {}

    # Update figure's plot_type from pre-analysis if it was OTHER
    if figure.plot_type == PlotType.OTHER and pre_analysis_result.get("plot_type"):
        try:
            detected_type = PlotType(pre_analysis_result["plot_type"])
            figure = figure.model_copy(update={"plot_type": detected_type})
            logger.info(f"Updated {figure.figure_id} plot_type to {detected_type.value} from pre-analysis")
        except ValueError:
            pass

    prompt = _get_prompt(figure, panel_label=panel_label, pre_analysis=pre_analysis_result)

    # Phase 2b: CV calibration — detect markers/bars at pixel level
    cv_context = ""
    try:
        from .cv_calibration import calibrate_image, format_calibration_prompt
        img_bytes = base64.b64decode(figure.image_base64)
        cal = calibrate_image(img_bytes)
        cv_context = format_calibration_prompt(cal)
        if cv_context:
            logger.info(
                f"CV calibration for {figure.figure_id}: "
                f"{len(cal.markers)} markers, {len(cal.bars)} bars"
            )
    except Exception as e:
        logger.warning(f"CV calibration failed for {figure.figure_id}: {e}")

    content: list[dict] = []

    # For multi-panel: send full figure first as context, then the cropped panel
    if full_figure_b64:
        content.append({"type": "text", "text": (
            "Here is the FULL multi-panel figure for context. "
            "Use it to read axis scales, tick values, color legends, "
            "and series labels. Then extract data ONLY from the cropped "
            "panel shown in the second image."
        )})
        content.append({
            "type": "image",
            "source": {
                "type": "base64",
                "media_type": "image/png",
                "data": full_figure_b64,
            },
        })
        content.append({"type": "text", "text": "Now extract data from THIS cropped panel:"})

    content.append({
        "type": "image",
        "source": {
            "type": "base64",
            "media_type": "image/png",
            "data": figure.image_base64,
        },
    })
    full_prompt = prompt
    if cv_context:
        full_prompt += "\n\n" + cv_context
    content.append({"type": "text", "text": full_prompt + "\n\nYou MUST call the report_extracted_data tool with your results."})

    response = await client.messages.create(
        model="claude-opus-4-6",
        max_tokens=8192,
        tools=[EXTRACT_DATA_TOOL],
        tool_choice={"type": "tool", "name": "report_extracted_data"},
        messages=[{"role": "user", "content": content}],
    )

    # Extract data from tool call (structured, no JSON parsing needed)
    data = _extract_tool_input(response, "report_extracted_data")

    if data:
        result = _parse_tool_result(data, figure.figure_id)
    else:
        # Fallback: try parsing raw text (shouldn't happen with tool_choice=tool)
        raw_text = ""
        for block in response.content:
            if hasattr(block, "text"):
                raw_text = block.text
                break
        logger.warning(f"Tool call not found for {figure.figure_id}, falling back to text parse")
        result = _parse_response(raw_text, figure.figure_id)

    # Attach legend and context metadata
    result.legend_text = figure.legend
    # text_mentions are attached by the caller (needs PDF access)

    return result, pre_analysis_result


def _parse_tool_result(data: dict, figure_id: str) -> ExtractedData:
    """Parse structured tool call output into ExtractedData.

    Unlike _parse_response which deals with raw JSON text, this receives
    already-parsed dict from Claude's tool_use response — no string parsing needed.
    """
    series = []
    for s in data.get("series", []):
        y_vals = s.get("y_values") or []
        x_vals = s.get("x_values") or []
        err_lo = s.get("error_bars_lower") or []
        err_hi = s.get("error_bars_upper") or []

        while len(err_lo) < len(y_vals):
            err_lo.append(None)
        while len(err_hi) < len(y_vals):
            err_hi.append(None)

        series.append(DataSeries(
            name=s.get("name", "Series"),
            x_values=x_vals,
            y_values=[v if v is not None else 0.0 for v in y_vals],
            error_bars_lower=err_lo,
            error_bars_upper=err_hi,
        ))

    try:
        plot_type = PlotType(data.get("plot_type", "other"))
    except ValueError:
        plot_type = PlotType.OTHER

    try:
        confidence = Confidence(data.get("confidence", "medium"))
    except ValueError:
        confidence = Confidence.MEDIUM

    return ExtractedData(
        figure_id=figure_id,
        plot_type=plot_type,
        title=data.get("title"),
        x_label=data.get("x_label"),
        y_label=data.get("y_label"),
        x_unit=data.get("x_unit"),
        y_unit=data.get("y_unit"),
        x_min=data.get("x_min"),
        x_max=data.get("x_max"),
        y_min=data.get("y_min"),
        y_max=data.get("y_max"),
        x_scale=data.get("x_scale"),
        y_scale=data.get("y_scale"),
        series=series,
        confidence=confidence,
        notes=data.get("notes"),
    )


def validate_extraction(result: ExtractedData, pre_analysis: dict | None = None) -> ExtractedData:
    """Lightweight heuristic validation — no API calls, no cost, instant.

    Checks for common extraction errors and flags them in notes.
    Returns the result with warnings appended to notes.
    """
    warnings: list[str] = []

    # 1. Empty series check
    if not result.series:
        warnings.append("No series extracted")
    else:
        for s in result.series:
            if not s.y_values:
                warnings.append(f"Series '{s.name}' has no y_values")
            if s.x_values and s.y_values and len(s.x_values) != len(s.y_values):
                warnings.append(f"Series '{s.name}': x_values ({len(s.x_values)}) and y_values ({len(s.y_values)}) length mismatch")

    # 2. Values within axis range
    if result.y_min is not None and result.y_max is not None:
        y_range = result.y_max - result.y_min
        tolerance = y_range * 0.15  # 15% tolerance for points near edges
        for s in result.series:
            for i, y in enumerate(s.y_values):
                if y is not None:
                    if y < result.y_min - tolerance or y > result.y_max + tolerance:
                        warnings.append(f"Series '{s.name}' point {i}: y={y} outside axis range [{result.y_min}, {result.y_max}]")
                        break  # one warning per series is enough

    # 3. Plot type consistency with pre-analysis
    if pre_analysis:
        pa_type = pre_analysis.get("plot_type")
        if pa_type and pa_type != "other" and pa_type != result.plot_type.value:
            warnings.append(f"Plot type mismatch: pre-analysis={pa_type}, extraction={result.plot_type.value}")

        # Series count check
        pa_n = pre_analysis.get("num_series")
        if pa_n and isinstance(pa_n, int) and len(result.series) > 0:
            if abs(len(result.series) - pa_n) > pa_n:  # more than 2x off
                warnings.append(f"Series count: pre-analysis expected ~{pa_n}, got {len(result.series)}")

        # Legend entry check
        pa_legend = pre_analysis.get("legend_entries", [])
        if pa_legend and len(pa_legend) > 0 and len(result.series) > 0:
            extracted_names = {s.name.lower().strip() for s in result.series}
            legend_names = {str(e).lower().strip() for e in pa_legend}
            missing = legend_names - extracted_names
            if missing and len(missing) <= 3:  # don't warn if totally different (may be subgroups)
                warnings.append(f"Legend entries not found in series: {', '.join(missing)}")

    # 4. Box plot sanity
    if result.plot_type == PlotType.BOX:
        for s in result.series:
            if s.y_values and len(s.y_values) == 1:
                # Check notes for box stats
                if result.notes:
                    import re
                    pat = re.escape(s.name) + r".*?min\s*=\s*([\d.\-]+).*?Q1\s*=\s*([\d.\-]+).*?median\s*=\s*([\d.\-]+).*?Q3\s*=\s*([\d.\-]+).*?max\s*=\s*([\d.\-]+)"
                    m = re.search(pat, result.notes, re.IGNORECASE)
                    if m:
                        mn, q1, med, q3, mx = [float(x) for x in m.groups()]
                        if not (mn <= q1 <= med <= q3 <= mx):
                            warnings.append(f"Box '{s.name}': stats not monotonic (min={mn} Q1={q1} median={med} Q3={q3} max={mx})")

    # 5. Error bar sanity (should be positive extents)
    for s in result.series:
        for i, v in enumerate(s.error_bars_lower):
            if v is not None and v < 0:
                warnings.append(f"Series '{s.name}': negative error_bars_lower at index {i} (should be positive extent)")
                break
        for i, v in enumerate(s.error_bars_upper):
            if v is not None and v < 0:
                warnings.append(f"Series '{s.name}': negative error_bars_upper at index {i} (should be positive extent)")
                break

    # 6. Duplicate series names
    names = [s.name for s in result.series]
    if len(names) != len(set(names)):
        from collections import Counter
        dupes = [n for n, c in Counter(names).items() if c > 1]
        warnings.append(f"Duplicate series names: {', '.join(dupes)}")

    # 7. Confidence downgrade if warnings found
    if warnings:
        note_text = "[Validation] " + "; ".join(warnings)
        existing_notes = result.notes or ""
        updates: dict = {"notes": existing_notes + "\n" + note_text if existing_notes else note_text}

        # Downgrade confidence if many warnings
        if len(warnings) >= 3 and result.confidence != Confidence.LOW:
            updates["confidence"] = Confidence.LOW
        elif len(warnings) >= 1 and result.confidence == Confidence.HIGH:
            updates["confidence"] = Confidence.MEDIUM

        result = result.model_copy(update=updates)
        logger.info(f"Validation for {result.figure_id}: {len(warnings)} warning(s): {'; '.join(warnings)}")

    return result


def _render_result_to_png(result: ExtractedData) -> str:
    """Render ExtractedData to a matplotlib PNG, return base64.

    This creates a plot that can be visually compared to the original figure.
    Handles all major plot types: box, bar, stacked_bar, histogram, violin,
    forest, kaplan_meier, scatter, line, dot_strip, paired, error_bar,
    waterfall, roc, heatmap, volcano, funnel, dose_response, area, bland_altman.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyBboxPatch
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for verification rendering")
        return ""

    try:
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = ["#3b82f6", "#ef4444", "#10b981", "#f59e0b", "#8b5cf6",
                  "#ec4899", "#14b8a6", "#f97316", "#6366f1", "#84cc16"]
        pt = result.plot_type

        # --- Helper: get unique ordered categories across all series ---
        def _unique_cats() -> list[str]:
            cats: list[str] = []
            seen: set[str] = set()
            for s in result.series:
                for v in s.x_values:
                    k = str(v)
                    if k not in seen:
                        seen.add(k)
                        cats.append(k)
            return cats

        # --- Helper: safe numeric x/y pairs (skips non-numeric x_values) ---
        def _xy_numeric(s):
            x, y = [], []
            for i, xv in enumerate(s.x_values):
                if xv is None:
                    continue
                try:
                    xf = float(xv)
                except (ValueError, TypeError):
                    continue
                if i < len(s.y_values) and s.y_values[i] is not None:
                    x.append(xf)
                    y.append(s.y_values[i])
            return x, y

        # --- Helper: build error arrays for a series against a category list ---
        def _err_arrays(s, unique_cats):
            vals, elo, ehi = [], [], []
            for cat in unique_cats:
                idx = next((i for i, x in enumerate(s.x_values) if str(x) == cat), None)
                vals.append(s.y_values[idx] if idx is not None and idx < len(s.y_values) else 0)
                elo.append(s.error_bars_lower[idx] if idx is not None and idx < len(s.error_bars_lower) and s.error_bars_lower[idx] else 0)
                ehi.append(s.error_bars_upper[idx] if idx is not None and idx < len(s.error_bars_upper) and s.error_bars_upper[idx] else 0)
            return vals, elo, ehi

        # ================================================================
        # 1. BOX
        # ================================================================
        if pt == PlotType.BOX:
            # Try to parse structured box stats from notes
            parsed_boxes = []
            outlier_map: dict[str, list[float]] = {}
            if result.notes:
                # Parse "GroupName: min=X, Q1=X, median=X, Q3=X, max=X"
                box_pat = re.compile(
                    r"(.+?):\s*min\s*=\s*([\d.\-eE+]+)\s*,\s*Q1\s*=\s*([\d.\-eE+]+)\s*,"
                    r"\s*median\s*=\s*([\d.\-eE+]+)\s*,\s*Q3\s*=\s*([\d.\-eE+]+)\s*,"
                    r"\s*max\s*=\s*([\d.\-eE+]+)",
                    re.IGNORECASE,
                )
                for m in box_pat.finditer(result.notes):
                    parsed_boxes.append({
                        "label": m.group(1).strip(),
                        "min": float(m.group(2)),
                        "q1": float(m.group(3)),
                        "med": float(m.group(4)),
                        "q3": float(m.group(5)),
                        "max": float(m.group(6)),
                    })
                # Parse "GroupName outliers: [val1, val2]"
                out_pat = re.compile(r"(.+?)\s+outliers:\s*\[([\d.,\s\-eE+]*)\]", re.IGNORECASE)
                for m in out_pat.finditer(result.notes):
                    name = m.group(1).strip()
                    vals_str = m.group(2).strip()
                    if vals_str:
                        outlier_map[name] = [float(v.strip()) for v in vals_str.split(",") if v.strip()]

            if parsed_boxes:
                positions = list(range(1, len(parsed_boxes) + 1))
                for i, b in enumerate(parsed_boxes):
                    bp = ax.bxp(
                        [{"whislo": b["min"], "q1": b["q1"], "med": b["med"],
                          "q3": b["q3"], "whishi": b["max"],
                          "fliers": outlier_map.get(b["label"], [])}],
                        positions=[positions[i]],
                        widths=0.5,
                        showfliers=True,
                        patch_artist=True,
                    )
                    for patch in bp.get("boxes", []):
                        patch.set_facecolor(colors[i % len(colors)])
                        patch.set_alpha(0.7)
                ax.set_xticks(positions)
                ax.set_xticklabels([b["label"] for b in parsed_boxes],
                                   rotation=35 if any(len(b["label"]) > 5 for b in parsed_boxes) else 0,
                                   ha="right")
            else:
                # Fallback: simple boxplot from y_values
                box_data = []
                labels = []
                for s in result.series:
                    labels.append(s.name)
                    box_data.append(s.y_values if s.y_values else [0])
                if box_data:
                    bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
                    for i, patch in enumerate(bp["boxes"]):
                        patch.set_facecolor(colors[i % len(colors)])
                        patch.set_alpha(0.7)

        # ================================================================
        # 2. BAR (with jittered-dot detection)
        # ================================================================
        elif pt == PlotType.BAR:
            unique_cats = _unique_cats()
            has_repeated = (
                unique_cats
                and any(len(s.y_values) > len(unique_cats) * 1.5 for s in result.series)
            )

            if has_repeated:
                # Jittered dot plot (individual observations)
                for si, s in enumerate(result.series):
                    col = colors[si % len(colors)]
                    by_cat: dict[str, list[float]] = {}
                    for i, y in enumerate(s.y_values):
                        cat = str(s.x_values[i]) if i < len(s.x_values) else ""
                        by_cat.setdefault(cat, []).append(y)
                    for ci, cat in enumerate(unique_cats):
                        pts = by_cat.get(cat, [])
                        if not pts:
                            continue
                        jitter = np.random.normal(0, 0.08, len(pts))
                        x_pos = ci + (si - len(result.series) / 2 + 0.5) * 0.25
                        ax.scatter([x_pos + j for j in jitter], pts, color=col, s=20,
                                   alpha=0.7, label=s.name if ci == 0 else None, zorder=3)
                ax.set_xticks(range(len(unique_cats)))
                ax.set_xticklabels(unique_cats,
                                   rotation=35 if unique_cats and max(len(c) for c in unique_cats) > 5 else 0,
                                   ha="right")
            else:
                # Standard grouped bar chart with error bars
                if unique_cats:
                    x_pos = np.arange(len(unique_cats))
                    width = 0.8 / max(1, len(result.series))
                    for si, s in enumerate(result.series):
                        col = colors[si % len(colors)]
                        offset = (si - len(result.series) / 2 + 0.5) * width
                        vals, elo, ehi = _err_arrays(s, unique_cats)
                        has_err = any(e for e in elo + ehi)
                        ax.bar(x_pos + offset, vals, width, color=col, alpha=0.8,
                               label=s.name,
                               yerr=[elo, ehi] if has_err else None, capsize=2)
                    ax.set_xticks(x_pos)
                    ax.set_xticklabels(unique_cats,
                                       rotation=35 if max(len(c) for c in unique_cats) > 5 else 0,
                                       ha="right")
                else:
                    # Numeric x-axis bar
                    for si, s in enumerate(result.series):
                        col = colors[si % len(colors)]
                        xv, yv = _xy_numeric(s)
                        if xv:
                            ax.bar(xv, yv, color=col, alpha=0.8, label=s.name)

        # ================================================================
        # 3. STACKED_BAR
        # ================================================================
        elif pt == PlotType.STACKED_BAR:
            unique_cats = _unique_cats()
            if unique_cats:
                x_pos = np.arange(len(unique_cats))
                bottoms = np.zeros(len(unique_cats))
                width = 0.6
                for si, s in enumerate(result.series):
                    col = colors[si % len(colors)]
                    vals, _, _ = _err_arrays(s, unique_cats)
                    vals_arr = np.array(vals, dtype=float)
                    ax.bar(x_pos, vals_arr, width, bottom=bottoms, color=col,
                           alpha=0.85, label=s.name)
                    bottoms += vals_arr
                ax.set_xticks(x_pos)
                ax.set_xticklabels(unique_cats,
                                   rotation=35 if max(len(c) for c in unique_cats) > 5 else 0,
                                   ha="right")

        # ================================================================
        # 4. HISTOGRAM
        # ================================================================
        elif pt == PlotType.HISTOGRAM:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    # x = bin centers, y = counts/frequencies
                    if len(xv) >= 2:
                        bar_width = abs(xv[1] - xv[0]) * 0.9
                    else:
                        bar_width = 1.0
                    ax.bar(xv, yv, width=bar_width, color=col, alpha=0.7,
                           edgecolor="white", label=s.name)
                elif yv:
                    # Only y-values: treat as raw data for histogram
                    ax.hist(yv, bins="auto", color=col, alpha=0.7, label=s.name)

        # ================================================================
        # 5. VIOLIN
        # ================================================================
        elif pt == PlotType.VIOLIN:
            # Attempt proper violin; fallback to box-like
            data_list = []
            labels = []
            for s in result.series:
                labels.append(s.name)
                data_list.append(s.y_values if s.y_values else [0])
            if data_list:
                # Violin needs at least 2 points per group for KDE
                can_violin = all(len(d) >= 2 for d in data_list)
                if can_violin:
                    parts = ax.violinplot(data_list, showmeans=True, showmedians=True)
                    for i, pc in enumerate(parts.get("bodies", [])):
                        pc.set_facecolor(colors[i % len(colors)])
                        pc.set_alpha(0.7)
                    ax.set_xticks(range(1, len(labels) + 1))
                    ax.set_xticklabels(labels,
                                       rotation=35 if any(len(l) > 5 for l in labels) else 0,
                                       ha="right")
                else:
                    bp = ax.boxplot(data_list, labels=labels, patch_artist=True)
                    for i, patch in enumerate(bp["boxes"]):
                        patch.set_facecolor(colors[i % len(colors)])
                        patch.set_alpha(0.7)

        # ================================================================
        # 6. FOREST
        # ================================================================
        elif pt == PlotType.FOREST:
            y_positions = []
            labels_list = []
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if not xv and s.x_values:
                    # x_values are study names, y_values are effect sizes
                    names = [str(v) for v in s.x_values]
                    n = min(len(names), len(s.y_values))
                    positions = list(range(n))
                    effects = s.y_values[:n]
                    elo = [abs(s.error_bars_lower[i]) if i < len(s.error_bars_lower) and s.error_bars_lower[i] is not None else 0 for i in range(n)]
                    ehi = [abs(s.error_bars_upper[i]) if i < len(s.error_bars_upper) and s.error_bars_upper[i] is not None else 0 for i in range(n)]
                    ax.errorbar(effects, positions, xerr=[elo, ehi],
                                fmt="o" if si == 0 else "s", color=col, capsize=3,
                                markersize=6, label=s.name, zorder=3)
                    if not labels_list:
                        labels_list = names[:n]
                        y_positions = positions
                elif xv and yv:
                    # x = effect size, y used for positioning
                    n = len(xv)
                    positions = list(range(len(y_positions), len(y_positions) + n))
                    elo = [abs(s.error_bars_lower[i]) if i < len(s.error_bars_lower) and s.error_bars_lower[i] is not None else 0 for i in range(n)]
                    ehi = [abs(s.error_bars_upper[i]) if i < len(s.error_bars_upper) and s.error_bars_upper[i] is not None else 0 for i in range(n)]
                    ax.errorbar(xv, positions, xerr=[elo, ehi],
                                fmt="o", color=col, capsize=3, markersize=6,
                                label=s.name, zorder=3)
                    y_positions = positions
            # Vertical reference line
            ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, zorder=1)
            # If mostly positive values (e.g., OR/RR), also add line at 1
            all_y = [v for s in result.series for v in s.y_values if v is not None]
            if all_y and min(all_y) > 0:
                ax.axvline(x=1, color="gray", linestyle="--", linewidth=0.8, zorder=1)
            if labels_list:
                ax.set_yticks(y_positions)
                ax.set_yticklabels(labels_list, fontsize=7)
            ax.invert_yaxis()

        # ================================================================
        # 7. KAPLAN_MEIER
        # ================================================================
        elif pt == PlotType.KAPLAN_MEIER:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.step(xv, yv, where="post", color=col, linewidth=1.5,
                            label=s.name)
            ax.set_ylim(-0.05, 1.05)

        # ================================================================
        # 8. SCATTER
        # ================================================================
        elif pt == PlotType.SCATTER:
            markers = ["o", "s", "^", "D", "v", "P", "*", "X", "p", "h"]
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                mk = markers[si % len(markers)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.scatter(xv, yv, color=col, s=20, alpha=0.8,
                               marker=mk, label=s.name)

        # ================================================================
        # 9. LINE
        # ================================================================
        elif pt == PlotType.LINE:
            markers = ["o", "s", "^", "D", "v", "P", "*", "X"]
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                mk = markers[si % len(markers)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    # Add error bands if available
                    elo = [s.error_bars_lower[i] if i < len(s.error_bars_lower) and s.error_bars_lower[i] is not None else 0 for i in range(len(yv))]
                    ehi = [s.error_bars_upper[i] if i < len(s.error_bars_upper) and s.error_bars_upper[i] is not None else 0 for i in range(len(yv))]
                    ax.plot(xv, yv, color=col, marker=mk, markersize=4,
                            linewidth=1.5, label=s.name)
                    if any(e for e in elo + ehi):
                        lower = [y - e for y, e in zip(yv, elo)]
                        upper = [y + e for y, e in zip(yv, ehi)]
                        ax.fill_between(xv, lower, upper, color=col, alpha=0.15)

        # ================================================================
        # 10. DOT_STRIP
        # ================================================================
        elif pt == PlotType.DOT_STRIP:
            unique_cats = _unique_cats()
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                by_cat: dict[str, list[float]] = {}
                for i, y in enumerate(s.y_values):
                    cat = str(s.x_values[i]) if i < len(s.x_values) else s.name
                    by_cat.setdefault(cat, []).append(y)
                for ci, cat in enumerate(unique_cats):
                    pts = by_cat.get(cat, [])
                    if not pts:
                        continue
                    jitter = np.random.normal(0, 0.1, len(pts))
                    x_pos = ci + (si - len(result.series) / 2 + 0.5) * 0.3
                    ax.scatter([x_pos + j for j in jitter], pts, color=col, s=18,
                               alpha=0.7, label=s.name if ci == 0 else None, zorder=3)
            if unique_cats:
                ax.set_xticks(range(len(unique_cats)))
                ax.set_xticklabels(unique_cats,
                                   rotation=35 if max(len(c) for c in unique_cats) > 5 else 0,
                                   ha="right")

        # ================================================================
        # 11. PAIRED
        # ================================================================
        elif pt == PlotType.PAIRED:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if len(xv) >= 2 and len(yv) >= 2:
                    # Assume pairs: (x[0],y[0]) connected to (x[1],y[1]), etc.
                    # Or two timepoints per subject
                    ax.plot(xv, yv, "o-", color=col, alpha=0.5, markersize=5,
                            label=s.name)
                elif yv and len(yv) % 2 == 0:
                    # Before-after pairs with x=0,1
                    half = len(yv) // 2
                    for i in range(half):
                        ax.plot([0, 1], [yv[i], yv[half + i]], "o-", color=col,
                                alpha=0.4, markersize=5)
                    ax.plot([], [], "o-", color=col, label=s.name)  # legend entry
                    ax.set_xticks([0, 1])
                    ax.set_xticklabels(["Before", "After"])

        # ================================================================
        # 12. ERROR_BAR
        # ================================================================
        elif pt == PlotType.ERROR_BAR:
            unique_cats = _unique_cats()
            is_cat = unique_cats and any(isinstance(s.x_values[0], str) for s in result.series if s.x_values)
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                if is_cat and unique_cats:
                    x_pos = np.arange(len(unique_cats))
                    offset = (si - len(result.series) / 2 + 0.5) * 0.2
                    vals, elo, ehi = _err_arrays(s, unique_cats)
                    has_err = any(e for e in elo + ehi)
                    ax.errorbar(x_pos + offset, vals,
                                yerr=[elo, ehi] if has_err else None,
                                fmt="o", color=col, capsize=4, markersize=6,
                                label=s.name)
                    ax.set_xticks(x_pos)
                    ax.set_xticklabels(unique_cats,
                                       rotation=35 if max(len(c) for c in unique_cats) > 5 else 0,
                                       ha="right")
                else:
                    xv, yv = _xy_numeric(s)
                    n = len(yv)
                    elo = [abs(s.error_bars_lower[i]) if i < len(s.error_bars_lower) and s.error_bars_lower[i] is not None else 0 for i in range(n)]
                    ehi = [abs(s.error_bars_upper[i]) if i < len(s.error_bars_upper) and s.error_bars_upper[i] is not None else 0 for i in range(n)]
                    has_err = any(e for e in elo + ehi)
                    ax.errorbar(xv, yv,
                                yerr=[elo, ehi] if has_err else None,
                                fmt="o", color=col, capsize=4, markersize=6,
                                label=s.name)

        # ================================================================
        # 13. WATERFALL
        # ================================================================
        elif pt == PlotType.WATERFALL:
            # Typically one series: ordered bars going up/down from baseline
            for si, s in enumerate(result.series):
                yv = [v for v in s.y_values if v is not None]
                if not yv:
                    continue
                x_indices = np.arange(len(yv))
                bar_colors = ["#ef4444" if v >= 0 else "#10b981" for v in yv]
                ax.bar(x_indices, yv, color=bar_colors, edgecolor="white",
                       linewidth=0.3, label=s.name)
                ax.axhline(y=0, color="black", linewidth=0.5)
            # Reference lines common in waterfall (e.g., +20%, -30%)
            if result.notes and "20" in result.notes:
                ax.axhline(y=20, color="gray", linestyle="--", linewidth=0.5)
                ax.axhline(y=-30, color="gray", linestyle="--", linewidth=0.5)

        # ================================================================
        # 14. ROC
        # ================================================================
        elif pt == PlotType.ROC:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.plot(xv, yv, color=col, linewidth=1.5, label=s.name)
            ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5, label="Reference")
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
            ax.set_aspect("equal", adjustable="box")

        # ================================================================
        # 15. HEATMAP
        # ================================================================
        elif pt == PlotType.HEATMAP:
            # Build matrix from series: each series = one row
            if result.series:
                matrix = []
                row_labels = []
                col_labels = []
                for s in result.series:
                    row_labels.append(s.name)
                    matrix.append(s.y_values if s.y_values else [])
                    if not col_labels and s.x_values:
                        col_labels = [str(v) for v in s.x_values]
                if matrix:
                    # Pad rows to same length
                    max_len = max(len(r) for r in matrix)
                    padded = [r + [0.0] * (max_len - len(r)) for r in matrix]
                    data = np.array(padded, dtype=float)
                    im = ax.imshow(data, cmap="RdBu_r", aspect="auto")
                    fig.colorbar(im, ax=ax, shrink=0.8)
                    if col_labels:
                        ax.set_xticks(range(min(len(col_labels), max_len)))
                        ax.set_xticklabels(col_labels[:max_len], rotation=45, ha="right", fontsize=7)
                    ax.set_yticks(range(len(row_labels)))
                    ax.set_yticklabels(row_labels, fontsize=7)

        # ================================================================
        # 16. VOLCANO
        # ================================================================
        elif pt == PlotType.VOLCANO:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.scatter(xv, yv, color=col, s=12, alpha=0.6, label=s.name)
            # Threshold lines
            ax.axhline(y=1.3, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)  # -log10(0.05)
            ax.axvline(x=-1, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)
            ax.axvline(x=1, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)

        # ================================================================
        # 17. FUNNEL
        # ================================================================
        elif pt == PlotType.FUNNEL:
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.scatter(xv, yv, color=col, s=25, alpha=0.7, label=s.name)
            # Reference line at pooled estimate (mean of x)
            all_x = [v for s in result.series for v in s.x_values if v is not None and not isinstance(v, str)]
            if all_x:
                mean_x = float(np.mean(all_x))
                ax.axvline(x=mean_x, color="black", linestyle="-", linewidth=0.8)
                # Approximate CI bounds (widening with SE)
                all_y = [v for s in result.series for v in s.y_values if v is not None]
                if all_y:
                    y_range = np.linspace(0, max(all_y), 50)
                    ax.plot(mean_x - 1.96 * y_range, y_range, "k--", linewidth=0.5, alpha=0.4)
                    ax.plot(mean_x + 1.96 * y_range, y_range, "k--", linewidth=0.5, alpha=0.4)
            ax.invert_yaxis()

        # ================================================================
        # 18. DOSE_RESPONSE
        # ================================================================
        elif pt == PlotType.DOSE_RESPONSE:
            markers = ["o", "s", "^", "D", "v"]
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                mk = markers[si % len(markers)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    n = len(yv)
                    elo = [abs(s.error_bars_lower[i]) if i < len(s.error_bars_lower) and s.error_bars_lower[i] is not None else 0 for i in range(n)]
                    ehi = [abs(s.error_bars_upper[i]) if i < len(s.error_bars_upper) and s.error_bars_upper[i] is not None else 0 for i in range(n)]
                    has_err = any(e for e in elo + ehi)
                    ax.errorbar(xv, yv,
                                yerr=[elo, ehi] if has_err else None,
                                fmt=f"-{mk}", color=col, capsize=3, markersize=5,
                                linewidth=1.5, label=s.name)
            ax.set_xscale("log")

        # ================================================================
        # 19. AREA
        # ================================================================
        elif pt == PlotType.AREA:
            # Stacked area if multiple series
            if len(result.series) > 1:
                all_x: list[float] = []
                all_y: list[list[float]] = []
                area_labels: list[str] = []
                for s in result.series:
                    xv, yv = _xy_numeric(s)
                    if xv and yv:
                        if not all_x:
                            all_x = xv
                        all_y.append(yv[:len(all_x)])
                        area_labels.append(s.name)
                if all_x and all_y:
                    # Pad to same length
                    min_len = min(len(all_x), *(len(y) for y in all_y))
                    ax.stackplot(all_x[:min_len], *[y[:min_len] for y in all_y],
                                 labels=area_labels,
                                 colors=colors[:len(all_y)], alpha=0.7)
            else:
                for si, s in enumerate(result.series):
                    col = colors[si % len(colors)]
                    xv, yv = _xy_numeric(s)
                    if xv and yv:
                        ax.fill_between(xv, yv, alpha=0.5, color=col, label=s.name)
                        ax.plot(xv, yv, color=col, linewidth=1)

        # ================================================================
        # 20. BLAND_ALTMAN
        # ================================================================
        elif pt == PlotType.BLAND_ALTMAN:
            all_x_vals: list[float] = []
            all_y_vals: list[float] = []
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.scatter(xv, yv, color=col, s=20, alpha=0.7, label=s.name)
                    all_x_vals.extend(xv)
                    all_y_vals.extend(yv)
            if all_y_vals:
                mean_diff = float(np.mean(all_y_vals))
                sd_diff = float(np.std(all_y_vals))
                ax.axhline(y=mean_diff, color="blue", linestyle="-", linewidth=1,
                            alpha=0.7, label=f"Bias ({mean_diff:.2f})")
                ax.axhline(y=mean_diff + 1.96 * sd_diff, color="red", linestyle="--",
                            linewidth=0.8, alpha=0.7, label=f"+1.96 SD")
                ax.axhline(y=mean_diff - 1.96 * sd_diff, color="red", linestyle="--",
                            linewidth=0.8, alpha=0.7, label=f"-1.96 SD")

        # ================================================================
        # FALLBACK: scatter
        # ================================================================
        else:
            markers = ["o", "s", "^", "D", "v", "P", "*", "X"]
            for si, s in enumerate(result.series):
                col = colors[si % len(colors)]
                mk = markers[si % len(markers)]
                xv, yv = _xy_numeric(s)
                if xv and yv:
                    ax.scatter(xv, yv, color=col, s=20, alpha=0.8,
                               marker=mk, label=s.name)

        # ── Axis scales ──
        if result.x_scale == "log" and pt != PlotType.DOSE_RESPONSE:
            try:
                ax.set_xscale("log")
            except Exception:
                pass
        if result.y_scale == "log":
            try:
                ax.set_yscale("log")
            except Exception:
                pass

        # ── Labels and title ──
        if result.x_label:
            xlabel = result.x_label + (f" ({result.x_unit})" if result.x_unit else "")
            ax.set_xlabel(xlabel, fontsize=9)
        if result.y_label:
            ylabel = result.y_label + (f" ({result.y_unit})" if result.y_unit else "")
            ax.set_ylabel(ylabel, fontsize=9)
        if result.title:
            ax.set_title(result.title, fontsize=10, pad=8)
        if len(result.series) > 1:
            ax.legend(fontsize=7, loc="best")

        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120)
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.getvalue()).decode()

    except Exception:
        logger.exception("Failed to render result to PNG")
        try:
            plt.close(fig)
        except Exception:
            pass
        return ""


async def digitize_figure(figure: Figure) -> tuple[list[ExtractedData], list[Figure]]:
    """Digitize a figure using a multi-phase pipeline:

    Phase 1 — Panel detection: identify sub-panels in the image
    Phase 2 — Pre-analysis: identify axes, scale, plot type, legend (per panel)
    Phase 3 — Data extraction: extract numerical data with full context
    Phase 4 — Heuristic validation: cheap sanity checks, no API calls

    Returns (results, panel_figures):
      - results: list of ExtractedData, one per panel
      - panel_figures: list of cropped Figure objects (empty for single-panel)
    """
    # ── Phase 1: Detect panels ──
    try:
        panels = await detect_panels(figure)
    except Exception as e:
        logger.warning(f"Panel detection failed for {figure.figure_id}: {e}, treating as single panel")
        panels = [{"label": "main", "plot_type": figure.plot_type.value, "bbox": [0, 0, figure.width, figure.height]}]

    # Single panel — run phases 2+3+4 directly
    if len(panels) <= 1:
        # Use panel detection's plot_type if it identified one
        if panels and panels[0].get("plot_type") and panels[0]["plot_type"] != "other":
            try:
                figure = figure.model_copy(update={"plot_type": PlotType(panels[0]["plot_type"])})
            except ValueError:
                pass
        result, pre_analysis = await _digitize_single(figure)
        result = validate_extraction(result, pre_analysis)
        return [result], []

    # ── Multi-panel — crop, pre-analyze, and extract each ──
    logger.info(f"Splitting {figure.figure_id} into {len(panels)} panels: {[p['label'] for p in panels]}")
    panel_figures = []
    panel_labels = []
    for p in panels:
        pf = _crop_panel(figure, p["bbox"], p["label"], p.get("plot_type", "other"))
        panel_figures.append(pf)
        panel_labels.append(p["label"])

    # Pass the full figure image as context so panels can reference
    # axes, legends, and labels from the complete figure
    full_b64 = figure.image_base64

    semaphore = asyncio.Semaphore(3)

    async def _extract_panel(pf: Figure, label: str) -> ExtractedData:
        async with semaphore:
            try:
                result, pre_analysis = await _digitize_single(pf, panel_label=label, full_figure_b64=full_b64)
                result = validate_extraction(result, pre_analysis)
                return result
            except Exception as e:
                logger.error(f"Panel extraction failed for {pf.figure_id}: {e}")
                return ExtractedData(
                    figure_id=pf.figure_id,
                    plot_type=pf.plot_type,
                    confidence=Confidence.LOW,
                    notes=f"Extraction failed: {str(e)}",
                )

    results = await asyncio.gather(*[_extract_panel(pf, lbl) for pf, lbl in zip(panel_figures, panel_labels)])
    return list(results), panel_figures


async def digitize_figures(
    figures: list[Figure], max_concurrent: int = 3
) -> tuple[list[ExtractedData], list[Figure]]:
    """Digitize multiple figures with bounded concurrency.

    Returns (results, panel_figures) — panel_figures contains cropped
    sub-panel Figure objects for multi-panel images so they can be served.
    """
    semaphore = asyncio.Semaphore(max_concurrent)
    all_panel_figures: list[Figure] = []

    async def _extract(fig: Figure) -> tuple[list[ExtractedData], list[Figure]]:
        async with semaphore:
            try:
                return await digitize_figure(fig)
            except Exception as e:
                logger.error(f"Extraction failed for {fig.figure_id}: {e}")
                return [ExtractedData(
                    figure_id=fig.figure_id,
                    plot_type=fig.plot_type,
                    confidence=Confidence.LOW,
                    notes=f"Extraction failed: {str(e)}",
                )], []

    nested = await asyncio.gather(*[_extract(fig) for fig in figures])
    results = []
    for data_list, panels in nested:
        results.extend(data_list)
        all_panel_figures.extend(panels)
    return results, all_panel_figures


# ── CSV / JSON export ──


def export_csv(results: list[ExtractedData], output_path: str) -> str:
    """Export extraction results to a CSV file. Returns the file path."""
    import csv
    from pathlib import Path

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "figure_id", "plot_type", "title", "series_name",
            "x_label", "y_label", "x_unit", "y_unit",
            "x_value", "y_value", "error_lower", "error_upper",
            "confidence", "notes",
        ])
        for r in results:
            for s in r.series:
                for i in range(len(s.y_values)):
                    x = s.x_values[i] if i < len(s.x_values) else ""
                    y = s.y_values[i]
                    el = s.error_bars_lower[i] if i < len(s.error_bars_lower) else None
                    eu = s.error_bars_upper[i] if i < len(s.error_bars_upper) else None
                    writer.writerow([
                        r.figure_id, r.plot_type.value, r.title or "", s.name,
                        r.x_label or "", r.y_label or "", r.x_unit or "", r.y_unit or "",
                        x, y,
                        el if el is not None else "",
                        eu if eu is not None else "",
                        r.confidence.value, r.notes or "",
                    ])
    return str(path)


def export_json(results: list[ExtractedData], output_path: str) -> str:
    """Export extraction results to a JSON file. Returns the file path."""
    from pathlib import Path

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [r.model_dump() for r in results]
    path.write_text(json.dumps(data, indent=2, default=str))
    return str(path)


async def extract_from_image(
    image_path: str | None = None,
    image_bytes: bytes | None = None,
    plot_type: str | None = None,
) -> tuple[list[ExtractedData], list[Figure]]:
    """Extract data from an image file or bytes.

    Convenience wrapper that handles loading, base64-encoding, and
    creating the Figure object before calling digitize_figure().
    """
    if image_bytes is None:
        if image_path is None:
            raise ValueError("Either image_path or image_bytes must be provided")
        from pathlib import Path
        image_bytes = Path(image_path).read_bytes()

    img = Image.open(io.BytesIO(image_bytes))
    if img.mode == "RGBA":
        bg = Image.new("RGB", img.size, (255, 255, 255))
        bg.paste(img, mask=img.split()[3])
        img = bg
    elif img.mode != "RGB":
        img = img.convert("RGB")

    # Resize if too large
    if img.width > 1500 or img.height > 1500:
        img.thumbnail((1500, 1500), Image.LANCZOS)

    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    b64 = base64.b64encode(buf.getvalue()).decode()

    pt = PlotType.OTHER
    if plot_type:
        try:
            pt = PlotType(plot_type)
        except ValueError:
            pass

    figure = Figure(
        figure_id="image",
        paper_id="standalone",
        page_number=1,
        image_index=0,
        width=img.width,
        height=img.height,
        image_base64=b64,
        plot_type=pt,
        plot_type_confidence=Confidence.MEDIUM,
    )

    return await digitize_figure(figure)
