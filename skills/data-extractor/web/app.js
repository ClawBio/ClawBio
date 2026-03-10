// Get Me The Data — Frontend (screenshot-only mode)

var currentImageId = null;
var extractedResults = [];

// ── Helpers ──

function show(id) { document.getElementById(id).classList.remove("hidden"); }
function hide(id) { document.getElementById(id).classList.add("hidden"); }

function setStatus(id, msg, type) {
    type = type || "info";
    var el = document.getElementById(id);
    if (!el) return;
    var colors = { info: "text-gray-500", success: "text-green-600", error: "text-red-600", loading: "text-blue-600" };
    el.className = "text-sm " + (colors[type] || colors.info);
    el.innerHTML = type === "loading" ? '<span class="spinner"></span> ' + msg : msg;
}

async function api(url, options) {
    options = options || {};
    var resp = await fetch(url, {
        headers: { "Content-Type": "application/json" },
        ...options,
    });
    if (!resp.ok) {
        var err = await resp.json().catch(function() { return { detail: resp.statusText }; });
        throw new Error(err.detail || "Request failed");
    }
    return resp;
}

// ── Image Upload ──

async function uploadImage(file) {
    setStatus("viewer-status", "Uploading image...", "loading");
    var formData = new FormData();
    formData.append("file", file);
    try {
        var resp = await fetch("/api/upload-image", { method: "POST", body: formData });
        if (!resp.ok) {
            var err = await resp.json();
            throw new Error(err.detail);
        }
        var data = await resp.json();
        currentImageId = data.image_id;
        document.getElementById("main-image").src = "/api/image/" + data.image_id;
        hide("drop-section");
        show("viewer-section");
        setStatus("viewer-status", "Image loaded. Draw boxes around plots to extract data.", "success");

        document.getElementById("main-image").onload = function() {
            _initDrawing();
        };
    } catch (e) {
        setStatus("viewer-status", "Upload failed: " + e.message, "error");
    }
}

function handleImageUpload(event) {
    var file = event.target.files[0];
    if (file) uploadImage(file);
}

async function loadExample(name) {
    try {
        var resp = await fetch("/static/examples/" + name + ".png");
        if (!resp.ok) throw new Error("Example not found");
        var blob = await resp.blob();
        var file = new File([blob], name + ".png", { type: "image/png" });
        uploadImage(file);
    } catch (e) {
        setStatus("viewer-status", "Failed to load example: " + e.message, "error");
    }
}

function clearImage() {
    currentImageId = null;
    extractedResults = [];
    pageBoxes = {};
    boxIdCounter = 0;
    _drawingInitialized = false;
    document.getElementById("main-image").src = "";
    document.getElementById("plot-boxes-layer").innerHTML = "";
    document.getElementById("page-results").innerHTML = "";
    document.getElementById("results-container").innerHTML = "";
    hide("viewer-section");
    hide("results-section");
    show("drop-section");
    document.getElementById("image-upload").value = "";
}

// ── Drag-and-drop + paste ──

// Prevent browser default file-open on ALL drag/drop events globally.
// Without this, dropping a file anywhere opens it in a new tab.
window.addEventListener("dragover", function(e) { e.preventDefault(); }, false);
window.addEventListener("drop", function(e) { e.preventDefault(); }, false);

document.addEventListener("DOMContentLoaded", function() {
    var zone = document.getElementById("drop-zone");
    zone.addEventListener("dragover", function(e) {
        e.preventDefault();
        e.stopPropagation();
        zone.classList.add("dragover");
    });
    zone.addEventListener("dragleave", function() { zone.classList.remove("dragover"); });
    zone.addEventListener("drop", function(e) {
        e.preventDefault();
        e.stopPropagation();
        zone.classList.remove("dragover");
        var file = e.dataTransfer.files[0];
        if (file && file.type.startsWith("image/")) {
            uploadImage(file);
        }
    });

    // Also allow dropping on the whole page (handles drops outside the zone)
    document.addEventListener("drop", function(e) {
        e.preventDefault();
        var file = e.dataTransfer && e.dataTransfer.files[0];
        if (file && file.type.startsWith("image/")) {
            uploadImage(file);
        }
    });

    // Paste from clipboard
    document.addEventListener("paste", function(e) {
        var items = e.clipboardData && e.clipboardData.items;
        if (!items) return;
        for (var i = 0; i < items.length; i++) {
            if (items[i].type.indexOf("image") !== -1) {
                var file = items[i].getAsFile();
                if (file) uploadImage(file);
                break;
            }
        }
    });
});

// ── Drawing rectangles ──

var pageBoxes = {};  // { 1: [ {id, num, x_pct, y_pct, w_pct, h_pct} ] }
var boxIdCounter = 0;
var _drawingInitialized = false;

function _initDrawing() {
    if (_drawingInitialized) return;
    _drawingInitialized = true;

    var drawLayer = document.getElementById("draw-layer");
    var drawing = false;
    var preview = null;
    var startX, startY;

    drawLayer.addEventListener("mousedown", function(e) {
        if (e.target !== drawLayer) return;
        e.preventDefault();
        drawing = true;
        var rect = drawLayer.getBoundingClientRect();
        startX = e.clientX - rect.left;
        startY = e.clientY - rect.top;

        preview = document.createElement("div");
        preview.className = "draw-preview";
        preview.style.left = startX + "px";
        preview.style.top = startY + "px";
        preview.style.width = "0px";
        preview.style.height = "0px";
        drawLayer.appendChild(preview);
    });

    document.addEventListener("mousemove", function(e) {
        if (!drawing || !preview) return;
        var rect = drawLayer.getBoundingClientRect();
        var curX = Math.max(0, Math.min(e.clientX - rect.left, rect.width));
        var curY = Math.max(0, Math.min(e.clientY - rect.top, rect.height));

        var x = Math.min(startX, curX), y = Math.min(startY, curY);
        var w = Math.abs(curX - startX), h = Math.abs(curY - startY);

        preview.style.left = x + "px";
        preview.style.top = y + "px";
        preview.style.width = w + "px";
        preview.style.height = h + "px";
    });

    document.addEventListener("mouseup", function(e) {
        if (!drawing || !preview) return;
        drawing = false;

        var rect = drawLayer.getBoundingClientRect();
        var px = parseFloat(preview.style.left);
        var py = parseFloat(preview.style.top);
        var pw = parseFloat(preview.style.width);
        var ph = parseFloat(preview.style.height);

        preview.remove();
        preview = null;

        if (pw < 20 || ph < 20) return;

        var x_pct = px / rect.width * 100;
        var y_pct = py / rect.height * 100;
        var w_pct = pw / rect.width * 100;
        var h_pct = ph / rect.height * 100;

        _addBox(x_pct, y_pct, w_pct, h_pct);
    });
}

function _addBox(x_pct, y_pct, w_pct, h_pct) {
    if (!pageBoxes[1]) pageBoxes[1] = [];
    boxIdCounter++;
    var num = pageBoxes[1].length + 1;
    var b = {
        id: boxIdCounter,
        num: num,
        x_pct: x_pct, y_pct: y_pct,
        w_pct: w_pct, h_pct: h_pct,
    };
    pageBoxes[1].push(b);
    _renderOneBox(b);
}

function _renumberBoxes() {
    var boxes = pageBoxes[1] || [];
    boxes.forEach(function(b, i) {
        b.num = i + 1;
        var label = document.querySelector("#plot-box-" + b.id + " .box-label");
        if (label) label.textContent = "#" + b.num;
    });
}

async function autoDetectPlots() {
    if (!currentImageId) return;
    setStatus("viewer-status", "Auto-detecting plots...", "loading");
    document.getElementById("auto-detect-btn").disabled = true;
    try {
        var resp = await api("/api/detect-plots/" + currentImageId);
        var data = await resp.json();
        var plots = data.plots || [];
        plots.forEach(function(p) {
            _addBox(p.x_pct, p.y_pct, p.w_pct, p.h_pct);
        });
        setStatus("viewer-status", plots.length + " plot(s) detected and added", "success");
    } catch (e) {
        setStatus("viewer-status", "Detection failed: " + e.message, "error");
    } finally {
        document.getElementById("auto-detect-btn").disabled = false;
    }
}

// ── Plot box rendering ──

function renderPlotBoxes() {
    var layer = document.getElementById("plot-boxes-layer");
    layer.innerHTML = "";
    var boxes = pageBoxes[1] || [];
    boxes.forEach(function(b) { _renderOneBox(b); });
}

function _renderOneBox(b) {
    var layer = document.getElementById("plot-boxes-layer");
    var el = document.createElement("div");
    el.className = "plot-box";
    el.id = "plot-box-" + b.id;
    el.style.cssText = "left:" + b.x_pct + "%;top:" + b.y_pct + "%;width:" + b.w_pct + "%;height:" + b.h_pct + "%";

    el.innerHTML =
        '<span class="box-label">#' + b.num + '</span>' +
        '<div class="box-controls">' +
            '<button class="box-btn btn-extract" onclick="extractBox(' + b.id + ')">Extract</button>' +
            '<button class="box-btn btn-delete" onclick="deleteBox(' + b.id + ')">&times;</button>' +
        '</div>' +
        '<div class="box-handle se" data-handle="se"></div>' +
        '<div class="box-handle sw" data-handle="sw"></div>' +
        '<div class="box-handle ne" data-handle="ne"></div>' +
        '<div class="box-handle nw" data-handle="nw"></div>';

    layer.appendChild(el);
    _makeBoxDraggable(b.id, el);
}

function _makeBoxDraggable(boxId, el) {
    function clamp(v, mn, mx) { return Math.max(mn, Math.min(mx, v)); }

    el.addEventListener("mousedown", function(e) {
        if (e.target.classList.contains("box-btn")) return;
        e.preventDefault();
        e.stopPropagation();
        var handle = e.target.dataset.handle;
        var dragging = handle || "move";
        var layer = document.getElementById("plot-boxes-layer");
        var lw = layer.clientWidth, lh = layer.clientHeight;
        var startX = e.clientX, startY = e.clientY;
        var startL = parseFloat(el.style.left) / 100 * lw;
        var startT = parseFloat(el.style.top) / 100 * lh;
        var startW = parseFloat(el.style.width) / 100 * lw;
        var startH = parseFloat(el.style.height) / 100 * lh;

        function onMove(e) {
            var dx = e.clientX - startX, dy = e.clientY - startY;
            var nl, nt, nw, nh;

            if (dragging === "move") {
                nl = clamp(startL + dx, 0, lw - startW); nt = clamp(startT + dy, 0, lh - startH);
                nw = startW; nh = startH;
            } else if (dragging === "se") {
                nl = startL; nt = startT;
                nw = clamp(startW + dx, 20, lw - startL); nh = clamp(startH + dy, 20, lh - startT);
            } else if (dragging === "nw") {
                nw = clamp(startW - dx, 20, startL + startW); nh = clamp(startH - dy, 20, startT + startH);
                nl = startL + startW - nw; nt = startT + startH - nh;
            } else if (dragging === "ne") {
                nl = startL; nw = clamp(startW + dx, 20, lw - startL);
                nh = clamp(startH - dy, 20, startT + startH); nt = startT + startH - nh;
            } else if (dragging === "sw") {
                nw = clamp(startW - dx, 20, startL + startW); nl = startL + startW - nw;
                nt = startT; nh = clamp(startH + dy, 20, lh - startT);
            }

            el.style.left = (nl / lw * 100) + "%";
            el.style.top = (nt / lh * 100) + "%";
            el.style.width = (nw / lw * 100) + "%";
            el.style.height = (nh / lh * 100) + "%";
        }

        function onUp() {
            document.removeEventListener("mousemove", onMove);
            document.removeEventListener("mouseup", onUp);
            var boxes = pageBoxes[1] || [];
            var b = boxes.find(function(x) { return x.id === boxId; });
            if (b) {
                b.x_pct = parseFloat(el.style.left);
                b.y_pct = parseFloat(el.style.top);
                b.w_pct = parseFloat(el.style.width);
                b.h_pct = parseFloat(el.style.height);
            }
        }

        document.addEventListener("mousemove", onMove);
        document.addEventListener("mouseup", onUp);
    });
}

function deleteBox(boxId) {
    var boxes = pageBoxes[1] || [];
    pageBoxes[1] = boxes.filter(function(b) { return b.id !== boxId; });
    var el = document.getElementById("plot-box-" + boxId);
    if (el) el.remove();
    var resultEl = document.getElementById("result-box-" + boxId);
    if (resultEl) resultEl.remove();
    _renumberBoxes();
}

// ── Extraction ──

async function extractBox(boxId) {
    var boxes = pageBoxes[1] || [];
    var b = boxes.find(function(x) { return x.id === boxId; });
    if (!b || !currentImageId) return;

    var el = document.getElementById("plot-box-" + boxId);
    if (el) el.classList.add("extracting");

    try {
        var resp = await api("/api/extract-image-region", {
            method: "POST",
            body: JSON.stringify({
                image_id: currentImageId,
                crop_x_pct: b.x_pct,
                crop_y_pct: b.y_pct,
                crop_w_pct: b.w_pct,
                crop_h_pct: b.h_pct,
            }),
        });
        var results = await resp.json();
        extractedResults = extractedResults.concat(results);

        if (el) { el.classList.remove("extracting"); el.classList.add("done"); }

        var container = document.getElementById("page-results");
        var oldResult = document.getElementById("result-box-" + boxId);
        if (oldResult) oldResult.remove();

        results.forEach(function(result, ri) {
            var rIdx = extractedResults.length - results.length + ri;
            var div = document.createElement("div");
            div.className = "border border-gray-200 rounded-lg p-4";
            div.id = "result-box-" + boxId;
            div.innerHTML = '<p class="text-xs font-semibold text-purple-600 mb-1">Region #' + b.num + '</p>' +
                _buildResultHtml(result, rIdx);
            container.appendChild(div);
        });

        if (extractedResults.length > 0) show("results-section");
    } catch (e) {
        if (el) el.classList.remove("extracting");
        setStatus("viewer-status", "Extraction failed: " + e.message, "error");
    }
}

// ── Result HTML ──

function _buildResultHtml(result, idx) {
    var confIcon = { high: "\u2713", medium: "~", low: "?" }[result.confidence] || "?";
    var pointCount = 0;
    if (result.series) result.series.forEach(function(s) { pointCount += s.y_values.length; });

    var displayTitle = result.title || result.figure_id;

    var html = '<div class="border-t border-gray-100 pt-3 mt-3">';
    html += '<div class="flex justify-between items-start mb-2"><div>' +
        '<p class="text-sm font-medium text-gray-700">' + displayTitle + '</p>' +
        '<p class="text-xs text-gray-500">' +
            '<span class="badge badge-' + result.plot_type + '">' + result.plot_type.replace("_", " ") + '</span>' +
            ' <span class="ml-1 confidence-' + result.confidence + '">' + confIcon + " " + result.confidence + '</span>' +
            ' <span class="ml-1 text-gray-400">' + (result.series ? result.series.length + " series, " + pointCount + " pts" : "") + '</span>' +
        '</p>' +
        (result.x_label ? '<p class="text-xs text-gray-400">X: ' + result.x_label + (result.x_unit ? " (" + result.x_unit + ")" : "") + (result.x_scale === "log" ? ' <span class="text-purple-500 font-medium">[log]</span>' : "") + " / Y: " + (result.y_label || "") + (result.y_unit ? " (" + result.y_unit + ")" : "") + (result.y_scale === "log" ? ' <span class="text-purple-500 font-medium">[log]</span>' : "") + '</p>' : "") +
        (result.notes ? '<p class="text-xs text-amber-600 mt-1">' + result.notes + '</p>' : "") +
    '</div>' +
    '<button onclick="downloadFigureCsv(' + idx + ')" class="text-xs bg-gray-100 px-2 py-1 rounded hover:bg-gray-200 flex-shrink-0">CSV</button></div>';

    // Context: legend text and text mentions
    if (result.legend_text || (result.text_mentions && result.text_mentions.length > 0)) {
        html += '<details class="mb-2"><summary class="text-xs text-blue-600 cursor-pointer hover:text-blue-800">Extraction context (legend &amp; mentions)</summary>';
        html += '<div class="mt-1 p-2 bg-blue-50 rounded text-xs text-gray-700 space-y-1">';
        if (result.legend_text) {
            html += '<p><span class="font-semibold text-gray-800">Legend:</span> ' + result.legend_text + '</p>';
        }
        if (result.text_mentions && result.text_mentions.length > 0) {
            html += '<div class="mt-1"><span class="font-semibold text-gray-800">Mentioned in text:</span><ul class="list-disc ml-4 mt-0.5">';
            result.text_mentions.forEach(function(m) {
                html += '<li>' + m + '</li>';
            });
            html += '</ul></div>';
        }
        html += '</div></details>';
    }

    // Check if ANY series has actual error bar values (not all null)
    var hasErrors = false;
    if (result.series) {
        result.series.forEach(function(s) {
            if (s.error_bars_lower && s.error_bars_lower.some(function(v) { return v != null && v !== 0; })) hasErrors = true;
            if (s.error_bars_upper && s.error_bars_upper.some(function(v) { return v != null && v !== 0; })) hasErrors = true;
        });
    }
    var ncols = hasErrors ? 4 : 2;

    // Data table — full width, no max-height
    html += '<div class="overflow-auto"><table class="data-table w-full text-sm border-collapse">';
    if (result.series && result.series.length > 0) {
        result.series.forEach(function(series) {
            html += '<tr class="bg-gray-50"><td colspan="' + ncols + '" class="font-medium py-1 px-2 text-gray-700 text-xs">' + series.name + '</td></tr>';
            html += '<tr class="border-b border-gray-200 text-gray-500 text-xs bg-white">' +
                '<th class="py-1 px-2 text-left">' + (result.x_label || "X") + '</th>' +
                '<th class="py-1 px-2 text-left">' + (result.y_label || "Y") + '</th>';
            if (hasErrors) html += '<th class="py-1 px-2 text-left">Err-</th><th class="py-1 px-2 text-left">Err+</th>';
            html += '</tr>';
            for (var i = 0; i < series.y_values.length; i++) {
                var x = i < series.x_values.length ? series.x_values[i] : "";
                var y = series.y_values[i];
                html += '<tr class="border-b border-gray-100">' +
                    '<td class="px-1"><input value="' + x + '" data-result="' + idx + '" data-series="' + series.name + '" data-row="' + i + '" data-col="x" onchange="updateCell(this)"></td>' +
                    '<td class="px-1"><input value="' + y + '" data-result="' + idx + '" data-series="' + series.name + '" data-row="' + i + '" data-col="y" onchange="updateCell(this)"></td>';
                if (hasErrors) {
                    var el = (series.error_bars_lower && i < series.error_bars_lower.length && series.error_bars_lower[i] != null) ? series.error_bars_lower[i] : "";
                    var eu = (series.error_bars_upper && i < series.error_bars_upper.length && series.error_bars_upper[i] != null) ? series.error_bars_upper[i] : "";
                    html += '<td class="px-1"><input value="' + el + '" data-result="' + idx + '" data-series="' + series.name + '" data-row="' + i + '" data-col="el" onchange="updateCell(this)"></td>' +
                        '<td class="px-1"><input value="' + eu + '" data-result="' + idx + '" data-series="' + series.name + '" data-row="' + i + '" data-col="eu" onchange="updateCell(this)"></td>';
                }
                html += '</tr>';
            }
        });
    } else {
        html += '<tr><td colspan="' + ncols + '" class="py-4 text-center text-gray-400">No data extracted</td></tr>';
    }
    html += '</table></div>';

    // SVG preview — below table, not side-by-side
    if (result.series && result.series.length > 0) {
        html += '<div class="mt-3">' + makeDataPlotSvg(result) + '</div>';
    }
    html += '</div>';
    return html;
}

// ── Cell editing ──

function updateCell(input) {
    var idx = parseInt(input.dataset.result);
    var seriesName = input.dataset.series;
    var row = parseInt(input.dataset.row);
    var col = input.dataset.col;
    var val = input.value === "" ? null : isNaN(input.value) ? input.value : parseFloat(input.value);

    var result = extractedResults[idx];
    if (!result) return;
    var seriesIdx = -1;
    var series = null;
    for (var i = 0; i < result.series.length; i++) {
        if (result.series[i].name === seriesName) { series = result.series[i]; seriesIdx = i; break; }
    }
    if (!series || seriesIdx < 0) return;

    // Update local data
    var field;
    if (col === "x") { series.x_values[row] = val; field = "x_values"; }
    else if (col === "y") { series.y_values[row] = val; field = "y_values"; }
    else if (col === "el") { series.error_bars_lower[row] = val; field = "error_bars_lower"; }
    else if (col === "eu") { series.error_bars_upper[row] = val; field = "error_bars_upper"; }

    // Save to server (fire-and-forget)
    if (currentImageId && field) {
        fetch("/api/edit-cell", {
            method: "PATCH",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify({
                image_id: currentImageId,
                result_index: idx,
                series_index: seriesIdx,
                field: field,
                point_index: row,
                value: val
            })
        }).catch(function() {});
    }

    // Refresh SVG preview
    var svgContainer = input.closest(".border-t").querySelector(".mt-3");
    if (svgContainer) {
        svgContainer.innerHTML = makeDataPlotSvg(result);
    }

    // Visual feedback
    input.style.backgroundColor = "#fef3c7";
    setTimeout(function() { input.style.backgroundColor = ""; }, 600);
}

// ── SVG Preview ──

function makeDataPlotSvg(result) {
    if (!result.series || result.series.length === 0) return "";

    var PADT = 20, PAD_LEFT = 55, PAD_RIGHT = 20;
    var colors = ["#3b82f6", "#ef4444", "#10b981", "#f59e0b", "#8b5cf6", "#ec4899"];

    // Collect all y values (including error extents)
    var allY = [];
    var isCategory = false;
    result.series.forEach(function(s) {
        if (!s.y_values) return;
        s.y_values.forEach(function(v) { if (v != null && !isNaN(v)) allY.push(v); });
        if (s.x_values && s.x_values.length > 0 && typeof s.x_values[0] === "string") isCategory = true;
        if (s.error_bars_lower) s.error_bars_lower.forEach(function(e, i) { if (e && s.y_values[i] != null) allY.push(s.y_values[i] - e); });
        if (s.error_bars_upper) s.error_bars_upper.forEach(function(e, i) { if (e && s.y_values[i] != null) allY.push(s.y_values[i] + e); });
    });
    // Force category routing for plot types handled in the categorical branch
    var categoricalTypes = ["histogram", "violin", "forest", "waterfall", "heatmap"];
    var continuousTypes = ["volcano", "funnel", "bland_altman", "bubble", "manhattan"];
    if (categoricalTypes.indexOf(result.plot_type) >= 0) isCategory = true;
    if (continuousTypes.indexOf(result.plot_type) >= 0) isCategory = false;
    // Stacked bar: compute cumulative totals per category for correct y range
    if (result.plot_type === "stacked_bar") {
        var stackTotals = {};
        result.series.forEach(function(s) {
            if (!s.x_values) return;
            for (var i = 0; i < s.x_values.length; i++) {
                var cat = String(s.x_values[i]);
                if (!stackTotals[cat]) stackTotals[cat] = 0;
                stackTotals[cat] += (s.y_values[i] || 0);
            }
        });
        Object.keys(stackTotals).forEach(function(k) { allY.push(stackTotals[k]); });
        allY.push(0); // stacked bars always start at 0
    }
    // Box plot whisker/outlier extents from notes
    if (result.plot_type === "box" && result.notes) {
        result.notes.split(/\n/).forEach(function(line) {
            var m = line.match(/min\s*=\s*([\d.\-]+)/i);
            if (m) allY.push(parseFloat(m[1]));
            m = line.match(/max\s*=\s*([\d.\-]+)/i);
            if (m) allY.push(parseFloat(m[1]));
            var om = line.match(/outliers?\s*:\s*\[([^\]]*)\]/i);
            if (om) om[1].split(",").forEach(function(v) { var n = parseFloat(v.trim()); if (!isNaN(n)) allY.push(n); });
        });
    }
    // Violin plot whisker extents from notes (same format as box)
    if (result.plot_type === "violin" && result.notes) {
        result.notes.split(/\n/).forEach(function(line) {
            var m = line.match(/min\s*=\s*([\d.\-]+)/i);
            if (m) allY.push(parseFloat(m[1]));
            m = line.match(/max\s*=\s*([\d.\-]+)/i);
            if (m) allY.push(parseFloat(m[1]));
        });
    }
    // Waterfall: ensure y range includes 0 and RECIST thresholds if close
    if (result.plot_type === "waterfall") {
        allY.push(0);
    }
    // Forest plot: y_values are effect sizes, include CI extents
    if (result.plot_type === "forest") {
        result.series.forEach(function(s) {
            if (!s.y_values) return;
            for (var i = 0; i < s.y_values.length; i++) {
                if (s.y_values[i] == null) continue;
                var eLo = (s.error_bars_lower && s.error_bars_lower[i]) || 0;
                var eHi = (s.error_bars_upper && s.error_bars_upper[i]) || 0;
                allY.push(s.y_values[i] - eLo);
                allY.push(s.y_values[i] + eHi);
            }
        });
    }
    allY = allY.filter(function(v) { return isFinite(v); });
    if (allY.length === 0) return "";

    // Collect x labels for category plots
    var catLabels = [];
    if (isCategory) {
        if (result.plot_type === "box") {
            result.series.forEach(function(s) { catLabels.push(s.name || ""); });
        } else if (result.series[0] && result.series[0].x_values) {
            result.series[0].x_values.forEach(function(v) { catLabels.push(String(v)); });
        }
    }

    // ── Adaptive dimensions based on data shape ──
    var nCatLabels = catLabels.length;
    // Unique categories for bar/box
    var uniqueCats = [];
    if (isCategory && result.series[0]) {
        var seen = {};
        result.series[0].x_values.forEach(function(v) {
            var k = String(v);
            if (!seen[k]) { seen[k] = true; uniqueCats.push(k); }
        });
    }
    var nUniqueCats = uniqueCats.length || nCatLabels;
    var nSeries = result.series.length;
    var totalPoints = 0;
    result.series.forEach(function(s) { totalPoints += s.y_values.length; });

    var W;
    if (result.plot_type === "box") {
        // Box: scale with number of groups
        W = Math.max(250, Math.min(600, 80 + nSeries * 70));
    } else if (isCategory) {
        // Bar: scale with unique categories × series
        var barsPerCat = nSeries;
        W = Math.max(250, Math.min(700, 80 + nUniqueCats * barsPerCat * 22 + nUniqueCats * 20));
    } else {
        // Scatter/line: roughly 4:3, wider if many points
        W = Math.max(300, Math.min(600, 350 + totalPoints * 0.5));
    }
    W = Math.round(W);

    var plotW = W - PAD_LEFT - PAD_RIGHT;

    // Label angle detection
    var maxLabelLen = 0;
    catLabels.forEach(function(l) { if (l.length > maxLabelLen) maxLabelLen = l.length; });
    var labelCharWidth = 4.5;
    var labelPixelWidth = maxLabelLen * labelCharWidth;
    var slotWidth = nCatLabels > 0 ? plotW / nCatLabels : plotW;
    var angleLabels = isCategory && nCatLabels > 0 && labelPixelWidth > slotWidth * 0.85;
    var PADB_LABEL = angleLabels ? Math.min(50, maxLabelLen * 3.2 + 8) : 18;
    var PADB_AXIS = (result.x_label ? 14 : 0);
    var PADB = PADB_LABEL + PADB_AXIS + 4;

    // Plot height: aim for roughly square plot area, clamped
    var targetPlotH = Math.max(150, Math.min(400, plotW * 0.75));
    var H = PADT + targetPlotH + PADB;
    var plotH = targetPlotH;

    // Scale detection
    function isLogScale(scaleStr) {
        if (!scaleStr) return false;
        var s = scaleStr.toLowerCase();
        return s === "log" || s === "log10" || s === "log2" || s === "logarithmic" || s.indexOf("log") === 0;
    }
    var yLog = isLogScale(result.y_scale);
    var xLog = isLogScale(result.x_scale);
    if (!yLog && result.y_min > 0 && result.y_max > 0 && result.y_max / result.y_min >= 1000) yLog = true;
    if (!xLog && result.x_min > 0 && result.x_max > 0 && result.x_max / result.x_min >= 1000) xLog = true;

    // Y range
    var yMin, yMax;
    if (result.y_min != null && result.y_max != null) {
        yMin = result.y_min; yMax = result.y_max;
    } else {
        yMin = Math.min.apply(null, allY); yMax = Math.max.apply(null, allY);
        if (yMin === yMax) { yMin -= 1; yMax += 1; }
        var yPad = (yMax - yMin) * 0.1;
        yMin -= yPad; yMax += yPad;
    }
    if (yLog && yMin <= 0) yMin = Math.min.apply(null, allY.filter(function(v) { return v > 0; })) * 0.5 || 0.01;
    var yRange = yMax - yMin;
    if (yRange === 0) { yRange = 2; yMin -= 1; yMax += 1; }

    function yToPixel(v) {
        if (yLog) {
            if (v <= 0) v = yMin;
            var logMin = Math.log10(yMin), logMax = Math.log10(yMax);
            var logRange = logMax - logMin;
            if (logRange === 0) logRange = 1;
            return PADT + plotH - ((Math.log10(v) - logMin) / logRange) * plotH;
        }
        return PADT + plotH - ((v - yMin) / yRange) * plotH;
    }

    function fmtTick(v, range) {
        if (range >= 100) return Math.round(v).toString();
        if (range >= 1) return v.toFixed(1);
        if (range >= 0.1) return v.toFixed(2);
        return v.toPrecision(3);
    }
    function fmtLogTick(v) {
        if (v >= 1 && v === Math.round(v)) return v.toString();
        if (v >= 0.01) return v.toPrecision(2);
        return v.toExponential(0);
    }

    // SVG header
    var parts = [];
    parts.push('<svg width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" xmlns="http://www.w3.org/2000/svg" style="max-width:100%">');
    parts.push('<rect width="' + W + '" height="' + H + '" fill="#fafafa" rx="4"/>');

    // Axes
    parts.push('<line x1="' + PAD_LEFT + '" y1="' + PADT + '" x2="' + PAD_LEFT + '" y2="' + (PADT + plotH) + '" stroke="#d1d5db"/>');
    parts.push('<line x1="' + PAD_LEFT + '" y1="' + (PADT + plotH) + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + (PADT + plotH) + '" stroke="#d1d5db"/>');

    // Y ticks
    if (yLog) {
        var logYMin = Math.log10(yMin), logYMax = Math.log10(yMax);
        var startPow = Math.floor(logYMin), endPow = Math.ceil(logYMax);
        for (var p = startPow; p <= endPow; p++) {
            var tv = Math.pow(10, p);
            if (tv < yMin * 0.9 || tv > yMax * 1.1) continue;
            var ty = yToPixel(tv);
            parts.push('<text x="' + (PAD_LEFT - 5) + '" y="' + (ty + 3) + '" text-anchor="end" font-size="8" fill="#9ca3af">' + fmtLogTick(tv) + '</text>');
            parts.push('<line x1="' + PAD_LEFT + '" y1="' + ty + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + ty + '" stroke="#f3f4f6"/>');
        }
    } else {
        for (var t = 0; t <= 4; t++) {
            var tv = yMin + (yRange * t / 4);
            var ty = yToPixel(tv);
            parts.push('<text x="' + (PAD_LEFT - 5) + '" y="' + (ty + 3) + '" text-anchor="end" font-size="8" fill="#9ca3af">' + fmtTick(tv, yRange) + '</text>');
            parts.push('<line x1="' + PAD_LEFT + '" y1="' + ty + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + ty + '" stroke="#f3f4f6"/>');
        }
    }

    // Category labels helper
    var xLabelBaseY = PADT + plotH + 3;
    function renderCatLabels(labels, positions) {
        for (var i = 0; i < labels.length; i++) {
            var lbl = labels[i];
            if (lbl.length > 18) lbl = lbl.substring(0, 16) + "\u2026";
            var lx = positions[i];
            if (angleLabels) {
                parts.push('<text x="' + lx + '" y="' + (xLabelBaseY + 6) + '" text-anchor="end" font-size="8" fill="#6b7280" transform="rotate(-35,' + lx + ',' + (xLabelBaseY + 6) + ')">' + lbl + '</text>');
            } else {
                parts.push('<text x="' + lx + '" y="' + (xLabelBaseY + 10) + '" text-anchor="middle" font-size="8" fill="#6b7280">' + lbl + '</text>');
            }
        }
    }

    // Plot type routing
    if (result.plot_type === "box") {
        // Box-and-whisker
        var nGroups = result.series.length;
        var groupW = plotW / nGroups;
        var boxW = Math.min(40, groupW - 10);

        var boxStats = {};
        var outlierData = {};
        if (result.notes) {
            result.notes.split(/\n/).forEach(function(line) {
                var statsMatch = line.match(/^(.+?):\s*min\s*=\s*([\d.\-]+)\s*,\s*Q1\s*=\s*([\d.\-]+)\s*,\s*median\s*=\s*([\d.\-]+)\s*,\s*Q3\s*=\s*([\d.\-]+)\s*,\s*max\s*=\s*([\d.\-]+)/i);
                if (statsMatch) {
                    boxStats[statsMatch[1].trim()] = {
                        min: parseFloat(statsMatch[2]), q1: parseFloat(statsMatch[3]),
                        median: parseFloat(statsMatch[4]), q3: parseFloat(statsMatch[5]),
                        max: parseFloat(statsMatch[6])
                    };
                }
                var outlierMatch = line.match(/^(.+?)\s*outliers?\s*:\s*\[([^\]]*)\]/i);
                if (outlierMatch) {
                    outlierData[outlierMatch[1].trim()] = outlierMatch[2].split(",").map(function(v) { return parseFloat(v.trim()); }).filter(function(v) { return !isNaN(v); });
                }
            });
        }

        var labelPositions = [];
        result.series.forEach(function(s, si) {
            var col = colors[si % colors.length];
            var cx = PAD_LEFT + si * groupW + groupW / 2;
            var bx = cx - boxW / 2;
            labelPositions.push(cx);

            var stats = boxStats[s.name];
            var median, q1, q3, wMin, wMax;
            if (stats) {
                median = stats.median; q1 = stats.q1; q3 = stats.q3;
                wMin = stats.min; wMax = stats.max;
            } else {
                median = s.y_values[0] || 0;
                var eLo = (s.error_bars_lower && s.error_bars_lower[0]) || 0;
                var eHi = (s.error_bars_upper && s.error_bars_upper[0]) || 0;
                wMin = median - eLo; wMax = median + eHi;
                q1 = median - eLo * 0.5; q3 = median + eHi * 0.5;
            }

            var yQ1 = yToPixel(q1), yQ3 = yToPixel(q3);
            var yMed = yToPixel(median), yWMin = yToPixel(wMin), yWMax = yToPixel(wMax);

            parts.push('<rect x="' + bx + '" y="' + yQ3 + '" width="' + boxW + '" height="' + Math.max(1, yQ1 - yQ3) + '" fill="' + col + '" opacity="0.2" stroke="' + col + '" stroke-width="1.5" rx="1"/>');
            parts.push('<line x1="' + bx + '" y1="' + yMed + '" x2="' + (bx + boxW) + '" y2="' + yMed + '" stroke="' + col + '" stroke-width="2"/>');
            parts.push('<line x1="' + cx + '" y1="' + yQ1 + '" x2="' + cx + '" y2="' + yWMin + '" stroke="' + col + '" stroke-width="1"/>');
            parts.push('<line x1="' + (cx - boxW * 0.3) + '" y1="' + yWMin + '" x2="' + (cx + boxW * 0.3) + '" y2="' + yWMin + '" stroke="' + col + '" stroke-width="1"/>');
            parts.push('<line x1="' + cx + '" y1="' + yQ3 + '" x2="' + cx + '" y2="' + yWMax + '" stroke="' + col + '" stroke-width="1"/>');
            parts.push('<line x1="' + (cx - boxW * 0.3) + '" y1="' + yWMax + '" x2="' + (cx + boxW * 0.3) + '" y2="' + yWMax + '" stroke="' + col + '" stroke-width="1"/>');

            (outlierData[s.name] || []).forEach(function(ov) {
                var oy = yToPixel(ov);
                parts.push('<circle cx="' + cx + '" cy="' + oy + '" r="2.5" fill="none" stroke="' + col + '" stroke-width="1" opacity="0.7"/>');
            });
        });
        var boxLabels = result.series.map(function(s) { return s.name || ""; });
        renderCatLabels(boxLabels, labelPositions);

    } else if (isCategory) {
        // Category-based plot routing
        var catOrder = [];
        var catSet = {};
        result.series.forEach(function(s) {
            if (!s.x_values) return;
            s.x_values.forEach(function(v) {
                var k = String(v);
                if (!catSet[k]) { catSet[k] = true; catOrder.push(k); }
            });
        });
        var nCats = catOrder.length;
        var nBarSeries = result.series.length;

        // Check if this is individual-dots data (many y-values per unique category)
        var hasRepeatedCats = false;
        result.series.forEach(function(s) {
            if (s.x_values && s.y_values && s.y_values.length > nCats * 1.5) hasRepeatedCats = true;
        });

        var groupW = plotW / Math.max(1, nCats);

        if (result.plot_type === "stacked_bar") {
            // ── Stacked bar chart ──
            var barW = Math.min(40, Math.max(8, groupW - 12));
            catOrder.forEach(function(cat, ci) {
                var cx = PAD_LEFT + ci * groupW + groupW / 2;
                var bx = cx - barW / 2;
                var yBase = 0; // stack from zero upward
                result.series.forEach(function(s, si) {
                    // Find the value for this category in this series
                    var val = 0;
                    for (var i = 0; i < s.x_values.length; i++) {
                        if (String(s.x_values[i]) === cat) { val = s.y_values[i] || 0; break; }
                    }
                    if (val === 0) { yBase += val; return; }
                    var col = colors[si % colors.length];
                    var segTop = yToPixel(yBase + val);
                    var segBot = yToPixel(yBase);
                    var segH = Math.max(1, segBot - segTop);
                    parts.push('<rect x="' + bx + '" y="' + segTop + '" width="' + barW + '" height="' + segH + '" fill="' + col + '" opacity="0.85" stroke="#fff" stroke-width="0.5"/>');
                    yBase += val;
                });
            });
            var barLabelPositions = [];
            for (var i = 0; i < nCats; i++) barLabelPositions.push(PAD_LEFT + i * groupW + groupW / 2);
            catLabels = catOrder;
            renderCatLabels(catLabels, barLabelPositions);

        } else if (result.plot_type === "error_bar") {
            // ── Error bar plot (means with error bars, no bars) ──
            var seriesSlotW = groupW / (nBarSeries + 1);
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var catIdx = catOrder.indexOf(String(s.x_values[i]));
                    if (catIdx < 0) catIdx = i;
                    var cx = PAD_LEFT + catIdx * groupW + seriesSlotW * (si + 1);
                    var py = yToPixel(s.y_values[i]);
                    // Mean marker
                    parts.push('<circle cx="' + cx + '" cy="' + py + '" r="4" fill="' + col + '"/>');
                    // Error bars
                    var eLo = (s.error_bars_lower && s.error_bars_lower[i]) || 0;
                    var eHi = (s.error_bars_upper && s.error_bars_upper[i]) || 0;
                    if (eHi > 0 || eLo > 0) {
                        var eTop = yToPixel(s.y_values[i] + eHi);
                        var eBot = yToPixel(s.y_values[i] - eLo);
                        parts.push('<line x1="' + cx + '" y1="' + eTop + '" x2="' + cx + '" y2="' + eBot + '" stroke="' + col + '" stroke-width="1.5"/>');
                        parts.push('<line x1="' + (cx - 4) + '" y1="' + eTop + '" x2="' + (cx + 4) + '" y2="' + eTop + '" stroke="' + col + '" stroke-width="1.5"/>');
                        parts.push('<line x1="' + (cx - 4) + '" y1="' + eBot + '" x2="' + (cx + 4) + '" y2="' + eBot + '" stroke="' + col + '" stroke-width="1.5"/>');
                    }
                }
            });
            var barLabelPositions = [];
            for (var i = 0; i < nCats; i++) barLabelPositions.push(PAD_LEFT + i * groupW + groupW / 2);
            catLabels = catOrder;
            renderCatLabels(catLabels, barLabelPositions);

        } else if (hasRepeatedCats || result.plot_type === "dot_strip") {
            // ── Individual data points mode — jittered dots per category ──
            var seriesSlotW = groupW / (nBarSeries + 1);
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                var bycat = {};
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var cat = i < s.x_values.length ? String(s.x_values[i]) : "";
                    if (!bycat[cat]) bycat[cat] = [];
                    bycat[cat].push(s.y_values[i]);
                }
                catOrder.forEach(function(cat, ci) {
                    var pts = bycat[cat] || [];
                    var cx = PAD_LEFT + ci * groupW + seriesSlotW * (si + 1);
                    pts.forEach(function(yv, pi) {
                        var jitter = (pi - (pts.length - 1) / 2) * 2.5;
                        var px = cx + jitter;
                        var py = yToPixel(yv);
                        parts.push('<circle cx="' + px + '" cy="' + py + '" r="3" fill="' + col + '" opacity="0.7"/>');
                    });
                });
            });
            var barLabelPositions = [];
            for (var i = 0; i < nCats; i++) barLabelPositions.push(PAD_LEFT + i * groupW + groupW / 2);
            catLabels = catOrder;
            renderCatLabels(catLabels, barLabelPositions);

        } else if (result.plot_type === "paired") {
            // ── Paired/spaghetti plot — connected before-after ──
            var seriesSlotW = groupW / 2;
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                // Draw connected line through each x position
                var pts = [];
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null) continue;
                    var catIdx = catOrder.indexOf(String(s.x_values[i]));
                    if (catIdx < 0) catIdx = i;
                    var px = PAD_LEFT + catIdx * groupW + groupW / 2;
                    var py = yToPixel(s.y_values[i]);
                    pts.push({ x: px, y: py });
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="3" fill="' + col + '" opacity="0.7"/>');
                }
                if (pts.length > 1) {
                    var pathD = "M" + pts[0].x + "," + pts[0].y;
                    for (var j = 1; j < pts.length; j++) pathD += "L" + pts[j].x + "," + pts[j].y;
                    parts.push('<path d="' + pathD + '" fill="none" stroke="' + col + '" stroke-width="1.2" opacity="0.5"/>');
                }
            });
            var barLabelPositions = [];
            for (var i = 0; i < nCats; i++) barLabelPositions.push(PAD_LEFT + i * groupW + groupW / 2);
            catLabels = catOrder;
            renderCatLabels(catLabels, barLabelPositions);

        } else if (result.plot_type === "histogram") {
            // ── Histogram — bars at bin centers ──
            // Treat each series independently; x_values are bin centers (as strings)
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                var numXs = s.x_values.map(Number).filter(function(v) { return isFinite(v); });
                // Calculate bar width from spacing between consecutive bins
                var barPixelW;
                if (numXs.length > 1) {
                    var gaps = [];
                    for (var g = 1; g < numXs.length; g++) gaps.push(Math.abs(numXs[g] - numXs[g - 1]));
                    gaps.sort(function(a, b) { return a - b; });
                    var binWidth = gaps[0]; // smallest gap = bin width
                    barPixelW = (binWidth / ((numXs[numXs.length - 1] - numXs[0]) || 1)) * plotW * 0.9;
                } else {
                    barPixelW = plotW / 3;
                }
                barPixelW = Math.max(4, Math.min(barPixelW, plotW / numXs.length - 2));
                // We need our own xToPixel for histogram since we're in categorical branch
                var hxMin = Math.min.apply(null, numXs);
                var hxMax = Math.max.apply(null, numXs);
                if (hxMin === hxMax) { hxMin -= 1; hxMax += 1; }
                var hxRange = hxMax - hxMin;
                var hxPad = hxRange * 0.1;
                hxMin -= hxPad; hxMax += hxPad; hxRange = hxMax - hxMin;
                function histXToPixel(v) { return PAD_LEFT + ((v - hxMin) / hxRange) * plotW; }
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var cx = histXToPixel(xv);
                    var by = yToPixel(s.y_values[i]);
                    var baseY = yToPixel(Math.max(yMin, 0));
                    var bh = baseY - by;
                    if (bh < 0) { by = by + bh; bh = -bh; }
                    parts.push('<rect x="' + (cx - barPixelW / 2) + '" y="' + by + '" width="' + barPixelW + '" height="' + Math.max(1, bh) + '" fill="' + col + '" opacity="0.75" stroke="' + col + '" stroke-width="0.5"/>');
                }
                // X tick labels for histogram
                var step = Math.max(1, Math.floor(numXs.length / 6));
                for (var i = 0; i < numXs.length; i += step) {
                    var tx = histXToPixel(numXs[i]);
                    var label = numXs[i] % 1 === 0 ? numXs[i].toString() : numXs[i].toFixed(1);
                    parts.push('<text x="' + tx + '" y="' + (xLabelBaseY + 10) + '" text-anchor="middle" font-size="8" fill="#9ca3af">' + label + '</text>');
                }
            });

        } else if (result.plot_type === "violin") {
            // ── Violin plot — rendered as box plots (box + whiskers + median) ──
            var nGroups = result.series.length;
            var vGroupW = plotW / Math.max(1, nGroups);
            var vBoxW = Math.min(50, vGroupW - 10);

            var violinStats = {};
            if (result.notes) {
                result.notes.split(/\n/).forEach(function(line) {
                    var statsMatch = line.match(/^(.+?):\s*min\s*=\s*([\d.\-]+)\s*,\s*Q1\s*=\s*([\d.\-]+)\s*,\s*median\s*=\s*([\d.\-]+)\s*,\s*Q3\s*=\s*([\d.\-]+)\s*,\s*max\s*=\s*([\d.\-]+)/i);
                    if (statsMatch) {
                        violinStats[statsMatch[1].trim()] = {
                            min: parseFloat(statsMatch[2]), q1: parseFloat(statsMatch[3]),
                            median: parseFloat(statsMatch[4]), q3: parseFloat(statsMatch[5]),
                            max: parseFloat(statsMatch[6])
                        };
                    }
                });
            }

            var vLabelPositions = [];
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                var cx = PAD_LEFT + si * vGroupW + vGroupW / 2;
                var bx = cx - vBoxW / 2;
                vLabelPositions.push(cx);

                var stats = violinStats[s.name];
                var median, q1, q3, wMin, wMax;
                if (stats) {
                    median = stats.median; q1 = stats.q1; q3 = stats.q3;
                    wMin = stats.min; wMax = stats.max;
                } else {
                    // Estimate from raw y_values
                    var sorted = (s.y_values || []).filter(function(v) { return v != null && !isNaN(v); }).sort(function(a, b) { return a - b; });
                    if (sorted.length === 0) return;
                    wMin = sorted[0]; wMax = sorted[sorted.length - 1];
                    q1 = sorted[Math.floor(sorted.length * 0.25)];
                    median = sorted[Math.floor(sorted.length * 0.5)];
                    q3 = sorted[Math.floor(sorted.length * 0.75)];
                }

                var yQ1 = yToPixel(q1), yQ3 = yToPixel(q3);
                var yMed = yToPixel(median), yWMin = yToPixel(wMin), yWMax = yToPixel(wMax);

                // Wider box shape with rounded corners for violin appearance
                parts.push('<rect x="' + bx + '" y="' + yQ3 + '" width="' + vBoxW + '" height="' + Math.max(1, yQ1 - yQ3) + '" fill="' + col + '" opacity="0.25" stroke="' + col + '" stroke-width="1.5" rx="' + Math.min(6, vBoxW / 4) + '"/>');
                // Median line
                parts.push('<line x1="' + bx + '" y1="' + yMed + '" x2="' + (bx + vBoxW) + '" y2="' + yMed + '" stroke="' + col + '" stroke-width="2.5"/>');
                // Whiskers (thinner for violin)
                parts.push('<line x1="' + cx + '" y1="' + yQ1 + '" x2="' + cx + '" y2="' + yWMin + '" stroke="' + col + '" stroke-width="1" stroke-dasharray="2,2"/>');
                parts.push('<line x1="' + cx + '" y1="' + yQ3 + '" x2="' + cx + '" y2="' + yWMax + '" stroke="' + col + '" stroke-width="1" stroke-dasharray="2,2"/>');
                // Whisker caps
                parts.push('<line x1="' + (cx - vBoxW * 0.2) + '" y1="' + yWMin + '" x2="' + (cx + vBoxW * 0.2) + '" y2="' + yWMin + '" stroke="' + col + '" stroke-width="1"/>');
                parts.push('<line x1="' + (cx - vBoxW * 0.2) + '" y1="' + yWMax + '" x2="' + (cx + vBoxW * 0.2) + '" y2="' + yWMax + '" stroke="' + col + '" stroke-width="1"/>');
            });
            var violinLabels = result.series.map(function(s) { return s.name || ""; });
            renderCatLabels(violinLabels, vLabelPositions);

        } else if (result.plot_type === "forest") {
            // ── Forest plot — horizontal: study names on left, effect sizes + CI ──
            // x_values = study names, y_values = effect sizes, error bars = CI
            var studies = [];
            result.series.forEach(function(s) {
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var name = i < s.x_values.length ? String(s.x_values[i]) : ("Study " + (i + 1));
                    var eLo = (s.error_bars_lower && s.error_bars_lower[i]) || 0;
                    var eHi = (s.error_bars_upper && s.error_bars_upper[i]) || 0;
                    var isOverall = /overall|pooled|summary|combined/i.test(name);
                    studies.push({ name: name, est: s.y_values[i], ciLo: s.y_values[i] - eLo, ciHi: s.y_values[i] + eHi, overall: isOverall });
                }
            });
            if (studies.length === 0) { parts.push('</svg>'); return parts.join(""); }

            // X axis = effect size (horizontal), Y axis = studies (vertical)
            var allEst = [];
            studies.forEach(function(st) { allEst.push(st.est); allEst.push(st.ciLo); allEst.push(st.ciHi); });
            var fxMin = Math.min.apply(null, allEst);
            var fxMax = Math.max.apply(null, allEst);
            if (fxMin === fxMax) { fxMin -= 1; fxMax += 1; }
            var fxPad = (fxMax - fxMin) * 0.15;
            fxMin -= fxPad; fxMax += fxPad;
            var fxRange = fxMax - fxMin;
            function forestXToPixel(v) { return PAD_LEFT + ((v - fxMin) / fxRange) * plotW; }

            var nStudies = studies.length;
            var rowH = plotH / Math.max(1, nStudies);

            // Vertical reference line at 0 (or 1 if all estimates > 0, suggesting OR/HR)
            var refVal = (fxMin > 0.5 && allEst.every(function(v) { return v > 0; })) ? 1 : 0;
            var refX = forestXToPixel(refVal);
            parts.push('<line x1="' + refX + '" y1="' + PADT + '" x2="' + refX + '" y2="' + (PADT + plotH) + '" stroke="#9ca3af" stroke-width="1" stroke-dasharray="4,3"/>');

            studies.forEach(function(st, si) {
                var cy = PADT + si * rowH + rowH / 2;
                var estX = forestXToPixel(st.est);
                var ciLoX = forestXToPixel(st.ciLo);
                var ciHiX = forestXToPixel(st.ciHi);
                var col = st.overall ? "#ef4444" : "#3b82f6";

                // CI horizontal line
                parts.push('<line x1="' + ciLoX + '" y1="' + cy + '" x2="' + ciHiX + '" y2="' + cy + '" stroke="' + col + '" stroke-width="1.5"/>');

                if (st.overall) {
                    // Diamond for overall estimate
                    var dW = Math.min(8, rowH * 0.4);
                    var dH = Math.min(6, rowH * 0.35);
                    parts.push('<polygon points="' + (estX - dW) + ',' + cy + ' ' + estX + ',' + (cy - dH) + ' ' + (estX + dW) + ',' + cy + ' ' + estX + ',' + (cy + dH) + '" fill="' + col + '" stroke="' + col + '"/>');
                } else {
                    // Square for individual studies
                    var sqSz = Math.min(7, rowH * 0.35);
                    parts.push('<rect x="' + (estX - sqSz / 2) + '" y="' + (cy - sqSz / 2) + '" width="' + sqSz + '" height="' + sqSz + '" fill="' + col + '"/>');
                }

                // Study name label on left
                var lbl = st.name;
                if (lbl.length > 22) lbl = lbl.substring(0, 20) + "\u2026";
                parts.push('<text x="' + (PAD_LEFT - 5) + '" y="' + (cy + 3) + '" text-anchor="end" font-size="7" fill="#374151">' + lbl + '</text>');
            });

            // X tick labels for effect size
            for (var t = 0; t <= 4; t++) {
                var txv = fxMin + (fxRange * t / 4);
                var txx = forestXToPixel(txv);
                parts.push('<text x="' + txx + '" y="' + (xLabelBaseY + 10) + '" text-anchor="middle" font-size="8" fill="#9ca3af">' + fmtTick(txv, fxRange) + '</text>');
            }

        } else if (result.plot_type === "waterfall") {
            // ── Waterfall — ordered bars from baseline, colored by sign ──
            var wfData = [];
            result.series.forEach(function(s) {
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var label = i < s.x_values.length ? String(s.x_values[i]) : "";
                    wfData.push({ label: label, value: s.y_values[i] });
                }
            });
            // Sort by value (most negative to most positive) for waterfall
            wfData.sort(function(a, b) { return a.value - b.value; });

            var wfN = wfData.length;
            var wfBarW = Math.max(2, Math.min(20, plotW / wfN - 1));
            var wfGroupW = plotW / Math.max(1, wfN);
            var baseY = yToPixel(0);

            // Reference lines at -30% and +20% (RECIST thresholds)
            if (yMin <= -30 && yMax >= -30) {
                var y30 = yToPixel(-30);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + y30 + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + y30 + '" stroke="#10b981" stroke-width="0.8" stroke-dasharray="4,3"/>');
                parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (y30 + 3) + '" font-size="7" fill="#10b981">-30%</text>');
            }
            if (yMin <= 20 && yMax >= 20) {
                var y20 = yToPixel(20);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + y20 + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + y20 + '" stroke="#ef4444" stroke-width="0.8" stroke-dasharray="4,3"/>');
                parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (y20 + 3) + '" font-size="7" fill="#ef4444">+20%</text>');
            }
            // Zero line
            parts.push('<line x1="' + PAD_LEFT + '" y1="' + baseY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + baseY + '" stroke="#6b7280" stroke-width="0.5"/>');

            wfData.forEach(function(d, i) {
                var cx = PAD_LEFT + i * wfGroupW + wfGroupW / 2;
                var by = yToPixel(d.value);
                var bh = baseY - by;
                var col = d.value < 0 ? "#10b981" : "#ef4444"; // green for negative, red for positive
                if (bh < 0) { by = by + bh; bh = -bh; }
                parts.push('<rect x="' + (cx - wfBarW / 2) + '" y="' + by + '" width="' + wfBarW + '" height="' + Math.max(1, bh) + '" fill="' + col + '" opacity="0.85"/>');
            });

        } else if (result.plot_type === "heatmap") {
            // ── Heatmap — grid of colored cells ──
            // Each series = one row, x_values = column labels, y_values = cell values
            var hmRows = [];
            var hmColLabels = [];
            result.series.forEach(function(s, si) {
                if (!s.y_values) return;
                hmRows.push({ name: s.name || ("Row " + (si + 1)), values: s.y_values });
                if (hmColLabels.length === 0 && s.x_values) {
                    hmColLabels = s.x_values.map(String);
                }
            });
            var nRows = hmRows.length;
            var nCols = hmColLabels.length || (hmRows[0] ? hmRows[0].values.length : 0);
            if (nCols === 0 && hmColLabels.length === 0) {
                for (var c = 0; c < nCols; c++) hmColLabels.push(String(c + 1));
            }

            // Find value range for color mapping
            var hmAllVals = [];
            hmRows.forEach(function(r) { r.values.forEach(function(v) { if (v != null && isFinite(v)) hmAllVals.push(v); }); });
            var hmMin = hmAllVals.length > 0 ? Math.min.apply(null, hmAllVals) : 0;
            var hmMax = hmAllVals.length > 0 ? Math.max.apply(null, hmAllVals) : 1;
            if (hmMin === hmMax) { hmMin -= 1; hmMax += 1; }

            function hmColor(v) {
                // Blue-white-red gradient: min=blue, mid=white, max=red
                var t = (v - hmMin) / (hmMax - hmMin);
                t = Math.max(0, Math.min(1, t));
                var r, g, b;
                if (t < 0.5) {
                    var s = t / 0.5; // 0..1 for blue-to-white
                    r = Math.round(59 + (255 - 59) * s);
                    g = Math.round(130 + (255 - 130) * s);
                    b = Math.round(246 + (255 - 246) * s);
                } else {
                    var s = (t - 0.5) / 0.5; // 0..1 for white-to-red
                    r = Math.round(255);
                    g = Math.round(255 - (255 - 68) * s);
                    b = Math.round(255 - (255 - 68) * s);
                }
                return "rgb(" + r + "," + g + "," + b + ")";
            }

            var cellW = plotW / Math.max(1, nCols);
            var cellH = plotH / Math.max(1, nRows);
            var showText = (nRows <= 10 && nCols <= 10);

            hmRows.forEach(function(row, ri) {
                // Row label on left
                var lbl = row.name;
                if (lbl.length > 14) lbl = lbl.substring(0, 12) + "\u2026";
                var ry = PADT + ri * cellH + cellH / 2;
                parts.push('<text x="' + (PAD_LEFT - 3) + '" y="' + (ry + 3) + '" text-anchor="end" font-size="7" fill="#6b7280">' + lbl + '</text>');

                for (var ci = 0; ci < nCols && ci < row.values.length; ci++) {
                    var v = row.values[ci];
                    var cx = PAD_LEFT + ci * cellW;
                    var cy = PADT + ri * cellH;
                    var fill = (v != null && isFinite(v)) ? hmColor(v) : "#f3f4f6";
                    parts.push('<rect x="' + cx + '" y="' + cy + '" width="' + cellW + '" height="' + cellH + '" fill="' + fill + '" stroke="#fff" stroke-width="0.5"/>');
                    if (showText && v != null && isFinite(v)) {
                        var txt = Math.abs(v) >= 100 ? Math.round(v).toString() : v.toFixed(1);
                        parts.push('<text x="' + (cx + cellW / 2) + '" y="' + (cy + cellH / 2 + 3) + '" text-anchor="middle" font-size="7" fill="#374151">' + txt + '</text>');
                    }
                }
            });

            // Column labels at bottom
            var colLabelPositions = [];
            for (var ci = 0; ci < nCols; ci++) colLabelPositions.push(PAD_LEFT + ci * cellW + cellW / 2);
            if (hmColLabels.length > 0) renderCatLabels(hmColLabels, colLabelPositions);

        } else {
            // ── Standard bar mode — one bar per x-value per series ──
            var barW = Math.min(28, Math.max(4, (groupW - 8) / nBarSeries));
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var catIdx = catOrder.indexOf(String(s.x_values[i]));
                    if (catIdx < 0) catIdx = i;
                    var bx = PAD_LEFT + catIdx * groupW + (groupW - nBarSeries * barW) / 2 + si * barW;
                    var by = yToPixel(s.y_values[i]);
                    var bh = yToPixel(Math.max(yMin, 0)) - by;
                    if (bh < 0) { by = by + bh; bh = -bh; }
                    parts.push('<rect x="' + bx + '" y="' + by + '" width="' + (barW - 1) + '" height="' + Math.max(1, bh) + '" fill="' + col + '" opacity="0.8" rx="1"/>');

                    var eLo = (s.error_bars_lower && s.error_bars_lower[i]) || 0;
                    var eHi = (s.error_bars_upper && s.error_bars_upper[i]) || 0;
                    if (eHi > 0 || eLo > 0) {
                        var ecx = bx + barW / 2;
                        var eTop = yToPixel(s.y_values[i] + eHi);
                        var eBot = yToPixel(s.y_values[i] - eLo);
                        parts.push('<line x1="' + ecx + '" y1="' + eTop + '" x2="' + ecx + '" y2="' + eBot + '" stroke="' + col + '" stroke-width="1"/>');
                        parts.push('<line x1="' + (ecx - 3) + '" y1="' + eTop + '" x2="' + (ecx + 3) + '" y2="' + eTop + '" stroke="' + col + '"/>');
                        parts.push('<line x1="' + (ecx - 3) + '" y1="' + eBot + '" x2="' + (ecx + 3) + '" y2="' + eBot + '" stroke="' + col + '"/>');
                    }
                }
            });
            var barLabelPositions = [];
            for (var i = 0; i < nCats; i++) barLabelPositions.push(PAD_LEFT + i * groupW + groupW / 2);
            catLabels = catOrder;
            renderCatLabels(catLabels, barLabelPositions);
        }

    } else {
        // Scatter / Line plot
        var allX = [];
        result.series.forEach(function(s) {
            if (!s.x_values) return;
            s.x_values.forEach(function(v) { var n = Number(v); if (isFinite(n)) allX.push(n); });
        });
        if (allX.length === 0) { parts.push('</svg>'); return parts.join(""); }

        var xMin, xMax;
        if (result.x_min != null && result.x_max != null) {
            xMin = result.x_min; xMax = result.x_max;
        } else {
            xMin = Math.min.apply(null, allX); xMax = Math.max.apply(null, allX);
            if (xMin === xMax) { xMin -= 1; xMax += 1; }
        }
        if (xLog && xMin <= 0) xMin = Math.min.apply(null, allX.filter(function(v) { return v > 0; })) * 0.5 || 0.01;
        var xRange = xMax - xMin;
        if (xRange === 0) { xRange = 2; xMin -= 1; xMax += 1; }

        function xToPixel(v) {
            if (xLog) {
                if (v <= 0) v = xMin;
                var logMin = Math.log10(xMin), logMax = Math.log10(xMax);
                var logRange = logMax - logMin;
                if (logRange === 0) logRange = 1;
                return PAD_LEFT + ((Math.log10(v) - logMin) / logRange) * plotW;
            }
            return PAD_LEFT + ((v - xMin) / xRange) * plotW;
        }

        // X ticks
        if (xLog) {
            var logXMin = Math.log10(xMin), logXMax = Math.log10(xMax);
            var xStartPow = Math.floor(logXMin), xEndPow = Math.ceil(logXMax);
            for (var p = xStartPow; p <= xEndPow; p++) {
                var txv = Math.pow(10, p);
                if (txv < xMin * 0.9 || txv > xMax * 1.1) continue;
                var txx = xToPixel(txv);
                parts.push('<text x="' + txx + '" y="' + (xLabelBaseY + 10) + '" text-anchor="middle" font-size="8" fill="#9ca3af">' + fmtLogTick(txv) + '</text>');
            }
        } else {
            for (var t = 0; t <= 4; t++) {
                var txv = xMin + (xRange * t / 4);
                var txx = PAD_LEFT + (t / 4) * plotW;
                parts.push('<text x="' + txx + '" y="' + (xLabelBaseY + 10) + '" text-anchor="middle" font-size="8" fill="#9ca3af">' + fmtTick(txv, xRange) + '</text>');
            }
        }

        var pt = result.plot_type;

        if (pt === "volcano") {
            // ── Volcano plot — scatter with threshold lines ──
            // Vertical dashed lines at x = ±1 (log2FC), horizontal at y = 1.3 (-log10(0.05))
            var vxThresh = 1; // log2FC threshold
            var vyThresh = 1.3; // -log10(0.05)
            // Parse thresholds from notes if available
            if (result.notes) {
                var fcMatch = result.notes.match(/log2?\s*FC\s*(?:threshold|cutoff)\s*[=:]\s*([\d.]+)/i);
                if (fcMatch) vxThresh = parseFloat(fcMatch[1]);
                var pvMatch = result.notes.match(/(?:p[\-_]?value|significance)\s*(?:threshold|cutoff)\s*[=:]\s*([\d.]+)/i);
                if (pvMatch) vyThresh = parseFloat(pvMatch[1]);
            }
            // Draw threshold lines
            if (xMin <= -vxThresh && xMax >= -vxThresh) {
                var vlx = xToPixel(-vxThresh);
                parts.push('<line x1="' + vlx + '" y1="' + PADT + '" x2="' + vlx + '" y2="' + (PADT + plotH) + '" stroke="#9ca3af" stroke-width="0.8" stroke-dasharray="4,3"/>');
            }
            if (xMin <= vxThresh && xMax >= vxThresh) {
                var vrx = xToPixel(vxThresh);
                parts.push('<line x1="' + vrx + '" y1="' + PADT + '" x2="' + vrx + '" y2="' + (PADT + plotH) + '" stroke="#9ca3af" stroke-width="0.8" stroke-dasharray="4,3"/>');
            }
            if (yMin <= vyThresh && yMax >= vyThresh) {
                var vhy = yToPixel(vyThresh);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + vhy + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + vhy + '" stroke="#9ca3af" stroke-width="0.8" stroke-dasharray="4,3"/>');
            }
            // Plot points colored by series (typically "up", "down", "non-significant")
            var volcanoColors = { "up": "#ef4444", "down": "#3b82f6", "non-significant": "#9ca3af", "ns": "#9ca3af", "not significant": "#9ca3af" };
            result.series.forEach(function(s, si) {
                var sName = (s.name || "").toLowerCase().trim();
                var col = volcanoColors[sName] || colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var px = xToPixel(xv);
                    var py = yToPixel(s.y_values[i]);
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="2.5" fill="' + col + '" opacity="0.7"/>');
                }
            });

        } else if (pt === "funnel") {
            // ── Funnel plot — scatter with inverted y-axis (SE), reference line ──
            // Vertical reference line at pooled estimate (mean of effect sizes or first x value)
            var pooledEst = 0;
            var allFx = [];
            result.series.forEach(function(s) {
                if (!s.x_values) return;
                s.x_values.forEach(function(v) { var n = Number(v); if (isFinite(n)) allFx.push(n); });
            });
            if (allFx.length > 0) pooledEst = allFx.reduce(function(a, b) { return a + b; }, 0) / allFx.length;
            // Parse pooled estimate from notes
            if (result.notes) {
                var poolMatch = result.notes.match(/(?:pooled|overall|summary)\s*(?:estimate|effect)?\s*[=:]\s*([\d.\-]+)/i);
                if (poolMatch) pooledEst = parseFloat(poolMatch[1]);
            }
            var refLineX = xToPixel(pooledEst);
            parts.push('<line x1="' + refLineX + '" y1="' + PADT + '" x2="' + refLineX + '" y2="' + (PADT + plotH) + '" stroke="#6b7280" stroke-width="0.8" stroke-dasharray="4,3"/>');

            // Pseudo-CI funnel lines (95% CI envelope)
            // SE on y-axis (inverted: higher SE = lower on plot = bottom)
            // For a standard funnel, CI bounds = pooled ± 1.96 * SE
            var funnelPts = 20;
            var seMax = yMax; // top of y range = max SE = bottom of plot
            var seMin = Math.max(0, yMin);
            for (var side = -1; side <= 1; side += 2) {
                var fPath = "";
                for (var fi = 0; fi <= funnelPts; fi++) {
                    var se = seMin + (seMax - seMin) * fi / funnelPts;
                    var bound = pooledEst + side * 1.96 * se;
                    if (bound < xMin || bound > xMax) continue;
                    var fx = xToPixel(bound);
                    var fy = yToPixel(se);
                    fPath += (fPath === "" ? "M" : "L") + fx + "," + fy;
                }
                if (fPath) parts.push('<path d="' + fPath + '" fill="none" stroke="#d1d5db" stroke-width="0.8" stroke-dasharray="3,2"/>');
            }

            // Plot points
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var px = xToPixel(xv);
                    var py = yToPixel(s.y_values[i]);
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="3" fill="' + col + '" opacity="0.7"/>');
                }
            });

        } else if (pt === "bland_altman") {
            // ── Bland-Altman — scatter with bias and LoA reference lines ──
            var baBias = 0, baSD = 0;
            // Try to parse bias and SD from notes
            if (result.notes) {
                var biasMatch = result.notes.match(/(?:bias|mean\s*diff(?:erence)?)\s*[=:]\s*([\d.\-]+)/i);
                if (biasMatch) baBias = parseFloat(biasMatch[1]);
                var sdMatch = result.notes.match(/(?:SD|std\s*dev)\s*[=:]\s*([\d.\-]+)/i);
                if (sdMatch) baSD = parseFloat(sdMatch[1]);
            }
            // If not in notes, compute from data
            if (baBias === 0 && baSD === 0) {
                var baYVals = [];
                result.series.forEach(function(s) {
                    if (!s.y_values) return;
                    s.y_values.forEach(function(v) { if (v != null && isFinite(v)) baYVals.push(v); });
                });
                if (baYVals.length > 0) {
                    baBias = baYVals.reduce(function(a, b) { return a + b; }, 0) / baYVals.length;
                    var baVariance = baYVals.reduce(function(a, b) { return a + (b - baBias) * (b - baBias); }, 0) / baYVals.length;
                    baSD = Math.sqrt(baVariance);
                }
            }
            // Draw reference lines
            var biasY = yToPixel(baBias);
            parts.push('<line x1="' + PAD_LEFT + '" y1="' + biasY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + biasY + '" stroke="#3b82f6" stroke-width="1"/>');
            parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (biasY + 3) + '" font-size="7" fill="#3b82f6">Bias</text>');
            var loaUp = baBias + 1.96 * baSD;
            var loaLo = baBias - 1.96 * baSD;
            if (loaUp >= yMin && loaUp <= yMax) {
                var loaUpY = yToPixel(loaUp);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + loaUpY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + loaUpY + '" stroke="#ef4444" stroke-width="0.8" stroke-dasharray="4,3"/>');
                parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (loaUpY + 3) + '" font-size="6" fill="#ef4444">+1.96SD</text>');
            }
            if (loaLo >= yMin && loaLo <= yMax) {
                var loaLoY = yToPixel(loaLo);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + loaLoY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + loaLoY + '" stroke="#ef4444" stroke-width="0.8" stroke-dasharray="4,3"/>');
                parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (loaLoY + 3) + '" font-size="6" fill="#ef4444">-1.96SD</text>');
            }
            // Plot points
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var px = xToPixel(xv);
                    var py = yToPixel(s.y_values[i]);
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="3" fill="' + col + '" opacity="0.7"/>');
                }
            });

        } else if (pt === "bubble") {
            // ── Bubble plot — scatter with variable circle size ──
            // If error_bars data exists, use it as size dimension; otherwise uniform
            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                // Determine bubble sizes: use error_bars_upper as size proxy if available
                var sizes = null;
                if (s.error_bars_upper && s.error_bars_upper.some(function(v) { return v != null && v > 0; })) {
                    sizes = s.error_bars_upper;
                }
                var maxSize = 0;
                if (sizes) sizes.forEach(function(v) { if (v != null && v > maxSize) maxSize = v; });
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var px = xToPixel(xv);
                    var py = yToPixel(s.y_values[i]);
                    var r = 4; // default radius
                    if (sizes && sizes[i] != null && maxSize > 0) {
                        r = 3 + (sizes[i] / maxSize) * 12; // scale between 3 and 15
                    }
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="' + r.toFixed(1) + '" fill="' + col + '" opacity="0.5" stroke="' + col + '" stroke-width="0.5"/>');
                }
            });

        } else if (pt === "manhattan") {
            // ── Manhattan plot — scatter with alternating chromosome colors + threshold ──
            // Significance threshold at y = -log10(5e-8) ≈ 7.3
            var mhThresh = 7.3;
            if (result.notes) {
                var threshMatch = result.notes.match(/(?:threshold|significance)\s*[=:]\s*([\d.]+)/i);
                if (threshMatch) mhThresh = parseFloat(threshMatch[1]);
            }
            if (mhThresh >= yMin && mhThresh <= yMax) {
                var mhThreshY = yToPixel(mhThresh);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + mhThreshY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + mhThreshY + '" stroke="#ef4444" stroke-width="0.8" stroke-dasharray="4,3"/>');
                parts.push('<text x="' + (PAD_LEFT + plotW + 2) + '" y="' + (mhThreshY + 3) + '" font-size="6" fill="#ef4444">5e-8</text>');
            }
            // Suggestive threshold at y = -log10(1e-5) ≈ 5
            var sugThresh = 5;
            if (sugThresh >= yMin && sugThresh <= yMax) {
                var sugThreshY = yToPixel(sugThresh);
                parts.push('<line x1="' + PAD_LEFT + '" y1="' + sugThreshY + '" x2="' + (PAD_LEFT + plotW) + '" y2="' + sugThreshY + '" stroke="#3b82f6" stroke-width="0.5" stroke-dasharray="3,3"/>');
            }
            // Alternating colors by chromosome group (use series or detect from x gaps)
            var mhColors = ["#3b82f6", "#6366f1", "#8b5cf6", "#a855f7", "#3b82f6", "#6366f1"];
            result.series.forEach(function(s, si) {
                var col = mhColors[si % mhColors.length];
                if (!s.x_values || !s.y_values) return;
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i])) continue;
                    var xv = Number(s.x_values[i]);
                    if (!isFinite(xv)) continue;
                    var px = xToPixel(xv);
                    var py = yToPixel(s.y_values[i]);
                    // Points above threshold are highlighted
                    var aboveThresh = s.y_values[i] >= mhThresh;
                    var ptCol = aboveThresh ? "#ef4444" : col;
                    var ptR = aboveThresh ? 3 : 2;
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="' + ptR + '" fill="' + ptCol + '" opacity="' + (aboveThresh ? 0.9 : 0.6) + '"/>');
                }
            });

        } else {
            // ── Default: Scatter / Line plot ──
            var isLine = (pt === "line" || pt === "dose_response" || pt === "roc" || pt === "area");
            var isStep = (pt === "kaplan_meier");

            result.series.forEach(function(s, si) {
                var col = colors[si % colors.length];
                if (!s.x_values || !s.y_values) return;
                var xVals = s.x_values.map(Number);

                var pts = [];
                for (var i = 0; i < s.y_values.length; i++) {
                    if (s.y_values[i] == null || isNaN(s.y_values[i]) || !isFinite(xVals[i])) continue;
                    pts.push({ x: xVals[i], y: s.y_values[i], idx: i });
                }
                pts.sort(function(a, b) { return a.x - b.x; });

                if ((isLine || isStep) && pts.length > 1) {
                    var pathD = "M" + xToPixel(pts[0].x) + "," + yToPixel(pts[0].y);
                    for (var j = 1; j < pts.length; j++) {
                        if (isStep) {
                            // Step function: horizontal then vertical (Kaplan-Meier style)
                            pathD += "H" + xToPixel(pts[j].x);
                            pathD += "V" + yToPixel(pts[j].y);
                        } else {
                            pathD += "L" + xToPixel(pts[j].x) + "," + yToPixel(pts[j].y);
                        }
                    }
                    parts.push('<path d="' + pathD + '" fill="none" stroke="' + col + '" stroke-width="1.5" opacity="0.8"/>');
                }

                pts.forEach(function(p) {
                    var px = xToPixel(p.x);
                    var py = yToPixel(p.y);
                    parts.push('<circle cx="' + px + '" cy="' + py + '" r="' + (isLine ? 2 : 3) + '" fill="' + col + '" opacity="0.8"/>');

                    var eLo = (s.error_bars_lower && s.error_bars_lower[p.idx]) || 0;
                    var eHi = (s.error_bars_upper && s.error_bars_upper[p.idx]) || 0;
                    if (eHi > 0 || eLo > 0) {
                        var eTop = yToPixel(p.y + eHi);
                        var eBot = yToPixel(p.y - eLo);
                        parts.push('<line x1="' + px + '" y1="' + eTop + '" x2="' + px + '" y2="' + eBot + '" stroke="' + col + '" stroke-width="1" opacity="0.5"/>');
                        parts.push('<line x1="' + (px - 2) + '" y1="' + eTop + '" x2="' + (px + 2) + '" y2="' + eTop + '" stroke="' + col + '" opacity="0.5"/>');
                        parts.push('<line x1="' + (px - 2) + '" y1="' + eBot + '" x2="' + (px + 2) + '" y2="' + eBot + '" stroke="' + col + '" opacity="0.5"/>');
                    }
                });
            });
        }
    }

    // Axis labels
    if (result.y_label) {
        var yLabelText = result.y_label + (result.y_unit ? ' (' + result.y_unit + ')' : '');
        if (yLabelText.length > 30) yLabelText = yLabelText.substring(0, 28) + "\u2026";
        parts.push('<text x="12" y="' + (PADT + plotH / 2) + '" text-anchor="middle" font-size="8" fill="#6b7280" transform="rotate(-90,12,' + (PADT + plotH / 2) + ')">' + yLabelText + '</text>');
    }
    if (result.x_label) {
        var xLabelText = result.x_label + (result.x_unit ? ' (' + result.x_unit + ')' : '');
        if (xLabelText.length > 50) xLabelText = xLabelText.substring(0, 48) + "\u2026";
        parts.push('<text x="' + (PAD_LEFT + plotW / 2) + '" y="' + (H - 2) + '" text-anchor="middle" font-size="8" fill="#6b7280">' + xLabelText + '</text>');
    }

    // Legend
    var legendX = PAD_LEFT, legendY = 5;
    result.series.forEach(function(s, si) {
        var name = s.name || ("Series " + (si + 1));
        if (name.length > 16) name = name.substring(0, 14) + "\u2026";
        var itemW = name.length * 5 + 16;
        if (legendX + itemW > W - 10) { legendX = PAD_LEFT; legendY += 12; }
        parts.push('<rect x="' + legendX + '" y="' + legendY + '" width="8" height="8" fill="' + colors[si % colors.length] + '" rx="1"/>');
        parts.push('<text x="' + (legendX + 11) + '" y="' + (legendY + 7) + '" font-size="8" fill="#6b7280">' + name + '</text>');
        legendX += itemW;
    });

    parts.push('</svg>');
    return parts.join("");
}

// ── Export ──

function downloadFigureCsv(idx) {
    var result = extractedResults[idx];
    if (!result) return;

    var hasErrors = false;
    result.series.forEach(function(s) {
        if (s.error_bars_lower && s.error_bars_lower.some(function(v) { return v != null && v !== 0; })) hasErrors = true;
        if (s.error_bars_upper && s.error_bars_upper.some(function(v) { return v != null && v !== 0; })) hasErrors = true;
    });

    var rows = [];
    var header = ["series", result.x_label || "x", result.y_label || "y"];
    if (hasErrors) header.push("error_lower", "error_upper");
    rows.push(header.join(","));

    result.series.forEach(function(s) {
        for (var i = 0; i < s.y_values.length; i++) {
            var x = i < s.x_values.length ? s.x_values[i] : "";
            var y = s.y_values[i];
            if (typeof x === "string" && x.indexOf(",") !== -1) x = '"' + x + '"';
            var row = [s.name, x, y];
            if (hasErrors) {
                var el = (s.error_bars_lower && i < s.error_bars_lower.length && s.error_bars_lower[i] != null) ? s.error_bars_lower[i] : "";
                var eu = (s.error_bars_upper && i < s.error_bars_upper.length && s.error_bars_upper[i] != null) ? s.error_bars_upper[i] : "";
                row.push(el, eu);
            }
            rows.push(row.join(","));
        }
    });
    var blob = new Blob([rows.join("\n")], { type: "text/csv" });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    var fname = (result.title || result.figure_id || "extracted").replace(/[^a-zA-Z0-9_\-]/g, "_").substring(0, 60);
    a.download = fname + ".csv";
    a.click();
    URL.revokeObjectURL(url);
}

function exportAllCsv() {
    if (extractedResults.length === 0) return;
    var rows = [];
    rows.push(["figure_id", "plot_type", "series", "x", "y", "error_lower", "error_upper"].join(","));

    extractedResults.forEach(function(result) {
        result.series.forEach(function(s) {
            for (var i = 0; i < s.y_values.length; i++) {
                var x = i < s.x_values.length ? s.x_values[i] : "";
                if (typeof x === "string" && x.indexOf(",") !== -1) x = '"' + x + '"';
                var el = (s.error_bars_lower && i < s.error_bars_lower.length && s.error_bars_lower[i] != null) ? s.error_bars_lower[i] : "";
                var eu = (s.error_bars_upper && i < s.error_bars_upper.length && s.error_bars_upper[i] != null) ? s.error_bars_upper[i] : "";
                rows.push([result.figure_id || "", result.plot_type || "", s.name, x, s.y_values[i], el, eu].join(","));
            }
        });
    });
    var blob = new Blob([rows.join("\n")], { type: "text/csv" });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = "extracted_data.csv";
    a.click();
    URL.revokeObjectURL(url);
}

function exportAllJson() {
    if (extractedResults.length === 0) return;
    var blob = new Blob([JSON.stringify(extractedResults, null, 2)], { type: "application/json" });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = "extracted_data.json";
    a.click();
    URL.revokeObjectURL(url);
}

function copyPandasCode() {
    if (extractedResults.length === 0) return;
    var code = "import pandas as pd\n\n";

    extractedResults.forEach(function(result, idx) {
        result.series.forEach(function(series) {
            var varName = "df_fig" + idx + "_" + series.name.replace(/\s+/g, "_").replace(/[^a-zA-Z0-9_]/g, "").toLowerCase();
            code += "# " + (result.title || result.figure_id || "Region " + idx) + " \u2014 " + series.name + "\n";
            code += varName + " = pd.DataFrame({\n";
            code += '    "' + (result.x_label || "x") + '": ' + JSON.stringify(series.x_values) + ",\n";
            code += '    "' + (result.y_label || "y") + '": ' + JSON.stringify(series.y_values) + ",\n";
            if (series.error_bars_lower && series.error_bars_lower.some(function(v) { return v !== null; })) {
                code += '    "error_lower": ' + JSON.stringify(series.error_bars_lower) + ",\n";
            }
            if (series.error_bars_upper && series.error_bars_upper.some(function(v) { return v !== null; })) {
                code += '    "error_upper": ' + JSON.stringify(series.error_bars_upper) + ",\n";
            }
            code += "})\n\n";
        });
    });

    navigator.clipboard.writeText(code).then(function() {
        alert("DataFrame code copied to clipboard!");
    });
}
