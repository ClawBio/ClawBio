"""The 170 historical comparisons with DnaSP 6.12 (Supplementary Table S1).

Two checks per row, both against recorded values rather than a new GUI run:
- test_recorded_python_value pins the skill's value recorded in S1 (a regression
  check of this code, not evidence of agreement with DnaSP);
- test_value_agrees_with_dnasp checks the current value against DnaSP's own
  printed figure at DnaSP's printed precision, except for the seven rows in
  DIVERGENT, whose recorded difference from DnaSP is documented there.
"""
from pathlib import Path
from functools import lru_cache
import json
import re
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dnasp as d

FIXTURES = Path(__file__).parent / 'fixtures'
ROWS = json.loads((FIXTURES / 'historical_comparisons.json').read_text(encoding='utf-8'))['rows']


@lru_cache(None)
def case_results(case):
    inputs = FIXTURES / 'inputs'
    options = {}
    files = {'A': 'rp49_ing34_outGUA.fas', 'B': 'rp49_36.nex',
             'C': 'COII_HsaPtr_outPpy.fas', 'D': 'COII_Apes.nex',
             'E': 'DmelOSRegion.nex', 'F': 'Ex_n1.fas', 'G': 'Ex_n1.fas'}
    if case.startswith('H'):
        a = d.load_alignment(inputs / 'rp49_5regions' / f'Region_{case[1:]}.fas')
    elif case == 'I':
        a = d.parse_vcf(inputs / 'vcf/Data_Example_DiploidPhased.vcf').alignments['Scaffold_2']
    else:
        a = d.load_alignment(inputs / files[case])
    if case in {'A', 'C'}:
        name = 'rp49.gua' if case == 'A' else 'Ppy1'
        i = a.names.index(name)
        options['outgroup'] = a.seqs[i]
        a = d.Alignment(a.names[:i] + a.names[i+1:], a.seqs[:i] + a.seqs[i+1:])
    if case == 'C':
        options['genetic_code'] = d.VERTEBRATE_MITOCHONDRIAL_CODE
    if case in {'B', 'D'}:
        name = 'rp49_36_pops.txt' if case == 'B' else 'COII_Apes_HsaPtr_pops.txt'
        options['pop_assignments'] = d.load_pop_file(inputs / name)
    if case == 'G':
        options['hka_loci'] = d.load_hka_file(inputs / 'hka_synthetic.tsv')
    return d.run_analysis(a, analyses=d.VALID_ANALYSES, **options)


def value_for(row):
    r = case_results(row['case']); label = row['label']
    rs = r['global']
    direct = {
        'S': 'S', 'Segregating sites S': 'S', 'η': 'Eta', 'Total mutations η': 'Eta',
        'Total mutations η (Eta)': 'Eta', 'Hap': 'H', 'Haplotypes h': 'H',
        'Hd': 'Hd', 'Haplotype diversity Hd': 'Hd', 'VarHd': 'VarHd', 'Variance of Hd': 'VarHd',
        'π': 'Pi', 'Nucleotide diversity π (per site)': 'Pi', 'θ_k': 'k', 'Avg pairwise differences k': 'k',
        'G+C': 'GC', 'G+C content': 'GC', 'Sites': 'L_total', 'NetSites': 'L_net', 'Net sites used': 'L_net',
        'Tajima D': 'TajimaD', "Tajima's D": 'TajimaD', 'Fu & Li D*': 'FuLiD_star',
        'Fu & Li D* (no outgroup)': 'FuLiD_star', 'Fu & Li F*': 'FuLiF_star',
        'Fu & Li F* (no outgroup)': 'FuLiF_star', 'R2': 'R2', 'Ramos-Onsins & Rozas R2': 'R2'}
    if label in direct: return getattr(rs, direct[label])
    if 'θ_W' in label: return rs.ThetaW_nuc
    if label.startswith('ZnS'): return r['ld'].ZnS
    if label in {'ZA', 'ZZ'}: return getattr(r['ld'], 'Za' if label == 'ZA' else 'ZZ')
    if label.startswith('Rm'): return r['recombination'].Rm
    if label.startswith('Raggedness'): return r['popsize'].raggedness
    if label == 'Mismatch: mean': return r['popsize'].mean
    if label == 'Mismatch: n pairwise comparisons': return r['popsize'].n_pairs
    if label == "Fu's Fs": return r['fufs'].Fs
    if label.startswith('Fu & Li D ('): return r['fuliout'].D
    if label.startswith('Fu & Li F ('): return r['fuliout'].F
    if label.startswith('SFS folded'): return sum(r['sfs'].folded.values())
    if label.startswith('SFS singletons'): return r['sfs'].folded.get(1, 0)
    div = r.get('divergence')
    if label.startswith('π within'):
        pop = label.removeprefix('π within ').replace('O(st)', 'Ost').replace('O(3+4)', 'O3_4')
        # Dataset labels are retained; the printed description uses conventional names.
        if row['case'] == 'B': return div.Pi2 if 'O(st)' in label else div.Pi1
        return div.Pi2 if pop == 'Human' else div.Pi1
    if label.startswith('Dxy'): return div.Dxy
    if label.startswith('Da '): return div.Da
    if label.startswith('Fixed-difference'): return div.n_fixed
    if label.startswith('Shared mutations'): return div.n_shared
    if '(Sx1)' in label: return div.n_private1
    if '(Sx2)' in label: return div.n_private2
    if label == 'reconciliation': return div.n_private1 + div.n_private2 + div.n_shared + div.n_fixed
    if label.startswith('Hudson Fst'): return r['fst'].fst_mean
    if label.startswith('MK '):
        field = next((x for x in ('Pn', 'Ps', 'Dn', 'Ds') if f' {x} ' in label), None)
        if field is None: field = 'NI' if 'Neutrality' in label else 'alpha' if 'α' in label else 'fisher_p'
        return getattr(r['mk'], field)
    if label == 'Ka (dN)': return r['kaks'].Ka
    if label == 'Ks (dS)': return r['kaks'].Ks
    if label == 'ω = Ka/Ks': return r['kaks'].omega
    if label.startswith('Synonymous sites'): return r['kaks'].S_sites
    if label.startswith('Nonsynonymous sites'): return r['kaks'].N_sites
    if 'locusA' in label: return r['hka'].loci_results[0]['theta_hat']
    if 'locusB' in label: return r['hka'].loci_results[1]['theta_hat']
    if label.startswith('T̂'): return r['hka'].T_hat
    if label == 'χ²': return r['hka'].chi2
    if label == 'df': return r['hka'].df
    if label == 'P-value': return r['hka'].p_value
    raise AssertionError(f'Unmapped historical statistic: {label}')


# Rows whose recorded DnaSP figure differs from the skill's, with the reason.
DIVERGENT = {
    ('A', 'Variance of Hd'): 'last displayed digit',
    ('C', 'ω = Ka/Ks'): 'DnaSP divides its two rounded 5-dp values (0.02561/0.83011)',
    ('C', 'Nonsynonymous sites (mean)'): 'last displayed digit (681 - 168.222 = 512.778)',
    ('H1', 'Fu & Li F*'): "DnaSP's multi-alignment output uses the Achaz (2009) variance; "
                          "the skill implements Simonsen et al. (1995)",
    ('H3', 'Fu & Li F*'): 'Achaz (2009) versus Simonsen et al. (1995) variance',
    ('H5', 'Fu & Li F*'): 'Achaz (2009) versus Simonsen et al. (1995) variance',
    ('I', 'Fu & Li F*'): 'Achaz (2009) versus Simonsen et al. (1995) variance',
}
# Cells where DnaSP printed both of its Fu and Li forms: the comparison target is
# the "v5-style" (Simonsen et al. 1995) figure the skill implements, never the
# biallelic-positions (Achaz 2009) figure printed beside it.
DNASP_TARGET = {
    ('A', 'Fu & Li D (with outgroup)'): '-1.77322',
    ('A', 'Fu & Li F (with outgroup)'): '-1.78504',
    ('F', 'Fu & Li D*'): '-0.52807',
    ('F', 'Fu & Li F*'): '-0.41685',
}
_NUMBER = re.compile(r'-?\d+(?:\.\d+)?')


def _dnasp_figures(display):
    """Numbers DnaSP printed for the row, ignoring parenthetical notes such as
    other settings ("13 under the standard code") or derivations."""
    text = display.replace('−', '-')
    previous = None
    while previous != text:
        previous, text = text, re.sub(r'\([^()]*\)', ' ', text)
    return _NUMBER.findall(text)


def test_multi_figure_cells_have_an_explicit_v5_target():
    multi = {(r['case'], r['label']) for r in ROWS
             if (r['case'], r['label']) not in DIVERGENT and len(_dnasp_figures(r['dnasp_display'])) > 1}
    assert multi == set(DNASP_TARGET)
    for row in ROWS:
        key = (row['case'], row['label'])
        if key in DNASP_TARGET:
            assert re.search(re.escape(DNASP_TARGET[key]) + r'\s*\(v5-style\)', row['dnasp_display']), key


def test_divergent_rows_are_exactly_the_documented_ones():
    recorded = {(r['case'], r['label']) for r in ROWS
                if r['historical_outcome'].startswith(('last displayed digit', 'Achaz'))}
    assert recorded == set(DIVERGENT)


@pytest.mark.parametrize('row', [r for r in ROWS if (r['case'], r['label']) not in DIVERGENT],
                         ids=lambda r: r['case'] + ':' + r['label'])
def test_value_agrees_with_dnasp(row):
    value = float(value_for(row))
    figures = _dnasp_figures(row['dnasp_display'])
    key = (row['case'], row['label'])
    target = DNASP_TARGET.get(key) or (figures[0] if len(figures) == 1 else None)
    assert target is not None, ('no single DnaSP figure to compare', row['dnasp_display'])
    decimals = len(target.split('.')[1]) if '.' in target else 0
    assert abs(value - float(target)) <= 0.5 * 10 ** -decimals + 1e-9, (value, row['dnasp_display'])


@pytest.mark.parametrize('row', ROWS, ids=lambda r: r['case'] + ':' + r['label'])
def test_recorded_python_value(row):
    # Regression of this code: the Python value recorded in S1 at its precision.
    shown = row['python_display']
    precision = row.get('python_precision', len(shown.split('.')[1]) if '.' in shown else 0)
    tolerance = 0.500001 * 10 ** -precision if precision else 1e-10
    assert value_for(row) == pytest.approx(float(shown), abs=tolerance, rel=0)
