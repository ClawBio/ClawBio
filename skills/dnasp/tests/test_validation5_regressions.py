"""Round-5 capture regressions; existing historical fixtures remain unchanged."""
import csv
import json
import math
from pathlib import Path
import re
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dnasp as d

FIXTURES = Path(__file__).resolve().parent / 'fixtures'


def ingroup(filename, outgroup=None, crop=None):
    aln = d.load_alignment(FIXTURES / 'inputs' / filename)
    selected = [(n, s[:crop]) for n, s in zip(aln.names, aln.seqs) if n != outgroup]
    return d.Alignment([n for n, _ in selected], [s for _, s in selected])


@pytest.mark.parametrize('filename,outgroup,window,step,count,last,midpoint,S', [
    ('rp49_ing34_outGUA.fas', 'rp49.gua', 100, 25, 69, '1701-1798', 1731, 6),
    ('COII_HsaPtr_outPpy.fas', 'Ppy1', 100, 50, 13, '601-684', 642, 9),
    ('DmelOSRegion.nex', None, 500, 250, 28, '6751-7107', 6926, 11),
])
def test_windows_capture_placement(filename, outgroup, window, step, count, last, midpoint, S):
    """A2_rp49_win100_25.out, C3_COII_win100_50.out, E2_Dmel_win500_250.out.

    CODIGO2.vb emits before testing To2 < nucw; CONTROLE.vb caps both ends.
    Inputs here are byte-identical to validation5/inputs.
    """
    windows = d.run_analysis(ingroup(filename, outgroup), window, step)['windows']
    assert len(windows) == count
    assert (windows[-1].region, windows[-1].midpoint, windows[-1].S) == (last, midpoint, S)
    start, end = map(int, last.split('-'))
    assert windows[-1].L_total == end - start + 1
    if outgroup == 'rp49.gua':
        assert windows[0].midpoint == 37
        assert round(windows[-1].Pi, 5) == 0.01834
        assert round(windows[-1].ThetaW_nuc, 5) == 0.02668
    if outgroup is None:
        assert [w.midpoint for w in windows[:7]] == [1, 251, 501, 751, 1001, 1251, 1501]
        assert all(w.S == w.L_net == 0 for w in windows[:7])
        assert windows[7].region == '1751-2250'
        assert windows[7].midpoint == 2204


@pytest.mark.parametrize('seq,window,step,regions,midpoints', [
    ('A' * 10, 3, 20, ['1-3', '10-10'], [2, 10]),
    ('A' * 10, 30, 2, ['1-10'], [5]),
    ('A' * 10, 4, 3, ['1-4', '4-7', '7-10'], [2, 5, 8]),
    ('A', 1, 1, ['1-1'], [1]),
    ('NNNRRY', 3, 3, ['1-3', '4-6'], [1, 4]),
])
def test_window_boundaries(seq, window, step, regions, midpoints):
    """A2/C3/E2 window captures above; synthetic boundaries for CONTROLE.vb's caps.

    E2_Dmel_win500_250.out also pins retention of windows with no usable sites.
    """
    result = d.run_analysis(d.Alignment(['a', 'b'], [seq, seq]), window, step)
    assert [w.region for w in result['windows']] == regions
    assert [w.midpoint for w in result['windows']] == midpoints


def test_window_exports_and_figure(tmp_path, monkeypatch):
    """A2_rp49_win100_25.out: midpoint 37 must survive TSV, JSON, report and plot."""
    captured = []
    if d.HAS_MPL:
        close = d.plt.close
        def capture(*args, **kwargs):
            if d.plt.get_fignums():
                for ax in d.plt.gcf().axes:
                    if ax.get_title() == 'Sliding Window Analysis':
                        captured.append(list(ax.lines[0].get_xdata()))
            return close(*args, **kwargs)
        monkeypatch.setattr(d.plt, 'close', capture)
    out = tmp_path / 'windows'
    assert d.main(['--input', str(FIXTURES / 'inputs/rp49_ing34_outGUA.fas'),
                   '--outgroup', 'rp49.gua', '--window', '100', '--step', '25',
                   '--output', str(out)]) == 0
    data = json.loads((out / 'summary.json').read_text())
    rows = list(csv.reader((out / 'results.tsv').read_text().splitlines(), delimiter='\t'))
    header = rows[2]
    assert header.count('Midpoint') == 1
    assert all(len(row) == len(header) for row in rows[3:])
    tsv = [dict(zip(header, row)) for row in rows[4:]]
    report = (out / 'report.md').read_text()
    assert '| Region | Midpoint | S |' in report
    assert len(tsv) == len(data['windows']) == 69
    for row, w in zip(tsv, data['windows']):
        assert int(row['Midpoint']) == w['midpoint']
        assert f"| {w['region']} | {w['midpoint']} | {w['S']} |" in report
    if d.HAS_MPL:
        assert captured == [[w['midpoint'] for w in data['windows']]]
    assert 'summary.json' in (out / 'reproducibility/checksums.sha256').read_text()


@pytest.mark.parametrize('filename,outgroup,variance,cv', [
    ('rp49_ing34_outGUA.fas', 'rp49.gua', 40.7780, 0.3987),
    ('DmelOSRegion.nex', None, 7298.4971, 1.3168),
])
def test_mismatch_captures(filename, outgroup, variance, cv):
    """A_rp49_all.out and E_Dmel_all.out, Pairwise No. of Differences panels."""
    result = d.run_analysis(ingroup(filename, outgroup), analyses={'popsize'})['popsize']
    assert round(result.variance, 4) == variance
    assert round(result.cv, 4) == cv


def test_mismatch_sequence_count_correction():
    """A_rp49_all.out formula; synthetic PairwiseDiff.vb 565-567/732 boundary.

    Four sequences produce six distances [1,2,3,1,2,1], mean 5/3, variance 2/3.
    The CV factor uses four sequences, not six pairs.
    """
    result = d.compute_mismatch(['AAA', 'AAT', 'ATT', 'TTT'])
    assert result.mean == pytest.approx(5 / 3)
    assert result.variance == pytest.approx(2 / 3)
    assert result.cv == pytest.approx((17 / 16) * math.sqrt(2 / 3) / (5 / 3))


@pytest.mark.parametrize('seqs,cv', [(['AAA', 'TTT'], 0.0), (['AAA', 'AAA'], None)])
def test_mismatch_two_sequences(seqs, cv):
    """A_rp49_all.out formula; PairwiseDiff.vb explicitly sets n=2 variance to 0."""
    result = d.compute_mismatch(seqs)
    assert result.variance == 0
    assert result.cv == cv


def test_fu_li_segregating_sites_setting():
    """A_rp49_all.out: second Fu/Li and outgroup panels use Segregating sites.

    FULI.vb also changes singleton and external-mutation capping (SingleMut and
    ExternaMut subtraction) in Eta mode. Substituting 92 is an rp49 reproduction,
    not a general S-to-Eta rule.
    """
    aln = ingroup('rp49_ing34_outGUA.fas', 'rp49.gua')
    original = d.load_alignment(FIXTURES / 'inputs/rp49_ing34_outGUA.fas')
    outgroup = original.seqs[original.names.index('rp49.gua')]
    result = d.run_analysis(aln, analyses={'fuliout'}, outgroup=outgroup)
    rs, out = result['global'], result['fuliout']
    assert [round(v, 5) for v in [rs.FuLiD_star, rs.FuLiF_star, out.D, out.F]] == [
        -1.52466, -1.58179, -1.77322, -1.78504]
    assert [round(v, 5) for v in d.fu_li_d_star_f_star(rs.k, 92, 40, 34)] == [-1.41418, -1.53274]


@pytest.mark.parametrize('code', ['standard', 'vertebrate-mitochondrial'])
def test_per_sequence_enc_transcriptions(code, tmp_path, monkeypatch):
    """M_ENC_standard_GUI_tables_transcribed.txt and
    M_ENC_vertebrate-mitochondrial_GUI_tables_transcribed.txt are labelled
    transcriptions, not DnaSP exports. All 11 + 11 values are checked at 3 dp.
    """
    table = (FIXTURES / 'validation5' / f'M_ENC_{code}_GUI_tables_transcribed.txt').read_text()
    expected = {n: float(v) for n, v in re.findall(r'^(Hsa\d+|Ptr\d+)\s+(\d+\.\d+)', table, re.M)}
    assert len(expected) == 11
    original = ingroup('COII_HsaPtr_outPpy.fas', crop=681)
    path = tmp_path / 'coding.fas'
    path.write_text(''.join(f'>{n}\n{s}\n' for n, s in zip(original.names, original.seqs)))
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out = tmp_path / 'enc'
    assert d.main(['--input', str(path), '--outgroup', 'Ppy1', '--analysis', 'codon',
                   '--genetic-code', code, '--output', str(out)]) == 0
    enc = json.loads((out / 'summary.json').read_text())['codon']['per_sequence_enc']
    assert set(enc) == set(expected)
    assert 'Ppy1' not in enc
    assert {n: round(v, 3) for n, v in enc.items()} == expected
    report = (out / 'report.md').read_text()
    assert '| Sequence | ENC |' in report
    for name, value in expected.items():
        assert f'| {name} | {value:.3f} |' in report


def test_enc_undefined_names_and_serialisation(tmp_path, monkeypatch):
    """L_stop_family.out: ENC undefined; preserve named entries as JSON null."""
    monkeypatch.setattr(d, 'HAS_MPL', False)
    path = tmp_path / 'stops.fas'
    path.write_text('>alpha\nAAATAAAAA\n>beta\nAAATAAAAA\n>out\nAAATAAAAA\n')
    out = tmp_path / 'out'
    assert d.main(['--input', str(path), '--outgroup', 'out', '--analysis', 'codon',
                   '--output', str(out)]) == 0
    codon = json.loads((out / 'summary.json').read_text())['codon']
    assert codon['per_sequence_enc'] == {'alpha': None, 'beta': None}
    assert codon['ENC'] is None
    assert '| alpha | n.a. |' in (out / 'report.md').read_text()


def test_vcf_window_plot_labels_snp_index(tmp_path, monkeypatch):
    """VCF windows slide over retained variant columns (DnaSP6_run_notes5.md, I).

    Their midpoints are SNP indices, so the figure axis must not read "bp".
    The alignment-based A2 plot keeps its base-pair label.
    """
    if not d.HAS_MPL:
        pytest.skip('matplotlib not available')
    labels = []
    close = d.plt.close

    def capture(*args, **kwargs):
        if d.plt.get_fignums():
            for ax in d.plt.gcf().axes:
                if ax.get_xlabel():
                    labels.append(ax.get_xlabel())
        return close(*args, **kwargs)
    monkeypatch.setattr(d.plt, 'close', capture)
    out = tmp_path / 'vcf'
    assert d.main(['--vcf', str(FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf'),
                   '--region', 'Scaffold_234', '--window', '2', '--step', '1',
                   '--output', str(out)]) == 0
    window_labels = [x for x in labels if 'variant index' in x or 'Position' in x]
    assert window_labels and all(x.startswith('Retained variant index') for x in window_labels)
    data = json.loads((out / 'summary.json').read_text())
    assert data['variant_sites_only'] is True
    assert all(1 <= w['midpoint'] <= data['global']['L_total'] for w in data['windows'])


@pytest.mark.parametrize('filename,outgroup,pairs', [
    ('ld_help.fas', None, {(1, 2): 0.125, (1, 18): 0.125, (2, 18): 0.250, (1, 4): -0.063, (4, 17): 0.188}),
    ('Ex_n1.fas', None, {(11, 27): -0.125, (11, 17): -0.063, (17, 21): 0.188, (17, 27): -0.125}),
    ('rp49_5regions/Region_2.fas', None, {(54, 63): 0.125}),
])
def test_ld_sign_follows_dnasp_allele_order_at_ties(filename, outgroup, pairs):
    """K_ld_help.out, F_Ex_n1.out and H_region_2.out LD pair grids (D column).

    At a 2/2 tie in these n = 4 files DnaSP's calculo_mas_freq1 keeps sequence
    1's allele as allele 1, so the sign of D depends on that ordering. |D|,
    |D'| and r^2 were already equal; this pins the sign on the captured rows.
    """
    r = d.run_analysis(ingroup(filename, outgroup), analyses={'ld'})['ld']
    got = {(p.site1, p.site2): p.D for p in r.pairs}
    for key, expected in pairs.items():
        assert key in got, key
        assert got[key] == pytest.approx(expected, abs=6e-4), (key, got[key], expected)


def test_ld_tie_rule_leaves_summaries_unchanged():
    """K_ld_help.out: ZnS 0.4667, Za 0.4667, ZZ 0.0000 and Rm 0 are sign-free."""
    r = d.run_analysis(ingroup('ld_help.fas'), analyses={'ld', 'recombination'})
    assert round(r['ld'].ZnS, 4) == 0.4667
    assert round(r['ld'].Za, 4) == 0.4667
    assert abs(r['ld'].ZZ) < 1e-9
    assert r['recombination'].Rm == 0
    # Reordering the sequences may flip signs but never magnitudes.
    aln = ingroup('ld_help.fas')
    reordered = d.Alignment(list(reversed(aln.names)), list(reversed(aln.seqs)))
    r2 = d.run_analysis(reordered, analyses={'ld'})['ld']
    a = {(p.site1, p.site2): abs(p.D) for p in r['ld'].pairs}
    b = {(p.site1, p.site2): abs(p.D) for p in r2.pairs}
    assert a == pytest.approx(b)


def test_result_json_envelope_file(tmp_path, monkeypatch):
    """result.json envelope written by every run (F_Ex_n1.out: n 4, net sites 46, S 10).

    Checked at file level so the standalone package, which has no ClawBio
    runner, exercises the same contract.
    """
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out = tmp_path / 'run'
    assert d.main(['--input', str(FIXTURES / 'inputs/Ex_n1.fas'), '--analysis', 'polymorphism,ld',
                   '--output', str(out)]) == 0
    payload = json.loads((out / 'result.json').read_text(encoding='utf-8'))
    assert payload['skill'] == 'dnasp' and payload['version'] == d.__version__
    assert payload['ok'] is True and payload['status'] == 'ok'
    assert payload['input_checksum'].startswith('sha256:') and payload['datasets'] == {'input': 'Ex_n1.fas'}
    assert (payload['summary']['n'], payload['summary']['L_net'], payload['summary']['S']) == (4, 46, 10)
    assert payload['summary']['analyses_completed'] == ['ld', 'polymorphism']
    summary_json = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    assert payload['data']['global'] == summary_json['global']
    assert {'report.md', 'summary.json', 'results.tsv', 'ld_pairs.tsv'} <= set(payload['data']['artifacts'])
    assert payload['preferred_artifacts'][0] == 'report.md'
    assert payload['chat_summary_lines'][0].startswith('DnaSP Ex_n1.fas: n = 4, net sites 46, S = 10')
    assert 'result.json' in (out / 'reproducibility/checksums.sha256').read_text(encoding='utf-8')


def test_multi_chrom_vcf_root_envelope_file(tmp_path, monkeypatch):
    """Data_Example_DiploidPhased.vcf (I): per-CHROM runs plus a root result.json."""
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out = tmp_path / 'vcf'
    assert d.main(['--vcf', str(FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf'),
                   '--analysis', 'polymorphism', '--output', str(out)]) == 0
    root = json.loads((out / 'result.json').read_text(encoding='utf-8'))
    assert root['summary']['chromosomes'] == 3
    assert {r['chrom'] for r in root['summary']['runs']} == {'Scaffold_2', 'Scaffold_34', 'Scaffold_234'}
    for run in root['data']['runs'].values():
        child = json.loads((out / run['directory'] / 'result.json').read_text(encoding='utf-8'))
        assert child['summary']['n'] == run['summary']['n']
    assert all((out / a).is_file() for a in root['data']['artifacts'])


def test_chat_lines_sanitise_user_controlled_names(tmp_path, monkeypatch):
    """Filenames feed chat_summary_lines; control characters and markup must not pass."""
    monkeypatch.setattr(d, 'HAS_MPL', False)
    src = FIXTURES / 'inputs/Ex_n1.fas'
    evil = tmp_path / 'ex\t**bold**`code`<b>[x](y).fas'
    evil.write_bytes(src.read_bytes())
    out = tmp_path / 'run'
    assert d.main(['--input', str(evil), '--analysis', 'polymorphism', '--output', str(out)]) == 0
    line = json.loads((out / 'result.json').read_text(encoding='utf-8'))['chat_summary_lines'][0]
    assert line.startswith('DnaSP exboldcodebx(y).fas: n = 4')
    assert '\t' not in line and '*' not in line and '`' not in line and '<' not in line
    assert d._display_label('a' * 200) == 'a' * 77 + '...'
    assert d._display_label('\n\x00') == 'input'


def test_multi_chrom_vcf_root_reproducibility_bundle(tmp_path, monkeypatch):
    """Split VCF runs: the root envelope is inside a root bundle whose checksums cover
    every child file, the root replay reruns the whole split run, and each child's
    replay names its own CHROM (Data_Example_DiploidPhased.vcf, three scaffolds)."""
    import hashlib
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out = tmp_path / 'vcf'
    assert d.main(['--vcf', str(FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf'),
                   '--analysis', 'polymorphism', '--output', str(out)]) == 0
    repro = out / 'reproducibility'
    for name in ('commands.sh', 'environment.yml', 'checksums.sha256', 'manifest.json'):
        assert (repro / name).is_file(), name
    checked = 0
    for line in (repro / 'checksums.sha256').read_text(encoding='utf-8').splitlines():
        digest, rel = line.split(maxsplit=1)
        rel = rel.lstrip('*')
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == digest, rel
        checked += 1
    listed = {line.split(maxsplit=1)[1].lstrip('*') for line in (repro / 'checksums.sha256').read_text(encoding='utf-8').splitlines()}
    assert 'result.json' in listed and 'Scaffold_2/result.json' in listed and 'Scaffold_2/report.md' in listed
    assert checked >= 10
    root_manifest = json.loads((repro / 'manifest.json').read_text(encoding='utf-8'))
    assert root_manifest['mode'] == 'multi-chrom' and set(root_manifest['chromosomes']) == {'Scaffold_2', 'Scaffold_34', 'Scaffold_234'}
    assert '--region' not in root_manifest['replay_arguments']
    assert 'result.json' in root_manifest['outputs'] and 'Scaffold_234/summary.json' in root_manifest['outputs']
    for chrom in ('Scaffold_2', 'Scaffold_34', 'Scaffold_234'):
        child = json.loads((out / chrom / 'reproducibility/manifest.json').read_text(encoding='utf-8'))
        args = child['replay_arguments']
        assert args[args.index('--region') + 1] == chrom
    commands = (repro / 'commands.sh').read_text(encoding='utf-8')
    assert '--region' not in commands and 'replay_output' in commands


def test_multi_chrom_directories_stay_inside_output_root(tmp_path, monkeypatch):
    """CHROM tokens "." and ".." are valid VCF values; their run directories must
    stay inside --output and never overwrite the root or its parent."""
    monkeypatch.setattr(d, 'HAS_MPL', False)
    src = (FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf').read_text(encoding='utf-8')
    lines = [('..\t' + l.split('\t', 1)[1]) if l.startswith('Scaffold_2\t')
             else ('.\t' + l.split('\t', 1)[1]) if l.startswith('Scaffold_34\t') else l
             for l in src.splitlines()]
    work = tmp_path / 'work'
    work.mkdir()
    vcf = work / 'dots.vcf'
    vcf.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    sentinel = tmp_path / 'report.md'
    sentinel.write_text('untouched', encoding='utf-8')
    out = work / 'out'
    assert d.main(['--vcf', str(vcf), '--analysis', 'polymorphism', '--output', str(out)]) == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == ['report.md', 'work']
    assert sorted(p.name for p in work.iterdir()) == ['dots.vcf', 'out']
    assert sentinel.read_text(encoding='utf-8') == 'untouched'
    root = json.loads((out / 'result.json').read_text(encoding='utf-8'))
    dirs = {r['chrom']: r['directory'] for r in root['summary']['runs']}
    assert set(dirs) == {'..', '.', 'Scaffold_234'}
    for chrom, rel in dirs.items():
        assert (out / rel).resolve().parent == out.resolve() and '/' not in rel and rel not in ('.', '..')
        assert (out / rel / 'result.json').is_file()
    assert dirs['..'].startswith('chrom_') and dirs['.'].startswith('chrom_') and dirs['..'] != dirs['.']


def test_multi_chrom_directory_names_never_collide(tmp_path, monkeypatch):
    """A CHROM whose literal name equals another CHROM's hash-suffixed directory
    must not share (and overwrite) that directory: three CHROMs, three runs."""
    import hashlib
    monkeypatch.setattr(d, 'HAS_MPL', False)
    suffix = hashlib.sha1(b'A/B').hexdigest()[:6]
    names = [f'A_B_{suffix}', 'A_B', 'A/B']          # ordered so the third collides twice
    src = (FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf').read_text(encoding='utf-8')
    mapping = {'Scaffold_2': names[0], 'Scaffold_34': names[1], 'Scaffold_234': names[2]}
    lines = []
    for line in src.splitlines():
        chrom = line.split('\t', 1)[0]
        lines.append(mapping[chrom] + '\t' + line.split('\t', 1)[1] if chrom in mapping else line)
    vcf = tmp_path / 'collide.vcf'
    vcf.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    out = tmp_path / 'out'
    assert d.main(['--vcf', str(vcf), '--analysis', 'polymorphism', '--output', str(out)]) == 0
    root = json.loads((out / 'result.json').read_text(encoding='utf-8'))
    dirs = {r['chrom']: r['directory'] for r in root['summary']['runs']}
    assert set(dirs) == set(names) and len(set(dirs.values())) == 3
    assert dirs[names[0]] == names[0] and dirs[names[1]] == 'A_B'
    assert dirs[names[2]].startswith(f'A_B_{suffix}') and dirs[names[2]] != names[0]
    sizes = {c: json.loads((out / rel / 'result.json').read_text(encoding='utf-8'))['summary']['n'] for c, rel in dirs.items()}
    assert sizes == {names[0]: 20, names[1]: 18, names[2]: 20}
    assert len(root['data']['runs']) == 3


def _vcf_with_chroms(tmp_path, mapping):
    src = (FIXTURES / 'inputs/vcf/Data_Example_DiploidPhased.vcf').read_text(encoding='utf-8')
    lines = []
    for line in src.splitlines():
        chrom = line.split('\t', 1)[0]
        lines.append(mapping[chrom] + '\t' + line.split('\t', 1)[1] if chrom in mapping else line)
    vcf = tmp_path / 'chroms.vcf'
    vcf.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return vcf


@pytest.mark.parametrize('names', [
    ['chr1', 'CHR1', 'CON'],                    # case-only variants; Windows reserved name
    ['Ω', 'Ω', 'Scaffold_234'],       # Greek capital omega vs Ohm sign (NFKC-equal)
    ['Scaffold_2', 'COM\u00b9', 'LPT\u00b3'],   # superscript-digit device names, after a normal CHROM
])
def test_multi_chrom_directories_distinct_on_case_insensitive_filesystems(tmp_path, monkeypatch, names):
    """CHROM names that differ only by case or Unicode form must get distinct
    directories on macOS/Windows; Windows device names are prefixed."""
    monkeypatch.setattr(d, 'HAS_MPL', False)
    vcf = _vcf_with_chroms(tmp_path, dict(zip(['Scaffold_2', 'Scaffold_34', 'Scaffold_234'], names)))
    out = tmp_path / 'out'
    assert d.main(['--vcf', str(vcf), '--analysis', 'polymorphism', '--output', str(out)]) == 0
    root = json.loads((out / 'result.json').read_text(encoding='utf-8'))
    dirs = {r['chrom']: r['directory'] for r in root['summary']['runs']}
    assert set(dirs) == set(names)
    keys = {d._fs_equivalence_key(rel) for rel in dirs.values()}
    assert len(keys) == 3
    for reserved in ('CON', 'COM\u00b9', 'LPT\u00b3'):
        if reserved in names:
            assert dirs[reserved].startswith('chrom_'), dirs[reserved]
    sizes = {c: json.loads((out / rel / 'result.json').read_text(encoding='utf-8'))['summary']['n'] for c, rel in dirs.items()}
    assert sizes == dict(zip(names, [20, 18, 20]))
