"""v0.5.1 regressions derived from the DnaSP help, VB source and input contract."""
from pathlib import Path
import json
import math
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dnasp as d


def alignment(seqs):
    return d.Alignment([f"s{i}" for i in range(len(seqs))], seqs)


def test_ambiguity_is_missing_consistently():
    r = d.run_analysis(alignment(['AAR', 'AAA', 'AAA', 'AAA']), analyses={'sfs'})
    assert r['global'].S == 0
    assert r['global'].L_net == 2
    assert sum(r['sfs'].folded.values()) == 0
    assert d.compute_fst({'a': ['RRR'] * 2, 'b': ['RRR'] * 2}).fst_pairwise[('a', 'b')] is None


@pytest.mark.parametrize('names,seqs', [(['a', 'a'], ['AAA', 'AAA']), (['a'], ['AAA', 'AAA']), (['a'], ['AAZ'])])
def test_invalid_alignment_rejected(names, seqs):
    with pytest.raises(ValueError):
        d.Alignment(names, seqs)


def test_nexus_quoted_labels_symbols_and_interleave(tmp_path):
    p = tmp_path / 'quoted.nex'
    p.write_text("#NEXUS\nbegin data; dimensions ntax=3 nchar=6;\nformat datatype=dna interleave=yes missing=X gap=~ matchchar=.;\nmatrix\n'a sample' AAA\nb AAT\nc AAC\n\n'a sample' AAX\nb ..~\nc ...\n; end;", encoding='utf-8')
    a = d.parse_nexus(p)
    assert a.names == ['a sample', 'b', 'c']
    assert a.seqs == ['AAAAAN', 'AATAA-', 'AACAAN']


@pytest.mark.parametrize('body', ['a AAA\nb AAT', 'a AAA\nb AA', 'a AAA\nb AAT\nb AAC'])
def test_nexus_dimensions_and_rows_checked(tmp_path, body):
    p = tmp_path / 'bad.nex'
    p.write_text('#NEXUS\nbegin data; dimensions ntax=3 nchar=3; matrix\n' + body + '\n; end;', encoding='utf-8')
    with pytest.raises(ValueError):
        d.parse_nexus(p)


@pytest.mark.parametrize('row', ['a 4 -1 100 20', 'a 4 1 nan 20', 'a 4 1 100 20 0', 'a 4 1 100 20 BAD', 'a 4 broken 100 20', 'a 4'])
def test_hka_rejects_bad_rows(tmp_path, row):
    p = tmp_path / 'loci.txt'
    p.write_text(row + '\nb 4 5 100 10\n', encoding='utf-8')
    with pytest.raises(ValueError, match='line 1'):
        d.load_hka_file(p)


def test_population_lengths_checked():
    with pytest.raises(ValueError):
        d.compute_divergence(['AAA'] * 2, ['AAAAAA'] * 2)


def test_indel_help_model_one():
    # DnaSP help: InDel (Insertion-Deletion) Polymorphism, 13 x 18 example.
    ref = 'AAAAAAGGGGGGGGGGGG'
    patterns = ['..................', '..................', '...C..--..........',
                '..................', '..---.............', '.........-------..',
                '..---....-------..', '.........-------..', '.C...........---..',
                '.........-------..', '.........-------..', '..................', '..................']
    seqs = [''.join(ref[i] if c == '.' else c for i, c in enumerate(p)) for p in patterns]
    r = d.compute_indel(seqs)
    assert r.n_events == 2
    assert r.net_sites == 11
    assert r.n_haplotypes == 3
    assert r.mean_event_length == 2.5
    assert r.mean_deletion_length == pytest.approx(8/3)
    assert r.k_indel == pytest.approx(0.435897435897)
    assert r.pi_indel == pytest.approx(0.039627039627)
    assert r.haplotype_diversity == pytest.approx(0.410256410256)
    assert r.theta_indel == pytest.approx(0.644493786401)
    assert r.tajima_d_indel == pytest.approx(-0.909202, abs=1e-6)
    assert r.n_excluded_events == 2


def test_ld_original_coordinates_and_vb_distance():
    # CODIGO2.vb inline LD loop; same example as linkage-disequilibrium help.
    seqs = ['ATATACGGGGTTA---TTAGA', 'CGATAC--GG-TA---TAACA',
            'AGATACGG-GATA---TAATA', 'ATAAACGGGGATA---GTAGT']
    r = d.run_analysis(alignment(seqs), analyses={'ld', 'recombination'})
    pair = next(p for p in r['ld'].pairs if (p.site1, p.site2) == (1, 18))
    assert pair.dist == 13


def test_indel_missing_fixed_and_adjacent_boundaries():
    # Mod34_BuscoNumEventosEnInDel skips missing columns inside fragments;
    # Mod34_Compute excludes fixed-gap sites from the denominator.
    r = d.compute_indel(['A--A', 'AN-A', 'AAAA', 'AAAA'])
    assert r.n_events == 1 and r.events[0].length == 1
    assert r.net_sites == 3
    r = d.compute_indel(['A-AA', 'AA-A', 'AAAA', 'AAAA'])
    assert r.n_events == 2  # adjacent events are not overlapping
    assert d.compute_indel(['----'] * 4).n_events == 0


def test_stop_family_rscu_and_missing_enc_class():
    # CodonUsage.vb MuestraRSCU includes AA 21; M23ENC interpolates missing F3.
    r = d.compute_codon_usage(['AAATAAAAA'])
    assert r.n_codons == 3
    assert r.rscu['TAA'] == 3
    assert r.rscu['TAG'] == r.rscu['TGA'] == 0
    s = 'TTTTTTGCTGCTCGTCGT'
    assert d.compute_codon_usage([s]).ENC == 20
    assert d.compute_codon_usage([s, s]).ENC == 20


def test_fs_log_odds_extreme_tail():
    # For H=n, upper tail is a single Ewens term; log lower tail tends to zero.
    n, theta = 500, 2.0
    expected = (n - 1) * math.log(theta) - sum(math.log(theta + i) for i in range(1, n))
    r = d.compute_fu_fs(['A'] * n, n, theta)
    assert r.Fs == pytest.approx(expected, abs=1e-9)


def test_faywu_excludes_multiallelic_columns():
    r = d.compute_fay_wu(['AA', 'AA', 'CA', 'GA'], 'AA')
    assert r.n_polarised == 0
    assert r.L_net == 1
    assert r.H is None


@pytest.mark.parametrize('args', [['--outgroup', 'absent'], ['--analysis', 'mk'], ['--analysis', 'hka'],
                                  ['--analysis', 'typo'], ['--window', '-1'], ['--step', '-1'],
                                  ['--analysis', 'codon'], ['--analysis', 'kaks']])
def test_explicit_cli_failure_is_nonzero(tmp_path, args):
    p = tmp_path / 'input.fas'
    p.write_text('>a\nAAAA\n>b\nAAAT\n>c\nAAAT\n>d\nAAAA\n', encoding='utf-8')
    assert d.main(['--input', str(p), '--output', str(tmp_path / 'out'), *args]) != 0


def test_demo_windows_outputs_and_provenance(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out = tmp_path / 'output with spaces'
    assert d.main(['--demo', '--window', '100', '--step', '25', '--output', str(out)]) == 0
    table = (out / 'results.tsv').read_text(encoding='utf-8')
    assert '1-100' in table and '201-300' in table and '1-25\t' not in table
    manifest = json.loads((out / 'reproducibility' / 'manifest.json').read_text(encoding='utf-8'))
    assert manifest['arguments'][:2] == ['--demo', '--window']
    assert 'codon' in manifest['completed']
    # A second run must preserve the first result rather than overwrite it.
    assert d.main(['--demo', '--output', str(out)]) != 0
    assert (out / 'results.tsv').read_text(encoding='utf-8') == table


def test_report_figures_code_and_missing_estimates(tmp_path, monkeypatch):
    if not d.HAS_MPL:
        pytest.skip('matplotlib not installed')
    code = d.VERTEBRATE_MITOCHONDRIAL_CODE
    aln = alignment(['TGATGG'] * 4)
    results = d.run_analysis(aln, analyses={'codon', 'fst'},
                             pop_assignments={'s0': 'a', 's1': 'a', 's2': 'b', 's3': 'b'},
                             genetic_code=code)
    seen = []
    close = d.plt.close
    def capture(*args, **kwargs):
        if d.plt.get_fignums():
            for ax in d.plt.gcf().axes:
                seen.extend(t.get_text() for t in ax.get_xticklabels())
                seen.extend(t.get_text() for t in ax.texts)
        return close(*args, **kwargs)
    monkeypatch.setattr(d.plt, 'close', capture)
    figures = d.make_figures(tmp_path, results)
    assert 'TGA\n(W)' in seen and 'ATA\n(M)' in seen
    assert 'n.a.' in seen
    report = d.write_report(tmp_path, 'test', aln, results, figures, genetic_code=code)
    import re
    for target in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', report.read_text(encoding='utf-8')):
        assert (report.parent / target).is_file()
    assert '| * | AGA |' in report.read_text(encoding='utf-8')
    vcf_report = d.write_report(tmp_path, 'test.vcf', aln, results, figures,
                                genetic_code=code, variant_sites_only=True)
    text = vcf_report.read_text(encoding='utf-8')
    assert 'unaffected' not in text and 'no invariant positions' not in text
    assert 'ascertainment' in text and 'monomorphic' in text


def test_neutrality_interpretation_requires_a_test():
    assert 'significance not assessed' in d._fu_li_interp(-3.0).lower()


def test_hka_only_cli(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    path = tmp_path / 'hka.tsv'
    path.write_text('a 11 10 1000 20\nb 11 5 1000 10\n', encoding='utf-8')
    assert d.main(['--hka-file', str(path), '--analysis', 'hka', '--output', str(tmp_path / 'out')]) == 0
    assert 'HKA' in (tmp_path / 'out/report.md').read_text(encoding='utf-8')


def test_replay_archives_inputs_and_handles_equals(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    path = tmp_path / 'input with spaces.fas'
    path.write_text('>a\nAAA\n>b\nAAT\n>c\nAAA\n>d\nAAT\n', encoding='utf-8')
    out = tmp_path / 'output with spaces'
    assert d.main([f'--input={path}', f'--output={out}']) == 0
    manifest = json.loads((out / 'reproducibility/manifest.json').read_text(encoding='utf-8'))
    assert str(path.resolve()) in manifest['inputs']
    replay = manifest['replay_arguments']
    path.unlink()
    assert d.main(replay) == 0
    assert (out / 'replay_output/results.tsv').is_file()


def test_clawbio_dispatch_accepts_real_dnasp_options(tmp_path, monkeypatch):
    from clawbio.cli import SKILLS
    allowed = SKILLS['dnasp']['allowed_extra_flags']
    assert {'--analysis', '--input2', '--pop-file', '--genetic-code', '--vcf', '--region', '--vcf-merge', '--hka-file'} <= allowed
    assert '--n-sim' not in allowed  # not an implemented DnaSP option
    from clawbio.cli import run_skill
    fixtures = Path(__file__).parent / 'fixtures/inputs'
    (tmp_path / 'input data.vcf').write_bytes((fixtures / 'vcf/Data_Example_DiploidPhased.vcf').read_bytes())
    (tmp_path / 'hka counts.tsv').write_bytes((fixtures / 'hka_synthetic.tsv').read_bytes())
    monkeypatch.chdir(tmp_path)
    for name, arguments in [
        ('vcf', ['--vcf', 'input data.vcf',
                 '--region', 'Scaffold_2', '--analysis', 'polymorphism']),
        ('hka', ['--hka-file=hka counts.tsv', '--analysis', 'hka']),
    ]:
        result = run_skill('dnasp', output_dir=str(tmp_path / name), extra_args=arguments)
        assert result['exit_code'] == 0, result


def test_result_json_envelope_and_runner_promotion(tmp_path):
    """AGENTS.md output contract: result.json envelope, promoted by clawbio.cli.run_skill.

    Ex_n1.fas (F_Ex_n1.out): n 4, net sites 46, S 10. The envelope carries the
    headline summary, the summary.json payload plus artifacts, and the fields
    the runner surfaces (chat_summary_lines, preferred_artifacts).
    """
    from clawbio.cli import run_skill
    fixtures = Path(__file__).parent / 'fixtures/inputs'
    out = tmp_path / 'run'
    result = run_skill('dnasp', input_path=str(fixtures / 'Ex_n1.fas'), output_dir=str(out),
                       extra_args=['--analysis', 'polymorphism,ld'])
    assert result['exit_code'] == 0, result
    payload = result['skill_result_json']
    assert payload['skill'] == 'dnasp' and payload['version'] == d.__version__
    assert payload['ok'] is True and payload['status'] == 'ok'
    assert payload['input_checksum'].startswith('sha256:') and payload['datasets'] == {'input': 'Ex_n1.fas'}
    assert (payload['summary']['n'], payload['summary']['L_net'], payload['summary']['S']) == (4, 46, 10)
    assert payload['summary']['analyses_completed'] == ['ld', 'polymorphism']
    summary_json = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    assert payload['data']['global'] == summary_json['global']
    assert payload['data']['ld']['ZnS'] == summary_json['ld']['ZnS']
    assert {'report.md', 'summary.json', 'results.tsv', 'ld_pairs.tsv'} <= set(payload['data']['artifacts'])
    assert 'result.json' in (out / 'reproducibility/checksums.sha256').read_text(encoding='utf-8')
    assert result['preferred_artifacts'][0] == 'report.md'
    assert result['chat_summary_lines'][0].startswith('DnaSP Ex_n1.fas: n = 4, net sites 46, S = 10')
    assert str(out / 'result.json') == result['result_json_path']


def test_multi_chrom_vcf_root_envelope_promoted(tmp_path):
    """Default multi-CHROM VCF runs write one envelope per CHROM directory plus a
    root result.json that clawbio.cli.run_skill promotes (Data_Example_DiploidPhased.vcf:
    Scaffold_2, Scaffold_34, Scaffold_234; sample sizes 20/18/20 as in DnaSP)."""
    from clawbio.cli import run_skill
    fixtures = Path(__file__).parent / 'fixtures/inputs'
    out = tmp_path / 'vcf'
    result = run_skill('dnasp', output_dir=str(out),
                       extra_args=['--vcf', str(fixtures / 'vcf/Data_Example_DiploidPhased.vcf'),
                                   '--analysis', 'polymorphism'])
    assert result['exit_code'] == 0, result
    payload = result['skill_result_json']
    assert payload['skill'] == 'dnasp' and payload['ok'] is True
    assert payload['summary']['chromosomes'] == 3 and payload['summary']['variant_sites_only'] is True
    sizes = {r['chrom']: r['n'] for r in payload['summary']['runs']}
    assert sizes == {'Scaffold_2': 20, 'Scaffold_34': 18, 'Scaffold_234': 20}
    assert set(result['preferred_artifacts']) == {'Scaffold_2/report.md', 'Scaffold_34/report.md', 'Scaffold_234/report.md'}
    assert all((out / a).is_file() for a in payload['data']['artifacts'])
    assert len(result['chat_summary_lines']) == 4 and result['chat_summary_lines'][0].endswith('3 chromosomes analysed separately.')
    assert (out / 'Scaffold_2/result.json').is_file()
