"""Windows encoding regressions for the dnasp CLI.

On Windows, Python writes redirected stdout and stderr (as when an agent
captures them) in the ANSI code page, usually cp1252, and opens text files in
that code page unless an encoding is given. The console summary uses symbols
such as pi, eta and theta that cp1252 cannot encode, so the run used to exit 1,
and results.tsv carried the input file name in the platform encoding.
"""
from pathlib import Path
import builtins
import io
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dnasp as d

SMALL_FASTA = '>a\nAAAT\n>b\nAATT\n>c\nAAAT\n>d\nATTT\n'


def _cp1252_stream():
    raw = io.BytesIO()
    return raw, io.TextIOWrapper(raw, encoding='cp1252', errors='strict', write_through=True)


def test_cli_completes_on_a_cp1252_console(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    out_raw, out = _cp1252_stream()
    err_raw, err = _cp1252_stream()
    monkeypatch.setattr(sys, 'stdout', out)
    monkeypatch.setattr(sys, 'stderr', err)

    code = d.main(['--demo', '--output', str(tmp_path / 'demo')])

    out.flush()
    err.flush()
    console = out_raw.getvalue().decode('cp1252')
    assert code == 0, err_raw.getvalue().decode('cp1252')
    assert 'Tajima D=' in console
    # Symbols the code page cannot show are escaped, not silently dropped.
    assert '\\u03c0=' in console


def test_console_symbols_are_unchanged_on_a_utf8_console(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    raw = io.BytesIO()
    out = io.TextIOWrapper(raw, encoding='utf-8', errors='strict', write_through=True)
    monkeypatch.setattr(sys, 'stdout', out)

    assert d.main(['--demo', '--output', str(tmp_path / 'demo')]) == 0

    out.flush()
    console = raw.getvalue().decode('utf-8')
    assert 'π=' in console
    assert '\\u03c0' not in console


def test_tsv_outputs_are_written_as_utf8(tmp_path, monkeypatch):
    monkeypatch.setattr(d, 'HAS_MPL', False)
    text_writes = []
    real_open = builtins.open

    def spy(file, mode='r', *args, **kwargs):
        if 'b' not in mode and any(flag in mode for flag in 'wax+'):
            text_writes.append((Path(file).name, kwargs.get('encoding')))
        return real_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(d, 'open', spy, raising=False)
    name = 'población_θ.fas'
    source = tmp_path / name
    source.write_text(SMALL_FASTA, encoding='utf-8')
    out = tmp_path / 'out'

    assert d.main(['--input', str(source), '--analysis', 'polymorphism,ld', '--output', str(out)]) == 0

    written = dict(text_writes)
    assert written.get('results.tsv') == 'utf-8'
    assert written.get('ld_pairs.tsv') == 'utf-8'
    assert all(encoding == 'utf-8' for _, encoding in text_writes), text_writes
    assert name in (out / 'results.tsv').read_text(encoding='utf-8')
