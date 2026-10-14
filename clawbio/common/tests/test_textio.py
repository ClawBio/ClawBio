"""Tests for clawbio.common.textio — cross-OS LF-only bundle writes."""

import os
from pathlib import Path

from clawbio.common.textio import write_text_lf, write_text_lf_atomic


def test_writes_lf_for_plain_content(tmp_path):
    p = tmp_path / "f.txt"
    write_text_lf(p, "a\nb\nc\n")
    assert p.read_bytes() == b"a\nb\nc\n"


def test_normalises_crlf_to_lf(tmp_path):
    p = tmp_path / "f.txt"
    write_text_lf(p, "a\r\nb\r\nc\r\n")
    raw = p.read_bytes()
    assert b"\r" not in raw
    assert raw == b"a\nb\nc\n"


def test_normalises_lone_cr_to_lf(tmp_path):
    p = tmp_path / "f.txt"
    write_text_lf(p, "a\rb\r")
    assert p.read_bytes() == b"a\nb\n"


def test_preserves_utf8(tmp_path):
    p = tmp_path / "f.txt"
    write_text_lf(p, "café — Ångström\n")
    assert p.read_bytes() == "café — Ångström\n".encode("utf-8")
    assert b"\r" not in p.read_bytes()


def test_returns_path(tmp_path):
    p = tmp_path / "f.txt"
    assert write_text_lf(p, "x\n") == p


def test_accepts_str_path(tmp_path):
    p = tmp_path / "f.txt"
    write_text_lf(str(p), "x\n")
    assert p.read_bytes() == b"x\n"


def test_atomic_writes_lf_and_normalises(tmp_path):
    from clawbio.common.textio import write_text_lf_atomic

    p = tmp_path / "f.txt"
    write_text_lf_atomic(p, "a\r\nb\rc\n")
    assert p.read_bytes() == b"a\nb\nc\n"


def test_atomic_replace_preserves_open_reader(tmp_path):
    """os.replace swaps in a new inode; a reader that already opened the file keeps
    reading the ORIGINAL content (bash executing commands.sh must not be corrupted when
    the wrapper regenerates it mid-run)."""
    from clawbio.common.textio import write_text_lf_atomic

    p = tmp_path / "script.sh"
    write_text_lf_atomic(p, "ORIGINAL-CONTENT\n")
    with p.open("r", encoding="utf-8") as held:
        first = held.read(4)
        write_text_lf_atomic(p, "REPLACED-CONTENT-DIFFERENT\n")
        rest = held.read()
    assert first + rest == "ORIGINAL-CONTENT\n"
    assert p.read_text(encoding="utf-8") == "REPLACED-CONTENT-DIFFERENT\n"


def test_atomic_preserves_executable_bit(tmp_path):
    import os
    from clawbio.common.textio import write_text_lf_atomic

    p = tmp_path / "script.sh"
    p.write_text("#!/bin/sh\n")
    p.chmod(0o755)
    write_text_lf_atomic(p, "#!/bin/sh\necho hi\n")
    assert p.stat().st_mode & 0o111, "executable bit must survive the atomic replace"


def test_atomic_new_file_honours_umask_not_mkstemp_0600(tmp_path):
    """mkstemp creates 0600. A new file must land on the umask mode a plain
    write would have given it, otherwise callers that OR in +x get 0711."""
    old = os.umask(0o022)
    try:
        p = write_text_lf_atomic(tmp_path / "fresh.txt", "hello\n")
        assert p.stat().st_mode & 0o777 == 0o644
    finally:
        os.umask(old)


def test_atomic_new_file_respects_a_restrictive_umask(tmp_path):
    old = os.umask(0o077)
    try:
        p = write_text_lf_atomic(tmp_path / "private.txt", "secret\n")
        assert p.stat().st_mode & 0o777 == 0o600
    finally:
        os.umask(old)
