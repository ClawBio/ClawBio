"""Permanent security regression tests for the RoboTerri bots.

These encode the two attacks from the security audit as tests that must always
pass: (1) a forged webhook sender must be rejected, and (2) one user's genomic
data must never be reachable from another user's session (identity isolation:
exactly one authenticated identity, no global/first/any fallback).

The functions under test live in bot/security.py (stdlib-only, no flask/openai),
so this runs without the bots' heavy runtime dependencies.
"""

import hashlib
import hmac
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bot/

from security import (  # noqa: E402
    is_sender_allowed,
    scoped_get,
    verify_whatsapp_signature,
)

SECRET = "test_app_secret"


def _sign(body: bytes, secret: str = SECRET) -> str:
    return "sha256=" + hmac.new(secret.encode(), body, hashlib.sha256).hexdigest()


# --- Attack 1: a forged webhook sender must be rejected ---

def test_forged_webhook_signature_rejected():
    body = b'{"from":"attacker","text":{"body":"run pharmgx on the owner genome"}}'
    assert verify_whatsapp_signature(body, "sha256=deadbeef", SECRET) is False
    assert verify_whatsapp_signature(body, None, SECRET) is False
    assert verify_whatsapp_signature(body, "garbage", SECRET) is False
    assert verify_whatsapp_signature(body, _sign(body, "wrong_secret"), SECRET) is False


def test_valid_webhook_signature_accepted():
    body = b'{"from":"+441234567","text":{"body":"hi"}}'
    assert verify_whatsapp_signature(body, _sign(body), SECRET) is True


def test_signature_fails_closed_without_secret():
    # If no app secret is configured, every payload is rejected (fail closed),
    # never silently accepted.
    body = b"{}"
    assert verify_whatsapp_signature(body, _sign(body), "") is False


# --- Attack 2: cross-user genome/data bleed must be impossible ---

def test_no_cross_user_file_bleed():
    store = {"+userA": {"path": "/tmp/A_genome.txt"}}
    # User B requests a skill run; must NOT receive user A's uploaded genome.
    assert scoped_get(store, "+userB") is None
    # User A still reaches their own file.
    assert scoped_get(store, "+userA") == {"path": "/tmp/A_genome.txt"}


def test_scoped_get_has_no_first_or_any_fallback():
    store = {"+userA": {"path": "/tmp/A.txt"}}
    # The patterns the audit flagged (first item / any item / next(iter)) must be
    # gone: an empty or missing identity yields nothing, never an arbitrary entry.
    assert scoped_get(store, "") is None
    assert scoped_get(store, None) is None


# --- Sender allow-list: deny by default ---

def test_sender_allowlist_denies_by_default():
    assert is_sender_allowed("+stranger", allowed=set(), admin="+owner") is False
    assert is_sender_allowed("+owner", allowed=set(), admin="+owner") is True
    assert is_sender_allowed("+friend", allowed={"+friend"}, admin="+owner") is True
    assert is_sender_allowed("", allowed={"+friend"}, admin="+owner") is False


def test_sender_allowlist_explicit_public_optin():
    # Public demo mode must be an explicit choice, never the silent default.
    assert is_sender_allowed("+anyone", allowed=set(), admin=None, allow_all=True) is True
    assert is_sender_allowed("", allowed=set(), admin=None, allow_all=True) is False


# --------------------------------------------------------------------------- #
# Filesystem containment for the chat-facing save_file / write_file tools.
#
# Writes stay inside the data directory, keep an allowed extension when the
# model renames an upload, and never touch a protected name. One
# implementation, one set of tests, all three adapters.
# --------------------------------------------------------------------------- #

import pytest  # noqa: E402

from security import (  # noqa: E402
    UnsafePath,
    is_allowed_extension,
    safe_write_path,
    sanitize_filename,
)


def test_write_stays_inside_root_when_folder_points_at_source_tree(tmp_path):
    root = tmp_path / "data"
    outside = tmp_path / "skills" / "nutrigx"
    outside.mkdir(parents=True)

    path = safe_write_path("../skills/nutrigx", "report.txt", root=root)

    assert path.resolve().is_relative_to(root.resolve())


def test_absolute_destination_outside_root_is_refused(tmp_path):
    root = tmp_path / "data"

    path = safe_write_path(str(tmp_path / "elsewhere"), "report.txt", root=root)

    assert path.resolve().is_relative_to(root.resolve())


@pytest.mark.parametrize("name", [".env", "SOUL.md", "claude.md", "CLAWBIO.PY"])
def test_protected_names_are_refused(tmp_path, name):
    with pytest.raises(UnsafePath):
        safe_write_path(".", name, root=tmp_path)


def test_upload_cannot_be_renamed_to_a_python_file(tmp_path):
    with pytest.raises(UnsafePath):
        safe_write_path(".", "nutrigx.py", root=tmp_path, require_allowed_extension=True)


def test_data_file_is_accepted(tmp_path):
    path = safe_write_path(".", "genome.vcf.gz", root=tmp_path, require_allowed_extension=True)

    assert path.name == "genome.vcf.gz"
    assert path.parent.resolve() == tmp_path.resolve()


def test_write_file_may_use_any_extension_by_default(tmp_path):
    path = safe_write_path(".", "notes.md", root=tmp_path)

    assert path.name == "notes.md"


def test_traversal_in_filename_is_flattened(tmp_path):
    path = safe_write_path(".", "../../etc/passwd", root=tmp_path)

    assert path.resolve().is_relative_to(tmp_path.resolve())
    assert "/" not in path.name


def test_bare_gz_is_not_an_allowed_extension():
    assert is_allowed_extension("genome.vcf.gz") is True
    assert is_allowed_extension("payload.gz") is False
    assert is_allowed_extension("payload.py") is False


def test_sanitize_filename_never_returns_empty():
    assert sanitize_filename("../..") == "unnamed_file"


def test_adapters_do_not_define_their_own_copies():
    """One implementation: adapters import the shared helpers, never copy them."""
    import ast

    shared = {
        "_resolve_dest",
        "_validate_path",
        "_sanitize_filename",
        "_is_allowed_extension",
    }
    bot_dir = Path(__file__).resolve().parents[1]
    offenders = []
    for adapter in sorted(bot_dir.glob("roboterri*.py")):
        tree = ast.parse(adapter.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in shared:
                offenders.append(f"{adapter.name}:{node.lineno} {node.name}")
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name) and t.id == "_PROTECTED_NAMES":
                        offenders.append(f"{adapter.name}:{node.lineno} _PROTECTED_NAMES")

    assert offenders == [], "adapters must import these from security.py, not redefine them: " + "; ".join(offenders)


# --- Uploads are stored per sender, so two senders cannot collide ---

from hypothesis import given  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

from security import upload_tmp_path  # noqa: E402


def test_same_filename_from_two_senders_does_not_collide():
    assert upload_tmp_path("+447700900001", "genome.txt") != upload_tmp_path("+447700900002", "genome.txt")


@given(owner=st.text(), filename=st.text())
def test_upload_tmp_path_stays_in_tempdir(owner, filename):
    import tempfile

    path = upload_tmp_path(owner, filename)
    assert path.parent == Path(tempfile.gettempdir())
    assert path.name.startswith("roboterri_")


def test_owner_id_cannot_inject_a_path():
    assert upload_tmp_path("../../etc", "x.vcf").parent == upload_tmp_path("a", "x.vcf").parent


def test_adapters_build_upload_paths_through_the_shared_helper():
    bot_dir = Path(__file__).resolve().parents[1]
    offenders = [
        f"{adapter.name}:{n}"
        for adapter in sorted(bot_dir.glob("roboterri*.py"))
        for n, line in enumerate(adapter.read_text(encoding="utf-8").splitlines(), 1)
        if 'f"roboterri_' in line
    ]
    assert offenders == [], "use security.upload_tmp_path: " + "; ".join(offenders)
