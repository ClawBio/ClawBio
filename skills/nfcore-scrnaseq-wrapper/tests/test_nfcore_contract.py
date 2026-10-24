from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from nfcore_4_1_0_contract import (
    NFCORE_SCRNASEQ_VERSION,
    OFFICIAL_PARAMS,
    WRAPPER_SUPPORTED_UPSTREAM_PARAMS,
    INTENTIONALLY_UNSUPPORTED_PARAMS,
    DEPRECATED_PARAMS,
    WRAPPER_DEPRECATED_ALIAS_PARAMS,
    PROTOCOLS_JSON_4_1_0,
    PRESETS_SUPPORTING_CUSTOM_PROTOCOL,
    PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL,
    KNOWN_PROTOCOL_TOKENS,
)


def test_contract_version_is_4_1_0():
    assert NFCORE_SCRNASEQ_VERSION == "4.1.0"


def test_official_params_include_known_4_1_0_parameters():
    assert OFFICIAL_PARAMS["aligner"]["default"] == "simpleaf"
    assert OFFICIAL_PARAMS["protocol"]["default"] == "auto"
    # nf-core/scrnaseq 4.1.0 declares no default for save_align_intermeds (it is an
    # opt-in boolean, false when unset — same as save_reference). The contract must
    # not invent a `default: True`, which would misrepresent the upstream schema.
    assert "default" not in OFFICIAL_PARAMS["save_align_intermeds"]
    assert OFFICIAL_PARAMS["save_align_intermeds"]["type"] == "boolean"
    assert OFFICIAL_PARAMS["skip_emptydrops"]["deprecated"] is True
    assert "cellrangermulti" in OFFICIAL_PARAMS["aligner"]["enum"]
    assert (
        OFFICIAL_PARAMS["email"]["pattern"]
        == r"^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$"
    )
    assert (
        OFFICIAL_PARAMS["email_on_fail"]["pattern"]
        == OFFICIAL_PARAMS["email"]["pattern"]
    )
    assert (
        OFFICIAL_PARAMS["max_multiqc_email_size"]["pattern"]
        == r"^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$"
    )


def test_every_official_param_is_supported_or_intentionally_unsupported():
    overlap = WRAPPER_SUPPORTED_UPSTREAM_PARAMS & INTENTIONALLY_UNSUPPORTED_PARAMS
    classified = WRAPPER_SUPPORTED_UPSTREAM_PARAMS | INTENTIONALLY_UNSUPPORTED_PARAMS

    assert overlap == set()
    assert classified == set(OFFICIAL_PARAMS)


def test_deprecated_params_are_not_primary_wrapper_surface():
    assert "skip_emptydrops" in DEPRECATED_PARAMS
    assert "skip_emptydrops" in WRAPPER_DEPRECATED_ALIAS_PARAMS
    assert "skip_emptydrops" in WRAPPER_SUPPORTED_UPSTREAM_PARAMS
    assert "skip_emptydrops" not in INTENTIONALLY_UNSUPPORTED_PARAMS


def test_skill_docs_list_intentionally_unsupported_params():
    # Anchor on this test file (not the process cwd) so the check is robust no
    # matter where pytest is invoked from, matching every other test in the suite.
    skill_md = Path(__file__).resolve().parent.parent / "SKILL.md"
    content = skill_md.read_text(encoding="utf-8")
    for name in sorted(INTENTIONALLY_UNSUPPORTED_PARAMS):
        assert f"`{name}`" in content


# ── Protocol matrix: faithful mirror of assets/protocols.json @ 4.1.0 ─────────
# Source of truth:
# https://github.com/nf-core/scrnaseq/blob/4.1.0/assets/protocols.json
# These tests freeze the upstream truth so the preflight protocol rules can never
# silently drift from the pipeline they wrap.


def test_protocol_matrix_mirrors_upstream_protocols_json():
    # Tokens are normalised (lowercase, separators stripped) for matching against
    # preflight._normalize_protocol_token.
    assert PROTOCOLS_JSON_4_1_0 == {
        "standard": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq"),
        "star": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
        "kallisto": ("10xv1", "10xv2", "10xv3", "10xv4", "dropseq", "smartseq"),
        "cellranger": ("auto", "10xv1", "10xv2", "10xv3", "10xv4"),
        "cellrangerarc": ("auto",),
    }


def test_only_star_and_kallisto_define_smartseq_protocol():
    assert PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL == {"star", "kallisto"}


def test_known_protocol_tokens_are_union_of_matrix():
    assert KNOWN_PROTOCOL_TOKENS == {
        "auto",
        "10xv1",
        "10xv2",
        "10xv3",
        "10xv4",
        "dropseq",
        "smartseq",
    }


def test_presets_supporting_custom_protocol_are_the_non_cellranger_map_presets():
    # nf-core/scrnaseq 4.1.0 usage docs document custom/unknown protocol strings
    # only for simpleaf/star/kallisto. Centralised in the contract (derived from
    # the routing matrix minus the Cell Ranger family) so preflight cannot drift
    # from CELLRANGER_FAMILY_PRESETS (audit F-8).
    assert PRESETS_SUPPORTING_CUSTOM_PROTOCOL == {"standard", "star", "kallisto"}


def test_cellrangermulti_is_protocol_agnostic():
    # cellrangermulti is driven by the multi samplesheet, not --protocol, so it
    # must NOT appear in the routing matrix (the wrapper defers it to upstream).
    assert "cellrangermulti" not in PROTOCOLS_JSON_4_1_0
