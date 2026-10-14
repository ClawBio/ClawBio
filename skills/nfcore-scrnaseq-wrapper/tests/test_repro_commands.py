import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

import repro_commands


REMOTE_SOURCE = {
    "source_kind": "remote_repo",
    "source_ref": "nf-core/scrnaseq",
    "resolved_version": "4.1.0",
}


def _build(**overrides):
    kwargs = dict(
        pipeline_source=REMOTE_SOURCE,
        profile="docker",
        demo=False,
        macos_docker_config=False,
        nextflow_version=None,
        generated_at="2026-05-31 00:00 UTC",
    )
    kwargs.update(overrides)
    return repro_commands.build_nextflow_commands_sh(**kwargs)


def test_remote_replay_is_nextflow_direct_and_self_anchoring():
    script = _build()
    assert 'BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"' in script
    assert 'OUTPUT_DIR="$(cd "$BUNDLE_DIR/.." && pwd)"' in script
    assert 'cd "$OUTPUT_DIR"' in script
    assert "nextflow run nf-core/scrnaseq -r 4.1.0" in script
    assert "-profile docker" in script
    assert "-params-file reproducibility/params.yaml" in script
    assert "-work-dir upstream/work" in script
    assert "nfcore_scrnaseq_wrapper.py" not in script
    assert "skills/" not in script


def test_remote_replay_preserves_object_store_work_dir():
    script = _build(work_dir="s3://bucket/scrnaseq/work")
    assert "-work-dir s3://bucket/scrnaseq/work" in script
    assert "-work-dir upstream/work" not in script


def test_demo_uses_test_profile_prefix():
    script = _build(demo=True)
    assert "-profile test,docker" in script


def test_resource_limits_guard_emitted_when_shipped():
    # A non-macOS docker real run ships reproducibility/resource_limits.config; the
    # replay must apply it on non-macOS hosts so a from-scratch reproduction on the
    # generating (memory-constrained) machine does not re-abort STAR_GENOMEGENERATE.
    script = _build(resource_limits_config=True)
    assert 'docker info --format' in script
    assert 'RESOURCE_CONFIG="-c reproducibility/resource_limits.config"' in script
    assert "$RESOURCE_CONFIG" in script


def test_resource_limits_guard_absent_by_default():
    script = _build()
    assert "resource_limits.config" not in script
    assert "RESOURCE_CONFIG" not in script


def test_idempotent_resume_guard_always_emitted():
    # Replay must be idempotent: reuse the previous run's cache when replaying against
    # a completed output dir, but start fresh on a first run. The guard adds -resume
    # only when a prior Nextflow session (.nextflow/) exists in the launch dir.
    script = _build()
    assert 'if [[ -d ".nextflow" ]]; then' in script
    assert 'RESUME="-resume"' in script
    assert "$RESUME" in script


def test_local_checkout_emits_portability_warning():
    local_source = {
        "source_kind": "local_checkout",
        "source_ref": "/abs/path/to/scrnaseq",
        "resolved_version": "4.1.0",
    }
    script = _build(pipeline_source=local_source)
    assert "NOT portable" in script
    assert "/abs/path/to/scrnaseq" in script


# ── macOS config is gated on the REPLAY host, not the generation host ─────────


def test_macos_guard_is_bash32_safe_uses_assigned_default():
    # macOS ships bash 3.2; an empty array under `set -u` errors. The generator
    # must use an always-assigned plain string variable, not an array expansion.
    script = _build(macos_docker_config=True)
    assert 'EXTRA_CONFIG=""' in script  # default empty assignment before the guard
    assert "EXTRA_CONFIG[@]" not in script


# ── Nextflow engine version pin (reproducibility) ─────────────────────────────


def test_nextflow_version_is_pinned_when_provided():
    script = _build(nextflow_version="24.10.1")
    assert 'export NXF_VER="24.10.1"' in script
    # The pin must precede the run command so Nextflow honours it.
    assert script.index('export NXF_VER="24.10.1"') < script.index("nextflow run")
