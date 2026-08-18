#!/bin/sh
# Everything that must pass before anything is pushed.
#
# This exists because a `cmd | tail -2 && git push` chain once masked seven test
# failures behind a pipe and pushed a repo whose demo did not work. Pipes discard
# exit codes; this script does not.
#
# Usage: sh skills/abstention-ledger/verify.sh [repo_root]
set -eu

ROOT="${1:-.}"
cd "$ROOT"
FAILED=0

step() {
    printf '\n== %s\n' "$1"
}

fail() {
    printf '   FAILED: %s\n' "$1"
    FAILED=$((FAILED + 1))
}

step "unit tests"
if python3 skills/abstention-ledger/tests/test_abstention_ledger.py; then
    :
else
    fail "test suite"
fi

step "demo runs and is self-consistent"
if python3 skills/abstention-ledger/abstention_ledger.py --demo --output /tmp/al_verify >/dev/null; then
    python3 - <<'PY' || exit 1
import json
d = json.load(open("/tmp/al_verify/result.json"))
assert d["records"] == 9, d["records"]
assert d["withheld"] == 5, d["withheld"]
assert d["reviewable"] == 4, d["reviewable"]
# The secondary-findings list must have loaded by SOME path. A zero here means
# every record was withheld for want of a list, which is not a real result.
assert d["screening"]["sf_list_size"] > 0, d["screening"]
assert "not a medical device" in open("/tmp/al_verify/report.md").read()
print(f"   ok  9 records, 5 withheld, 4 reviewable, SF list {d['screening']['sf_list_size']} genes")
print(f"   ok  source: {d['screening']['sf_source'][:80]}...")
PY
else
    fail "demo run"
fi

step "prohibited-claim lint"
# Only our own prose. In a ClawBio checkout the root README.md belongs to the
# library and legitimately discusses variant classification; linting it would
# report their vocabulary as our violation. Ours is identified by its title.
LINT_TARGETS=""
for f in out/run1/report.md out/demo/report.md out/gap/report.md \
         SUBMISSION.md skills/abstention-ledger/SKILL.md; do
    [ -f "$f" ] && LINT_TARGETS="$LINT_TARGETS $f"
done
if [ -f README.md ] && head -5 README.md | grep -q 'Abstention Ledger'; then
    LINT_TARGETS="$LINT_TARGETS README.md"
fi
if [ -n "$LINT_TARGETS" ]; then
    # shellcheck disable=SC2086
    python3 skills/abstention-ledger/lint_claims.py $LINT_TARGETS || fail "claims lint"
else
    fail "no lint targets found"
fi

step "no secrets outside .env"
# The pattern is split so this script does not contain the literal it searches
# for. Without that, the check matches itself and reports a leak on every run —
# which it did, and a self-match that always fires is a check nobody trusts.
KEY_PREFIX='tvly''-'
if grep -rI "$KEY_PREFIX" --exclude-dir=.git --exclude='.env' . >/dev/null 2>&1; then
    printf '   files:\n'
    grep -rIl "$KEY_PREFIX" --exclude-dir=.git --exclude='.env' . 2>/dev/null | sed 's/^/     /'
    fail "an API key appears outside .env"
else
    printf '   ok  no key material outside .env\n'
fi

printf '\n'
if [ "$FAILED" -ne 0 ]; then
    printf '%s check group(s) FAILED — do not push.\n' "$FAILED"
    exit 1
fi
printf 'all checks passed — safe to push.\n'
