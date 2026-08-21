"""A malformed 2xx from the GI API is refused and named, never substituted.

Pure-local: every response here is a fake, so these need no key and no network.

Two defects are pinned. The first is the envelope: ``predict`` and
``wait_for_job`` returned whatever the API sent straight to the report writer,
so a 200 carrying no usable ``data`` produced a zero-valued report and an OK
line on stderr. Only ``submit_async`` checked anything, via its ``job_id``
read, and even that accepted ``{"data": {}}``.

The second is what is *inside* ``data``. The envelope check stops at "``data``
is a non-empty object", which is the right scope for a helper shared by six
tasks, so the report writer still meets whatever ``data.summary`` actually
contains — and it runs after ``run_skill``'s try/except has closed. Two wrong
answers were available there: ``x or {}`` let a truthy wrong type through to an
AttributeError traceback, and coercing it to ``{}`` instead turned that into a
zero-valued report printed as a success, which is worse because it is quiet.
The third option is a typed refusal, and that is what these assert.

These drive the client methods and ``run_skill`` rather than asserting on their
source. The last defect of this shape was a source-level edit that matched
nothing and still read as correct in review.
"""

from __future__ import annotations

import copy
import json
import sys

import pytest

from clawbio.gi import gi_client, gi_runner


class _Resp:
    """Enough of a ``requests.Response`` for the paths under test."""

    def __init__(self, payload=None, status_code=200, text=""):
        self._payload = payload
        self.status_code = status_code
        self.text = text
        self.headers: dict = {}

    @property
    def ok(self):
        return self.status_code < 400

    def json(self):
        if self._payload is _NOT_JSON:
            raise ValueError("not json")
        return self._payload


_NOT_JSON = object()


class _Session:
    """Returns a canned response to whatever the client asks for."""

    def __init__(self, resp):
        self._resp = resp
        self.headers: dict = {}

    def post(self, url, json=None, headers=None, timeout=None):
        return self._resp

    def get(self, url, timeout=None):
        return self._resp


def _client(resp):
    c = gi_client.Client.__new__(gi_client.Client)
    c.api_key = "gi_test"
    c.base_url = "https://api.genomicintelligence.ai"
    c.timeout = 5.0
    c._session = _Session(resp)
    return c


MALFORMED_ENVELOPES = [
    {},                                    # no data key at all
    {"meta": {}},                          # meta only
    {"data": None, "meta": {}},            # null data
    {"data": "summary", "meta": {}},       # non-object data
    {"data": [1, 2], "meta": {}},          # array data
    {"data": {}, "meta": {"request_id": "req-1"}},  # object, but no result in it
    [],                                    # not an object at all
    "text",
    None,
]


class TestTheEnvelopeIsCheckedOnEveryPath:
    """All three paths read content out of ``data``, so all three check it.

    ``submit_async`` was the only one guarded, and its ``job_id`` read let
    ``{"data": {}}`` through as far as the check itself. The sync and
    job-result paths had nothing.
    """

    @pytest.mark.parametrize("payload", MALFORMED_ENVELOPES)
    def test_predict_refuses_it(self, payload):
        with pytest.raises(gi_client.GIError):
            _client(_Resp(payload)).predict("promoter", "ACGT" * 100)

    @pytest.mark.parametrize("payload", MALFORMED_ENVELOPES)
    def test_a_finished_job_refuses_it(self, payload):
        with pytest.raises(gi_client.GIError):
            _client(_Resp(payload)).wait_for_job("job-1")

    @pytest.mark.parametrize("payload", MALFORMED_ENVELOPES)
    def test_submit_async_refuses_it(self, payload):
        with pytest.raises(gi_client.GIError):
            _client(_Resp(payload)).submit_async("promoter", "ACGT" * 100)

    def test_an_empty_data_never_reaches_the_job_id_read(self):
        """``{"data": {}}`` is refused by the envelope check itself.

        It used to reach ``data.job_id`` on this path only, and nothing at all
        caught it on the other two.
        """
        with pytest.raises(gi_client.GIError) as exc:
            _client(_Resp({"data": {}, "meta": {}})).submit_async("promoter", "ACGT")
        assert "empty object" in str(exc.value)

    def test_a_data_without_job_id_is_still_refused(self):
        """The job_id check stays: a non-empty ``data`` can still lack it."""
        resp = _Resp({"data": {"status": "queued"}, "meta": {}})
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).submit_async("promoter", "ACGT")
        assert "job_id" in str(exc.value)

    def test_a_non_json_200_is_not_returned_as_a_result(self):
        """``_check`` builds an ``{"error": ...}`` dict for an unparseable body.

        It has no ``data`` key, so the envelope check refuses it rather than
        handing the caller an error envelope as a successful prediction. On the
        job path the same body used to raise ValueError out of ``.json()``.
        """
        for call in (
            lambda c: c.predict("promoter", "ACGT" * 100),
            lambda c: c.wait_for_job("job-1"),
        ):
            with pytest.raises(gi_client.GIError):
                call(_client(_Resp(_NOT_JSON, text="<html>gateway</html>")))

    def test_a_well_formed_envelope_passes_through(self):
        body = {"data": {"summary": {"total_windows": 3}}, "meta": {"request_id": "r"}}
        assert _client(_Resp(body)).predict("promoter", "ACGT" * 100) == body
        assert _client(_Resp(body)).wait_for_job("job-1") == body

    def test_health_is_not_enveloped(self):
        """``/health`` is deliberately un-enveloped and must not be checked."""
        assert _client(_Resp({"status": "ok"})).health() == {"status": "ok"}


class TestNestedFieldsOfTheWrongType:
    """The envelope passing is not a promise about what is inside ``data``."""

    _MALFORMED = [
        {"data": {"summary": "all good"}},            # summary as a string
        {"data": {"summary": ["a", "b"]}},            # summary as an array
        {"data": {"summary": 0.94}},                  # summary as a float
        {"data": {"summary": {}, "prediction": "high"}},   # prediction as a string
        {"data": {"summary": {}, "input": "9198bp"}},      # input as a string
        {"data": {"summary": {}}, "meta": "req-1"},        # meta as a string
        {"data": {"summary": {}, "regions": "chr1"}},      # array field as a string
        {"data": {"summary": {}, "sites": [1, 2, 3]}},     # non-object elements
        {"data": {"summary": {}, "transcripts": "ENST1"}},
    ]

    @pytest.mark.parametrize("body", _MALFORMED)
    @pytest.mark.parametrize(
        "task", ["promoter", "splice", "enhancer", "chromatin", "expression", "annotation"]
    )
    def test_a_wrong_typed_field_is_refused_not_coerced(self, task, body):
        """Whichever task reads the offending field must refuse the body.

        A task that never reads it is entitled to succeed — ``enhancer`` does
        not touch ``data.transcripts`` — so the assertion is on the failure
        mode, not on every combination failing: either a typed refusal naming
        the field, or a clean summary. Never an AttributeError, and never a
        summary built out of a substituted empty value.
        """
        try:
            out = gi_runner._summarize(task, body)
        except gi_runner.ResponseShapeError as e:
            assert "should be an" in str(e)
            return
        assert isinstance(out, dict)
        assert isinstance(out["raw_summary"], dict)

    @pytest.mark.parametrize(
        "body,task,field",
        [
            ({"data": {"summary": "all good"}}, "promoter", "data.summary"),
            ({"data": {"summary": ["a", "b"]}}, "promoter", "data.summary"),
            ({"data": {"summary": 0.94}}, "promoter", "data.summary"),
            ({"data": {"summary": {}, "regions": "chr1"}}, "promoter", "data.regions"),
            ({"data": {"summary": {}, "sites": [1, 2, 3]}}, "splice", "data.sites"),
            ({"data": {"summary": {}, "prediction": "high"}}, "expression", "data.prediction"),
            ({"data": {"summary": {}, "input": "9198bp"}}, "expression", "data.input"),
            ({"data": {"summary": {}, "transcripts": "ENST1"}}, "annotation", "data.transcripts"),
        ],
    )
    def test_the_offending_field_is_named(self, body, task, field):
        with pytest.raises(gi_runner.ResponseShapeError) as exc:
            gi_runner._summarize(task, body)
        assert field in str(exc.value)

    @pytest.mark.parametrize("body", _MALFORMED)
    @pytest.mark.parametrize("task", ["promoter", "splice", "expression", "annotation"])
    def test_the_report_never_half_writes(self, tmp_path, task, body):
        """A refusal may happen; a traceback or a silent zero report may not."""
        try:
            summary = gi_runner._summarize(task, body)
            gi_runner._write_report(
                task, summary, body, tmp_path, tmp_path / "in.fa", "seq", 9198, 12.0,
            )
        except gi_runner.ResponseShapeError:
            return
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    def test_absent_and_null_are_still_legitimate(self):
        """Only a *present, wrong-typed* field is malformed.

        A task with no ``prediction`` omits it; treating that as ``{}`` is
        correct and must not become a refusal, or every sparse-but-valid
        response breaks.
        """
        body = {"data": {"summary": {"total_windows": 5}, "regions": None},
                "meta": None}
        out = gi_runner._summarize("promoter", body)
        assert out["regions"] == []
        assert out["raw_summary"] == {"total_windows": 5}

    def test_a_well_formed_body_still_reports_its_rows(self, tmp_path):
        body = {
            "data": {
                "summary": {"promoter_windows": 2, "total_windows": 5},
                "regions": [{"window_index": 0, "start": 10, "end": 20, "probability": 0.9}],
            },
            "meta": {"request_id": "req-1", "inference_time_ms": 42},
        }
        summary = gi_runner._summarize("promoter", body)
        assert summary["regions"] == body["data"]["regions"]
        gi_runner._write_report(
            "promoter", summary, body, tmp_path, tmp_path / "in.fa", "seq", 800, 12.0,
        )
        report = (tmp_path / "report.md").read_text()
        assert "req-1" in report and "0.900" in report


# Recorded from a live POST /v1/tasks/expression/predict on 2026-08-21 against
# PROD (contract revision 5), submitting skills/gi-expression/example_data/
# expression_hbb_k562.fa with the skill's own default options. Trimmed to the
# keys under test; the values are verbatim.
#
# The windowing pair is echoed in BOTH places, with identical values:
# data.input, which also carries submitted_sequence_length, and
# meta.task_specific_counts, which does not. gi-expression/SKILL.md named only
# the meta path while the report reads data.input, so the two disagreed on
# paper and nothing pinned either. Both are real; this pins both, so a future
# response that drops one fails here rather than in a user's report.
LIVE_EXPRESSION_BODY = {
    "data": {
        "task": "expression",
        "model": "g0-expression",
        "input": {
            "sequence_length": 9198,
            "tss_index": 4599,
            "scored_window": [0, 9198],
            "submitted_sequence_length": 9198,
        },
        "summary": {"expression_log_tpm": 0.9492, "expression_tpm": 1.5837},
        "prediction": {
            "expression_log_tpm": 0.9492,
            "expression_tpm": 1.5837,
            "unit": "log(TPM+1)",
        },
    },
    "meta": {
        "request_id": "req-live",
        "inference_time_ms": 88,
        "task_specific_counts": {
            "task": "expression",
            "tss_index": 4599,
            "scored_window": [0, 9198],
        },
    },
}


class TestTheExpressionWindowingProvenanceSurvivesToTheReport:
    """The "Scored window" line is the only defence against a wrong --tss-index.

    SKILL.md says an offset that is merely wrong "does not error, it lies", and
    points the reader at this line. Nothing asserted it end to end, so the path
    the runner reads could have gone stale silently -- which is exactly the
    failure 4077ee4 fixed for the splice table, where the report read field
    names the API had never returned and rendered as "-" for months.
    """

    def test_both_documented_paths_are_present_in_a_live_response(self):
        counts = LIVE_EXPRESSION_BODY["meta"]["task_specific_counts"]
        inp = LIVE_EXPRESSION_BODY["data"]["input"]
        assert inp["scored_window"] == counts["scored_window"] == [0, 9198]
        assert inp["tss_index"] == counts["tss_index"] == 4599
        # Only data.input carries this, which is why the runner reads there.
        assert "submitted_sequence_length" not in counts

    def test_summarize_picks_up_the_window(self):
        out = gi_runner._summarize("expression", LIVE_EXPRESSION_BODY)
        assert out["tss_index"] == 4599
        assert out["scored_window"] == [0, 9198]
        assert out["submitted_sequence_length"] == 9198
        assert out["log_tpm"] == 0.9492

    def test_the_report_states_the_scored_window(self, tmp_path):
        summary = gi_runner._summarize("expression", LIVE_EXPRESSION_BODY)
        gi_runner._write_report(
            "expression", summary, LIVE_EXPRESSION_BODY, tmp_path,
            tmp_path / "in.fa", "hbb", 9198, 12.0,
        )
        report = (tmp_path / "report.md").read_text()
        assert "Scored window" in report
        assert "[0, 9198)" in report
        assert "TSS index 4599" in report

    def test_a_wrong_typed_window_is_refused_not_indexed(self, tmp_path):
        """A dict here used to raise KeyError(0) out of the report writer.

        ``window[0]`` on ``{"start": ..., "end": ...}`` raises, and it raises
        inside the one block whose stated purpose is refusing wrong-typed
        response fields -- so the runner's ResponseShapeError diagnostic was
        bypassed by an uncaught KeyError.
        """
        body = copy.deepcopy(LIVE_EXPRESSION_BODY)
        body["data"]["input"]["scored_window"] = {"start": 0, "end": 9198}
        summary = gi_runner._summarize("expression", body)
        with pytest.raises(gi_runner.ResponseShapeError, match="scored_window"):
            gi_runner._write_report(
                "expression", summary, body, tmp_path,
                tmp_path / "in.fa", "hbb", 9198, 12.0,
            )


class TestTheRunnerExitsTwoRatherThanReportingNothing:
    """The whole point: a diagnostic and exit 2, not a zero-valued report.

    Drives ``run_skill`` end to end rather than grepping it for a handler.
    """

    @staticmethod
    def _fasta(tmp_path):
        fa = tmp_path / "in.fa"
        fa.write_text(">seq\n" + "ACGT" * 100 + "\n")  # 400 bp, above the promoter floor
        return fa

    def _run(self, tmp_path, monkeypatch, payload):
        fa = self._fasta(tmp_path)

        class _FakeClient:
            def __init__(self, *a, **kw):
                pass

            def predict(self, *a, **kw):
                return payload

        monkeypatch.setattr(gi_runner, "Client", _FakeClient)
        monkeypatch.setenv("GI_API_KEY", "gi_test")
        monkeypatch.setattr(
            sys, "argv",
            ["gi_promoter.py", "--input", str(fa), "--output", str(tmp_path / "out")],
        )
        return gi_runner.run_skill(task="promoter", demo_path=tmp_path / "demo.fa")

    def test_a_wrong_typed_summary_exits_two_and_names_the_field(
        self, tmp_path, monkeypatch, capsys
    ):
        code = self._run(tmp_path, monkeypatch, {"data": {"summary": "all good"}, "meta": {}})
        err = capsys.readouterr().err
        assert code == 2
        assert "unexpected API response shape" in err
        assert "data.summary" in err
        assert "OK — wrote" not in err
        assert not (tmp_path / "out" / "report.md").exists()

    def test_a_well_formed_response_still_writes_its_report(self, tmp_path, monkeypatch):
        code = self._run(
            tmp_path, monkeypatch,
            {"data": {"summary": {"promoter_windows": 1, "total_windows": 4}},
             "meta": {"request_id": "req-1"}},
        )
        assert code == 0
        assert (tmp_path / "out" / "report.md").exists()


class TestGIErrorsOwnFields:
    """``GIError``'s two load-bearing attributes, which nothing was holding.

    ``code`` is the one the class docstring tells callers to branch on, and it
    was renamed from ``non_json`` to ``http_error`` on this branch -- a
    breaking change for anyone who had branched on the old value, and the
    suite would not have noticed either the rename or a revert of it. The
    ``X-Request-Id`` fallback exists for exactly the non-JSON case, where the
    envelope that normally carries ``request_id`` is not there to read, and
    that is the case a support ticket needs a correlation id most.
    """

    def test_a_non_json_body_is_reported_as_http_error(self):
        resp = _Resp(_NOT_JSON, status_code=502, text="<html>gateway</html>")
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.code == "http_error"
        assert exc.value.status == 502
        assert "gateway" in exc.value.message

    def test_an_error_body_that_is_not_our_envelope_still_gets_a_code(self):
        """A different line from the one above.

        ``_check`` writes ``http_error`` into a synthetic envelope when the
        body will not parse; ``GIError.__init__`` defaults to it when the body
        parsed but carries no ``error`` key -- a proxy or framework error page
        is often perfectly good JSON. Both spell the same published enum value
        and each needs its own test, or one can be renamed without the other.
        """
        resp = _Resp({"detail": "Not Found"}, status_code=404)
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.code == "http_error"

    def test_a_code_the_api_sent_is_carried_through_unchanged(self):
        """``http_error`` is the default, not an override."""
        resp = _Resp(
            {"error": {"code": "validation_failed", "message": "sequence is 4 bp"}},
            status_code=422,
        )
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT")
        assert exc.value.code == "validation_failed"

    def test_the_request_id_falls_back_to_the_header(self):
        resp = _Resp(_NOT_JSON, status_code=503, text="upstream unavailable")
        resp.headers = {"X-Request-Id": "req-from-header"}
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.request_id == "req-from-header"
        assert "req-from-header" in str(exc.value)

    def test_the_envelope_request_id_wins_over_the_header(self):
        """Every error envelope carries one; the header is only the fallback."""
        resp = _Resp(
            {"error": {"code": "rate_limited", "message": "slow down",
                       "request_id": "req-from-body"}},
            status_code=429,
        )
        resp.headers = {"X-Request-Id": "req-from-header"}
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.request_id == "req-from-body"

    def test_no_request_id_anywhere_still_produces_a_readable_message(self):
        resp = _Resp(_NOT_JSON, status_code=500, text="boom")
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.request_id is None
        assert "request_id=unset" in str(exc.value)


class TestTheRateLimitHeadersAreNotDiscarded:
    """Six SKILL.md files tell the reader to consult these; nothing kept them.

    The caps are per-key and are retuned server-side rather than published as
    a fixed number, so the headers are the only statement of the live
    allowance. ``GIError`` captured ``X-Request-Id`` and nothing else, and
    ``result.json`` holds the parsed body, which does not carry headers -- so
    a user following the documentation had nowhere to look.

    Header names here are the lowercase wire form measured against PROD on
    2026-08-21, not the canonical capitalisation, because that is what a
    caller actually receives.
    """

    LIVE_HEADERS = {
        "ratelimit-limit": "200",
        "ratelimit-policy": "200;w=60",
        "ratelimit-remaining": "199",
        "ratelimit-reset": "1",
        "x-request-id": "d6bc8156-0eb4-4451-843d-84ba7c6eebcf",
    }

    def test_the_lowercase_wire_form_is_normalised(self):
        out = gi_client.rate_limit_from(self.LIVE_HEADERS)
        assert out == {
            "RateLimit-Limit": "200",
            "RateLimit-Remaining": "199",
            "RateLimit-Reset": "1",
            "RateLimit-Policy": "200;w=60",
        }

    def test_a_header_the_response_did_not_send_is_omitted_not_nulled(self):
        """``/health`` sends none of these, and no 2xx observed sends
        ``Retry-After``. An invented null would read as a measured zero."""
        assert "Retry-After" not in gi_client.rate_limit_from(self.LIVE_HEADERS)
        assert gi_client.rate_limit_from({}) == {}
        assert gi_client.rate_limit_from(None) == {}

    def test_a_429_states_the_allowance_in_its_message(self):
        """The one status where the numbers are the message rather than noise.

        gi-promoter/SKILL.md says `Retry-After` on a 429 is the wait; before
        this the runner printed the error line with neither.
        """
        resp = _Resp({"error": {"code": "rate_limited", "message": "slow down"}},
                     status_code=429)
        resp.headers = dict(self.LIVE_HEADERS, **{"retry-after": "30", "ratelimit-remaining": "0"})
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT" * 100)
        assert exc.value.rate_limit["RateLimit-Remaining"] == "0"
        assert exc.value.rate_limit["Retry-After"] == "30"
        assert "Retry-After: 30" in str(exc.value)

    def test_a_422_carries_them_without_burying_its_own_message(self):
        """A validation error's message says what to fix; the numbers do not."""
        resp = _Resp(
            {"error": {"code": "validation_failed", "message": "sequence is 4 bp"}},
            status_code=422,
        )
        resp.headers = dict(self.LIVE_HEADERS)
        with pytest.raises(gi_client.GIError) as exc:
            _client(resp).predict("promoter", "ACGT")
        assert exc.value.rate_limit["RateLimit-Limit"] == "200"
        assert "RateLimit-Limit" not in str(exc.value)
        assert "sequence is 4 bp" in str(exc.value)

    def test_a_successful_call_leaves_the_allowance_on_the_client(self):
        body = {"data": {"summary": {"total_windows": 3}}, "meta": {"request_id": "r"}}
        resp = _Resp(body)
        resp.headers = dict(self.LIVE_HEADERS)
        c = _client(resp)
        c.predict("promoter", "ACGT" * 100)
        assert c.last_rate_limit["RateLimit-Remaining"] == "199"

    def test_health_sends_none_and_that_is_recorded_as_none(self):
        """Measured: /health carries no rate-limit headers at all."""
        c = _client(_Resp({"status": "ok"}))
        c.health()
        assert c.last_rate_limit == {}

    def test_result_json_records_them_end_to_end(self, tmp_path, monkeypatch):
        """Through run_skill and a real Client, which is where it has to work.

        Stubbing gi_runner.Client would intercept above the client attribute
        this reads, so the seam is the session.
        """
        resp = _Resp(
            {"data": {"summary": {"promoter_windows": 1, "total_windows": 4}},
             "meta": {"request_id": "req-1"}},
        )
        resp.headers = dict(self.LIVE_HEADERS)

        def _make_client(*a, **kw):
            return _client(resp)

        fa = tmp_path / "in.fa"
        fa.write_text(">seq\n" + "ACGT" * 100 + "\n")
        monkeypatch.setattr(gi_runner, "Client", _make_client)
        monkeypatch.setenv("GI_API_KEY", "gi_test")
        monkeypatch.setattr(
            sys, "argv",
            ["gi_promoter.py", "--input", str(fa), "--output", str(tmp_path / "out")],
        )
        assert gi_runner.run_skill(task="promoter", demo_path=tmp_path / "demo.fa") == 0
        result = json.loads((tmp_path / "out" / "result.json").read_text())
        assert result["rate_limit"] == {
            "RateLimit-Limit": "200",
            "RateLimit-Remaining": "199",
            "RateLimit-Reset": "1",
            "RateLimit-Policy": "200;w=60",
        }
        assert "RateLimit-*" in (tmp_path / "out" / "report.md").read_text()


class TestTheSpliceTableReadsTheFieldNamesTheApiSends:
    """Every shipped splice report rendered three of its four columns as "-".

    The renderer read ``position`` / ``kind`` / ``probability``; the API
    returns ``name`` / ``start`` / ``end`` / ``site_type`` / ``score``, and
    only ``strand`` overlapped. A missing key renders as a dash rather than
    raising, so nothing failed -- and nothing tested it, which means reverting
    the fix would leave the suite green and the reports empty again.

    The rendered row is the only place the field names are observable, so
    that is what these assert.
    """

    BODY = {
        "data": {
            "summary": {"total_sites": 2, "donor_sites": 1, "acceptor_sites": 1},
            "sites": [
                {"name": "site_1", "start": 2324, "end": 2330,
                 "site_type": "donor", "strand": "+", "score": 0.9871},
                {"name": "site_2", "start": 4110, "end": 4117,
                 "site_type": "acceptor", "strand": "+", "score": 0.4432},
            ],
        },
        "meta": {"request_id": "req-splice"},
    }

    def _report(self, tmp_path, body):
        summary = gi_runner._summarize("splice", body)
        gi_runner._write_report(
            "splice", summary, body, tmp_path, tmp_path / "in.fa", "hbb", 5000, 12.0,
        )
        return (tmp_path / "report.md").read_text()

    def _row(self, report, needle):
        rows = [l for l in report.splitlines() if needle in l]
        assert len(rows) == 1, f"expected one row containing {needle!r}, got {len(rows)}"
        return rows[0]

    def test_every_column_carries_its_value(self, tmp_path):
        row = self._row(self._report(tmp_path, self.BODY), "site_1")
        assert "2324-2330" in row      # start/end, not the old `position`
        assert "donor" in row          # site_type, not the old `kind`
        assert "0.9871" in row         # score, not the old `probability`
        assert "| + |" in row          # strand, the one name that never changed

    def test_the_span_is_a_span_and_not_one_endpoint(self, tmp_path):
        """start/end bound one variable-width token with the junction inside it
        (GI-054), so printing a single endpoint states a boundary the response
        does not give."""
        row = self._row(self._report(tmp_path, self.BODY), "site_2")
        assert "4110-4117" in row

    def test_the_pre_fix_field_names_would_render_nothing(self, tmp_path):
        """The shape every shipped report actually saw, pinned as the failure.

        Without this, a revert of the renamed columns is invisible: the table
        still renders, just with no data in it.
        """
        body = copy.deepcopy(self.BODY)
        body["data"]["sites"] = [
            {"position": 2324, "kind": "donor", "strand": "+", "probability": 0.9871}
        ]
        report = self._report(tmp_path, body)
        assert "| - | - | - | + | - |" in report
