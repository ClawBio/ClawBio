"""What the gi-* client actually PUTS ON THE WIRE.

The existing gi tests are thorough about responses -- malformed envelopes,
wrong-typed fields, refusals -- and asserted nothing about the request. That
left the headline feature the worst covered path: ``--tss-index`` reaches the
API only through the two lines in ``gi_client._build_body`` that copy it into
the body, and every ``--tss-index`` test drove
``gi_runner.validate_expression_input``, a pure local function. Deleting those
two lines kept the whole suite green while ``--tss-index`` became a silent
no-op: the user submits a 50 kbp locus, the API scores some default window,
and ``report.md`` states a prediction for a TSS the model never saw.

That is the same shape as the two response-side defects already fixed on this
branch -- a client reading or writing a field name nobody checked -- so these
assert the POSTed JSON rather than the source.

Pure-local: the session is a recorder, so these need no key and no network.
"""

from __future__ import annotations

import sys

import pytest

from clawbio.gi import gi_client, gi_runner


class _Resp:
    status_code = 200
    ok = True
    text = ""

    def __init__(self, payload):
        self._payload = payload
        self.headers: dict = {}

    def json(self):
        return self._payload


class _Recorder:
    """A session that records the request instead of sending it."""

    def __init__(self, payload):
        self._payload = payload
        self.headers: dict = {}
        self.calls: list[dict] = []

    def post(self, url, json=None, headers=None, timeout=None):
        self.calls.append({"url": url, "json": json, "headers": headers})
        return _Resp(self._payload)

    def get(self, url, timeout=None):
        self.calls.append({"url": url, "json": None, "headers": None})
        return _Resp(self._payload)

    @property
    def body(self) -> dict:
        assert len(self.calls) == 1, f"expected exactly one request, got {len(self.calls)}"
        return self.calls[0]["json"]


OK_PREDICT = {"data": {"summary": {"expression_log_tpm": 0.9}}, "meta": {"request_id": "r"}}
OK_SUBMIT = {"data": {"job_id": "job-1"}, "meta": {"request_id": "r"}}


def _client(payload):
    c = gi_client.Client.__new__(gi_client.Client)
    c.api_key = "gi_test"
    c.base_url = "https://api.genomicintelligence.ai"
    c.timeout = 5.0
    c._session = _Recorder(payload)
    return c


class TestTssIndexReachesTheWire:
    """The two lines in _build_body that nothing was holding."""

    def test_predict_sends_tss_index(self):
        c = _client(OK_PREDICT)
        c.predict("expression", sequence="ACGT" * 10, tss_index=4599)
        assert c._session.body["tss_index"] == 4599

    def test_submit_async_sends_tss_index(self):
        c = _client(OK_SUBMIT)
        c.submit_async("expression", sequence="ACGT" * 10, tss_index=4599)
        assert c._session.body["tss_index"] == 4599
        assert c._session.calls[0]["headers"] == {"Prefer": "respond-async"}

    def test_a_zero_tss_index_is_still_sent(self):
        """0 is a legal offset and falsy; `if tss_index:` would drop it."""
        c = _client(OK_PREDICT)
        c.predict("expression", sequence="ACGT" * 10, tss_index=0)
        assert c._session.body["tss_index"] == 0

    def test_tss_index_is_omitted_when_unset(self):
        """Every per-task request model is additionalProperties: false, so a
        null tss_index on a non-expression task is a 422 extra_forbidden
        rather than a silent ignore. Absent must mean absent."""
        c = _client(OK_PREDICT)
        c.predict("promoter", sequence="ACGT" * 10)
        assert "tss_index" not in c._session.body

    def test_the_other_fields_are_forwarded_too(self):
        c = _client(OK_PREDICT)
        c.predict(
            "expression",
            sequence="ACGTACGT",
            sequence_name="hbb",
            model="g0-expression",
            options={"description": "K562"},
            tss_index=4599,
        )
        body = c._session.body
        assert body["sequence"] == "ACGTACGT"
        assert body["sequence_name"] == "hbb"
        assert body["model"] == "g0-expression"
        assert body["options"] == {"description": "K562"}
        assert c._session.calls[0]["url"].endswith("/v1/tasks/expression/predict")

    def test_model_and_options_are_omitted_when_unset(self):
        c = _client(OK_PREDICT)
        c.predict("promoter", sequence="ACGT" * 10)
        assert set(c._session.body) == {"sequence", "sequence_name"}


class TestTheCliOptionReachesTheWire:
    """End to end through run_skill, which is where the no-op would bite.

    The client-level tests above pass an argument directly; these prove the
    CLI flag is wired to it, so neither half can be deleted quietly.
    """

    @staticmethod
    def _locus(tmp_path, bases):
        fa = tmp_path / "in.fa"
        fa.write_text(">locus\n" + ("ACGT" * (bases // 4))[:bases] + "\n")
        return fa

    def _run(self, tmp_path, monkeypatch, argv, task, bases):
        """Drives a *real* Client with a recording session underneath it.

        Stubbing gi_runner.Client instead would intercept above _build_body,
        which is the code under test -- and it does not fail when the
        tss_index lines are deleted. Asserting the recorded body is the whole
        point, so the seam has to be the session, not the client.
        """
        fa = self._locus(tmp_path, bases)
        recorder = _Recorder({
            "data": {
                "summary": {"expression_log_tpm": 0.9, "promoter_windows": 1,
                            "total_windows": 4},
                "input": {"tss_index": 12000, "scored_window": [7401, 16599],
                          "submitted_sequence_length": bases},
            },
            "meta": {"request_id": "req-1"},
        })

        def _make_client(*a, **kw):
            c = gi_client.Client.__new__(gi_client.Client)
            c.api_key = "gi_test"
            c.base_url = "https://api.genomicintelligence.ai"
            c.timeout = 5.0
            c._session = recorder
            return c

        monkeypatch.setattr(gi_runner, "Client", _make_client)
        monkeypatch.setenv("GI_API_KEY", "gi_test")
        monkeypatch.setattr(
            sys, "argv",
            ["gi_x.py", "--input", str(fa), "--output", str(tmp_path / "out")] + argv,
        )
        code = gi_runner.run_skill(task=task, demo_path=tmp_path / "demo.fa")
        return code, recorder

    def test_the_flag_reaches_the_posted_body(self, tmp_path, monkeypatch):
        code, recorder = self._run(
            tmp_path, monkeypatch, ["--tss-index", "12000"], "expression", 30000,
        )
        assert code == 0
        assert recorder.body["tss_index"] == 12000

    def test_a_non_expression_task_never_sends_one(self, tmp_path, monkeypatch):
        """tss_index is declared on ExpressionPredictRequest alone, and every
        per-task request model is additionalProperties: false."""
        code, recorder = self._run(tmp_path, monkeypatch, [], "promoter", 4000)
        assert code == 0
        assert "tss_index" not in recorder.body

    def test_the_report_states_the_window_the_flag_asked_for(self, tmp_path, monkeypatch):
        """The user-visible end of the same wire.

        A --tss-index that never reached the API would still render a
        confident report; this ties the flag to what the report prints.
        """
        code, _ = self._run(
            tmp_path, monkeypatch, ["--tss-index", "12000"], "expression", 30000,
        )
        assert code == 0
        report = (tmp_path / "out" / "report.md").read_text()
        assert "TSS index 12000" in report
