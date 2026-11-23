from __future__ import annotations

import importlib.util
import json
from io import BytesIO
from pathlib import Path
from types import ModuleType
from urllib.error import HTTPError

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "capture_ncbi_fixtures.py"
EMAIL = "tester@example.org"
ESEARCH = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">\n'
    b"<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart>"
    b"<IdList><Id>NC_001604.1</Id><Id>V01146.1</Id></IdList>"
    b"<TranslationSet/><QueryTranslation>x</QueryTranslation></eSearchResult>\n"
)


@pytest.fixture
def capture(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    spec = importlib.util.spec_from_file_location("capture_ncbi_fixtures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_missing_email_exits_without_calling_ncbi(capture, monkeypatch, tmp_path) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    called = []
    monkeypatch.setattr(capture.Entrez, "esummary", lambda **kw: called.append(kw))
    assert capture.main(["--out", str(tmp_path)]) == 2
    assert called == []


def test_records_replies_errors_and_manifest_without_email(
    capture, monkeypatch, tmp_path, capsys
) -> None:
    monkeypatch.setenv("NCBI_EMAIL", EMAIL)
    requests: list[dict] = []

    def fake_esummary(**params):
        requests.append(params)
        if params["id"] == "NC_999999999":
            raise HTTPError("url", 400, "Bad Request", {}, BytesIO(b"<ERROR>bad id</ERROR>"))
        return BytesIO(b"<eSummaryResult/>")

    def fake_esearch(**params):
        return BytesIO(ESEARCH)

    monkeypatch.setattr(capture.Entrez, "esummary", fake_esummary)
    monkeypatch.setattr(capture.Entrez, "esearch", fake_esearch)

    assert capture.main(["--out", str(tmp_path)]) == 0

    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["esummary_missing.xml"]["status"] == 400
    assert (tmp_path / "esummary_missing.xml").read_bytes() == b"<ERROR>bad id</ERROR>"
    assert manifest["esummary_lambda.xml"]["status"] == 200
    assert requests[-1] == {"db": "nuccore", "id": "NC_001604.1,V01146.1"}
    assert "esummary_search.xml" in manifest
    written = "".join(p.read_text(errors="replace") for p in tmp_path.iterdir())
    out = capsys.readouterr()
    assert EMAIL not in written + out.out + out.err
