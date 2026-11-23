from __future__ import annotations

import importlib.util
from io import StringIO
from pathlib import Path
from types import ModuleType

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "fetch_demo.py"
EMAIL = "tester@example.org"


@pytest.fixture
def fd(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    spec = importlib.util.spec_from_file_location("fetch_demo", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "DEMOS", (module.Demo("mini", "NC_000001.1", 1_200),))
    return module


def genbank_text(accession: str, length: int) -> str:
    rec = SeqRecord(Seq("ACGT" * (length // 4)), id=accession, name=accession.split(".")[0])
    rec.annotations["molecule_type"] = "DNA"
    return rec.format("genbank")


def install_fake_efetch(module: ModuleType, monkeypatch: pytest.MonkeyPatch, text: str) -> list:
    calls: list[dict] = []

    def fake_efetch(**kwargs: str) -> StringIO:
        calls.append(kwargs)
        return StringIO(text)

    monkeypatch.setattr(module.Entrez, "efetch", fake_efetch)
    return calls


def test_missing_email_exits_without_calling_ncbi(fd, monkeypatch, tmp_path) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    calls = install_fake_efetch(fd, monkeypatch, "")
    assert fd.main(["--out", str(tmp_path)]) == 2
    assert calls == []


def test_downloads_and_writes_lf_file(fd, monkeypatch, tmp_path, capsys) -> None:
    monkeypatch.setenv("NCBI_EMAIL", EMAIL)
    calls = install_fake_efetch(fd, monkeypatch, genbank_text("NC_000001.1", 1_200))

    assert fd.main(["--out", str(tmp_path)]) == 0

    assert calls == [
        {"db": "nuccore", "id": "NC_000001.1", "rettype": "gbwithparts", "retmode": "text"}
    ]
    data = (tmp_path / "mini.gbk").read_bytes()
    assert data.startswith(b"LOCUS")
    assert b"\r\n" not in data
    out = capsys.readouterr()
    assert EMAIL not in out.out + out.err


def test_existing_file_is_skipped_without_force(fd, monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("NCBI_EMAIL", EMAIL)
    calls = install_fake_efetch(fd, monkeypatch, genbank_text("NC_000001.1", 1_200))
    (tmp_path / "mini.gbk").write_text("old")

    assert fd.main(["--out", str(tmp_path)]) == 0
    assert calls == []
    assert (tmp_path / "mini.gbk").read_text() == "old"

    assert fd.main(["--out", str(tmp_path), "--force"]) == 0
    assert (tmp_path / "mini.gbk").read_text().startswith("LOCUS")


@pytest.mark.parametrize(
    ("accession", "length", "message"),
    [("NC_000001.1", 800, "expected 1,200 bp"), ("NC_999999.1", 1_200, "expected NC_000001.1")],
)
def test_unexpected_record_is_rejected(
    fd, monkeypatch, tmp_path, capsys, accession, length, message
) -> None:
    monkeypatch.setenv("NCBI_EMAIL", EMAIL)
    install_fake_efetch(fd, monkeypatch, genbank_text(accession, length))

    assert fd.main(["--out", str(tmp_path)]) == 1
    assert not (tmp_path / "mini.gbk").exists()
    assert message in capsys.readouterr().err
