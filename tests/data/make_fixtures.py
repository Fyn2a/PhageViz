"""Regenerate the synthetic GenBank fixtures in this directory.

    python tests/data/make_fixtures.py

Sequences are seeded, so the output is identical on every run. Designed open reading frames
are built from ``CODONS`` (a fixed codon -> amino acid table), so their expected proteins are
known without running a translator; they are stored in a ``/note="expected_aa=..."``
qualifier that the tests read.
"""

from __future__ import annotations

import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent

# Standard-code sense codons used to build stop-free ORFs.
CODONS = {
    "GCT": "A", "AAA": "K", "GAA": "E", "CTG": "L", "TCT": "S", "GGT": "G", "CGT": "R",
    "ACT": "T", "GTT": "V", "CAG": "Q", "GAT": "D", "TTC": "F", "ATC": "I", "AAC": "N",
    "CCG": "P", "TGG": "W", "CAT": "H", "TAC": "Y", "TGC": "C", "ATG": "M",
}  # fmt: skip


def background(length: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    return [rng.choice("ACGT") for _ in range(length)]


def orf(start_codon: str, n_codons: int, seed: int) -> tuple[str, str]:
    """A complete ORF: start codon, ``n_codons`` sense codons, TAA. Returns (dna, protein)."""
    rng = random.Random(seed)
    body = [rng.choice(sorted(CODONS)) for _ in range(n_codons)]
    dna = start_codon + "".join(body) + "TAA"
    return dna, "M" + "".join(CODONS[c] for c in body)  # any start codon initiates with Met


def place(seq: list[str], start: int, bases: str) -> None:
    """Write ``bases`` at 1-based ``start``."""
    seq[start - 1 : start - 1 + len(bases)] = bases


def loc(start: int, end: int, strand: int) -> SimpleLocation:
    """1-based inclusive coordinates to a Biopython location."""
    return SimpleLocation(start - 1, end, strand)


def circular(record_id: str, seq: list[str], description: str) -> SeqRecord:
    return SeqRecord(
        Seq("".join(seq)),
        id=record_id,
        name=record_id.split(".")[0],
        description=description,
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )


def wraps_origin() -> list[SeqRecord]:
    """Two circular 3,000 bp genomes, each with one designed ORF crossing the origin."""
    # Record 1: + strand 300 nt ORF split 200/100 across the origin, join(2801..3000,1..100).
    dna1, aa1 = orf("ATG", 98, seed=11)
    seq1 = background(3_000, seed=1)
    seq1[2_800:3_000] = dna1[:200]
    seq1[0:100] = dna1[200:]
    first = circular("WRAP_1.1", seq1, "Synthetic circular phage, + strand origin crossing")
    first.features = [
        # join(2801..3000,1..100) -> start 2801, end 3100, wraps
        SeqFeature(
            CompoundLocation([loc(2_801, 3_000, 1), loc(1, 100, 1)]),
            type="CDS",
            qualifiers={
                "locus_tag": ["wrapF"],
                "product": ["portal protein"],
                "note": [f"expected_aa={aa1}"],
            },
        ),
        # Ends exactly on the last base: not wrapped.
        SeqFeature(
            loc(2_401, 3_000, -1),
            type="CDS",
            qualifiers={"locus_tag": ["edge"], "product": ["tail fiber"]},
        ),
        SeqFeature(
            loc(1_001, 1_300, 1),
            type="CDS",
            qualifiers={"locus_tag": ["mid"], "product": ["integrase"]},
        ),
        SeqFeature(
            loc(1_501, 1_576, 1),
            type="tRNA",
            qualifiers={"locus_tag": ["trna1"], "product": ["tRNA-Gly"]},
        ),
    ]

    # Record 2: - strand ORF with a GTG start, complement(join(2801..3000,1..151)).
    # Reading order runs 151 -> 1 -> 3000 -> 2801, so the + strand text over
    # [2801..3000][1..151] is the reverse complement of the ORF.
    seq2 = background(3_000, seed=2)
    dna2, aa2 = orf("GTG", 115, seed=22)  # 351 nt
    genomic = str(Seq(dna2).reverse_complement())
    seq2[2_800:3_000] = genomic[:200]
    seq2[0:151] = genomic[200:]
    second = circular("WRAP_2.1", seq2, "Synthetic circular phage, - strand origin crossing")
    second.features = [
        # -> start 2801, end 3151, strand -1, wraps; GTG start gives Met
        SeqFeature(
            CompoundLocation([loc(1, 151, -1), loc(2_801, 3_000, -1)]),
            type="CDS",
            qualifiers={
                "locus_tag": ["wrapR"],
                "product": ["holin"],
                "note": [f"expected_aa={aa2}"],
            },
        ),
    ]
    return [first, second]


def pharokka_sample() -> SeqRecord:
    """Linear 2,000 bp genome with Pharokka-style qualifiers and ID edge cases."""
    seq = background(2_000, seed=2)
    place(seq, 101, "ATGAAATAA")  # + strand: M K stop
    place(seq, 201, "TTAGGCCAT")  # reverse complement of ATG GCC TAA: M A stop
    place(seq, 1_011, "GTGAAATAA")  # alternative start, complete gene: M K
    place(seq, 1_021, "GTGGCCTAA")  # alternative start, 5'-partial gene: V A
    record = SeqRecord(
        Seq("".join(seq)),
        id="PHAGE_1.1",
        name="PHAGE_1",
        description="Synthetic phage annotated by Pharokka",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    record.features = [
        SeqFeature(loc(1, 2_000, 1), type="source"),
        SeqFeature(loc(101, 109, 1), type="gene", qualifiers={"locus_tag": ["pk_0001"]}),
        # No translation -> "MK"; function wins over keywords.
        SeqFeature(
            loc(101, 109, 1),
            type="CDS",
            qualifiers={
                "locus_tag": ["pk_0001"],
                "product": ["portal protein"],
                "phrog": ["phrog_15"],
                "function": ["head and packaging"],
            },
        ),
        # Minus strand, no translation -> "MA"; function in different case; bare PHROG number.
        SeqFeature(
            loc(201, 209, -1),
            type="CDS",
            qualifiers={
                "locus_tag": ["pk_0002"],
                "product": ["hypothetical protein"],
                "phrog": ["136"],
                "function": ["Tail"],
            },
        ),
        # Trailing '*' in translation is stripped; No_PHROG -> missing.
        SeqFeature(
            loc(301, 600, 1),
            type="CDS",
            qualifiers={
                "locus_tag": ["pk_0003"],
                "product": ["hypothetical protein"],
                "phrog": ["No_PHROG"],
                "function": ["unknown function"],
                "translation": ["MSTV*"],
            },
        ),
        # No locus_tag -> "PHAGE_1.1_0004"; no function -> keyword "lysis".
        SeqFeature(
            loc(701, 1_000, 1),
            type="CDS",
            qualifiers={"product": ["endolysin"], "translation": ["MEND"]},
        ),
        # GTG start on a complete gene -> "MK"; on a 5'-partial gene (<1021) -> "VA".
        SeqFeature(loc(1_011, 1_019, 1), type="CDS", qualifiers={"locus_tag": ["pk_gtg"]}),
        SeqFeature(
            SimpleLocation(BeforePosition(1_020), 1_029, 1),
            type="CDS",
            qualifiers={"locus_tag": ["pk_partial"]},
        ),
        # Duplicate locus tags -> "dup", "dup_2".
        SeqFeature(
            loc(1_101, 1_400, 1),
            type="CDS",
            qualifiers={"locus_tag": ["dup"], "translation": ["MD"]},
        ),
        SeqFeature(
            loc(1_401, 1_700, -1),
            type="CDS",
            qualifiers={"locus_tag": ["dup"], "translation": ["MDD"]},
        ),
        # RNA gene without product -> product "tRNA", category other/keyword.
        SeqFeature(
            loc(1_801, 1_876, 1),
            type="tRNA",
            qualifiers={"locus_tag": ["pk_trna"], "function": ["tRNA"]},
        ),
        SeqFeature(loc(1_900, 1_950, 1), type="misc_feature"),
    ]
    return record


def write(records: list[SeqRecord], name: str) -> None:
    with (HERE / name).open("w", encoding="utf-8", newline="\n") as handle:
        SeqIO.write(records, handle, "genbank")


if __name__ == "__main__":
    write(wraps_origin(), "wraps_origin.gbk")
    write([pharokka_sample()], "pharokka_sample.gbk")
    print("wrote wraps_origin.gbk, pharokka_sample.gbk")
