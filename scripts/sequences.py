"""
sequences.py
------------
Load the candidate set and tie it to the structure.

The screen asks, of every candidate: what do you have at the positions that
5ICE says matter? Answering that needs a correspondence between each
candidate's numbering and the reference's.

That is done by pairwise alignment to the reference rather than one multiple
alignment of everything. The reason is accountability. In a multiple
alignment, a badly aligned candidate quietly corrupts shared columns and
nothing says which ones. Aligned one at a time against the reference, every
candidate carries its own identity and coverage, so a mapping that cannot be
trusted can be flagged and set aside instead of silently scored.

Before any of that, the reference sequence is checked against the structure's
own residues. If UniProt Q5C9L7 and the model in 5ICE disagree at the
positions being scored, every number downstream is wrong, so the check runs
first and refuses to continue on a mismatch.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, index_to_one, three_to_index

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "output"

REFERENCE_ACC = "Q5C9L7"      # Thalictrum flavum 6OMT, the 5ICE sequence
NON_PROTEIN = {"HOH", "EDO", "K", "SAH", "2H4"}


@dataclass
class Record:
    acc: str
    name: str
    organism: str
    seq: str
    source: str                     # which file it came from
    reviewed: bool
    evidence: int | None            # UniProt PE level, 1 = protein evidence
    gene: str = ""
    flags: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.seq)


HEADER = re.compile(
    r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<name>.*?)"
    r"(?:\s+OS=(?P<os>.*?))?(?:\s+OX=(?P<ox>\d+))?(?:\s+GN=(?P<gn>\S+))?"
    r"(?:\s+PE=(?P<pe>\d+))?(?:\s+SV=\d+)?$"
)


def parse_header(desc: str) -> dict:
    m = HEADER.match(desc.strip())
    if not m:
        return {"acc": desc.split()[0], "name": desc, "organism": "",
                "reviewed": False, "evidence": None, "gene": ""}
    g = m.groupdict()
    return {
        "acc": g["acc"],
        "name": (g["name"] or "").strip(),
        "organism": (g["os"] or "").strip(),
        "reviewed": g["db"] == "sp",
        "evidence": int(g["pe"]) if g["pe"] else None,
        "gene": g["gn"] or "",
    }


def load_fasta(path: Path, source: str) -> list[Record]:
    out = []
    for rec in SeqIO.parse(str(path), "fasta"):
        meta = parse_header(rec.description)
        out.append(Record(
            acc=meta["acc"], name=meta["name"], organism=meta["organism"],
            seq=str(rec.seq).upper(), source=source,
            reviewed=meta["reviewed"], evidence=meta["evidence"],
            gene=meta["gene"],
        ))
    return out


def structure_sequence(pdb: Path) -> dict[int, str]:
    """Residue number -> one-letter code, for the modelled chain."""
    chain = None
    for model in PDBParser(QUIET=True).get_structure("s", str(pdb)):
        for c in model:
            chain = c
            break
        break
    out = {}
    for res in chain:
        if res.get_resname().strip() in NON_PROTEIN:
            continue
        if not is_aa(res, standard=True):
            continue
        try:
            out[res.get_id()[1]] = index_to_one(three_to_index(
                res.get_resname().strip()))
        except (KeyError, ValueError):
            continue
    return out


def check_reference(ref: Record, pdb: Path, positions: list[int]) -> dict:
    """
    Confirm the reference sequence and the deposited model agree.

    PDB numbering here is supposed to equal UniProt numbering (the DBREF says
    residues 3-350 map 1:1). This verifies that at every scored position
    rather than trusting it.
    """
    model = structure_sequence(pdb)
    rows, mismatches = [], []
    for num in positions:
        in_struct = model.get(num)
        in_seq = ref.seq[num - 1] if 0 < num <= len(ref.seq) else None
        ok = in_struct is not None and in_struct == in_seq
        rows.append({"position": num, "structure": in_struct,
                     "sequence": in_seq, "match": ok})
        if not ok:
            mismatches.append(num)

    overlap = [n for n in model if 0 < n <= len(ref.seq)]
    agree = sum(1 for n in overlap if model[n] == ref.seq[n - 1])
    return {
        "scored_positions": rows,
        "mismatched_positions": mismatches,
        "modelled_residues": len(model),
        "overlap_checked": len(overlap),
        "overlap_identical": agree,
        "overlap_percent": round(100.0 * agree / len(overlap), 1) if overlap else 0.0,
        "ok": not mismatches,
    }


def make_aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.end_insertion_score = 0.0      # don't punish ragged termini
    a.end_deletion_score = 0.0
    return a


def map_positions(ref_seq: str, cand_seq: str, positions: list[int],
                  aligner: Align.PairwiseAligner) -> dict:
    """
    Align one candidate to the reference and read off what it has at each
    reference position.

    Returns the residue per position (None where the candidate has a gap),
    plus identity and coverage so the mapping's reliability can be judged.
    """
    aln = aligner.align(ref_seq, cand_seq)[0]
    ref_to_cand: dict[int, int] = {}
    identical = aligned = 0

    for (r0, r1), (c0, c1) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(r1 - r0):
            ri, ci = r0 + off, c0 + off
            ref_to_cand[ri + 1] = ci          # 1-based ref -> 0-based cand
            aligned += 1
            if ref_seq[ri] == cand_seq[ci]:
                identical += 1

    residues = {}
    for pos in positions:
        ci = ref_to_cand.get(pos)
        residues[pos] = cand_seq[ci] if ci is not None else None

    return {
        "residues": residues,
        "identity": round(100.0 * identical / aligned, 1) if aligned else 0.0,
        "coverage": round(100.0 * aligned / len(ref_seq), 1),
        "aligned_positions": aligned,
        "gapped_at": sorted(p for p, r in residues.items() if r is None),
        "score": float(aln.score),
    }


def load_all() -> tuple[list[Record], Record]:
    known = load_fasta(DATA / "known_6omt.fasta", "known_6omt")
    pool = load_fasta(DATA / "plant_omts.fasta", "plant_omts")

    by_acc: dict[str, Record] = {}
    for rec in known + pool:          # known first, so its label wins
        if rec.acc in by_acc:
            by_acc[rec.acc].flags.append("in both input files")
            continue
        by_acc[rec.acc] = rec

    records = list(by_acc.values())
    ref = by_acc.get(REFERENCE_ACC)
    if ref is None:
        raise SystemExit(f"reference {REFERENCE_ACC} not present in the input")
    return records, ref


def main() -> None:
    OUT.mkdir(exist_ok=True)
    site = json.loads((OUT / "active_site.json").read_text())
    positions = [p["resnum"] for p in site["scored_positions"]]

    records, ref = load_all()
    print(f"loaded {len(records)} unique sequences")
    print(f"reference: {ref.acc}  {ref.name}  ({ref.organism}), "
          f"{ref.length} aa")

    print("\nChecking the reference against the deposited model")
    print("-" * 50)
    chk = check_reference(ref, DATA / "5ICE.pdb", positions)
    print(f"  modelled residues        {chk['modelled_residues']}")
    print(f"  overlap checked          {chk['overlap_checked']}")
    print(f"  identical over overlap   {chk['overlap_identical']} "
          f"({chk['overlap_percent']}%)")
    for row in chk["scored_positions"]:
        mark = "ok" if row["match"] else "MISMATCH"
        print(f"    {row['position']:>4}  structure {row['structure']}  "
              f"sequence {row['sequence']}   {mark}")
    if not chk["ok"]:
        raise SystemExit(
            f"\nrefusing to continue: positions {chk['mismatched_positions']} "
            "differ between the structure and the reference sequence")
    print("\n  reference and structure agree at every scored position")

    (OUT / "reference_check.json").write_text(json.dumps(chk, indent=2))
    print(f"wrote {OUT / 'reference_check.json'}")


if __name__ == "__main__":
    main()
