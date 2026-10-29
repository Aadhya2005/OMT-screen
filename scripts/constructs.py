"""
constructs.py
-------------
Turn the shortlist into expression constructs you could actually order.

The screen ends with a ranked list, which is an answer to "what should we
test" but not something a lab can act on. This step closes that gap: each
chosen candidate is codon-optimised for E. coli, cleared of the restriction
sites used for cloning, given primers, and assembled in silico into
pET-28a(+), with the finished plasmid written as GenBank.

The whole construct-design pipeline is carried over unchanged from the
earlier RsTAL project, which is the point: it was built to be reused, and
reusing it on a different protein is a test of whether it actually was.

Four constructs are built, chosen to make the experiment interpretable
rather than to maximise the count:

  Q9LEL6  Coptis japonica 6OMT      the enzyme the published E. coli
                                    platform uses; the benchmark
  Q5C9L7  Thalictrum flavum 6OMT    the sequence the structure came from
  Q6WUC1  Papaver somniferum 6OMT   a third characterised enzyme
  A0A0B4VG62  Sinopodophyllum       annotated as a 6OMT but never
                                    characterised; the one genuine
                                    prediction in the set

Three knowns and one unknown. If the three knowns express and the unknown
does not, the screen was wrong about it, and that is a result worth having.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"
CONSTRUCTS = OUT / "constructs"
# pET-28a(+) is not committed. Point this at your own copy; SETUP.md says
# where to download it. A path relative to the project root works.
VECTOR = ROOT / "data" / "pET-28a.dna"

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent / "pipeline"))

from sequences import load_all                                   # noqa: E402
from pipeline.assemble import (annotate, build_construct,        # noqa: E402
                               digest_insert, find_expression_product,
                               pcr_product, verify)
from pipeline.codon_tools import summarise, translate_cds        # noqa: E402
from pipeline.optimise import find_sites, optimise_cds           # noqa: E402
from pipeline.primers import (design_cloning_primers,            # noqa: E402
                              design_internal_seq_primers,
                              pair_heterodimer_dg, primer_table_rows)

CHOSEN = {
    "Q9LEL6": "Coptis japonica 6OMT, the published platform enzyme",
    "Q5C9L7": "Thalictrum flavum 6OMT, the sequence 5ICE was solved from",
    "Q6WUC1": "Papaver somniferum 6OMT",
    "A0A0B4VG62": "Sinopodophyllum hexandrum, annotated 6OMT, uncharacterised",
}

T7_PROMOTER_CORE = "TAATACGACTCACTATAG"
SEED = 1


def load_vector(path: Path):
    """
    Load pET-28a(+), flipping it if the cassette is stored on the bottom
    strand. The SnapGene distribution is, which the earlier project found
    the hard way.
    """
    fmt = {"dna": "snapgene", "gb": "genbank", "gbk": "genbank"}.get(
        path.suffix.lower().lstrip("."), "fasta")
    rec = next(SeqIO.parse(str(path), fmt))
    rec.annotations.setdefault("molecule_type", "DNA")

    if T7_PROMOTER_CORE not in str(rec.seq).upper():
        if T7_PROMOTER_CORE in str(rec.seq.reverse_complement()).upper():
            rc = rec.reverse_complement()
            rc.id, rc.name = path.stem, path.stem
            rc.annotations["molecule_type"] = "DNA"
            rc.annotations["topology"] = "circular"
            rec = rc
        else:
            raise SystemExit("T7 promoter not found on either strand")
    rec.annotations.setdefault("topology", "circular")
    if rec.id in ("", "<unknown id>"):
        rec.id = path.stem
    return rec


def build_one(rec, vector, note: str) -> dict:
    name = f"pET28a_{rec.acc}"
    print(f"\n{'=' * 66}\n{rec.acc}  {rec.organism}\n{note}\n{'=' * 66}")
    print(f"  protein {len(rec.seq)} aa")

    opt = optimise_cds(rec.seq, seed=SEED)
    cds = opt["optimised_cds"]
    if translate_cds(cds) != rec.seq:
        raise SystemExit(f"{rec.acc}: translation changed during optimisation")

    naive_stats = summarise(opt["naive_cds"], "naive")
    stats = summarise(cds, "optimised")
    print(f"  CAI {naive_stats['CAI']:.3f} -> {stats['CAI']:.3f}   "
          f"GC {naive_stats['GC_percent']:.1f}% -> {stats['GC_percent']:.1f}%")
    dropped = opt.get("dropped_constraints") or []
    if dropped:
        print(f"  solver fell back, dropped: {'; '.join(dropped)}")
    leftover = find_sites(cds)
    print(f"  forbidden sites remaining: {leftover or 'none'}")

    fwd, rev = design_cloning_primers(cds, label=rec.acc)
    internals = design_internal_seq_primers(cds, label=rec.acc)
    het = pair_heterodimer_dg(fwd, rev)
    warn = [p.name for p in (fwd, rev, *internals) if p.warnings]
    print(f"  primers: {2 + len(internals)} designed, "
          f"heterodimer dG {het:.1f} kcal/mol, "
          f"{'warnings on ' + ', '.join(warn) if warn else 'all clean'}")

    amplicon = pcr_product(cds, fwd.sequence, rev.sequence)
    fragment = digest_insert(amplicon)
    plasmid = annotate(build_construct(vector, fragment, name=name),
                       len(rec.seq))

    checks = verify(plasmid, rec.seq, vector)
    passed = sum(c.passed for c in checks)
    for c in checks:
        if not c.passed:
            print(f"    FAILED: {c.line()}")
    print(f"  verification: {passed}/{len(checks)} checks passed, "
          f"plasmid {len(plasmid.seq)} bp")

    CONSTRUCTS.mkdir(parents=True, exist_ok=True)
    SeqIO.write(plasmid, str(CONSTRUCTS / f"{name}.gb"), "genbank")
    (CONSTRUCTS / f"{name}_CDS.fasta").write_text(
        f">{rec.acc}_codon_optimised len={len(cds)}nt seed={SEED}\n"
        + "\n".join(cds[i:i + 60] for i in range(0, len(cds), 60)) + "\n")
    pd.DataFrame(primer_table_rows([fwd, rev, *internals])).to_csv(
        CONSTRUCTS / f"{name}_primers.csv", index=False)

    expr = find_expression_product(str(plasmid.seq))
    return {
        "accession": rec.acc,
        "organism": rec.organism,
        "note": note,
        "protein_aa": len(rec.seq),
        "cds_nt": len(cds),
        "cai_naive": naive_stats["CAI"],
        "cai_optimised": stats["CAI"],
        "gc_naive": naive_stats["GC_percent"],
        "gc_optimised": stats["GC_percent"],
        "forbidden_sites_remaining": leftover,
        "dropped_constraints": dropped,
        "plasmid_bp": len(plasmid.seq),
        "checks_passed": passed,
        "checks_total": len(checks),
        "all_passed": passed == len(checks),
        "heterodimer_dG": round(het, 2),
        "primer_warnings": warn,
        "expressed_aa": expr.get("length_aa"),
        "n_terminal_leader": expr.get("leader"),
        "files": [f"{name}.gb", f"{name}_CDS.fasta", f"{name}_primers.csv"],
    }


def main() -> None:
    if not VECTOR.exists():
        raise SystemExit(
            f"vector not found at {VECTOR}.\n"
            "Download pET-28a(+) from SnapGene's plasmid library or Addgene "
            "and set VECTOR to its path.")

    records, _ = load_all()
    by_acc = {r.acc: r for r in records}
    vector = load_vector(VECTOR)
    print(f"vector {vector.id}, {len(vector.seq)} bp")

    ranked = pd.read_csv(OUT / "screen_ranked.csv").set_index("accession")

    results = []
    for acc, note in CHOSEN.items():
        rec = by_acc.get(acc)
        if rec is None:
            print(f"\n{acc}: not in the candidate set, skipped")
            continue
        row = build_one(rec, vector, note)
        if acc in ranked.index:
            row["screen_rank"] = int(ranked.loc[acc, "rank"])
            row["site_score"] = float(ranked.loc[acc, "site_score"])
        results.append(row)

    (OUT / "constructs.json").write_text(json.dumps(results, indent=2))

    print(f"\n{'=' * 66}\nSummary\n{'=' * 66}")
    print(f"  {'accession':<12} {'rank':>5} {'score':>6} {'aa':>5} "
          f"{'plasmid':>8} {'checks':>8}")
    for r in results:
        print(f"  {r['accession']:<12} {r.get('screen_rank', '-'):>5} "
              f"{r.get('site_score', 0):>6.1f} {r['protein_aa']:>5} "
              f"{r['plasmid_bp']:>8} "
              f"{r['checks_passed']}/{r['checks_total']:>6}")
    ok = all(r["all_passed"] for r in results)
    print(f"\n  {'every construct verified' if ok else 'SOME CHECKS FAILED'}")
    print(f"  GenBank files in {CONSTRUCTS}")


if __name__ == "__main__":
    main()
