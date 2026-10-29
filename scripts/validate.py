"""
validate.py
-----------
Does the ranking actually separate known 6OMTs from the rest?

A scoring scheme can look sensible and rank nothing correctly. The way to
find out is to include sequences whose behaviour is already published and
check where they land:

  positives  - reviewed, protein-level characterised (RS)-norcoclaurine
               6-O-methyltransferases. These should rank at the top.
  near miss  - the ipecac OMTs, listed under the same EC number but acting
               on a different position of a different alkaloid. A scheme
               that cannot separate these from the real thing is not
               discriminating, it is just measuring family membership.
  putative   - annotated as 6OMT but only inferred from homology. No claim
               either way; reported for interest.

If the positives do not rank highly, the honest conclusion is that the
scoring is wrong, and this script says so rather than letting a ranking be
published that its own controls contradict.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"


def banner(t: str) -> None:
    print(f"\n{t}\n" + "-" * len(t))


def main() -> int:
    df = pd.read_csv(OUT / "screen.csv")
    meta = json.loads((OUT / "screen_meta.json").read_text())
    n = len(df)

    df = df.sort_values(["site_score", "identity_to_ref"],
                        ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    df["percentile"] = 100.0 * (1 - df.index / (n - 1))

    ok = True

    banner(f"Controls among {n} candidates")
    controls = df[df["control"].notna()].sort_values("rank")
    print(f"  {'rank':>5}  {'score':>6}  {'ident':>6}  {'class':<10} "
          f"{'accession':<12} organism")
    for _, r in controls.iterrows():
        print(f"  {r['rank']:>5}  {r['site_score']:>6.1f}  "
              f"{r['identity_to_ref']:>6.1f}  {r['control']:<10} "
              f"{r['accession']:<12} {str(r['organism'])[:38]}")

    pos = df[df["control"] == "positive"]
    near = df[df["control"] == "near miss"]
    putative = df[df["control"] == "putative"]

    # --- test 1: do the positives rank at the top? ---------------------
    banner("Test 1: positives should rank near the top")
    if pos.empty:
        print("  no positives in the set; cannot test")
        ok = False
    else:
        worst = int(pos["rank"].max())
        frac = 100.0 * worst / n
        print(f"  {len(pos)} positives, worst rank {worst} of {n} "
              f"(top {frac:.1f}%)")
        for _, r in pos.iterrows():
            print(f"    rank {int(r['rank']):>4}  score {r['site_score']:>5.1f}  "
                  f"{r['accession']}  {str(r['organism'])[:40]}")
        if frac <= 10.0:
            print("  PASS: every characterised 6OMT lands in the top 10%")
        else:
            print("  FAIL: a characterised 6OMT ranks outside the top 10%")
            ok = False

    # --- test 2: are the near misses separated? ------------------------
    banner("Test 2: near misses should score below the positives")
    if near.empty or pos.empty:
        print("  not testable with this control set")
    else:
        lo_pos = float(pos["site_score"].min())
        hi_near = float(near["site_score"].max())
        print(f"  lowest positive score   {lo_pos:.1f}")
        print(f"  highest near-miss score {hi_near:.1f}")
        for _, r in near.iterrows():
            print(f"    rank {int(r['rank']):>4}  score {r['site_score']:>5.1f}  "
                  f"{r['accession']}  {str(r['name'])[:48]}")
        if hi_near < lo_pos:
            print("  PASS: the two groups are separated")
        else:
            print("  PARTIAL: the near misses overlap the positives.")
            print("  The active site alone does not distinguish which ring")
            print("  position gets methylated, which is a real limit of a")
            print("  residue-identity score and is reported as such.")

    # --- test 3: is the catalytic pair doing work? ---------------------
    banner("Test 3: does the catalytic pair separate anything?")
    intact = df["catalytic_intact"].sum()
    print(f"  {intact} of {n} candidates keep both catalytic residues")
    if not pos.empty:
        print(f"  positives keeping both: "
              f"{int(pos['catalytic_intact'].sum())} of {len(pos)}")
    if intact == n:
        print("  NOTE: every candidate has them, so this criterion separates")
        print("  nothing here. It is retained because losing them would be")
        print("  disqualifying, not because it discriminates in this set.")

    # --- test 4: mapping reliability -----------------------------------
    banner("Test 4: how much of the set can be scored at all?")
    bad = int((~df["mapping_reliable"]).sum())
    print(f"  {bad} of {n} have identity below "
          f"{meta['thresholds']['min_identity']}% or coverage below "
          f"{meta['thresholds']['min_coverage']}%")
    print("  Their position mapping is not trustworthy, so their scores are")
    print("  reported but should not be read as evidence.")
    if not pos.empty:
        print(f"  all positives reliable: "
              f"{bool(pos['mapping_reliable'].all())}")

    # --- test 5: does this beat plain sequence identity? ---------------
    banner("Test 5: does the active-site score add anything over identity?")
    rel = df[df["mapping_reliable"]].copy()
    corr = float(rel["site_score"].corr(rel["identity_to_ref"], method="spearman"))
    print(f"  Spearman correlation with global identity: {corr:.2f} "
          f"(over {len(rel)} reliably mapped candidates)")

    # Where would the positives rank if we sorted on identity alone?
    by_id = rel.sort_values("identity_to_ref", ascending=False).reset_index(drop=True)
    by_id["id_rank"] = by_id.index + 1
    by_site = rel.sort_values("site_score", ascending=False).reset_index(drop=True)
    by_site["site_rank"] = by_site.index + 1
    merged = by_id[["accession", "id_rank"]].merge(
        by_site[["accession", "site_rank", "control"]], on="accession")

    pos_m = merged[merged["control"] == "positive"]
    near_m = merged[merged["control"] == "near miss"]
    print(f"  {'accession':<12} {'by identity':>12} {'by site score':>14}  class")
    for _, r in pd.concat([pos_m, near_m]).iterrows():
        print(f"  {r['accession']:<12} {int(r['id_rank']):>12} "
              f"{int(r['site_rank']):>14}  {r['control']}")

    if corr > 0.9:
        print("\n  The two orderings agree almost completely. The structural")
        print("  weighting is not adding information here, and the honest")
        print("  reading is that this ranking is close to a similarity search.")
        ok = False
    else:
        print(f"\n  The orderings differ (rho = {corr:.2f}), so the site score")
        print("  is not simply reproducing overall similarity.")

    # --- the shortlist --------------------------------------------------
    banner("Shortlist: top candidates that are not already known 6OMTs")
    fresh = df[(df["control"].isna()) & (df["mapping_reliable"])].head(12)
    print(f"  {'rank':>5} {'score':>6} {'ident':>6} {'flags':>6}  "
          f"{'accession':<12} name")
    for _, r in fresh.iterrows():
        print(f"  {r['rank']:>5} {r['site_score']:>6.1f} "
              f"{r['identity_to_ref']:>6.1f} {int(r['n_flags']):>6}  "
              f"{r['accession']:<12} {str(r['name'])[:44]}")

    df.to_csv(OUT / "screen_ranked.csv", index=False)
    summary = {
        "candidates": n,
        "positives": [
            {"accession": r["accession"], "rank": int(r["rank"]),
             "score": float(r["site_score"]), "organism": r["organism"]}
            for _, r in pos.iterrows()
        ],
        "near_misses": [
            {"accession": r["accession"], "rank": int(r["rank"]),
             "score": float(r["site_score"]), "name": r["name"]}
            for _, r in near.iterrows()
        ],
        "putative": [
            {"accession": r["accession"], "rank": int(r["rank"]),
             "score": float(r["site_score"])} for _, r in putative.iterrows()
        ],
        "worst_positive_rank": int(pos["rank"].max()) if not pos.empty else None,
        "unreliable_mappings": bad,
        "controls_pass": ok,
    }
    (OUT / "validation.json").write_text(json.dumps(summary, indent=2))

    banner("Verdict")
    print("  CONTROLS PASS" if ok else "  CONTROLS FAIL")
    print(f"  wrote {OUT / 'screen_ranked.csv'} and {OUT / 'validation.json'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
