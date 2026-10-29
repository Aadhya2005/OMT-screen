"""
screen.py
---------
Rank plant O-methyltransferases as candidate 6OMTs for expression in E. coli.

Two independent judgements, kept separate because they answer different
questions and a candidate can pass one and fail the other:

  CATALYTIC PLAUSIBILITY - does it have the right residues where 5ICE says
  the chemistry happens? Weighted by what each position does. The His/Asp
  pair at the acceptor hydroxyl is the catalysis and counts most; residues
  that merely line the pocket count least.

  EXPRESSIBILITY - would it fold in a bacterium? Plant enzymes routinely
  fail here, and this is the filter that the Fermelanta literature keeps
  running into. These checks are crude heuristics, labelled as such.

Two positions are deliberately NOT scored on identity. Gly165 and Cys253
touch the substrate through their backbone carbonyl, not their side chains,
so any residue there presents the same atom and scoring the amino acid would
be measuring nothing. They are reported and excluded.

Nothing here predicts activity. The output is a ranked shortlist with the
reasoning attached, which is what a screen is for.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

import sys
sys.path.insert(0, str(Path(__file__).parent))
from sequences import load_all, make_aligner, map_positions

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"

BACKBONE = {"N", "CA", "C", "O", "OXT"}
BLOSUM = substitution_matrices.load("BLOSUM62")

# Weights by what the position does, set from the structure rather than guessed.
W_CATALYTIC = 3.0     # polar atom within 3.5 A of the acceptor hydroxyl
W_POLAR = 2.0         # side-chain polar contact with the substrate
W_SHAPE = 1.0         # side chain lines the pocket

# Kyte-Doolittle hydropathy, for the crude membrane-segment check.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
POLAR_AA = set("STNQCYHKRDE")

# Identity below this and the position mapping is not trustworthy.
MIN_IDENTITY = 25.0
MIN_COVERAGE = 60.0


def position_roles(site: dict) -> dict[int, dict]:
    """Decide each position's weight and whether its side chain is involved."""
    by_num = {}
    acceptor = {r["resnum"] for r in site["acceptor_neighbours"]}

    for row in site["substrate_pocket"]:
        num = row["resnum"]
        atom = row["protein_atom"]
        backbone_only = atom in BACKBONE
        # Whether a contact is polar depends on the ATOM doing the touching,
        # not on the residue's type. Thr311 reaches the substrate with its
        # CG2 methyl while its hydroxyl points away; counting that as a
        # hydrogen-bonding contact would score the wrong chemistry.
        polar_atom = (atom[0].upper() in {"N", "O", "S"}) and not backbone_only

        if num in acceptor and not backbone_only:
            weight, role = W_CATALYTIC, "catalytic"
        elif backbone_only:
            weight, role = 0.0, "backbone contact only"
        elif polar_atom:
            weight, role = W_POLAR, "polar contact"
        else:
            weight, role = W_SHAPE, "pocket lining"
        by_num[num] = {
            "reference_aa": row["aa"],
            "weight": weight,
            "role": role,
            "distance": row["min_dist"],
            "via": row["protein_atom"],
        }

    # acceptor neighbours contacting through a side chain but outside the
    # 4.5 A substrate list still matter
    for row in site["acceptor_neighbours"]:
        num = row["resnum"]
        if num not in by_num and row["atom"] not in BACKBONE:
            by_num[num] = {
                "reference_aa": row["aa"], "weight": W_CATALYTIC,
                "role": "catalytic", "distance": row["dist"],
                "via": row["atom"],
            }
    return by_num


def similarity(a: str, b: str) -> float:
    """1.0 identical, 0.5 conservative, 0.0 otherwise."""
    if a == b:
        return 1.0
    try:
        s = BLOSUM[a, b]
    except (KeyError, IndexError):
        return 0.0
    return 0.5 if s >= 1 else 0.0


def hydropathy_max(seq: str, window: int = 19) -> float:
    """Highest windowed hydropathy: a crude stand-in for a membrane segment."""
    if len(seq) < window:
        return float("nan")
    vals = np.array([KD.get(c, 0.0) for c in seq])
    kernel = np.ones(window) / window
    return float(np.convolve(vals, kernel, mode="valid").max())


def expressibility(seq: str, ref_len: int) -> dict:
    """
    Crude checks on whether this would behave in E. coli.

    Every one of these is a heuristic, not a prediction. They are here to
    flag candidates worth a second look, not to rule anything out on their
    own, and the report says so.
    """
    flags = []
    n = len(seq)

    if n < 0.6 * ref_len:
        flags.append(f"short ({n} aa vs reference {ref_len})")
    if n > 1.6 * ref_len:
        flags.append(f"long ({n} aa vs reference {ref_len})")

    hyd = hydropathy_max(seq)
    if hyd >= 1.6:
        flags.append(f"hydrophobic stretch (max window {hyd:.2f})")

    cys = seq.count("C")
    cys_frac = cys / n if n else 0.0
    if cys >= 8 and cys_frac > 0.025:
        flags.append(f"cysteine-rich ({cys} Cys); disulfides do not form "
                     "in the E. coli cytoplasm")

    # a crude signal-peptide shape: hydrophobic run early on
    if n > 30 and hydropathy_max(seq[:30], window=12) >= 2.0:
        flags.append("hydrophobic N-terminus, possible targeting peptide")

    return {"flags": flags, "max_hydropathy": round(hyd, 2) if hyd == hyd else None,
            "cys_count": cys, "length": n}


def classify_control(rec) -> str | None:
    """
    Which sequences have published behaviour we can test the ranking against.

    Positives are the three reviewed, protein-level 6OMTs. Near misses are
    the ipecac enzymes: same EC listing, but they methylate a different
    position of a different alkaloid, so a scoring scheme that cannot tell
    them from the real thing is not discriminating.
    """
    name = rec.name.lower()
    if "norcoclaurine 6-o-methyltransferase" in name and rec.reviewed \
            and rec.evidence == 1:
        return "positive"
    if "cephaeline" in name or "demethylcephaeline" in name:
        return "near miss"
    if "norcoclaurine 6-o-methyltransferase" in name:
        return "putative"
    return None


def main() -> None:
    site = json.loads((OUT / "active_site.json").read_text())
    roles = position_roles(site)
    scored = {n: r for n, r in roles.items() if r["weight"] > 0}
    excluded = {n: r for n, r in roles.items() if r["weight"] == 0}

    print("Positions and what each is worth")
    print("-" * 64)
    for num in sorted(roles):
        r = roles[num]
        print(f"  {r['reference_aa']}{num:<5} weight {r['weight']:<4} "
              f"{r['role']:<22} {r['distance']:>5.2f} A via {r['via']}")
    print(f"\n  scoring on {len(scored)} positions; "
          f"{len(excluded)} excluded as backbone-only "
          f"({', '.join(f'{roles[n]['reference_aa']}{n}' for n in sorted(excluded))})")

    records, ref = load_all()
    aligner = make_aligner()
    total_weight = sum(r["weight"] for r in scored.values())
    positions = sorted(scored)

    rows = []
    for rec in records:
        m = map_positions(ref.seq, rec.seq, positions, aligner)
        reliable = (m["identity"] >= MIN_IDENTITY
                    and m["coverage"] >= MIN_COVERAGE)

        got = 0.0
        matches, per_pos = 0, {}
        for pos in positions:
            aa = m["residues"][pos]
            w = scored[pos]["weight"]
            ref_aa = scored[pos]["reference_aa"]
            s = similarity(ref_aa, aa) if aa else 0.0
            got += w * s
            if aa == ref_aa:
                matches += 1
            per_pos[pos] = aa or "-"

        exp = expressibility(rec.seq, ref.length)
        rows.append({
            "accession": rec.acc,
            "name": rec.name,
            "organism": rec.organism,
            "reviewed": rec.reviewed,
            "evidence": rec.evidence,
            "control": classify_control(rec),
            "identity_to_ref": m["identity"],
            "coverage": m["coverage"],
            "mapping_reliable": reliable,
            "site_score": round(100.0 * got / total_weight, 1),
            "exact_matches": matches,
            "positions_scored": len(positions),
            "catalytic_intact": all(
                m["residues"][p] == scored[p]["reference_aa"]
                for p in positions if scored[p]["role"] == "catalytic"),
            "gapped_at": ",".join(str(p) for p in m["gapped_at"]),
            "expressibility_flags": "; ".join(exp["flags"]),
            "n_flags": len(exp["flags"]),
            "length": exp["length"],
            "max_hydropathy": exp["max_hydropathy"],
            "cys_count": exp["cys_count"],
            "residues": "".join(per_pos[p] for p in positions),
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(["site_score", "identity_to_ref"], ascending=False)
    df.to_csv(OUT / "screen.csv", index=False)

    meta = {
        "reference": {"accession": ref.acc, "organism": ref.organism,
                      "length": ref.length},
        "positions_scored": positions,
        "positions_excluded_backbone": sorted(excluded),
        "weights": {str(n): scored[n]["weight"] for n in positions},
        "total_weight": total_weight,
        "candidates": len(df),
        "thresholds": {"min_identity": MIN_IDENTITY, "min_coverage": MIN_COVERAGE},
    }
    (OUT / "screen_meta.json").write_text(json.dumps(meta, indent=2))

    print(f"\nscreened {len(df)} candidates; wrote {OUT / 'screen.csv'}")
    unreliable = int((~df["mapping_reliable"]).sum())
    print(f"  {unreliable} candidates have an unreliable mapping and are "
          "flagged, not scored on merit")


if __name__ == "__main__":
    main()
