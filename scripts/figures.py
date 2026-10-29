"""
figures.py
----------
Figures for the screen, drawn from screen_ranked.csv and active_site.json.
No number is typed here; everything is read from the computed output.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output"

INK = "#1a1a1a"
MUTED = "#6b6b76"
TEAL = "#1f6f8b"
WARM = "#c75b39"
GREEN = "#4a7c59"
PURPLE = "#8a6fa8"
LIGHT = "#e8e6e1"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
})

CONTROL_STYLE = {
    "positive":  (GREEN,  "characterised 6OMT", 120, "o"),
    "near miss": (WARM,   "same EC, different reaction", 110, "s"),
    "putative":  (PURPLE, "annotated, not characterised", 95, "^"),
}


def load():
    df = pd.read_csv(OUT / "screen_ranked.csv")
    site = json.loads((OUT / "active_site.json").read_text())
    val = json.loads((OUT / "validation.json").read_text())
    return df, site, val


def fig_ranking(df: pd.DataFrame) -> Path:
    """Where every candidate scores, and where the controls landed."""
    fig, ax = plt.subplots(figsize=(9.6, 4.4))

    rel = df[df["mapping_reliable"]]
    unrel = df[~df["mapping_reliable"]]

    ax.scatter(unrel["rank"], unrel["site_score"], s=22, color=LIGHT,
               edgecolors="none", zorder=1, label="mapping unreliable")
    plain = rel[rel["control"].isna()]
    ax.scatter(plain["rank"], plain["site_score"], s=26, color=MUTED,
               alpha=0.55, edgecolors="none", zorder=2,
               label="candidate")

    for kind, (colour, label, size, marker) in CONTROL_STYLE.items():
        sub = df[df["control"] == kind]
        if sub.empty:
            continue
        ax.scatter(sub["rank"], sub["site_score"], s=size, color=colour,
                   marker=marker, edgecolors="white", linewidths=1.4,
                   zorder=4, label=label)

    pos = df[df["control"] == "positive"]
    if not pos.empty:
        worst = int(pos["rank"].max())
        ax.axvline(worst + 0.5, color=GREEN, lw=1.1, ls="--", alpha=0.8, zorder=3)
        ax.annotate(f"all characterised 6OMTs\nrank {worst} or better\n"
                    f"(top {100.0 * worst / len(df):.1f}%)",
                    (worst + 16, 86), fontsize=9.5, color=GREEN, va="top")

    ax.set_xlabel("rank")
    ax.set_ylabel("active-site score")
    ax.set_ylim(-3, 105)
    ax.set_xlim(-4, len(df) + 4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right", ncol=2)
    ax.set_title("Every candidate scored, with the controls marked",
                 fontsize=11.5, color=INK, pad=10)
    fig.tight_layout()
    p = OUT / "fig_ranking.png"
    fig.savefig(p, dpi=220, facecolor="white")
    plt.close(fig)
    return p


def fig_vs_identity(df: pd.DataFrame, val: dict) -> Path:
    """Does the site score just reproduce sequence identity?"""
    rel = df[df["mapping_reliable"]]
    rho = float(rel["site_score"].corr(rel["identity_to_ref"], method="spearman"))

    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    plain = rel[rel["control"].isna()]
    ax.scatter(plain["identity_to_ref"], plain["site_score"], s=30,
               color=MUTED, alpha=0.5, edgecolors="none", zorder=2)

    for kind, (colour, label, size, marker) in CONTROL_STYLE.items():
        sub = rel[rel["control"] == kind]
        if sub.empty:
            continue
        ax.scatter(sub["identity_to_ref"], sub["site_score"], s=size + 20,
                   color=colour, marker=marker, edgecolors="white",
                   linewidths=1.4, zorder=4, label=label)

    ax.set_xlabel("global sequence identity to the reference (%)")
    ax.set_ylabel("active-site score")
    ax.set_title(f"Spearman rho = {rho:.2f}\nthe two measures disagree enough "
                 "to be worth computing separately",
                 fontsize=10.5, color=INK, pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5, loc="lower right")
    fig.tight_layout()
    p = OUT / "fig_score_vs_identity.png"
    fig.savefig(p, dpi=220, facecolor="white")
    plt.close(fig)
    return p


def fig_positions(df: pd.DataFrame, site: dict) -> Path:
    """What the top candidates actually have at each scored position."""
    meta = json.loads((OUT / "screen_meta.json").read_text())
    positions = meta["positions_scored"]
    weights = meta["weights"]
    ref_aa = {p["resnum"]: p["aa"] for p in site["substrate_pocket"]}

    top = df[df["mapping_reliable"]].head(16).copy()
    rows = [list(str(r)) for r in top["residues"]]
    labels = [f"{a[:28]}" for a in top["name"].fillna("")]
    accs = list(top["accession"])
    ctrl = list(top["control"].fillna(""))

    fig, ax = plt.subplots(figsize=(10.2, 6.6))
    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -1.9)
    ax.axis("off")

    # header: position, reference residue, weight
    for j, pos in enumerate(positions):
        w = weights[str(pos)]
        colour = GREEN if w >= 3 else (TEAL if w >= 2 else MUTED)
        ax.text(j, -1.55, f"{ref_aa.get(pos, '?')}{pos}", ha="center",
                va="center", fontsize=9.5, color=colour, rotation=90,
                fontweight="bold" if w >= 3 else "normal")
        ax.text(j, -0.78, f"w{w:g}", ha="center", va="center",
                fontsize=8, color=MUTED)

    for i, row in enumerate(rows):
        for j, aa in enumerate(row):
            ref = ref_aa.get(positions[j], "?")
            if aa == ref:
                fc, tc = TEAL, "white"
            elif aa == "-":
                fc, tc = "#f0efec", MUTED
            else:
                fc, tc = "#f6f1ec", INK
            ax.add_patch(plt.Rectangle((j - 0.46, i - 0.42), 0.92, 0.84,
                                       facecolor=fc, edgecolor="white",
                                       linewidth=1.2, zorder=2))
            ax.text(j, i, aa, ha="center", va="center", fontsize=10,
                    color=tc, zorder=3,
                    fontweight="bold" if aa == ref else "normal")

        c = ctrl[i]
        colour = CONTROL_STYLE.get(c, (INK,))[0] if c else INK
        ax.text(-0.75, i, f"{accs[i]}", ha="right", va="center",
                fontsize=9.5, color=colour,
                fontweight="bold" if c == "positive" else "normal")
        ax.text(len(positions) - 0.25, i, labels[i], ha="left", va="center",
                fontsize=9, color=MUTED)

    ax.text(-0.75, -1.55, "position", ha="right", va="center",
            fontsize=9.5, color=INK, fontweight="bold")
    ax.text(-0.75, -0.78, "weight", ha="right", va="center",
            fontsize=8.5, color=MUTED)

    ax.set_title("Residues at the scored positions, top 16 reliably mapped "
                 "candidates\nfilled = identical to the reference",
                 fontsize=11, color=INK, pad=16, loc="left")
    fig.tight_layout()
    p = OUT / "fig_positions.png"
    fig.savefig(p, dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    df, site, val = load()
    for fn, args in ((fig_ranking, (df,)),
                     (fig_vs_identity, (df, val)),
                     (fig_positions, (df, site))):
        print(f"wrote {fn(*args)}")


if __name__ == "__main__":
    main()
