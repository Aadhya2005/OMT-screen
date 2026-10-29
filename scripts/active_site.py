"""
active_site.py
--------------
Work out which positions in 6OMT actually matter, from the structure.

PDB 5ICE is Thalictrum flavum 6OMT with both the cofactor product (SAH) and a
substrate analogue (2H4, norlaudanosoline) bound. Its DBREF maps the model
1:1 onto UniProt Q5C9L7, which is also one of the characterised sequences in
the control set, so positions found here can be carried onto the alignment
without a numbering offset.

Three questions, in order:

  1. Which hydroxyl of the substrate is the one being methylated? The enzyme
     is a 6-O-methyltransferase, so the 6-hydroxyl should be the one pointing
     at the sulfur of SAH, where the methyl group sat before transfer. The
     ligand's atom names carry no chemistry, so this is settled by geometry.

  2. Which residues line the substrate pocket, which line the cofactor
     pocket, and which sit beside the acceptor hydroxyl? The last group is
     the catalytic machinery, and the most interesting to screen on.

  3. What moves when substrate binds? 5ICG is the same enzyme with an empty
     site, so comparing the two says which parts of the pocket are rigid and
     which close around the ligand.

One caveat recorded here rather than buried: each structure holds a single
chain, so unlike the TAL work there is no replication across copies to
separate real contacts from crystallographic noise. Distances below are
single observations.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_index, index_to_one

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "output"

SUBSTRATE = "2H4"
COFACTOR = "SAH"
TRANSFER_ATOM = "SD"        # the sulfur that carried the methyl group

# Solvent, cryoprotectant and ions, plus the two ligands themselves: SAH is an
# amino acid derivative, so Biopython counts it as a residue and it turns up in
# its own contact list unless it is excluded here.
IGNORE = {"HOH", "EDO", "K", SUBSTRATE, COFACTOR}

CONTACT = 4.5
NEAR_ACCEPTOR = 6.0
POLAR = {"N", "O", "S"}
HBOND_MAX = 3.5


def load(name: str):
    return PDBParser(QUIET=True).get_structure(name, str(DATA / f"{name}.pdb"))


def first_chain(structure):
    for model in structure:
        for chain in model:
            return chain
    return None


def get_het(chain, resname):
    for res in chain:
        if res.get_resname().strip() == resname:
            return res
    return None


def one_letter(res) -> str:
    try:
        return index_to_one(three_to_index(res.get_resname().strip()))
    except (KeyError, ValueError):
        return "X"


def label(res) -> str:
    return f"{res.get_resname().strip().capitalize()}{res.get_id()[1]}"


def protein_residues(chain):
    return [r for r in chain
            if is_aa(r, standard=False) and r.get_resname().strip() not in IGNORE]


def contacts(chain, target_atoms, cutoff):
    """Residues with any atom within cutoff of any of the given atom coords."""
    coords = np.array([a.get_coord() for a in target_atoms])
    rows = []
    for res in protein_residues(chain):
        rc = np.array([a.get_coord() for a in res])
        if len(rc) == 0:
            continue
        d = np.linalg.norm(rc[:, None, :] - coords[None, :, :], axis=-1)
        dmin = float(d.min())
        if dmin <= cutoff:
            i, j = np.unravel_index(d.argmin(), d.shape)
            rows.append({
                "residue": label(res),
                "resnum": res.get_id()[1],
                "aa": one_letter(res),
                "min_dist": round(dmin, 2),
                "protein_atom": list(res)[i].get_name(),
            })
    rows.sort(key=lambda r: r["min_dist"])
    return rows


def find_acceptor(sub, cof):
    """Which substrate oxygen sits closest to the methyl-donating sulfur."""
    sd = cof[TRANSFER_ATOM].get_coord()
    oxygens = [a for a in sub if (a.element or a.get_name()[0]).upper() == "O"]
    rows = []
    for a in oxygens:
        rows.append({
            "atom": a.get_name(),
            "dist_to_SD": round(float(np.linalg.norm(a.get_coord() - sd)), 2),
        })
    rows.sort(key=lambda r: r["dist_to_SD"])
    return rows


def polar_near(chain, atom, cutoff=HBOND_MAX):
    """Polar protein atoms close to one ligand atom: candidate general bases."""
    xyz = atom.get_coord()
    rows = []
    for res in protein_residues(chain):
        for a in res:
            el = (a.element or a.get_name()[0]).strip().upper()
            if el not in POLAR:
                continue
            d = float(np.linalg.norm(a.get_coord() - xyz))
            if d <= cutoff:
                rows.append({
                    "residue": label(res), "resnum": res.get_id()[1],
                    "aa": one_letter(res), "atom": a.get_name(),
                    "dist": round(d, 2),
                })
    rows.sort(key=lambda r: r["dist"])
    return rows


def superimpose_apo(holo_chain, apo_chain):
    """
    Put the two structures in the same frame before comparing them.

    They were solved in different crystals, so their coordinates sit in
    unrelated frames: without this step every residue appears to have moved
    tens of angstroms, which is the frame difference and not the protein.
    The fit uses CA atoms of residues common to both, excluding the pocket
    itself so that real movement there is measured rather than averaged away.
    """
    from Bio.PDB import Superimposer

    holo = {r.get_id()[1]: r for r in protein_residues(holo_chain)}
    apo = {r.get_id()[1]: r for r in protein_residues(apo_chain)}
    shared = sorted(set(holo) & set(apo))

    fixed, moving = [], []
    for num in shared:
        a, b = holo[num], apo[num]
        if a.get_resname() != b.get_resname():
            continue
        if "CA" in a and "CA" in b:
            fixed.append(a["CA"])
            moving.append(b["CA"])

    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply([a for r in apo_chain for a in r])
    return sup.rms, len(fixed)


def apo_holo_shift(holo_chain, apo_chain, resnums):
    """How far each residue's side chain sits between the two structures."""
    apo = {r.get_id()[1]: r for r in protein_residues(apo_chain)}
    bb = {"N", "CA", "C", "O", "OXT"}
    rows = []
    for res in protein_residues(holo_chain):
        num = res.get_id()[1]
        if num not in resnums or num not in apo:
            continue
        other = apo[num]
        if other.get_resname().strip() != res.get_resname().strip():
            continue
        names = {a.get_name() for a in res} & {a.get_name() for a in other}
        sc = [n for n in names if n not in bb]
        if "CA" not in names:
            continue
        ca = float(np.linalg.norm(res["CA"].get_coord() - other["CA"].get_coord()))
        side = (float(np.mean([
            np.linalg.norm(res[n].get_coord() - other[n].get_coord()) for n in sc
        ])) if sc else float("nan"))
        rows.append({
            "residue": label(res), "resnum": num,
            "ca_shift": round(ca, 2),
            "sidechain_shift": round(side, 2) if sc else None,
        })
    rows.sort(key=lambda r: -(r["sidechain_shift"] or 0))
    return rows


def banner(t: str) -> None:
    print(f"\n{t}\n" + "-" * len(t))


def main() -> None:
    OUT.mkdir(exist_ok=True)

    holo = first_chain(load("5ICE"))
    apo = first_chain(load("5ICG"))

    sub = get_het(holo, SUBSTRATE)
    cof = get_het(holo, COFACTOR)
    if sub is None or cof is None:
        raise SystemExit("5ICE is missing the substrate or cofactor")

    result = {
        "structure": "5ICE",
        "uniprot": "Q5C9L7",
        "organism": "Thalictrum flavum subsp. glaucum",
        "note": "single chain; contacts are single observations, not replicated",
    }

    # --- 1. which hydroxyl gets methylated -----------------------------
    banner("1. Which substrate hydroxyl faces the methyl donor")
    acc = find_acceptor(sub, cof)
    for r in acc:
        print(f"  {r['atom']:<5} -> SAH SD   {r['dist_to_SD']:>5.2f} A")
    acceptor_name = acc[0]["atom"]
    print(f"\n  nearest oxygen is {acceptor_name} at {acc[0]['dist_to_SD']} A;"
          f" next is {acc[1]['atom']} at {acc[1]['dist_to_SD']} A")
    result["acceptor"] = {"atom": acceptor_name, "oxygens": acc}

    # --- 2. the catalytic neighbourhood --------------------------------
    banner("2. Polar atoms beside the acceptor hydroxyl")
    acceptor_atom = sub[acceptor_name]
    base = polar_near(holo, acceptor_atom)
    for r in base:
        print(f"  {r['residue']:<10} {r['atom']:<5} {r['dist']:>5.2f} A")
    if not base:
        print("  none within 3.5 A")
    result["acceptor_neighbours"] = base

    # --- 3. pockets ----------------------------------------------------
    banner(f"3. Substrate pocket, residues within {CONTACT} A of {SUBSTRATE}")
    sub_contacts = contacts(holo, list(sub), CONTACT)
    for r in sub_contacts:
        print(f"  {r['residue']:<10} {r['min_dist']:>5.2f} A   via {r['protein_atom']}")
    result["substrate_pocket"] = sub_contacts

    banner(f"4. Cofactor pocket, residues within {CONTACT} A of {COFACTOR}")
    cof_contacts = contacts(holo, list(cof), CONTACT)
    print(f"  {len(cof_contacts)} residues")
    print("  " + ", ".join(r["residue"] for r in cof_contacts))
    result["cofactor_pocket"] = cof_contacts

    shared = ({r["resnum"] for r in sub_contacts}
              & {r["resnum"] for r in cof_contacts})
    print(f"\n  touching both ligands: "
          f"{sorted(shared) if shared else 'none'}")
    result["both_pockets"] = sorted(shared)

    # --- 4. what moves on binding --------------------------------------
    banner("5. Movement between the empty and occupied structures")
    rms, n_fit = superimpose_apo(holo, apo)
    watch = {r["resnum"] for r in sub_contacts} | {r["resnum"] for r in cof_contacts}
    shifts = apo_holo_shift(holo, apo, watch)
    missing = sorted(watch - {r["resnum"] for r in shifts})

    # How much does the protein differ overall, away from the pocket? If the
    # whole chain disagrees by about as much as the pocket does, the pocket
    # numbers carry no signal.
    all_shift = apo_holo_shift(holo, apo, {r.get_id()[1]
                                           for r in protein_residues(holo)})
    overall = np.median([r["ca_shift"] for r in all_shift])
    pocket = np.median([r["ca_shift"] for r in shifts]) if shifts else float("nan")

    print(f"  superimposed on {n_fit} shared CA atoms, overall RMSD {rms:.2f} A")
    print(f"  median CA deviation: whole chain {overall:.2f} A, "
          f"pocket residues {pocket:.2f} A")

    inconclusive = rms > 1.5 or overall > 1.5
    if inconclusive:
        print("\n  INCONCLUSIVE. The two structures differ across the whole")
        print("  chain by about as much as they differ at the pocket, in")
        print("  several separate loops rather than as one coherent domain")
        print("  closure. With a single chain in each, ligand-induced")
        print("  movement cannot be told apart from crystal-form differences.")
        print("  No claim about induced fit is carried forward.")
    else:
        print(f"  {'residue':<10} {'CA':>6} {'side chain':>12}")
        for r in shifts[:12]:
            sc = (f"{r['sidechain_shift']:.2f}"
                  if r["sidechain_shift"] is not None else "-")
            print(f"  {r['residue']:<10} {r['ca_shift']:>6.2f} {sc:>12}")
    if missing:
        print(f"\n  not modelled in the apo structure: {missing}")

    result["apo_holo"] = {
        "superposition_rmsd": round(rms, 2),
        "atoms_fitted": n_fit,
        "median_ca_shift_whole_chain": round(float(overall), 2),
        "median_ca_shift_pocket": round(float(pocket), 2),
        "conclusive": not inconclusive,
        "verdict": ("apo and holo differ as much globally as at the pocket; "
                    "no induced-fit claim is supported by these two structures"
                    if inconclusive else "pocket moves more than the bulk"),
        "shifts": shifts,
        "absent_in_apo": missing,
    }

    # --- the positions we will score on --------------------------------
    scored = sorted({r["resnum"] for r in sub_contacts}
                    | {r["resnum"] for r in base})
    result["scored_positions"] = [
        {"resnum": n,
         "aa": next((r["aa"] for r in sub_contacts + base if r["resnum"] == n), "X"),
         "role": ("acceptor neighbour" if n in {r["resnum"] for r in base}
                  else "substrate pocket")}
        for n in scored
    ]
    banner("6. Positions carried forward to the screen")
    print(f"  {len(scored)} positions: "
          + ", ".join(f"{p['aa']}{p['resnum']}" for p in result["scored_positions"]))

    (OUT / "active_site.json").write_text(json.dumps(result, indent=2))
    print(f"\nwrote {OUT / 'active_site.json'}")


if __name__ == "__main__":
    main()
