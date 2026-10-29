# Running this on your own machine

## 1. Check Python

```bash
python3 --version
```

Needs 3.10 or newer.

## 2. Open it in VS Code

Install [VS Code](https://code.visualstudio.com) and the Python extension by
Microsoft, then File -> Open Folder on this project.

`.vscode/` is already set up: the interpreter points at the virtual
environment you make below, and the Run and Debug panel has all six steps
listed in order, so you can run them with a click instead of typing.

Terminal -> New Terminal opens a shell inside the editor for the commands
that follow.

## 3. Virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Your prompt should now start with `(.venv)`.

On Windows the `source` line will not work, because that is a bash command
and the VS Code terminal runs PowerShell. Use the `Activate.ps1` version. If
PowerShell refuses to run it, run this once and answer Y:

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

## 4. Download the inputs

None of these are committed, because none of them are ours to redistribute.

| File | Where | Save as |
|---|---|---|
| 6OMT with substrate and cofactor | `files.rcsb.org/download/5ICE.pdb` | `data/5ICE.pdb` |
| the same enzyme, empty site | `files.rcsb.org/download/5ICG.pdb` | `data/5ICG.pdb` |
| characterised 6OMTs | uniprot.org, search `ec:2.1.1.128`, Download -> FASTA | `data/known_6omt.fasta` |
| plant O-methyltransferases | uniprot.org, search `(o-methyltransferase) AND (taxonomy_id:33090)`, Download -> FASTA | `data/plant_omts.fasta` |
| pET-28a(+) | SnapGene plasmid library or Addgene | anywhere; see step 6 |

The two PDB links download directly if you paste them into a browser.

The README's numbers come from a download of 39,967 sequences that included
unreviewed (TrEMBL) entries. UniProt updates constantly, so your pool will
differ a little and your numbers will too. The scripts do not assume a count.
Adding `AND (reviewed:true)` to the search gives a much smaller pool of
curated sequences, which runs faster but will not reproduce the README's
ranks. Note the date you downloaded if you write anything up.

## 5. Run the analysis

```bash
python3 scripts/active_site.py
python3 scripts/sequences.py
python3 scripts/screen.py
python3 scripts/validate.py
python3 scripts/figures.py
```

Run them in that order; each reads what the previous one wrote.

`validate.py` is the one that matters. It exits non-zero if the controls
fail, which is the signal that the scoring is wrong and the ranking should
Expect `CONTROLS PASS`, with all three characterised 6OMTs
scoring 93 or above and landing in the top few dozen. Several unlabelled
sequences tie with them at the top, so the exact ranks 1 to 19 can shuffle.

## 6. Build the constructs

This step needs pET-28a(+). Open `scripts/constructs.py` and set `VECTOR` to
wherever you saved it:

```python
VECTOR = Path("data/pET-28a.dna")
```

Then:

```bash
python3 scripts/constructs.py
```

Expect four plasmids, each passing 8 of 8 checks. They land in
`output/constructs/` as GenBank files.

Codon optimisation is a real computation and takes a few minutes per
construct. It is not stuck.

## 7. Look at the results

In VS Code, open `output/` in the sidebar. The `.png` figures preview
inline; `screen_ranked.csv` opens as a table.

For the plasmids, load `output/constructs/*.gb` into **Benchling**
(benchling.com, free) or SnapGene Viewer. You should see the circular map
with the T7 promoter, the His tags, your insert and single NdeI and XhoI
sites.

Worth doing by eye at least once: open `data/5ICE.pdb` in the 3D viewer at
`rcsb.org/structure/5ICE`, find His256, Asp257 and the ligand 2H4, and check
that the histidine really is sitting on the substrate hydroxyl. That is the
observation the whole screen is built on.

## 8. If something breaks

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError` | venv not active | re-run the activate line |
| `refusing to continue: positions ... differ` | wrong structure or wrong reference | check you downloaded 5ICE, not another entry |
| `reference Q5C9L7 not present` | the UniProt search returned a different set | make sure `known_6omt.fasta` came from the `ec:2.1.1.128` search |
| `vector not found` | step 6 path not set | edit `VECTOR` in `scripts/constructs.py` |
| `CONTROLS FAIL` | the scoring is not separating the controls | that is the script working; do not use the ranking, find out why |

## 9. Change something

- In `screen.py`, change `W_CATALYTIC`, `W_POLAR` and `W_SHAPE` and re-run
  `validate.py`. Do the controls still pass? How much do the weights have to
  move before they stop?
- Raise `MIN_IDENTITY` from 25% and watch how many candidates become
  unscoreable.
- In `active_site.py`, widen `CONTACT` from 4.5 A and see which extra
  residues join the scored set, and whether that helps or hurts the controls.

Write down what you tried and what happened. That part is yours.
