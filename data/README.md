# Inputs

Nothing in this folder is committed. None of it is ours to redistribute, and
the sequence databases change over time, so the honest thing is to fetch your
own copy and record when you did.

Download these four and save them here under exactly these names:

| File | Where |
|---|---|
| `5ICE.pdb` | `https://files.rcsb.org/download/5ICE.pdb` |
| `5ICG.pdb` | `https://files.rcsb.org/download/5ICG.pdb` |
| `known_6omt.fasta` | uniprot.org, search `ec:2.1.1.128`, Download -> FASTA (canonical) |
| `plant_omts.fasta` | uniprot.org, search `(o-methyltransferase) AND (taxonomy_id:33090)`, Download -> FASTA (canonical) |

`scripts/constructs.py` additionally needs pET-28a(+), from SnapGene's plasmid
library or Addgene. Save it here as `pET-28a.dna` (or `.gb`), or point the
`VECTOR` path in that script wherever you put it.

The README's numbers come from 39,967 sequences downloaded on 29 October
2026, which included unreviewed (TrEMBL) entries. UniProt is updated
continuously, so your pool will differ and your numbers will drift from the
published ones. The scripts do not assume a count and the controls are what
decide whether a run is trustworthy. Note the date you downloaded if you write
anything up.

`SETUP.md` in the project root has the full walkthrough.
