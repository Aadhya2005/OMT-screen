# Screening plant O-methyltransferases for a 6OMT that works in E. coli

A structure-guided screen of 39,967 plant O-methyltransferase sequences from
UniProt, asking which look like (S)-norcoclaurine 6-O-methyltransferases and
which would be worth trying in *E. coli*. All 39,967 are scored, but only
13,987 of them map onto the reference structure well enough for the score to
be trusted; the other 25,980 are scored and flagged. Ends with four verified
expression constructs.

**Computational only. Nothing here has been tested at the bench. The ranking
says which sequences are worth testing; it does not say that any of them have
activity.**

## Why this enzyme

Benzylisoquinoline alkaloids (morphine, codeine, berberine) are made in
plants and can be produced instead by rebuilding the plant pathway inside
*E. coli*. The published route runs glucose to L-tyrosine to dopamine to
(S)-norcoclaurine, then three methylation steps to (S)-reticuline, which is
the branch point to everything downstream.

6OMT is the first of those methylations, and Robin et al. put "rate-limiting"
in the title of the structure paper. It is also a good test case for
screening, because plant O-methyltransferases are a large, promiscuous
family: many enzymes look like 6OMT and do something else entirely.

## What the screen does

| Stage | Script | Output |
|---|---|---|
| Find the functional positions | `active_site.py` | which residues matter, and why |
| Tie sequences to the structure | `sequences.py` | reference check, position mapping |
| Score and filter the candidates | `screen.py` | `screen.csv` |
| Test the ranking against controls | `validate.py` | `screen_ranked.csv`, pass or fail |
| Draw the result | `figures.py` | three figures |
| Build constructs for the shortlist | `constructs.py` | four verified plasmids |

Run them in that order.

## What came out

Numbers below are from one run on 29 October 2026, on 39,967 sequences
(UniProt, Viridiplantae, Swiss-Prot and TrEMBL together). UniProt changes
constantly, so a fresh download will give slightly different numbers.

**The structure identifies its own chemistry.** The substrate analogue has
four hydroxyls. The one 3.40 A from the cofactor's sulfur is the one being
methylated; the next nearest is 4.59 A away. Sitting on that hydroxyl are
His256 at 2.99 A and Asp257 at 2.64 A, which is the His/Asp general-base
pair that plant O-methyltransferases use. None of that was assumed: it falls
out of distance measurements.

**The controls pass.** The three characterised 6OMTs score 100.0, 93.2 and
93.2 and land at ranks 1, 7 and 19 out of 39,967. The 74 near-miss sequences
(ipecac and related enzymes, which carry the same EC number but methylate a
different position of a different alkaloid) score between 29.5 and 65.9. The
lowest positive scores 93.2, so the two groups do not overlap.

**Ranks 7 and 19 are partly a tie effect.** The top 23 sequences all score
93.2 or higher, and several of them are unlabelled Papaver, Thalictrum and
Cinnamomum sequences that are probably 6OMTs nobody has characterised. The
score cannot separate things inside that group, so the order within it is
decided by sort order and the difference between rank 7 and rank 19 means
nothing.

**It is not a similarity search in disguise.** Spearman correlation between
the active-site score and plain sequence identity is 0.49 over the 13,987
candidates whose mapping is reliable. Had that come out above 0.9, the honest
conclusion would have been that the structural work added nothing, and
`validate.py` fails the run on exactly that condition. Ranked by identity
alone the three positives would sit at 1, 5 and 28.

**The catalytic pair filters hard.** Only 8,898 of the 39,967 candidates keep
both His256 and Asp257 equivalents. All three positives do.

## Four things done deliberately

**The reference is checked against the structure before anything else.**
Position numbering is how this kind of analysis goes wrong without anyone
noticing. The model in 5ICE is compared residue by residue against UniProt
Q5C9L7 at every scored position, and `sequences.py` refuses to run on a
mismatch. They agree at all 17 positions, 99.7% over the whole overlap.

**Two positions are excluded from scoring on purpose.** Gly165 and Cys253
touch the substrate through their backbone carbonyl, not their side chains.
Any amino acid there presents the same atom, so scoring the residue identity
would measure nothing. The same logic demoted Thr113 and Thr311: both reach
the substrate with a carbon and not their hydroxyl, so they count as pocket
lining. Classification follows the contacting atom, not the residue type.

**Each candidate is aligned to the reference on its own.** The usual
approach is one multiple alignment of everything, but then a badly aligned
sequence corrupts shared columns invisibly. Aligned one at a time, every
candidate carries its own identity and coverage, so an untrustworthy mapping
gets flagged instead of being scored as though it were sound. 25,980 of the
39,967 (65%) fall below the thresholds (25% identity, 60% coverage) and are
marked. That is a lot, and it is expected: the pool includes fragments and
distant members of the family that share little with a 6OMT.

**The apo-versus-holo comparison was abandoned.** 5ICG is the same enzyme
with an empty active site, and comparing the two should say what closes on
binding. It does not: the structures differ across the whole chain by about
as much as they differ at the pocket, and the differences sit in several
separate loops instead of adding up to one coherent motion. With a single
chain in each there is no way to separate ligand-induced movement from
crystal-form differences, so `active_site.py` reports the comparison as
inconclusive and carries no claim forward.

## Two bugs this work exposed

Both were in the construct-design pipeline carried over from the earlier
RsTAL project, and neither was visible when that pipeline had only ever been
run on one protein.

**A latent crash in the optimiser.** Two of the hairpin constraints trip a
bug in DNA Chisel where localising a constraint returns `None`, which the
solver then dereferences. Whether it fires depends on the sequence. RsTAL
never triggered it; two of these four 6OMTs do. The fix tries the full
constraint set first and falls back to a reduced one only on failure. It
also prints which constraints were dropped, so a weakened run shows up in
the output instead of passing unnoticed.

**A primer check that flagged the wrong thing.** Every reverse primer came
back with a self-dimer around -7.1 kcal/mol, identical across four unrelated
proteins. The cause is that both cloning sites are palindromes: `CTCGAG` is
its own reverse complement, so any primer carrying it pairs with a copy of
itself. That is inherent to restriction cloning and no clamp removes it. The
dimer's melting temperature is about 5 C, so at a 60 C anneal it does not
exist. The check now asks whether a dimer is stable near the annealing
temperature, not merely whether it is energetically strong.

That second one is the mirror of a lesson from the first project, where
hairpin warnings had to be moved off Tm and onto dG. Neither number is
enough on its own. The question is whether the structure holds together at
the temperature the reaction actually runs at.

## The constructs

Four of them, and the number is small on purpose: a small set keeps the
experiment interpretable. Each is codon-optimised for *E. coli*, cleared of
the cloning sites, given primers, assembled into pET-28a(+) and verified by
the same eight checks used in the earlier project.

| Accession | Organism | Rank | Score | Why |
|---|---|---|---|---|
| Q9LEL6 | *Coptis japonica* | 7 | 93.2 | the enzyme the published platform uses |
| Q5C9L7 | *Thalictrum flavum* | 1 | 100.0 | the sequence 5ICE was solved from |
| Q6WUC1 | *Papaver somniferum* | 19 | 93.2 | a third characterised enzyme |
| A0A0B4VG62 | *Sinopodophyllum hexandrum* | 31 | 88.6 | annotated as a 6OMT, never characterised |

Three knowns and one unknown. If the three express and the fourth does not,
the screen was wrong about it, which is a result worth having.

All four pass 8 of 8 verification checks. The two NdeI forward primers for
Q9LEL6 and Q5C9L7 carry a low-GC warning (37 to 39% against a 40% floor).
That comes from the deliberately AT-rich start of the gene, and the Q9LEL6
primer also has a run of five A's. Neither is a design fault, but both are
worth knowing before ordering.

The GenBank files open in SnapGene, Benchling or any sequence viewer.

**How the fourth construct was chosen.** It was picked during development on
a smaller, reviewed-only subset of this data, where it ranked 4th. On the
full 39,967 it ranks 31st, and the validation output lists twelve
higher-ranked candidates that are not known 6OMTs (for example J7MCT6 at rank
3 and A0A2G5CX73 at rank 6). A fresh pick from the full-pool shortlist would
be a different, and probably better, choice. It was kept as built.

## Running it

Requires Python 3.10 or newer.

```bash
pip install -r requirements.txt

python3 scripts/active_site.py
python3 scripts/sequences.py
python3 scripts/screen.py
python3 scripts/validate.py
python3 scripts/figures.py
python3 scripts/constructs.py      # needs pET-28a(+), see SETUP.md
```

Inputs are not committed. `SETUP.md` says where to get each of them.

## Limitations

Sequence similarity predicts function weakly. Everything here is a
hypothesis about what is worth testing.

The pool is unfiltered. It includes unreviewed entries, fragments and
sequences labelled only "uncharacterized protein". In the printed control
table the class "putative" just means a sequence is not in either control
set, so it says nothing about what the enzyme does.

The expressibility filters are crude heuristics: a hydrophobicity window
standing in for a membrane-segment predictor, a cysteine count standing in
for disulfide prediction. They flag candidates for a second look. They are
nowhere near good enough to rule anything out, and a real version would use
dedicated predictors.

The active site was mapped from a single chain of a single structure, so
unlike the earlier TAL work there is no replication across copies to
separate real contacts from crystallographic noise.

Only one family of near misses is tested (the ipecac enzymes). The score has
not been checked against other close relatives such as 4'OMTs, so whether it
can separate those is untested here. Settling which ring position gets
methylated would need either many more characterised enzymes with known
regiochemistry, or modelling the substrate in each candidate instead of
reading residue identities off an alignment.

## References

- Minami H et al. Microbial production of plant benzylisoquinoline alkaloids. *PNAS* 2008;105:7393
- Nakagawa A et al. A bacterial platform for fermentative production of plant alkaloids. *Nat Commun* 2011;2:326
- Robin AY et al. Crystal structure of norcoclaurine-6-O-methyltransferase, a key rate-limiting step in the synthesis of benzylisoquinoline alkaloids. *Plant J* 2016;87:641
- Sequences from UniProt; coordinates from the RCSB Protein Data Bank

## How this was built

The code was written with AI assistance (Claude) under my direction. I chose
the target enzyme, selected the structures and the candidate set, directed
the analysis and checked the output. Every number in this README came out of
these scripts; none of them were copied from anywhere else.
