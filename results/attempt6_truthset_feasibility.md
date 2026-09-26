# Attempt-6 truth-set feasibility (pre-evaluation, 2026-09-26 20:16 IST)
Checked BEFORE locking any truth set or running any enrichment evaluation.

Overlap of the 56 natural panel variants (33 distinct MSA columns) with candidate
"known functional" position sets, all in the 6PTZ reference frame:
- Trp-chain columns (359, 410, 413, 436): 0 variant overlap (fully conserved across panel).
- FAD-contact columns (<=4.5 A, 15 cols): 0 variant overlap.
- ET-pathway-adjacent columns (<=6 A of a chain Trp, 23 cols): 1 variant.

Consequence: a truth set built from chain/FAD-contact POSITIONS makes the locked
enrichment benchmark unwinnable by construction (nothing to recover at the top).
Benchmark design options, to be resolved and locked BEFORE evaluation:
1. Literature variant-level truth set: published mutagenesis substitutions with
   magnetosensitivity/photochemical phenotypes (DmCry, ClCry4, AtCry1), mapped to
   panel columns - overlap to be measured after literature grounding.
2. UniProt FT-annotated functional positions on the reference sequences, mapped to
   panel columns - overlap to be measured.
3. If both overlaps are empty/near-empty: invoke amendment "OR match PLUS a named
   proven plus point" arm, or ChatGPT redirection consult (rule 6) - recorded here
   pre-outcome, not after.

## Update 20:17 IST — UniProt FT route also empty
UniProt functional features on all 4 reference sequences (Q5IZC5 erCry1,
A0A219TAG4, A0A2I4SZI9 erCry4, A0A386QUR4 ClCry4): 20 distinct annotated MSA
columns (FAD binding sites, electron-transfer Trp sites 318/372/395).
Overlap with the 33 variant columns: ZERO. All annotated functional positions
are fully conserved across the 121-row avian panel.
The DmCry mutagenesis literature (e.g. PMC9977682) maps onto the same conserved
core, so route 1 is near-certainly empty too.
Standing observation (candidate scientific finding): natural avian Cry1 variation
is purged from every annotated functional position - purifying selection on the
magnetosensory core. Quantified: 0/33 variant columns overlap 20/20 annotated.
DECISION (pre-evaluation): the locked enrichment benchmark is unwinnable by
construction. Per amendment item 3 + rule 6: ChatGPT redirection consult
(paste route, fresh conversation, own lease) to redesign the attempt-6 win
condition - leading candidate: constraint-signature contrast (migratory-lineage
private variants vs sedentary under the G3 battery) using this purging signal.
To be locked as a new dated amendment BEFORE any evaluation.
