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
