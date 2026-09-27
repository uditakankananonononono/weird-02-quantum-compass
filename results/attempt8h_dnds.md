# W02 queue #11 — Episodic-selection screen on avian Cry4 (10-species 5v5)

**Amendment:** 13:12 IST + addendum 37010e4 (13:15 IST), ref `8c424cb`. Delegated executor, **report-only, UNCOMMITTED**, returned to lane-29 for review/commit.
**Interpretation lock (verbatim):** report-only screen; a significant elevation would be a screen signal for discussion, never an adaptive-evolution claim; no causal residues, no migratory-prediction claim.

## Headline result — NEGATIVE screen

| Contrast | median omega | mean omega | n pairs |
|---|---|---|---|
| cross-class (migratory x sedentary) | 0.111674 | 0.121736 | 25 |
| within-class | 0.136581 | 0.130531 | 20 |

**Mann-Whitney U = 211.0, two-sided p = 0.3792 — no screen signal.** Cross-class dN/dS is, if anything, *lower* than within-class (opposite of the elevation direction). All omega values sit in the 0.079-0.190 range typical of purifying selection on Cry4. Honest negative, reported as such.

## Pavo cristatus exclusion (disclosed verbatim per addendum)
Pavo cristatus is excluded from this screen: no committed label exists anywhere in the repo (absent from results/h1_scoring.json, n_partial=0) and no protein-MSA row exists in cry1_panel_msa.fasta. This exclusion is disclosed here and must appear in the paper.

## DEVIATION requiring lane-29 review/veto — back-map target MSA
The amendment text says to back-map onto **cry1_panel_msa.fasta**. The ARM-2 identity check (below) shows all 10 CDS translations match the **committed Cry4 subpanel MSA (data/mafft_ebi_crosscheck.fasta, accession-keyed) at 0.994-1.000**, and match Cry1 panel rows at only 0.03-0.10 (different accessions, Cry1-family proteins). Back-mapping Cry4 CDS onto Cry1 rows would be invalid, so the screen was executed against the Cry4 subpanel MSA. This deviation is flagged for lane-29 to accept or veto.

Identity check (CDS translation vs each MSA):

| species | vs Cry4 subpanel MSA | vs Cry1 panel row | CDS aa |
|---|---|---|---|
| Columba_livia | 0.9943 | 0.0984 | 526 |
| Erithacus_rubecula | 0.9982 | 0.0774 | 557 |
| Gallus_gallus | 0.9981 | 0.1014 | 530 |
| Passer_domesticus | 1.0000 | 0.0274 | 359* |
| Catharus_ustulatus | 0.9981 | 0.0952 | 528 |
| Ficedula_albicollis | 0.9981 | 0.1032 | 528 |
| Hirundo_rustica | 0.9981 | 0.0968 | 528 |
| Serinus_canaria | 0.9981 | 0.0987 | 528 |
| Struthio_camelus | 0.9981 | 0.1047 | 530 |
| Zonotrichia_albicollis | 0.9981 | 0.0952 | 528 |

## ARM 1 — CDS retrieval ledger (10/10, zero failures; sha256 per accession in results/attempt8h_cds_ledger.json)
- 4 REFSEQ nuccore efetch fasta_cds_na: KX168611 (Columba), MN709784 (Erithacus), NM_001039596 (Gallus), AY494987 (Passer)
- 3 XP protein efetch fasta_cds_na: XP_032936297.1 (Catharus), XP_039941661.1 (Hirundo), XP_068774195.1 (Struthio)
- 3 ENS via Ensembl REST: lookup/id (Translation -> Parent transcript) then sequence/id(transcript)?type=cds: ENSFALP00000002532.1 (Ficedula, via ENSFALT00000002544), ENSSCAP00000009299.1 (Serinus, via ENSSCAT00000010504), ENSZALP00000018676.1 (Zonotrichia, via ENSZALT00000024665). Direct sequence/id on the Translation id silently returns protein; alphabet was verified nucleotide for all three.
- Labels sourced ONLY from results/h1_scoring.json rows: 5 migratory (Ficedula, Erithacus, Zonotrichia, Catharus, Hirundo), 5 sedentary (Passer, Serinus, Columba, Struthio, Gallus).

## ARM 2 — back-map
CDS back-mapped onto the accession-keyed Cry4 subpanel protein MSA rows; N-containing codons gap-masked (see caveats); all 10 alignments equal length.

## ARM 3 — pairwise Nei-Gojobori dN/dS + Mann-Whitney (locked fallback)
PAML/HyPhy were NOT available in this environment (`which yn00 codeml paml hyphy HYPHYMP` all empty), so the amendment's locked fallback was used: pairwise Nei-Gojobori (Biopython 1.88 Bio.codonalign.codonseq.cal_dn_ds, method NG86) over all 45 pairs + two-sided MWU. This is a heuristic screen, not a branch-site test.

## Caveats (compact limitations, per regime)
1. **Power:** 5v5 classes; lane-29's floor at label level is min two-sided p = 0.00794. Weakly powered; a real effect could be missed. This negative does not rule out episodic selection.
2. **Pairwise non-independence:** the 45 pair values are non-independent (each species enters 9 pairs); MWU treats them as independent. Screen-level heuristic, exactly as pre-registered.
3. **Passer fragment:** AY494987 covers 359 of ~528 aa; pairs involving Passer use only the covered columns (identity 1.000 over the fragment).
4. **KX168611 ambiguity:** 2 IUPAC codes (R,Y) masked to N, the affected codons gap-masked; only pairs with Columba skip those columns.
5. **Pairwise NG, not a branch model:** cannot localize episodes to lineages or residues; no causal-residue or migratory-prediction claim is made.

## Reproducibility
Determinism verified: two consecutive runs produce identical values (45 pairs, U=211.0, p=0.37918336397400576). Compute script: src/attempt8h_dnds.py (reads data/cds/*.fna pinned by the sha256 ledger). CDS bytes re-derivable via the retrieval routes above.
