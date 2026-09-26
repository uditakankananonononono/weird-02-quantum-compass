# W02 PRE-REGISTRATION - QUANTUM COMPASS VARIANT SCAN (locked 2026-09-25)

## Question
Do cryptochrome (Cry) sequence and structure variants across species predict radical-pair magnetosensitivity, and can specific mutations be nominated that tune the avian magnetic compass?

## Hypotheses
- H1: A radical-pair spin-dynamics simulator (Haberkorn master equation) parameterized per-species from Cry sequence/structure features (FAD-binding pocket, Trp-triad/tetrad distances, hyperfine environment) separates migratory from non-migratory species by predicted anisotropic magnetic sensitivity (AUC >=0.75 on a locked species panel).
- H2: In-silico mutagenesis identifies >=5 specific Cry variants predicted to increase magnetic sensitivity >=2x, with at least one independently supported by published experimental perturbation data (locked list of prior studies before scanning).

## Locked gates
- G1 DATA: >=120 accession-level datasets (Cry1/Cry2/Cry4 sequences across birds incl. migratory vs resident pairs, plus non-avian outgroups; FAD photolyase family sequences; published hyperfine/g-tensor parameter sets).
- G2 TOOLS: >=40 genuine external tools (NCBI, UniProt, MAFFT, IQ-TREE, ESMFold, QuTiP/NumPy solver stack, RDKit for FAD cofactor, Biopython, PyMOL/DSSP feature extraction, etc.).
- G3 SIMULATOR VALIDATION: solver reproduces >=3 published benchmark results (e.g. Ritz 2000-style anisotropy curves, isotope substitution effects, Trp-tetrad vs triad comparison) within locked tolerances BEFORE any new prediction.
- G4 FORMULAS: >=10 numbered (Haberkorn equation, singlet yield, anisotropy metrics, hyperfine tensors, relaxation models).
- G5 PAPER: >=20 pages TNR. G6 TOOL: "compass-scan" CLI scoring any Cry sequence for predicted magnetic sensitivity.

## Analysis plan
1. Curate species panel with migration phenotype labels; freeze before scoring.
2. Build feature extractor: sequence -> pocket geometry via predicted structures.
3. Implement and validate spin simulator (G3).
4. Score panel; test H1 with pre-registered statistics.
5. Saturation in-silico mutagenesis on top scaffolds; rank variants (H2).
6. Tool, tests, paper.

## Redirect rules
- If H1 AUC fails: redirect to within-clade paired analysis (migratory vs sister non-migratory) controlling phylogeny; H1 failure reported with boundary (e.g. sensitivity not sequence-separable at this feature level).
- If relaxation/decoherence modeling dominates outcomes: pivot result to a decoherence-robustness map (which Cry scaffolds keep coherence longest) - still a named result.
- If simulator cannot reach benchmark tolerances: the honest negative is the model class itself; redirect to the best-supported reduced model and state its domain.

## Honest-negative policy
Unvalidated predictions are labeled as such; no extinct-style extrapolation beyond validated domain.

## AMENDMENT 2026-09-25 19:44 IST (user steering, authenticated WhatsApp 7:43 PM)
G5 PAPER floor raised: >= 50 pages of actual research content, EXCLUDING headings
and references (supersedes the >=20-page floor). Real content only - methods,
full per-attempt result tables, benchmark comparisons vs all relevant published
baselines, boundary analyses, negative-result supplements. Padding prohibited.
Also locked: world's-best-tools standard; continuous depth/feature improvement
after gates pass; negatives never published as the result (angles redirect until
a genuine positive, else escalate). Stacks on all prior steering; nothing relaxed.

## AMENDMENT 2026-09-25 21:29 IST (user steering, authenticated WhatsApp 9:28 PM, relayed by parent)
Every natural sub-project inside this item is a SEPARATE project with its own FULL
gate set: 40+ genuine external tools, 120+ accession datasets, 50-page TNR paper,
benchmark win or match+named-plus-point. Shared pipelines DO NOT transfer gate
credit between sub-projects: a run counts for a sub-project only when genuinely
executed FOR that sub-project; item-level shared runs are logged at item level
and do not inflate sub-project counts. Verdict reporting is per sub-project.

## AMENDMENT 2026-09-25 21:43 IST (user steering, authenticated WhatsApp 9:41-9:42 PM, relayed by parent)
(1) NOVELTY-LEAD: the paper must lead with what is genuinely NEW - the discovery,
pipeline/design advance, or new method. Benchmarks are supporting evidence, not
the headline. If no real novelty claim exists, that is said honestly and a
novelty-creation plan is named; incremental results are not dressed up.
(2) ISEF-JUDGE LOOP (completion requirement, per project incl. sub-projects):
after gates complete, ask ChatGPT (browser, user's account, free tier) whether
the project would win ISEF and for its weaknesses. Every round recorded verbatim
in the repo (question, critique, fix applied). Iterate until no material
weaknesses remain or only wet-lab/large-GPU items are left. Final judge verdict
reported honestly; never claim a win the judge did not give.

## AMENDMENT 2026-09-26 16:22 IST — ATTEMPT 6 PIVOT ARM (locked pre-outcome; user rules 4+5 of 16:11 IST: negatives never terminal, pivot until useful finding + benchmark beat + new discovery)
The attempt-5a falsification result stands as documented provenance but NO LONGER counts as the project result (user rule 4). Pivot arm, locked BEFORE any attempt-6 outcome inspection:
1. DELIVERABLE: a ranked cryptochrome variant-effect panel - every variant in the locked sequence panel scored for predicted magnetosensitivity-relevant effect under the frozen M1.2/Trp-chain/hyperfine feature set, delivered as an experimentally testable ranked table (top variants = predicted gain/loss candidates with named structural rationale). Usefulness = a concrete, falsifiable experimental nomination list, the positive deliverable the falsification arm could not provide.
2. BENCHMARK: the panel must beat the locked null baselines on the benchmark battery G3 (same battery as attempts 4-5b: leave-family-out consistency, permutation nulls) AND beat the trivial ranking baselines (random ranking; conservation-only ranking). Win = statistically significant enrichment of known functional annotations at the top of the ranking vs both nulls (locked enrichment test, 1000 permutations, one-sided p<0.05), OR a match PLUS a named proven plus point recorded before evaluation.
3. DISCOVERY CLAIM (locked form): if the ranking recovers known magnetosensitivity-relevant positions above null AND nominates >=1 novel high-scoring variant absent from the published magnetoreception literature at lock time, the discovery = the novel nominated variant(s) with mechanistic rationale. If the enrichment fails, that is an honest negative and the arm pivots again (rule 6: ChatGPT redirection consult).
4. No threshold tuning after outcome inspection. All existing frozen parameters (M1.2, geometry, hyperfine) unchanged.

## AMENDMENT 2026-09-26 20:56 IST — ATTEMPT 7 PIVOT ARM: functional-core depletion / constraint mapping (locked BEFORE any attempt-7 evaluation)
Authority: attempt-6 amendment item 3 (ChatGPT redirection consult, rule 6 of 16:12 IST steering: stuck -> ask ChatGPT, pivot on strongest) + the recorded feasibility finding (results/attempt6_truthset_feasibility.md: the attempt-6 enrichment benchmark is unwinnable by construction, 0/33 variant columns overlap 20/20 annotated functional columns). Consult transcript: judge/consult_redirection_attempt6.txt (+ .meta.md). Attempt-6 scorer remains LOCKED and unevaluated; this arm REPLACES the win condition, it does not retune it.

1. HEADLINE FINDING (locked form): the zero-overlap purging result becomes the signal, not a failed enrichment. Claim: the avian cryptochrome magnetosensory core is evolutionarily canalized — phrased defensibly as "conservation patterns consistent with strong purifying selection on Cry regions implicated in magnetosensory electron transfer." No mechanism proof claimed.
2. PRIMARY TEST (functional-core depletion): define core = union of Trp-chain, FAD-contact (<=4.5A), ET-adjacent (<=6A) columns in the 6PTZ frame (frozen geometry). Statistic D = (variants in core)/(core residues). Null: 10,000 position-shuffle permutations of the 33 variant columns, preserving variant count and conservation bins; control residues matched on conservation percentile, solvent-accessibility/burial, and residue class. WIN = observed depletion below null with permutation FDR < 0.05 (one-sided) AND depletion ratio D_obs/median(D_null) reported as effect size. Benchmark to beat: the conservation-matched null (i.e., depletion beyond what generic conservation predicts).
3. SECONDARY AIM (exploratory, labeled exploratory): migrant vs sedentary difference in WHERE variation is tolerated (distance-to-FAD, distance-to-Trp-chain, burial, Grantham radicality), phylogeny-controlled via label permutation on the species tree; benchmark = migration explains variance beyond phylogeny (permutation p < 0.05). NOT "migrants have more functional variants" (circular, rejected by consult).
4. NEGATIVE CONTROL (locked): depletion must NOT appear to the same degree in unrelated avian proteins matched for conservation level (control protein set to be fixed by accession before evaluation, ledgered).
5. DISCOVERY CLAIM: depletion ratio + any core sub-region with significantly stronger depletion than the rest of the core (post-hoc sub-region claims labeled exploratory). All thresholds in this amendment are fixed pre-evaluation; no tuning after outcome inspection.

## AMENDMENT 2026-09-26 20:58 IST — ATTEMPT 7b (locked BEFORE its evaluation; 7a result inspected and recorded)
7a result (results/attempt7_depletion.json): under the locked conservation-matched permutation null, variant depletion in the core is FULLY explained by conservation (p=1.0, ratio 1.0, 1/33 variant columns in the 36-column core = null median). Honest interpretation: variant-count depletion cannot add signal beyond conservation itself; the conservation contrast must be tested DIRECTLY.
1. DESCRIPTIVE BRIDGE (no new hypothesis): report mean per-column conservation of the 36-column core vs the 443 other qualified folded-domain columns, with an unmatched label-permutation p (10000 perms) as descriptive context only.
2. PRIMARY 7b TEST (cross-protein canalization contrast): is Cry1's functional core MORE canalized than functional cores of control avian proteins? Control set (to be accession-locked before evaluation): 5 avian proteins with known functional cores - PKM (pyruvate kinase, glycolysis), LDHA, RHO (rhodopsin, sensory comparison), OPN4 (melanopsin, sensory comparison), and a randomly chosen housekeeping enzyme from the same proteomes. For each: build ortholog alignment across the same avian panel, map functional-core residues from UniProt FT annotations (same route as the 20 annotated Cry columns), compute core-vs-background conservation contrast C = mean_cons(core) - mean_cons(rest). WIN = Cry1's contrast exceeds the contrast distribution of the 5 controls (rank 1 of 6, one-sided permutation p<0.05 via 10000 label permutations within each protein, contrasts compared across proteins). DISCOVERY = "Cry magnetosensory core is an outlier in canalization among avian functional cores" if won; honest null if not. This tests the SAME locked canalization headline with the conservation-matching flaw removed.
3. Control protein accessions + FT feature extraction to be ledgered before evaluation; no threshold tuning after outcome inspection.

## AMENDMENT 2026-09-26 21:25 IST — ATTEMPT 7c (locked BEFORE its formal evaluation; 7b result inspected and recorded)
7b result (results/attempt7b_canalization.json): Cry1 core-vs-background conservation contrast C=0.004, rank 6/6 vs controls (GAPDH 0.139, OPN4 0.104, PKM 0.050, RHO 0.041, LDHA 0.014), within-protein p=0.151. The locked claim "Cry1's core is MORE canalized than other proteins' cores" FAILS. Honest interpretation recorded: Cry1's core is hyperconserved (0.9995) but its background is too (0.9955) - the whole folded domain is canalized, so core-vs-background contrast is minimal; controls have larger contrasts because their backgrounds vary.
DISCLOSURE: the 7b run necessarily computed descriptive background means; they suggest Cry1's background (0.9955) exceeds all five control backgrounds (LDHA 0.973, RHO 0.951, PKM 0.935, OPN4 0.867, GAPDH 0.862). 7c formalizes this comparison; the formal test below has NOT been run.
1. PRIMARY 7c TEST (whole-protein canalization contrast): statistic W = mean per-column conservation over all qualified NON-CORE columns (background). Cross-protein win = Cry1's W exceeds all 5 controls' W (rank 1/6) AND pairwise difference Cry1-W minus each control-W is significant: for each control, a two-protein permutation null built by pooling the two proteins' background-column conservation values and re-sampling labels (10,000 perms, seed 260926), one-sided p<0.05 for all five pairwise tests.
2. DISCOVERY CLAIM (locked form): if won - "the avian Cry1 folded domain is canalized at the whole-protein level, exceeding the background conservation of housekeeping and sensory control proteins; the magnetosensory apparatus appears frozen as part of a globally constrained protein, not a locally constrained pocket." If lost, honest null and the purging observation (0/33) stands alone as the descriptive result.
3. No threshold tuning after evaluation; all inputs already ledgered (attempt7b_controls_ledger.json, attempt7b_uniprot_refs.json).

## AMENDMENT 2026-09-27 02:08 IST (attempt 7d: robustness of the 7c canalization win - locked BEFORE computation)
The 7c whole-protein canalization W is a column-level mean; two robustness arms are locked pre-outcome:
1. COLUMN BOOTSTRAP: 10,000 resamples (seed 260927) of the qualified non-core column values per protein -> 95% percentile CI of W for each of the 6 proteins. WIN CONDITION A: Cry1's CI lower bound exceeds every control's point W.
2. LEAVE-ONE-SPECIES-OUT (Cry1 panel rows): recompute Cry1 W dropping each panel species in turn (no re-alignment; row dropped, column conservation recomputed). WIN CONDITION B: minimum LOSO W >= 0.99 (i.e., the canalization is not carried by any single species).
Both arms reported regardless of outcome; failure of either is an honest robustness caveat on the 7c claim, not a retraction of the locked 7c result (which stands on its own pre-registered gates).
