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
