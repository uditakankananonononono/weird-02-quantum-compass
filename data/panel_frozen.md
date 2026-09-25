# W02 SPECIES PANEL - FROZEN 2026-09-25 19:26 IST (BEFORE any H1 scoring)
## Data grounding (accession-level, logged toward G1)
- 997 avian cryptochrome nuccore accessions -> 943 main-ORF protein translations
  (54 records lack /translation: pseudogene/unannotated CDS; excluded, logged).
- Species resolved via NCBI esummary: 187 species (data/cry_aves_species.json).
- Migration labels: AVONET Supplementary dataset 1, BirdLife taxonomy
  (Tobias et al. 2022, Ecology Letters, doi:10.1111/ele.13898; figshare 16586228).
  173/187 mapped (14 excluded: NCBI-vs-BirdLife name mismatch, listed in
  data/migration_labels_raw.json provenance).
- Paralog labels: reference-seeded phylogenetic classification (50 title-named
  Cry1/Cry2/Cry4 + 349 DASH/photolyase outgroup seeds; FAMSA + FastTree WAG,
  3-nearest-reference vote): Cry1 296, Cry2 279, Cry4 12, DASH 313
  (data/paralog_labels.json).

## FROZEN RULES
R1. Primary H1 contrast: species with a classified CRY1 sequence AND AVONET
    Migration=3 (obligate migrant, n=36) vs Migration=1 (sedentary, n=80).
    Cry1 chosen as primary scaffold because it is the only paralog with
    panel-scale coverage; the avian-magnetoreception literature supports a
    retinal Cry1a role alongside Cry4.
R2. Migration=2 (partial migrants, n=37 class-Cry1 species where available)
    are HELD OUT of the primary contrast; reported separately.
R3. CRY4 COVERAGE GAP (honest limitation, logged pre-scoring): only 4 species
    carry a classified Cry4 accession - avian Cry4 is severely under-annotated
    in RefSeq (consistent with its recent characterization). Cry4-specific
    analyses are restricted to those 4 species and labeled low-power.
R4. REDIRECT for the Cry4 gap (pre-scoring, in-project): de novo Cry4 mining -
    Cry4-seed profile-HMM (hmmbuild on the 12 classified Cry4 sequences) scanned
    against additional avian proteomes (NCBI RefSeq/GenBank + Ensembl), with
    tree-verified orthology; any species gained join the Cry4 sub-panel under
    the same migration-label rules. Executed and reported before H2 mutagenesis.
R5. Per species, ONE representative sequence per paralog: longest classified
    accession; tie-break lowest accession number. Locked.
R6. H1 AUC >= 0.75 gate and redirect rules from PRE-REGISTRATION.md unchanged.
