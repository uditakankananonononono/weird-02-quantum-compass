# W02 PREREG ATTEMPT-5B - locked 09:09 IST 2026-09-26, BEFORE any 5b metric evaluation
Source: ChatGPT consult round 3 rank-2; parent approved as separate lock 09:07.

## Hypothesis (tetrad-vs-triad tradeoff)
Adding the fourth tryptophan shifts the radical-pair response from maximum instantaneous sensitivity toward a prolonged magnetically responsive state: a persistence/dynamic-range advantage despite reduced peak sensitivity.

## Locked gates (from erCry4 locked M1.2 lifetimes: tetrad tau=0.1495us, triad tau=0.0034us; anisotropy-vs-field A(B) curves from the attempt-4 amended-resolution engine)
- B1: persistence integral ratio I_tetrad/I_triad >= 1.25, where I = A_peak x tau x (1-exp(-T/tau)), T=10us locked. HONEST NOTE (pre-registered): B1 is expected-by-construction (longer tau inflates the time integral) - it is recorded for completeness and carries no evidential weight alone.
- B2 (discriminating): tetrad peak anisotropy <= triad peak anisotropy (the tradeoff direction).
- B3 (independent anti-circularity metric, per consult): field discrimination width (FWHM of A(B) in uT) is preserved within 2x: |log2(FWHM_tetrad/FWHM_triad)| <= 1 - the persistence gain must not destroy field selectivity.
- Triad A(B) curve: triad tau 0.0034us lies below the amended grid minimum; locked rule = nearest computed grid point (0.1us) as triad proxy, noted as a limitation. Success = B1 AND B2 AND B3.

## AMENDMENT 1 (09:09, BEFORE extended computation): both A(B) curves are still rising at the amended grid edge (20mT), so FWHM is undefined on-grid (grid artifact, not a physics result). Locked extension: recompute A(B) for tau 0.1us (triad proxy) and 0.15us (tetrad) on an extended log-spaced grid 1e-6..0.5 T, 24 points, 9 orientations (same amended resolution); FWHM computed on the extended curve. If a curve still never falls to half-max, FWHM is reported as >grid-range and B3 is scored INCONCLUSIVE (attempt-5b verdict then rests on B2 alone, reported honestly). All other gates unchanged.
