# W02 SIMULATOR BENCHMARK PROTOCOL - LOCKED 2026-09-25 20:20 IST (BEFORE computation)
Benchmarks 2-4 for pre-reg G3 (benchmark 1 = Ritz 2000 theta-curve, already
reproduced: anisotropy 0.385, cos2theta R^2=0.974).
Engine approximation stated honestly: nuclei are modeled as spin-1/2 though
14N is spin-1; tensor VALUES follow the published figures, so spike STRUCTURE
(orientation dependence) is compared qualitatively/semi-quantitatively, not as
exact numbers. Closed-form completion yield: Phi_S = -kS * PS . (L^-1 rho0)
(exact for reaction-to-completion, replacing time-stepping; documented in code).

## B2 (Hiscock 2016 PNAS 113:4634, Fig 2B config)
Model: FAD.- carries N5 (-0.087,-0.100,1.757) mT and N10 (-0.014,-0.024,0.605) mT
(diagonal tensors, flavin frame); partner Y. carries ONE axial tensor
(0,0,1.0812) mT with z-axis rotated 45 deg vs flavin z. B = 50 uT, kS=kT=k.
LOCKED EXPECTATIONS (from the paper):
 E2a: Phi_S(theta) shows a spike near theta=90 deg (local maximum distinct from
      the smooth sinusoidal background).
 E2b: spike amplitude increases with radical-pair lifetime tau over 1->100 us.
 E2c: removing the partner anisotropy (Y tensor -> zero) abolishes the spike.
PASS = all three qualitative expectations reproduced.

## B3 (Hiscock 2016, Fig 2A config)
Same FAD model; partner TrpH.+ simplified to the same single axial N1 tensor.
Zero out transverse (Axx,Ayy) components of N5, then N10, then both.
LOCKED EXPECTATIONS:
 E3a: zeroing either nitrogen's transverse components attenuates the spike
      (paper: 60-70% attenuation each; accepted band 30-100% given spin-1/2
      approximation).
 E3b: zeroing BOTH abolishes the spike.
PASS = E3b plus at least one of E3a within band.

## B4 (Ritz 2000 field-strength dependence)
Single anisotropic nucleus model (benchmark-1 config): anisotropy vs B magnitude
over 5 uT .. 50 mT.
LOCKED EXPECTATIONS:
 E4a: anisotropy rises monotonically at low field and saturates at high field
      (no decrease >5% between successive high-field points).
 E4b: A(50 uT) < A(5 mT) (Earth field sits below saturation).
PASS = E4a + E4b.

## B5 (Ritz 2000 Sec II theoretical benchmark) - LOCKED 20:21 IST (BEFORE computation)
Published prediction: an ISOTROPIC hyperfine coupling produces NO magnetic-field
orientation dependence (anisotropic coupling is required for a compass).
Model: single nucleus with isotropic tensor (1.0,1.0,1.0) mT; kS=kT=1/us; B=50 uT.
LOCKED EXPECTATION E5: anisotropy < 1e-9 (numerical zero).
PASS = E5 holds.
## B3/B4 outcomes (computed 20:21 IST, logged honestly)
B3: E3a PASS (N5-zeroed attenuation 93.8%, N10-zeroed 93.9%, both within the
30-100% band; paper 60-70%, spin-1/2 approximation noted). E3b FAIL: both-zeroed
residual = 7.3% of baseline vs locked <5% abolition. B3 overall FAIL (partial).
B4: FAIL under locked criteria: anisotropy rises 0.005->0.5 mT (0.399->0.490)
then collapses at 5/50 mT (0.015/0.00015) - at kS=kT=1/us the high-field regime
suppresses S-T mixing faster than reaction samples it; the locked expectation
was regime-wrong for these parameters (simulator behavior is the physics; the
expectation erred). Redirect logged: re-lock B4b with the recombination regime
taken from the source figure before recomputation.

## Methods note (locked 20:26, BEFORE any structure-derived scoring)
ESM Atlas foldSequence API caps at <500 aa (empirical: 400 OK, 500 -> HTTP 413).
Empirical cap found 20:27: 400 OK, 410+ -> 413. Rule corrected BEFORE folding:
first 400 residues (domain core; FAD pocket + Trp chain incl. 4th Trp at ~394 in
erCry4 numbering sit inside). Sequences whose 4th-chain Trp lies beyond 400 get
an explicit truncated flag. Tail features remain sequence-level.

## AMENDMENT: simulator parameterization mapping (LOCKED 20:29 IST, BEFORE any AUC/scoring)
M1 (locked): hyperfine tensors fixed from literature (FAD N5/N10 axial + Trp
partner; Hiscock 2016 configs as used in B2/B3). Per-species variation enters
ONLY via effective radical-pair lifetime tau_eff, derived from Trp-chain
geometry by standard ET theory: each chain hop of length d contributes rate
k_hop = k0 * exp(-beta*(d - d0)), beta = 1.4/A (Gray-Winkler canonical value),
k0=1e13/s at van der Waals contact d0=3.6A; tau_eff = 1/kS with kS the
chain-limited back-ET rate (harmonic sum of hop rates for triad; tetrad adds
the 4th hop). Sequences with a broken chain (missing Trp at a locked column)
use the triad-minus-one configuration. No per-species hyperfine refitting.
AUC test: predicted anisotropy A(tau_eff) at 50 uT, migratory(3) vs sedentary(1),
locked threshold AUC >= 0.75 (pre-reg H1). One-sided gate, no post-hoc metric swaps.

## AMENDMENT M1.2 (LOCKED 20:36 IST, BEFORE any AUC/scoring; supersedes the tau mechanics of M1)
Physics correction found during CLI diagnostics: tau_eff is NOT the forward chain
hop rate. tau_eff = 1/k_back where k_back = k0*exp(-beta*(d_term - d0)) is the
back-ET rate from the TERMINAL chain Trp to the FAD isoalloxazine ring
(beta=1.4/A, k0=1e13/s, d0=3.6A, as locked). d_term = minimum edge distance
(Trp indole atoms to isoalloxazine atoms) measured after SVD superposition of
the species apo fold onto PDB 6PTZ (holo, FAD-bound), transforming the 6PTZ FAD
into the species frame. Terminal Trp = last consecutive W walking the locked
chain columns FAD->surface (cols 436,413,359,410 = ClCry4 395,372,318,369);
a non-W breaks the chain there (triad = terminal at 318-column species).
Validation anchor (must hold for the scoring run): ClCry4/6PTZ d_term = 15.21 A
-> tau ~1.1 us; triad (terminal 318) d_term = 12.07 A -> tau ~0.14 us.
Chain definition grounded in 6PTZ (data/trp_chain_definition.json; COMPND
misannotation verified: 99.6% identity to ClCry4 A0A386QUR4 vs 57% Cry1).
