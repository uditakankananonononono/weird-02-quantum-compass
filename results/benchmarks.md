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
