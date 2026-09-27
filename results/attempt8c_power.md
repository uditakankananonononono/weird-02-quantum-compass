# Queue #15: minimum detectable effects for the migration-specific-signature nulls
Amendment 2026-09-27 11:17 IST locked BEFORE computation. Retrospective MDE at power 0.80 (not observed power); one-sided alpha 0.05 for 7f-b and the classifier (matching their locked tests), two-sided for 7g-b. Null SDs recovered by re-running the exact locked permutation streams (seed 260927, 10,000 permutations); every observed statistic reproduced its committed JSON value (sanity locks in the JSON).

## Results
| test | n (mig/sed) | observed | MDE (native units) |
|---|---|---|---|
| 7f-b gxp Cry1 DeltaW | 36/80 | +0.00142 | 0.0055 |
| 7f-b gxp Cry2 DeltaW | 35/77 | +0.00038 | 0.0042 |
| 7f-b intersection Cry1 | 35/77 | +0.00063 | 0.0061 |
| 7f-b intersection Cry2 | 35/77 | +0.00038 | 0.0042 |
| 7g-b tail net charge Cry1 | 36/80 | -0.00000 | 1.1e-05 |
| 7g-b tail net charge Cry2 | 35/73 | -0.00170 | 0.0022 |
| 7g-b charged density Cry1 | 36/80 | -0.00002 | 9.7e-05 |
| 7g-b charged density Cry2 | 35/73 | +0.00091 | 0.0026 |
| LOO classifier AUC (best: attempt 3a) | 36/80 | 0.402 | 0.642 |

## Locked framing
At the locked panel sizes the nulls exclude migration-specific effects at or above each MDE; smaller effects cannot be excluded. Concretely: a migratory-specific conservation elevation above ~0.5 percentage points of background conservation (DeltaW >= 0.0055) is excluded - the observed +0.0014 is fourfold below the detection floor; tail-charge separations above ~0.002 charge/residue are excluded (the marginal Cry2 net-charge signal, -0.0017 at p=0.024, sits BELOW its own MDE of 0.0022, i.e. even that signal is in the underpowered regime, consistent with the wrong-gene interpretation); and a sequence/structure classifier would need AUC >= 0.64 to be detectable at this panel size - the locked win thresholds were near or below the detectable range, which bounds what the classification nulls can say. No gate is retuned; this is interpretive context for the existing negatives.
