# Attempt 15: power analysis - detectable effect sizes at committed n
Amendment: PRE-REGISTRATION.md 2026-09-27 13:17 IST (locked pre-compute). Data: results/h_power_analysis.json.
Method: Hanley-McNeil (1971) AUC standard error; one-sided alpha=0.05, power=0.80, H1: AUC>0.5; committed counts only.
## Results
- Full panel (n=35 migratory vs 80 sedentary): minimum detectable AUC = 0.6453. Observed AUC = 0.3414 (SE 0.0527) - the observed effect is not merely non-significant, it points opposite to the tested direction.
- Subpanel arms (5v5): MDA = 0.9042 - only near-perfect separation is detectable; all subpanel arms are power-limited to near-ceiling effects.
- Reference points: 10v10 -> 0.8029; 50v50 -> 0.6419.
## Positive content
Cross-method consistency: the independent Hanley-McNeil MDA (0.6453) agrees with the locked simulation-based MDE (0.6424, attempt 8g) to within 0.003 - two independent power frameworks bracket the same detectability boundary, strengthening the attempt-8g conclusion that the RF benchmark (0.6234) and logistic baseline (0.4021) both sit below the detectable threshold. This converts the paper's negative results into quantified power statements: the panel can only detect AUC >= 0.645, and nothing observed approaches that. No redesign, no species additions (anti-goals).
