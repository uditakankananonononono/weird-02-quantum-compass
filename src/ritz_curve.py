"""W02 G3 quantitative step: theta-dependence curve for the FAD-Trp compass model.
Model: FAD radical with axial N5 hyperfine (A_perp=-0.5145, A_par=+1.0243 mT,
Ritz 2000 Table 1 N5 values) + Trp radical with no hyperfine (minimal model).
Expectation (Ritz 2000): singlet yield varies with field direction; curve well
fit by a + b*cos(2*theta). Saves results/ritz2000_curve.json.
"""
import sys, json, math
sys.path.insert(0, "src")
import numpy as np
from radical_pair import RadicalPair, anisotropy

A_N5 = np.diag([-0.5145, -0.5145, 1.0243])  # mT, axial N5 of FAD (Ritz 2000)
rp = RadicalPair(hyperfine_A=[A_N5], hyperfine_B=[], kS=1.0, kT=1.0)

thetas = np.linspace(0, math.pi / 2, 10)
yields = [rp.singlet_yield(0.05, float(th), 0.0, t_max=5.0, dt=0.05) for th in thetas]

# fit y = a + b*cos(2 theta)
X = np.column_stack([np.ones_like(thetas), np.cos(2 * thetas)])
coef, res, *_ = np.linalg.lstsq(X, np.array(yields), rcond=None)
pred = X @ coef
ss_res = float(np.sum((np.array(yields) - pred) ** 2))
ss_tot = float(np.sum((np.array(yields) - np.mean(yields)) ** 2))
r2 = 1 - ss_res / ss_tot

out = {
    "model": "FAD(N5 axial) + Trp(no hf), kS=kT=1/us, B=0.05 mT",
    "thetas": [round(float(t), 4) for t in thetas],
    "singlet_yields": [round(y, 6) for y in yields],
    "anisotropy": round(anisotropy(yields), 6),
    "cos2theta_fit": {"a": round(float(coef[0]), 6), "b": round(float(coef[1]), 6), "R2": round(r2, 6)},
    "benchmark": "Ritz, Adem, Schulten 2000 Biophys J 78:707 - compass anisotropy from axial hyperfine",
}
open("results/ritz2000_curve.json", "w").write(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
