"""Benchmark tests for the W02 radical-pair simulator (pre-reg G3 anchors).
1) Isotropic hyperfine coupling -> zero anisotropy (symmetry argument).
2) Axial anisotropic hyperfine -> nonzero anisotropy, extremal at theta=0 vs pi/2.
3) Field-strength sanity: yield differs between 0 mT and geomagnetic 0.05 mT.
"""
import sys, math
sys.path.insert(0, "src")
import numpy as np
from radical_pair import RadicalPair, anisotropy

ANGLES = [0.0, math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2]

def yields_for(rp, B=0.05, t_max=5.0, dt=0.05):
    return [rp.singlet_yield(B, th, 0.0, t_max=t_max, dt=dt) for th in ANGLES]

def test_isotropic_zero_anisotropy():
    rp = RadicalPair(hyperfine_A=[np.eye(3) * 1.0], kS=1.0, kT=1.0)
    a = anisotropy(yields_for(rp))
    print(f"isotropic anisotropy: {a:.6f}")
    assert a < 1e-6, f"isotropic case must have ~zero anisotropy, got {a}"

def test_axial_nonzero_anisotropy():
    # axially symmetric tensor (N5-like): A_perp small, A_parallel large
    A = np.diag([-0.5, -0.5, 2.0])
    rp = RadicalPair(hyperfine_A=[A], kS=1.0, kT=1.0)
    ys = yields_for(rp)
    a = anisotropy(ys)
    print(f"axial anisotropy: {a:.6f}; yields: {[round(y,4) for y in ys]}")
    assert a > 1e-4, "axial tensor must produce nonzero anisotropy"
    assert ys[0] != ys[-1], "theta=0 and theta=pi/2 must differ"

def test_field_strength_sensitivity():
    A = np.diag([-0.5, -0.5, 2.0])
    rp = RadicalPair(hyperfine_A=[A], kS=1.0, kT=1.0)
    y0 = rp.singlet_yield(0.0, math.pi / 2, 0.0, t_max=5.0, dt=0.05)
    yG = rp.singlet_yield(0.05, math.pi / 2, 0.0, t_max=5.0, dt=0.05)
    print(f"yield 0mT={y0:.4f} vs 0.05mT={yG:.4f}")
    assert abs(yG - y0) > 1e-6, "geomagnetic field must shift the yield"

if __name__ == "__main__":
    test_isotropic_zero_anisotropy()
    test_axial_nonzero_anisotropy()
    test_field_strength_sensitivity()
    print("ALL BENCHMARK TESTS PASSED")
