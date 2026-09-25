#!/usr/bin/env python3
"""Side-effect-free shared benchmark/simulator pieces (extracted 20:31)."""
import json, numpy as np
from scipy.linalg import solve
from radical_pair import RadicalPair, G_E

def completion_yield(rp, B_mT, theta, phi=0.0):
    """Phi_S = -kS * PS . (L^-1 rho0), exact reaction-to-completion."""
    H = rp.hamiltonian(B_mT, theta, phi)
    d = rp.dim
    K = 0.5 * (rp.kS * rp.PS + rp.kT * rp.PT)
    L = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d))) \
        - (np.kron(np.eye(d), K) + np.kron(K.T, np.eye(d)))
    rho0 = (rp.PS / np.trace(rp.PS)).reshape(-1, order='F')
    X = solve(L, rho0)
    PS_flat = rp.PS.reshape(-1, order='F')
    return float(np.real(-rp.kS * (PS_flat @ X)))

def rot_y_45(T):
    c = s = np.sqrt(0.5)
    R = np.array([[c,0,s],[0,1,0],[-s,0,c]])
    return R @ T @ R.T

N5 = np.diag([-0.087, -0.100, 1.757])
N10 = np.diag([-0.014, -0.024, 0.605])
Y45 = rot_y_45(np.diag([0.0, 0.0, 1.0812]))
THETAS = np.linspace(0, np.pi, 37)

def curve(rp, B=0.05):
    return np.array([completion_yield(rp, B, t) for t in THETAS])
