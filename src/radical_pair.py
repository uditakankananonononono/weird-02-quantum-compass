"""W02 radical-pair spin dynamics simulator (Haberkorn master equation).

Model: two electron spins-1/2 (radicals A, B) coupled to nuclear spins via
hyperfine tensors. Recombination from singlet (kS) and triplet (kT) channels.
Singlet yield Phi_S(B-direction) is the compass observable.

Reference model: Ritz, Adem, Schulten (2000) Biophys J 78:707-718 - one
anisotropically coupled nucleus (e.g. N5 of FAD) gives orientation-dependent
yield; isotropic hyperfine gives zero anisotropy.
Units: angular frequencies (rad/us) unless noted; 1 mT ~ 1.76e8 rad/s for electron.
"""
import numpy as np
from scipy.linalg import expm
from scipy.sparse import kron as sp_kron, eye as sp_eye, csc_matrix

G_E = 1.76086e5  # electron gyromagnetic ratio, rad / (us . mT)

def _spin_ops():
    sx = np.array([[0, 0.5], [0.5, 0]], complex)
    sy = np.array([[0, -0.5j], [0.5j, 0]], complex)
    sz = np.array([[0.5, 0], [0, -0.5]], complex)
    return sx, sy, sz

class RadicalPair:
    """Hilbert space: electron A x electron B x nuclei (each spin-1/2)."""
    def __init__(self, hyperfine_A=None, hyperfine_B=None, kS=1.0, kT=1.0):
        """hyperfine_*: list of 3x3 tensors (mT) for nuclei on each radical."""
        self.hfA = hyperfine_A or []
        self.hfB = hyperfine_B or []
        self.kS, self.kT = kS, kT
        self.n_nuc = len(self.hfA) + len(self.hfB)
        self.dim = 4 * (2 ** self.n_nuc)
        sx, sy, sz = _spin_ops()
        self._build_spin_ops(sx, sy, sz)
        self._build_projectors()

    def _kron_all(self, factors):
        out = factors[0]
        for f in factors[1:]:
            out = np.kron(out, f)
        return out

    def _build_spin_ops(self, sx, sy, sz):
        n = self.n_nuc
        eye2 = np.eye(2)
        def place(op, idx):
            return self._kron_all([op if i == idx else eye2 for i in range(2 + n)])
        self.SA = [place(sx, 0), place(sy, 0), place(sz, 0)]
        self.SB = [place(sx, 1), place(sy, 1), place(sz, 1)]
        self.nuclei = []
        for j in range(n):
            self.nuclei.append([place(sx, 2 + j), place(sy, 2 + j), place(sz, 2 + j)])

    def _build_projectors(self):
        # Singlet/triplet projectors on the two-electron subspace, lifted to full space
        d = self.dim
        SA2 = sum(S @ S for S in self.SA)
        SB2 = sum(S @ S for S in self.SB)
        SAdotB = sum(self.SA[i] @ self.SB[i] for i in range(3))
        S2 = SA2 + SB2 + 2 * SAdotB
        # S(S+1) eigenvalues: 0 (singlet), 2 (triplet)
        self.PS = (S2 - 2 * np.eye(d)) / (0 - 2) * -1.0  # = (2I - S2)/2... normalized below
        self.PS = (2 * np.eye(d) - S2) / 2.0
        # verify idempotent-ish (PS^2 ~ PS on spin subspace)
        self.PT = np.eye(d) - self.PS

    def hamiltonian(self, B_mT, theta, phi):
        """H = g * B . (SA + SB) + sum_i SA.A_i.I_i (+ SB.A_j.I_j), rad/us."""
        Bvec = B_mT * np.array([np.sin(theta) * np.cos(phi),
                                np.sin(theta) * np.sin(phi),
                                np.cos(theta)])
        H = G_E * sum(Bvec[i] * (self.SA[i] + self.SB[i]) for i in range(3))
        for A, I in zip(self.hfA, self.nuclei[:len(self.hfA)]):
            for a in range(3):
                for b in range(3):
                    H = H + G_E * A[a, b] * (self.SA[a] @ I[b])
        for A, I in zip(self.hfB, self.nuclei[len(self.hfA):]):
            for a in range(3):
                for b in range(3):
                    H = H + G_E * A[a, b] * (self.SB[a] @ I[b])
        return H

    def singlet_yield(self, B_mT, theta, phi, t_max=10.0, dt=0.02):
        """Phi_S via trace formula: kS * integral_0^tmax Tr[PS rho(t)] dt,
        rho(0) = PS / Tr(PS) (singlet-born). Haberkorn propagator."""
        H = self.hamiltonian(B_mT, theta, phi)
        d = self.dim
        K = 0.5 * (self.kS * self.PS + self.kT * self.PT)
        L = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d))) \
            - (np.kron(np.eye(d), K) + np.kron(K.T, np.eye(d)))
        rho0 = (self.PS / np.trace(self.PS)).reshape(-1, order="F")
        prop = expm(L * dt)
        steps = int(t_max / dt)
        acc = 0.0
        rho = rho0
        PS_flat = self.PS.reshape(-1, order="F")
        for _ in range(steps):
            rho = prop @ rho
            acc += float(np.real(PS_flat @ rho))
        return self.kS * acc * dt

def anisotropy(yields):
    yields = np.asarray(yields, float)
    return (yields.max() - yields.min()) / yields.mean()
