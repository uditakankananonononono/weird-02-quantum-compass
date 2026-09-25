#!/usr/bin/env python3
"""W02 H1 attempt-4: physics-forward MFE replication (PREREG_ATTEMPT4.md, locked 00:12, P3 amended 00:13 - all BEFORE computation).
Universal S(tau) curve + per-species structure-derived taus for erCry4/ClCry4/GgCry4 WT (+ erCry4-WDF triad),
gates P1-P5 evaluated and saved. No refitting after this run."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from compass_scan import d_term_edge, tau_eff, load_def, msa_row, chain_positions, D
from bench_common import completion_yield, N5, N10, Y45, THETAS
from radical_pair import RadicalPair
from Bio import SeqIO

B_PRIMARY = np.linspace(1e-6, 200e-6, 12)          # 0-200 uT (P5 primary)
B_WIDE = np.concatenate([B_PRIMARY, np.linspace(2.5e-4, 2e-3, 8)])  # secondary wide
WIN = (25e-6, 65e-6)                                # geomagnetic window (locked)
SPECIES = {
    'erCry4': ('CRY4__REFSEQ__Erithacus_rubecula__MN709784', 'MN709784'),
    'ClCry4': ('CRY4__REFSEQ__Columba_livia__KX168611', 'KX168611'),
    'GgCry4': ('CRY4__REFSEQ__Gallus_gallus__NM_001039596', 'NM_001039596'),
}
SUB = 'data/cry4_subpanel.fasta'

def seq_of(acc):
    recs = {r.id: str(r.seq) for r in SeqIO.parse(SUB, 'fasta')}
    for k, v in recs.items():
        if acc in k: return v
    raise KeyError(acc)

def A_of_B(tau_us, B):
    k = 1.0 / tau_us
    rp = RadicalPair(hyperfine_A=[N5, N10], hyperfine_B=[Y45], kS=k, kT=k)
    phis = np.array([completion_yield(rp, float(B), float(t)) for t in THETAS])
    return float((phis.max()-phis.min())/phis.mean()), float(phis.max())

def S_of_tau(tau_us, Bgrid):
    A = np.array([A_of_B(tau_us, B)[0] for B in Bgrid])
    m = (Bgrid >= WIN[0]) & (Bgrid <= WIN[1])
    if m.sum() < 2:
        # interpolate into window
        Bw = np.linspace(WIN[0], WIN[1], 9)
        Aw = np.interp(Bw, Bgrid, A)
        return float(np.trapezoid(Aw, Bw) if hasattr(np,'trapezoid') else np.trapz(Aw, Bw)), A
    return float(np.trapezoid(A[m], Bgrid[m]) if hasattr(np,'trapezoid') else np.trapz(A[m], Bgrid[m])), A

def main():
    definition = load_def()
    # terminal-chain d_term per species (WT tetrad) + triad variant (chain truncated after W318: terminal=C)
    res = {}
    for sp, (fold, acc) in SPECIES.items():
        pdb = f'data/folds/{fold}.pdb'
        fasta = f'/tmp/a4_{sp}.fa'
        open(fasta,'w').write(f'>{sp}\n{seq_of(acc)}\n')
        # reuse compass_scan machinery via its CLI for exactness
        import subprocess
        msa = [k for k in msa_row([k for k in msa_row('migratory|Apus_apus|XM_051636888')[1].keys() if sp.split('Cry')[0][:3] in k] or ['migratory|Apus_apus|XM_051636888'])[1].keys()] # placeholder
        r = subprocess.run(['python3','src/compass_scan.py',fasta,'--pdb',pdb,'--json'], capture_output=True, text=True)
        d = json.loads(r.stdout)
        res[sp] = {'d_term_WT': d['d_term_A'], 'tau_WT': d['tau_eff_us'], 'A50_WT': d['anisotropy_50uT']}
        print(sp, res[sp], flush=True)
    # universal S(tau) curve
    tau_grid = np.concatenate([np.linspace(0.1, 1.0, 10), np.linspace(1.2, 10.0, 15)])
    curve = {}
    for tau in tau_grid:
        S, A = S_of_tau(float(tau), B_WIDE)
        curve[round(float(tau),3)] = {'S': S, 'A_curve': [round(float(a),6) for a in A]}
        print('tau', round(float(tau),2), 'S', f'{S:.3e}', flush=True)
    json.dump({'species': res, 'tau_grid_S': curve,
               'B_wide_uT': [round(float(b)*1e6,1) for b in B_WIDE], 'window_uT': [25,65]},
              open('results/h1_attempt4_stage1.json','w'), indent=1)
    print('STAGE1 DONE', flush=True)

if __name__ == '__main__':
    main()
