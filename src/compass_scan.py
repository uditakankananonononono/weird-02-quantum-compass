#!/usr/bin/env python3
"""compass-scan v0.2 - score a cryptochrome sequence for predicted magnetic sensitivity.
Locked M1 (20:29 IST): chain defined by erCry4-reference MSA columns
(data/trp_chain_definition.json); hops = consecutive chain-residue C-alpha
distances in the species' own fold; k_hop = 1e13*exp(-1.4*(d-3.6)); k_eff =
harmonic sum; tau_eff = 1/k_eff; anisotropy at 50 uT from the B1/B2/B5-validated
solver. Broken chain (non-W at a chain column) -> triad-minus-one default.
Usage: python3 compass_scan.py sequence.fasta --pdb fold.pdb [--msa-row NAME] [--json]
"""
import sys, json, argparse, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from bench_common import completion_yield, N5, N10, Y45, THETAS

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
BETA, D0, K0, B50 = 1.4, 3.6, 1e13, 50e-6

def load_def():
    return json.load(open(os.path.join(D, 'trp_chain_definition.json')))

def msa_row(name):
    msa = {}
    for line in open(os.path.join(D, 'cry1_panel_msa.fasta')):
        if line.startswith('>'): k = line[1:].strip(); msa[k] = ''
        else: msa[k] += line.strip()
    if name in msa: return msa[name], msa
    return None, msa

def chain_positions(seq, pdb_path, definition, aligned_row=None):
    """Return (aa_at_chain_columns, [(residx, ca_xyz)]) for the locked chain columns."""
    cols = definition['msa_cols']
    refkey_aln = None
    if aligned_row is None:
        # align novel sequence to erCry4 reference sequence, walk columns
        refseq = ''.join(c for c in open(os.path.join(D,'ref_sequences.fasta')).read()
                         .split('>tr|A0A2I4SZI9|')[1].split('>')[0].splitlines()[1:])
        al = PairwiseAligner().align(refseq, seq)[0]
        # build ref-residue -> query-residue map
        r2q = {}
        for (rs, re_), (qs, qe_) in zip(*al.aligned):
            for i in range(re_-rs): r2q[rs+i+1] = qs+i+1
        residx = [r2q.get(r) for r in definition['ref_chain_residx']]
    else:
        # column -> residue index in this row
        col2res = {}
        ri = 0
        for c, ch in enumerate(aligned_row):
            if ch != '-': ri += 1; col2res[c] = ri
        residx = [col2res.get(c) for c in cols]
        aas = [aligned_row[c] if c < len(aligned_row) else '-' for c in cols]
    # C-alpha coords from fold by residue index
    cas = {}
    for line in open(pdb_path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            cas[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    pts = [(i, cas[i]) for i in residx if i in cas]
    if aligned_row is None:
        aas = None
    return aas, pts

def k_hop(d): return K0 * np.exp(-BETA * max(d - D0, 0.0))

def tau_eff(hops, default=(15.4, 18.1)):
    gs = [g for g in hops if g <= 30.0]
    if len(gs) < 2: gs = list(default)
    k_eff = 1.0 / sum(1.0/k_hop(g) for g in gs)
    return 1.0 / max(k_eff * 1e-6, 1e-12)

def anisotropy_50uT(tau_us):
    k = 1.0 / tau_us
    from radical_pair import RadicalPair
    rp = RadicalPair(hyperfine_A=[N5, N10], hyperfine_B=[Y45], kS=k, kT=k)
    phis = np.array([completion_yield(rp, B50, float(t)) for t in THETAS])
    return float((phis.max()-phis.min())/phis.mean()), float(phis.max())

def main():
    ap = argparse.ArgumentParser(prog='compass-scan')
    ap.add_argument('fasta'); ap.add_argument('--pdb', required=True)
    ap.add_argument('--msa-row'); ap.add_argument('--json', action='store_true')
    a = ap.parse_args()
    rec = next(SeqIO.parse(a.fasta, 'fasta'))
    definition = load_def()
    row = None
    if a.msa_row:
        row, _ = msa_row(a.msa_row)
    aas, pts = chain_positions(str(rec.seq), a.pdb, definition, aligned_row=row)
    broken = aas is not None and any(x != 'W' for x in aas)
    hops = [float(np.linalg.norm(p1[1]-p0[1])) for p0, p1 in zip(pts, pts[1:])]
    tau = tau_eff([] if broken else hops)
    A, phimax = anisotropy_50uT(tau)
    out = {'id': rec.id, 'chain_aas': aas, 'chain_broken': broken,
           'hops_A': [round(h,2) for h in hops], 'tau_eff_us': round(tau,4),
           'anisotropy_50uT': round(A,5), 'phi_s_max': round(phimax,5),
           'model': 'M1 locked 2026-09-25 20:29 IST; solver B1/B2/B5-validated'}
    print(json.dumps(out, indent=1) if a.json else
          f"{rec.id}: A@50uT={A:.4f} tau={tau:.3f}us broken={broken} hops={out['hops_A']}")

if __name__ == '__main__':
    main()
