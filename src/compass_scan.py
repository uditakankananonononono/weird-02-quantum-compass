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

TRP_RING = ('CD1','CD2','CE2','CE3','CZ2','CZ3','CH2','NE1','CG')
import re
ISO = re.compile(r'^(N1|C2|O2|N3|C4|O4|C4A|C4X|N5|C5A|C5X|C6|C7|C7M|C8|C8M|C9|C9A|C9X|N10|C10)$')

def _load_6ptz():
    cas, ring, fad = {}, {}, []
    for line in open('/tmp/6ptz.pdb'):
        if line.startswith('HETATM') and line[17:20].strip() == 'FAD' and ISO.match(line[12:16].strip()):
            fad.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        if line.startswith('ATOM'):
            r = int(line[22:26]); atom = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if atom == 'CA': cas[r] = xyz
            if line[17:20].strip() == 'TRP' and atom in TRP_RING: ring.setdefault(r, []).append(xyz)
    return cas, ring, np.array(fad)

def d_term_edge(species_pdb, chain_residx_terminal, pairs):
    """pairs = [(species_residx, clcry4/6PTZ_residx)] from MSA columns; SVD superpose."""
    cas6, ring6, fad6 = _load_6ptz()
    casS = {}
    ringS = {}
    for line in open(species_pdb):
        if line.startswith('ATOM'):
            r = int(line[22:26]); atom = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if atom == 'CA': casS[r] = xyz
            if line[17:20].strip() == 'TRP' and atom in TRP_RING: ringS.setdefault(r, []).append(xyz)
    pairs = [(vs, u6) for vs, u6 in pairs if vs in casS and u6 in cas6]
    if len(pairs) < 100: raise ValueError(f'too few CA pairs for superposition: {len(pairs)}')
    A = np.array([casS[vs] for vs, _ in pairs]); B = np.array([cas6[u6] for _, u6 in pairs])
    # SVD superposition A(species) -> B(6PTZ): map species Trp ring into 6PTZ frame
    Ac, Bc = A.mean(0), B.mean(0)
    H = (A - Ac).T @ (B - Bc)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0: Vt[-1] *= -1; R = Vt.T @ U.T
    term = chain_residx_terminal
    if term not in ringS: raise ValueError(f'terminal Trp {term} missing ring atoms')
    ringS_6 = [(R @ (x - Ac)) + Bc for x in ringS[term]]
    return min(float(np.linalg.norm(x - f)) for x in ringS_6 for f in fad6)

def tau_eff(d_term):
    k_back = K0 * np.exp(-BETA * max(d_term - D0, 0.0))
    return 1.0 / max(k_back * 1e-6, 1e-12)   # microseconds

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
    # terminal Trp: last consecutive W walking FAD->surface chain columns
    chain_res = definition['ref_chain_residx']  # erCry4 numbering == 6PTZ numbering
    pairs = []
    if row is not None:
        col2res = {}
        ri = 0
        for c, ch in enumerate(row):
            if ch != '-': ri += 1; col2res[c] = ri
        species_res = [col2res.get(c) for c in definition['msa_cols']]
        # full correspondence set: species residx <-> ClCry4(=6PTZ) residx per column
        _, msa_all = msa_row(a.msa_row)
        clkey = [k for k in msa_all if 'A0A386QUR4' in k][0]
        clrow = msa_all[clkey]
        ri_cl = 0; cl_col2res = {}
        for c, ch in enumerate(clrow):
            if ch != '-': ri_cl += 1; cl_col2res[c] = ri_cl
        for c in range(min(len(row), len(clrow))):
            if c in col2res and c in cl_col2res:
                pairs.append((col2res[c], cl_col2res[c]))
    else:
        species_res = [p[0] for p in pts]
    terminal_idx = 0
    if aas is not None:
        for j, x in enumerate(aas):
            if x == 'W': terminal_idx = j
            else: break
    broken = aas is not None and aas[0] != 'W'
    if broken:
        d_term, tau = 3.6, tau_eff(3.6)   # proximal W absent -> contact recombination
    else:
        d_term = d_term_edge(a.pdb, species_res[terminal_idx], pairs)
        tau = tau_eff(d_term)
    A, phimax = anisotropy_50uT(tau)
    out = {'id': rec.id, 'chain_aas': aas, 'chain_broken': broken,
           'terminal_chain_pos': terminal_idx, 'd_term_A': round(d_term,2),
           'tau_eff_us': round(tau,4),
           'anisotropy_50uT': round(A,5), 'phi_s_max': round(phimax,5),
           'model': 'M1 locked 2026-09-25 20:29 IST; solver B1/B2/B5-validated'}
    print(json.dumps(out, indent=1) if a.json else
          f"{rec.id}: A@50uT={A:.4f} tau={tau:.3f}us d_term={d_term:.2f}A terminal_pos={terminal_idx}")

if __name__ == '__main__':
    main()
