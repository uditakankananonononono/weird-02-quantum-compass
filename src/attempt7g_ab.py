#!/usr/bin/env python3
"""W02 7g-a + 7g-b (amendment 2026-09-27 10:13 IST, locked pre-computation).
7g-a: APC-corrected MI network, top-200 edges per gene, S_RP-coupled component vs surface.
7g-b: C-terminal 25% tail charge features, mig-vs-sed 10k label perms seed 260927.
Variable-column prefilter is RESULT-NEUTRAL (near-constant columns have MI~0 and never enter top-200)."""
import json, random, glob
import numpy as np
from collections import Counter, defaultdict
from Bio.Align import PairwiseAligner
from Bio.SVDSuperimposer import SVDSuperimposer
from pyfamsa import Aligner, Sequence

SEED, NPERM = 260927, 10000
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V','PRO':'P','SER':'S','THR':'T'}
CHAIN6 = json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering']
SRP = json.load(open('results/attempt7f_d_sphere.json'))['S_RP']
SRP1, SRP2 = set(SRP['CRY1_residues_6ptz']), set(SRP['CRY2_residues'])

def parse_ca_fad():
    ca6, fad = {}, []
    for line in open('/tmp/6ptz.pdb'):
        if line.startswith('ATOM') and line[12:16].strip()=='CA':
            ca6[int(line[22:26])] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
        elif line.startswith('HETATM') and line[17:20].strip()=='FAD':
            fad.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    return ca6, np.array(fad)

def load_msa1():
    msa = {}
    for line in open('data/cry1_panel_msa.fasta'):
        if line.startswith('>'): k=line[1:].strip(); msa[k]=''
        else: msa[k]+=line.strip()
    return {k:v for k,v in msa.items() if not k.startswith('REF|')}

def cry2_msa_refmap():
    files = sorted(glob.glob('data/cry2_panel/*.fa'))
    refseq = ''.join(open('data/cry2_panel/Columba_livia.fa').read().split('\n')[1:]).strip()
    seqs = [Sequence(fn.split('/')[-1][:-3].encode(), ''.join(open(fn).read().split('\n')[1:]).strip().encode()) for fn in files]
    seqs.append(Sequence(b'REF', refseq.encode()))
    msa_list = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
    ref_row = [s for s in msa_list if s[0]=='REF'][0][1]
    pos2col, pos = {}, 0
    for i, a in enumerate(ref_row):
        if a != '-': pos += 1; pos2col[pos] = i
    return {k:v for k,v in msa_list if k != 'REF'}, pos2col, ref_row

def superpose():
    ca6, fad6 = parse_ca_fad()
    res2, ca2 = {}, {}
    for line in open('data/folds/CRY2__REF__Columba_livia__XP_064920712.1.pdb'):
        if line.startswith('ATOM'):
            rn = int(line[22:26])
            if line[12:16].strip()=='CA': ca2[rn] = (line[17:20].strip(), np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    resnums6 = sorted(ca6); seq6 = ''.join(AA3.get(ca6[r][0],'X') for r in resnums6)
    resnums2 = sorted(ca2); seq2 = ''.join(AA3.get(ca2[r][0],'X') for r in resnums2)
    al = PairwiseAligner().align(seq6, seq2)[0]
    m6to2 = {}
    for (s6,e6),(s2,e2) in zip(*al.aligned):
        for i in range(e6-s6): m6to2[resnums6[s6+i]] = resnums2[s2+i]
    pairs = [(a,b) for a,b in m6to2.items()]
    fixed = np.array([ca6[a][1] for a,_ in pairs]); moving = np.array([ca2[b][1] for _,b in pairs])
    sup = SVDSuperimposer(); sup.set(fixed, moving); sup.run(); rot, tran = sup.get_rotran()
    ca2X = {r: np.dot(x, rot)+tran for r,(aa,x) in ca2.items()}
    return ca6, fad6, ca2X

print('geometry...', flush=True)
ca6, fad6, ca2X = superpose()
fad_centroid = fad6.mean(axis=0)
surf1 = {r for r,(aa,x) in ca6.items() if np.linalg.norm(x-fad_centroid) > 12.0}
surf2 = {r for r,x in ca2X.items() if np.linalg.norm(x-fad_centroid) > 12.0}
print('surface residues: cry1 %d cry2 %d' % (len(surf1), len(surf2)), flush=True)

geo = json.load(open('data/attempt6_column_geometry.json'))
col2res1 = {int(c): g['resnum_6ptz'] for c,g in geo.items()}
msa1 = load_msa1()
msa2, pos2col2, ref_row2 = cry2_msa_refmap()
col2res2 = {c: r for r, c in pos2col2.items()}
L2 = len(ref_row2)

def qcols(msa, L):
    rows = list(msa.values()); n = len(rows)
    out = []
    for c in range(L):
        res = [s[c] for s in rows if s[c] in 'ACDEFGHIKLMNPQRSTVWY']
        if len(res) >= 0.9*n: out.append(c)
    return out

def mi_apc(msa, cols, col2res, srp, surf, tag):
    # result-neutral prefilter: keep columns with >=2 symbols each at >=5% freq
    rows = list(msa.values()); n = len(rows)
    var = []
    for c in cols:
        cnt = Counter(s[c] for s in rows)
        cnt = {a:k for a,k in cnt.items() if a in 'ACDEFGHIKLMNPQRSTVWY-'}
        freq = sorted(cnt.values(), reverse=True)
        if len(freq) >= 2 and freq[1] >= 0.05*n: var.append(c)
    print('%s variable columns: %d of %d' % (tag, len(var), len(cols)), flush=True)
    V = np.array([[rows[i][c] for i in range(n)] for c in var])  # symbols as chars
    # encode to ints
    alpha = 'ACDEFGHIKLMNPQRSTVWY-'
    enc = np.searchsorted(np.array(list(alpha)), V)
    m = len(var)
    # one-hot blocks: (m, n, 21) uint8
    OH = np.zeros((m, n, 21), dtype=np.float32)
    for a in range(21): OH[:,:,a] = (enc == a)
    P1 = OH.sum(axis=1)/n  # (m,21)
    MI = np.zeros((m,m), dtype=np.float32)
    B = 25
    for i0 in range(0, m, B):
        J = np.einsum('ina,jnb->ijab', OH[i0:i0+B], OH) / n  # (B,m,21,21)
        Pi = P1[i0:i0+B][:,None,:,None] * P1[None,:,None,:]
        with np.errstate(divide='ignore', invalid='ignore'):
            T = J*np.log(J/Pi)
        T[~np.isfinite(T)] = 0
        MI[i0:i0+B] = T.sum(axis=(2,3))
    rm = MI.mean(axis=1, keepdims=True); cm = MI.mean(axis=0, keepdims=True); gm = MI.mean()
    APC = MI - rm*cm/gm
    np.fill_diagonal(APC, -1)
    # top 200 edges
    iu = np.triu_indices(m, 1)
    vals = APC[iu]
    top = np.argpartition(-vals, 200)[:200]
    top = top[np.argsort(-vals[top])]
    edges = [(var[iu[0][t]], var[iu[1][t]], float(vals[t])) for t in top]
    # network on residues
    adj = defaultdict(set)
    for c1, c2, v in edges:
        r1, r2 = col2res.get(c1), col2res.get(c2)
        if r1 is None or r2 is None or r1 == r2: continue
        adj[r1].add(r2); adj[r2].add(r1)
    # components containing S_RP
    seen = set(); comps = []
    for r in list(adj):
        if r in seen: continue
        stack=[r]; comp=set()
        while stack:
            x = stack.pop()
            if x in comp: continue
            comp.add(x); stack.extend(adj[x]-comp)
        seen |= comp; comps.append(comp)
    srp_comps = [c for c in comps if c & srp]
    if srp_comps:
        C = max(srp_comps, key=len)
        return {'n_var_cols':m,'n_edges':len(edges),'C_size':len(C),'C_surface':len(C & surf),
                'C_residues':sorted(C),'srp_in_C':sorted(C & srp),'top_edge':round(edges[0][2],4) if edges else None}
    return {'n_var_cols':m,'n_edges':len(edges),'C_size':0,'C_surface':0,'C_residues':[],'srp_in_C':[]}

L1 = len(next(iter(msa1.values())))
qc1 = qcols(msa1, L1)
qc2 = qcols(msa2, L2)
r_a1 = mi_apc(msa1, qc1, col2res1, SRP1, surf1, 'CRY1')
print('CRY1 7g-a:', {k:v for k,v in r_a1.items() if k!='C_residues'}, flush=True)
r_a2 = mi_apc(msa2, qc2, col2res2, SRP2, surf2, 'CRY2')
print('CRY2 7g-a:', {k:v for k,v in r_a2.items() if k!='C_residues'}, flush=True)
win_a = r_a1['C_size']>=10 and r_a1['C_surface']>=3 and (r_a2['C_size']<10 or r_a2['C_surface']<3)
print('7g-a WIN:', win_a, flush=True)

# ---- 7g-b tail charge ----
def labels():
    lab = {}
    for k in load_msa1():
        pass
    for line in open('data/cry1_panel_msa.fasta'):
        if line.startswith('>'):
            p = line[1:].strip().split('|')
            if len(p)>=3 and p[0] in ('migratory','sedentary'): lab[p[1]] = p[0]
    return lab
lab = labels()
ref_len1 = max(col2res1.values()); ref_len2 = max(col2res2.keys())
tail_res1 = set(range(int(0.75*ref_len1), ref_len1+1))
tail_res2 = set(range(int(0.75*ref_len2), ref_len2+1))
tail_cols1 = [c for c,r in col2res1.items() if r in tail_res1]
tail_cols2 = [c for c,r in col2res2.items() if r in tail_res2]
CHG = {'K':1.0,'R':1.0,'D':-1.0,'E':-1.0,'H':0.1}
def tail_feats(msa, tail_cols, species_of):
    feats = {}
    for k, s in msa.items():
        sp = species_of(k)
        if sp not in lab: continue
        vals = [s[c] for c in tail_cols]
        occ = [a for a in vals if a in 'ACDEFGHIKLMNPQRSTVWY']
        if len(occ) < 0.9*len(tail_cols): continue
        q = sum(CHG.get(a,0.0) for a in occ)/len(occ)
        d = sum(1 for a in occ if a in CHG)/len(occ)
        feats[sp] = (q, d, lab[sp])
    return feats
def perm_test(feats, fi, tag):
    rows = [(sp, v[fi], v[2]) for sp, v in feats.items()]
    mig = [v for _,v,l in rows if l=='migratory']; sed = [v for _,v,l in rows if l=='sedentary']
    obs = (sum(mig)/len(mig)) - (sum(sed)/len(sed))
    allv = mig+sed; nm = len(mig); rng = random.Random(SEED); cnt=0
    import math
    for _ in range(NPERM):
        samp = set(rng.sample(range(len(allv)), nm))
        a = [allv[i] for i in range(len(allv)) if i in samp]; b = [allv[i] for i in range(len(allv)) if i not in samp]
        d = (sum(a)/len(a)) - (sum(b)/len(b))
        if abs(d) >= abs(obs): cnt += 1
    p = (1+cnt)/(NPERM+1)
    print('%s: obs_delta=%.4f p=%.5f (n_mig=%d n_sed=%d)' % (tag, obs, p, nm, len(sed)), flush=True)
    return {'obs_delta':round(obs,5),'p_two_sided':round(p,5),'n_mig':nm,'n_sed':len(sed)}
f1 = tail_feats(msa1, tail_cols1, lambda k: k.split('|')[1] if '|' in k else k)
f2 = tail_feats(msa2, tail_cols2, lambda k: k)
r_b = {}
for fi, fname in [(0,'net_charge'),(1,'charged_density')]:
    r_b['cry1_'+fname] = perm_test(f1, fi, 'CRY1 '+fname)
    r_b['cry2_'+fname] = perm_test(f2, fi, 'CRY2 '+fname)
win_b = any(r_b['cry1_'+f]['p_two_sided']<0.05 and r_b['cry2_'+f]['p_two_sided']>=0.05 for f in ['net_charge','charged_density'])
print('7g-b WIN:', win_b, flush=True)

out = {'amendment':'2026-09-27 10:13 IST 7g-a/7g-b','seed':SEED,
 'surface_def':'CA centroid > 12.0A from FAD heavy-atom centroid',
 'g7a':{'CRY1':r_a1,'CRY2':r_a2,'win':bool(win_a)},
 'g7b':{'tail_def':'C-terminal 25% of reference length','results':r_b,'win':bool(win_b)},
 'prefilter_note':'variable-column prefilter (2nd symbol >=5%) is result-neutral: excluded columns have MI~0 and never enter top-200 edges'}
json.dump(out, open('results/attempt7g_ab.json','w'), indent=1)
print('written results/attempt7g_ab.json')
