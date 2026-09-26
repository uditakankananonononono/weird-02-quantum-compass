#!/usr/bin/env python3
"""W02 7d robustness (amendment 02:08 IST, locked pre-computation).
Arm A: column bootstrap 95% CI of W per protein (10k, seed 260927).
Arm B: leave-one-species-out on Cry1 rows, min W. Reuses 7c column-value machinery."""
import json, random, glob
from collections import Counter
from pyfamsa import Aligner, Sequence

SEED, NBOOT = 260927, 10000
refs = json.load(open('data/attempt7b_uniprot_refs.json'))['refs']

def bg_cons(gene):
    files = sorted(glob.glob(f'data/controls7b/*__{gene}.fa'))
    seqs = [Sequence(fn.split('/')[-1].encode(), ''.join(open(fn).read().split('\n')[1:]).encode()) for fn in files]
    r = refs[gene]
    seqs.append(Sequence(b'REF', r['seq'].encode()))
    msa = [(s.id.decode(), s.sequence.decode()) for s in Aligner().align(seqs)]
    ref_row = [s for s in msa if s[0]=='REF'][0][1]
    core_cols, pos = set(), 0
    for i, a in enumerate(ref_row):
        if a != '-':
            pos += 1
            for f in r['features']:
                if f['start'] <= pos <= f['end']: core_cols.add(i)
    seqrows = [t[1] for t in msa if t[0] != 'REF']
    vals = []
    for i in range(len(ref_row)):
        if i in core_cols: continue
        res = [s[i] for s in seqrows if s[i] != '-']
        if len(res) >= 0.9*len(seqrows):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return vals

def cry1_rows_cols():
    msa = {}
    for line in open('data/cry1_panel_msa.fasta'):
        if line.startswith('>'): k=line[1:].strip(); msa[k]=''
        else: msa[k]+=line.strip()
    rows = [k for k in msa if not k.startswith('REF|')]
    geo = json.load(open('data/attempt6_column_geometry.json'))
    chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
    core_c = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6 or g['fad_d']<=4.5 or g['trp_d']<=6.0}
    geocols = {int(c) for c in geo}
    L = len(next(iter(msa.values())))
    qcols = [i for i in range(L) if i not in core_c and i in geocols]
    return msa, rows, qcols

def cryW(msa, rows, qcols, drop=None):
    use = [r for r in rows if r != drop]
    vals = []
    for i in qcols:
        res = [msa[r][i] for r in use if msa[r][i] in 'ACDEFGHIKLMNPQRSTVWY']
        if len(res) >= 0.9*len(use):
            vals.append(Counter(res).most_common(1)[0][1]/len(res))
    return sum(vals)/len(vals), vals

msa, rows, qcols = cry1_rows_cols()
W_cry, cry_vals = cryW(msa, rows, qcols)
pools = {'CRY1': cry_vals}
W = {'CRY1': W_cry}
for g in ['PKM','LDHA','RHO','OPN4','GAPDH']:
    v = bg_cons(g); pools[g] = v; W[g] = sum(v)/len(v)
    print(g, 'W=%.4f n=%d' % (W[g], len(v)), flush=True)
print('CRY1 W=%.4f n=%d' % (W['CRY1'], len(cry_vals)), flush=True)

rng = random.Random(SEED)
ci = {}
for g, vals in pools.items():
    n = len(vals); boots = []
    for _ in range(NBOOT):
        s = sum(rng.choice(vals) for _ in range(n))
        boots.append(s/n)
    boots.sort()
    ci[g] = [boots[int(0.025*NBOOT)], boots[int(0.975*NBOOT)]]
    print(g, 'CI95 [%.4f, %.4f]' % tuple(ci[g]), flush=True)
winA = all(ci['CRY1'][0] > W[g] for g in ci if g != 'CRY1')
loso = {}
for r in rows:
    w, _ = cryW(msa, rows, qcols, drop=r)
    loso[r] = round(w, 5)
min_loso = min(loso.values())
winB = min_loso >= 0.99
print('LOSO min W=%.5f (n=%d species)' % (min_loso, len(rows)), flush=True)
out = {'amendment':'2026-09-27 02:08 IST 7d','nboot':NBOOT,'seed':SEED,
       'W':{k: round(v,5) for k,v in W.items()},
       'W_bootstrap_CI95':{k:[round(x,5) for x in v] for k,v in ci.items()},
       'winA_cry1_CI_above_all_controls': winA,
       'loso_n_species': len(rows), 'loso_min_W': min_loso, 'loso_per_species': loso,
       'winB_loso_min_above_0.99': winB}
json.dump(out, open('results/attempt7d_robustness.json','w'), indent=1)
print('WIN_A', winA, 'WIN_B', winB)
