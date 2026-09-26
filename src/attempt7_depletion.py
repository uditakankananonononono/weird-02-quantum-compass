#!/usr/bin/env python3
"""W02 attempt-7: functional-core depletion / constraint mapping.
LOCKED 2026-09-26 20:56 IST amendment (PRE-REGISTRATION.md) BEFORE this evaluation.
Primary test: D = (# variant columns in core)/(core size).
Null: 10,000 stratified shuffles of the 33 observed variant columns over all
gap-qualified folded-domain columns, preserving conservation-quartile
distribution (bins from the observed set). Strict arm adds consensus-residue
volume class matching. Burial matching via per-residue SASA NOT possible
(sasa_reference_folds.json holds totals only) - logged as limitation.
Win: one-sided permutation p < 0.05 AND depletion ratio reported.
RNG seed locked 260926 (same as attempt-6).
"""
import json, math, random
from collections import Counter

D = 'data/'
GAPCAP = 0.90; NPERM = 10000; RNG_SEED = 260926
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,
       'H':153.2,'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,
       'T':116.1,'W':227.8,'Y':193.6,'V':140.0}

msa = {}
for line in open(D+'cry1_panel_msa.fasta'):
    if line.startswith('>'): k=line[1:].strip(); msa[k]=''
    else: msa[k]+=line.strip()
rows = [k for k in msa if not k.startswith('REF|')]
L = len(next(iter(msa.values())))

geo = json.load(open(D+'attempt6_column_geometry.json'))
chain6 = set(json.load(open(D+'trp_chain_definition.json'))['chain_6ptz_numbering'])
chain_cols = {int(c) for c,g in geo.items() if g['resnum_6ptz'] in chain6}
fad_d = {int(c): g['fad_d'] for c,g in geo.items()}
trp_d = {int(c): g['trp_d'] for c,g in geo.items()}
core = {c for c in fad_d if c in chain_cols or fad_d[c]<=4.5 or trp_d[c]<=6.0}

# qualified columns + per-column conservation + consensus residue
qual, consv, consaa = [], {}, {}
for col in range(L):
    res = [msa[r][col] for r in rows if msa[r][col] in VOL]
    if len(res) < GAPCAP*len(rows) or col not in fad_d: continue
    c, f = Counter(res).most_common(1)[0]
    qual.append(col); consv[col] = f/len(res); consaa[col] = c

sv = sorted(consv.values())
def quartile(v):
    import bisect; return bisect.bisect_right(sv, v)*4//(len(sv)+1)
def vclass(a): return 0 if VOL[a]<120 else (1 if VOL[a]<170 else 2)

var_cols = [c for c in qual if len({msa[r][c] for r in rows if msa[r][c] in VOL}) > 1]
D_obs = sum(c in core for c in var_cols)/len(core)

rng = random.Random(RNG_SEED)
byq = {}
for c in qual: byq.setdefault(quartile(consv[c]), []).append(c)
obs_q = Counter(quartile(consv[c]) for c in var_cols)
byqv = {}
for c in qual: byqv.setdefault((quartile(consv[c]), vclass(consaa[c])), []).append(c)
obs_qv = Counter((quartile(consv[c]), vclass(consaa[c])) for c in var_cols)

def draw_null(strata, obs_counts):
    cols = []
    for s, n in obs_counts.items():
        pool = strata.get(s, [])
        if len(pool) < n: return None
        cols += rng.sample(pool, n)
    return sum(c in core for c in cols)/len(core)

nulls, nulls_strict = [], []
tries = 0
while len(nulls) < NPERM and tries < 20*NPERM:
    tries += 1
    v = draw_null(byq, obs_q)
    if v is not None: nulls.append(v)
    vs = draw_null(byqv, obs_qv)
    if vs is not None: nulls_strict.append(vs)

def summarize(nulls):
    nulls.sort()
    p = (1 + sum(x <= D_obs for x in nulls))/(len(nulls)+1)
    med = nulls[len(nulls)//2]
    return {'n': len(nulls), 'p_one_sided_depletion': round(p,5),
            'median_null': round(med,4), 'depletion_ratio': round(D_obs/med,3) if med>0 else None}

out = {'locked_amendment': '2026-09-26 20:56 IST attempt-7', 'nperm': NPERM, 'seed': RNG_SEED,
 'core_size': len(core), 'n_variant_cols': len(var_cols), 'variant_cols_in_core': sum(c in core for c in var_cols),
 'D_obs': round(D_obs,4), 'null_conservation_quartile': summarize(nulls),
 'null_strict_quartile_x_volclass': summarize(nulls_strict) if len(nulls_strict)>=1000 else 'infeasible strata',
 'win_rule': 'p<0.05 one-sided + depletion ratio', 'p_primary': None}
out['p_primary'] = out['null_conservation_quartile']['p_one_sided_depletion']
out['win'] = out['p_primary'] < 0.05
json.dump(out, open('results/attempt7_depletion.json','w'), indent=1)
print(json.dumps(out, indent=1)[:900])
