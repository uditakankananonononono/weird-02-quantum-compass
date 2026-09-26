#!/usr/bin/env python3
"""W02 post-hoc robustness (judge round-5 weakness 5): leave-one-FAMILY-out
re-evaluation of the locked attempt-2 classifier. Same features, same estimator
(L2 logistic C=1.0, numpy Newton), same labels. Family = AVONET Family3.
Pipeline check: plain LOO must reproduce locked AUC 0.2066.
Labeled POST-HOC ROBUSTNESS - locked result already seen; no gate changes."""
import json, numpy as np, openpyxl

feats = json.load(open('data/cry1_features_v1.json'))
namemap = json.load(open('data/avonet_name_map.json'))
wb = openpyxl.load_workbook('data/avonet_supp1.xlsx', read_only=True)
ws = wb['AVONET3_BirdTree']
rows = ws.iter_rows(values_only=True)
hdr = list(next(rows))
i_sp = hdr.index('Species3'); i_fam = hdr.index('Family3')
fam = {}; gen2fam = {}
for r in rows:
    if r[i_sp]:
        fam[str(r[i_sp])] = str(r[i_fam])
        gen2fam[str(r[i_sp]).split()[0]] = str(r[i_fam])


# Taxonomic synonyms for AVONET3 (BirdTree) matching - recorded for audit
SYNONYMS = {
 'Ammospiza caudacuta': 'FAMILY:Passerellidae',      # genus absent in BirdTree; Ammospiza is Passerellidae
 'Ammospiza nelsoni': 'Ammodramus nelsoni',
 'Antrostomus carolinensis': 'Caprimulgus carolinensis',
 'Calidris pugnax': 'GENUS:Calidris',                # Philomachus pugnax absent; genus Calidris fixes family Scolopacidae
 'Calonectris borealis': 'Calonectris diomedea',
 'Motacilla alba alba': 'Motacilla alba',
 'Zonotrichia leucophrys gambelii': 'Zonotrichia leucophrys',
 'Apteryx rowi': 'Apteryx australis',
 'Balearica regulorum gibbericeps': 'Balearica regulorum',
 'Buceros rhinoceros silvestris': 'Buceros rhinoceros',
 'Chlamydotis macqueenii': 'Chlamydotis undulata',
 'Cyanistes caeruleus': 'Parus caeruleus',
 'Dryobates pubescens': 'Picoides pubescens',
 'Melozone crissalis': 'GENUS:Melozone',
 'Pithys albifrons albifrons': 'Pithys albifrons',
 'Pseudopipra pipra': 'Pipra pipra',
}

keys = sorted(feats)
X, y, F = [], [], []
missing = []
for k in keys:
    sp = k.split('|')[1].replace('_', ' ')
    av = namemap.get(k.split('|')[1], k.split('|')[1])
    syn = SYNONYMS.get(sp)
    if syn and syn.startswith('FAMILY:'):
        f = syn.split(':',1)[1]
    elif syn and syn.startswith('GENUS:'):
        f = gen2fam.get(syn.split(':',1)[1])
    elif syn:
        f = fam.get(syn)
    else:
        f = fam.get(av) or fam.get(sp)
    if f is None: missing.append(k); continue
    v = feats[k]
    g = list(v['trp_gaps']) + [0, 0, 0]
    X.append([v['trp_total'], v['trp_cterm150'], v['ident_erCry4'], g[0], g[1], g[2]])
    y.append(1 if v['group'] == 'migratory' else 0)
    F.append(f)
X = np.array(X, float); y = np.array(y, float); F = np.array(F)

def fit_logreg(Xtr, ytr, C=1.0, iters=200):
    n, p = Xtr.shape
    Xd = np.hstack([Xtr, np.ones((n, 1))])
    w = np.zeros(p + 1)
    lam = 1.0 / (C * n)
    for _ in range(iters):
        z = Xd @ w
        pr = 1 / (1 + np.exp(-z))
        g = Xd.T @ (pr - ytr) / n + lam * np.r_[w[:-1], 0.0]
        W = pr * (1 - pr) + 1e-9
        H = (Xd.T * W) @ Xd / n + lam * np.eye(p + 1)
        step = np.linalg.solve(H, g)
        w -= step
        if np.max(np.abs(step)) < 1e-8: break
    return w

def predict(Xte, w):
    Xd = np.hstack([Xte, np.ones((len(Xte), 1))])
    return 1 / (1 + np.exp(-(Xd @ w)))

def auc(sc, yy):
    pos = sc[yy == 1]; neg = sc[yy == 0]
    return float(np.mean([(p > neg).mean() + 0.5 * (p == neg).mean() for p in pos]))

def loo(X, y, groups=None):
    # groups: if None, each observation its own fold (LOO); else leave-one-group-out
    uniq = np.arange(len(y)) if groups is None else np.unique(groups)
    sc = np.full(len(y), np.nan)
    for u in uniq:
        te = (np.arange(len(y)) == u) if groups is None else (groups == u)
        tr = ~te
        if y[tr].sum() == 0 or y[tr].sum() == tr.sum():  # single-class train
            sc[te] = 0.5; continue
        w = fit_logreg(X[tr], y[tr])
        sc[te] = predict(X[te], w)
    return auc(sc, y), sc

# standardize features (attempt-2 Newton logistic on raw features; check LOO repro)
res = {'label': 'POST-HOC ROBUSTNESS (judge round-5 weakness 5); locked attempt-2 result unchanged',
       'n': int(len(y)), 'n_mig': int(y.sum()), 'missing_family_rows': missing}
# pipeline check on FULL 116-row panel (no family requirement)
Xa, ya = [], []
for k in sorted(feats):
    v = feats[k]
    g = list(v['trp_gaps']) + [0,0,0]
    Xa.append([v['trp_total'], v['trp_cterm150'], v['ident_erCry4'], g[0], g[1], g[2]])
    ya.append(1 if v['group']=='migratory' else 0)
Xa = np.array(Xa,float); ya = np.array(ya,float)
loo_full_auc, _ = loo(Xa, ya)
res['pipeline_check_LOO_auc_full116'] = round(loo_full_auc, 4)
res['n_full116_mig'] = int(ya.sum())
loo_auc, _ = loo(X, y)
res['pipeline_check_LOO_auc'] = round(loo_auc, 4)
res['locked_LOO_auc'] = 0.2066
lofo_auc, sc = loo(X, y, F)
res['LOFO_auc'] = round(lofo_auc, 4)
res['n_families'] = int(len(np.unique(F)))
res['families'] = sorted(set(F.tolist()))
res['conclusion'] = ('null survives phylogenetic blocking' if lofo_auc < 0.75 else 'NULL DOES NOT SURVIVE LOFO - investigate')
json.dump(res, open('results/h1_attempt2_lofo.json', 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('n','n_mig','pipeline_check_LOO_auc','locked_LOO_auc','LOFO_auc','n_families','missing_family_rows','conclusion')}, indent=1))
