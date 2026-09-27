#!/usr/bin/env python3
"""W02 queue #6 confidence comparison (amendment 2026-09-27 11:27 IST, locked pre-computation)."""
import json, glob, sys
import numpy as np
sys.path.insert(0, 'src')

geo = json.load(open('data/attempt6_column_geometry.json'))
chain6 = set(json.load(open('data/trp_chain_definition.json'))['chain_6ptz_numbering'])
core_cols = sorted(int(c) for c, g in geo.items()
                   if g['resnum_6ptz'] in chain6 or g['fad_d'] <= 4.5 or g['trp_d'] <= 6.0)

def load_msa(path):
    msa = {}
    for line in open(path):
        if line.startswith('>'): k = line[1:].strip(); msa[k] = ''
        else: msa[k] += line.strip()
    return msa

def col_to_residx(row, cols):
    m = {}
    idx = 0
    cols = set(cols)
    for i, ch in enumerate(row):
        if ch != '-':
            idx += 1
            if i in cols: m[i] = idx
    return m  # column -> 1-based residue index

def parse_pdb(path):
    ca = {}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            ca[int(line[22:26])] = float(line[60:66])
    return ca

def stats(vals, thr):
    vals = list(vals)
    if not vals: return None
    a = np.array(vals)
    return {'n': len(vals), 'mean': round(float(a.mean()), 3), 'frac_ge_thr': round(float((a >= thr).mean()), 3)}

msa1 = load_msa('data/cry1_panel_msa.fasta')
msa4 = load_msa('data/mafft_ebi_crosscheck.fasta')

folds = {}
for p in sorted(glob.glob('data/folds/*.pdb')):
    name = p.split('/')[-1][:-4]
    folds[name] = parse_pdb(p)

def scale_thr(ca):  # per-file detection (clarification 11:28)
    return 0.70 if float(np.median(list(ca.values()))) <= 1.5 else 70.0

def region_stats(ca, core_idx):
    if not ca: return None
    thr = scale_thr(ca)
    L = max(ca)
    tail_lo = int(0.75 * L) + 1
    core_b = [b for r, b in ca.items() if r in core_idx]
    tail_b = [b for r, b in ca.items() if r >= tail_lo and r not in core_idx]
    rest_b = [b for r, b in ca.items() if r < tail_lo and r not in core_idx]
    return {'len': L, 'thr': thr, 'overall_mean': round(float(np.mean(list(ca.values()))), 3),
            'CORE': stats(core_b, thr), 'TAIL': stats(tail_b, thr), 'REST': stats(rest_b, thr)}

out = {'amendments': ['2026-09-27 11:27 IST queue #6', '2026-09-27 11:28 IST clarification'],
       'scale_note': 'per-file detection; ESMFold 0-1 (thr 0.70), AFDB 0-100 (thr 70)',
       'n_folds': len(folds), 'classes': {}}
classes = {}
for name, ca in folds.items():
    if name.startswith(('migratory__', 'sedentary__', 'REF__')):
        rowkey = name.replace('__', '|')
        row = msa1.get(rowkey)
        if row is None:
            # REF rows: fold REF__sp__Q5IZC5__CRY1_ERIRU -> row REF|sp|Q5IZC5|CRY1_ERIRU
            classes.setdefault('unmatched', []).append(name); continue
        m = col_to_residx(row, core_cols)
        core_idx = set(m.values())
        cls = 'cry1_panel' if name.startswith(('migratory__', 'sedentary__')) else 'ref_rows'
    elif name.startswith('CRY4__'):
        acc = name.split('__')[3]
        row = msa4.get(acc)
        if row is None: classes.setdefault('unmatched', []).append(name); continue
        m = col_to_residx(row, core_cols)
        core_idx = set(m.values())
        cls = 'cry4'
    elif name.startswith('CRY2__'):
        import attempt8_redirect3 as r3
        core6 = sorted({geo[c]['resnum_6ptz'] for c in map(str, core_cols)})
        core_idx = {r3.m6to2[r] for r in core6 if r in r3.m6to2}
        cls = 'cry2_ref'
    elif name.startswith('AFDB__'):
        uid = name.split('__')[1]
        partner = 'REF__sp__' + uid + '__' + None if False else None
        # find REF fold with same UniProt id
        pref = [n for n in folds if n.startswith('REF__') and n.split('__')[2] == uid]
        if pref:
            rowkey = pref[0].replace('__', '|')
            row = msa1.get(rowkey)
            if row is not None and len([c for c in row if c != '-']) == max(folds[name]):
                m = col_to_residx(row, core_cols)
                core_idx = set(m.values())
            else:
                core_idx = set()
        else:
            core_idx = set()
        cls = 'afdb'
    else:
        classes.setdefault('unmatched', []).append(name); continue
    st = region_stats(ca, core_idx)
    st['core_mapped_residues'] = len(core_idx)
    classes.setdefault(cls, {})[name] = st

for cls, d in classes.items():
    if cls == 'unmatched':
        out['classes'][cls] = d; continue
    oms = [v['overall_mean'] for v in d.values()]
    cm = [v['CORE']['mean'] for v in d.values() if v['CORE']]
    cf = [v['CORE']['frac_ge_thr'] for v in d.values() if v['CORE']]
    tm = [v['TAIL']['mean'] for v in d.values() if v['TAIL']]
    out['classes'][cls] = {
        'n_folds': len(d),
        'overall_mean_median': round(float(np.median(oms)), 3),
        'overall_mean_min': round(float(np.min(oms)), 3),
        'core_mean_median': round(float(np.median(cm)), 3) if cm else None,
        'core_mean_min': round(float(np.min(cm)), 3) if cm else None,
        'core_frac_ge_thr_median': round(float(np.median(cf)), 3) if cf else None,
        'tail_mean_median': round(float(np.median(tm)), 3) if tm else None,
        'per_fold': {k: v for k, v in d.items() if cls in ('cry4', 'cry2_ref', 'afdb', 'ref_rows')},
    }
# paired REF vs AFDB
pairs = []
for name in folds:
    if name.startswith('AFDB__'):
        uid = name.split('__')[1]
        pref = [n for n in folds if n.startswith('REF__') and n.split('__')[2] == uid]
        if pref and pref[0] in classes.get('ref_rows', {}) and name in classes.get('afdb', {}):
            r, a = classes['ref_rows'][pref[0]], classes['afdb'][name]
            if r['CORE'] and a['CORE']:
                pairs.append({'uid': uid, 'ref_core_mean': r['CORE']['mean'], 'afdb_core_mean': a['CORE']['mean'],
                              'ref_overall': r['overall_mean'], 'afdb_overall': a['overall_mean']})
out['ref_vs_afdb_pairs'] = pairs
census = json.load(open('data/rcsb_cry_census_details.json'))
res = [v['resolution'][0] for v in census.values() if v.get('resolution')]
out['census'] = {'n_solved': len(census), 'methods': sorted({v.get('method', '?') for v in census.values()}),
                 'resolution_median': round(float(np.median(res)), 2) if res else None,
                 'resolution_range': [min(res), max(res)] if res else None,
                 '6ptz': census.get('6PTZ', {})}
json.dump(out, open('results/attempt8d_confidence.json', 'w'), indent=1)
print(json.dumps({k: (v if k != 'classes' else {c: {kk: vv for kk, vv in cl.items() if kk != 'per_fold'} for c, cl in v.items()}) for k, v in out.items()}, indent=1, default=str)[:2600])
