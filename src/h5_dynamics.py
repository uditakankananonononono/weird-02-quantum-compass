#!/usr/bin/env python3
"""W02 queue #5 dynamics arms (amendment locked 16:20 IST, commit e81f26b).
ARM 5a: metapredict per-residue disorder on committed Cry4/Cry1/Cry2 panels.
ARM 5b: ProDy GNM flexibility on committed AFDB models. Report-only; MD deferred."""
import sys, os, json, glob, hashlib, subprocess, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import re
from Bio import SeqIO

def load_seqs():
    panels = {}
    panels['Cry4'] = [(r.id, str(r.seq)) for r in SeqIO.parse('data/cry4_subpanel.fasta', 'fasta')]
    panels['Cry1'] = [(r.id, str(r.seq).replace('-', '')) for r in SeqIO.parse('data/cry1_panel_msa.fasta', 'fasta')]
    panels['Cry2'] = []
    for f in sorted(glob.glob('data/cry2_panel/*.fa')):
        for r in SeqIO.parse(f, 'fasta'):
            panels['Cry2'].append((r.id, str(r.seq).replace('-', '')))
    return panels

def arm5a():
    import metapredict as meta
    ver = subprocess.run(['pip', 'show', 'metapredict'], capture_output=True, text=True).stdout
    ver = [l.split(':')[1].strip() for l in ver.splitlines() if l.startswith('Version')][0]
    panels = load_seqs()
    out = {'tool': f'metapredict {ver}', 'x_handling': 'clarification 16:29 IST: X-free segments >=10aa, X=NaN, affected: sedentary|Nothoprocta_perdicaria', 'per_sequence': {}, 'summaries': {}}
    for pan, recs in panels.items():
        rows = []
        for rid, seq in recs:
            # clarification 16:29 IST: X-free segments, X -> NaN, segments <10 skipped
            full = np.full(len(seq), np.nan)
            for m in re.finditer(r'[^Xx]+', seq):
                seg = m.group(0)
                if len(seg) >= 10:
                    full[m.start():m.end()] = meta.predict_disorder(seg)
            s = full
            tail = s[-60:] if len(s) >= 60 else s
            row = {'id': rid, 'len': len(s), 'disordered_fraction': round(float(np.nanmean((s > 0.5).astype(float))), 4),
                   'longest_stretch': int(max((len(list(g)) for k, g in __import__('itertools').groupby(np.nan_to_num(s, nan=0) > 0.5) if k), default=0)),
                   'tail60_disorder_fraction': round(float(np.nanmean((tail > 0.5).astype(float))), 4)}
            out['per_sequence'][f'{pan}|{rid}'] = row
            rows.append(row)
        out['summaries'][pan] = {
            'n': len(rows),
            'disordered_fraction_median': round(float(np.median([r['disordered_fraction'] for r in rows])), 4),
            'longest_stretch_median': float(np.median([r['longest_stretch'] for r in rows])),
            'tail60_median': round(float(np.median([r['tail60_disorder_fraction'] for r in rows])), 4)}
    from scipy.stats import mannwhitneyu
    out['group_comparisons'] = {}
    for met in ('disordered_fraction', 'longest_stretch', 'tail60_disorder_fraction'):
        groups = {pan: [out['per_sequence'][f'{pan}|{rid}'][met] for rid, _ in recs] for pan, recs in panels.items()}
        for other in ('Cry1', 'Cry2'):
            u = mannwhitneyu(groups['Cry4'], groups[other], alternative='two-sided')
            out['group_comparisons'][f'Cry4_vs_{other}|{met}'] = {'U': float(u.statistic), 'p_report_only': float(u.pvalue)}
    return out

def arm5b():
    import prody
    models = json.load(open('data/afdb_reference_models.json'))
    os.makedirs('data/afdb_models_v6', exist_ok=True)
    ledger = {}
    out = {'tool': f'prody {prody.__version__}', 'per_model': {}, 'failed_downloads': []}
    for uid, m in models.items():
        fn = f'data/afdb_models_v6/AF-{uid}-F1-model_v6.pdb'
        if not os.path.exists(fn):
            url = f'https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.pdb'
            r = subprocess.run(['curl', '-sf', '-o', fn, url])
            if r.returncode != 0:
                out['failed_downloads'].append(uid); continue
            ledger[uid] = hashlib.sha256(open(fn, 'rb').read()).hexdigest()
        ag = prody.parsePDB(fn).select('calpha')
        gnm = prody.GNM('gnm'); gnm.buildKirchhoff(ag, cutoff=10.0)
        gnm.calcModes(n_modes=10)
        msf = prody.calcSqFlucts(gnm[:10])
        n = len(msf); tail = msf[-60:] if n >= 60 else msf
        core = msf[:-60] if n >= 60 else msf
        out['per_model'][uid] = {'name': m.get('name'), 'n_ca': n,
                                 'tail60_msf_mean': round(float(tail.mean()), 4),
                                 'core_msf_mean': round(float(core.mean()), 4),
                                 'tail_core_ratio': round(float(tail.mean() / core.mean()), 4) if core.mean() > 0 else None}
    json.dump(ledger, open('data/afdb_models_v6/sha256_ledger.json', 'w'), indent=1)
    return out

if __name__ == '__main__':
    which = sys.argv[1]
    if which == '5a':
        res = arm5a(); fn = 'results/h5a_disorder.json'
    else:
        res = arm5b(); fn = 'results/h5b_flexibility.json'
    res['amendment_commit'] = 'e81f26b'
    json.dump(res, open(fn, 'w'), indent=1)
    print(which, '->', fn, flush=True)
