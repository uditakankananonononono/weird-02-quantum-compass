#!/usr/bin/env python3
"""W02 7g-a SENSITIVITY (locked in 2026-09-27 10:13 amendment: matched-intersection
species sub-MSAs, identical pipeline). Reuses the executed 7g module verbatim so the
pipeline cannot diverge; only the row set changes (species present in BOTH panels)."""
import sys, json
sys.path.insert(0, 'src')
import attempt7g_ab as g  # executes the locked full-panel run identically (deterministic, same seed)

sp1 = {k.split('|')[1] for k in g.msa1}
sp2 = set(g.msa2.keys())
inter = sp1 & sp2
print('intersection species:', len(inter), flush=True)
m1 = {k: v for k, v in g.msa1.items() if k.split('|')[1] in inter}
m2 = {k: v for k, v in g.msa2.items() if k in inter}
qc1 = g.qcols(m1, g.L1)
qc2 = g.qcols(m2, g.L2)
r1 = g.mi_apc(m1, qc1, g.col2res1, g.SRP1, g.surf1, 'CRY1-intersect')
print('CRY1 7g-a sens:', {k: v for k, v in r1.items() if k != 'C_residues'}, flush=True)
r2 = g.mi_apc(m2, qc2, g.col2res2, g.SRP2, g.surf2, 'CRY2-intersect')
print('CRY2 7g-a sens:', {k: v for k, v in r2.items() if k != 'C_residues'}, flush=True)
win = r1['C_size'] >= 10 and r1['C_surface'] >= 3 and (r2['C_size'] < 10 or r2['C_surface'] < 3)
out = {'amendment': '2026-09-27 10:13 IST 7g-a sensitivity (matched-intersection sub-MSAs, identical pipeline)',
       'n_intersection_species': len(inter), 'CRY1': r1, 'CRY2': r2, 'win': bool(win)}
json.dump(out, open('results/attempt7g_a_sensitivity.json', 'w'), indent=1)
print('7g-a sensitivity WIN:', win, flush=True)
