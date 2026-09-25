#!/usr/bin/env python3
"""H1 scoring per locked pre-reg + M1/M1.2. Fires only when folds complete.
Locked: AUC >= 0.75, migratory(3)=positive vs sedentary(1), anisotropy@50uT from
M1.2-parameterized validated solver. Held-out partial(2) reported separately (R2).
"""
import json, os, sys, subprocess
sys.path.insert(0, 'src')
os.chdir('/home/sandbox/weird10/projects/W02')

def run_scoring():
    aln = {}
    for line in open('data/cry1_panel_msa.fasta'):
        if line.startswith('>'): k = line[1:].strip(); aln[k] = ''
        else: aln[k] += line.strip()
    rows = []
    fails = []
    for rowname in aln:
        if rowname.startswith('REF|'): continue
        grp = rowname.split('|')[0]
        pdb = 'data/folds/' + rowname.replace('|','__') + '.pdb'
        if not os.path.exists(pdb):
            fails.append((rowname, 'no fold')); continue
        open('/tmp/hs.fasta','w').write(f'>{rowname}\n{aln[rowname].replace("-","")}\n')
        r = subprocess.run(['python3','src/compass_scan.py','/tmp/hs.fasta','--pdb',pdb,
                            '--msa-row',rowname,'--json'], capture_output=True, text=True, timeout=300)
        try:
            out = json.loads(r.stdout)
            rows.append({'row': rowname, 'group': grp, 'A': out['anisotropy_50uT'],
                         'tau': out['tau_eff_us'], 'd_term': out['d_term_A'],
                         'terminal_pos': out['terminal_chain_pos'], 'broken': out['chain_broken']})
        except Exception:
            fails.append((rowname, (r.stderr or r.stdout)[-160:]))
    def auc(pos, neg):
        allx = [(x,1) for x in pos] + [(x,0) for x in neg]
        allx.sort(key=lambda t: t[0])
        rank_sum = sum(i+1 for i,(x,y) in enumerate(allx) if y==1)
        n1, n0 = len(pos), len(neg)
        return (rank_sum - n1*(n1+1)/2) / (n1*n0)
    mig = [r['A'] for r in rows if r['group']=='migratory']
    sed = [r['A'] for r in rows if r['group']=='sedentary']
    par = [r['A'] for r in rows if r['group']=='partial']
    result = {'n_mig': len(mig), 'n_sed': len(sed), 'n_partial': len(par),
              'AUC_mig_vs_sed': round(auc(mig, sed), 4) if mig and sed else None,
              'gate': 'AUC >= 0.75 (locked)',
              'median_A_mig': sorted(mig)[len(mig)//2] if mig else None,
              'median_A_sed': sorted(sed)[len(sed)//2] if sed else None,
              'partial_median_A': sorted(par)[len(par)//2] if par else None,
              'fails': fails, 'rows': rows}
    json.dump(result, open('results/h1_scoring.json','w'), indent=1)
    print('AUC:', result['AUC_mig_vs_sed'], 'n:', len(mig), len(sed), 'fails:', len(fails))
    print('medians mig/sed:', result['median_A_mig'], result['median_A_sed'])

if __name__ == '__main__':
    run_scoring()
