#!/usr/bin/env python3
"""W02 queue #8 (12:31 IST amendment): generic ML vs committed physics descriptors.
LR + RF, seed 260927, LOO CV, AUC vs attempt8c MDE floor 0.6424. Report-only.
Feature set: the committed v1 physics descriptors used by the attempt-3a logistic benchmark
(tail_after_ref, charge_pH7, pI, aromatic_frac, kd_mean, length, ident_erCry1, trp_total),
the set the 8c classifier_auc MDE (0.6424) is calibrated on."""
import json, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score

FEATS=['tail_after_ref','charge_pH7','pI','aromatic_frac','kd_mean','length','ident_erCry1','trp_total']
SEED=260927; MDE=0.6424
d=json.load(open('data/cry1_features_v1.json'))
keys=sorted(d)
X=np.array([[d[k][f] for f in FEATS] for k in keys],float)
y=np.array([1 if d[k]['group']=='migratory' else 0 for k in keys])
loo=LeaveOneOut()
def loo_probs(clf_factory):
    p=np.zeros(len(y))
    for tr,te in loo.split(X):
        clf=clf_factory(); clf.fit(X[tr],y[tr])
        p[te]=clf.predict_proba(X[te])[:,1]
    return p
res={}
for name,factory in [
  ('logistic_regression', lambda: LogisticRegression(C=1.0, max_iter=5000, random_state=SEED)),
  ('random_forest', lambda: RandomForestClassifier(n_estimators=100, random_state=SEED))]:
    p=loo_probs(factory)
    auc=float(roc_auc_score(y,p))
    res[name]={'loo_auc':round(auc,4),
      'vs_mde_floor_0.6424':'BELOW detectable-effect bound' if auc<MDE else 'at/above MDE floor',
      'delta_vs_floor':round(auc-MDE,4)}
out={
 'amendment':'2026-09-27 12:31 IST queue #8 (delegated executor, report-only)',
 'executor_note':'computed by delegated executor at ref 4845c5e; uncommitted, returned to lane-29 for review',
 'features':FEATS,
 'feature_source':'data/cry1_features_v1.json (frozen v1 physics descriptor matrix; the attempt-3a feature set the attempt8c classifier_auc MDE is calibrated on)',
 'panel':{'n':int(len(y)),'n_migratory':int(y.sum()),'n_sedentary':int((1-y).sum()),
   'note':'116-species locked Cry1 panel, labels from locked panel list (same n_mig=36/n_sed=80 as 7f-b/8c)'},
 'seed':SEED,'evaluation':'leave-one-out CV, held-out probabilities pooled, single AUC',
 'mde_floor':MDE,'mde_source':'results/attempt8c_power.json classifier_auc.MDE_AUC (attempt 3a LOO AUC vs 0.5, observed_best_auc 0.402)',
 'reference_committed_result':{'attempt_3a_logistic_loo_auc':0.402},
 'results':res,
 'plm_embedding_arm':'UNAVAILABLE - no local PLM embeddings tool in executor env (no esm/torch local weights); multi-GB model downloads barred by protocol. Arm disclosed, not run.',
 'interpretation_lock':'REPORT-ONLY. This benchmark bounds descriptor information content only. Any AUC below 0.6424 is below the panel detectable-effect bound. LOCKED ANTI-GOAL: no migratory-vs-sedentary prediction is claimed as a project result.'}
json.dump(out, open('results/attempt8g_ml_benchmark.json','w'), indent=1)
print(json.dumps(res,indent=1))
print('n=',len(y),'mig=',y.sum())
