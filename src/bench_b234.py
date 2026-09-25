"""W02 benchmarks B2/B3/B4 per results/benchmarks.md (locked 20:20 IST pre-computation)."""
from bench_common import *

def spike_amp(y, smooth):
    return float((y - smooth).max())

# --- B2 ---
res = {'B2': {}, 'B3': {}, 'B4': {}}
def smooth_bg(y):
    # sinusoidal background estimate: fit a + b*cos(2t)
    A = np.vstack([np.ones_like(THETAS), np.cos(2*THETAS)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return A @ coef

for tau_us, k in [(1.0, 1.0), (10.0, 0.1), (100.0, 0.01)]:
    rp = RadicalPair(hyperfine_A=[N5, N10], hyperfine_B=[Y45], kS=k, kT=k)
    y = curve(rp)
    res['B2'][f'tau_{tau_us}us'] = {'phi_max': float(y.max()), 'phi_min': float(y.min()),
        'spike_amp_over_background': spike_amp(y, smooth_bg(y))}
    print('B2 tau', tau_us, res['B2'][f'tau_{tau_us}us'], flush=True)
rp0 = RadicalPair(hyperfine_A=[N5, N10], hyperfine_B=[], kS=1.0, kT=1.0)
y0 = curve(rp0)
res['B2']['no_partner_anisotropy'] = {'spike_amp_over_background': spike_amp(y0, smooth_bg(y0))}
res['B2']['E2a_spike_present'] = res['B2']['tau_1.0us']['spike_amp_over_background'] > 3 * res['B2']['no_partner_anisotropy']['spike_amp_over_background']
res['B2']['E2b_grows_with_tau'] = (res['B2']['tau_100.0us']['spike_amp_over_background'] > res['B2']['tau_1.0us']['spike_amp_over_background'])
res['B2']['E2c_no_spike_without_partner'] = res['B2']['no_partner_anisotropy']['spike_amp_over_background'] < 0.2 * res['B2']['tau_1.0us']['spike_amp_over_background']
res['B2']['pass'] = all([res['B2']['E2a_spike_present'], res['B2']['E2b_grows_with_tau'], res['B2']['E2c_no_spike_without_partner']])
print('B2 pass:', res['B2']['pass'], flush=True)

# --- B3 ---
def zeroT(T):
    T2 = T.copy(); T2[0,0] = 0; T2[1,1] = 0; return T2
def b3_amp(tens):
    y = curve(RadicalPair(hyperfine_A=tens, hyperfine_B=[Y45], kS=1.0, kT=1.0))
    return spike_amp(y, smooth_bg(y))
base = b3_amp([N5,N10])
att = {}
for name, tens in [('N5_zeroed', [zeroT(N5), N10]), ('N10_zeroed', [N5, zeroT(N10)]), ('both_zeroed', [zeroT(N5), zeroT(N10)])]:
    att[name] = b3_amp(tens)
res['B3'] = {'baseline_spike_amp': base, 'attenuated': att,
    'E3a_N5_attenuation_pct': round(100*(1-att['N5_zeroed']/base),1) if base>0 else None,
    'E3a_N10_attenuation_pct': round(100*(1-att['N10_zeroed']/base),1) if base>0 else None,
    'E3b_both_abolishes': att['both_zeroed'] < 0.05*base if base>0 else None}
res['B3']['pass'] = bool(res['B3']['E3b_both_abolishes']) and any(
    30 <= res['B3'][k] <= 100 for k in ('E3a_N5_attenuation_pct','E3a_N10_attenuation_pct') if res['B3'][k] is not None)
print('B3:', res['B3'], flush=True)

# --- B4 ---
A1 = np.diag([-0.0989, -0.0989, 1.7569])  # Ritz-style single anisotropic nucleus
rp = RadicalPair(hyperfine_A=[A1], hyperfine_B=[], kS=1.0, kT=1.0)
anis = {}
for B in [0.005, 0.05, 0.5, 5.0, 50.0]:
    y = curve(rp, B=B)
    anis[str(B)] = float((y.max()-y.min())/y.mean())
    print('B4', B, anis[str(B)], flush=True)
vals = list(anis.values())
high_field = vals[2:]
res['B4'] = {'anisotropy_vs_B_mT': anis,
    'E4a_high_field_monotone_sat': all(high_field[i+1] >= high_field[i] - 0.05*high_field[i] for i in range(len(high_field)-1)),
    'E4b_earth_below_sat': anis['0.05'] < anis['5.0']}
res['B4']['pass'] = res['B4']['E4a_high_field_monotone_sat'] and res['B4']['E4b_earth_below_sat']
res['summary'] = {'B2': res['B2']['pass'], 'B3': res['B3']['pass'], 'B4': res['B4']['pass']}
json.dump(res, open('/home/sandbox/weird10/projects/W02/results/bench_b234.json','w'), indent=1)
print(json.dumps(res['summary']))
