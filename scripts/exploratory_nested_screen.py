"""Exploratory (11/18): nested honest screen of ~40 candidate acoustic features added to log dose.
Outer: hold one mouse out. Inner: pick the best feature by LOMO RMSE among the other 5 mice only.
Reports outer RMSE of 'dose only' vs 'dose + inner-selected feature' (selection included in the score),
plus a within-mouse permutation null of the whole procedure. Reads CSVs only."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'burst_ml')
from scipy.special import logit
from bml.data import load_labels, load_bursts, recordings_for
lab = load_labels(); B = load_bursts()
rows = []
for _, r in lab.iterrows():
    m, t = int(r.mouse), int(r.target)
    b = pd.concat([B[(B.mouse == m) & (B.recording == x)] for x in recordings_for(m, t)])
    f = dict(mouse=m, target=t, y=logit((r.k + .5) / (r.n + 1)), ld=np.log(r.dose), nb=len(b))
    for band in ['bb', 'sub', 'ultra', 'wb']:
        x = b[band].values
        for th in (1.1, 1.2, 1.5):
            f[f'frac_{band}>{th}'] = np.mean(x > th)
        f[f'logsumex_{band}'] = np.log1p(np.sum(np.clip(x - 1, 0, None)))
        f[f'logmax_{band}'] = np.log(x.max())
        f[f'logmean_{band}'] = np.log(x.mean())
    h = b.h2.values; n = len(h)
    f['log_h2_first30'] = np.log(h[:30].sum() / 1e4 + 1e-6); f['log_h2_last30'] = np.log(h[-30:].sum() / 1e4 + 1e-6)
    f['log_h2_max'] = np.log(h.max()); f['log_h2_cv'] = np.log(h.std() / h.mean())
    f['log_h2_top10'] = np.log(np.sort(h)[-10:].sum() / 1e4)
    f['h2_wf_ratio'] = np.log(b.h2_wf.mean()); f['logvolt'] = np.log(b.volt.mean() + 1e-6)
    f['logvoltmax'] = np.log(b.volt.max() + 1e-6); f['lognb'] = np.log(n)
    f['log_f0'] = np.log(b.f0.mean())
    rows.append(f)
d = pd.DataFrame(rows)
cands = [c for c in d.columns if c not in ('mouse', 'target', 'y', 'ld', 'nb')]
cands = [c for c in cands if d[c].std() > 1e-9]
cen = lambda df, c: df[c] - df.groupby('mouse')[c].transform('mean')
def fit_pred(tr, te, cols):
    X = np.column_stack([cen(tr, c) for c in cols]); b = np.linalg.lstsq(X, cen(tr, 'y'), rcond=None)[0]
    Xt = np.column_stack([te[c] - te[c].mean() for c in cols]); return (te.y - te.y.mean()) - Xt @ b
def rmse_lomo(df, cols):
    e = np.concatenate([fit_pred(df[df.mouse != m], df[df.mouse == m], cols) for m in df.mouse.unique()])
    return np.sqrt(np.mean(e ** 2))
def nested(df):
    base, sel = [], []; picks = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        scores = {c: rmse_lomo(tr, ['ld', c]) for c in cands}
        best = min(scores, key=scores.get); picks.append(best)
        base.append(fit_pred(tr, te, ['ld'])); sel.append(fit_pred(tr, te, ['ld', best]))
    r = lambda e: np.sqrt(np.mean(np.concatenate(e) ** 2))
    return r(base), r(sel), picks
b, s, picks = nested(d)
print('n candidates', len(cands)); print('outer RMSE dose only %.3f, dose+inner-selected %.3f' % (b, s)); print('picks', picks)
full = sorted(((rmse_lomo(d, ['ld', c]), c) for c in cands))[:8]; print('top single (non-nested):', [(round(x, 3), c) for x, c in full])
rng = np.random.default_rng(1); gains = []
for i in range(100):
    dd = d.copy()
    for c in cands: dd[c] = dd.groupby('mouse')[c].transform(lambda s: rng.permutation(s.values))
    bb, ss, _ = nested(dd); gains.append(bb - ss)
print('observed gain %.3f; null mean %.3f, p = %.3f' % (b - s, np.mean(gains), np.mean(np.array(gains) >= b - s)))
