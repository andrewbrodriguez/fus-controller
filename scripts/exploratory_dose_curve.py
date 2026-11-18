"""Exploratory (11/18). Pre-stated dose-curve shapes vs straight line in log dose, within-mouse LOMO RMSE:
 Q: quadratic in log dose; L: linear in dose; I: isotonic (monotone) in log dose. Null: shuffle y within mouse, best of 3."""
import numpy as np, pandas as pd
from scipy.special import logit
from sklearn.isotonic import IsotonicRegression
d = pd.read_csv('results/histology/slices/cumulative_dose_delivery.csv'); d = d[d.dose > 0].reset_index(drop=True)
d['y'] = logit((d.gfp_cells + .5) / (d.cells + 1)); d['ld'] = np.log(d.dose)
def cen(df, c): return df[c] - df.groupby('mouse')[c].transform('mean')
def lomo(df, kind):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        yc, yt = cen(tr, 'y'), te.y - te.y.mean()
        if kind == 'I':
            x = cen(tr, 'ld'); iso = IsotonicRegression(out_of_bounds='clip').fit(x, yc)
            e += list(yt - iso.predict(te.ld - te.ld.mean())); continue
        f = {'line': lambda t: [cen(t, 'ld')], 'Q': lambda t: [cen(t, 'ld'), cen(t.assign(q=t.ld ** 2), 'q')],
             'L': lambda t: [cen(t, 'dose')]}[kind]
        b = np.linalg.lstsq(np.column_stack(f(tr)), yc, rcond=None)[0]
        Xt = np.column_stack([c - c.mean() if False else c for c in [ (te.ld - te.ld.mean()), (te.ld**2 - te.ld.pow(2).mean())]][:len(b)]) if kind in ('line', 'Q') else np.column_stack([te.dose - te.dose.mean()])
        e += list(yt - Xt @ b)
    return np.sqrt(np.mean(np.square(e)))
sc = {k: lomo(d, k) for k in ['line', 'Q', 'L', 'I']}; print({k: round(v, 3) for k, v in sc.items()})
rng = np.random.default_rng(0); null = []
for _ in range(500):
    dd = d.copy(); dd['y'] = dd.groupby('mouse').y.transform(lambda s: rng.permutation(s.values))
    b = lomo(dd, 'line'); null.append(max(b - lomo(dd, k) for k in 'QLI'))
null = np.array(null)
for k in 'QLI': print(k, 'gain %.3f, p(best of 3) %.3f' % (sc['line'] - sc[k], (null >= sc['line'] - sc[k]).mean()))
