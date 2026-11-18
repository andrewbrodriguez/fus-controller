"""Exploratory (11/18): does a per-position offset (learned from other mice) add to log dose
in held-out mice? Within-mouse LOMO RMSE on logit GFP+ fraction, plus permutation null
(position labels shuffled within mouse). Reads CSVs only."""
import numpy as np, pandas as pd
from scipy.special import logit
d = pd.read_csv('results/histology/slices/cumulative_dose_delivery.csv'); d = d[d.dose > 0].copy()
d['y'] = logit((d.gfp_cells + .5) / (d.cells + 1)); d['ld'] = np.log(d.dose)
def cen(df, c): return df[c] - df.groupby('mouse')[c].transform('mean')
def design(df, mode, pos=None):
    p = df.target if pos is None else pos
    cols = [cen(df, 'ld')]
    if mode == 'pos':
        for k in range(2, 7): cols.append(cen(df.assign(z=(p == k).astype(float)), 'z'))   # T1 baseline
    if mode == 'T1': cols.append(cen(df.assign(z=(p == 1).astype(float)), 'z'))
    if mode == 'side': cols.append(cen(df.assign(z=p.isin([1, 2, 3]).astype(float)), 'z'))
    return np.column_stack(cols)
def lomo(df, mode, pos=None):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        ptr = None if pos is None else pos[tr.index]; pte = None if pos is None else pos[te.index]
        X = design(tr, mode, ptr); b = np.linalg.lstsq(X, cen(tr, 'y'), rcond=None)[0]
        Xt = design(te, mode, pte); e += list((te.y - te.y.mean()) - Xt @ b)
    return np.sqrt(np.mean(np.square(e)))
base = lomo(d, 'ld')
res = {m: lomo(d, m) for m in ['T1', 'side', 'pos']}
print('dose only %.3f' % base, {k: round(v, 3) for k, v in res.items()})
rng = np.random.default_rng(0)
for mode in res:
    gains = []
    for _ in range(500):
        pos = d.groupby('mouse').target.transform(lambda s: pd.Series(rng.permutation(s.values), index=s.index))
        gains.append(base - lomo(d, mode, pos))
    print(mode, 'gain %.3f, perm p %.3f' % (base - res[mode], np.mean(np.array(gains) >= base - res[mode])))
