"""Pre-registered test in burst_ml/plan_collapse.md. Run once."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'burst_ml')
from scipy.special import logit
from bml.data import load_labels, recordings_for
lab = load_labels(); C = pd.read_csv('burst_ml/cache/collapse.csv')
rows = []
for _, r in lab.iterrows():
    m, t = int(r.mouse), int(r.target)
    b = pd.concat([C[(C.mouse == m) & (C.recording == x)] for x in recordings_for(m, t)])
    rows.append(dict(mouse=m, target=t, y=logit((r.k + .5) / (r.n + 1)), ld=np.log(r.dose), c2=b.c2.mean(), r2=b.r2.mean()))
d = pd.DataFrame(rows); print(d[['c2', 'r2']].describe().round(3))
cen = lambda df, c: df[c] - df.groupby('mouse')[c].transform('mean')
def lomo(df, cols):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        b = np.linalg.lstsq(np.column_stack([cen(tr, c) for c in cols]), cen(tr, 'y'), rcond=None)[0]
        e += list((te.y - te.y.mean()) - np.column_stack([te[c] - te[c].mean() for c in cols]) @ b)
    return np.sqrt(np.mean(np.square(e)))
sets = {'ld': ['ld'], 'ld+c2': ['ld', 'c2'], 'ld+r2': ['ld', 'r2'], 'ld+c2+r2': ['ld', 'c2', 'r2']}
sc = {k: lomo(d, v) for k, v in sets.items()}; print({k: round(v, 3) for k, v in sc.items()})
rng = np.random.default_rng(0); null = []
for _ in range(1000):
    dd = d.copy()
    for c in ('c2', 'r2'): dd[c] = dd.groupby('mouse')[c].transform(lambda s: rng.permutation(s.values))
    b = lomo(dd, ['ld']); null.append(max(b - lomo(dd, v) for k, v in sets.items() if k != 'ld'))
null = np.array(null)
for k in list(sets)[1:]: print(k, 'gain %.3f, p(best of 3) %.3f' % (sc['ld'] - sc[k], (null >= sc['ld'] - sc[k]).mean()))
