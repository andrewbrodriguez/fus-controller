"""Pre-registered test in burst_ml/plan_events.md. Run once."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'burst_ml')
from scipy.special import logit
from bml.data import load_labels, recordings_for
lab = load_labels(); C = pd.read_csv('burst_ml/cache/events.csv'); F = ['k2', 'f2', 'p2', 'ku', 'fu']
rows = []
for _, r in lab.iterrows():
    m, t = int(r.mouse), int(r.target)
    b = pd.concat([C[(C.mouse == m) & (C.recording == x)] for x in recordings_for(m, t)])
    rows.append(dict(mouse=m, target=t, y=logit((r.k + .5) / (r.n + 1)), ld=np.log(r.dose), **{f: b[f].mean() for f in F}))
d = pd.DataFrame(rows)
cen = lambda df, c: df[c] - df.groupby('mouse')[c].transform('mean')
def lomo(df, cols):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        b = np.linalg.lstsq(np.column_stack([cen(tr, c) for c in cols]), cen(tr, 'y'), rcond=None)[0]
        e += list((te.y - te.y.mean()) - np.column_stack([te[c] - te[c].mean() for c in cols]) @ b)
    return np.sqrt(np.mean(np.square(e)))
base = lomo(d, ['ld']); sc = {f: lomo(d, ['ld', f]) for f in F}; print('dose %.3f' % base, {k: round(v, 3) for k, v in sc.items()})
rng = np.random.default_rng(0); null = []
for _ in range(1000):
    dd = d.copy()
    for f in F: dd[f] = dd.groupby('mouse')[f].transform(lambda s: rng.permutation(s.values))
    b0 = lomo(dd, ['ld']); null.append(max(b0 - lomo(dd, ['ld', f]) for f in F))
null = np.array(null)
for f in F: print(f, 'gain %.3f, p(best of 5) %.3f' % (base - sc[f], (null >= base - sc[f]).mean()))
