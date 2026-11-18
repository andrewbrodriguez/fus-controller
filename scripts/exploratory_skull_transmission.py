"""Exploratory (11/18). Three PRE-SPECIFIED candidates (written before running):
 A: log(mean f0 / mean volt) added to log dose   (skull transmission proxy)
 B: log dose replaced by log(sum h2 / mean f0)   (receive-path-corrected dose)
 C: log dose replaced by log(sum h2 / mean f0^2)
Within-mouse LOMO RMSE on logit GFP+ fraction; null shuffles the candidate within mouse,
min over the three. Reads the burst cache + labels only."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'burst_ml')
from scipy.special import logit
from bml.data import load_labels, load_bursts, recordings_for
lab = load_labels(); Bu = load_bursts(); Bu = Bu[Bu.volt > 0.01]  # drop near-zero-drive baseline rows
rows = []
for _, r in lab.iterrows():
    m, t = int(r.mouse), int(r.target)
    b = pd.concat([Bu[(Bu.mouse == m) & (Bu.recording == x)] for x in recordings_for(m, t)])
    rows.append(dict(mouse=m, target=t, y=logit((r.k + .5) / (r.n + 1)), ld=np.log(r.dose),
                     A=np.log(b.f0.mean() / b.volt.mean()), B=np.log(r.dose / b.f0.mean()), C=np.log(r.dose / b.f0.mean() ** 2)))
d = pd.DataFrame(rows)
cen = lambda df, c: df[c] - df.groupby('mouse')[c].transform('mean')
def lomo(df, cols):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        b = np.linalg.lstsq(np.column_stack([cen(tr, c) for c in cols]), cen(tr, 'y'), rcond=None)[0]
        e += list((te.y - te.y.mean()) - np.column_stack([te[c] - te[c].mean() for c in cols]) @ b)
    return np.sqrt(np.mean(np.square(e)))
sets = {'dose': ['ld'], 'A: ld+f0/volt': ['ld', 'A'], 'B: dose/f0': ['B'], 'C: dose/f0^2': ['C']}
sc = {k: lomo(d, v) for k, v in sets.items()}; print({k: round(v, 3) for k, v in sc.items()})
rng = np.random.default_rng(0); base = sc['dose']; null = []
for _ in range(1000):
    dd = d.copy()
    for c in 'ABC': dd[c] = dd.groupby('mouse')[c].transform(lambda s: rng.permutation(s.values))
    null.append(max(base - lomo(dd, v) for k, v in sets.items() if k != 'dose'))
null = np.array(null)
for k in list(sets)[1:]: print(k, 'gain %.3f, p(best of 3) %.3f' % (base - sc[k], (null >= base - sc[k]).mean()))
