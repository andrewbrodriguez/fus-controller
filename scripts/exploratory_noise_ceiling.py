"""Exploratory (11/18): split-half noise ceiling of the delivery measure vs what log dose explains.
Per mouse, sections split into two halves (all splits for mice with >=4 sections, odd/even for
others); per target, logit GFP+ fraction summed over each half; within-mouse correlation between
halves -> Spearman-Brown reliability. Compare with within-mouse R^2 of log dose."""
import numpy as np, pandas as pd, itertools
from scipy.special import logit
s = pd.read_csv('results/histology/slices/slices_dose_delivery.csv')
print(s.columns.tolist())
s = s[(s.dose > 0) & (s.dropped.isna() | (s.dropped == False)) & s.n_cells.notna()]
rng = np.random.default_rng(0)
def half_logit(g):
    return logit((g.n_gfp_pos.sum() + .5) / (g.n_cells.sum() + 1))
rel, r2dose = [], []
for m, gm in s.groupby('mouse'):
    secs = gm.label.astype(str) + gm.section.astype(str)
    gm = gm.assign(sec=secs); u = gm.sec.unique()
    if len(u) < 2: continue
    rs = []
    for _ in range(200):
        p = rng.permutation(u); A, B = set(p[:len(u)//2]), set(p[len(u)//2:])
        a = gm[gm.sec.isin(A)].groupby('target').apply(half_logit); b = gm[gm.sec.isin(B)].groupby('target').apply(half_logit)
        t = a.index.intersection(b.index)
        if len(t) >= 4: rs.append(np.corrcoef(a[t], b[t])[0, 1])
    r = np.nanmean(rs); sb = 2 * r / (1 + r)
    full = gm.groupby('target').apply(half_logit); dose = gm.groupby('target').dose.first()
    t = full.index.intersection(dose.index); rd = np.corrcoef(np.log(dose[t]), full[t])[0, 1]
    print(f'mouse {m}: {len(u)} sections, split-half r {r:.2f}, Spearman-Brown reliability {sb:.2f}, dose r^2 {rd**2:.2f}')
    rel.append(sb); r2dose.append(rd**2)
print('mean reliability %.2f; mean dose r^2 %.2f' % (np.mean(rel), np.mean(r2dose)))
