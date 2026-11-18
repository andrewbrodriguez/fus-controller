"""Exploratory (11/18), ONE pre-stated candidate: dt = minutes from each target's (last) sonication to the
last sonication of its session (AAV goes in after the last target); added to log dose, within-mouse LOMO RMSE;
null = shuffle dt within mouse. Motivated by PubMed 41397344 (AAV timing after FUS mattered; harmonic dose did not)."""
import numpy as np, pandas as pd
from scipy.special import logit
d = pd.read_csv('results/histology/slices/cumulative_dose_delivery.csv'); d = d[d.dose > 0].copy()
rec = pd.read_csv('results/acoustic_eda/recordings.csv', parse_dates=['time']); rec = rec[rec.mouse <= 6]
g = rec.groupby(['mouse', 'position']).time.max().reset_index().rename(columns={'position': 'target'})
g['dt'] = (g.groupby('mouse').time.transform('max') - g.time).dt.total_seconds() / 60
d = d.merge(g[['mouse', 'target', 'dt']], on=['mouse', 'target']).reset_index(drop=True)
d['y'] = logit((d.gfp_cells + .5) / (d.cells + 1)); d['ld'] = np.log(d.dose)
cen = lambda df, c: df[c] - df.groupby('mouse')[c].transform('mean')
def lomo(df, cols):
    e = []
    for m in df.mouse.unique():
        tr, te = df[df.mouse != m], df[df.mouse == m]
        b = np.linalg.lstsq(np.column_stack([cen(tr, c) for c in cols]), cen(tr, 'y'), rcond=None)[0]
        e += list((te.y - te.y.mean()) - np.column_stack([te[c] - te[c].mean() for c in cols]) @ b)
    return np.sqrt(np.mean(np.square(e)))
a, b = lomo(d, ['ld']), lomo(d, ['ld', 'dt']); rng = np.random.default_rng(0); n = []
for _ in range(2000):
    dd = d.copy(); dd['dt'] = dd.groupby('mouse').dt.transform(lambda s: rng.permutation(s.values)); n.append(a - lomo(dd, ['ld', 'dt']))
print('dt range min', d.dt.min(), d.dt.max(), '| dose %.3f, dose+dt %.3f, gain %.3f, p %.3f' % (a, b, a - b, np.mean(np.array(n) >= a - b)))
c = lambda k: cen(d, k).values; print('within-mouse coef on dt (logit per min):', np.linalg.lstsq(np.column_stack([c('ld'), c('dt')]), c('y'), rcond=None)[0].round(3))
