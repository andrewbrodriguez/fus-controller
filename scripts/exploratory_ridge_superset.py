"""Exploratory (11/18): fair superset test. Dose (log) is unpenalised; ALL candidate features of
exploratory_nested_screen.py enter with a ridge penalty chosen by inner leave-one-mouse-out on the
training mice only (grid includes a penalty so large that the model is dose alone).
Also in-sample fit (must be >= dose alone) and outer held-out RMSE at several fixed penalties."""
import sys, io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    sys.argv = ['x']; exec(open('scripts/exploratory_nested_screen.py').read().split("b, s, picks = nested(d)")[0])
F = cands + ['logvolt']; F = list(dict.fromkeys(F))
def prep(tr, te):
    mu = tr[F].groupby(tr.mouse).transform('mean'); sd = (tr[F] - mu).std().replace(0, 1)
    Ztr = ((tr[F] - mu) / sd).values; Zte = ((te[F] - te[F].mean()) / sd).values
    return Ztr, Zte
def fit_ridge(tr, te, lam):
    Ztr, Zte = prep(tr, te)
    x = cen(tr, 'ld').values[:, None]; y = cen(tr, 'y').values
    X = np.hstack([x, Ztr]); P = np.diag([0] + [lam] * len(F))
    b = np.linalg.solve(X.T @ X + P + 1e-9 * np.eye(X.shape[1]), X.T @ y)
    pred_te = np.hstack([(te.ld - te.ld.mean()).values[:, None], Zte]) @ b
    return (te.y - te.y.mean()).values - pred_te, y - X @ b
lams = [1e0, 1e1, 1e2, 1e3, 1e4, 1e9]
def lomo(df, lam):
    return np.concatenate([fit_ridge(df[df.mouse != m], df[df.mouse == m], lam)[0] for m in df.mouse.unique()])
print('in-sample RMSE dose only %.3f' % np.sqrt(np.mean(fit_ridge(d, d, 1e9)[1] ** 2)),
      '| all features lam=1: %.3f' % np.sqrt(np.mean(fit_ridge(d, d, 1.)[1] ** 2)))
print('outer held-out RMSE by fixed penalty:', {l: round(float(np.sqrt(np.mean(lomo(d, l) ** 2))), 3) for l in lams})
# nested: choose lambda on training mice only
e = []
picks = []
for m in d.mouse.unique():
    tr, te = d[d.mouse != m], d[d.mouse == m]
    inner = {l: np.sqrt(np.mean(lomo(tr, l) ** 2)) for l in lams}; l = min(inner, key=inner.get); picks.append(l)
    e += list(fit_ridge(tr, te, l)[0])
print('nested-chosen penalty per held-out mouse:', picks, '-> outer RMSE %.3f (dose only 0.921)' % np.sqrt(np.mean(np.square(e))))
