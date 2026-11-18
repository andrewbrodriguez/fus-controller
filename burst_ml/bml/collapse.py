"""Within-burst 2nd-harmonic collapse features (plan: burst_ml/plan_collapse.md). Mice 1-6, label-free.
c2 = log(E2nd second half / first half) of the 30 ms record; r2 = log(E_t / E_{t-1}) of whole-record 2f0 energy.
Writes burst_ml/cache/collapse.csv (per sonicating burst). Reads the zip via fus.archive only."""
import warnings, numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt, hilbert
import fus.archive as A
from fus.extract import extract
from . import CACHE
warnings.simplefilter('ignore')
rows = []
for m in A.recordings():
    if m.mouse > 6: continue
    with A.local_path(m) as p: r = extract(p).recording
    fs, f2 = r.sample_rate_hz, 2 * r.carrier_hz
    sos = butter(4, [f2 * .9, f2 * 1.1], btype='band', fs=fs, output='sos')
    wf = [np.asarray(b.data0, float) for b in np.atleast_1d(r.raw.niscope.wf)]
    E = []; C = []
    for w in wf:
        e = np.abs(hilbert(sosfiltfilt(sos, w - w.mean()))) ** 2; h = len(e) // 2
        E.append(e.sum()); C.append(np.log((e[h:].sum() + 1e-12) / (e[:h].sum() + 1e-12)))
    E = np.array(E); nb = r.n_baseline
    r2 = np.r_[np.nan, np.log((E[1:] + 1e-12) / (E[:-1] + 1e-12))]
    for t in range(nb, len(wf)):
        rows.append(dict(mouse=m.mouse, recording=m.recording if hasattr(m, 'recording') else str(m), t=t - nb, c2=C[t], r2=r2[t], E=E[t] / E[:nb].mean()))
    print(m, len(wf), flush=True)
pd.DataFrame(rows).to_csv(CACHE / 'collapse.csv', index=False)
