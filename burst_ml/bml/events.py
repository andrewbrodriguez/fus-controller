"""Event statistics of the within-record envelope (plan: burst_ml/plan_events.md). Mice 1-6, label-free."""
import warnings, numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt, hilbert
from scipy.stats import kurtosis
import fus.archive as A
from fus.extract import extract
from . import CACHE
warnings.simplefilter('ignore')
def stats(e):
    med = np.median(e) + 1e-12
    return np.log(kurtosis(e, fisher=False) + 1e-12), np.mean(e > 2 * med), np.log(e.max() / med)
rows = []
for m in A.recordings():
    if m.mouse > 6: continue
    with A.local_path(m) as p: r = extract(p).recording
    fs = r.sample_rate_hz
    sos2 = butter(4, [2 * r.carrier_hz * .9, 2 * r.carrier_hz * 1.1], btype='band', fs=fs, output='sos')
    sosu = butter(4, [1.5 * r.carrier_hz * .9, 1.5 * r.carrier_hz * 1.1], btype='band', fs=fs, output='sos')
    for t, b in enumerate(np.atleast_1d(r.raw.niscope.wf)):
        if t < r.n_baseline: continue
        w = np.asarray(b.data0, float); w = w - w.mean()
        e2 = np.abs(hilbert(sosfiltfilt(sos2, w))) ** 2; eu = np.abs(hilbert(sosfiltfilt(sosu, w))) ** 2
        k2, f2, p2 = stats(e2); ku, fu, _ = stats(eu)
        rows.append(dict(mouse=m.mouse, recording=m.recording, t=t - r.n_baseline, k2=k2, f2=f2, p2=p2, ku=ku, fu=fu))
    print(m.mouse, m.recording, flush=True)
pd.DataFrame(rows).to_csv(CACHE / 'events.csv', index=False)
