import numpy as np, sys, json
from multiprocessing import Pool
sys.path.insert(0, '/home/user/arrhythmia-detection-pipeline')
from experiments.phase9_waveform.morph import ktz
def task(P):
    le, apd = ktz.lyapunov_periodic(P, tm=400_000)
    a = apd[~np.isnan(apd)][-1000:]
    return {"P": int(P), "le_per_ts": float(le), "le_per_beat": float(le * P), "n_distinct_apd_last1000": int(len(np.unique(a))),
            "apd_mean": float(a.mean()), "apd_min": float(a.min()), "apd_max": float(a.max()), "n_nan": int(np.isnan(apd).sum())}
if __name__ == "__main__":
    with Pool(4) as pool:
        for r in pool.imap(task, range(80, 302, 2)):
            print(json.dumps(r), flush=True)
