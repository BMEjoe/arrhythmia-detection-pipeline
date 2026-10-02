import numpy as np, sys, json, time
from multiprocessing import Pool
sys.path.insert(0, '/home/user/arrhythmia-detection-pipeline')
from experiments.phase9_waveform.morph import lr1
def task(args):
    g, pcl = args
    apd, nead, y, _ = lr1.run(np.full(120, pcl), lr1.Y0, g, 1.0, 1.0, 2.5)
    a = apd[-60:]
    apd2, nead2, y2, ls = lr1.run(np.full(150, pcl), y, g, 1.0, 1.0, 2.5, 0.01, True)
    le = float(np.sum(ls) / (150 * pcl / 1000.0))
    a2 = apd2[~np.isnan(apd2)]
    return {"gsi": g, "pcl": pcl, "le_per_s": le, "n_nan": int(np.isnan(apd2).sum()),
            "n_distinct": int(len(np.unique(np.round(a2, 1)))), "apd_min": float(a2.min()) if len(a2) else None,
            "apd_max": float(a2.max()) if len(a2) else None, "ead_mean": float(nead2.mean())}
if __name__ == "__main__":
    gs = [float(x) for x in sys.argv[1].split(',')]
    pcls = np.arange(float(sys.argv[2]), float(sys.argv[3]) + 1e-9, float(sys.argv[4]))
    T = [(g, float(p)) for g in gs for p in pcls]
    t0 = time.time()
    with Pool(4) as pool:
        for r in pool.imap(task, T):
            print(json.dumps(r), flush=True)
    print("elapsed", time.time() - t0, file=sys.stderr)
