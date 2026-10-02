"""
Logistic KTz map model of a paced cardiac cell (DEV morphology family).

Source (primary, open): Gall WM et al., "Cardiac reentry modeled by spatiotemporal chaos in a
coupled map lattice", arXiv:2202.12406, single-cell model Eqs. (1)-(3), Jacobian Eq. (7),
parameters Sect. 2.1 / Fig. 2-3: K = 0.6, Ie = 0, delta = lambda = 0.001, T = 0.154,
xr = -0.48; stimulus pulse I = 0.1 for 10 time steps (ts) every P ts (P from pulse start);
APD between upward and downward zero crossings of x (Sect. 3.1).
    x(t+1) = f((x(t) - K y(t) + z(t) + Ie + I(t)) / T),  f(u) = u / (1 + |u|)
    y(t+1) = x(t)
    z(t+1) = (1 - delta) z(t) - lambda (x(t) - xr)
Driven here by an arbitrary sequence of integer pacing intervals P_n (ts): pulse n starts at
sum_{m<n} P_m.
"""
from __future__ import annotations

import numpy as np
from numba import njit

K, IE, DELTA, LAM, T_PAR, XR = 0.6, 0.0, 0.001, 0.001, 0.154, -0.48
I_STIM, STIM_LEN = 0.1, 10


@njit(cache=True)
def run(periods, x0, y0, z0, T=T_PAR, xr=XR, want_le=False):
    """Iterate the map driven by integer pacing periods (ts).  Returns
    (apd[n_beats], dep_time[n_beats], le_sum, n_steps, x_final, y_final, z_final, le_per_beat);
    apd[n] = duration of the first x > 0 excursion starting after pulse n onset (NaN if none
    before the next pulse).  If want_le, the largest LE (per ts) is accumulated by QR of the
    Jacobian product (Eckmann-Ruelle / Benettin)."""
    nb = periods.shape[0]
    apd = np.full(nb, np.nan)
    dep = np.full(nb, np.nan)
    x, y, z = x0, y0, z0
    Q = np.eye(3)
    le = 0.0
    le_beat = np.zeros(nb)
    nsteps = 0
    t = 0
    for n in range(nb):
        P = periods[n]
        up = -1
        for s in range(P):
            I = I_STIM if s < STIM_LEN else 0.0
            u = (x - K * y + z + IE + I) / T
            xn = u / (1.0 + abs(u))
            if want_le:
                fp = 1.0 / ((1.0 + abs(u)) ** 2) / T
                J = np.array([[fp, -K * fp, fp], [1.0, 0.0, 0.0], [-LAM, 0.0, 1.0 - DELTA]])
                M = J @ np.ascontiguousarray(Q)
                Q, R = np.linalg.qr(M)
                for i in range(3):
                    if R[i, i] < 0:
                        R[i, :] = -R[i, :]
                        Q[:, i] = -Q[:, i]
                le += np.log(abs(R[0, 0]))
                le_beat[n] += np.log(abs(R[0, 0]))
            zn = (1.0 - DELTA) * z - LAM * (x - xr)
            yn = x
            # crossings (Sect. 3.1): depolarization x(t)<0<x(t+1), repolarization x(t)>0>x(t+1)
            if up < 0 and x < 0.0 and xn > 0.0 and np.isnan(apd[n]):
                up = t + 1
            elif up >= 0 and x > 0.0 and xn < 0.0 and np.isnan(apd[n]):
                apd[n] = (t + 1) - up
                dep[n] = up
            x, y, z = xn, yn, zn
            t += 1
            nsteps += 1
    return apd, dep, le, nsteps, x, y, z, le_beat


def lyapunov_periodic(P, n_beats=None, tm=1_000_000, tt=95_000):
    """Largest LE (per ts) at constant pacing P, IC x=y=z=0, transient tt, total tm (Fig. 3)."""
    nb_t = int(np.ceil(tt / P))
    nb = int(np.ceil(tm / P)) if n_beats is None else n_beats
    per = np.full(nb_t, P, dtype=np.int64)
    _, _, _, _, x, y, z, _ = run(per, 0.0, 0.0, 0.0)
    per = np.full(nb, P, dtype=np.int64)
    apd, dep, le, ns, *_ = run(per, x, y, z, want_le=True)
    return le / ns, apd
