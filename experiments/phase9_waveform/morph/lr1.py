"""
Modified Luo-Rudy I (LR1) ventricular action-potential model with early afterdepolarizations
(EADs) (TEST morphology family).

Sources:
- Model modification and the published behaviour reproduced in verification (primary, open):
  Tran DX, Sato D, Yochelis A, Weiss JN, Garfinkel A, Qu Z, "Bifurcation and chaos in a model
  of cardiac early afterdepolarizations", Phys Rev Lett 102:258103 (2009), NIH author
  manuscript PMC2726623: Eq. (1); E_si = 80 mV and E_K = -77 mV fixed; tau_d -> alpha tau_d,
  tau_f -> beta tau_f, tau_x -> gamma tau_x; stimulus 2 ms, 30 uA/cm^2; RK4 with dt = 0.01 ms;
  APD = duration of V > -72 mV; Fig. 4: gamma = 2.5, chaos at intermediate PCL,
  LE = 0.38 s^-1 at PCL = 0.907 s.
- LR1 equations and parameters (the original Luo & Rudy, Circ Res 68:1501 (1991), is not
  accessible): the CellML encoding luo_rudy_1991.cellml from the Physiome Model Repository
  (models.physiomeproject.org/workspace/luo_rudy_1991), transcribed term by term below.
With E_si fixed, [Ca]i no longer affects any current and is not integrated.
Scaling a gate's time constant by c with fixed steady state is alpha_y, beta_y -> alpha_y/c,
beta_y/c.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

R_GAS, T_K, F_C = 8314.0, 310.0, 96484.6
NAO, NAI, KI, KO = 140.0, 18.0, 145.0, 5.4
PR_NAK = 0.01833
G_NA, G_KP, G_B, E_B = 23.0, 0.0183, 0.03921, -59.87
RTF = R_GAS * T_K / F_C
E_NA = RTF * math.log(NAO / NAI)
G_K = 0.282 * math.sqrt(KO / 5.4)
G_K1 = 0.6047 * math.sqrt(KO / 5.4)
E_K1 = RTF * math.log(KO / KI)
E_SI_TRAN, E_K_TRAN = 80.0, -77.0          # Tran et al. 2009
STIM_AMP, STIM_DUR = 30.0, 2.0             # uA/cm^2, ms (Tran et al. 2009)
APD_THRESH = -72.0                         # mV (Tran et al. 2009, Fig. 4)
# CellML initial state (V, m, h, j, d, f, X)
Y0 = np.array([-84.3801107371, 0.00171338077730188, 0.982660523699656, 0.989108212766685,
               0.00302126301779861, 0.999967936476325, 0.0417603108167287])


@njit(cache=True)
def rhs(y, dy, istim, gsi, ad, bf, gx):
    V, m, h, j, d, f, X = y[0], y[1], y[2], y[3], y[4], y[5], y[6]
    # fast sodium
    am = 0.32 * (V + 47.13) / (1.0 - math.exp(-0.1 * (V + 47.13))) if abs(V + 47.13) > 1e-9 else 3.2
    bm = 0.08 * math.exp(-V / 11.0)
    if V < -40.0:
        ah = 0.135 * math.exp((80.0 + V) / -6.8)
        bh = 3.56 * math.exp(0.079 * V) + 310000.0 * math.exp(0.35 * V)
        aj = (-127140.0 * math.exp(0.2444 * V) - 0.00003474 * math.exp(-0.04391 * V)) * (V + 37.78) / \
            (1.0 + math.exp(0.311 * (V + 79.23)))
        bj = 0.1212 * math.exp(-0.01052 * V) / (1.0 + math.exp(-0.1378 * (V + 40.14)))
    else:
        ah = 0.0
        bh = 1.0 / (0.13 * (1.0 + math.exp((V + 10.66) / -11.1)))
        aj = 0.0
        bj = 0.3 * math.exp(-0.0000002535 * V) / (1.0 + math.exp(-0.1 * (V + 32.0)))
    i_na = G_NA * m * m * m * h * j * (V - E_NA)
    # slow inward (Tran: E_si fixed)
    a_d = 0.095 * math.exp(-0.01 * (V - 5.0)) / (1.0 + math.exp(-0.072 * (V - 5.0))) / ad
    b_d = 0.07 * math.exp(-0.017 * (V + 44.0)) / (1.0 + math.exp(0.05 * (V + 44.0))) / ad
    a_f = 0.012 * math.exp(-0.008 * (V + 28.0)) / (1.0 + math.exp(0.15 * (V + 28.0))) / bf
    b_f = 0.0065 * math.exp(-0.02 * (V + 30.0)) / (1.0 + math.exp(-0.2 * (V + 30.0))) / bf
    i_si = gsi * d * f * (V - E_SI_TRAN)
    # time-dependent K (Tran: E_K fixed)
    a_x = 0.0005 * math.exp(0.083 * (V + 50.0)) / (1.0 + math.exp(0.057 * (V + 50.0))) / gx
    b_x = 0.0013 * math.exp(-0.06 * (V + 20.0)) / (1.0 + math.exp(-0.04 * (V + 20.0))) / gx
    if V > -100.0:
        if abs(V + 77.0) > 1e-9:
            xi = 2.837 * (math.exp(0.04 * (V + 77.0)) - 1.0) / ((V + 77.0) * math.exp(0.04 * (V + 35.0)))
        else:
            xi = 2.837 * 0.04 / math.exp(0.04 * (V + 35.0))
    else:
        xi = 1.0
    i_k = G_K * X * xi * (V - E_K_TRAN)
    # time-independent K
    ak1 = 1.02 / (1.0 + math.exp(0.2385 * (V - E_K1 - 59.215)))
    bk1 = (0.49124 * math.exp(0.08032 * (V + 5.476 - E_K1)) + math.exp(0.06175 * (V - (E_K1 + 594.31)))) / \
        (1.0 + math.exp(-0.5143 * (V - E_K1 + 4.753)))
    i_k1 = G_K1 * ak1 / (ak1 + bk1) * (V - E_K1)
    kp = 1.0 / (1.0 + math.exp((7.488 - V) / 5.98))
    i_kp = G_KP * kp * (V - E_K1)
    i_b = G_B * (V - E_B)
    dy[0] = -(i_na + i_si + i_k + i_k1 + i_kp + i_b) + istim
    dy[1] = am * (1.0 - m) - bm * m
    dy[2] = ah * (1.0 - h) - bh * h
    dy[3] = aj * (1.0 - j) - bj * j
    dy[4] = a_d * (1.0 - d) - b_d * d
    dy[5] = a_f * (1.0 - f) - b_f * f
    dy[6] = a_x * (1.0 - X) - b_x * X


@njit(cache=True)
def _step(y, h, istim, gsi, ad, bf, gx, k1, k2, k3, k4, yt):
    rhs(y, k1, istim, gsi, ad, bf, gx)
    for i in range(7):
        yt[i] = y[i] + 0.5 * h * k1[i]
    rhs(yt, k2, istim, gsi, ad, bf, gx)
    for i in range(7):
        yt[i] = y[i] + 0.5 * h * k2[i]
    rhs(yt, k3, istim, gsi, ad, bf, gx)
    for i in range(7):
        yt[i] = y[i] + h * k3[i]
    rhs(yt, k4, istim, gsi, ad, bf, gx)
    for i in range(7):
        y[i] = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])


@njit(cache=True)
def run(cl_ms, y0, gsi, ad, bf, gx, dt=0.01, le=False, d0=1e-8):
    """Pace with cycle lengths cl_ms[n] (ms; stimulus at the start of each cycle).  Returns
    (apd[n] ms, n_ead[n], y_end, le_sum_per_beat[n]).  apd[n] = time from the first upward
    crossing of APD_THRESH after stimulus n to the next downward crossing (NaN if the cell does
    not repolarize within the cycle; the AP then continues into the next cycle and that beat's
    APD is NaN).  n_ead[n] = number of local voltage maxima above -40 mV after the upstroke
    (EAD count).  If le, a second trajectory at relative distance d0 (in the scaled state
    space V/100, gates) is integrated and renormalized after every beat (Benettin); the log
    growth per beat is returned."""
    nb = cl_ms.shape[0]
    apd = np.full(nb, np.nan)
    nead = np.zeros(nb, dtype=np.int64)
    lsum = np.zeros(nb)
    y = y0.copy()
    k1 = np.empty(7); k2 = np.empty(7); k3 = np.empty(7); k4 = np.empty(7); yt = np.empty(7)
    z = y0.copy()
    scale = np.array([100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    if le:
        v = np.ones(7) / math.sqrt(7.0)
        for i in range(7):
            z[i] = y[i] + d0 * v[i] * scale[i]
    above = y[0] > APD_THRESH
    t_up = -1.0
    for n in range(nb):
        nsteps = int(round(cl_ms[n] / dt))
        vprev2 = y[0]
        vprev = y[0]
        for s in range(nsteps):
            t = s * dt
            istim = STIM_AMP if t < STIM_DUR else 0.0
            _step(y, dt, istim, gsi, ad, bf, gx, k1, k2, k3, k4, yt)
            if le:
                _step(z, dt, istim, gsi, ad, bf, gx, k1, k2, k3, k4, yt)
            V = y[0]
            if (not above) and V > APD_THRESH:
                above = True
                t_up = t
            elif above and V <= APD_THRESH:
                above = False
                if t_up >= 0.0 and np.isnan(apd[n]):
                    apd[n] = t - t_up
                t_up = -1.0
            if t_up >= 0.0 and vprev > vprev2 and vprev > V and vprev > -40.0 and t - t_up > 20.0:
                nead[n] += 1
            vprev2 = vprev
            vprev = V
        if t_up >= 0.0:
            t_up = t_up - cl_ms[n]       # AP continues into the next cycle
        if le:
            dist = 0.0
            for i in range(7):
                dd = (z[i] - y[i]) / scale[i]
                dist += dd * dd
            dist = math.sqrt(dist)
            lsum[n] = math.log(dist / d0)
            for i in range(7):
                z[i] = y[i] + (z[i] - y[i]) * (d0 / dist)
    return apd, nead, y, lsum
