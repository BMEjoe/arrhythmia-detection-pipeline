"""
Seidel-Herzel (SH) model of the baroreceptor-cardiac reflex.

SOURCE: Seidel H, Herzel H, "Bifurcations in a nonlinear model of the baroreceptor-
cardiac reflex", Physica D 115:145-160, 1998 (not accessible here).  Equations (1)-(14)
and the parameter values as restated in Dudkowska A, Makowiec D, "Influence of
stochastic delays in Seidel-Herzel model of human cardiorespiratory system",
arXiv:q-bio/0603016 (2006), Section II, cross-checked against the restatement in
Hung K et al., PMC13305918 (2026), Eqs. (1)-(9) (identical forms; systole "of fixed
duration tau_sys", beat when the sinus phase reaches 1).

  (1)  nu_b = k1 (p - p0) + k2 dp/dt
  (2)  nu_s = max(0, nu_s0 - ksb nu_b + ksr |sin(pi fr t + dphi_s)|)
  (3)  nu_p = max(0, nu_p0 + kpb nu_b + kpr |sin(pi fr t + dphi_p)|)
  (4)  dc_cNa/dt = -c_cNa/tau_cNa + kc nu_s(t - theta_cNa)
  (5)  dphi/dt = f_s f_p / T0
  (6)  f_s = 1 + kphi_c [c_cNa + (chat_c - c_cNa) c_cNa^n_c / (chat_c^n_c + c_cNa^n_c)]
  (7)  f_p = 1 - kphi_p [nu_pd + (nuhat_p - nu_pd) nu_pd^n_p / (nuhat_p^n_p + nu_pd^n_p)] F(phi),
       nu_pd = nu_p(t - theta_p)
  (8)  F(phi) = phi^1.3 (phi - 0.45) (1 - phi)^3 / ((1 - 0.8)^3 + (1 - phi)^3)
  (9)  S'_i = S0 + kSc c_cNa + kSt T_{i-1}
  (10) S_i = S'_i + (Shat - S'_i) S'_i^nS / (S'_i^nS + Shat^nS)
  (11) tau_v = tau_v0 - dtau_v [c_vNa + (chat_v - c_vNa) c_vNa^n_v / (chat_v^n_v + c_vNa^n_v)]
  (12) dc_vNa/dt = -c_vNa/tau_vNa + kv nu_s(t - theta_vNa)
  (13) systole (t_i <= t < t_i + tau_sys): p = d_{i-1} + S_i x exp(1 - x), x = (t - t_i)/tau_sys
  (14) diastole: dp/dt = -p / tau_v(t)
  Beat: phi reaches 1 -> phi reset to 0, t_i = beat time, T_{i-1} = previous heart period,
  d_{i-1} = pressure at the beat, new S_i from (9)-(10).

Parameters: Dudkowska & Makowiec Section II; chat_v = 1.0 (their proposed modification
of 10.0, which lets tau_v become negative; they report it does not change the regular
regime); theta_p = 0.5 s (their Section III).  "Constant respiration" (their Figs 3b-5,
used for the published bifurcation diagrams): |sin| replaced by its mean 2/pi.
Dynamical noise (their Section IV, Eq. 15): stochastic delays theta + xi, xi ~ U(-xihat,
xihat), drawn once per cardiac cycle.

Beat events: the model's own beats (phase reset); RR = heart periods T_i in seconds.
Integration: RK4, h = 0.001 s (as Dudkowska & Makowiec), delays on the grid, half-step
delayed values interpolated.
"""
from __future__ import annotations

import numba
import numpy as np

PARAMS = dict(k1=0.02, k2=0.00125, p0=50.0,
              nus0=0.8, ksb=0.7, ksr=0.1, dphis=0.0,
              nup0=0.0, kpb=0.3, kpr=0.1, dphip=0.0, fr=0.2,
              tau_c=2.0, kc=1.2, tau_vna=2.0, kv=1.2,
              T0=1.1, kphic=1.6, chatc=2.0, nc=2.0,
              kphip=5.8, nuhatp=2.5, np_=2.0,
              S0=25.0, kSc=40.0, kSt=10.0, Shat=70.0, nS=2.5,
              tauv0=2.2, dtauv=1.2, chatv=1.0, nv=1.5,
              tau_sys=0.125, theta_p=0.5,
              sys_dpdt=1.0)   # 1: nu_b uses dp/dt of Eq. (13) during systole (default); 0: variant check only

# state vector layout
T_, P_, CC_, CV_, PHI_, TI_, D_, S_, TPREV_, NB_, HEAD_, DPDT_ = range(12)
NSTATE = 12


@numba.njit(cache=True)
def _sat(x, xhat, n):
    if x <= 0.0:
        return 0.0
    xn = x ** n
    return x + (xhat - x) * xn / (xhat ** n + xn)


@numba.njit(cache=True)
def _F(phi):
    if phi <= 0.0:
        return 0.0
    om = 1.0 - phi
    if om < 0.0:
        om = 0.0
    return phi ** 1.3 * (phi - 0.45) * om ** 3 / (0.2 ** 3 + om ** 3)


@numba.njit(cache=True)
def _resp(t, fr, dph, const):
    if const:
        return 2.0 / np.pi
    return abs(np.sin(np.pi * fr * t + dph))


@numba.njit(cache=True)
def _derivs(t, cc, cv, phi, p, nus_dc, nus_dv, nup_d, pr):
    # pr: packed params array (see pack())
    dcc = -cc / pr[12] + pr[13] * nus_dc
    dcv = -cv / pr[14] + pr[15] * nus_dv
    fs = 1.0 + pr[17] * _sat(cc, pr[18], pr[19])
    fp = 1.0 - pr[20] * _sat(nup_d, pr[21], pr[22]) * _F(phi)
    dphi = fs * fp / pr[16]
    tauv = pr[28] - pr[29] * _sat(cv, pr[30], pr[31])
    dp = -p / tauv
    return dcc, dcv, dphi, dp


@numba.njit(cache=True)
def _nus_nup(t, p, dpdt, pr, const):
    nub = pr[0] * (p - pr[2]) + pr[1] * dpdt
    nus = pr[3] - pr[4] * nub + pr[5] * _resp(t, pr[11], pr[6], const)
    nup = pr[7] + pr[8] * nub + pr[9] * _resp(t, pr[11], pr[10], const)
    return max(0.0, nus), max(0.0, nup)


@numba.njit(cache=True)
def _buf(buf, head, back, L):
    return buf[(head - back) % L]


@numba.njit(cache=True)
def advance(s, nus_buf, nup_buf, pr, lag_c, lag_v, lag_p, offs_c, offs_v, dt, const, max_steps):
    """Integrate until the next beat (or max_steps).  Returns 1 if a beat occurred.
    nus_buf / nup_buf: ring buffers of nu_s, nu_p on this copy's own time grid, head = last
    written index.  offs_c/offs_v: per-beat delay offsets (steps) for stochastic delays."""
    L = nus_buf.shape[0]
    tau_sys = pr[32]
    for _ in range(max_steps):
        t = s[T_]; p = s[P_]; cc = s[CC_]; cv = s[CV_]; phi = s[PHI_]
        head = int(s[HEAD_])
        nb = int(s[NB_])
        lc = lag_c + offs_c[nb % offs_c.shape[0]]
        lv = lag_v + offs_v[nb % offs_v.shape[0]]
        # delayed inputs at t, t+dt/2, t+dt (grid values back from head; head = value at t)
        a0c = _buf(nus_buf, head, lc, L); a1c = _buf(nus_buf, head, lc - 1, L)
        a0v = _buf(nus_buf, head, lv, L); a1v = _buf(nus_buf, head, lv - 1, L)
        a0p = _buf(nup_buf, head, lag_p, L); a1p = _buf(nup_buf, head, lag_p - 1, L)
        ahc = 0.5 * (a0c + a1c); ahv = 0.5 * (a0v + a1v); ahp = 0.5 * (a0p + a1p)
        insys = (t - s[TI_]) < tau_sys
        k1 = _derivs(t, cc, cv, phi, p, a0c, a0v, a0p, pr)
        k2 = _derivs(t + 0.5 * dt, cc + 0.5 * dt * k1[0], cv + 0.5 * dt * k1[1], phi + 0.5 * dt * k1[2],
                     p + 0.5 * dt * k1[3], ahc, ahv, ahp, pr)
        k3 = _derivs(t + 0.5 * dt, cc + 0.5 * dt * k2[0], cv + 0.5 * dt * k2[1], phi + 0.5 * dt * k2[2],
                     p + 0.5 * dt * k2[3], ahc, ahv, ahp, pr)
        k4 = _derivs(t + dt, cc + dt * k3[0], cv + dt * k3[1], phi + dt * k3[2], p + dt * k3[3], a1c, a1v, a1p, pr)
        ccn = cc + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        cvn = cv + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        phin = phi + dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        tn = t + dt
        # pressure
        if insys:
            xe = (tn - s[TI_]) / tau_sys
            if xe <= 1.0:
                pn = s[D_] + s[S_] * xe * np.exp(1.0 - xe)
                dpn = s[S_] / tau_sys * np.exp(1.0 - xe) * (1.0 - xe)
            else:                                   # systole ends inside the step
                pend = s[D_] + s[S_]
                rem = tn - (s[TI_] + tau_sys)
                tauv = pr[28] - pr[29] * _sat(cvn, pr[30], pr[31])
                pn = pend * np.exp(-rem / tauv)
                dpn = -pn / tauv
        else:
            pn = p + dt / 6.0 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
            tauv = pr[28] - pr[29] * _sat(cvn, pr[30], pr[31])
            dpn = -pn / tauv
        beat = False
        if phin >= 1.0:
            frac = (1.0 - phi) / (phin - phi) if phin > phi else 1.0
            tb = t + frac * dt
            pb = p + frac * (pn - p)
            ccb = cc + frac * (ccn - cc)
            s[TPREV_] = tb - s[TI_]
            Sp = pr[24] + pr[25] * ccb + pr[26] * s[TPREV_]
            if Sp < 0.0:
                Sp = 0.0
            Spn = Sp ** pr[27]
            s[S_] = Sp + (pr[23] - Sp) * Spn / (Spn + pr[23] ** pr[27])
            s[TI_] = tb
            s[D_] = pb
            xe = (tn - tb) / tau_sys
            pn = pb + s[S_] * xe * np.exp(1.0 - xe)
            dpn = s[S_] / tau_sys * np.exp(1.0 - xe) * (1.0 - xe)
            phin = phin - 1.0
            s[NB_] = nb + 1
            beat = True
        s[T_] = tn; s[P_] = pn; s[CC_] = ccn; s[CV_] = cvn; s[PHI_] = phin; s[DPDT_] = dpn
        if pr[33] == 0.0 and (tn - s[TI_]) < tau_sys:   # variant: no systolic dp/dt in nu_b
            nus, nup = _nus_nup(tn, pn, 0.0, pr, const)
        else:
            nus, nup = _nus_nup(tn, pn, dpn, pr, const)
        head = (head + 1) % L
        nus_buf[head] = nus
        nup_buf[head] = nup
        s[HEAD_] = head
        if beat:
            return 1
    return 0


def pack(P):
    keys = ["k1", "k2", "p0", "nus0", "ksb", "ksr", "dphis", "nup0", "kpb", "kpr", "dphip", "fr",
            "tau_c", "kc", "tau_vna", "kv", "T0", "kphic", "chatc", "nc", "kphip", "nuhatp", "np_",
            "Shat", "S0", "kSc", "kSt", "nS", "tauv0", "dtauv", "chatv", "nv", "tau_sys", "sys_dpdt"]
    return np.array([float(P[k]) for k in keys])


class SH:
    """One model copy (state + ring buffers)."""

    def __init__(self, theta_c, theta_v, P=None, dt=0.001, const_resp=True, xihat=0.0, rng=None,
                 n_offsets=200_000):
        self.P = dict(PARAMS, **(P or {}))
        self.pr = pack(self.P)
        self.dt = dt
        self.const = bool(const_resp)
        self.lag_c = int(round(theta_c / dt)); self.lag_v = int(round(theta_v / dt))
        self.lag_p = int(round(self.P["theta_p"] / dt))
        xi_steps = int(round(xihat / dt))
        self.L = max(self.lag_c, self.lag_v, self.lag_p) + xi_steps + 4
        if xi_steps > 0:
            self.offs_c = rng.integers(-xi_steps, xi_steps + 1, n_offsets).astype(np.int64)
            self.offs_v = rng.integers(-xi_steps, xi_steps + 1, n_offsets).astype(np.int64)
        else:
            self.offs_c = np.zeros(1, np.int64); self.offs_v = np.zeros(1, np.int64)
        self.s = np.zeros(NSTATE)
        self.s[P_] = 80.0; self.s[CC_] = 0.5; self.s[CV_] = 0.5; self.s[PHI_] = 0.0
        self.s[TI_] = -1.0; self.s[D_] = 80.0; self.s[S_] = 60.0; self.s[TPREV_] = 0.8
        self.nus = np.full(self.L, 0.3); self.nup = np.full(self.L, 0.3)
        self.s[HEAD_] = 0

    def copy(self):
        o = object.__new__(SH)
        o.__dict__ = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.__dict__.items()}
        return o

    def next_beat(self, max_steps=20_000):
        ok = advance(self.s, self.nus, self.nup, self.pr, self.lag_c, self.lag_v, self.lag_p,
                     self.offs_c, self.offs_v, self.dt, self.const, max_steps)
        if not ok:
            raise RuntimeError("no beat within max_steps")
        return self.s[TI_]

    def beats(self, n):
        out = np.empty(n)
        for k in range(n):
            out[k] = self.next_beat()
        return out

    def vector(self):
        """State at a beat for the LE distance: continuous variables + delay histories
        (ordered back from the head), weighted to comparable scale."""
        s = self.s
        L = self.L
        h = int(s[HEAD_])
        idx = (h - np.arange(L - 2)) % L
        wv = 1.0 / np.sqrt(L)
        return np.concatenate([[s[P_] / 100.0, s[CC_], s[CV_], s[PHI_], s[D_] / 100.0, s[S_] / 100.0,
                                s[TPREV_], s[T_] - s[TI_]], self.nus[idx] * wv, self.nup[idx] * wv]), idx


def simulate_rr(theta_c, theta_v, n_beats, transient_beats=600, const_resp=True, xihat=0.0, rng=None, P=None):
    m = SH(theta_c, theta_v, P=P, const_resp=const_resp, xihat=xihat, rng=rng)
    m.beats(transient_beats)
    tb = m.beats(n_beats + 1)
    return np.diff(tb)


def lyapunov(theta_c, theta_v, n_beats=4000, transient_beats=800, eps=1e-7, const_resp=True, P=None,
             seed=0, n_blocks=20):
    """Largest LE of the beat-to-beat (Poincare section at phi = 0) dynamics by the
    two-copy Benettin method on the full state (continuous variables + delay histories
    + time offset for the respiratory forcing), renormalized at every beat.  Returns
    (lambda per beat, SE, lambda per second, convergence curve)."""
    rng = np.random.default_rng(seed)
    ref = SH(theta_c, theta_v, P=P, const_resp=const_resp)
    ref.beats(transient_beats)
    pert = ref.copy()
    v, idx = ref.vector()
    # random initial perturbation of continuous state and histories
    d0 = rng.normal(size=6)
    pert.s[[P_, CC_, CV_, D_, S_, TPREV_]] += eps * d0 * np.array([100, 1, 1, 100, 100, 1])
    pert.nus += eps * rng.normal(size=pert.L)
    pert.nup += eps * rng.normal(size=pert.L)
    logs = np.empty(n_beats)
    tsum = 0.0
    for k in range(n_beats):
        t0 = ref.s[TI_]
        ref.next_beat(); pert.next_beat()
        tsum += ref.s[TI_] - t0
        vr, ir = ref.vector(); vp, ip = pert.vector()
        dvec = vp - vr
        dvec = np.append(dvec, (pert.s[TI_] - ref.s[TI_]) if not ref.const else 0.0)
        dist = float(np.sqrt(np.sum(dvec ** 2)))
        logs[k] = np.log(dist / eps)
        sc = eps / dist
        # rescale the perturbed copy toward the reference (aligned by "steps back from head")
        for j in (P_, CC_, CV_, PHI_, D_, S_, TPREV_, DPDT_):
            pert.s[j] = ref.s[j] + (pert.s[j] - ref.s[j]) * sc
        dt_off = (pert.s[TI_] - ref.s[TI_]) * sc
        pert.s[T_] = ref.s[T_] + ((pert.s[T_] - pert.s[TI_]) - (ref.s[T_] - ref.s[TI_])) * sc + dt_off
        pert.s[TI_] = ref.s[TI_] + dt_off
        pert.s[NB_] = ref.s[NB_]
        L = ref.L
        hr, hp = int(ref.s[HEAD_]), int(pert.s[HEAD_])
        back = np.arange(L)
        for buf_r, buf_p in ((ref.nus, pert.nus), (ref.nup, pert.nup)):
            br = buf_r[(hr - back) % L]
            bp = buf_p[(hp - back) % L]
            buf_p[(hp - back) % L] = br + (bp - br) * sc
    lam = float(logs.mean())
    blocks = logs[: n_beats - n_beats % n_blocks].reshape(n_blocks, -1).mean(1)
    conv = {int(k): float(logs[:k].mean()) for k in (n_beats // 8, n_beats // 4, n_beats // 2, n_beats)}
    return lam, float(blocks.std(ddof=1) / np.sqrt(n_blocks)), lam / (tsum / n_beats), conv
