"""Pan-Tompkins detector regression tests (deterministic, fixed seeds)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from dataclasses import replace
import numpy as np
import pytest
import final_pipeline as fp
from synth_ecg import make_ecg, score

REG = np.arange(1.0, 40.0, 0.85)


def run(x, fs, cfg=fp.CFG):
    return fp.detect_r_peaks(x, fs, cfg)


def report(name, m):
    print(f"[{name}] TP={m['TP']} FP={m['FP']} FN={m['FN']} Se={m['Se']:.3f} "
          f"PPV={m['PPV']:.3f} F1={m['F1']:.3f} err_mean={m['mean_err_ms']:.1f}ms "
          f"|err|max={m['max_abs_err_ms']:.1f}ms")


# --------------------------- basic / conditions -----------------------------
@pytest.mark.parametrize("fs", [200, 250, 360, 500, 1000])
def test_clean_multiple_sampling_rates(fs):
    x, tr = make_ecg(fs, REG)
    m = score(run(x, fs), tr, fs); report(f"clean fs={fs}", m)
    assert (m["FP"], m["FN"]) == (0, 0)
    assert m["max_abs_err_ms"] <= 20


def test_baseline_wander_noise_and_polarity():
    x, tr = make_ecg(360, REG, wander=1.0, noise=0.03, seed=11)
    m = score(run(x, 360), tr, 360); report("wander+noise", m)
    assert m["F1"] >= 0.99
    m = score(run(-x, 360), tr, 360); report("negative polarity", m)
    assert m["F1"] >= 0.99


@pytest.mark.parametrize("fs", [200, 360])
@pytest.mark.parametrize("offset", [0.0, 0.1, 1.0, 2.0, 5.0, 10.0, -1.0, -5.0, 1000.0])
def test_dc_offset_robustness(fs, offset):
    x, tr = make_ecg(fs, REG, offset=offset)
    m = score(run(x, fs), tr, fs); report(f"offset={offset} fs={fs}", m)
    assert (m["FP"], m["FN"]) == (0, 0)


def test_dc_offset_does_not_change_result():
    x, tr = make_ecg(200, REG)
    base = run(x, 200)
    for off in (1.0, 10.0, -7.0):
        np.testing.assert_array_equal(run(x + off, 200), base)


def test_first_beat_near_record_start_and_no_wraparound():
    """A beat in the last samples must not create/affect a detection at the start."""
    bt = np.arange(0.6, 20, 0.9)
    x, tr = make_ecg(200, bt)
    x[-30:] += 50.0                                # huge spike at the very end
    det = run(x, 200)
    assert det[0] > 0.4 * 200                      # nothing spurious near t=0
    m = score(det, tr, 200)
    assert m["FN"] == 0


# ------------------------------ rhythm cases -------------------------------
@pytest.mark.parametrize("name,rr", [("tachy", 0.40), ("brady", 1.6), ("regular", 0.8)])
def test_rate_range(name, rr):
    bt = np.arange(1.0, 1 + 40 * rr, rr)
    x, tr = make_ecg(360, bt, t_delay=0.25 if rr > 0.5 else 0.19)
    m = score(run(x, 360), tr, 360); report(name, m)
    assert m["F1"] >= 0.97


def test_irregular_rr():
    rng = np.random.default_rng(5)
    rr = rng.uniform(0.55, 1.3, 45)
    bt = 1.0 + np.cumsum(rr)
    x, tr = make_ecg(360, bt)
    m = score(run(x, 360), tr, 360); report("irregular", m)
    assert m["F1"] >= 0.97


def test_premature_beats_detected():
    rr = [0.85] * 12
    rr[5] = 0.42; rr[6] = 1.28          # premature beat, compensatory pause
    rr[9] = 0.45; rr[10] = 1.25
    bt = 1.0 + np.cumsum(rr * 2)
    bt = 1.0 + np.cumsum(np.array(rr * 2))
    x, tr = make_ecg(360, bt)
    m = score(run(x, 360), tr, 360); report("premature", m)
    assert m["FN"] == 0 and m["FP"] == 0


def test_pause():
    bt = np.concatenate([np.arange(1.0, 12, 0.85), np.arange(15.0, 30, 0.85)])
    x, tr = make_ecg(360, bt)
    m = score(run(x, 360), tr, 360); report("pause", m)
    assert m["FN"] == 0 and m["FP"] == 0


def test_broad_qrs():
    x, tr = make_ecg(360, REG, sigmas=[3.0] * len(REG), amps=[1.0] * len(REG))
    m = score(run(x, 360), tr, 360); report("broad QRS", m)
    assert m["F1"] >= 0.99


def test_noisy_beats():
    x, tr = make_ecg(360, REG, noise=0.08, seed=21)
    m = score(run(x, 360), tr, 360); report("noisy", m)
    assert m["F1"] >= 0.98


# ---------------------------- T-wave discrimination -------------------------
@pytest.mark.parametrize("t_amp", [0.3, 0.6, 0.9])
def test_tall_t_wave_is_not_a_qrs(t_amp):
    x, tr = make_ecg(360, REG, t_amp=t_amp, t_delay=0.26)
    m = score(run(x, 360), tr, 360); report(f"tall T {t_amp}", m)
    assert m["FP"] == 0 and m["FN"] == 0


def test_twave_slope_window_covers_the_qrs_rise():
    """The maximal QRS slope lies inside the integration window that ENDS at the
    integrated peak (the samples that produced it), for narrow and broad QRS."""
    for sig in (1.0, 3.0):
        x, _ = make_ecg(200, [2.0, 3.0], duration=4.0, sigmas=[sig, sig])
        hp, d, integ = fp._pan_tompkins_filter_chain_200hz(x)
        q = int(np.argmax(integ[:500]))
        lo = q - 60
        k = lo + int(np.argmax(np.abs(d[lo:q + 20])))
        assert q - 29 <= k <= q, (sig, q, k)


def test_slope_rule_rejects_t_wave_that_passes_amplitude_threshold():
    """T wave whose integrated peak exceeds THRESHOLD I1 but whose (band-passed)
    slope is < half the QRS slope (ratios ~0.38 energy / ~0.42 slope): ONLY the
    360-ms slope test can reject it."""
    kw = dict(t_amp=2.0, t_sigma=0.06, t_delay=0.28)
    x, tr = make_ecg(360, REG, **kw)
    with_rule = score(run(x, 360), tr, 360)
    without = score(run(x, 360, replace(fp.CFG, pt_twave_window_seconds=0.0)), tr, 360)
    report("slope-rule T wave, rule ON", with_rule); report("rule OFF (control)", without)
    assert without["FP"] > 10                     # T waves really are accepted without the rule
    assert with_rule["FP"] == 0 and with_rule["FN"] == 0


def test_premature_qrs_inside_360ms_not_rejected_as_twave():
    rr = [0.85] * 8 + [0.33] + [0.85] * 10       # coupling < 360 ms, same morphology
    bt = 1.0 + np.cumsum(rr)
    x, tr = make_ecg(360, bt, t_amp=0.2)
    m = score(run(x, 360), tr, 360); report("premature <360ms", m)
    assert m["FN"] == 0 and m["FP"] == 0


# --------------------------------- searchback --------------------------------
def weak_beat_ecg(fs=360, weak=0.45):
    bt = np.arange(1.0, 30.0, 0.85)
    amps = np.ones(len(bt)); k = 18; amps[k] = weak
    x, _ = make_ecg(fs, bt, amps=amps, t_amps=amps)
    return x, bt, k


def test_searchback_recovers_weak_beat():
    x, bt, k = weak_beat_ecg()
    tr = np.rint(bt * 360).astype(int)
    on = run(x, 360)
    off = run(x, 360, replace(fp.CFG, pt_searchback_seconds=1e9))
    m_on, m_off = score(on, tr, 360), score(off, tr, 360)
    report("weak beat, searchback ON", m_on); report("weak beat, searchback OFF", m_off)
    assert m_off["FN"] == 1                        # the weak beat really needs searchback
    assert m_on["FN"] == 0 and m_on["FP"] == 0


def test_searchback_uses_quarter_weighting_and_trace():
    x, bt, k = weak_beat_ecg(fs=200)
    hp, d, integ = fp._pan_tompkins_filter_chain_200hz(x)
    trace = []
    fp._pan_tompkins_decide(integ, hp, d, fp.CFG, 200.0, 30, trace=trace)
    sb = [t for t in trace if t["searchback"]]
    assert len(sb) == 1
    t = sb[0]
    assert t["spki_after"] == pytest.approx(0.25 * t["height"] + 0.75 * t["spki_before"])
    for t in trace:
        if not t["searchback"]:
            assert t["spki_after"] == pytest.approx(0.125 * t["height"] + 0.875 * t["spki_before"])


def test_missing_beat_is_not_invented():
    bt = np.delete(np.arange(1.0, 30.0, 0.85), 15)       # truly absent (no waveform)
    x, tr = make_ecg(360, bt)
    m = score(run(x, 360), tr, 360); report("missing beat", m)
    assert m["FP"] == 0 and m["FN"] == 0


def test_rejected_candidates_persist_until_next_qrs():
    """Searchback recovers a beat that was rejected *earlier* than the trigger."""
    x, bt, k = weak_beat_ecg(fs=200)
    hp, d, integ = fp._pan_tompkins_filter_chain_200hz(x)
    trace = []
    acc = fp._pan_tompkins_decide(integ, hp, d, fp.CFG, 200.0, 30, trace=trace)
    sb = [t for t in trace if t["searchback"]][0]
    later = [a for a in acc if a > sb["index"]]
    assert later and later[0] - sb["index"] > 100    # trigger came from a later beat


# ------------------------------ RR statistics -------------------------------
def tracker():
    return fp._RRTracker(8, 0.92, 1.16, 1.66)


def test_rr_limits_use_average2_and_paper_factors():
    t = tracker()
    for _ in range(8):
        t.add(200)
    lo, hi, miss = t.limits()
    assert (lo, hi, miss) == pytest.approx((184.0, 232.0, 332.0))


def test_rr_average1_vs_average2_diverge_on_outliers():
    t = tracker()
    for _ in range(8):
        t.add(200)
    t.add(100)        # premature: outside [184, 232] -> in AVERAGE1 only
    t.add(400)        # pause
    assert t.avg2 == pytest.approx(200.0)
    assert t.avg1 == pytest.approx(np.mean([200] * 6 + [100, 400]))
    assert t.avg1 != t.avg2
    assert t.irregular()


def test_regular_is_not_irregular_and_needs_eight_beats():
    t = tracker()
    for _ in range(7):
        t.add(200); assert not t.irregular()
    t.add(203)
    assert not t.irregular()
    t.add(300)
    assert t.irregular()


@pytest.mark.parametrize("rr", [100, 300])          # tachycardic / bradycardic @200 Hz
def test_rr_limits_scale_with_rhythm(rr):
    t = tracker()
    for _ in range(8):
        t.add(rr)
    lo, hi, miss = t.limits()
    assert (lo, hi, miss) == pytest.approx((0.92 * rr, 1.16 * rr, 1.66 * rr))
    assert not t.irregular()


def test_single_previous_rr_is_not_used():
    """One premature RR must not shift RR_AVERAGE2 (the baseline used only RR[-1])."""
    t = tracker()
    for _ in range(8):
        t.add(200)
    t.add(90)
    assert t.avg2 == pytest.approx(200.0)


def _weak_after(rr_before, weak=0.45, fs=360):
    rr = list(rr_before) + [0.85, 0.85]
    bt = 1.0 + np.cumsum(rr)
    amps = np.ones(len(bt)); amps[len(rr_before)] = weak
    x, _ = make_ecg(fs, bt, amps=amps, t_amps=amps)
    return x, bt, len(rr_before)


def test_irregular_rhythm_halves_first_threshold_only():
    """Same weak beat (integrated peak between THR_I1/2 and THR_I1), searchback
    disabled so only the first-pass threshold matters: missed in a regular
    rhythm, detected in an irregular rhythm (THRESHOLD I1 halved, Eq.22)."""
    cfg = replace(fp.CFG, pt_searchback_seconds=1e9)
    rng = np.random.default_rng(8)
    x_reg, bt_reg, k = _weak_after([0.85] * 20)
    x_irr, bt_irr, k2 = _weak_after(list(rng.uniform(0.55, 1.25, 20)))
    det_reg = score(run(x_reg, 360, cfg), np.rint(bt_reg * 360).astype(int), 360)
    det_irr = score(run(x_irr, 360, cfg), np.rint(bt_irr * 360).astype(int), 360)
    report("regular, weak beat", det_reg); report("irregular, weak beat", det_irr)
    assert det_reg["FN"] == 1
    assert det_irr["FN"] == 0 and det_irr["FP"] == 0


def test_detections_respect_refractory_period():
    for seed in (31, 32):
        x, _ = make_ecg(360, REG, noise=0.15, seed=seed)
        det = run(x, 360)
        assert np.all(np.diff(det) >= round(0.2 * 360))


def test_candidates_use_full_refractory_spacing():
    x, _ = make_ecg(200, REG, noise=0.05, seed=3)
    hp, d, integ = fp._pan_tompkins_filter_chain_200hz(x)
    seen = []
    orig = fp.find_peaks
    def spy(sig, *a, **k):
        seen.append(k.get("distance"))
        return orig(sig, *a, **k)
    fp.find_peaks = spy
    try:
        fp._pan_tompkins_decide(integ, hp, d, fp.CFG, 200.0, 30)
    finally:
        fp.find_peaks = orig
    assert seen and all(v == 40 for v in seen)      # 0.2 s * 200 Hz, not 20


# ------------------------------- startup -----------------------------------
def test_short_record_returns_array():
    out = run(np.zeros(5), 200)
    assert isinstance(out, np.ndarray) and out.size == 0
