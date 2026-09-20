"""Pan-Tompkins filter regression tests (1985 paper Eqs. 1-9)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from scipy.signal import freqz
import final_pipeline as fp

FS = 200.0


# ---- reference recursive realisations, with explicit zero-state guards ----
def x_at(x, i):
    return x[i] if i >= 0 else 0.0


def lp_recursive(x):
    """Paper Eq.3 with unscaled recursion; the 1/36 gain is applied to the OUTPUT."""
    y = np.zeros(len(x))
    for n in range(len(x)):
        y[n] = ((2 * y[n - 1] if n >= 1 else 0) - (y[n - 2] if n >= 2 else 0)
                + x_at(x, n) - 2 * x_at(x, n - 6) + x_at(x, n - 12))
    return y / 36.0


def lp_recursive_normalised_state(x):
    """Same filter, gain applied to the recursion inputs only (feedback un-scaled)."""
    y = np.zeros(len(x))
    for n in range(len(x)):
        y[n] = ((2 * y[n - 1] if n >= 1 else 0) - (y[n - 2] if n >= 2 else 0)
                + (x_at(x, n) - 2 * x_at(x, n - 6) + x_at(x, n - 12)) / 36.0)
    return y


def lp_old_buggy(x):
    """The baseline implementation: whole recurrence (incl. feedback) / 36."""
    y = np.zeros(len(x))
    for n in range(len(x)):
        v = x[n]
        if n >= 1: v += 2 * y[n - 1]
        if n >= 2: v -= y[n - 2]
        if n >= 6: v -= 2 * x[n - 6]
        if n >= 12: v += x[n - 12]
        y[n] = v / 36.0
    return y


def hp_recursive(x):
    y = np.zeros(len(x))
    for n in range(len(x)):
        y[n] = ((y[n - 1] if n >= 1 else 0) - x_at(x, n) / 32.0 + x_at(x, n - 16)
                - x_at(x, n - 17) + x_at(x, n - 32) / 32.0)
    return y


def hp_printed_eq6(x):
    """Paper's printed Eq.6 read literally: y = 32x[n-16] - [y[n-1] + x[n] - x[n-32]] (/32)."""
    y = np.zeros(len(x))
    for n in range(len(x)):
        y[n] = 32 * x_at(x, n - 16) - ((y[n - 1] if n >= 1 else 0) + x_at(x, n) - x_at(x, n - 32))
    return y / 32.0


def impulse(n=400, at=0):
    x = np.zeros(n); x[at] = 1.0
    return x


# ------------------------------- low-pass ---------------------------------
def test_lp_fir_equals_recursion_impulse():
    x = impulse()
    fir = np.convolve(x, fp._PT_LP_TAPS)[:len(x)]
    np.testing.assert_allclose(fir, lp_recursive(x), atol=1e-12)
    np.testing.assert_allclose(fir, lp_recursive_normalised_state(x), atol=1e-12)


def test_lp_impulse_is_unit_gain_triangle():
    h = fp._PT_LP_TAPS
    np.testing.assert_allclose(h * 36, [1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1])
    assert np.argmax(h) == 5            # linear-phase group delay (paper quotes 6)
    np.testing.assert_allclose(h, h[::-1])
    assert h.sum() == pytest.approx(1.0)


def test_lp_dc_gain_and_cutoff():
    w, H = freqz(fp._PT_LP_TAPS, [1.0], worN=np.linspace(0, np.pi, 20001))
    f = w * FS / (2 * np.pi)
    mag = np.abs(H)
    assert mag[0] == pytest.approx(1.0, abs=1e-12)
    f3 = f[np.argmax(mag < 1 / np.sqrt(2))]
    assert f3 == pytest.approx(10.77, abs=0.02)
    # closed form Eq.2: |H| = sin^2(6 wT/2) / (36 sin^2(wT/2))
    wt = w[1:]
    ref = np.sin(3 * wt) ** 2 / (36 * np.sin(wt / 2) ** 2)
    np.testing.assert_allclose(mag[1:], np.abs(ref), atol=1e-10)


def test_lp_baseline_bug_is_real():
    """Dividing the whole recurrence by 36 is NOT the same filter."""
    x = impulse()
    good = lp_recursive(x)
    bad = lp_old_buggy(x)
    assert np.max(np.abs(good - bad)) > 1e-3
    dc = np.ones(400)
    assert lp_recursive(dc)[-1] == pytest.approx(1.0, abs=1e-9)
    assert abs(lp_old_buggy(dc)[-1] - 1.0) > 0.05


# ------------------------------- high-pass --------------------------------
def test_hp_fir_equals_recursion_impulse_and_random():
    rng = np.random.default_rng(1)
    for x in (impulse(), rng.standard_normal(500)):
        fir = np.convolve(x, fp._PT_HP_TAPS)[:len(x)]
        np.testing.assert_allclose(fir, hp_recursive(x), atol=1e-12)


def test_hp_impulse_response_shape_and_finite():
    h = fp._PT_HP_TAPS
    assert len(h) == 32
    assert h[16] == pytest.approx(1 - 1 / 32)
    assert np.all(np.delete(h, 16) == -1 / 32)
    # recursive impulse response terminates (pole at z=1 fully cancelled)
    r = hp_recursive(impulse(300))
    assert np.max(np.abs(r[32:])) < 1e-12


def test_hp_dc_zero_and_low_frequency_rejection():
    assert fp._PT_HP_TAPS.sum() == pytest.approx(0.0, abs=1e-14)
    w, H = freqz(fp._PT_HP_TAPS, [1.0], worN=np.linspace(0, np.pi, 20001))
    f = w * FS / (2 * np.pi)
    mag = np.abs(H)
    assert mag[0] < 1e-12
    assert mag[np.argmin(abs(f - 0.5))] < 0.05
    # transfer function check: z^-16 - (1/32)(1 - z^-32)/(1 - z^-1)
    z = np.exp(1j * w[1:])
    ref = z ** -16 - (1 / 32) * (1 - z ** -32) / (1 - 1 / z)
    np.testing.assert_allclose(H[1:], ref, atol=1e-10)
    print('HP |H| at 5/15/40 Hz', [float(mag[np.argmin(abs(f - q))]) for q in (5, 15, 40)])


def test_hp_printed_eq6_read_literally_has_pole_at_minus_one():
    """Documents the 1985 typo: literal Eq.6 does not decay; intended filter does."""
    lit = hp_printed_eq6(impulse(400))
    assert np.max(np.abs(lit[300:])) > 1e-3          # alternating, non-decaying tail
    ok = hp_recursive(impulse(400))
    assert np.max(np.abs(ok[300:])) < 1e-12


def test_combined_bandpass_response():
    b = np.convolve(fp._PT_LP_TAPS, fp._PT_HP_TAPS)
    w, H = freqz(b, [1.0], worN=np.linspace(0, np.pi, 40001))
    f = w * FS / (2 * np.pi)
    mag = np.abs(H)
    assert mag[0] < 1e-12
    band = f[mag >= mag.max() / np.sqrt(2)]
    assert band[0] == pytest.approx(4.9, abs=0.3)
    assert band[-1] == pytest.approx(11.8, abs=0.3)
    assert 5 < f[np.argmax(mag)] < 12


# ------------------------------- derivative --------------------------------
def test_derivative_impulse_matches_eq9_with_two_sample_delay():
    x = impulse(50, at=10)
    d = np.convolve(x, fp._PT_DERIV_TAPS)[:50]
    # Eq.9 (non-causal, T=1): y[n] = (1/8)(-x[n-2] - 2x[n-1] + 2x[n+1] + x[n+2])
    xe = np.concatenate([x, np.zeros(4)])
    eq9 = np.array([(-xe[n - 2] * (n >= 2) - 2 * xe[n - 1] * (n >= 1) + 2 * xe[n + 1] + xe[n + 2]) / 8
                    for n in range(50)])
    np.testing.assert_allclose(d[2:], eq9[:-2], atol=1e-14)   # exactly two samples later


def test_derivative_ramp_and_sinusoid():
    n = np.arange(200)
    d = np.convolve(3.0 * n, fp._PT_DERIV_TAPS)[:200]
    np.testing.assert_allclose(d[10:], 3.0, atol=1e-12)              # slope per sample
    for f0 in (1.0, 5.0, 20.0):
        w = 2 * np.pi * f0 / FS
        x = np.sin(w * n)
        d = np.convolve(x, fp._PT_DERIV_TAPS)[:200]
        gain_eq = (np.sin(2 * w) + 2 * np.sin(w)) / 4          # Eq.8 with T=1
        expected = gain_eq * np.cos(w * (n - 2))               # exact 2-sample delay
        np.testing.assert_allclose(d[10:], expected[10:], atol=1e-12)
    # near-linear (ideal-derivative) at 5 Hz
    w = 2 * np.pi * 5 / FS
    assert (np.sin(2 * w) + 2 * np.sin(w)) / 4 == pytest.approx(w, rel=0.02)


def test_derivative_baseline_coefficients_were_wrong():
    old = np.array([2.0, 1.0, 0.0, -1.0, -2.0]) / 8
    assert not np.allclose(old, fp._PT_DERIV_TAPS)


# ---------------------------- startup / chain ------------------------------
def test_chain_no_negative_index_wraparound():
    rng = np.random.default_rng(2)
    x = rng.standard_normal(600)
    y = x.copy(); y[-150:] += 1e6                # would leak into the start via x[-k]
    a = fp._pan_tompkins_filter_chain_200hz(x)
    b = fp._pan_tompkins_filter_chain_200hz(y)
    for u, v in zip(a, b):
        np.testing.assert_array_equal(u[:400], v[:400])    # causal: future cannot affect past


@pytest.mark.parametrize("c", [0.0, 1.0, -5.0, 1000.0])
def test_chain_constant_input_has_no_startup_transient(c):
    hp, d, integ = fp._pan_tompkins_filter_chain_200hz(np.full(400, c))
    assert np.max(np.abs(hp)) < 1e-9 * max(1, abs(c))
    assert np.max(np.abs(d)) < 1e-9 * max(1, abs(c))
    assert np.max(integ) < 1e-12 * max(1, c * c)


def test_chain_equals_unpadded_steady_state_and_preserves_input():
    rng = np.random.default_rng(3)
    x = rng.standard_normal(500) + 3.0
    x0 = x.copy()
    hp, d, integ = fp._pan_tompkins_filter_chain_200hz(x)
    np.testing.assert_array_equal(x, x0)              # input samples untouched
    # samples far from the start equal plain (zero-state) filtering of the input
    lp = np.convolve(x, fp._PT_LP_TAPS)[:500]
    ref = np.convolve(lp, fp._PT_HP_TAPS)[:500]
    np.testing.assert_allclose(hp[100:], ref[100:], atol=1e-10)


def test_processing_delay_measured():
    hp_d, int_d = fp._pan_tompkins_processing_delay_200hz()
    assert hp_d == 21            # LP 5 + HP 16
    assert int_d == 37           # + 2 (derivative) + (30-1)/2 window centre
