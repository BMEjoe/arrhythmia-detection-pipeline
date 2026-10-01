"""Phase 7: behaviour of the existing MIT-BIH helpers in final_pipeline.py that the
Phase 7 evaluation reuses unchanged (match_detected_peaks_to_annotations,
make_rr_windows, validate_r_peaks_against_annotations, _annotation_label) and of
the thin Phase 7 bookkeeping in experiments/phase7_mitbih/data.py.  Synthetic
inputs only; no MIT-BIH file is needed."""
import numpy as np
import pytest

import final_pipeline as fp
from experiments.phase7_mitbih import data as D

FS = 360.0
TOL = int(round(fp.ANNOTATION_MATCH_TOLERANCE_MS * FS / 1000.0))  # 27 samples


def test_tolerance_is_27_samples_at_360hz():
    assert TOL == 27


@pytest.mark.parametrize("sym,label", [
    ("N", 0), ("L", 0), ("R", 0), ("e", 0), ("j", 0),
    ("A", 1), ("a", 1), ("J", 1), ("S", 1), ("V", 1), ("E", 1), ("F", 1), ("f", 1),
    ("Q", 1),                      # in ABNORMAL and EXCLUDED; ABNORMAL takes effect
    ("/", -1), ("P", -1), ("B", -1), ("?", -1), ("x", -1),
    ("+", -1), ("~", -1), ("|", -1), ('"', -1), ("!", -1), ("[", -1), ("]", -1)])
def test_effective_label_of_every_symbol(sym, label):
    assert D.effective_label(sym) == label
    got = fp._annotation_label(sym)
    assert (got if got is not None else -1) == label


def test_q_is_in_both_sets_and_abnormal_wins():
    assert "Q" in fp.ABNORMAL_BEAT_SYMBOLS and "Q" in fp.EXCLUDED_BEAT_SYMBOLS
    assert fp._annotation_label("Q") == 1


def test_match_labels_perfect_and_one_to_one():
    ref = np.array([100, 400, 700, 1000])
    sym = np.array(["N", "V", "/", "N"])
    det = np.array([102, 395, 700, 1000, 1010])       # last is a second peak near 1000
    lab = fp.match_detected_peaks_to_annotations(det, ref, sym, FS)
    assert lab.tolist() == [0, 1, -1, 0, -1]          # paced -> -1; duplicate unmatched -> -1


def test_match_outside_tolerance_is_unmatched():
    ref = np.array([1000])
    lab = fp.match_detected_peaks_to_annotations(np.array([1000 + TOL + 1]), ref, np.array(["N"]), FS)
    assert lab.tolist() == [-1]
    lab = fp.match_detected_peaks_to_annotations(np.array([1000 + TOL]), ref, np.array(["N"]), FS)
    assert lab.tolist() == [0]


def test_non_beat_annotation_can_steal_a_match():
    """With all annotations passed (as analyze_mitbih_record_windows does), a detected
    peak nearer to a non-beat annotation than to its beat is matched to the non-beat
    one and labelled -1; with beat annotations only it gets the beat's label."""
    ref = np.array([105, 120])
    sym = np.array(["+", "N"])
    det = np.array([100])
    assert fp.match_detected_peaks_to_annotations(det, ref, sym, FS).tolist() == [-1]
    beat = np.isin(sym, list(D.BEAT_SYMBOLS))
    assert fp.match_detected_peaks_to_annotations(det, ref[beat], sym[beat], FS).tolist() == [0]


def test_greedy_match_agrees_with_fp_labels():
    rng = np.random.default_rng(0)
    ref = np.cumsum(rng.integers(150, 400, 300))
    sym = rng.choice(["N", "V", "A", "/", "+"], 300)
    det = np.sort(np.concatenate([ref[rng.random(300) > 0.05], rng.integers(0, ref[-1], 10)]))
    det = det + rng.integers(-30, 31, len(det))
    m = D.greedy_match(det, ref, FS)
    lab = fp.match_detected_peaks_to_annotations(det, ref, sym, FS)
    expect = np.array([-1 if j < 0 else D.effective_label(sym[j]) for j in m])
    assert np.array_equal(lab, expect)


def test_validate_counts():
    ref = np.array([100, 400, 700, 1000])
    det = np.array([100, 400, 1000, 1300])           # one missed (700), one extra (1300)
    v = fp.validate_r_peaks_against_annotations(det, ref, FS, 75.0)
    assert (v["TP"], v["FP"], v["FN"]) == (3, 1, 1)
    assert v["sensitivity"] == pytest.approx(0.75) and v["PPV"] == pytest.approx(0.75)


def _labels(n_beats, abnormal=(), bad=()):
    lab = np.zeros(n_beats, dtype=int)
    lab[list(abnormal)] = 1
    lab[list(bad)] = -1
    return lab


def test_make_rr_windows_alignment_and_threshold():
    n_rr = 512
    rr = np.arange(n_rr, dtype=float)
    # window 0: intervals 0..255 -> beats 1..256; 26 abnormal beats (>= 10 %)
    # window 1: intervals 256..511 -> beats 257..512; 25 abnormal beats (< 10 %)
    lab = _labels(n_rr + 1, abnormal=list(range(1, 27)) + list(range(257, 282)))
    w = fp.make_rr_windows(rr, lab)
    assert [x["start_rr"] for x in w] == [0, 256]
    assert w[0]["label"] == 1 and w[0]["abnormal_fraction"] == pytest.approx(26 / 256)
    assert w[1]["label"] == 0 and w[1]["abnormal_fraction"] == pytest.approx(25 / 256)
    assert np.array_equal(w[1]["rr"], rr[256:512])


def test_make_rr_windows_drops_any_unusable_ending_beat():
    rr = np.ones(512)
    w = fp.make_rr_windows(rr, _labels(513, bad=[300]))
    assert [x["start_rr"] for x in w] == [0]
    w = fp.make_rr_windows(rr, _labels(513, bad=[256]))   # last ending beat of window 0
    assert [x["start_rr"] for x in w] == [256]


def test_make_rr_windows_does_not_check_window_starting_beat():
    """Beat `start` begins interval rr[start] but its label is not checked: an
    unmatched peak at a window's first beat does not drop the window."""
    rr = np.ones(512)
    w = fp.make_rr_windows(rr, _labels(513, bad=[0]))      # starting beat of window 0
    assert [x["start_rr"] for x in w] == [0, 256]
    # beat 256 ends window 0 (dropped) and starts window 1 (kept)
    w = fp.make_rr_windows(rr, _labels(513, bad=[256]))
    assert [x["start_rr"] for x in w] == [256]


def test_missed_beat_gives_merged_interval_and_window_is_kept():
    """A beat the detector missed leaves no -1 label: the two RR intervals around it
    merge into one long interval labelled by the next detected beat."""
    ref = np.arange(1, 600) * 300                           # 599 annotated beats, 300 samples apart
    sym = np.array(["N"] * len(ref))
    sym[100] = "V"
    det = np.delete(ref, 100)                               # the V beat is missed
    lab = fp.match_detected_peaks_to_annotations(det, ref, sym, FS)
    assert np.all(lab >= 0) and np.sum(lab == 1) == 0      # the abnormal label is lost
    rr, _, _ = fp.extract_rr_intervals(det, FS, fp.CFG.rr_min_seconds, fp.CFG.rr_max_seconds)
    assert rr[99] == pytest.approx(600 / FS)                # merged interval
    w = fp.make_rr_windows(rr, lab)
    assert w[0]["start_rr"] == 0 and w[0]["label"] == 0


def test_make_rr_windows_length_check():
    with pytest.raises(ValueError):
        fp.make_rr_windows(np.ones(300), np.zeros(300, dtype=int))


def test_nn_touched_marks_interval_ending_at_and_following_abnormal():
    lab = np.array([0, 0, 1, 0, 0, 1, 1, 0, 0])            # beats 2, 5, 6 abnormal
    t = D.nn_touched(lab)                                   # 8 intervals
    assert t.tolist() == [False, True, True, False, True, True, True, False]


def test_edit_matches_phase6_rule():
    rr = np.array([1.0, 1.0, 0.6, 1.4, 1.0, 1.2])
    mask = np.array([False, False, True, True, False, False])
    out, info = D.edit_nn(rr, mask)
    assert out.tolist() == pytest.approx([1.0, 1.0, 1.0, 1.0, 1.0, 1.2])
    assert info["nn_fraction"] == pytest.approx(4 / 6)


def test_annotation_series_excludes_beats_after_vf_annotations():
    s = np.array([100, 400, 500, 600, 900, 1200])
    y = np.array(["N", "[", "!", "]", "N", "V"])
    bs, by, lab = D.annotation_beat_series(s, y)
    assert bs.tolist() == [100, 900, 1200]
    assert lab.tolist() == [0, -1, 1]


def test_hrv_values():
    rr = np.array([0.8, 0.9, 0.8, 0.84])
    h = D.hrv(rr)
    assert h["sdnn_ms"] == pytest.approx(np.std(rr, ddof=1) * 1000)
    assert h["rmssd_ms"] == pytest.approx(np.sqrt(np.mean(np.array([0.1, -0.1, 0.04]) ** 2)) * 1000)
    assert h["pnn50"] == pytest.approx(200 / 3)


def test_subjects_201_202_shared():
    assert D.SUBJECT["201"] == D.SUBJECT["202"]
    assert len(set(D.SUBJECT.values())) == 47 and len(D.SUBJECT) == 48
