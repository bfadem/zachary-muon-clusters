"""
Muon Shower Coincidence Analysis
Muhlenberg College 8-detector plastic scintillator array
Data: data_acquisition_2-27.csv (February 27, 2026 run)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import poisson, norm
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

DATA_FILE   = "data_acquisition_2-27.csv"
N_NANOS     = 8
VETO_US     = 20_000          # ±20 ms sync veto window
DROP_FIRST  = 2               # anomalous startup ST entries to drop per Nano
WINDOWS_MS  = [1, 2, 5, 10, 20, 50]
MIN_FOLD    = 2               # minimum multiplicity to record

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — PARSE
# ─────────────────────────────────────────────────────────────────────────────
def parse_csv(path):
    records      = []   # (nano_id, event_type, timestamp_us)
    master_syncs = []
    master_start = None

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if parts[0] == "MasterStart":
                master_start = int(parts[1])
            elif parts[0] == "MasterSync":
                master_syncs.append(int(parts[1]))
            elif len(parts) == 3:
                records.append((int(parts[0]), parts[1], int(parts[2])))

    return records, master_syncs, master_start


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — BUILD SYNC TABLE
# ─────────────────────────────────────────────────────────────────────────────
def build_sync_table(records):
    """
    Returns sync_array shape (n_syncs, 8):
      sync_array[k, i] = Nano (i+1) local timestamp at sync event k
    Also returns raw_st_per_nano dict {nano_id: sorted array of raw ST timestamps}
    for use in the veto step.
    """
    st_per_nano = {i: [] for i in range(1, N_NANOS + 1)}

    for (nano_id, etype, ts) in records:
        if etype == "ST":
            st_per_nano[nano_id].append(ts)

    # Sort and drop first DROP_FIRST anomalous entries
    for nano_id in st_per_nano:
        st_per_nano[nano_id] = sorted(st_per_nano[nano_id])

    raw_st_per_nano = {nano_id: np.array(st_per_nano[nano_id])
                       for nano_id in st_per_nano}

    for nano_id in st_per_nano:
        st_per_nano[nano_id] = st_per_nano[nano_id][DROP_FIRST:]

    # Trim all Nanos to same length
    n_syncs = min(len(v) for v in st_per_nano.values())
    sync_array = np.zeros((n_syncs, N_NANOS), dtype=np.int64)
    for i, nano_id in enumerate(range(1, N_NANOS + 1)):
        sync_array[:, i] = st_per_nano[nano_id][:n_syncs]

    return sync_array, raw_st_per_nano


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — INTERPOLATION CALIBRATION
# ─────────────────────────────────────────────────────────────────────────────
def build_offset_table(sync_array):
    """
    offset_table[k, i] = sync_array[k, i] - sync_array[k, 0]
    (clock offset of Nano i+1 relative to Nano 1 at sync k)
    """
    return sync_array - sync_array[:, 0:1]


def calibrate_timestamp(t_nano, nano_idx, sync_array, offset_table):
    """
    Convert a single raw Nano timestamp to Nano-1 reference time via
    piecewise-linear interpolation between consecutive sync pairs.
    nano_idx: 0-based (0 = Nano 1)
    """
    st = sync_array[:, nano_idx]
    ot = offset_table[:, nano_idx]

    # Find bracketing sync interval
    k = np.searchsorted(st, t_nano, side="right") - 1

    # Clamp to valid range (extrapolate with first/last pair at edges)
    k = np.clip(k, 0, len(st) - 2)

    dt_st = st[k + 1] - st[k]
    if dt_st == 0:
        frac = 0.0
    else:
        frac = (t_nano - st[k]) / dt_st

    offset = ot[k] + frac * (ot[k + 1] - ot[k])
    return t_nano - offset


def calibrate_all_pulses(records, sync_array, offset_table):
    """
    Returns a list of (t_corr_us, nano_id) for all PT events.
    """
    calibrated = []
    for (nano_id, etype, ts) in records:
        if etype != "PT":
            continue
        nano_idx = nano_id - 1
        t_corr = calibrate_timestamp(ts, nano_idx, sync_array, offset_table)
        calibrated.append((t_corr, nano_id))
    return calibrated


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — PER-NANO SYNC VETO
# ─────────────────────────────────────────────────────────────────────────────
def apply_sync_veto(calibrated, raw_st_per_nano, veto_us=VETO_US):
    """
    Remove PT events within ±veto_us of any raw ST timestamp for that Nano.
    Returns (kept_events, n_vetoed).
    """
    kept    = []
    vetoed  = 0

    for (t_corr, nano_id) in calibrated:
        raw_ts_for_nano = raw_st_per_nano[nano_id]
        # Find the closest ST timestamp
        idx = np.searchsorted(raw_ts_for_nano, t_corr)
        near = False
        for j in [idx - 1, idx]:
            if 0 <= j < len(raw_ts_for_nano):
                if abs(t_corr - raw_ts_for_nano[j]) < veto_us:
                    near = True
                    break
        if near:
            vetoed += 1
        else:
            kept.append((t_corr, nano_id))

    return kept, vetoed


def effective_run_duration(raw_st_per_nano, veto_us=VETO_US):
    """
    T_eff per Nano = raw run duration − (n_ST × 2 × veto_us)
    Returns dict {nano_id: T_eff_s} and overall T_eff_s (mean across Nanos).
    """
    t_eff = {}
    for nano_id, sts in raw_st_per_nano.items():
        if len(sts) < 2:
            continue
        raw_dur_s   = (sts[-1] - sts[0]) / 1e6
        veto_total_s = len(sts) * 2 * veto_us / 1e6
        t_eff[nano_id] = max(raw_dur_s - veto_total_s, 0.0)
    return t_eff


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — COINCIDENCE COUNTING
# ─────────────────────────────────────────────────────────────────────────────
def find_coincidences(kept_events, window_us, min_fold=MIN_FOLD):
    """
    Greedy sliding-window coincidence search.
    Each pulse used at most once (earliest-first).
    Returns list of clusters: each cluster is a list of (t_corr, nano_id).
    """
    # Sort by calibrated time
    events = sorted(kept_events, key=lambda x: x[0])
    t_arr  = np.array([e[0] for e in events])
    n_arr  = np.array([e[1] for e in events])

    used     = np.zeros(len(events), dtype=bool)
    clusters = []

    for i in range(len(events)):
        if used[i]:
            continue
        t0 = t_arr[i]

        # Find all events within window
        j_end = np.searchsorted(t_arr, t0 + window_us, side="right")
        candidates = list(range(i, j_end))

        # Pick at most one per Nano (first occurrence)
        seen_nanos = {}
        for j in candidates:
            if not used[j]:
                nid = n_arr[j]
                if nid not in seen_nanos:
                    seen_nanos[nid] = j

        if len(seen_nanos) >= min_fold:
            cluster_indices = list(seen_nanos.values())
            for j in cluster_indices:
                used[j] = True
            clusters.append([(t_arr[j], n_arr[j]) for j in cluster_indices])

    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — EXPECTED ACCIDENTALS
# ─────────────────────────────────────────────────────────────────────────────
def expected_accidentals(window_us, fold, rates, run_duration_s):
    """
    E_k = Σ_{k-tuples} Π r_i × (2τ)^(k-1) × T_eff
    rates: dict {nano_id: rate_hz}
    """
    tau     = window_us / 1e6
    nano_ids = sorted(rates.keys())
    total    = 0.0
    for combo in combinations(nano_ids, fold):
        total += np.prod([rates[n] for n in combo]) * (2 * tau) ** (fold - 1) * run_duration_s
    return total


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7 — STATISTICAL SIGNIFICANCE
# ─────────────────────────────────────────────────────────────────────────────
def significance(observed, expected):
    """
    p = P(X >= observed | Poisson(expected))
    sigma = norm.ppf(1 - p)
    """
    if expected <= 0:
        return np.nan
    p = 1.0 - poisson.cdf(observed - 1, expected)
    p = max(p, 1e-15)   # floor to avoid -inf
    return norm.ppf(1.0 - p)


# ─────────────────────────────────────────────────────────────────────────────
# HELPER — PER-DETECTOR RATES
# ─────────────────────────────────────────────────────────────────────────────
def compute_rates(kept_events, t_eff_per_nano):
    rates = {}
    counts = {}
    for (_, nano_id) in kept_events:
        counts[nano_id] = counts.get(nano_id, 0) + 1
    for nano_id, cnt in counts.items():
        t = t_eff_per_nano.get(nano_id, 1.0)
        rates[nano_id] = cnt / t if t > 0 else 0.0
    return rates, counts


# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────
def plot_detector_rates(kept_events, t_eff_per_nano, rates):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))

    # --- bar chart ---
    nano_ids = sorted(rates.keys())
    ax = axes[0]
    ax.bar([str(n) for n in nano_ids], [rates[n] for n in nano_ids],
           color="steelblue", edgecolor="k")
    ax.set_xlabel("Nano ID")
    ax.set_ylabel("Rate (Hz)")
    ax.set_title("Per-Detector Pulse Rates")
    ax.axhline(np.mean(list(rates.values())), color="red",
               linestyle="--", label="mean")
    ax.legend()

    # --- rate vs time (60 s bins) ---
    ax = axes[1]
    all_t = sorted([t for (t, _) in kept_events])
    t_min, t_max = all_t[0], all_t[-1]
    bin_s = 60
    bins  = np.arange(t_min, t_max + bin_s * 1e6, bin_s * 1e6)

    colors = plt.cm.tab10(np.linspace(0, 1, N_NANOS))
    for i, nano_id in enumerate(nano_ids):
        t_nano = np.array([t for (t, n) in kept_events if n == nano_id])
        counts, edges = np.histogram(t_nano, bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2 / 1e6 / 60  # convert to minutes
        ax.plot(centers, counts / bin_s, color=colors[i],
                label=f"Nano {nano_id}", alpha=0.8)

    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Rate (Hz)")
    ax.set_title("Pulse Rate vs. Time (60 s bins)")
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig("plot_01_detector_rates.png", dpi=150)
    plt.close()
    print("Saved plot_01_detector_rates.png")


def plot_clock_offsets(sync_array, offset_table):
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    axes = axes.flatten()

    sync_indices = np.arange(sync_array.shape[0])
    # Use Nano-1 times as x-axis (minutes)
    x_min = sync_array[:, 0] / 1e6 / 60

    for i in range(1, N_NANOS):   # skip Nano 1 (reference)
        ax = axes[i - 1]
        offsets_us = offset_table[:, i]

        ax.scatter(x_min, offsets_us / 1e3, s=6, alpha=0.6, label=f"Nano {i+1}")

        # Linear fit
        coeffs = np.polyfit(x_min, offsets_us / 1e3, 1)
        fit    = np.polyval(coeffs, x_min)
        ax.plot(x_min, fit, "r-", linewidth=1.5,
                label=f"slope={coeffs[0]:.1f} ms/min")

        ax.set_title(f"Nano {i+1} − Nano 1")
        ax.set_ylabel("Offset (ms)")
        ax.legend(fontsize=7)

    axes[-1].set_visible(False)   # 8 Nanos → 7 pairs → one blank
    for ax in axes[5:7]:
        ax.set_xlabel("Time (min)")

    plt.suptitle("Pairwise Clock Offsets (Nano_i − Nano_1)", fontsize=13)
    plt.tight_layout()
    plt.savefig("plot_02_clock_offsets.png", dpi=150)
    plt.close()
    print("Saved plot_02_clock_offsets.png")


def plot_calibration_residuals(sync_array, offset_table):
    """
    Residual = change in offset over consecutive 10-s sync intervals.
    This sets the timing resolution floor.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    for i in range(1, N_NANOS):
        diff = np.diff(offset_table[:, i])   # µs change per 10 s interval
        ax.plot(diff / 1e3, alpha=0.7, label=f"Nano {i+1}")

    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xlabel("Sync interval index")
    ax.set_ylabel("Δoffset per 10 s interval (ms)")
    ax.set_title("Calibration Residuals — Timing Resolution Floor")
    ax.legend(fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig("plot_03_calibration_residuals.png", dpi=150)
    plt.close()
    print("Saved plot_03_calibration_residuals.png")


def plot_obs_vs_exp(results_by_window, windows_ms):
    """
    results_by_window: {window_ms: {fold: (observed, expected)}}
    """
    folds   = [2, 3, 4, 5, 6]
    colors  = plt.cm.plasma(np.linspace(0.1, 0.9, len(folds)))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # --- Observed vs Expected (log-log) ---
    ax = axes[0]
    for fi, fold in enumerate(folds):
        obs_list = []
        exp_list = []
        for w in windows_ms:
            if fold in results_by_window[w]:
                o, e = results_by_window[w][fold]
                obs_list.append(o)
                exp_list.append(e)
        if exp_list:
            ax.scatter(exp_list, obs_list, color=colors[fi],
                       label=f"{fold}-fold", s=60, zorder=3)

    lim_min = min(0.01, min(e for w in results_by_window.values()
                             for (o, e) in w.values() if e > 0))
    lim_max = max(max(o for w in results_by_window.values()
                       for (o, e) in w.values()) * 2, 10)
    ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--",
            linewidth=1, label="expected = observed")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Expected accidentals")
    ax.set_ylabel("Observed coincidences")
    ax.set_title("Observed vs. Expected Coincidences")
    ax.legend(fontsize=8)

    # --- Significance vs window size ---
    ax = axes[1]
    for fi, fold in enumerate(folds):
        sigs = []
        for w in windows_ms:
            if fold in results_by_window[w]:
                o, e = results_by_window[w][fold]
                sigs.append(significance(o, e))
            else:
                sigs.append(np.nan)
        ax.plot(windows_ms, sigs, "o-", color=colors[fi],
                label=f"{fold}-fold", linewidth=1.5)

    ax.axhline(3, color="gray", linestyle="--", linewidth=0.8, label="3σ")
    ax.axhline(5, color="red",  linestyle="--", linewidth=0.8, label="5σ")
    ax.set_xscale("log")
    ax.set_xlabel("Window size (ms)")
    ax.set_ylabel("Significance (σ)")
    ax.set_title("Significance vs. Window Size")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("plot_04_obs_vs_exp_and_significance.png", dpi=150)
    plt.close()
    print("Saved plot_04_obs_vs_exp_and_significance.png")


def plot_coincidence_timeline(clusters_by_window, window_ms=1):
    clusters = clusters_by_window[window_ms]
    if not clusters:
        print(f"No coincidences found at {window_ms} ms window.")
        return

    fig, ax = plt.subplots(figsize=(14, 5))

    fold_colors = {2: "blue", 3: "green", 4: "orange", 5: "red",
                   6: "purple", 7: "brown", 8: "black"}

    for cluster in clusters:
        t_min_min = min(t for (t, _) in cluster) / 1e6 / 60
        fold      = len(cluster)
        color     = fold_colors.get(fold, "gray")
        ax.scatter([t_min_min], [fold], color=color, alpha=0.6,
                   s=20, zorder=2)

    # Legend
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=c, markersize=8, label=f"{f}-fold")
               for f, c in fold_colors.items() if f >= 2]
    ax.legend(handles=handles, fontsize=8)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Multiplicity (fold)")
    ax.set_title(f"Coincidence Timeline ({window_ms} ms window)")

    plt.tight_layout()
    plt.savefig("plot_05_coincidence_timeline.png", dpi=150)
    plt.close()
    print("Saved plot_05_coincidence_timeline.png")


def plot_shower_candidate(clusters_by_window, window_ms=1, top_n=3):
    """
    Timestamp ladder for the highest-multiplicity clusters.
    """
    clusters = clusters_by_window[window_ms]
    if not clusters:
        return

    # Sort by multiplicity desc, then span asc
    def sort_key(c):
        fold = len(c)
        span = max(t for (t, _) in c) - min(t for (t, _) in c)
        return (-fold, span)

    top = sorted(clusters, key=sort_key)[:top_n]

    fig, axes = plt.subplots(1, len(top), figsize=(5 * len(top), 5))
    if len(top) == 1:
        axes = [axes]

    for ax, cluster in zip(axes, top):
        t0   = min(t for (t, _) in cluster)
        fold = len(cluster)
        span = max(t for (t, _) in cluster) - t0
        t0_min = t0 / 1e6 / 60

        for (t, nano_id) in sorted(cluster, key=lambda x: x[1]):
            dt_us = t - t0
            ax.plot([dt_us, dt_us], [nano_id - 0.4, nano_id + 0.4],
                    "b-", linewidth=3)
            ax.scatter([dt_us], [nano_id], color="red", zorder=5, s=60)
            ax.text(dt_us + 0.1, nano_id + 0.1, f"{dt_us:.1f} µs",
                    fontsize=8)

        ax.set_yticks(sorted(nano_id for (_, nano_id) in cluster))
        ax.set_yticklabels([f"Nano {n}" for (_, n) in
                            sorted(cluster, key=lambda x: x[1])])
        ax.set_xlabel("Δt from first hit (µs)")
        ax.set_title(f"{fold}-fold @ {t0_min:.2f} min\nspan={span:.1f} µs")
        ax.grid(axis="x", alpha=0.3)

    plt.suptitle(f"Top {top_n} High-Multiplicity Events ({window_ms} ms window)",
                 fontsize=12)
    plt.tight_layout()
    plt.savefig("plot_06_shower_candidates.png", dpi=150)
    plt.close()
    print("Saved plot_06_shower_candidates.png")


# ─────────────────────────────────────────────────────────────────────────────
# DIAGNOSTIC PLOTS — TIMING VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def plot_diag_A_sync_intervals(sync_array):
    """
    Plot A: Sync interval regularity.
    ST[k+1] - ST[k] per Nano, as a time series and histogram.
    Proves the sync pulses are evenly spaced and that no pulses are missed.
    """
    intervals = np.diff(sync_array, axis=0)          # shape (n_syncs-1, 8), µs
    deviation_ms = (intervals / 1e6 - 10.0) * 1e3    # deviation from 10 s, in ms

    colors = plt.cm.tab10(np.linspace(0, 1, N_NANOS))
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Time series
    ax = axes[0]
    x = np.arange(intervals.shape[0])
    for i in range(N_NANOS):
        ax.plot(x, intervals[:, i] / 1e6, color=colors[i],
                alpha=0.7, label=f"Nano {i+1}", linewidth=0.8)
    ax.axhline(10.0, color="k", linestyle="--", linewidth=1, label="10 s nominal")
    ax.set_xlabel("Sync interval index")
    ax.set_ylabel("Interval (s)")
    ax.set_title("Sync Interval Duration per Nano  [ST(k+1) − ST(k)]")
    ax.legend(fontsize=7, ncol=4)

    # Histogram of deviations from 10 s
    ax = axes[1]
    all_devs = deviation_ms.flatten()
    lo, hi = np.percentile(all_devs, 0.5), np.percentile(all_devs, 99.5)
    bins = np.linspace(lo - 1, hi + 1, 60)
    for i in range(N_NANOS):
        ax.hist(deviation_ms[:, i], bins=bins, alpha=0.4,
                color=colors[i], label=f"Nano {i+1}")
    ax.axvline(0, color="k", linestyle="--", linewidth=1)
    ax.set_xlabel("Deviation from 10 s nominal (ms)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Sync Interval Deviations")
    ax.legend(fontsize=7, ncol=4)

    plt.tight_layout()
    plt.savefig("diag_A_sync_intervals.png", dpi=150)
    plt.close()
    rms_all = np.sqrt(np.mean(all_devs**2))
    print(f"Saved diag_A_sync_intervals.png  (RMS deviation = {rms_all:.2f} ms)")


def plot_diag_B_sync_jitter(sync_array, offset_table):
    """
    Plot B: Jitter in simultaneous sync receipt.
    Removes the smooth clock-drift trend from each Nano's offset series;
    the residual is the jitter in sync detection.  Its RMS is the irreducible
    timing floor that limits calibration accuracy regardless of interpolation.
    """
    x = np.arange(sync_array.shape[0])
    detrended = np.zeros_like(offset_table, dtype=float)
    for i in range(1, N_NANOS):
        coeffs = np.polyfit(x, offset_table[:, i], 1)
        detrended[:, i] = offset_table[:, i] - np.polyval(coeffs, x)

    colors = plt.cm.tab10(np.linspace(0, 1, N_NANOS - 1))
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Time series of detrended offsets
    ax = axes[0]
    for i in range(1, N_NANOS):
        ax.plot(x, detrended[:, i] / 1e3, color=colors[i - 1],
                alpha=0.7, label=f"Nano {i+1}", linewidth=0.8)
    ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Sync index")
    ax.set_ylabel("Detrended offset (ms)")
    ax.set_title("Sync Detection Jitter — Offset After Removing Linear Drift Trend")
    ax.legend(fontsize=7, ncol=4)

    # Histogram
    ax = axes[1]
    all_ms = detrended[:, 1:].flatten() / 1e3
    rms_all = np.sqrt(np.mean(all_ms**2))
    lo, hi = np.percentile(all_ms, 0.5), np.percentile(all_ms, 99.5)
    bins = np.linspace(lo, hi, 50)
    for i in range(1, N_NANOS):
        vals = detrended[:, i] / 1e3
        rms_i = np.sqrt(np.mean(vals**2))
        ax.hist(vals, bins=bins, alpha=0.4, color=colors[i - 1],
                label=f"Nano {i+1} (σ={rms_i:.2f} ms)")
    ax.axvline(0, color="k", linestyle="--")
    ax.set_xlabel("Detrended offset (ms)")
    ax.set_ylabel("Count")
    ax.set_title(f"Sync Jitter Distribution  (RMS all Nanos = {rms_all:.2f} ms)")
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig("diag_B_sync_jitter.png", dpi=150)
    plt.close()
    print(f"Saved diag_B_sync_jitter.png  (RMS jitter = {rms_all:.2f} ms)")


def plot_diag_C_oos_residuals(sync_array, offset_table):
    """
    Plot C: Out-of-sample calibration residuals (even/odd sync split).
    Calibration is built from even-indexed syncs only; odd-indexed syncs are
    held out and their offsets are predicted by interpolation.
    Residual = predicted - actual.  RMS is the true achieved timing accuracy
    for PT events falling between sync pulses.
    """
    even_idx = np.arange(0, sync_array.shape[0], 2)
    odd_idx  = np.arange(1, sync_array.shape[0], 2)

    if len(even_idx) < 3 or len(odd_idx) < 1:
        print("Not enough sync events for even/odd split; skipping Plot C.")
        return

    sync_train   = sync_array[even_idx]
    offset_train = offset_table[even_idx]
    sync_test    = sync_array[odd_idx]
    offset_test  = offset_table[odd_idx]

    n_test    = sync_test.shape[0]
    residuals = np.zeros((n_test, N_NANOS))   # µs

    for i in range(1, N_NANOS):
        st_tr = sync_train[:, i]
        ot_tr = offset_train[:, i]
        for idx in range(n_test):
            t  = sync_test[idx, i]
            k  = np.clip(np.searchsorted(st_tr, t, side="right") - 1,
                         0, len(st_tr) - 2)
            dt = st_tr[k + 1] - st_tr[k]
            frac = (t - st_tr[k]) / dt if dt != 0 else 0.0
            predicted = ot_tr[k] + frac * (ot_tr[k + 1] - ot_tr[k])
            residuals[idx, i] = predicted - offset_test[idx, i]

    colors = plt.cm.tab10(np.linspace(0, 1, N_NANOS - 1))
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Time series
    ax = axes[0]
    for i in range(1, N_NANOS):
        ax.plot(odd_idx, residuals[:, i] / 1e3, ".",
                color=colors[i - 1], alpha=0.6, markersize=4,
                label=f"Nano {i+1}")
    ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Sync index (odd = held-out test points)")
    ax.set_ylabel("Predicted − Actual offset (ms)")
    ax.set_title("Out-of-Sample Calibration Residuals  (train = even syncs, test = odd syncs)")
    ax.legend(fontsize=7, ncol=4)

    # Histogram
    ax = axes[1]
    all_ms = residuals[:, 1:].flatten() / 1e3
    overall_rms = np.sqrt(np.mean(all_ms**2))
    lo, hi = np.percentile(all_ms, 0.5), np.percentile(all_ms, 99.5)
    bins = np.linspace(lo, hi, 50)
    for i in range(1, N_NANOS):
        vals = residuals[:, i] / 1e3
        rms_i = np.sqrt(np.mean(vals**2))
        ax.hist(vals, bins=bins, alpha=0.4, color=colors[i - 1],
                label=f"Nano {i+1} (σ={rms_i:.1f} ms)")
    ax.axvline(0, color="k", linestyle="--")
    ax.set_xlabel("Predicted − Actual offset (ms)")
    ax.set_ylabel("Count")
    ax.set_title(f"Out-of-Sample Residual Distribution  (RMS all = {overall_rms:.2f} ms)")
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig("diag_C_oos_residuals.png", dpi=150)
    plt.close()
    print(f"Saved diag_C_oos_residuals.png  (overall OOS RMS = {overall_rms:.2f} ms)")


def plot_diag_D_max_interp_error(sync_array, offset_table):
    """
    Plot D: Worst-case interpolation error per 10-second sync interval.
    For a PT event at the midpoint of interval k, the maximum error
    from piecewise-linear interpolation is |Δoffset[k]| / 4.
    """
    d_offset    = np.diff(offset_table, axis=0)   # µs change per interval
    max_err_ms  = np.abs(d_offset) / 4.0 / 1e3    # ms, worst case at midpoint
    x_min       = sync_array[:-1, 0] / 1e6 / 60   # minutes

    colors = plt.cm.tab10(np.linspace(0, 1, N_NANOS - 1))
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))

    # Time series
    ax = axes[0]
    for i in range(1, N_NANOS):
        ax.plot(x_min, max_err_ms[:, i], color=colors[i - 1],
                alpha=0.7, label=f"Nano {i+1}", linewidth=0.8)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Max interpolation error (ms)")
    ax.set_title("Worst-Case Timing Error for a PT Event at Interval Midpoint  [|Δoffset| / 4]")
    ax.legend(fontsize=7, ncol=4)

    # Histogram with median and 95th percentile marked
    ax = axes[1]
    all_ms = max_err_ms[:, 1:].flatten()
    med    = np.median(all_ms)
    p95    = np.percentile(all_ms, 95)
    bins   = np.linspace(0, np.percentile(all_ms, 99.5), 50)
    for i in range(1, N_NANOS):
        ax.hist(max_err_ms[:, i], bins=bins, alpha=0.4,
                color=colors[i - 1], label=f"Nano {i+1}")
    ax.axvline(med, color="k",   linestyle="--", label=f"median = {med:.2f} ms")
    ax.axvline(p95, color="red", linestyle="--", label=f"95th pct = {p95:.2f} ms")
    ax.set_xlabel("Max interpolation error (ms)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Worst-Case Interpolation Error per Interval")
    ax.legend(fontsize=7, ncol=2)

    plt.tight_layout()
    plt.savefig("diag_D_max_interp_error.png", dpi=150)
    plt.close()
    print(f"Saved diag_D_max_interp_error.png  "
          f"(median={med:.2f} ms, 95th pct={p95:.2f} ms)")


def plot_diag_E_dt_distribution(kept_events, window_us=1000):
    """
    Plot E: Distribution of time differences within all cross-Nano pairs
    found within the coincidence window.  A genuine signal from correlated
    events (showers) piles up near Δt = 0; accidental pairs are uniformly
    distributed across the window.
    """
    events = sorted(kept_events, key=lambda x: x[0])
    t_arr  = np.array([e[0] for e in events])
    n_arr  = np.array([e[1] for e in events])

    dt_values = []
    for i in range(len(events)):
        t0    = t_arr[i]
        j_end = np.searchsorted(t_arr, t0 + window_us, side="right")
        for j in range(i + 1, j_end):
            if n_arr[j] != n_arr[i]:
                dt_values.append(t_arr[j] - t_arr[i])

    dt_values = np.array(dt_values)
    if len(dt_values) == 0:
        print("No cross-Nano pairs found; skipping Plot E.")
        return

    # Expected flat level per bin if all pairs were accidental
    n_bins      = 50
    bin_width   = window_us / n_bins
    flat_level  = len(dt_values) * bin_width / window_us

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Full window
    ax    = axes[0]
    bins  = np.linspace(0, window_us, n_bins + 1)
    cnts, edges = np.histogram(dt_values, bins=bins)
    ctrs  = (edges[:-1] + edges[1:]) / 2
    ax.bar(ctrs / 1e3, cnts, width=bin_width / 1e3,
           color="steelblue", edgecolor="none", alpha=0.8)
    ax.axhline(flat_level, color="red", linestyle="--",
               label=f"flat background ({flat_level:.0f}/bin)")
    ax.set_xlabel("Δt (ms)")
    ax.set_ylabel("Pair count")
    ax.set_title(f"Δt Distribution — All Cross-Nano Pairs  ({window_us/1e3:.0f} ms window)")
    ax.legend()

    # Zoomed to first 10% of window
    ax       = axes[1]
    zoom_us  = window_us * 0.10
    mask     = dt_values <= zoom_us
    bins_z   = np.linspace(0, zoom_us, n_bins + 1)
    cnts_z, edges_z = np.histogram(dt_values[mask], bins=bins_z)
    ctrs_z   = (edges_z[:-1] + edges_z[1:]) / 2
    bw_z     = zoom_us / n_bins
    flat_z   = len(dt_values) * bw_z / window_us
    ax.bar(ctrs_z / 1e3, cnts_z, width=bw_z / 1e3,
           color="steelblue", edgecolor="none", alpha=0.8)
    ax.axhline(flat_z, color="red", linestyle="--",
               label=f"flat ({flat_z:.1f}/bin)")
    ax.set_xlabel("Δt (ms)")
    ax.set_ylabel("Pair count")
    ax.set_title(f"Δt Distribution — Zoomed to First {zoom_us/1e3:.2f} ms")
    ax.legend()

    plt.tight_layout()
    plt.savefig("diag_E_dt_distribution.png", dpi=150)
    plt.close()
    print(f"Saved diag_E_dt_distribution.png  (total pairs = {len(dt_values)})")


def plot_diag_F_null_test(kept_events, rates, T_eff_s, window_us=1000,
                          n_trials=200):
    """
    Plot F: Null test via independent random time-shifting of each Nano.
    Each trial shifts every Nano's events by a different random offset
    (uniformly drawn, wrapping modulo run duration) and re-runs the
    coincidence search.  The distribution of shuffle counts is compared
    to the real data and to the Poisson accidental expectation.
    """
    print(f"  Running {n_trials} time-shift trials  "
          f"(window = {window_us/1e3:.0f} ms) ...")

    events_by_nano = {}
    for (t, nano_id) in kept_events:
        events_by_nano.setdefault(nano_id, []).append(t)

    t_min_g       = min(t for (t, _) in kept_events)
    t_max_g       = max(t for (t, _) in kept_events)
    run_dur_us    = t_max_g - t_min_g

    # Real coincidence counts at this window
    real_clusters = find_coincidences(kept_events, window_us)
    real_fold = {}
    for c in real_clusters:
        f = len(c)
        real_fold[f] = real_fold.get(f, 0) + 1

    folds_to_track = [2, 3, 4, 5]
    null_counts    = {f: [] for f in folds_to_track}
    rng = np.random.default_rng(seed=42)

    for trial in range(n_trials):
        if (trial + 1) % 50 == 0:
            print(f"    trial {trial + 1}/{n_trials}")
        shifted = []
        for nano_id, times in events_by_nano.items():
            # Shift in (10%, 90%) of run to avoid edge effects
            shift = rng.uniform(run_dur_us * 0.1, run_dur_us * 0.9)
            for t in times:
                t_shift = t_min_g + (t - t_min_g + shift) % run_dur_us
                shifted.append((t_shift, nano_id))
        clusters = find_coincidences(shifted, window_us)
        fc = {}
        for c in clusters:
            f = len(c)
            fc[f] = fc.get(f, 0) + 1
        for f in folds_to_track:
            null_counts[f].append(fc.get(f, 0))

    fig, axes = plt.subplots(1, len(folds_to_track), figsize=(16, 5))
    for ax, fold in zip(axes, folds_to_track):
        null     = np.array(null_counts[fold])
        real_val = real_fold.get(fold, 0)
        exp_val  = expected_accidentals(window_us, fold, rates, T_eff_s)
        sig_val  = significance(real_val, exp_val)
        p_shuffle = np.mean(null >= real_val)

        lo = min(null.min(), real_val) - 0.5
        hi = max(null.max(), real_val) + 1.5
        bins = np.arange(lo, hi, 1.0)
        ax.hist(null, bins=bins, color="gray", edgecolor="k",
                alpha=0.7, label=f"Shuffled (n={n_trials})")
        ax.axvline(real_val, color="red",  linewidth=2,
                   label=f"Real data = {real_val}")
        ax.axvline(exp_val,  color="blue", linewidth=1.5, linestyle="--",
                   label=f"Poisson exp = {exp_val:.1f}")
        ax.set_xlabel("Coincidence count")
        ax.set_ylabel("Trials")
        ax.set_title(
            f"{fold}-fold\n"
            f"Real={real_val},  Exp={exp_val:.1f},  {sig_val:+.1f}σ\n"
            f"p(shuffle ≥ real) = {p_shuffle:.3f}"
        )
        ax.legend(fontsize=7)

    plt.suptitle(
        f"Null Test: Time-Shifted Coincidence Counts  "
        f"({window_us/1e3:.0f} ms window, {n_trials} trials)",
        fontsize=12
    )
    plt.tight_layout()
    plt.savefig("diag_F_null_test.png", dpi=150)
    plt.close()
    print("Saved diag_F_null_test.png")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Muon Shower Coincidence Analysis")
    print("=" * 60)

    # ── Step 1: Parse ────────────────────────────────────────────────────────
    records, master_syncs, master_start = parse_csv(DATA_FILE)
    print(f"\n[Parse] {len(records)} records, {len(master_syncs)} MasterSync events")

    # ── Step 2: Sync table ───────────────────────────────────────────────────
    sync_array, raw_st_per_nano = build_sync_table(records)
    print(f"[Sync]  sync_array shape: {sync_array.shape}  "
          f"(dropped first {DROP_FIRST} ST per Nano)")

    # ── Step 3: Calibration ──────────────────────────────────────────────────
    offset_table = build_offset_table(sync_array)
    calibrated   = calibrate_all_pulses(records, sync_array, offset_table)
    print(f"[Calib] {len(calibrated)} PT events calibrated to Nano-1 reference")

    # ── Step 4: Sync veto ────────────────────────────────────────────────────
    kept_events, n_vetoed = apply_sync_veto(calibrated, raw_st_per_nano)
    t_eff_per_nano = effective_run_duration(raw_st_per_nano)
    T_eff_s  = np.mean(list(t_eff_per_nano.values()))
    print(f"[Veto]  {n_vetoed} events vetoed  |  {len(kept_events)} events kept")
    print(f"        Effective run duration: {T_eff_s/60:.2f} min  "
          f"({T_eff_s:.1f} s)")

    # ── Per-detector rates ───────────────────────────────────────────────────
    rates, counts = compute_rates(kept_events, t_eff_per_nano)
    print("\n[Rates]")
    print(f"  {'Nano':>6}  {'Count':>7}  {'T_eff (s)':>10}  {'Rate (Hz)':>10}")
    for nano_id in sorted(rates.keys()):
        print(f"  {nano_id:>6}  {counts.get(nano_id,0):>7}  "
              f"{t_eff_per_nano.get(nano_id, 0):>10.1f}  "
              f"{rates[nano_id]:>10.4f}")

    # Clock drift summary (ppm relative to Nano 1)
    print("\n[Clock drift relative to Nano 1]")
    x_min = sync_array[:, 0] / 1e6 / 60
    duration_s = (sync_array[-1, 0] - sync_array[0, 0]) / 1e6
    for i in range(1, N_NANOS):
        total_drift_us = offset_table[-1, i] - offset_table[0, i]
        ppm = total_drift_us / duration_s
        print(f"  Nano {i+1}: total drift = {total_drift_us/1e3:.1f} ms  "
              f"({ppm:+.1f} ppm)")

    # ── Steps 5–7: Coincidences, accidentals, significance ──────────────────
    print("\n[Coincidences]")
    print(f"  {'Window':>8}  "
          + "  ".join(f"{f}-fold obs/exp/σ" for f in range(2, 7)))

    results_by_window = {}
    clusters_by_window = {}

    for w_ms in WINDOWS_MS:
        w_us = w_ms * 1000
        clusters = find_coincidences(kept_events, w_us)
        clusters_by_window[w_ms] = clusters

        # Count by fold
        fold_counts = {}
        for cluster in clusters:
            f = len(cluster)
            fold_counts[f] = fold_counts.get(f, 0) + 1

        results_by_window[w_ms] = {}
        row = f"  {w_ms:>5} ms  "
        for fold in range(2, 7):
            obs = fold_counts.get(fold, 0)
            exp = expected_accidentals(w_us, fold, rates, T_eff_s)
            sig = significance(obs, exp)
            results_by_window[w_ms][fold] = (obs, exp)
            row += f"  {obs}/{exp:.1f}/{sig:+.1f}σ"
        print(row)

    # ── Primary shower candidate detail ─────────────────────────────────────
    print("\n[Shower Candidates — top events at 1 ms window]")
    clusters_1ms = clusters_by_window[1]

    def sort_key(c):
        return (-len(c), max(t for (t, _) in c) - min(t for (t, _) in c))

    top_clusters = sorted(clusters_1ms, key=sort_key)[:5]
    for cluster in top_clusters:
        fold = len(cluster)
        t0   = min(t for (t, _) in cluster)
        span = max(t for (t, _) in cluster) - t0
        nids = sorted(n for (_, n) in cluster)
        print(f"  {fold}-fold @ {t0/1e6/60:.4f} min  span={span:.1f} µs  "
              f"Nanos={nids}")

    # ── Plots ────────────────────────────────────────────────────────────────
    print("\n[Plots]")
    plot_detector_rates(kept_events, t_eff_per_nano, rates)
    plot_clock_offsets(sync_array, offset_table)
    plot_calibration_residuals(sync_array, offset_table)
    plot_obs_vs_exp(results_by_window, WINDOWS_MS)
    plot_coincidence_timeline(clusters_by_window, window_ms=1)
    plot_shower_candidate(clusters_by_window, window_ms=1, top_n=3)

    # ── Diagnostic timing validation plots ───────────────────────────────────
    print("\n[Diagnostic Plots — Timing Validation]")
    plot_diag_A_sync_intervals(sync_array)
    plot_diag_B_sync_jitter(sync_array, offset_table)
    plot_diag_C_oos_residuals(sync_array, offset_table)
    plot_diag_D_max_interp_error(sync_array, offset_table)
    plot_diag_E_dt_distribution(kept_events, window_us=1000)
    plot_diag_F_null_test(kept_events, rates, T_eff_s, window_us=1000,
                          n_trials=200)

    print("\nDone.")


if __name__ == "__main__":
    main()
