# Claude Code Prompt: Muon Shower Coincidence Analysis

## Context

You are analyzing cosmic ray muon detector data from an 8-detector plastic scintillator array at Muhlenberg College. Each detector has an Arduino Nano that timestamps muon pulses using `micros()`. A master Arduino Uno polls the Nanos via I2C and broadcasts sync pulses every 10 seconds to all Nanos simultaneously via a shared wire. The result is a CSV log file with pulse timestamps (PT) and sync timestamps (ST) per detector, plus MasterSync times from the Uno.

## Data Format

The CSV log file has this structure:
```
# Logging started at <datetime>
MasterStart,<timestamp_us>
MasterSync,<timestamp_us>
<nano_id>,PT,<timestamp_us>
<nano_id>,ST,<timestamp_us>
...
# Done
```

- `MasterStart` / `MasterSync`: Uno clock times in µs
- `<nano_id>,PT,...`: pulse timestamp from Nano `nano_id` (1–8) in that Nano's local `micros()` clock
- `<nano_id>,ST,...`: sync timestamp — the time Nano `nano_id` recorded when it received the sync pulse
- All 8 Nanos receive each sync pulse simultaneously (shared wire), so `ST` differences between Nanos reflect only clock offsets, not propagation delay

## Key Findings from Initial Analysis (February 27, 2026 run)

- 8 Nanos active, IDs 1–8
- Run duration: ~59.85 minutes, 28,935 total PT events
- Per-detector rates: 0.81–1.25 Hz (healthy, stable)
- 360 sync events per Nano; first 2 ST entries per Nano are anomalous (startup artifact) and must be dropped
- Clock drift rates relative to Nano 1: −1,266 to +456 ppm over the hour
- Best calibration: pairwise interpolation (Nano_i − Nano_1 offset interpolated between consecutive sync pairs)
- Residual timing uncertainty after calibration: σ ≈ 1,500–2,000 µs per detector pair
- Sync veto needed: ±20 ms around raw ST timestamps per Nano removes sync-pulse bleedthrough into PT channel (129 events, 0.4%)
- Primary shower candidate: 5-fold coincidence at t = 49.97 min, Nanos 1,3,5,6,7, span = 3.6 µs, significance 4.5σ at 1 ms window
- 3-fold excess: 4.3–4.6σ at 1–2 ms windows vs. Poisson accidental expectation

## Analysis Pipeline to Implement

### Step 1: Parse

```python
records, master_syncs, master_start = [], [], None
with open("data_file.csv") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"): continue
        parts = line.split(",")
        if parts[0] == "MasterStart": master_start = int(parts[1])
        elif parts[0] == "MasterSync": master_syncs.append(int(parts[1]))
        elif len(parts) == 3: records.append((int(parts[0]), parts[1], int(parts[2])))
```

### Step 2: Build Sync Table

```python
# For each Nano, sort ST events and drop the first 2 (anomalous startup entries)
# Result: sync_array shape (n_syncs, 8), where sync_array[k, i] = Nano (i+1)'s
# local timestamp at sync event k
# Trim all Nanos to same length (min across Nanos)
```

### Step 3: Interpolation Calibration (Nano-to-Nano, no master needed)

```python
# Reference: Nano 1 (index 0)
# offset_table[k, i] = sync_array[k, i] - sync_array[k, 0]
# For PT event from Nano i at local time t_nano:
#   Find k such that sync_array[k, i] <= t_nano <= sync_array[k+1, i]
#   Interpolate: offset = offset_table[k,i] + frac * (offset_table[k+1,i] - offset_table[k,i])
#   t_corr = t_nano - offset
# Edge cases: extrapolate using first/last two sync pairs
```

### Step 4: Per-Nano Sync Veto

```python
# For each Nano, build array of raw ST timestamps
# Veto any PT event where |raw_PT_timestamp - any_raw_ST_timestamp| < VETO_US (= 20000)
# Recalculate effective run duration: subtract (n_ST * 2 * VETO_US / 1e6) per Nano
```

### Step 5: Coincidence Counting

```python
# Sort calibrated PT events by t_corr
# Greedy sliding-window: for each unused pulse i, find all pulses within window_us
# from t_corr[i] that are from distinct Nanos and not yet used
# Mark all members of cluster as used; record if multiplicity >= 2
# Key: each pulse used at most once (greedy earliest-first)
```

### Step 6: Expected Accidentals

```python
# For k-fold coincidence with window tau (seconds):
# E_k = sum over all k-tuples of detectors: prod(r_i) * (2*tau)^(k-1) * T_eff
# where r_i = rate of detector i, T_eff = effective run duration
from itertools import combinations
import numpy as np
def expected_accidentals(window_us, fold, rates, run_duration_s):
    tau = window_us / 1e6
    return sum(np.prod([rates[n] for n in combo]) * (2*tau)**(fold-1) * run_duration_s
               for combo in combinations(range(1,9), fold))
```

### Step 7: Statistical Significance

```python
from scipy.stats import poisson, norm
# p_value = 1 - poisson.cdf(observed - 1, expected)  # P(X >= observed)
# sigma = norm.ppf(1 - p_value)
```

## Suggested Output Plots

1. **Per-detector pulse rates** (bar chart + rate vs. time in 60s bins)
2. **Pairwise clock offsets** (Nano_i − Nano_1 vs. sync index, with linear fit)
3. **Calibration residuals** (offset change per 10s sync interval — sets timing resolution floor)
4. **Observed vs. expected coincidences** (log-log, per fold, vs. window size)
5. **Significance vs. window size** (for ≥3, ≥4, ≥5-fold)
6. **Coincidence timeline** (scatter plot of all coincidence events vs. time, colored by multiplicity)
7. **High-multiplicity event detail** (timestamp ladder for the primary shower candidate)

## Coincidence Windows to Explore
`[1, 2, 5, 10, 20, 50]` ms

## Known Issues / Watch Points

- **First 2 ST entries per Nano**: anomalous (startup timing), must be dropped before building sync table
- **8-fold artifact at t ≈ 9.58 s**: sync pulse bleedthrough — removed by ST veto
- **Clock drift is nonlinear over 1 hour**: linear fit is a good approximation, but piecewise or spline interpolation between consecutive sync pairs is more accurate (already implemented in Step 3 above)
- **Significance falls at large windows**: expected — accidentals dominate; peak significance at tight windows is the physically meaningful result
- **4-fold event at t = 24.32 min (Nanos 4,5,7,8, span 9.9 ms)**: not an artifact but ambiguous — span is wide and ~2 accidentals expected at 10 ms window

## Environment

Recommended: Python 3, conda environment with:
```
conda create -n muon-analysis python=3.11 numpy pandas scipy matplotlib jupyter
conda activate muon-analysis
```

## Future Extensions

1. **Nonlinear calibration**: Replace linear fits with spline interpolation between consecutive sync pairs for drift correction
2. **Shower direction**: If detector (x,y,z) positions are known, fit arrival time differences to a plane wave to reconstruct shower angle
3. **Energy proxy**: High-multiplicity events (6, 7, 8-fold) suggest higher-energy primaries; build a multiplicity spectrum
4. **Rate stability**: Monitor per-detector rates in short time bins to look for correlated rate increases during shower episodes
5. **Hardware improvement**: Trigger master sync timestamp at interrupt level (not in software loop) to reduce sync jitter and improve calibration accuracy
