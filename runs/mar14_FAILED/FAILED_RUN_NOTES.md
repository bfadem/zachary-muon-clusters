# March 14, 2026 Run — FAILED

**This run is NOT usable. All outputs from `muon_analysis_2hr.py` on this
data should be ignored.**

---

## Root Cause

The Phase 02 follower firmware (`firmware/phase_02/FOLLOWER-NANO_phase02.ino`)
configures the sync pin with `INPUT` instead of `INPUT_PULLUP`:

```cpp
// BUG — current code:
pinMode(SYNC_PIN, INPUT);

// FIX — one-line change:
pinMode(SYNC_PIN, INPUT_PULLUP);
```

Without an internal pull-up resistor, the sync line floats slightly between
pulses. At 10-second sync intervals this is harmless — the line has plenty
of time to settle and the 5 ms debounce window comfortably rejects any
bounce. But at **1-second sync intervals**, the bounce from the rising edge
of one sync pulse can retrigger the ISR within the 5 ms debounce window of
the *next* pulse, causing the Nano to record two ST events (one genuine, one
bounce) for every sync.

---

## Symptom

Sync timestamps show alternating short and long intervals instead of a
uniform ~1,000,000 µs:

- Nano 1: alternating ~480,000 µs and ~580,000 µs
- Nano 2: alternating ~113,000 µs and ~887,000 µs

The ratio between short and long intervals varies by Nano depending on
exactly when the bounce occurs relative to the true sync edge.

---

## Effect on Analysis

The jump-interval veto in `muon_analysis_2hr.py` correctly identifies these
pairs of short/long intervals as corrupted (the d_delta anomaly is enormous
— often 200–400 ms, far above the 5 ms threshold). Approximately **68% of
all sync intervals are flagged**, and the PT events falling in those
intervals are discarded. With only ~32% of the data usable:

- Per-detector rates are underestimated
- Expected accidental counts are wrong
- Final significance: ~2σ at 2-fold, 1 ms (not significant)

---

## Required Fix Before Next Run

1. Open `firmware/phase_02/FOLLOWER-NANO_phase02.ino`
2. Change `pinMode(SYNC_PIN, INPUT)` to `pinMode(SYNC_PIN, INPUT_PULLUP)`
3. Reflash all 8 Nano followers
4. Verify sync intervals are uniform (~1,000,000 µs ± a few hundred µs)
   before starting a long run — check `diag_A_sync_intervals_2hr.png`
   after a short test run

No other firmware changes are needed. The Phase 02 master firmware
(`firmware/phase_02/MASTER-UNO_phase02.ino`) is not affected.

---

## What the 2-Hour Run Was Supposed to Achieve

- 1-second sync intervals → timing resolution floor ~150–200 µs (vs ~1.3 ms
  at 10-second intervals)
- 2-hour run → ~4× more statistics than Feb 27
- Sub-millisecond coincidence windows (0.5 ms) now physically meaningful
- Expected: push 2-fold significance from ~4.8σ toward or past 5σ

Once the firmware is fixed, `muon_analysis_2hr.py` is ready to use. Update
`DATA_FILE` to point at the new CSV and run normally.
