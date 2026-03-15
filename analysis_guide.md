# Muon Shower Coincidence Analysis — Step-by-Step Guide

This document describes the analysis pipeline in plain English, without code.
It is intended as a reference for reproducing the analysis independently,
step by step, and for understanding the purpose and logic of each stage.

The corresponding Python implementation is in `muon_analysis.py`.

---

## The Big Picture

We have 8 muon detectors running simultaneously. Each detector records a
timestamp every time it sees a muon. If muons arrive independently and at
random (a Poisson process), then two detectors firing close together in time
is purely accidental. But if a cosmic ray air shower passes through the array,
it sends a burst of secondary particles to several detectors nearly
simultaneously — a genuine coincidence.

The goal is to count how many multi-detector coincidences we observe, compare
that to how many we would expect from pure chance, and test whether the excess
is statistically significant.

The challenge is that each detector has its own independent clock, and those
clocks drift relative to each other. Before we can compare timestamps across
detectors, we have to calibrate them onto a common time reference.

---

## The Data Format

The data file is a CSV log with four types of rows:

- `MasterStart, <timestamp>` — the Uno's clock time when logging began. Not used in the analysis.
- `MasterSync, <timestamp>` — the Uno's clock time when it broadcast a sync pulse. Not used in the analysis (the Uno's clock is too unreliable).
- `<nano_id>, PT, <timestamp>` — a pulse timestamp from detector nano_id. This is the primary data.
- `<nano_id>, ST, <timestamp>` — a sync timestamp from detector nano_id. This records the moment that Nano's local clock read when it received the sync pulse.

All timestamps are in microseconds, from each device's own `micros()` clock.
The clocks are not synchronized at startup — they each start from zero
independently, and they drift at different rates over time.

---

## Step 1: Parse the File

Read through the file line by line and sort every row into one of three
buckets:

1. **records** — all the PT and ST events, stored as (nano_id, event_type, timestamp)
2. **master_syncs** — the Uno's sync timestamps (kept but not used for calibration)
3. **master_start** — the Uno's start timestamp (kept but not used)

After parsing, you should have roughly 29,000 records for a one-hour run.

---

## Step 2: Build the Sync Table

The sync pulse is broadcast every 10 seconds over a shared wire connected to
all 8 Nanos simultaneously. Every Nano records the moment it received that
pulse, in its own local clock. Because the pulse arrives at all Nanos at the
same physical instant, the differences between their ST timestamps at each sync
event reflect nothing but clock offsets and drift — not any real time
difference.

**The sync table** is an array of shape (n_syncs × 8), where entry [k, i] is
the local timestamp that Nano i recorded for sync event k.

To build it:

1. Collect all ST events, grouped by Nano.
2. Sort each Nano's ST list by timestamp.
3. **Drop the first 2 ST entries for every Nano.** The very first sync pulses
   after startup are anomalous — the Nanos are still initializing, and these
   timestamps are unreliable. This was confirmed by inspection of the raw data.
4. Trim all Nanos to the same number of syncs (use the minimum across all 8),
   so the table is rectangular.

You should end up with approximately 357 sync events per Nano after dropping
the first two.

---

## Step 3: Build the Offset Table and Calibrate

**The offset table** captures, for each sync event k and each Nano i, how far
Nano i's clock is ahead of (or behind) Nano 1's clock at that moment:

```
offset[k, i] = ST[k, i] − ST[k, 0]
```

(Nano 1 is the reference; its offset is zero by definition.)

These offsets change slowly over time as clocks drift. A typical Nano might
drift by −1,000 to −4,000 ms over an hour relative to Nano 1.

**Calibrating a PT event** means converting a raw timestamp from Nano i into
the equivalent time on Nano 1's clock. For a PT event at local time t from
Nano i:

1. Find the two sync events that bracket t — the last sync before t, and the
   first sync after t.
2. Linearly interpolate the offset between those two sync points to estimate
   what the offset was at the moment t was recorded.
3. Subtract that interpolated offset from t to get the corrected time in Nano
   1's reference frame.

This is called piecewise-linear interpolation. It assumes the clock drift is
approximately linear within each 10-second interval, which is a good
approximation for a crystal oscillator.

**Edge cases:** For PT events before the first sync or after the last sync,
extrapolate using the first or last sync pair respectively.

After calibration, all PT events are expressed in Nano 1's time reference and
can be directly compared.

---

## Step 4: The Sync Veto

A practical problem: the sync pulse occasionally bleeds through into the PT
(pulse) channel, causing a Nano to record a false muon hit at the exact moment
of the sync. These fake hits cluster tightly around the ST timestamps.

To remove them, for every PT event from Nano i, check whether its raw
(uncorrected) timestamp falls within 20 ms of any ST timestamp for that same
Nano. If it does, discard it.

The veto window is ±20 ms around each sync timestamp. This is much larger than
any real muon pulse width, so it safely removes bleedthrough without
significantly affecting the physics data.

After vetoing, recalculate the **effective run duration** for each Nano:

```
T_eff = raw_run_duration − (number_of_syncs × 2 × veto_window)
```

This accounts for the dead time introduced by the veto. Use T_eff (not the raw
run duration) in all subsequent rate and accidental calculations.

About 130–140 events are typically removed by the veto, representing ~0.5% of
the data.

---

## Step 5: Per-Detector Rates

For each Nano, count how many PT events survived the veto and divide by that
Nano's T_eff. This gives the rate in Hz.

Typical rates are 0.8–1.3 Hz per detector for a plastic scintillator of this
size.

These rates are the foundation of the accidental coincidence calculation —
they characterize how often each detector fires due to random (unrelated) muon
passages.

---

## Step 6: Find Coincidences

A coincidence is a group of PT events from different detectors that all fall
within a chosen time window of each other. We search for these using a greedy
sliding-window algorithm:

1. Sort all calibrated PT events by corrected timestamp.
2. Walk through the list from earliest to latest. For each event that has not
   yet been assigned to a cluster:
   - Call this event the "seed." Open a window of width W starting at its timestamp.
   - Collect all other events within that window that are from a *different* Nano
     and have not yet been used.
   - If there is more than one Nano represented (including the seed), this is a
     coincidence cluster. Mark all members as used.
3. Record every cluster with multiplicity ≥ 2 (i.e., at least 2 different Nanos).

Key rules:
- Each PT event can belong to at most one cluster (greedy, earliest-first).
- Only one hit per Nano per cluster (the earliest available one within the window).

The window width W is a free parameter. We test W = 1, 2, 5, 10, 20, 50 ms.

The **multiplicity** (or fold) of a cluster is the number of distinct Nanos
represented. A 3-fold coincidence means 3 different detectors all fired within
the window.

---

## Step 7: Expected Accidental Coincidences

If all detectors fire independently at random (Poisson processes), what is the
expected number of k-fold coincidences purely by chance?

For a k-fold coincidence with window half-width τ (in seconds), rates r₁
through r₈ (in Hz), and effective run duration T (in seconds):

```
E_k = T × Σ (over all k-tuples of detectors) [ r_i1 × r_i2 × ... × r_ik × (2τ)^(k−1) ]
```

The sum is over all possible combinations of k detectors chosen from 8.

The logic behind this formula:
- The first detector fires at some rate r_i1.
- For each such firing, the probability that detector i2 also fires within the
  window is approximately r_i2 × 2τ.
- Similarly for each additional detector in the tuple.
- Summing over all k-tuples gives the total expected rate of k-fold accidentals.
- Multiplying by T gives the expected count over the run.

This formula is derived rigorously in `accidental_coincidence_derivation.pdf`.

For a 2-fold coincidence at 1 ms window with typical rates (~1 Hz each), you
expect roughly 200 accidental pairs per hour. For 3-fold, roughly 1 per hour.
For 5-fold, essentially zero.

---

## Step 8: Statistical Significance

For each (window, fold) combination, we now have an observed count O and an
expected accidental count E. The question is: how unlikely is it to observe O
or more events if the true rate is E?

We model the accidental counts as Poisson-distributed with mean E. The
one-sided p-value is:

```
p = P(X ≥ O | Poisson(E)) = 1 − P(X ≤ O−1 | Poisson(E))
```

We convert this to a number of standard deviations (σ) using the inverse of
the normal CDF:

```
σ = Φ⁻¹(1 − p)
```

A result of +4σ means the probability of seeing this many events (or more) by
chance is about 1 in 30,000. The convention for claiming a discovery in
physics is 5σ.

This calculation is derived in `statistical_significance_derivation.pdf`.

---

## Diagnostic Checks (Timing Validation)

Beyond the main analysis, six diagnostic plots are produced to validate that
the calibration is working correctly. These are described briefly here:

**Plot A — Sync interval regularity:** Checks that ST[k+1] − ST[k] is close to
10 seconds for every Nano. Any missing sync pulses or large jitter will appear
here immediately.

**Plot B — Sync detection jitter:** After removing the smooth clock drift trend,
plots the residual scatter in sync detection timing. The residuals contain two
components: (1) genuine ISR jitter (~4–8 µs) reflecting real variability in
when each Nano registers the sync pulse, and (2) isolated spike artifacts caused
by interrupt deferral during I²C callbacks (~10–20 ms), affecting ~1.4% of sync
intervals. The spike artifacts dominate the RMS of this plot. They are not
irreducible — firmware v2 eliminates them by removing Serial output from the
onRequest callback.

**Plot C — Out-of-sample calibration residuals:** The most rigorous test of
calibration accuracy. The calibration is built using only even-numbered sync
events, then used to predict the offset at odd-numbered sync events (which were
not used in building the calibration). The RMS of the prediction errors characterizes the timing accuracy for PT events
falling between sync pulses. However, this RMS is dominated by the ~1.4% of
intervals that contain spike artifacts. Excluding those intervals, the combined
timing floor is ≈ 75 µs, making the 1 ms coincidence window well-justified.
Spike intervals are handled separately by the jump veto.

**Plot D — Maximum interpolation error:** For each 10-second sync interval,
computes the maximum timing error that could result from piecewise-linear
interpolation for an event at the midpoint of that interval. This is
|Δoffset| / 4, where Δoffset is the change in clock offset across the interval.
Note that the spike intervals (~1.4% of all intervals) produce anomalously large
values in this plot. Excluding those, the median maximum interpolation error is
well below 0.5 ms, consistent with the 75 µs combined timing floor.

**Plot E — Δt distribution:** For all pairs of events from different Nanos that
fall within the coincidence window, histograms the time difference |t_i − t_j|.
If real correlated events (showers) are present, there should be an excess of
pairs near Δt = 0 above a flat accidental background. This is the most visually
intuitive evidence for real coincidences.

**Plot F — Null test (time-shift):** Repeats the coincidence search 200 times
after randomly shifting each Nano's event list by an independent random offset.
The shifted data contains only accidentals by construction. Comparing the real
observed counts to the distribution of shuffled counts provides an
analysis-independent confirmation of the statistical significance.

---

## Key Numbers from the February 27, 2026 Run

| Quantity | Value |
|----------|-------|
| Run duration | ~59.9 min |
| Total PT events | 28,935 |
| Events after veto | 28,798 |
| Sync events per Nano | 357 (after dropping first 2) |
| Per-detector rates | 0.81 – 1.25 Hz |
| Clock drift range | −1,264 to +0 ppm (relative to Nano 1) |
| Out-of-sample calibration RMS (all intervals) | ~1.26 ms |
| Out-of-sample calibration RMS (spike intervals excluded) | ~0.075 ms |
| Fraction of spike intervals | ~1.4% |
| 2-fold coincidences at 1 ms | 274 observed, 202 expected, +4.8σ |
| 3-fold coincidences at 1 ms | 7 observed, 0.8 expected, +4.1σ |
| 5-fold coincidence at 1 ms | 1 observed, ~0 expected, +4.5σ |
| Best shower candidate | 5-fold at t = 49.97 min, span = 3.6 µs, Nanos 1,3,5,6,7 |

---

## What the Results Mean and Their Limitations

**What can be claimed:** There is a statistically significant excess of
multi-detector coincidences above the accidental expectation. This excess is
not an artifact of the analysis pipeline (confirmed by the null test). It is
inconsistent with the hypothesis that all detectors fire independently.

**What cannot be claimed from this data alone:**
- That the coincidences are definitively from cosmic ray air showers (as
  opposed to correlated electrical noise).
- That any individual event (including the 5-fold candidate) is certainly a
  shower rather than a noise artifact.
- A discovery-level result (5σ) — the current significance is ~4σ, which is
  strong evidence but below the conventional threshold.

**What would strengthen the result:**
- A second independent data run showing the same excess.
- Hardware investigation to rule out common-mode electrical noise as an
  alternative explanation for multi-detector coincidences.

**Note on timing resolution:** The calibration timing floor of ~75 µs
(excluding spike intervals) means the 1 ms coincidence window is well-supported
by the data quality. The spike intervals (~1.4% of the run) are handled by the
jump veto. The firmware v2 rewrite eliminates the source of the spikes entirely.
