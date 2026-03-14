# Muon Shower Analysis — Run Log

## Run 1: February 27, 2026 (branch: main)
- Duration: ~60 minutes
- Sync interval: 10 seconds
- Analysis: `muon_analysis.py`
- Jump veto: `muon_analysis_jump_veto.py`
- Key results: 2-fold +4.8σ, 3-fold +4.1σ, 5-fold +4.5σ at 1 ms window
- Notes: Timing floor ~1,500–2,000 µs. Jump artifacts identified in
  multiple Nanos at sync indices k=89, 133, 286 (correlated across all Nanos,
  suggesting Uno timing glitches rather than individual Nano interrupt deferral).
  Results robust to jump veto (max significance change 0.2σ across all
  key window/fold combinations).

## Run 2: [date TBD] (branch: analysis-2hr-1s-sync)
- Duration: 2 hours
- Sync interval: 1 second
- Analysis: `muon_analysis_2hr.py`
- Expected timing floor: ~150–200 µs
- Notes: Shorter sync interval improves calibration resolution by ~10x,
  making the 1 ms coincidence window genuinely meaningful and enabling
  sub-millisecond windows (0.5 ms included). Longer run duration improves
  statistics. Update DATA_FILE in muon_analysis_2hr.py to match the actual
  CSV filename before running.
