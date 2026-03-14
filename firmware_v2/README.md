# Firmware v2 — Muon Shower Detector Array

Muhlenberg College, B. Fadem

This directory contains the rewritten Arduino firmware for the 8-detector
muon shower array. It replaces the original firmware in the repository root
(`FOLLOWER-NANO_synced_code_8-28.ino` and `MASTER-UNO_synced_code_8-28.ino`).

The CSV output format is **identical to v1** — `muon_analysis.py` requires
no changes.

---

## What Changed and Why

The original firmware was analyzed and found to have 14 issues. They are
grouped below by severity.

### Critical (would cause data corruption or incorrect timestamps)

1. **Serial.println() inside sendData() on the Nano** — The I2C `onRequest`
   callback executes with interrupts disabled. Any Serial output inside it
   blocks for 10–20 ms, during which all ISRs are frozen. This caused sync
   pulse timestamps (ST events) to be deferred by one full poll cycle,
   introducing systematic timing errors in the calibration data.
   *Fix: removed all Serial output from sendData().*

2. **Single-slot buffer on the Nano** — v1 stored only the most recent event.
   If a second event arrived before the master polled, the first was silently
   overwritten. Under normal conditions (~1 Hz rate) this was infrequent, but
   during sync pulses (where PT and ST can arrive close together) it caused
   event loss.
   *Fix: replaced with an 8-event circular FIFO queue.*

3. **PT-over-ST priority on the Nano** — v1 would overwrite an ST event in
   the buffer if a PT arrived afterward. Because the master polls Nanos
   sequentially and there is latency between the sync pulse and the poll, this
   caused ST events to be displaced by PT events and delivered one cycle late.
   *Fix: removed priority logic; events are enqueued in true arrival order.*

4. **lastSyncTime updated before pulse fires on the Uno** — v1 sampled `now`
   at the top of `loop()` and used that as `lastSyncTime` after calling
   `sendSyncPulse()`. Because `sendSyncPulse()` takes ~10 ms to execute, the
   interval to the next sync was slightly shorter than intended.
   *Fix: `lastSyncTime = micros()` is called after `sendSyncPulse()` returns.*

5. **Extra delay(10) after sendSyncPulse() in loop() on the Uno** — v1 had a
   `delay(10)` following the `sendSyncPulse()` call, adding 10 ms of dead time
   on top of the 10 ms pulse width. During this 20 ms window, any PT events
   arriving at the Nanos could displace the ST event in the single-slot buffer.
   *Fix: removed the extra delay.*

### Significant (reliability or robustness issues)

6. **No I2C timeout on the Uno** — if a Nano was unresponsive, `Wire.requestFrom()`
   could block indefinitely, halting the entire acquisition loop.
   *Fix: added a 5 ms timeout; unresponsive Nanos are skipped.*

7. **Wire.available() == 6 instead of >= 6** — an exact equality check could
   fail spuriously if extra bytes were present.
   *Fix: changed to >= 6.*

8. **Baud rate 9600 on both devices** — at 9600 baud, transmitting a single
   CSV line (~20 bytes) takes ~20 ms, during which the Serial TX buffer backs
   up and `Serial.println()` blocks. At 115200 baud this drops to ~1.7 ms.
   *Fix: increased to 115200 on both Uno and Nanos.*

9. **No FIFO drain loop on the Uno** — v1 retrieved at most one event per Nano
   per `loop()` iteration. With the v2 FIFO queue on the Nanos, the master
   must loop until it receives the empty sentinel (0xFFFFFFFF) to fully drain
   each Nano's queue.
   *Fix: added drain loop in requestTimestamps().*

10. **INPUT instead of INPUT_PULLUP on interrupt pins** — floating input pins
    can generate spurious interrupts when nothing is connected or when signal
    lines are briefly disconnected.
    *Fix: changed to INPUT_PULLUP on PULSE_PIN and SYNC_PIN.*

### Minor (code quality)

11. **Removed `lastEventType` variable on the Nano** — declared and assigned
    but never read. Dead code removed.

12. **Removed "Sync Pulse Sent" debug print on the Uno** — added no information
    to the data file and introduced unnecessary Serial delay before polling.

13. **lastSyncTime initialized to startTime (not 0) on the Uno** — removes the
    `|| lastSyncTime == 0` special-case logic that was needed to handle the
    uninitialized state.

14. **Named constants replace magic numbers** — pin numbers, addresses, timing
    parameters, and baud rates are now defined as named constants at the top of
    each file for easy configuration.

---

## Configuring Each Nano

Before uploading `FOLLOWER-NANO_v2.ino` to each Nano, set the two constants
at the top of the file:

```cpp
#define SLAVE_ADDRESS  0x08   // I2C address (see table below)
#define SLAVE_ID       1      // Unique ID shown in CSV output (1–8)
```

### I2C Address Map

| Nano ID | SLAVE_ID | SLAVE_ADDRESS |
|---------|----------|---------------|
| Nano 1  | 1        | 0x08          |
| Nano 2  | 2        | 0x09          |
| Nano 3  | 3        | 0x0A          |
| Nano 4  | 4        | 0x0B          |
| Nano 5  | 5        | 0x0C          |
| Nano 6  | 6        | 0x0D          |
| Nano 7  | 7        | 0x0E          |
| Nano 8  | 8        | 0x0F          |

The Uno master (`MASTER-UNO_v2.ino`) requires no per-device configuration.

---

## Compatibility

The CSV output format is identical to v1:

```
MasterStart,<timestamp_us>
MasterSync,<timestamp_us>
<slaveID>,PT,<timestamp_us>
<slaveID>,ST,<timestamp_us>
```

`muon_analysis.py` is fully compatible with v2 data and requires no changes.
