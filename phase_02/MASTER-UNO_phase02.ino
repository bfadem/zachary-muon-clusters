// ============================================================================
//  MASTER UNO — Phase 02
//  Muhlenberg College Muon Detector Array
// ============================================================================
//
//  KEY IMPROVEMENTS OVER PHASE 01
//  ─────────────────────────────────────────────────────────────────────────
//  1. Dedicated PAUSE line (pin 7) asserted before every sync event.
//     All Nanos disable particle-capture while paused, guaranteeing
//     the ST is ready when polled — no more ST/PT delivery-priority jumps.
//  2. Sync interval reduced to 1 second (from 10 s), reducing the
//     inter-sync calibration uncertainty by ~10×.
//  3. I²C at 400 kHz Fast mode — per-Nano poll takes ~150 µs.
//  4. Serial at 115200 baud — eliminates logging bottleneck.
//  5. 8-byte I²C receive packets: 4 bytes ticks_hi + 2 bytes ticks_lo
//     + 1 byte Nano ID + 1 byte event type.
//  6. ST collection separated into its own doSync() function, which
//     only runs during the pause window.  Normal requestTimestamps()
//     handles PT collection.
//
//  DATA FORMAT (CSV)
//  ─────────────────────────────────────────────────────────────────────────
//  # Phase02 timestamp unit: 500 ns per tick (Timer1 prescaler 8)
//  MasterStart,<micros>
//  MasterSync,<micros>             ← Uno's micros() at end of sync window
//  <nano_id>,ST,<ticks_hi>,<ticks_lo>   ← guaranteed ST from paused Nano
//  <nano_id>,PT,<ticks_hi>,<ticks_lo>   ← normal particle pulse
//
//  Timestamp reconstruction in analysis:
//    ticks    = (uint64_t(ticks_hi) << 16) | ticks_lo
//    time_us  = ticks * 0.5   (1 tick = 500 ns)
//
//  PIN ASSIGNMENTS (Arduino Uno)
//  ─────────────────────────────────────────────────────────────────────────
//  Pin 6  (SYNC_PIN)   → sync pulse wire  → all Nanos pin 3
//  Pin 7  (PAUSE_PIN)  → pause/resume wire → all Nanos pin 2
//  Pin A4 (SDA)        ↔ I²C SDA bus
//  Pin A5 (SCL)        ↔ I²C SCL bus
//  USB                 → host PC (serial data logging)
//
// ============================================================================

#include <Wire.h>

// ── Run configuration ─────────────────────────────────────────────────────────
const unsigned long RUN_DURATION_US  = 60UL * 60UL * 1000000UL;  // 60 minutes
const unsigned long SYNC_INTERVAL_US = 1000000UL;                  // 1 second

// ── Hardware pins ─────────────────────────────────────────────────────────────
const int SYNC_PIN  = 6;   // Sync pulse output  → all Nanos pin 3
const int PAUSE_PIN = 7;   // Pause/resume output → all Nanos pin 2

// ── I²C configuration ────────────────────────────────────────────────────────
const uint8_t NANO_ADDRESSES[]  = {0x08, 0x09, 0x0A, 0x0B,
                                    0x0C, 0x0D, 0x0E, 0x0F};
const int     NUM_NANOS          = 8;
const int     BYTES_PER_PACKET   = 8;   // 4 ticks_hi + 2 ticks_lo + 1 ID + 1 type

// ── Timing state ─────────────────────────────────────────────────────────────
unsigned long startTime    = 0;
unsigned long lastSyncTime = 0;

// ─────────────────────────────────────────────────────────────────────────────
void setup() {
    Wire.begin();
    Wire.setClock(400000);       // Fast mode (400 kHz)

    Serial.begin(115200);

    pinMode(SYNC_PIN,  OUTPUT);  digitalWrite(SYNC_PIN,  LOW);
    pinMode(PAUSE_PIN, OUTPUT);  digitalWrite(PAUSE_PIN, LOW);

    startTime = micros();
    Serial.print("MasterStart,");
    Serial.println(startTime);
    Serial.println("# Phase02 timestamp unit: 500 ns per tick (Timer1 prescaler 8)");
}

// ─────────────────────────────────────────────────────────────────────────────
void loop() {
    unsigned long now = micros();

    if (now - startTime >= RUN_DURATION_US) {
        Serial.println("Done");
        while (true) {}
    }

    // Fire a sync event at the configured interval
    if (lastSyncTime == 0 || (now - lastSyncTime) >= SYNC_INTERVAL_US) {
        doSync();
        lastSyncTime = micros();   // record after sync completes
    }

    // Normal particle-pulse polling
    requestTimestamps();
}

// ── doSync: pause → pulse → collect STs → resume ─────────────────────────────
//
// The sequence guarantees that when the Uno polls each Nano, no PT event
// can have arrived to preempt the ST in the Nano's output buffer.
//
void doSync() {

    // ── Step 1: Assert PAUSE line ─────────────────────────────────────────
    // Nanos respond in their INT0 ISR (any-edge), which runs within a few
    // microseconds of the rising edge on their pin 2.
    digitalWrite(PAUSE_PIN, HIGH);
    delayMicroseconds(300);   // Allow all eight Nanos' INT0 ISRs to complete

    // ── Step 2: Send sync pulse ───────────────────────────────────────────
    // The rising edge on each Nano's pin 3 (INT1) triggers the syncISR,
    // which reads TCNT1 and packs the ST packet.  The 5 ms pulse width
    // ensures the edge is registered even under marginal signal conditions.
    digitalWrite(SYNC_PIN, HIGH);
    delay(5);                  // 5 ms pulse width
    digitalWrite(SYNC_PIN, LOW);

    // ── Step 3: Record Uno's own time and collect all STs ────────────────
    Serial.print("MasterSync,");
    Serial.println(micros());

    for (int i = 0; i < NUM_NANOS; i++) {
        Wire.requestFrom(NANO_ADDRESSES[i], (uint8_t)BYTES_PER_PACKET);

        if (Wire.available() == BYTES_PER_PACKET) {
            uint8_t buf[BYTES_PER_PACKET];
            for (int j = 0; j < BYTES_PER_PACKET; j++) {
                buf[j] = Wire.read();
            }

            uint8_t slaveID  = buf[6];
            uint8_t dataType = buf[7];

            if (dataType == 2) {
                // Expected: ST packet — print ticks_hi and ticks_lo
                uint32_t hi;
                uint16_t lo;
                memcpy(&hi, buf,     4);
                memcpy(&lo, buf + 4, 2);

                Serial.print(slaveID);
                Serial.print(",ST,");
                Serial.print(hi);
                Serial.print(",");
                Serial.println(lo);

            } else if (dataType == 1) {
                // A PT slipped through before pause took effect — log it
                uint32_t hi;
                uint16_t lo;
                memcpy(&hi, buf,     4);
                memcpy(&lo, buf + 4, 2);
                Serial.print(slaveID);
                Serial.print(",PT,");
                Serial.print(hi);
                Serial.print(",");
                Serial.println(lo);
                // This Nano's ST will be missing for this sync event;
                // the analysis should flag and skip this sync pair.

            } else {
                // type == 0: no event returned — ST not captured
                Serial.print("WARN_NO_ST,");
                Serial.println(slaveID);
            }
        }
    }

    // ── Step 4: Deassert PAUSE line — resume normal data-taking ──────────
    digitalWrite(PAUSE_PIN, LOW);
}

// ── requestTimestamps: normal PT polling ──────────────────────────────────────
void requestTimestamps() {
    for (int i = 0; i < NUM_NANOS; i++) {
        Wire.requestFrom(NANO_ADDRESSES[i], (uint8_t)BYTES_PER_PACKET);

        if (Wire.available() == BYTES_PER_PACKET) {
            uint8_t buf[BYTES_PER_PACKET];
            for (int j = 0; j < BYTES_PER_PACKET; j++) {
                buf[j] = Wire.read();
            }

            uint8_t slaveID  = buf[6];
            uint8_t dataType = buf[7];

            if (dataType == 1) {
                uint32_t hi;
                uint16_t lo;
                memcpy(&hi, buf,     4);
                memcpy(&lo, buf + 4, 2);

                Serial.print(slaveID);
                Serial.print(",PT,");
                Serial.print(hi);
                Serial.print(",");
                Serial.println(lo);
            }
            // type == 0: no pending event — nothing to log
            // type == 2: ST during normal polling (should not occur; ignored)
        }
    }
}
