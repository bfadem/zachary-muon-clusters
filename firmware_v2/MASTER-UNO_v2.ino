// ============================================================
// MASTER-UNO_v2.ino
// Muon Shower Detector Array — Uno Master Firmware v2
// Muhlenberg College, B. Fadem
// ============================================================

#include <Wire.h>

// ── Configuration constants ──────────────────────────────────
const unsigned long RUN_DURATION_US  = 60UL * 60UL * 1000000UL; // 60 min
const unsigned long SYNC_INTERVAL_US = 10UL * 1000000UL;         // 10 s
const unsigned long SYNC_PULSE_MS    = 10UL;                      // 10 ms pulse width
const int           SYNC_PIN         = 6;
const int           NUM_NANOS        = 8;
const long          SERIAL_BAUD      = 115200;
const unsigned long I2C_TIMEOUT_US   = 5000;  // 5 ms timeout per Nano poll

const byte NANO_ADDRESSES[NUM_NANOS] = {
  0x08, 0x09, 0x0A, 0x0B,
  0x0C, 0x0D, 0x0E, 0x0F
};

// ── Global state ─────────────────────────────────────────────
unsigned long startTime    = 0;
unsigned long lastSyncTime = 0;

// ── setup ────────────────────────────────────────────────────
void setup() {
  Wire.begin();
  Serial.begin(SERIAL_BAUD);

  pinMode(SYNC_PIN, OUTPUT);
  digitalWrite(SYNC_PIN, LOW);

  startTime    = micros();
  lastSyncTime = startTime;   // Initialized to startTime, not 0 —
                               // prevents immediate sync pulse at startup

  Serial.print("MasterStart,");
  Serial.println(startTime);
}

// ── loop ─────────────────────────────────────────────────────
void loop() {
  unsigned long now = micros();

  // End of run
  if (now - startTime >= RUN_DURATION_US) {
    Serial.println("Done");
    while (true);
  }

  // Sync pulse
  if (now - lastSyncTime >= SYNC_INTERVAL_US) {
    sendSyncPulse();
    lastSyncTime = micros();  // CRITICAL: recorded AFTER pulse fires,
                               // so the next interval is measured from
                               // the actual post-pulse time
  }

  // Poll all Nanos
  requestTimestamps();
}

// ── sendSyncPulse ────────────────────────────────────────────
void sendSyncPulse() {
  // Optional: log MasterSync timestamp before raising the pin.
  // Remove these three lines if you want to minimize the delay
  // between the Serial.print and the rising edge.
  unsigned long syncTime = micros();
  Serial.print("MasterSync,");
  Serial.println(syncTime);

  digitalWrite(SYNC_PIN, HIGH);
  delay(SYNC_PULSE_MS);
  digitalWrite(SYNC_PIN, LOW);
  // Removed: Serial.println("Sync Pulse Sent") — adds no information
  // to data file and delays start of polling
}

// ── requestTimestamps ────────────────────────────────────────
// Polls each Nano in turn, draining its FIFO completely before
// moving to the next Nano.
void requestTimestamps() {
  for (int i = 0; i < NUM_NANOS; i++) {
    byte addr = NANO_ADDRESSES[i];

    // Drain this Nano's FIFO completely before moving on
    while (true) {
      Wire.requestFrom(addr, (uint8_t)6);

      // Wait for 6 bytes with timeout
      unsigned long t0 = micros();
      while (Wire.available() < 6) {
        if (micros() - t0 > I2C_TIMEOUT_US) break;
      }

      if (Wire.available() < 6) break;  // Nano unresponsive, skip

      unsigned long timestamp = 0;
      for (int j = 0; j < 4; j++) {
        ((byte*)&timestamp)[j] = Wire.read();
      }
      byte slaveID  = Wire.read();
      byte dataType = Wire.read();

      // Empty sentinel — Nano has no more events
      if (timestamp == 0xFFFFFFFF) break;

      // Log to Serial in CSV format
      Serial.print(slaveID);
      Serial.print(",");
      if (dataType == 1) {
        Serial.print("PT,");
      } else if (dataType == 2) {
        Serial.print("ST,");
      } else {
        Serial.print("UT,");
      }
      Serial.println(timestamp);
    }
  }
}
