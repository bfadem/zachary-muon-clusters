// ============================================================
// FOLLOWER-NANO_v2.ino
// Muon Shower Detector Array — Nano Follower Firmware v2
// Muhlenberg College, B. Fadem
//
// Change SLAVE_ADDRESS and SLAVE_ID for each Nano before uploading.
// ============================================================

// ── Per-Nano configuration ───────────────────────────────────
#define SLAVE_ADDRESS  0x08   // I2C address: 0x08-0x0F for Nanos 1-8
#define SLAVE_ID       1      // Unique ID: 1-8

// ── Pin assignments ──────────────────────────────────────────
const int PULSE_PIN = 2;   // Particle pulse input (hardware interrupt 0)
const int SYNC_PIN  = 3;   // Sync pulse input (hardware interrupt 1)

// ── FIFO event queue ─────────────────────────────────────────
const byte QUEUE_SIZE = 8;

struct Event {
  unsigned long timestamp;
  byte eventType;   // 1 = PT (particle), 2 = ST (sync)
};

volatile Event eventQueue[QUEUE_SIZE];
volatile byte qHead  = 0;   // index of next event to send
volatile byte qTail  = 0;   // index of next empty slot
volatile byte qCount = 0;   // number of events currently in queue

// ── Debounce state ───────────────────────────────────────────
volatile unsigned long lastPulseTime = 0;
volatile unsigned long lastSyncTime  = 0;

// ── Enqueue (called from ISRs only) ─────────────────────────
void enqueue(unsigned long ts, byte etype) {
  if (qCount >= QUEUE_SIZE) return;  // queue full, drop event
  eventQueue[qTail].timestamp = ts;
  eventQueue[qTail].eventType = etype;
  qTail = (qTail + 1) % QUEUE_SIZE;
  qCount++;
}

// ── ISRs ─────────────────────────────────────────────────────
void pulseISR() {
  unsigned long now = micros();
  if (now - lastPulseTime > 5000) {
    enqueue(now, 1);
    lastPulseTime = now;
  }
}

void syncISR() {
  unsigned long now = micros();
  if (now - lastSyncTime > 5000) {
    enqueue(now, 2);
    lastSyncTime = now;
  }
}

// ── I2C onRequest callback ───────────────────────────────────
// CRITICAL: no Serial output here — executes with interrupts disabled.
// Sends one event from the front of the FIFO.
// Sends sentinel 0xFFFFFFFF / type 0 when queue is empty.
void sendData() {
  unsigned long tsToSend   = 0xFFFFFFFF;
  byte          typeToSend = 0;

  if (qCount > 0) {
    tsToSend   = eventQueue[qHead].timestamp;
    typeToSend = eventQueue[qHead].eventType;
    qHead = (qHead + 1) % QUEUE_SIZE;
    qCount--;
  }

  Wire.write((byte*)&tsToSend, 4);
  Wire.write(SLAVE_ID);
  Wire.write(typeToSend);
}

// ── setup ────────────────────────────────────────────────────
void setup() {
  Wire.begin(SLAVE_ADDRESS);
  Wire.onRequest(sendData);

  pinMode(PULSE_PIN, INPUT_PULLUP);
  pinMode(SYNC_PIN,  INPUT_PULLUP);

  attachInterrupt(digitalPinToInterrupt(PULSE_PIN), pulseISR, FALLING);
  attachInterrupt(digitalPinToInterrupt(SYNC_PIN),  syncISR,  RISING);

  // Serial used only for startup confirmation — never again after setup
  Serial.begin(115200);
  Serial.print("Nano ");
  Serial.print(SLAVE_ID);
  Serial.println(" ready.");
}

// ── loop ─────────────────────────────────────────────────────
void loop() {
  // All work is done in ISRs and sendData().
  // Nothing to do here.
}
