// WORKING TIME STAMP CODE 8/28, just get starting timestamps working

#include <Wire.h>

// === CONFIGURATION ===
// Set these per Nano
#define SLAVE_ADDRESS 0x08 // I2C address
#define SLAVE_ID      1     // Unique ID

volatile bool pulseDetected = false;
volatile unsigned long pulseTimestamp = 0;

volatile bool syncDetected = false;
volatile unsigned long syncTimestamp = 0;

// Last event flag: 1 = pulse, 2 = sync
volatile byte lastEventType = 0;

volatile unsigned long lastPulseTime = 0;
volatile unsigned long lastSyncTime = 0;

void sendData();
void pulseISR();
void syncISR();

void setup() {
  Wire.begin(SLAVE_ADDRESS);
  Wire.onRequest(sendData);

  pinMode(2, INPUT);  // Pulse input
  pinMode(3, INPUT);  // Sync input

  attachInterrupt(digitalPinToInterrupt(2), pulseISR, FALLING);
  attachInterrupt(digitalPinToInterrupt(3), syncISR, RISING);

  Serial.begin(9600);

}

void loop() {
  // Nothing to do here
}

void pulseISR() {
  unsigned long now = micros();
  if (!pulseDetected && (now - lastPulseTime > 5000)) {
    pulseTimestamp = now;
    pulseDetected = true;
    lastEventType = 1;
    lastPulseTime = now;
  }
}

void syncISR() {
  unsigned long now = micros();
  if (!syncDetected && (now - lastSyncTime > 5000)) {
    syncTimestamp = now;
    syncDetected = true;
    lastEventType = 2;
    lastSyncTime = now;
  }
}

void sendData() {
  unsigned long timestampToSend = 0xFFFFFFFF;
  byte typeToSend = 0;

  if (pulseDetected) {
    timestampToSend = pulseTimestamp;
    typeToSend = 1;
    pulseDetected = false;
  } else if (syncDetected) {
    timestampToSend = syncTimestamp;
    typeToSend = 2;
    syncDetected = false;
  }

  Wire.write((byte *)&timestampToSend, sizeof(timestampToSend));
  Wire.write(SLAVE_ID);
  Wire.write(typeToSend);

  // Debug output
  if (typeToSend == 1) {
    Serial.println("SPT"); // Sent Pulse Timestamp
  } else if (typeToSend == 2) {
    Serial.println("SST"); // Sent Sync Timestamp
  } else {
    Serial.println("NNE"); // No New Event
  }
}

