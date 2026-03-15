// WORKING CODE WITH TIME STAMP 8/28

#include <Wire.h>

// ===== CONFIGURATION =====
const unsigned long runDuration = 10UL * 60UL * 1000000UL;  // 10 minutes
const unsigned long syncInterval = 10UL * 1000000UL;        // 10 seconds
const int SYNC_PIN = 6;

const byte nanoAddresses[] = {
  0x08, 0x09, 0x0A, 0x0B,
  0x0C, 0x0D, 0x0E, 0x0F
};
const int numNanos = sizeof(nanoAddresses) / sizeof(nanoAddresses[0]);

unsigned long startTime = 0;
unsigned long lastSyncTime = 0;

void setup() {
  Wire.begin();         
  Serial.begin(9600);   
  pinMode(SYNC_PIN, OUTPUT);
  digitalWrite(SYNC_PIN, LOW);
  startTime = micros(); 

  Serial.print("MST: "); // Master Start Time
  Serial.println(startTime);
}

void loop() {
  unsigned long now = micros();

  if (now - startTime >= runDuration) {
    Serial.println("Done");
    while (true);
  }

  if (now - lastSyncTime >= syncInterval || lastSyncTime == 0) {
    sendSyncPulse();
    delay(10);
    lastSyncTime = now;
  }

  // Always poll the Nanos for data
  requestTimestamps();
}

void sendSyncPulse() {
  digitalWrite(SYNC_PIN, HIGH);
  delay(10);  // 10ms sync pulse
  digitalWrite(SYNC_PIN, LOW);
  Serial.println("Sync Pulse Sent");
}

void requestTimestamps() {
  for (int i = 0; i < numNanos; i++) {
    byte addr = nanoAddresses[i];
    Wire.requestFrom(addr, (uint8_t)6);  // 4 bytes timestamp, 1 byte ID, 1 byte type

    if (Wire.available() == 6) {
      unsigned long timestamp = 0;
      for (int j = 0; j < 4; j++) {
        ((byte *)&timestamp)[j] = Wire.read();
      }
      byte slaveID = Wire.read();
      byte dataType = Wire.read();  // 1 = normal pulse, 2 = sync pulse

      if (timestamp != 0xFFFFFFFF) {
        Serial.print("Nano ");
        Serial.print(slaveID);
        Serial.print(" - ");
        if (dataType == 1) {
          Serial.print("PT: ");
        } else if (dataType == 2) {
          Serial.print("ST: ");
        } else {
          Serial.print("UT: ");
        }
        Serial.println(timestamp);
      }
    }
  }
}
