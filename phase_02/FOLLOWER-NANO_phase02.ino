// ============================================================================
//  FOLLOWER NANO — Phase 02
//  Muhlenberg College Muon Detector Array
// ============================================================================
//
//  KEY IMPROVEMENTS OVER PHASE 01
//  ─────────────────────────────────────────────────────────────────────────
//  1. Particle pulse timestamps via Timer1 Input Capture Unit (ICP1, pin 8)
//     — hardware latches the timer value with zero software latency,
//       eliminating the 3–4 µs random ISR-entry jitter of the old approach.
//  2. Timer1 prescaler 8  →  1 tick = 500 ns  (8× finer than micros()).
//  3. Sync timestamps read via TCNT1 directly in the INT1 ISR — same
//     500 ns resolution, and the Nano-to-Nano difference cancels the small
//     fixed ISR latency.
//  4. Dedicated pause line (pin 2) driven HIGH by the Uno before every sync
//     event, guaranteeing that no PT event can preempt ST delivery.
//  5. I²C bus at 400 kHz (Fast mode) — reduces per-Nano poll time to ~150 µs.
//  6. Serial at 115200 baud — eliminates logging bottleneck.
//  7. 8-byte I²C payload: 4 bytes ticks_hi + 2 bytes ticks_lo + 1 ID + 1 type.
//
//  TIMESTAMP FORMAT
//  ─────────────────────────────────────────────────────────────────────────
//  A 48-bit tick count is transmitted as two fields:
//    ticks_hi  (uint32_t) = Timer1 overflow count
//    ticks_lo  (uint16_t) = Timer1 TCNT1 value at capture moment
//  Reconstruction: ticks = ((uint64_t)ticks_hi << 16) | ticks_lo
//  Conversion:     time_us = ticks × 0.5   (since 1 tick = 500 ns)
//  Range:          2^48 × 500 ns ≈ 1.4 × 10^8 seconds — no overflow concern.
//
//  PIN ASSIGNMENTS (ATmega328P / Arduino Nano)
//  ─────────────────────────────────────────────────────────────────────────
//  Pin 8  (PB0, ICP1)  ←  Particle pulse from discriminator (FALLING edge)
//  Pin 3  (PD3, INT1)  ←  Sync pulse from Uno pin 6         (RISING  edge)
//  Pin 2  (PD2, INT0)  ←  Pause/resume from Uno pin 7       (ANY     edge)
//  Pin A4 (PC4, SDA)   ↔  I²C SDA bus
//  Pin A5 (PC5, SCL)   ↔  I²C SCL bus
//
//  PER-NANO CONFIGURATION — change these two lines on each Nano
//  ─────────────────────────────────────────────────────────────────────────
#define SLAVE_ADDRESS  0x08   // I²C address: 0x08→Nano1, 0x09→Nano2, … 0x0F→Nano8
#define SLAVE_ID       1      // Numeric ID:  1 through 8
//  ─────────────────────────────────────────────────────────────────────────

#include <Wire.h>

// ── Timer1 state ─────────────────────────────────────────────────────────────
// timer1_ovf counts how many times TCNT1 has rolled over (every 32.768 ms).
// Full 48-bit timestamp = (timer1_ovf << 16) | TCNT1_at_capture.
volatile uint32_t timer1_ovf = 0;

// ── Pause flag ───────────────────────────────────────────────────────────────
// Set HIGH by Uno before sync pulse, LOW after.  While paused, ICP captures
// are discarded so that no PT event can preempt the ST in the send buffer.
volatile bool paused = false;

// ── Pre-assembled 8-byte I²C send packets ────────────────────────────────────
// Each packet: [ticks_hi 4B][ticks_lo 2B][SLAVE_ID 1B][type 1B]
// type: 0 = no event, 1 = particle pulse (PT), 2 = sync timestamp (ST)
volatile uint8_t pt_packet[8] = {0xFF,0xFF,0xFF,0xFF, 0xFF,0xFF, SLAVE_ID, 0};
volatile bool    pt_ready     = false;

volatile uint8_t st_packet[8] = {0xFF,0xFF,0xFF,0xFF, 0xFF,0xFF, SLAVE_ID, 0};
volatile bool    st_ready     = false;

// Forward declaration
void i2c_send_callback();

// ── ISR: Timer1 overflow ─────────────────────────────────────────────────────
ISR(TIMER1_OVF_vect) {
    timer1_ovf++;
}

// ── ISR: ICP1 — hardware-captured particle pulse (pin 8, falling edge) ───────
//
// ICR1 is latched by hardware at the exact moment the edge arrives, with no
// software latency.  We read it immediately (first instruction), then read
// timer1_ovf, and correct for the rare case where an overflow occurred between
// the hardware latch and our software read.
//
ISR(TIMER1_CAPT_vect) {
    // Discard if we are in the sync window or a PT is already pending
    if (paused || pt_ready) return;

    uint16_t cap = ICR1;          // ← read latched value FIRST
    uint32_t ovf = timer1_ovf;   // ← then read overflow count

    // Overflow boundary correction:
    // If TOV1 is set (overflow pending, OVF ISR not yet serviced) AND the
    // captured value is in the lower half of the range (just rolled past 0),
    // the overflow happened BEFORE the capture — count it.
    if ((TIFR1 & _BV(TOV1)) && !(cap & 0x8000)) {
        ovf++;
    }

    // Pack entire packet atomically inside this ISR
    pt_packet[0] = (uint8_t)(ovf      );
    pt_packet[1] = (uint8_t)(ovf >>  8);
    pt_packet[2] = (uint8_t)(ovf >> 16);
    pt_packet[3] = (uint8_t)(ovf >> 24);
    pt_packet[4] = (uint8_t)(cap      );
    pt_packet[5] = (uint8_t)(cap >>  8);
    pt_packet[6] = SLAVE_ID;
    pt_packet[7] = 1;   // event type: PT
    pt_ready     = true;
}

// ── ISR: INT1 — sync pulse (pin 3, rising edge) ───────────────────────────────
//
// We read TCNT1 directly (not ICR1) because ICP1 is already assigned to the
// particle pulse.  The INT1 ISR adds ~2–4 instruction-cycle latency (~125–250 ns)
// relative to the rising edge, but this fixed offset cancels in Nano-to-Nano
// differences and does not affect calibration quality.
//
ISR(INT1_vect) {
    if (st_ready) return;   // discard if previous ST not yet transmitted

    uint16_t t   = TCNT1;        // read Timer1 directly
    uint32_t ovf = timer1_ovf;
    if ((TIFR1 & _BV(TOV1)) && !(t & 0x8000)) {
        ovf++;
    }

    st_packet[0] = (uint8_t)(ovf      );
    st_packet[1] = (uint8_t)(ovf >>  8);
    st_packet[2] = (uint8_t)(ovf >> 16);
    st_packet[3] = (uint8_t)(ovf >> 24);
    st_packet[4] = (uint8_t)(t       );
    st_packet[5] = (uint8_t)(t  >>  8);
    st_packet[6] = SLAVE_ID;
    st_packet[7] = 2;   // event type: ST
    st_ready     = true;
}

// ── ISR: INT0 — pause/resume line (pin 2, any edge) ──────────────────────────
ISR(INT0_vect) {
    paused = (PIND & _BV(PD2)) ? true : false;
}

// ── I²C request callback (invoked from Wire's TWI ISR) ───────────────────────
void i2c_send_callback() {
    if (pt_ready) {
        Wire.write((const uint8_t*)pt_packet, 8);
        pt_ready = false;
    } else if (st_ready) {
        Wire.write((const uint8_t*)st_packet, 8);
        st_ready = false;
    } else {
        // No pending event — send all-0xFF "no event" packet
        const uint8_t empty[8] = {0xFF,0xFF,0xFF,0xFF, 0xFF,0xFF, SLAVE_ID, 0};
        Wire.write(empty, 8);
    }
}

// ── setup ─────────────────────────────────────────────────────────────────────
void setup() {

    Serial.begin(115200);

    Wire.begin(SLAVE_ADDRESS);
    Wire.setClock(400000);       // Fast mode I²C (400 kHz)
    Wire.onRequest(i2c_send_callback);

    // ── Timer1 configuration ─────────────────────────────────────────────
    //
    //  Mode   : Normal (WGM13:0 = 0000) — free-running up-count, no PWM output
    //  Prescaler: 8  →  tick = 8 / 16 MHz = 500 ns
    //  ICP edge : Falling (ICES1 = 0)  — discriminator output is active-low
    //  Noise canceller: ON (ICNC1 = 1) — 4 × 62.5 ns = 250 ns fixed latency,
    //                                     filters sub-250 ns glitches,
    //                                     cancels in Nano-to-Nano differences
    //
    TCCR1A = 0;                    // Normal mode, no compare-output pins
    TCCR1B = _BV(CS11)            // Prescaler = 8 (CS12:CS11:CS10 = 010)
           | _BV(ICNC1);          // Noise canceller ON; ICES1=0 → falling edge
    TIMSK1 = _BV(TOIE1)           // Overflow interrupt
           | _BV(ICIE1);          // Input Capture interrupt

    // ── INT0 (pin 2): pause/resume, any logic change ──────────────────────
    EICRA = (EICRA & ~(_BV(ISC01) | _BV(ISC00))) | _BV(ISC00);  // ISC01:00=01
    EIMSK |= _BV(INT0);

    // ── INT1 (pin 3): sync pulse, rising edge ─────────────────────────────
    EICRA |= _BV(ISC11) | _BV(ISC10);   // ISC11:10 = 11 → rising edge
    EIMSK |= _BV(INT1);

    sei();   // Enable global interrupts
}

// ── loop: all work is interrupt-driven ───────────────────────────────────────
void loop() {}
