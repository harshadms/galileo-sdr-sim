# S02: Diagnose and Fix I/NAV Message Generation and Timing

**Goal:** Ensure I/NAV messages are generated with correct CRC, bit alignment, and timing synchronization so GNSS-SDR can decode them reliably.

**Demo:** GNSS-SDR logs show successful I/NAV message decoding with valid CRC; no bit errors in message structure.

## Must-Haves

- I/NAV CRC computation verified (bitwise comparison against reference implementation)
- Bit alignment correct per Galileo ICD (no off-by-one errors in spare bits, page structure)
- Message timing aligned with code phase boundaries (4ms symbol boundaries)
- Navigation message page regeneration timing synchronized with satellite transmission time
- GNSS-SDR receiver successfully decodes messages from generated signal (verified by log inspection)

## Tasks

- [ ] **T01: Audit I/NAV generation code for known issues**
  Review inav-msg.cpp, check CRC algorithm, bit packing, spare bit counts against GitHub issues #6, #8

- [ ] **T02: Diagnose timing alignment between message generation and symbol boundaries**
  Trace when messages are regenerated (currently every 30s?) vs when symbols update (every 4ms)

- [ ] **T03: Fix CRC and bit-packing bugs (if found)**
  Implement correct CRC24Q computation, verify spare bit counts per ICD

- [ ] **T04: Synchronize message timing with code phase transitions**
  Ensure page boundaries align with actual transmission time, not arbitrary 30s intervals

- [ ] **T05: Validate against GNSS-SDR decoding**
  Generate test signal, check GNSS-SDR logs for successful I/NAV decode; verify no CRC errors

## Files Likely Touched

- `src/inav-msg.cpp` — CRC computation, bit packing, spare bits
- `src/galileo-sdr.cpp` — Lines 575-585 (page generation timing)
- `include/structures.h` — May need to verify page structure alignment
