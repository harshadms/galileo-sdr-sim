# S05: Validate TTFF and Overall Receiver Performance

**Goal:** Confirm all prior slices work together to deliver fast acquisition and reliable tracking in a realistic scenario.

**Demo:** GNSS-SDR acquires satellites within 40 seconds, achieves C/NO >30 dB-Hz on 4+ satellites, computes position within 2 minutes.

## Must-Haves

- TTFF <40 seconds (measured from first sample to first valid position fix)
- C/NO >30 dB-Hz on minimum 4 satellites (per GNSS-SDR tracking output)
- Position convergence <2 minutes after first lock
- RMS position error <5 meters (vs ground truth)
- Continuous tracking for full 5-minute scenario (no loss of lock, no jumps)
- Reproducibility: same ephemeris/location produces identical signal (bit-for-bit, except timestamp differences)

## Tasks

- [ ] **T01: Set up automated test harness**
  Script to generate signal, run GNSS-SDR, parse results, measure TTFF, C/NO, position error

- [ ] **T02: Test on reference ephemeris (week171.rnx)**
  Generate 5-minute signal at known location; run GNSS-SDR; capture metrics

- [ ] **T03: Analyze results against acceptance criteria**
  If TTFF >40s or C/NO <30 dB-Hz, backtrack to prior slices; identify remaining issue

- [ ] **T04: Stress test: multiple ephemeris files and locations**
  Verify slices work across different scenarios, not just the test case

- [ ] **T05: Final integration verification**
  Confirm no regressions from S01-S04; measure end-to-end metrics

- [ ] **T06: Document results and acceptance**
  Write test report; confirm M001 acceptance criteria are met

## Files Likely Touched

- None (verification slice — tests existing code)
- May create test scripts in a new `tests/` directory

## Verification Approach

- **GNSS-SDR console output parsing** — Extract lock status, C/NO, position, time tags
- **Position error calculation** — Compare computed solution to ground truth (input location)
- **TTFF measurement** — Time from scenario start to first valid position fix
- **Tracking continuity** — Check for gaps in solution output (indicates loss of lock)
