# M001: Fix Core Signal Stability and Precision

**Gathered:** 2026-03-19
**Status:** Ready for planning

## Project Description

Galileo E1 SDR signal simulator currently generates signals that GNSS receivers can detect (PRN acquisition) but cannot reliably lock onto or use for positioning. Receivers show: 1 km pseudorange offsets, frequent tracking loss, 70+ minute TTFF (should be <40s), position jumps, and low C/NO.

## Why This Milestone

The simulator is fundamentally broken for its primary use case — testing GNSS receivers. M001 fixes the core signal generation to enable reliable receiver lock, accurate positioning, and fast acquisition. This is the foundation all downstream work depends on.

## User-Visible Outcome

### When this milestone is complete, the user can:

- Generate a 5-minute Galileo E1 signal file
- Feed it to GNSS-SDR and achieve receiver lock within 40 seconds
- See computed position within 2-5 meters of the simulated ground truth
- Maintain stable tracking and position for the entire 5-minute duration without jumps or loss of lock

### Entry point / environment

- Entry point: `./usrp_galileo -e <rinex> -l <lat,lon,hgt> -U 1 -b 1 -d 300 -o output.bin`
- Environment: Local dev (Ubuntu Linux), GNSS-SDR post-processing
- Live dependencies: none (offline file generation)

## Completion Class

- Contract complete means: All 9 active requirements have verification passing (C/NO >30 dB-Hz, TTFF <40s, position <5m, no jumps for 5 min, 16-bit output, I/NAV decodes, pseudorange <5m error, smooth tracking)
- Integration complete means: Generated file passes end-to-end test with GNSS-SDR; receiver locks and positions correctly
- Operational complete means: Simulator generates file consistently; no crashes or edge cases on different ephemeris/location inputs

## Final Integrated Acceptance

To call this milestone complete, we must prove:

1. GNSS-SDR successfully locks onto 4+ satellites within 40 seconds (C/NO >30 dB-Hz for each)
2. GNSS-SDR computes position with <5m error compared to ground truth (within first 2 minutes of lock)
3. Position remains stable for full 5-minute scenario — no jumps >1m, no loss of lock
4. Generated output file is 16-bit I/Q (verified by file format inspection)
5. Navigation messages decode correctly in GNSS-SDR log (I/NAV CRC passes, bit alignment correct)

## Risks and Unknowns

- **Pseudorange offset root cause unclear** — Could be epoch/week calculation, ionospheric model, satellite position interpolation, or relativistic correction. Need systematic diagnosis.
- **I/NAV message timing alignment** — Regenerated every 30s but updates every symbol (4ms). Timing mismatch could cause decoding failures even with correct CRC.
- **Code phase discontinuities on channel allocation** — `set_code_phase` logic may not initialize correctly; could cause tracking jumps when new channels open.
- **Frequency update rate too coarse** — 0.1s updates to Doppler may cause drift; real receivers track continuous phase.
- **8-bit quantization masking real issues** — Fixing to 16-bit will reveal signal quality problems currently hidden.

## Existing Codebase / Prior Art

- `channel.cpp::allocateChannel()` — Channel allocation, pseudorange computation, visibility checks
- `channel.cpp::computeRange()` — Pseudorange calculation; calls geodesy, ephemeris lookup, ionospheric correction
- `gal-sig.cpp::computeCodePhase()` — Code phase and Doppler update from pseudorange rate
- `inav-msg.cpp::generateINavMsg()` — I/NAV message generation; prior fixes in commits dd22226, 0d8d2bc, fa59d4c
- `galileo-sdr.cpp` (main loop, lines 460-730) — Sample generation, channel updates, file output
- `geodesy.cpp` — Satellite position, range, elevation/azimuth computation
- `gnss-time.cpp` — Galileo time arithmetic (epoch fixed in commit 377a5e7)
- Test file: `rinex_files/week171.rnx` — Real Galileo ephemeris for 2021-06-20

## Relevant Requirements

- R001 — Receiver can acquire and lock (primary acceptance)
- R002 — TTFF <40s (acceptance criterion)
- R003 — Position accuracy <5m (acceptance criterion)
- R004 — Stable tracking 5+ min (acceptance criterion)
- R005 — 16-bit I/Q precision (enabling requirement)
- R006 — I/NAV message integrity (enabling requirement)
- R007 — Pseudorange accuracy (core requirement, likely root cause of offset)
- R008 — Code phase / carrier sync (likely root cause of jumps)
- R009 — Navigation message timing (enabling requirement)

## Scope

### In Scope

- Diagnose pseudorange offset root cause and fix
- Validate/fix I/NAV message generation and timing
- Upgrade to 16-bit I/Q output
- Ensure code phase and carrier frequency transitions are smooth (no discontinuities)
- Validate against GNSS-SDR on test ephemeris
- Improve TTFF and tracking robustness

### Out of Scope / Non-Goals

- Real-time USRP transmission (file output only; USRP testing deferred to M002)
- Real receiver validation (u-Blox, Septentrio deferred to M002)
- GPS L1 C/A signal generation (future work, M003)
- Multi-SDR support (future work, M003)
- Performance optimization beyond what's needed for correctness

## Technical Constraints

- Must remain compatible with existing RINEX ephemeris loader
- Must preserve channel allocation and visibility logic (it works correctly)
- Simulation step size is 0.1 seconds (fixed; don't change)
- Output format can change (e.g., 8-bit to 16-bit) as long as GNSS-SDR can read it
- No external dependencies beyond what's already required (UHD, Boost, Curses, glib2.0)

## Integration Points

- **GNSS-SDR** — Reads generated binary I/Q file; must decode I/NAV, track satellites, compute position
- **RINEX ephemeris files** — Input; no changes to parser needed
- **Receiver position (ground truth)** — Injected via command-line; used for pseudorange computation and verification

## Open Questions

- What exact ionospheric model (NequickG) parameters are we using? Is the model correctly applied?
- How sensitive is pseudorange to satellite position interpolation? Should we compute position more frequently than 0.1s?
- Are there edge cases in week/epoch arithmetic that still cause time errors despite prior fixes?
