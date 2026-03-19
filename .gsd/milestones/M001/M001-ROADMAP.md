# M001: Fix Core Signal Stability and Precision

**Vision:** Generate reliable Galileo E1 signals that GNSS receivers can lock onto, decode, and use for accurate positioning with fast convergence.

## Success Criteria

- GNSS-SDR achieves receiver lock within 40 seconds on 4+ satellites
- Computed position is within 2-5 meters of simulated ground truth
- Position remains stable for 5+ minutes with no jumps or loss of lock
- I/NAV messages decode correctly (CRC validation passes)
- Signal output is 16-bit I/Q with no quantization clipping

## Key Risks / Unknowns

- **Pseudorange offset (1 km)** — Root cause unknown; could be in ephemeris timing, interpolation, ionospheric model, or relativistic correction
- **I/NAV message timing** — Regenerated every 30s but symbol rate is 4ms; timing mismatch may prevent decoding
- **Code phase discontinuities** — `set_code_phase` logic may cause tracking jumps when channels are allocated
- **Frequency update coarseness** — 0.1s updates to Doppler may drift; could contribute to TTFF delay

## Proof Strategy

- Pseudorange offset → retire in S03 by proving computed position matches ground truth (<5m error)
- I/NAV timing/CRC → retire in S02 by decoding messages successfully in GNSS-SDR receiver
- Code phase smoothness → retire in S04 by observing no position jumps during continuous 5-min tracking
- 16-bit precision → retire in S01 by inspection and verifying no clipping in generated samples

## Verification Classes

- Contract verification: GNSS-SDR receiver lock test (C/NO, TTFF, position error, tracking stability)
- Integration verification: End-to-end signal generation → GNSS-SDR processing → position output
- Operational verification: None beyond correctness (simulator is offline file generator)
- UAT / human verification: Visual inspection of GNSS-SDR console output (lock status, solution convergence)

## Milestone Definition of Done

This milestone is complete only when all are true:

- All 5 slices have verified outputs (16-bit samples, I/NAV decoded, pseudorange accurate, smooth tracking, fast TTFF)
- GNSS-SDR successfully locks within 40 seconds on test ephemeris (week171.rnx)
- Position converges to within 5m of ground truth
- Continuous 5-minute scenario runs without loss of lock or position jumps
- Generated signal file is reproducible (same ephemeris/location produces same output)

## Requirement Coverage

- Covers: R001, R002, R003, R004, R005, R006, R007, R008, R009
- Partially covers: none
- Leaves for later: R010 (real receiver validation, M002)

## Slices

- [ ] **S01: Upgrade to 16-bit I/Q precision** `risk:low` `depends:[]`
  > After this: Signal output file contains 16-bit signed I/Q samples instead of 8-bit; no quantization clipping observed in hex dump

- [ ] **S02: Diagnose and fix I/NAV message generation and timing** `risk:high` `depends:[S01]`
  > After this: GNSS-SDR successfully decodes I/NAV messages from generated signals; CRC validation passes in receiver log

- [ ] **S03: Diagnose and fix pseudorange computation** `risk:high` `depends:[S02]`
  > After this: Computed position matches ground truth within 5m; no 1 km offsets observed in GNSS-SDR solution

- [ ] **S04: Ensure smooth code phase and carrier frequency transitions** `risk:medium` `depends:[S03]`
  > After this: Position remains stable during entire 5-min scenario; no jumps >1m observed; tracking does not drop

- [ ] **S05: Validate TTFF and overall receiver performance** `risk:medium` `depends:[S04]`
  > After this: GNSS-SDR achieves lock within 40 seconds; C/NO >30 dB-Hz for 4+ satellites; position converges in <2 min

## Boundary Map

### S01 → S02
Produces:
  - 16-bit signed I/Q samples in binary file (I1 Q1 I2 Q2 ... In Qn format)
  - Removal of quantization clipping logic; raw signal scaled appropriately

Consumes: nothing (first slice, builds on existing output code)

### S01 → S03
Produces:
  - Same as S01 (16-bit samples enable diagnosis of pseudorange issues)

Consumes: nothing

### S02 → S03
Produces:
  - Correct I/NAV message bits, properly aligned and timed
  - CRC that validates in GNSS-SDR receiver
  - Page boundaries correctly marked with transmission time

Consumes from S01:
  - 16-bit I/Q samples (better signal fidelity for message decoding)

### S03 → S04
Produces:
  - Pseudorange computed with <5m RMS error
  - Ionospheric correction applied correctly
  - Satellite position accurately interpolated

Consumes from S02:
  - Correct I/NAV messages (receiver uses these to decode ephemeris and position)

### S04 → S05
Produces:
  - Code phase initialized smoothly on channel allocation (no discontinuities)
  - Carrier frequency updated without jumps
  - Doppler tracking stable across entire scenario

Consumes from S03:
  - Accurate pseudorange (enables smooth Doppler estimation)

### S05 (final integration)
Produces:
  - Complete, validated signal file that GNSS-SDR can process end-to-end

Consumes from all prior slices:
  - 16-bit samples, correct I/NAV, accurate pseudorange, smooth tracking
