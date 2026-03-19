# Requirements

This file is the explicit capability and coverage contract for the project.

## Active

### R001 — Receiver can acquire and lock onto satellites
- Class: primary-user-loop
- Status: active
- Description: GNSS receiver must achieve lock on Galileo E1 signals with stable C/NO (>30 dB-Hz for at least 4 satellites visible)
- Why it matters: Without lock, receiver cannot extract navigation data or compute position
- Source: user
- Primary owning slice: M001/S02
- Supporting slices: M001/S01, M001/S03
- Validation: unmapped
- Notes: Currently fails — receivers detect PRN but lose lock within seconds

### R002 — Time-to-First-Fix under 40 seconds
- Class: quality-attribute
- Status: active
- Description: From cold start, receiver must compute first position fix in <40 seconds
- Why it matters: Affects usability and practical testability
- Source: user
- Primary owning slice: M001/S04
- Supporting slices: M001/S02, M001/S03
- Validation: unmapped
- Notes: Currently 70+ minutes; root cause likely pseudorange offset or tracking instability

### R003 — Position accuracy within 2-5 meters of ground truth
- Class: quality-attribute
- Status: active
- Description: Computed position must be within 2-5m of simulated reference location (after lock achieved)
- Why it matters: Validates pseudorange computation and ionospheric correction accuracy
- Source: user
- Primary owning slice: M001/S03
- Supporting slices: M001/S02
- Validation: unmapped
- Notes: Currently shows 1 km offsets; suggests systematic pseudorange error

### R004 — Stable continuous tracking for 5+ minutes
- Class: quality-attribute
- Status: active
- Description: Once locked, receiver must maintain lock and position fix for entire 5+ minute scenario without jumps or loss of lock
- Why it matters: Demonstrates signal consistency and navigation message integrity
- Source: user
- Primary owning slice: M001/S03
- Supporting slices: M001/S02, M001/S04
- Validation: unmapped
- Notes: Currently sees frequent tracking loss and position jumps

### R005 — 16-bit I/Q precision
- Class: quality-attribute
- Status: active
- Description: Signal samples must be 16-bit signed integers, not 8-bit (removes quantization noise source)
- Why it matters: Improves SNR, enables reliable tracking and message decoding
- Source: user
- Primary owning slice: M001/S01
- Supporting slices: none
- Validation: unmapped
- Notes: Currently 8-bit; adds quantization noise that masks actual signal quality issues

### R006 — I/NAV message integrity
- Class: core-capability
- Status: active
- Description: Generated I/NAV messages must have correct CRC, bit packing, spare bit counts, and bit alignment per Galileo IS-GPS-200 (ICD)
- Why it matters: Receivers cannot decode corrupted messages; causes tracking loss and wrong position
- Source: user (and GitHub issues #6, #8)
- Primary owning slice: M001/S02
- Supporting slices: none
- Validation: unmapped
- Notes: Prior fixes applied (commits dd22226, 0d8d2bc, fa59d4c); needs validation against live decoding

### R007 — Pseudorange computation accuracy
- Class: core-capability
- Status: active
- Description: Computed pseudorange must account for: satellite position (ephemeris), receiver position, ionospheric delay, and relativistic corrections
- Why it matters: Systematic pseudorange error directly causes 1 km position offset and slow convergence
- Source: user
- Primary owning slice: M001/S03
- Supporting slices: M001/S02
- Validation: unmapped
- Notes: Suspected issues: epoch time, interpolation between 0.1s updates, or ionospheric model

### R008 — Code phase and carrier frequency synchronization
- Class: core-capability
- Status: active
- Description: Code phase and carrier frequency must transition smoothly without discontinuities; correct initialization on channel allocation
- Why it matters: Discontinuities cause tracking jumps and loss of lock
- Source: user (position jumps observed)
- Primary owning slice: M001/S04
- Supporting slices: M001/S02, M001/S03
- Validation: unmapped
- Notes: Potential issues: set_code_phase logic, rho0 initialization, frequency update timing

### R009 — Navigation message timing correctness
- Class: core-capability
- Status: active
- Description: I/NAV messages must have correct transmission time alignment with code phase and satellite ephemeris
- Why it matters: Receiver must correlate message bits with precise timing to decode correctly
- Source: user
- Primary owning slice: M001/S02
- Supporting slices: none
- Validation: unmapped
- Notes: Currently regenerated every 30s; may not align correctly with 0.1s step timing

### R010 — Signal works with real receivers
- Class: differentiator
- Status: deferred
- Description: Generated signals must pass validation against u-Blox, Septentrio, or other real GNSS receivers (not just GNSS-SDR)
- Why it matters: Proves signal is truly realistic and production-grade
- Source: user (optional but valuable)
- Primary owning slice: M002/S01
- Supporting slices: none
- Validation: unmapped
- Notes: Deferred to M002; M001 focuses on GNSS-SDR + accuracy

## Traceability

| ID | Class | Status | Primary owner | Supporting | Proof |
|---|---|---|---|---|---|
| R001 | primary-user-loop | active | M001/S02 | S01, S03 | unmapped |
| R002 | quality-attribute | active | M001/S04 | S02, S03 | unmapped |
| R003 | quality-attribute | active | M001/S03 | S02 | unmapped |
| R004 | quality-attribute | active | M001/S03 | S02, S04 | unmapped |
| R005 | quality-attribute | active | M001/S01 | none | unmapped |
| R006 | core-capability | active | M001/S02 | none | unmapped |
| R007 | core-capability | active | M001/S03 | S02 | unmapped |
| R008 | core-capability | active | M001/S04 | S02, S03 | unmapped |
| R009 | core-capability | active | M001/S02 | none | unmapped |
| R010 | differentiator | deferred | M002/S01 | none | unmapped |

## Coverage Summary

- Active requirements: 9
- Mapped to slices: 9 (all in M001)
- Validated: 0
- Unmapped active requirements: 0
