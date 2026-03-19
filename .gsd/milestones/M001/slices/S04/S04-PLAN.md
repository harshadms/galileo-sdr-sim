# S04: Ensure Smooth Code Phase and Carrier Frequency Transitions

**Goal:** Eliminate position jumps and tracking discontinuities by ensuring code phase and carrier frequency update smoothly without step changes.

**Demo:** Position output remains within 1m of mean value during entire 5-minute scenario; no sudden jumps or tracking loss observed.

## Must-Haves

- Code phase initialized correctly on channel allocation (no discontinuities when set_code_phase flag is used)
- Carrier frequency updates smoothly based on pseudorange rate (no step changes)
- Position jump threshold: no single fix more than 1m away from previous fix
- Tracking loss: 0 occurrences during 5-minute scenario
- Doppler estimation stable (frequency updates follow smooth pseudorange rate curve)

## Tasks

- [ ] **T01: Review code phase initialization logic (set_code_phase in galileo-sdr.cpp)**
  Verify calculation at lines 493-515; check boundary conditions and rounding

- [ ] **T02: Trace code phase and frequency updates through a scenario**
  Log code_phase, f_code, f_carr, and rho0 for all channels during full run; plot to identify jumps

- [ ] **T03: Check rho0 state variable correctness**
  Verify that rho0 (previous pseudorange) is updated consistently; ensure no stale values cause discontinuities

- [ ] **T04: Review Doppler calculation (f_carr = -rhorate / LAMBDA_E1)**
  Verify pseudorange rate computation; check for division errors or sign mistakes

- [ ] **T05: Fix any identified discontinuities**
  Smooth transitions, fix rounding errors, ensure state consistency

- [ ] **T06: Validate tracking stability against GNSS-SDR**
  Generate test signal, check GNSS-SDR tracking logs for continuous lock and smooth position (no jumps)

## Files Likely Touched

- `src/galileo-sdr.cpp` — Lines 493-515 (set_code_phase initialization), 486-492 (code phase/frequency update), 612-618 (phase advance per sample)
- `src/gal-sig.cpp` — computeCodePhase() function
- `include/structures.h` — channel_t structure (rho0, code_phase, f_carr, set_code_phase fields)
