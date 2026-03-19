# S03: Diagnose and Fix Pseudorange Computation

**Goal:** Identify and eliminate the 1 km systematic pseudorange offset so computed position matches ground truth within 2-5 meters.

**Demo:** GNSS-SDR solution shows position within 5m of simulated ground truth location; no 1 km offsets observed.

## Must-Haves

- Pseudorange error RMS <5 meters (verified via GNSS-SDR position solution comparison)
- Systematic offset (1 km bias) eliminated
- Ionospheric delay applied correctly per NequickG model
- Satellite position interpolation verified (no gaps between 0.1s updates)
- Ephemeris validity and epoch correctness confirmed

## Tasks

- [ ] **T01: Audit pseudorange computation (channel.cpp::computeRange)**
  Review satellite position calculation, ionospheric correction, relativistic term; compare against reference (e.g., RTKLIB)

- [ ] **T02: Verify ephemeris timing (epoch, week, TOE/TOC)**
  Ensure satellite position is computed at correct transmission time; check for lingering epoch bugs despite prior fixes

- [ ] **T03: Check ionospheric model application (NequickG)**
  Verify parameters, grid lookup, interpolation; compare output to known reference values

- [ ] **T04: Trace pseudorange through a full scenario**
  Generate detailed logs showing pseudorange vs time for one satellite; overlay with GNSS-SDR measurements to spot divergence

- [ ] **T05: Fix root cause (if identified)**
  Implement correction; regenerate signals

- [ ] **T06: Validate position accuracy against GNSS-SDR**
  Generate test signal, check GNSS-SDR solution matches ground truth within 5m

## Files Likely Touched

- `src/channel.cpp` — computeRange(), lines ~100-110
- `src/geodesy.cpp` — satellite_position(), range calculation
- `src/iono.cpp` — NequickG ionospheric delay
- `src/gnss-time.cpp` — epoch/week/TOE arithmetic
- `include/constants.h` — Physical constants (speed of light, relativity, etc.)
