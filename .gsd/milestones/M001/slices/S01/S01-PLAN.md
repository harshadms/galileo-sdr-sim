# S01: Upgrade to 16-bit I/Q Precision

**Goal:** Change output format from 8-bit to 16-bit signed I/Q samples to improve signal fidelity and remove quantization noise as a confounding variable in debugging.

**Demo:** Generated signal file contains 16-bit signed integers; hex dump shows full range of values (-32768 to 32767) with no clipping.

## Must-Haves

- 16-bit signed I/Q output format confirmed in generated file
- No quantization clipping to ±127 range observed
- File size approximately double the 8-bit version (for same sample count)
- GNSS-SDR can successfully read and process the 16-bit file
- Backward compatibility check: other slices don't break

## Tasks

- [ ] **T01: Change output data type from signed char (8-bit) to short (16-bit)**
  Modify galileo-sdr.cpp buffer allocation and output logic

- [ ] **T02: Update scaling/normalization for 16-bit range**
  Ensure signal amplitude uses full 16-bit range without clipping

- [ ] **T03: Verify file format and test with GNSS-SDR**
  Generate test file, inspect hex dump, run through GNSS-SDR to confirm it reads correctly

## Files Likely Touched

- `src/galileo-sdr.cpp` — Lines 350-353 (buffer allocation), 554-620 (clipping logic), 637 (fwrite)
- `include/structures.h` — May need buffer type changes if using typedef
