# Galileo SDR Simulator

## What This Is

A Galileo E1B/C baseband signal generator that produces realistic IQ samples suitable for real-time or offline processing by GNSS receivers. Currently generates signals that receivers can detect (PRN acquisition works) but cannot reliably lock onto or decode (low C/NO, frequent tracking loss, poor positioning).

## Core Value

Receivers can generate stable Galileo E1 signals with accurate pseudorange, reliable I/NAV messages, and smooth tracking characteristics — enabling offline validation, real-time testing, and research on Galileo signal processing.

## Current State

**What works:**
- Ephemeris loading from RINEX files
- PRN code generation (E1B/E1C)
- Basic I/NAV message generation and transmission
- Channel allocation and signal modulation
- File output (binary I/Q samples)
- GNSS-SDR receiver can detect PRN numbers

**What's broken:**
- Receivers cannot lock (frequent loss of tracking)
- Pseudorange offsets up to 1 km (should be <5m)
- TTFF 70+ minutes (should be <40s)
- Position jumps during scenario
- 8-bit I/Q quantization reduces signal quality

**Known issues documented:**
- I/NAV message CRC and bit-packing bugs (partially fixed in recent commits)
- Epoch/week time calculation bugs (fixed in commit 377a5e7)
- Code phase initialization ambiguity (addressed in commit e4c6654)

## Architecture / Key Patterns

**Signal generation pipeline:**
1. Load ephemeris from RINEX (includes epoch, clock correction, orbital elements)
2. For each 0.1-second time step:
   - Allocate channels for visible satellites (elevation mask 10°)
   - Compute pseudorange for each satellite (geodesy + ionospheric corrections)
   - Compute Doppler shift and code phase from pseudorange rate
   - Generate I/NAV navigation messages
   - Modulate baseband signal: E1B (data + pilot) at 2.6 MHz, BOC(1,1)
   - Write I/Q samples to file (currently 8-bit, will upgrade to 16-bit)

**Key modules:**
- `channel.cpp` — channel allocation, visibility checks, pseudorange computation
- `gal-sig.cpp` — code generation, phase computation, signal modulation
- `inav-msg.cpp` — I/NAV message generation and bit packing
- `galileo-sdr.cpp` — main loop, sample generation, time management
- `geodesy.cpp` — Earth geometry, satellite position, range computation
- `rinex.cpp` — ephemeris parsing
- `iono.cpp` — ionospheric delay (NequickG model, recently refactored)

**Established patterns:**
- Galileo system time represented as (week, second) tuple
- Pseudorange computed using iterative SV transmission time calculation
- Code phase tracked per-symbol (4ms periods)
- Navigation messages updated every 2 seconds (page boundary)

## Capability Contract

See `.gsd/REQUIREMENTS.md` for explicit capability contract and coverage mapping.

## Milestone Sequence

- [ ] M001: Fix core signal stability and precision — Enable reliable receiver lock, accurate positioning, fast TTFF
- [ ] M002: Validate against live receivers — Prove signals work with real GNSS-SDR, u-Blox, Septentrio

Current commit: `0c06823f84f14b9a036c7f8a9b54fca5b2f1b201`
