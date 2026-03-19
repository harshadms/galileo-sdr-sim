# Decisions Register

<!-- Append-only. Never edit or remove existing rows.
     To reverse a decision, add a new row that supersedes it.
     Read this file at the start of any planning or research phase. -->

| # | When | Scope | Decision | Choice | Rationale | Revisable? |
|---|------|-------|----------|--------|-----------|------------|
| D001 | M001 | scope | M001 focus | Core signal stability only (no USRP, no real receivers yet) | USRP and real receiver validation blocked by signal quality; defer to M002 | Yes — if signal fixes are simpler than expected |
| D002 | M001 | precision | Output format | Upgrade from 8-bit to 16-bit I/Q | 8-bit quantization adds noise that masks real signal issues; 16-bit removes this confound | No |
| D003 | M001 | strategy | Diagnosis order | Fix I/NAV → pseudorange → code phase/tracking | I/NAV integrity is prerequisite for receiver to decode anything; pseudorange affects position; tracking affects lock | Yes — if diagnosis reveals different root cause order |
| D004 | M001 | dependency | S01 independence | S01 (16-bit) has no dependencies | Can change output format independently; enables all downstream diagnosis | No |
| D005 | M001 | requirement | TTFF threshold | <40 seconds | User requirement; current 70+ minutes unacceptable | No |
| D006 | M001 | requirement | Position accuracy | 2-5 meters | User requirement; current 1 km systematic error unacceptable | No |
| D007 | M001 | requirement | Tracking stability | 5+ minutes, no jumps | User requirement; enables practical testing scenario | No |
