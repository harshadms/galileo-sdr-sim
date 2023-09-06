# GALILEO-SDR-SIM

GALILEO-SDR-SIM generates Galileo baseband signal data streams, which can be converted to RF using software-defined radio (SDR) platforms.

### To Do List:
- [ ] Phase 1: Implementing the galileo modulation scheme and spreading. Galileo signals are modulated by a new modulation technique;
which is called Binary Offset Carrier (BOC) modulation. The Navigation data in this phase will be random bits. 
- [ ] Testing the the signal modulation part. The goal is to verifying the modulation section. The gnss-sdr-sim should be able to track the generated signal. 
- [ ] Phase 2: Generating the navigation bits. Reading from the Galileo rinex file or RTCM stream; Then generating the raw bits for Galileo Navigation messages. 
- [ ] Testing phase 2 implementation 
- [ ] Phase 3: adding pseudoranges for one satellite. 
- [ ] Phase 4: adding pseudoranges for all visible satellites.
- [ ] Final test 

### Some notes:

gal-sdr-sim.c:190 - readTestvectToEph()
- This function reads and extracts ephemeris data from the test vector file. (called just once)

testvect_update_navm() - readTestVectors.c
- Accepts a TOW as an argument and populates navigation message for each satellite from the test vector file

Our task:

- Read RINEX file for current TOW
- Get ephemeris data 
- Start from wordtype 1 (IODNav values) (Refer to galileo_inav_message.cc:839 Galileo_Inav_Message::page_jk_decoder(const char* data_jk) from GNSS-SDR)
- Refer to chartoU8(), data_extract(), and data_extract_tU32() functions from conversion.c. These functions are responsible for converting test vector data to ephemeris and to bits. 
- **Use - rinex_files/galileo_2022-02-20-0800hrs.rnx**

**Task**

- [ ] Modify these functions to convert from RINEX to nav message bits. read_rinex() -> chartoU8() -> nav bits 