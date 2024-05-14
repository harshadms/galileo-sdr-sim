#ifndef GAlSIM_H
#define GAlSIM_H

#include "constants.h"
#include "structures.h"
#include <vector>
#include <pthread.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <signal.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <stdio.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/shm.h>
#include <queue>
#include <map>

using namespace std;

const int WordAllocationE1[30] = {
	2, 4, 6, 7, 8, 17, 19, 16, 0, 0, 1, 3, 5, 0, 16,
	2, 4, 6, 9, 10,17, 19, 16, 0, 0, 1, 3, 5, 0, 16
};

void sigint_handler(int code);
void date2gal(const datetime_t *t, galtime_t *g);

/* Nav message reader function declarations*/
uint32_t data_extract_uint32_t(const std::vector<uint8_t> *str_in, int offset_in, int len);
void data_extract(vector<uint8_t> *str_out, const vector<uint8_t> str_in, int offset_out, int offset_in, int len);
void chartouint8_t(vector<uint8_t> *out, char *in, int len);
int data_extract_signedInt(int len, int data);
void extract_ephemeris(vector<uint8_t> *nav_msg_uint8, ephem_t *eph, uint8_t wordtype, int svid, int *eph_index);

vector<ephem_t> load_ephemeris(const char *fname);
int epoch_matcher(galtime_t obsTime, vector<ephem_t> eph, int index);

vector<uint8_t> load_page(FILE *TV_ptr, int prn, bool advance_fptr, int tow, galtime_t *g);
void cnv_encd(long input_len, int *in_array, int *out_array);

/*! \brief GNSS-Time functions */
void date2gal(const datetime_t *t, galtime_t *g);
void gal2date(const galtime_t *g, datetime_t *t);
double subGalTime(galtime_t g1, galtime_t g0);
galtime_t incGalTime(galtime_t g, double dt);
long get_nanos();
double gps_time(datetime_t *t);

/*! \brief Coord conversion functions - geodesy.cpp */
void neu2azel(double *azel, const double *neu);
void ecef2neu(const double *xyz, double t[3][3], double *neu);
void ltcmat(const double *llh, double t[3][3]);
void llh2xyz(const double *llh, double *xyz);
void xyz2llh(const double *xyz, double *llh);
void satpos(ephem_t eph, galtime_t g, double *pos, double *vel, double *clk);
double dotProd(const double *x1, const double *x2);
double normVect(const double *x);
void subVect(double *y, const double *x1, const double *x2);
double ionosphericDelay(const ionoutc_t *ionoutc, galtime_t g, double *llh, double *azel);
int checkSatVisibility(ephem_t eph, galtime_t g, double *xyz, double elvMask, double *azel, int prn);

/*! \brief Galileo Signal functions - gal-sig.cpp */
void hex_to_binary_converter(short *dest, bool c1, int prn);
void sboc(short *dest, short *in_prn_ca, int len, int m, int n);
void codegen_E1B(short *ca, int prn);
void codegen_E1C(short *ca, int prn);
void computeRange(range_t *rho, ephem_t eph, ionoutc_t *ionoutc, galtime_t g, double xyz[], int prn);
void computeCodePhase(channel_t *chan, range_t rho1, double dt, galtime_t grx);

/*! \brief Galileo Signal functions - inav-msg  .cpp */
void generate_page(galtime_t g, ephem_t *eph, ionoutc_t *ionoutc, int *even_page, int *odd_page, int word_type);
void shift_and_insert(int *arr, int index, uint32_t newElement1, uint32_t newElement2, int arr_size);
void cnv_encd(long input_len, int *in_array, int *out_array);
int generateINavMsg(galtime_t g, channel_t *chan, ephem_t *eph, ionoutc_t *iono);
void generateFrame(int *half_page, int *frame);
void generate_word(ephem_t *Ephemeris, ionoutc_t *IonoGal, unsigned int *EphData, int word_type);
unsigned char GalConvolutionEncode(unsigned char &ConvEncodeBits, unsigned int EncodeWord);
int generateINavMsg_v2(galtime_t g, channel_t *chan, ephem_t *eph, ionoutc_t *iono);

/*! \brief Data conversion datatype.cpp */
#define BIT_ISSET(var, n) !!((long)(var) & (1 << (n)))
#define COMPOSE_BITS(data, start, width) (((data) & ((1UL << width) - 1)) << start)

void encode_int_to_bits(int *page, int *offset, long value, int num_bits);
void encode_double_to_bits(int *page, int *offset, double value, int num_bits);
uint32_t data_extract_uint32_t(const vector<uint8_t> str_in, int offset_in, int len);
void data_extract(vector<uint8_t> *str_out, const vector<uint8_t> str_in, int offset_out, int offset_in, int len);
void chartouint8_t(vector<uint8_t> *out, char *in, int len);
int data_extract_signedInt(int len, int data);

/*! \brief Debug functions - debug.cpp */
int writecsv(int prn, double tow, double values[3], int type);
void print_eph2(ephem_t *eph, int prn);

/*! \brief Galileo functions - galileo-sdr.cpp */
void set_scenario_start_time(galtime_t *g0, galtime_t gmin, galtime_t gmax, datetime_t *t0, datetime_t *tmin, datetime_t *tmax, bool timeoverwrite, ionoutc_t *ionoutc, int neph, vector<ephem_t> eph1[MAX_SAT]);

/*! \brief Channel functions - channel.cpp */
void init_channel(channel_t *chan, int *allocatedSat);
int allocateChannel(channel_t *chan,
                    // map<int, vector<Rinex3Nav::DataGAL>> *navGAL,
                    vector<ephem_t> *eph_vector,
                    ionoutc_t ionoutc, galtime_t grx, double *xyz,
                    double elvMask, map<int, int> *sm,
                    vector<int> current_eph,
                    int *allocatedSat);

/*! \brief Rinex and ephemeris function - rinex.cpp */
//int epoch_matcher(double obsTime, vector<ephem_t> eph, int index);
int readRinexV3(vector<ephem_t> eph_vector[MAX_SAT], ionoutc_t *ionoutc, char *fname);
int readContentsData(char *str, double *data, datetime_t *time, bool read_time);
unsigned char getGalileoUra(double data);
void convertD2E(char *str);

/*! \brief USRP functions */
int init_usrp(usrp_conf_t usrpconf, sim_t *s);

/*! \brief FIFO functions */
int is_fifo_write_ready(sim_t *s);
bool is_finished_generation(sim_t *s);
size_t fifo_read(short *buffer, size_t samples, sim_t *s);
size_t get_sample_length(sim_t *s);
extern void *galileo_task(void *arg);
extern int is_fifo_write_ready(sim_t *s);

/*! \brief datatypes_v2 */
int roundi(double data);
int roundu(double data);
double UnscaleDouble(double value, int scale);
int UnscaleInt(double value, int scale);
unsigned int UnscaleUint(double value, int scale);
long long int UnscaleLong(double value, int scale);
unsigned long long int UnscaleULong(double value, int scale);
int AssignBits(unsigned int Data, int BitNumber, int BitStream[]);
unsigned char ConvolutionEncode(unsigned char EncodeBits);
unsigned int Crc24qEncode(unsigned int *BitStream, int Length);

// /*! \brief FIFO functions - socket.cpp */
// int sockinit(short port);
// void sockclose(int s);
// void socksend(int s, void *data, int siz);
// long int timem();
// int udpinit(short port);
// int udprecv(int s, void *data, int size);
// static long get_nanos2(void);
// void *bitstreamer_thread(vector<queue<int>> &queues, map<int, int> &sm, bool *exit_flag);
// void *dt_thread(void *);
// void *locations_thread(double llhr_init[3]);

#ifndef CRC_CODE_H
#define CRC_CODE_H


const unsigned int Crc24q[256] = {
	0x00000000, 0x864CFB00, 0x8AD50D00, 0x0C99F600,
	0x93E6E100, 0x15AA1A00, 0x1933EC00, 0x9F7F1700,
	0xA1813900, 0x27CDC200, 0x2B543400, 0xAD18CF00,
	0x3267D800, 0xB42B2300, 0xB8B2D500, 0x3EFE2E00,
	0xC54E8900, 0x43027200, 0x4F9B8400, 0xC9D77F00,
	0x56A86800, 0xD0E49300, 0xDC7D6500, 0x5A319E00,
	0x64CFB000, 0xE2834B00, 0xEE1ABD00, 0x68564600,
	0xF7295100, 0x7165AA00, 0x7DFC5C00, 0xFBB0A700,
	0x0CD1E900, 0x8A9D1200, 0x8604E400, 0x00481F00,
	0x9F370800, 0x197BF300, 0x15E20500, 0x93AEFE00,
	0xAD50D000, 0x2B1C2B00, 0x2785DD00, 0xA1C92600,
	0x3EB63100, 0xB8FACA00, 0xB4633C00, 0x322FC700,
	0xC99F6000, 0x4FD39B00, 0x434A6D00, 0xC5069600,
	0x5A798100, 0xDC357A00, 0xD0AC8C00, 0x56E07700,
	0x681E5900, 0xEE52A200, 0xE2CB5400, 0x6487AF00,
	0xFBF8B800, 0x7DB44300, 0x712DB500, 0xF7614E00,
	0x19A3D200, 0x9FEF2900, 0x9376DF00, 0x153A2400,
	0x8A453300, 0x0C09C800, 0x00903E00, 0x86DCC500,
	0xB822EB00, 0x3E6E1000, 0x32F7E600, 0xB4BB1D00,
	0x2BC40A00, 0xAD88F100, 0xA1110700, 0x275DFC00,
	0xDCED5B00, 0x5AA1A000, 0x56385600, 0xD074AD00,
	0x4F0BBA00, 0xC9474100, 0xC5DEB700, 0x43924C00,
	0x7D6C6200, 0xFB209900, 0xF7B96F00, 0x71F59400,
	0xEE8A8300, 0x68C67800, 0x645F8E00, 0xE2137500,
	0x15723B00, 0x933EC000, 0x9FA73600, 0x19EBCD00,
	0x8694DA00, 0x00D82100, 0x0C41D700, 0x8A0D2C00,
	0xB4F30200, 0x32BFF900, 0x3E260F00, 0xB86AF400,
	0x2715E300, 0xA1591800, 0xADC0EE00, 0x2B8C1500,
	0xD03CB200, 0x56704900, 0x5AE9BF00, 0xDCA54400,
	0x43DA5300, 0xC596A800, 0xC90F5E00, 0x4F43A500,
	0x71BD8B00, 0xF7F17000, 0xFB688600, 0x7D247D00,
	0xE25B6A00, 0x64179100, 0x688E6700, 0xEEC29C00,
	0x3347A400, 0xB50B5F00, 0xB992A900, 0x3FDE5200,
	0xA0A14500, 0x26EDBE00, 0x2A744800, 0xAC38B300,
	0x92C69D00, 0x148A6600, 0x18139000, 0x9E5F6B00,
	0x01207C00, 0x876C8700, 0x8BF57100, 0x0DB98A00,
	0xF6092D00, 0x7045D600, 0x7CDC2000, 0xFA90DB00,
	0x65EFCC00, 0xE3A33700, 0xEF3AC100, 0x69763A00,
	0x57881400, 0xD1C4EF00, 0xDD5D1900, 0x5B11E200,
	0xC46EF500, 0x42220E00, 0x4EBBF800, 0xC8F70300,
	0x3F964D00, 0xB9DAB600, 0xB5434000, 0x330FBB00,
	0xAC70AC00, 0x2A3C5700, 0x26A5A100, 0xA0E95A00,
	0x9E177400, 0x185B8F00, 0x14C27900, 0x928E8200,
	0x0DF19500, 0x8BBD6E00, 0x87249800, 0x01686300,
	0xFAD8C400, 0x7C943F00, 0x700DC900, 0xF6413200,
	0x693E2500, 0xEF72DE00, 0xE3EB2800, 0x65A7D300,
	0x5B59FD00, 0xDD150600, 0xD18CF000, 0x57C00B00,
	0xC8BF1C00, 0x4EF3E700, 0x426A1100, 0xC426EA00,
	0x2AE47600, 0xACA88D00, 0xA0317B00, 0x267D8000,
	0xB9029700, 0x3F4E6C00, 0x33D79A00, 0xB59B6100,
	0x8B654F00, 0x0D29B400, 0x01B04200, 0x87FCB900,
	0x1883AE00, 0x9ECF5500, 0x9256A300, 0x141A5800,
	0xEFAAFF00, 0x69E60400, 0x657FF200, 0xE3330900,
	0x7C4C1E00, 0xFA00E500, 0xF6991300, 0x70D5E800,
	0x4E2BC600, 0xC8673D00, 0xC4FECB00, 0x42B23000,
	0xDDCD2700, 0x5B81DC00, 0x57182A00, 0xD154D100,
	0x26359F00, 0xA0796400, 0xACE09200, 0x2AAC6900,
	0xB5D37E00, 0x339F8500, 0x3F067300, 0xB94A8800,
	0x87B4A600, 0x01F85D00, 0x0D61AB00, 0x8B2D5000,
	0x14524700, 0x921EBC00, 0x9E874A00, 0x18CBB100,
	0xE37B1600, 0x6537ED00, 0x69AE1B00, 0xEFE2E000,
	0x709DF700, 0xF6D10C00, 0xFA48FA00, 0x7C040100,
	0x42FA2F00, 0xC4B6D400, 0xC82F2200, 0x4E63D900,
	0xD11CCE00, 0x57503500, 0x5BC9C300, 0xDD853800
};

#endif
#endif