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

const int WordAllocationE1[15] = {
	2, 4, 6, 7, 8, 17, 19, 16, 0, 0, 1, 3, 5, 0, 16,
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
void init_usrp(usrp_conf_t usrpconf, sim_t *s);

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

// static const unsigned int Crc24q[] = {
//     0x000000, 0x864CFB, 0x8AD50D, 0x0C99F6, 0x93E6E1, 0x15AA1A, 0x1933EC, 0x9F7F17,
//     0xA18139, 0x27CDC2, 0x2B5434, 0xAD18CF, 0x3267D8, 0xB42B23, 0xB8B2D5, 0x3EFE2E,
//     0xC54E89, 0x430272, 0x4F9B84, 0xC9D77F, 0x56A868, 0xD0E493, 0xDC7D65, 0x5A319E,
//     0x64CFB0, 0xE2834B, 0xEE1ABD, 0x685646, 0xF72951, 0x7165AA, 0x7DFC5C, 0xFBB0A7,
//     0x0CD1E9, 0x8A9D12, 0x8604E4, 0x00481F, 0x9F3708, 0x197BF3, 0x15E205, 0x93AEFE,
//     0xAD50D0, 0x2B1C2B, 0x2785DD, 0xA1C926, 0x3EB631, 0xB8FACA, 0xB4633C, 0x322FC7,
//     0xC99F60, 0x4FD39B, 0x434A6D, 0xC50696, 0x5A7981, 0xDC357A, 0xD0AC8C, 0x56E077,
//     0x681E59, 0xEE52A2, 0xE2CB54, 0x6487AF, 0xFBF8B8, 0x7DB443, 0x712DB5, 0xF7614E,
//     0x19A3D2, 0x9FEF29, 0x9376DF, 0x153A24, 0x8A4533, 0x0C09C8, 0x00903E, 0x86DCC5,
//     0xB822EB, 0x3E6E10, 0x32F7E6, 0xB4BB1D, 0x2BC40A, 0xAD88F1, 0xA11107, 0x275DFC,
//     0xDCED5B, 0x5AA1A0, 0x563856, 0xD074AD, 0x4F0BBA, 0xC94741, 0xC5DEB7, 0x43924C,
//     0x7D6C62, 0xFB2099, 0xF7B96F, 0x71F594, 0xEE8A83, 0x68C678, 0x645F8E, 0xE21375,
//     0x15723B, 0x933EC0, 0x9FA736, 0x19EBCD, 0x8694DA, 0x00D821, 0x0C41D7, 0x8A0D2C,
//     0xB4F302, 0x32BFF9, 0x3E260F, 0xB86AF4, 0x2715E3, 0xA15918, 0xADC0EE, 0x2B8C15,
//     0xD03CB2, 0x567049, 0x5AE9BF, 0xDCA544, 0x43DA53, 0xC596A8, 0xC90F5E, 0x4F43A5,
//     0x71BD8B, 0xF7F170, 0xFB6886, 0x7D247D, 0xE25B6A, 0x641791, 0x688E67, 0xEEC29C,
//     0x3347A4, 0xB50B5F, 0xB992A9, 0x3FDE52, 0xA0A145, 0x26EDBE, 0x2A7448, 0xAC38B3,
//     0x92C69D, 0x148A66, 0x181390, 0x9E5F6B, 0x01207C, 0x876C87, 0x8BF571, 0x0DB98A,
//     0xF6092D, 0x7045D6, 0x7CDC20, 0xFA90DB, 0x65EFCC, 0xE3A337, 0xEF3AC1, 0x69763A,
//     0x578814, 0xD1C4EF, 0xDD5D19, 0x5B11E2, 0xC46EF5, 0x42220E, 0x4EBBF8, 0xC8F703,
//     0x3F964D, 0xB9DAB6, 0xB54340, 0x330FBB, 0xAC70AC, 0x2A3C57, 0x26A5A1, 0xA0E95A,
//     0x9E1774, 0x185B8F, 0x14C279, 0x928E82, 0x0DF195, 0x8BBD6E, 0x872498, 0x016863,
//     0xFAD8C4, 0x7C943F, 0x700DC9, 0xF64132, 0x693E25, 0xEF72DE, 0xE3EB28, 0x65A7D3,
//     0x5B59FD, 0xDD1506, 0xD18CF0, 0x57C00B, 0xC8BF1C, 0x4EF3E7, 0x426A11, 0xC426EA,
//     0x2AE476, 0xACA88D, 0xA0317B, 0x267D80, 0xB90297, 0x3F4E6C, 0x33D79A, 0xB59B61,
//     0x8B654F, 0x0D29B4, 0x01B042, 0x87FCB9, 0x1883AE, 0x9ECF55, 0x9256A3, 0x141A58,
//     0xEFAAFF, 0x69E604, 0x657FF2, 0xE33309, 0x7C4C1E, 0xFA00E5, 0xF69913, 0x70D5E8,
//     0x4E2BC6, 0xC8673D, 0xC4FECB, 0x42B230, 0xDDCD27, 0x5B81DC, 0x57182A, 0xD154D1,
//     0x26359F, 0xA07964, 0xACE092, 0x2AAC69, 0xB5D37E, 0x339F85, 0x3F0673, 0xB94A88,
//     0x87B4A6, 0x01F85D, 0x0D61AB, 0x8B2D50, 0x145247, 0x921EBC, 0x9E874A, 0x18CBB1,
//     0xE37B16, 0x6537ED, 0x69AE1B, 0xEFE2E0, 0x709DF7, 0xF6D10C, 0xFA48FA, 0x7C0401,
//     0x42FA2F, 0xC4B6D4, 0xC82F22, 0x4E63D9, 0xD11CCE, 0x575035, 0x5BC9C3, 0xDD8538};

const unsigned int Crc24q[256] = {
	0x00000000u, 0x01864CFBu, 0x028AD50Du, 0x030C99F6u, 0x0493E6E1u, 0x0515AA1Au, 0x061933ECu, 0x079F7F17u,
	0x08A18139u, 0x0927CDC2u, 0x0A2B5434u, 0x0BAD18CFu, 0x0C3267D8u, 0x0DB42B23u, 0x0EB8B2D5u, 0x0F3EFE2Eu,
	0x10C54E89u, 0x11430272u, 0x124F9B84u, 0x13C9D77Fu, 0x1456A868u, 0x15D0E493u, 0x16DC7D65u, 0x175A319Eu,
	0x1864CFB0u, 0x19E2834Bu, 0x1AEE1ABDu, 0x1B685646u, 0x1CF72951u, 0x1D7165AAu, 0x1E7DFC5Cu, 0x1FFBB0A7u,
	0x200CD1E9u, 0x218A9D12u, 0x228604E4u, 0x2300481Fu, 0x249F3708u, 0x25197BF3u, 0x2615E205u, 0x2793AEFEu,
	0x28AD50D0u, 0x292B1C2Bu, 0x2A2785DDu, 0x2BA1C926u, 0x2C3EB631u, 0x2DB8FACAu, 0x2EB4633Cu, 0x2F322FC7u,
	0x30C99F60u, 0x314FD39Bu, 0x32434A6Du, 0x33C50696u, 0x345A7981u, 0x35DC357Au, 0x36D0AC8Cu, 0x3756E077u,
	0x38681E59u, 0x39EE52A2u, 0x3AE2CB54u, 0x3B6487AFu, 0x3CFBF8B8u, 0x3D7DB443u, 0x3E712DB5u, 0x3FF7614Eu,
	0x4019A3D2u, 0x419FEF29u, 0x429376DFu, 0x43153A24u, 0x448A4533u, 0x450C09C8u, 0x4600903Eu, 0x4786DCC5u,
	0x48B822EBu, 0x493E6E10u, 0x4A32F7E6u, 0x4BB4BB1Du, 0x4C2BC40Au, 0x4DAD88F1u, 0x4EA11107u, 0x4F275DFCu,
	0x50DCED5Bu, 0x515AA1A0u, 0x52563856u, 0x53D074ADu, 0x544F0BBAu, 0x55C94741u, 0x56C5DEB7u, 0x5743924Cu,
	0x587D6C62u, 0x59FB2099u, 0x5AF7B96Fu, 0x5B71F594u, 0x5CEE8A83u, 0x5D68C678u, 0x5E645F8Eu, 0x5FE21375u,
	0x6015723Bu, 0x61933EC0u, 0x629FA736u, 0x6319EBCDu, 0x648694DAu, 0x6500D821u, 0x660C41D7u, 0x678A0D2Cu,
	0x68B4F302u, 0x6932BFF9u, 0x6A3E260Fu, 0x6BB86AF4u, 0x6C2715E3u, 0x6DA15918u, 0x6EADC0EEu, 0x6F2B8C15u,
	0x70D03CB2u, 0x71567049u, 0x725AE9BFu, 0x73DCA544u, 0x7443DA53u, 0x75C596A8u, 0x76C90F5Eu, 0x774F43A5u,
	0x7871BD8Bu, 0x79F7F170u, 0x7AFB6886u, 0x7B7D247Du, 0x7CE25B6Au, 0x7D641791u, 0x7E688E67u, 0x7FEEC29Cu,
	0x803347A4u, 0x81B50B5Fu, 0x82B992A9u, 0x833FDE52u, 0x84A0A145u, 0x8526EDBEu, 0x862A7448u, 0x87AC38B3u,
	0x8892C69Du, 0x89148A66u, 0x8A181390u, 0x8B9E5F6Bu, 0x8C01207Cu, 0x8D876C87u, 0x8E8BF571u, 0x8F0DB98Au,
	0x90F6092Du, 0x917045D6u, 0x927CDC20u, 0x93FA90DBu, 0x9465EFCCu, 0x95E3A337u, 0x96EF3AC1u, 0x9769763Au,
	0x98578814u, 0x99D1C4EFu, 0x9ADD5D19u, 0x9B5B11E2u, 0x9CC46EF5u, 0x9D42220Eu, 0x9E4EBBF8u, 0x9FC8F703u,
	0xA03F964Du, 0xA1B9DAB6u, 0xA2B54340u, 0xA3330FBBu, 0xA4AC70ACu, 0xA52A3C57u, 0xA626A5A1u, 0xA7A0E95Au,
	0xA89E1774u, 0xA9185B8Fu, 0xAA14C279u, 0xAB928E82u, 0xAC0DF195u, 0xAD8BBD6Eu, 0xAE872498u, 0xAF016863u,
	0xB0FAD8C4u, 0xB17C943Fu, 0xB2700DC9u, 0xB3F64132u, 0xB4693E25u, 0xB5EF72DEu, 0xB6E3EB28u, 0xB765A7D3u,
	0xB85B59FDu, 0xB9DD1506u, 0xBAD18CF0u, 0xBB57C00Bu, 0xBCC8BF1Cu, 0xBD4EF3E7u, 0xBE426A11u, 0xBFC426EAu,
	0xC02AE476u, 0xC1ACA88Du, 0xC2A0317Bu, 0xC3267D80u, 0xC4B90297u, 0xC53F4E6Cu, 0xC633D79Au, 0xC7B59B61u,
	0xC88B654Fu, 0xC90D29B4u, 0xCA01B042u, 0xCB87FCB9u, 0xCC1883AEu, 0xCD9ECF55u, 0xCE9256A3u, 0xCF141A58u,
	0xD0EFAAFFu, 0xD169E604u, 0xD2657FF2u, 0xD3E33309u, 0xD47C4C1Eu, 0xD5FA00E5u, 0xD6F69913u, 0xD770D5E8u,
	0xD84E2BC6u, 0xD9C8673Du, 0xDAC4FECBu, 0xDB42B230u, 0xDCDDCD27u, 0xDD5B81DCu, 0xDE57182Au, 0xDFD154D1u,
	0xE026359Fu, 0xE1A07964u, 0xE2ACE092u, 0xE32AAC69u, 0xE4B5D37Eu, 0xE5339F85u, 0xE63F0673u, 0xE7B94A88u,
	0xE887B4A6u, 0xE901F85Du, 0xEA0D61ABu, 0xEB8B2D50u, 0xEC145247u, 0xED921EBCu, 0xEE9E874Au, 0xEF18CBB1u,
	0xF0E37B16u, 0xF16537EDu, 0xF269AE1Bu, 0xF3EFE2E0u, 0xF4709DF7u, 0xF5F6D10Cu, 0xF6FA48FAu, 0xF77C0401u,
	0xF842FA2Fu, 0xF9C4B6D4u, 0xFAC82F22u, 0xFB4E63D9u, 0xFCD11CCEu, 0xFD575035u, 0xFE5BC9C3u, 0xFFDD8538u,
};

#endif
#endif