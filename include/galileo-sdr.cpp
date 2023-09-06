#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include "getopt.h"
#else
#include <unistd.h>
#endif

#include "socket.h"
#include "usrp_galileo.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <ncurses.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <cmath>
#include <signal.h>

int samples_per_code;
std::vector<int> current_eph;
bool stop_signal_called_gal_task = false;
bool use_cursors = true;
bool advance_fptr = false;

FILE *TV_ptr;
FILE *TV_ptrs_sv[MAX_SAT];

void sigint_handler_gal_task(int code)
{
    (void)code;
    endwin();
    stop_signal_called_gal_task = true;
    fprintf(stderr, "Done\n");
}

#define COS_TAB_LENGTH 512
static const int cosTable512[COS_TAB_LENGTH] = {
    250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
    245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
    230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
    207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
    176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
    138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100, 97,
    94, 91, 89, 86, 83, 80, 77, 74, 71, 68, 65, 62, 59, 56, 53, 50,
    47, 44, 41, 38, 35, 32, 29, 26, 23, 20, 17, 14, 11, 8, 5, 2,
    -2, -5, -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
    -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
    -97, -100, -103, -105, -108, -111, -114, -116, -119, -122, -125, -127, -130, -132, -135, -138,
    -140, -143, -145, -148, -150, -153, -155, -157, -160, -162, -164, -167, -169, -171, -173, -176,
    -178, -180, -182, -184, -186, -188, -190, -192, -194, -196, -198, -200, -202, -204, -205, -207,
    -209, -210, -212, -214, -215, -217, -218, -220, -221, -223, -224, -225, -227, -228, -229, -230,
    -232, -233, -234, -235, -236, -237, -238, -239, -240, -241, -241, -242, -243, -244, -244, -245,
    -245, -246, -247, -247, -248, -248, -248, -249, -249, -249, -249, -250, -250, -250, -250, -250,
    -250, -250, -250, -250, -250, -249, -249, -249, -249, -248, -248, -248, -247, -247, -246, -245,
    -245, -244, -244, -243, -242, -241, -241, -240, -239, -238, -237, -236, -235, -234, -233, -232,
    -230, -229, -228, -227, -225, -224, -223, -221, -220, -218, -217, -215, -214, -212, -210, -209,
    -207, -205, -204, -202, -200, -198, -196, -194, -192, -190, -188, -186, -184, -182, -180, -178,
    -176, -173, -171, -169, -167, -164, -162, -160, -157, -155, -153, -150, -148, -145, -143, -140,
    -138, -135, -132, -130, -127, -125, -122, -119, -116, -114, -111, -108, -105, -103, -100, -97,
    -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
    -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11, -8, -5, -2,
    2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47,
    50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 91, 94,
    97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
    140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
    178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
    209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
    232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
    245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250};

static const int sinTable512[COS_TAB_LENGTH] = {
    2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47,
    50, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 91, 94,
    97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
    140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
    178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
    209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
    232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
    245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250,
    250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
    245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
    230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
    207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
    176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
    138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100, 97,
    94, 91, 89, 86, 83, 80, 77, 74, 71, 68, 65, 62, 59, 56, 53, 50,
    47, 44, 41, 38, 35, 32, 29, 26, 23, 20, 17, 14, 11, 8, 5, 2,
    -2, -5, -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
    -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
    -97, -100, -103, -105, -108, -111, -114, -116, -119, -122, -125, -127, -130, -132, -135, -138,
    -140, -143, -145, -148, -150, -153, -155, -157, -160, -162, -164, -167, -169, -171, -173, -176,
    -178, -180, -182, -184, -186, -188, -190, -192, -194, -196, -198, -200, -202, -204, -205, -207,
    -209, -210, -212, -214, -215, -217, -218, -220, -221, -223, -224, -225, -227, -228, -229, -230,
    -232, -233, -234, -235, -236, -237, -238, -239, -240, -241, -241, -242, -243, -244, -244, -245,
    -245, -246, -247, -247, -248, -248, -248, -249, -249, -249, -249, -250, -250, -250, -250, -250,
    -250, -250, -250, -250, -250, -249, -249, -249, -249, -248, -248, -248, -247, -247, -246, -245,
    -245, -244, -244, -243, -242, -241, -241, -240, -239, -238, -237, -236, -235, -234, -233, -232,
    -230, -229, -228, -227, -225, -224, -223, -221, -220, -218, -217, -215, -214, -212, -210, -209,
    -207, -205, -204, -202, -200, -198, -196, -194, -192, -190, -188, -186, -184, -182, -180, -178,
    -176, -173, -171, -169, -167, -164, -162, -160, -157, -155, -153, -150, -148, -145, -143, -140,
    -138, -135, -132, -130, -127, -125, -122, -119, -116, -114, -111, -108, -105, -103, -100, -97,
    -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
    -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11, -8, -5, -2};

// Receiver antenna attenuation in dB for boresight angle = 0:5:180 [deg]
double ant_pat_db[37] = {0.00, 0.00, 0.22, 0.44, 0.67, 1.11, 1.56, 2.00,
                         2.44, 2.89, 3.56, 4.22, 4.89, 5.56, 6.22, 6.89,
                         7.56, 8.22, 8.89, 9.78, 10.67, 11.56, 12.44, 13.33,
                         14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67,
                         24.00, 25.56, 27.33, 29.33, 31.56};

int allocatedSat[MAX_SAT];

/*! \brief Subtract two vectors of double
 *  \param[out] y Result of subtraction
 *  \param[in] x1 Minuend of subtraction
 *  \param[in] x2 Subtrahend of subtraction
 */
void subVect(double *y, const double *x1, const double *x2)
{
    y[0] = x1[0] - x2[0];
    y[1] = x1[1] - x2[1];
    y[2] = x1[2] - x2[2];

    return;
}

/*! \brief Compute Norm of Vector
 *  \param[in] x Input vector
 *  \returns Length (Norm) of the input vector
 */
double normVect(const double *x)
{
    return (sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
}

/*! \brief Used to debug - remove */
int writecsv(int prn, double tow, double values[3], int type)
{
    FILE *csvFile = fopen("data.csv", "a");

    if (csvFile == NULL)
    {
        perror("Error opening file");
        return -1;
    }

    // Write the values to the CSV file
    switch (type)
    {
    case 0:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, satpos", prn, values[0], values[1], values[2], tow);
        break;

    case 1:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, los", prn, values[0], values[1], values[2], tow);
        break;

    case 2:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, rxpos", prn, values[0], values[1], values[2], tow);
        break;

    case 3:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, range", prn, values[0], values[1], values[2], tow);
        break;
    }

    // Add a newline to separate rows
    fprintf(csvFile, "\n");

    // Close the CSV file
    fclose(csvFile);

    return 1;
}

/*! \brief Compute dot-product of two vectors
 *  \param[in] x1 First multiplicand
 *  \param[in] x2 Second multiplicand
 *  \returns Dot-product of both multiplicands
 */
double dotProd(const double *x1, const double *x2)
{
    return (x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
}

// Converts hex E1 PRN code to binary (-1,1)
void hex_to_binary_converter(short *dest, bool c1, int prn)
{
    int index = 0;

    for (int i = 0; i < 1024; i++)
    {
        char from;
        if (!c1)
        {
            from = GALILEO_E1_B_PRIMARY_CODE[prn][i];
        }
        else
        {
            from = GALILEO_E1_C_PRIMARY_CODE[prn][i];
        }

        switch (from)
        {
        case '0':
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case '1':
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case '2':
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case '3':
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case '4':
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case '5':
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case '6':
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case '7':
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case '8':
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case '9':
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case 'A':
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case 'B':
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case 'C':
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case 'D':
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            dest[index] = -1;
            index++;
            break;
        case 'E':
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = 1;
            index++;
            break;
        case 'F':
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            dest[index] = -1;
            index++;
            break;
        default:
            break;
        }
    }
}

void sboc(short *dest, short *prn, int len, int m, int n)
{
    constexpr uint32_t length_in = CA_SEQ_LEN_E1;
    const auto period = static_cast<uint32_t>(CA_SEQ_LEN_E1 / length_in);

    int i, j, N = 2 * m / n;
    for (i = 0; i < len; i++)
        for (j = 0; j < N; j++)
            dest[N * i + j] = prn[i];

    // Mix sub carrier
    for (i = 0; i < N * len / 2; i++)
    {
        dest[2 * i] = -dest[2 * i];
    }

    // for (uint32_t i = 0; i < length_in; i++)
    // {
    //     for (uint32_t j = 0; j < (period / 2); j++)
    //     {
    //         dest[i * period + j] = prn[i];
    //     }
    //     for (uint32_t j = (period / 2); j < period; j++)
    //     {
    //         dest[i * period + j] = -prn[i];
    //     }
    // }
}

/* !\brief generate the BOC E1B code sequence for a given Satellite Vehicle PRN
 *  \param[in] prn PRN number of the Satellite Vehicle
 *  \param[out] ca Caller-allocated integer array of 1023 bytes
 */
void codegen_E1B(short *ca, int prn)
{
    // TODO - figure out the delay for each satellite? I think there is a standard
    // somefor, for now keep it at 0 int delay[] = { 	  5,   6,   7,   8,  17,  18,
    // 139, 140, 141, 251, 	252, 254, 255, 256, 257, 258, 469, 470, 471, 472, 	473,
    // 474, 509, 512, 513, 514, 515, 516, 859, 860, 	861, 862};

    // Get the PRN code using the hex to binary converter
    short *tmp_ca = (short *)malloc(CA_SEQ_LEN_E1 * sizeof(short));
    hex_to_binary_converter(tmp_ca, false, prn - 1);
    sboc(ca, tmp_ca, CA_SEQ_LEN_E1, 1, 1);

    // for (int i=0; i<CA_SEQ_LEN_E1*2; i++)
    //     fprintf(stderr, "%d", ca[i]);

    // exit(1);
}

void codegen_E1C(short *ca, int prn)
{
    // TODO - figure out the delay for each satellite? I think there is a standard
    // somefor, for now keep it at 0 int delay[] = { 	  5,   6,   7,   8,  17,  18,
    // 139, 140, 141, 251, 	252, 254, 255, 256, 257, 258, 469, 470, 471, 472, 	473,
    // 474, 509, 512, 513, 514, 515, 516, 859, 860, 	861, 862};

    // Get the PRN code using the hex to binary converter
    short *tmp_ca = (short *)malloc(CA_SEQ_LEN_E1 * sizeof(short));
    hex_to_binary_converter(tmp_ca, true, prn - 1);
    sboc(ca, tmp_ca, CA_SEQ_LEN_E1, 1, 1);
}

/*! \brief Convert a UTC date into a GAL date
 *  \param[in] t input date in UTC form
 *  \param[out] g output date GAL form
 */
void date2gal(const datetime_t *t, galtime_t *g)
{
    int doy[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int ye;
    int de;
    int lpdays;

    ye = t->y - 1980;

    // Compute the number of leap days since Jan 5/Jan 6, 1980.
    lpdays = ye / 4 + 1;
    if ((ye % 4) == 0 && t->m <= 2)
        lpdays--;

    // Compute the number of days elapsed since Jan 5/Jan 6, 1980.
    de = ye * 365 + doy[t->m - 1] + t->d + lpdays - 6;

    // Convert time to GPS weeks and seconds.
    g->week = de / 7;
    g->sec = (double)(de % 7) * SECONDS_IN_DAY + t->hh * SECONDS_IN_HOUR + t->mm * SECONDS_IN_MINUTE + t->sec;

    return;
}

void gal2date(const galtime_t *g, datetime_t *t) // checked
{
    // Convert Julian day number to calendar date
    int c = (int)(7 * g->week + floor(g->sec / 86400.0) + 2444245.0) + 1537;
    int d = (int)((c - 122.1) / 365.25);
    int e = 365 * d + d / 4;
    int f = (int)((c - e) / 30.6001);

    t->d = c - e - (int)(30.6001 * f);
    t->m = f - 1 - 12 * (f / 14);
    t->y = d - 4715 - ((7 + t->m) / 10);

    t->hh = ((int)(g->sec / 3600.0)) % 24;
    t->mm = ((int)(g->sec / 60.0)) % 60;
    t->sec = g->sec - 60.0 * floor(g->sec / 60.0);

    return;
}

/*! \brief Convert Earth-centered Earth-fixed (ECEF) into Lat/Long/Height
 *  \param[in] xyz Input Array of X, Y and Z ECEF coordinates
 *  \param[out] llh Output Array of Latitude, Longitude and Height
 */
void xyz2llh(const double *xyz, double *llh)
{
    double a, eps, e, e2;
    double x, y, z;
    double rho2, dz, zdz, nh, slat, n, dz_new;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;

    eps = 1.0e-3;
    e2 = e * e;

    if (normVect(xyz) < eps)
    {
        // Invalid ECEF vector
        llh[0] = 0.0;
        llh[1] = 0.0;
        llh[2] = -a;

        return;
    }

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];

    rho2 = x * x + y * y;
    dz = e2 * z;

    while (1)
    {
        zdz = z + dz;
        nh = sqrt(rho2 + zdz * zdz);
        slat = zdz / nh;
        n = a / sqrt(1.0 - e2 * slat * slat);
        dz_new = n * e2 * slat;

        if (fabs(dz - dz_new) < eps)
            break;

        dz = dz_new;
    }

    llh[0] = atan2(zdz, sqrt(rho2));
    llh[1] = atan2(y, x);
    llh[2] = nh - n;

    return;
}

/*! \brief Convert Lat/Long/Height into Earth-centered Earth-fixed (ECEF)
 *  \param[in] llh Input Array of Latitude, Longitude and Height
 *  \param[out] xyz Output Array of X, Y and Z ECEF coordinates
 */
void llh2xyz(const double *llh, double *xyz)
{
    double n;
    double a;
    double e;
    double e2;
    double clat;
    double slat;
    double clon;
    double slon;
    double d, nph;
    double tmp;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;
    e2 = e * e;

    clat = cos(llh[0]);
    slat = sin(llh[0]);
    clon = cos(llh[1]);
    slon = sin(llh[1]);
    d = e * slat;

    n = a / sqrt(1.0 - d * d);
    nph = n + llh[2];

    tmp = nph * clat;
    xyz[0] = tmp * clon;
    xyz[1] = tmp * slon;
    xyz[2] = ((1.0 - e2) * n + llh[2]) * slat;

    return;
}

/*! \brief Compute the intermediate matrix for LLH to ECEF
 *  \param[in] llh Input position in Latitude-Longitude-Height format
 *  \param[out] t Three-by-Three output matrix
 */
void ltcmat(const double *llh, double t[3][3])
{
    double slat, clat;
    double slon, clon;

    slat = sin(llh[0]);
    clat = cos(llh[0]);
    slon = sin(llh[1]);
    clon = cos(llh[1]);

    t[0][0] = -slat * clon;
    t[0][1] = -slat * slon;
    t[0][2] = clat;
    t[1][0] = -slon;
    t[1][1] = clon;
    t[1][2] = 0.0;
    t[2][0] = clat * clon;
    t[2][1] = clat * slon;
    t[2][2] = slat;

    return;
}

/*! \brief Convert Earth-centered Earth-Fixed to ?
 *  \param[in] xyz Input position as vector in ECEF format
 *  \param[in] t Intermediate matrix computed by \ref ltcmat
 *  \param[out] neu Output position as North-East-Up format
 */
void ecef2neu(const double *xyz, double t[3][3], double *neu)
{
    neu[0] = t[0][0] * xyz[0] + t[0][1] * xyz[1] + t[0][2] * xyz[2];
    neu[1] = t[1][0] * xyz[0] + t[1][1] * xyz[1] + t[1][2] * xyz[2];
    neu[2] = t[2][0] * xyz[0] + t[2][1] * xyz[1] + t[2][2] * xyz[2];

    return;
}

/*! \brief Convert North-East-Up to Azimuth + Elevation
 *  \param[in] neu Input position in North-East-Up format
 *  \param[out] azel Output array of azimuth + elevation as double
 */
void neu2azel(double *azel, const double *neu)
{
    double ne;

    azel[0] = atan2(neu[1], neu[0]);
    if (azel[0] < 0.0)
        azel[0] += (2.0 * PI);

    ne = sqrt(neu[0] * neu[0] + neu[1] * neu[1]);
    azel[1] = atan2(neu[2], ne);

    return;
}

/*! \brief Compute Satellite position, velocity and clock at given time
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at which position is to be computed
 *  \param[out] pos Computed position (vector)
 *  \param[out] vel Computed velocity (vector)
 *  \param[clk] clk Computed clock
 */
void satpos(ephem_t eph, galtime_t g, double *pos, double *vel, double *clk)
{
    // Computing Satellite Velocity using the Broadcast Ephemeris
    // http://www.ngs.noaa.gov/gps-toolbox/bc_velo.htm

    double tk;
    double mk;
    double ek;
    double ekold;
    double ekdot;
    double cek, sek;
    double pk;
    double pkdot;
    double c2pk, s2pk;
    double uk;
    double ukdot;
    double cuk, suk;
    double ok;
    double sok, cok;
    double ik;
    double ikdot;
    double sik, cik;
    double rk;
    double rkdot;
    double xpk, ypk;
    double xpkdot, ypkdot;

    double relativistic, OneMinusecosE = 0, tmp;

    tk = g.sec - eph.toe.sec;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    mk = eph.m0 + eph.n * tk;
    ek = mk;
    ekold = ek + 1.0;

    int cpt = 0;

    while ((fabs(ek - ekold) > 1.0E-14) && cpt < 500)
    {

        cpt++;
        ekold = ek;
        OneMinusecosE = 1.0 - eph.ecc * cos(ekold);
        ek = ek + (mk - ekold + eph.ecc * sin(ekold)) / OneMinusecosE;
    }

    sek = sin(ek);
    cek = cos(ek);

    ekdot = eph.n / OneMinusecosE;

    relativistic = -4.442807633E-10 * eph.ecc * eph.sqrta * sek;

    pk = atan2(eph.sq1e2 * sek, cek - eph.ecc) + eph.aop;
    pkdot = eph.sq1e2 * ekdot / OneMinusecosE;

    s2pk = sin(2.0 * pk);
    c2pk = cos(2.0 * pk);

    uk = pk + eph.cus * s2pk + eph.cuc * c2pk;
    suk = sin(uk);
    cuk = cos(uk);
    ukdot = pkdot * (1.0 + 2.0 * (eph.cus * c2pk - eph.cuc * s2pk));

    rk = eph.A * OneMinusecosE + eph.crc * c2pk + eph.crs * s2pk;
    rkdot = eph.A * eph.ecc * sek * ekdot + 2.0 * pkdot * (eph.crs * c2pk - eph.crc * s2pk);

    ik = eph.inc0 + eph.idot * tk + eph.cic * c2pk + eph.cis * s2pk;

    sik = sin(ik);
    cik = cos(ik);
    ikdot = eph.idot + 2.0 * pkdot * (eph.cis * c2pk - eph.cic * s2pk);

    xpk = rk * cuk;
    ypk = rk * suk;
    xpkdot = rkdot * cuk - ypk * ukdot;
    ypkdot = rkdot * suk + xpk * ukdot;

    ok = eph.omg0 + tk * eph.omgkdot - OMEGA_EARTH * eph.toe.sec;
    sok = sin(ok);
    cok = cos(ok);

    // Compute position
    pos[0] = xpk * cok - ypk * cik * sok;
    pos[1] = xpk * sok + ypk * cik * cok;
    pos[2] = ypk * sik;

    tmp = ypkdot * cik - ypk * sik * ikdot;

    // Compute velocity
    vel[0] = -eph.omgkdot * pos[1] + xpkdot * cok - tmp * sok;
    vel[1] = eph.omgkdot * pos[0] + xpkdot * sok + tmp * cok;
    vel[2] = ypk * cik * ikdot + ypkdot * sik;

    // Satellite clock correction
    tk = g.sec - eph.toc.sec;

    if (tk > SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if (tk < -SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    clk[0] = eph.af0 + tk * (eph.af1 + tk * eph.af2) + relativistic - eph.bgde5b;

    clk[1] = eph.af1 + 2.0 * tk * eph.af2;

    //*dts = eph->f0 + eph->f1 * tk + eph->f2 * tk * tk;

    // clk[1] = eph.af1 + 2.0 * tk * eph.af2;

    // char s[25];

    // sprintf(s, "/tmp/sim_pos_debug_%d", eph.PRN);

    // FILE *fp = fopen(s,"a");
    // fwrite(pos, sizeof(double), 3, fp);
    // fwrite(&g.sec, sizeof(double), 1, fp);
    // fclose(fp);

    return;
}

double subGalTime(galtime_t g1, galtime_t g0)
{
    double dt;

    dt = g1.sec - g0.sec;
    dt += (double)(g1.week - g0.week) * SECONDS_IN_WEEK;

    return (dt);
}

galtime_t incGalTime(galtime_t g, double dt)
{
    g.sec = g.sec + dt;
    return g;
}

int readRinexNavAll(ephem_t eph[][MAX_SAT], ionoutc_t *ionoutc,
                    const char *fname) // checked
{
    FILE *fp;
    int ieph;

    int sv;
    char str[MAX_CHAR];
    char tmp[30];

    datetime_t t;
    galtime_t g;
    galtime_t g0;
    double dt;

    int flags = 0x0;

    if (NULL == (fp = fopen(fname, "rt")))
        return (-1);

    // Clear valid flag
    for (ieph = 0; ieph < EPHEM_ARRAY_SIZE; ieph++)
        for (sv = 0; sv < MAX_SAT; sv++)
            eph[ieph][sv].vflg = 0;

    // Read header lines
    while (1)
    {
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        if (strncmp(str + 60, "END OF HEADER", 13) == 0)
            break;
        else if (strncmp(str + 60, "IONOSPHERIC CORR", 16) == 0)
        {
            strncpy(tmp, str + 6, 11);
            tmp[11] = 0;
            ionoutc->ai0 = atof(tmp);
            // cout << "ion_corr0"  << ionoutc->ion_corr0 << endl;

            strncpy(tmp, str + 18, 11);
            tmp[11] = 0;
            ionoutc->ai1 = atof(tmp);
            // cout << ionoutc->ion_corr1 << endl;

            strncpy(tmp, str + 31, 11);
            tmp[11] = 0;
            ionoutc->ai2 = atof(tmp);
            // cout << ionoutc->ion_corr2 << endl;

            strncpy(tmp, str + 43, 11);
            tmp[11] = 0;
            ionoutc->ai3 = atof(tmp);
            // cout << ionoutc->ion_corr3 << endl;

            // flags |= 0x1;
        }
        else if (strncmp(str + 60, "TIME SYSTEM CORR", 16) == 0)
        {
            strncpy(tmp, str + 7, 17);
            tmp[17] = 0;
            // cout << tmp << endl;
            ionoutc->A0 = atof(tmp);
            // cout << ionoutc->time_corr0 << endl;

            strncpy(tmp, str + 24, 15);
            tmp[15] = 0;
            // cout << tmp << endl;
            ionoutc->A1 = atof(tmp);
            // cout << ionoutc->time_corr1 << endl;

            strncpy(tmp, str + 39, 6);
            tmp[6] = 0;
            // cout << tmp << endl;
            ionoutc->tot = atof(tmp);
            // cout << ionoutc->time_corr2 << endl;

            strncpy(tmp, str + 46, 4);
            tmp[4] = 0;
            // cout << tmp << endl;
            ionoutc->wnt = atof(tmp);
            // cout << ionoutc->time_corr3 << endl;

            // flags |= 0x1 << 1;
        }

        else if (strncmp(str + 60, "LEAP SECONDS", 12) == 0)
        {
            strncpy(tmp, str, 6);
            tmp[6] = 0;
            ionoutc->dtls = atoi(tmp);
            // cout << ionoutc->dtls << endl;

            // flags |= 0x1 << 3;
        }
    }

    // ionoutc->vflg = false;
    // if (flags == 0xF) // Read all Iono/UTC lines
    // 	ionoutc->vflg = true;
    ionoutc->vflg = TRUE;
    // Read ephemeris blocks
    g0.week = -1;
    ieph = 0;

    std::vector<int> v;
    int in = 0;

    while (1)
    // for (int i = 1; i <= 20; ++i)
    {
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        // PRN
        strncpy(tmp, str + 1, 2);
        tmp[2] = 0;
        sv = atoi(tmp) - 1;
        // std::cout << "sv : " << sv << std::endl;

        // EPOCH
        strncpy(tmp, str + 4, 4);
        tmp[4] = 0;
        t.y = atoi(tmp);
        // cout << "y : " << t.y << endl;

        strncpy(tmp, str + 9, 2);
        tmp[2] = 0;
        t.m = atoi(tmp);
        // cout << "m : " << t.m << endl;

        strncpy(tmp, str + 12, 2);
        tmp[2] = 0;
        t.d = atoi(tmp);
        // std::cout << "d : " << t.d << std::endl;

        strncpy(tmp, str + 15, 2);
        tmp[2] = 0;
        t.hh = atoi(tmp);

        strncpy(tmp, str + 18, 2);
        tmp[2] = 0;
        t.mm = atoi(tmp);
        // cout << "mm : " << t.mm << endl;

        strncpy(tmp, str + 21, 2);
        tmp[2] = 0;
        t.sec = atof(tmp);
        // cout << "sec : " << t.sec << endl;

        date2gal(&t, &g);

        if (g0.week == -1)
            g0 = g;

        // Check current time of clock
        dt = subGalTime(g, g0);

        if (dt > SECONDS_IN_HOUR)
        {
            g0 = g;
            ieph++; // a new set of ephemerides

            if (ieph >= EPHEM_ARRAY_SIZE)
                break;
        }

        // g0 = g;
        // ieph++;

        // Date and time
        eph[ieph][sv].t = t;

        // SV CLK
        eph[ieph][sv].toc = g;

        // 0
        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].af0 = atof(tmp);
        // cout << "af0 : " << eph[ieph][sv].af0 << endl;

        // 1
        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].af1 = atof(tmp);
        // cout << "af1 : " << eph[ieph][sv].af1 << endl;

        // 2
        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].af2 = atof(tmp);
        // cout << "af2 : " << eph[ieph][sv].af2 << endl;

        // BROADCAST ORBIT - 1 // 3
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].iode = (int)atof(tmp);
        // cout << "iode : " << eph[ieph][sv].iode << endl;

        // 4
        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].crs = atof(tmp);
        // cout << "crs : " << eph[ieph][sv].crs << endl;

        // 5
        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].deltan = atof(tmp);
        // cout << "deltan : " << eph[ieph][sv].deltan << endl;

        // 6
        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].m0 = atof(tmp);
        // cout << "m0 : " << eph[ieph][sv].m0 << endl;

        // BROADCAST ORBIT - 2
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].cuc = atof(tmp);

        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].ecc = atof(tmp);

        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].cus = atof(tmp);

        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].sqrta = atof(tmp);

        // BROADCAST ORBIT - 3
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].toe.sec = atof(tmp);
        // eph[ieph][sv].toc.sec;

        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].cic = atof(tmp);

        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].omg0 = atof(tmp);

        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].cis = atof(tmp);

        // BROADCAST ORBIT - 4
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].inc0 = atof(tmp);

        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].crc = atof(tmp);

        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].aop = atof(tmp);

        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].omgdot = atof(tmp);
        // cout << "orbit 4 :  " << eph[ieph][sv].omgdot << endl;

        // BROADCAST ORBIT - 5
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].idot = atof(tmp);

        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].codeL2 = (int)atof(tmp);

        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].toe.week = (int)atof(tmp);

        // BROADCAST ORBIT - 6
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 23, 19);
        tmp[19] = 0;
        eph[ieph][sv].svhlth = (int)atof(tmp);
        if ((eph[ieph][sv].svhlth > 0) && (eph[ieph][sv].svhlth < 32))
            eph[ieph][sv].svhlth += 32; // Set MSB to 1

        strncpy(tmp, str + 42, 19);
        tmp[19] = 0;
        eph[ieph][sv].bgde5a = atof(tmp);

        strncpy(tmp, str + 61, 19);
        tmp[19] = 0;
        eph[ieph][sv].bgde5b = (int)atof(tmp);

        // BROADCAST ORBIT - 7
        if (NULL == fgets(str, MAX_CHAR, fp))
            break;

        strncpy(tmp, str + 4, 19);
        tmp[19] = 0;
        eph[ieph][sv].toc.sec = atof(tmp);

        // Set valid flag
        eph[ieph][sv].vflg = 1;

        if (!(find(v.begin(), v.end(), sv) != v.end()))
        {
            v.push_back(sv);
        }

        // Update the working variables
        eph[ieph][sv].A = eph[ieph][sv].sqrta * eph[ieph][sv].sqrta;
        eph[ieph][sv].n =
            sqrt(GM_EARTH / (eph[ieph][sv].A * eph[ieph][sv].A * eph[ieph][sv].A)) +
            eph[ieph][sv].deltan;
        eph[ieph][sv].sq1e2 = sqrt(1.0 - eph[ieph][sv].ecc * eph[ieph][sv].ecc);
        eph[ieph][sv].omgkdot = eph[ieph][sv].omgdot - OMEGA_EARTH;
    }

    fclose(fp);

    if (g0.week >= 0)
        ieph += 1; // Number of sets of ephemerides

    return (ieph);
}

double ionosphericDelay(const ionoutc_t *ionoutc, galtime_t g, double *llh, double *azel)
{
    double iono_delay = 0.0;
    double E, phi_u, lam_u, F;

    if (ionoutc->enable == FALSE)
        return (0.0); // No ionospheric delay

    E = azel[1] / PI;
    phi_u = llh[0] / PI;
    lam_u = llh[1] / PI;

    // Obliquity factor
    F = 1.0 + 16.0 * pow((0.53 - E), 3.0);

    if (ionoutc->vflg == FALSE)
        iono_delay = F * 5.0e-9 * SPEED_OF_LIGHT;
    else
    {
        double t, psi, phi_i, lam_i, phi_m, phi_m2, phi_m3;
        double AMP, PER, X, X2, X4;

        // Earth's central angle between the user position and the earth projection
        // of ionospheric intersection point (semi-circles)
        psi = 0.0137 / (E + 0.11) - 0.022;

        // Geodetic latitude of the earth projection of the ionospheric intersection
        // point (semi-circles)
        phi_i = phi_u + psi * cos(azel[0]);
        if (phi_i > 0.416)
            phi_i = 0.416;
        else if (phi_i < -0.416)
            phi_i = -0.416;

        // Geodetic longitude of the earth projection of the ionospheric
        // intersection point (semi-circles)
        lam_i = lam_u + psi * sin(azel[0]) / cos(phi_i * PI);

        // Geomagnetic latitude of the earth projection of the ionospheric
        // intersection point (mean ionospheric height assumed 350 km)
        // (semi-circles)
        phi_m = phi_i + 0.064 * cos((lam_i - 1.617) * PI);
        phi_m2 = phi_m * phi_m;
        phi_m3 = phi_m2 * phi_m;

        AMP = ionoutc->ai0 + ionoutc->ai1 * phi_m + ionoutc->ai2 * phi_m2 +
              ionoutc->ai3 * phi_m3;
        if (AMP < 0.0)
            AMP = 0.0;

        // PER = ionoutc->beta0 + ionoutc->beta1 * phi_m + ionoutc->beta2 * phi_m2 +
        // ionoutc->beta3 * phi_m3; if (PER < 72000.0) 	PER = 72000.0;

        PER = 72000.0;

        // Local time (sec)
        t = SECONDS_IN_DAY / 2.0 * lam_i + g.sec;
        while (t >= SECONDS_IN_DAY)
            t -= SECONDS_IN_DAY;
        while (t < 0)
            t += SECONDS_IN_DAY;

        // Phase (radians)
        X = 2.0 * PI * (t - 50400.0) / PER;

        if (fabs(X) < 1.57)
        {
            X2 = X * X;
            X4 = X2 * X2;
            iono_delay =
                F * (5.0e-9 + AMP * (1.0 - X2 / 2.0 + X4 / 24.0)) * SPEED_OF_LIGHT;
        }
        else
            iono_delay = F * 5.0e-9 * SPEED_OF_LIGHT;
    }

    return (iono_delay);
}

/*! \brief Compute range between a satellite and the receiver
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] xyz position of the receiver
 */
void computeRange(range_t *rho, ephem_t eph, ionoutc_t *ionoutc, galtime_t g, double xyz[], int prn)
{
    double pos[3] = {0.0}, vel[3] = {0.0}, clk[2] = {0.0};
    double los[3] = {0.0};
    double tau = 0.0;
    double range = 0.0, rate = 0.0;
    double xrot = 0.0, yrot = 0.0;

    double llh[3] = {0.0}, neu[3] = {0.0};
    double tmat[3][3] = {0};

    // SV position at time of the pseudorange observation.
    satpos(eph, g, pos, vel, clk);

    // Receiver to satellite vector and light-time.
    subVect(los, pos, xyz);
    tau = normVect(los) / SPEED_OF_LIGHT;

    // Extrapolate the satellite position backwards to the transmission time.
    pos[0] -= vel[0] * tau;
    pos[1] -= vel[1] * tau;
    pos[2] -= vel[2] * tau;

    // Earth rotation correction. The change in velocity can be neglected.
    xrot = pos[0] + pos[1] * GNSS_OMEGA_EARTH_DOT * tau;
    yrot = pos[1] - pos[0] * GNSS_OMEGA_EARTH_DOT * tau;
    pos[0] = xrot;
    pos[1] = yrot;

    // New observer to satellite vector and satellite range.
    subVect(los, pos, xyz);
    range = normVect(los);

    // range = sqrt(std::pow((xyz[0] - pos[0]), 2) + std::pow((xyz[1] - pos[1]),
    // 2) + std::pow((xyz[2] - pos[2]), 2));

    rho->d = range;

    // Pseudorange.
    rho->range = range - SPEED_OF_LIGHT * clk[0];

    double r[3] = {range, range - SPEED_OF_LIGHT * clk[0], 0};

    // Azimuth and elevation angles.
    xyz2llh(xyz, llh);
    ltcmat(llh, tmat);
    ecef2neu(los, tmat, neu);
    neu2azel(rho->azel, neu);

    return;
}

/*! \brief Compute the code phase for a given channel (satellite)
 *  \param chan Channel on which we operate (is updated)
 *  \param[in] rho1 Current range, after \a dt has expired
 *  \param[in dt delta-t (time difference) in seconds
 */
void computeCodePhase(channel_t *chan, range_t rho1, double dt, galtime_t grx) // checked
{
    double ms;
    int ims;
    double rhorate;

    // Pseudorange rate.
    rhorate = (rho1.range - chan->rho0.range) / dt;

    // Carrier and code frequency.
    chan->f_carr = (-rhorate / LAMBDA_E1); // + GALILEO_E1_SUB_CARRIER_A_RATE_HZ;

    chan->f_code = CODE_FREQ_E1 + chan->f_carr * CARR_TO_CODE_E1;

    ms = (grx.sec - rho1.range / SPEED_OF_LIGHT) * 1000.0;

    int ipage = ms / 2000.0; // 1 word = 250 symbols = 1000ms = 1s

    ms -= ipage * 2000;

    int ibit = (unsigned int)ms / 4; // 1 symbol = 1 code = 4ms
    ms -= ibit * 4;
    double code_phase = ms / 4 * CA_SEQ_LEN_E1;

    ms -= ipage * PAGE_TRANS_TIME_ms;

    ibit = (ibit + (N_SYM_PAGE / 2)) % N_SYM_PAGE;

    chan->code_phase = code_phase;

    chan->ibit = ibit;         // gal: 1 bit = 4 ms
    chan->ipage = ipage % 360; // Number of half pages in a frame
    chan->icode = ims;         // 1 code = 4 ms

    // Save current pseudorange
    chan->g0 = grx;

    chan->rho0 = rho1;
    return;
}

int checkSatVisibility(ephem_t eph, galtime_t g, double *xyz, double elvMask, double *azel, int prn) // checked
{
    double llh[3], neu[3];
    double pos[3], vel[3], clk[3], los[3];
    double tmat[3][3];

    if (eph.vflg != 1)
    {
        return (-1); // Invalid}
    }

    xyz2llh(xyz, llh);
    ltcmat(llh, tmat);

    satpos(eph, g, pos, vel, clk);

    subVect(los, pos, xyz);
    ecef2neu(los, tmat, neu);
    neu2azel(azel, neu);

    if (azel[1] * R2D > elvMask)
        return (1); // Visible
    // else

    return (0); // Invisible
}

/* Apply FEC, Interleaving and add Sync pattern (gal-osnma-sdr)*/
void generateFrame(int *half_page, int *frame)
{
    int fec_tmp_sbf[240] = {0};
    // FEC Encoding
    cnv_encd(N_BIT_PAGE, half_page, fec_tmp_sbf);

    int int_tmp_sbf[240] = {0};
    // Interleaving
    for (int r = 0; r < 8; r++)
    {
        for (int c = 0; c < 30; c++)
        {
            int_tmp_sbf[r * 30 + c] = fec_tmp_sbf[c * 8 + r];
        }
    }

    // Add sync
    const int sync_pattern[10] = {0, 1, 0, 1, 1, 0, 0, 0, 0, 0};

    for (int i = 0; i < 10; i++)
        frame[i] = sync_pattern[i];
    for (int i = 0; i < 240; i++)
        frame[i + 10] = int_tmp_sbf[i];
}

#define BIT_ISSET(var, n) !!((uint32_t)(var) & (1 << (n)))

int generateNavMsg(galtime_t g, channel_t *chan, bool advance_fptr)
{
    int iwrd, isbf;
    galtime_t g0;
    unsigned long wn, tow;
    unsigned sbfwrd;
    unsigned long prevwrd;
    int nib, j;
    int *page = (int *)malloc(PAGE_SIZE * sizeof(int));
    int *odd_page = (int *)malloc(120 * sizeof(int));
    int *even_page = (int *)malloc(120 * sizeof(int));
    int index = 0;

    g0.week = g.week;
    g0.sec = (double)(((unsigned long)(g.sec))) + 2; // Align with the full page length = 2 sec

    wn = (unsigned long)(g0.week);
    tow = ((unsigned long)g0.sec);

    // Get page for the current satellite from the
    vector<uint8_t> nav_message = load_page(TV_ptrs_sv[chan->prn - 1], chan->prn, advance_fptr, tow, &g0);

    /* Separate even and odd parts */
    for (int i = 0; i < 15; i++)
    {
        for (j = 7; j >= 0; j--)
        {
            even_page[index] = BIT_ISSET(nav_message[i], j);
            odd_page[index] = BIT_ISSET(nav_message[i + 15], j);
            index++;
        }
    }

    generateFrame(even_page, &page[0]);
    generateFrame(odd_page, &page[250]);

    chan->page = page;

    return (1);
}

int allocateChannel(channel_t *chan,
                    // map<int, vector<Rinex3Nav::DataGAL>> *navGAL,
                    vector<ephem_t> *eph_vector,
                    ionoutc_t ionoutc, galtime_t grx, double *xyz,
                    double elvMask, map<int, int> *sm)
{
    int nsat = 0;
    int i, sv;
    double azel[2];

    range_t rho;
    double ref[3] = {0.0};
    double r_ref, r_xyz;
    double phase_ini;

    ephem_t eph;

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (!eph_vector[sv][0].vflg)
            continue;

        current_eph[sv] = epoch_matcher(grx.sec, eph_vector[sv], current_eph[sv]);

        eph = eph_vector[sv][current_eph[sv]];

        if (checkSatVisibility(eph, grx, xyz, 0, azel, sv + 1) == 1)
        {
            nsat++; // Number of visible satellites

            if (allocatedSat[sv] == -1) // Visible but not allocated
            {
                // Allocated new satellite
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn == 0)
                    {
                        // Initialize channel
                        chan[i].prn = sv + 1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];
                        // chan[i].g0 = grx;

                        // Insert latest channel assignment to the map
                        sm->insert({chan[i].prn, i});

                        // C/A code generation
                        codegen_E1B(chan[i].ca_E1B, chan[i].prn);
                        codegen_E1C(chan[i].ca_E1C, chan[i].prn);

                        // Generate navigation message
                        generateNavMsg(grx, &chan[i], advance_fptr);

                        // Initialize pseudorange
                        computeRange(&rho, eph, &ionoutc, grx, xyz, chan[i].prn);
                        chan[i].rho0 = rho;

                        // Initialize carrier phase
                        r_xyz = rho.range;

                        computeRange(&rho, eph, &ionoutc, grx, ref, chan[i].prn);
                        r_ref = rho.range;

                        phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_L1;
                        chan[i].carr_phase = phase_ini - floor(phase_ini);

                        break;
                    }
                }

                // Set satellite allocation channel
                if (i < MAX_CHAN)
                    allocatedSat[sv] = i;
            }
        }
        else if (allocatedSat[sv] >= 0) // Not visible but allocated
        {
            // Clear channel
            chan[allocatedSat[sv]].prn = 0;

            // Clear satellite allocation flag
            allocatedSat[sv] = -1;
        }
    }
    // advance_fptr = true;
    return (nsat);
}

static long get_nanos(void)
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}

void *galileo_task(void *arg)
{
    struct hash_queue hq;

    vector<queue<int>> queues;
    map<int, int> chn_prn_map;
    int queue_index;

    // hq.queue_bits = queues;
    // hq.chn_prn_map = chn_prn_map;

    sim_t *s = (sim_t *)arg;

    // queue<int> q2;

    // for (int i = 0; i < MAX_SAT; i++)
    //     queues.push_back(queue<int>());

    // pthread_t th_bits;
    bool exit_flag = false;

    // Stores current data bit from all channels
    int dataBit[MAX_CHAN];

    // Values for TOW fix - these determine whether corrections received over
    // socket are applied
    bool tow_fixed;
    bool local_fix = false;
    bool fixed_prn = false;

    clock_t tstart, tend;

    FILE *fp;

    int sv;
    int neph, ieph;
    ephem_t eph1[EPHEM_ARRAY_SIZE][MAX_SAT];
    galtime_t g0;

    double llh[3];

    int i;
    channel_t chan[MAX_CHAN];
    double elvmask = 10; // in degree

    int ip, qp;
    int iTable;
    short *iq_buff = NULL;
    signed char *iq8_buff = NULL;

    galtime_t grx;
    double delt;
    int isamp;

    int iumd;
    int numd;

    double xyz[USER_MOTION_SIZE][3];

    char navfile[MAX_CHAR];
    char outfile[MAX_CHAR];
    char tv_file[MAX_CHAR];

    double samp_freq;
    int iq_buff_size;
    int data_format;

    int result;

    float gain[MAX_CHAN];
    double path_loss;
    double ant_gain;
    double ant_pat[37];
    int ibs; // boresight angle index

    datetime_t t0, tmin, tmax;
    galtime_t gmin, gmax;
    double dt;
    int igrx;

    double duration;
    int iduration;
    int verb;

    bool use_usrp;
    bool use_bits_from_streamer;

    int timeoverwrite = FALSE; // Overwrite the TOC and TOE in the RINEX file

    ionoutc_t ionoutc;

    ////////////////////////////////////////////////////////////
    // Read options
    ////////////////////////////////////////////////////////////

    strcpy(navfile, s->opt.navfile);
    strcpy(tv_file, s->opt.tvfile);
    strcpy(outfile, s->opt.outfile);

    g0.week = s->opt.g0.week;
    g0.sec = s->opt.g0.sec;

    gal2date(&g0, &t0);

    iduration = USER_MOTION_SIZE;

    samp_freq = (double)TX_SAMPLERATE;

    verb = s->opt.verb;

    iq_buff_size = NUM_IQ_SAMPLES;

    delt = 1.0 / (double)TX_SAMPLERATE;

    int interactive = s->opt.interactive;

    timeoverwrite = s->opt.timeoverwrite;

    ionoutc.enable = s->opt.iono_enable;

    use_usrp = s->opt.use_usrp;

    duration = s->opt.iduration / 10.0;

    use_bits_from_streamer = s->opt.use_bit_stream;

    llh[0] = s->opt.llh[0];
    llh[1] = s->opt.llh[1];
    llh[2] = s->opt.llh[2];

    ////////////////////////////////////////////////////////////
    // Start auxiliary task threads
    ////////////////////////////////////////////////////////////

    // pthread_create(&th_bits, NULL, &bitstreamer_thread, (void *)&hq);
    // std::thread thread_obj(bitstreamer_thread, std::ref(queues), std::ref(chn_prn_map), &exit_flag);

    // pthread_t th_range;
    // pthread_create(&th_range, NULL, &dt_thread, NULL);

    // Update location via UDP Socket
    std::thread th_loc(&locations_thread, llh);

    ////////////////////////////////////////////////////////////
    // Load navigation messages and satellite ephemeris
    ////////////////////////////////////////////////////////////

    vector<ephem_t> eph_vector[MAX_SAT];
    ephem_t eph;

    // Load file pointers to individual satellites
    std::cerr << s->opt.tvfile[0];
    if (s->opt.navfile[0] == 0 && s->opt.tvfile[0] != 0)
    {
        for (int sv = 0; sv < MAX_SAT; sv++)
        {
            char tfname[MAX_CHAR];
            snprintf(tfname, sizeof(tfname), "%s/%d.csv", tv_file, sv + 1);
            TV_ptrs_sv[sv] = fopen(tfname, "r");
            eph_vector[sv] = load_ephemeris(tfname);
        }
    }

    // Keep track of current ephemeris index for each satellite
    for (int i = 0; i < MAX_SAT; i++)
        current_eph.push_back(0);

    llh[0] = llh[0] / R2D; // convert to RAD
    llh[1] = llh[1] / R2D; // convert to RAD
    llh2xyz(llh, xyz[0]);  // Convert llh to xyz

    iduration = (int)(duration * 10.0 + 0.5);

    ////////////////////////////////////////////////////////////
    // Receiver position
    ////////////////////////////////////////////////////////////

    fprintf(stderr, "\nUsing static location mode.\n");
    numd = iduration;

    fprintf(stderr, "xyz = %11.1f, %11.1f, %11.1f\n", xyz[0][0], xyz[0][1],
            xyz[0][2]);
    fprintf(stderr, "llh = %11.6f, %11.6f, %11.1f\n", llh[0] * R2D, llh[1] * R2D,
            llh[2]);

    /* sv = PRN-1 */
    for (sv = 0; sv < MAX_SAT; sv++)
    {
        eph = eph_vector[sv][0];

        if (eph.vflg == 1)
        {
            gmin = eph.toc;
            gal2date(&gmin, &tmin);
            break;
        }
    }

    gmax.sec = 0;
    gmax.week = 0;
    tmax.sec = 0;
    tmax.mm = 0;
    tmax.hh = 0;
    tmax.d = 0;
    tmax.m = 0;
    tmax.y = 0;

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        eph = eph_vector[sv][eph_vector[sv].size() - 2];

        if (eph.vflg == 1)
        {
            if (eph.toc.sec > gmax.sec)
                gmax = eph.toc;
        }
    }

    gal2date(&gmax, &tmax);

    if (g0.week >= 0) // Scenario start time has been set.
    {
        if (timeoverwrite == TRUE)
        {
            galtime_t gtmp;
            datetime_t ttmp;
            double dsec;

            gtmp.week = g0.week;
            gtmp.sec = (double)(((int)(g0.sec)) / 7200) * 7200.0;

            dsec = subGalTime(gtmp, gmin);

            // Overwrite the UTC reference week number
            ionoutc.wnt = gtmp.week;
            ionoutc.tot = (int)gtmp.sec;

            // Overwrite the TOC and TOE to the scenario start time
            for (sv = 0; sv < MAX_SAT; sv++)
            {
                for (i = 0; i < neph; i++)
                {
                    if (eph1[i][sv].vflg == 1)
                    {
                        gtmp = incGalTime(eph1[i][sv].toc, dsec);
                        gal2date(&gtmp, &ttmp);
                        eph1[i][sv].toc = gtmp;
                        eph1[i][sv].t = ttmp;

                        gtmp = incGalTime(eph1[i][sv].toe, dsec);
                        eph1[i][sv].toe = gtmp;
                    }
                }
            }
        }
        else
        {
            if (subGalTime(g0, gmin) < 0.0 || subGalTime(gmax, g0) < 0.0)
            {
                datetime_t tl;
                gal2date(&g0, &tl);
                fprintf(stderr, "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec, g0.week,
                        g0.sec);

                fprintf(stderr, "ERROR: Invalid start time.\n");
                fprintf(stderr, "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmin.y, tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec, gmin.week,
                        gmin.sec);
                fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmax.y, tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec, gmax.week,
                        gmax.sec);
                exit(1);
            }
        }
    }
    else
    {
        g0 = gmin;
        t0 = tmin;
    }

    datetime_t tl;
    gal2date(&g0, &tl);

    fprintf(stderr, "\ntmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", tmin.y,
            tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec, gmin.week, gmin.sec);

    fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", tmax.y,
            tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec, gmax.week, gmax.sec);

    fprintf(stderr, "\nDuration = %.1f [sec]\n", ((double)numd) / 10.0);

    fprintf(stderr, "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
            tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec, g0.week,
            g0.sec);

    ////////////////////////////////////////////////////////////
    // Read data from file and setup bit and TOW correction queues
    ////////////////////////////////////////////////////////////

    // Initial reception time
    grx = g0;

    double dt1, correction;

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        current_eph[sv] = epoch_matcher(grx.sec, eph_vector[sv], current_eph[sv]);
        if (current_eph[sv] == -1)
        {
            fprintf(stderr, "ERROR: No current set of ephemerides has been found.\n");
            exit(1);
        }
    }

    ////////////////////////////////////////////////////////////
    // Baseband signal buffer and output file
    ////////////////////////////////////////////////////////////

    // Allocate I/Q buffer
    iq_buff = (short *)calloc(2 * iq_buff_size, sizeof(short));

    // Open output file
    // "-" can be used as name for stdout
    if (strcmp("-", outfile))
    {
        if (NULL == (fp = fopen(outfile, "wb")))
        {
            fprintf(stderr, "ERROR: Failed to open output file.\n");
            exit(1);
        }
    }
    else
    {
        fp = stdout;
    }

    ////////////////////////////////////////////////////////////
    // Initialize channels
    ////////////////////////////////////////////////////////////

    // Clear all channels
    for (i = 0; i < MAX_CHAN; i++)
    {
        // Set length of code array at runtime based on samples
        chan[i].set_code_phase = true;
        chan[i].prn = 0;
        chan[i].ca_E1B = (short *)malloc(2 * CA_SEQ_LEN_E1 * sizeof(short));
        chan[i].ca_E1C = (short *)malloc(2 * CA_SEQ_LEN_E1 * sizeof(short));
        chan[i].page = (int *)malloc(PAGE_SIZE * sizeof(int));
    }

    // Clear satellite allocation flag
    for (sv = 0; sv < MAX_SAT; sv++)
        allocatedSat[sv] = -1;

    dt = 0.10000002314200000;

    grx = incGalTime(grx, dt);

    allocateChannel(chan, eph_vector, ionoutc, grx, xyz[0], elvmask, &chn_prn_map);

    for (i = 0; i < MAX_CHAN; i++)
    {
        if (chan[i].prn > 0)
            fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.5f\n", chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.range, grx.sec);
    }

    ////////////////////////////////////////////////////////////
    // Receiver antenna gain pattern
    ////////////////////////////////////////////////////////////

    for (i = 0; i < 37; i++)
        ant_pat[i] = pow(10.0, -ant_pat_db[i] / 20.0);

    tstart = clock();

    // Update receiver time

    int outbuf_idx = 0;
    float cosPhase[1];
    float sinPhase[1];
    double temp_phase = 0;

    bool debug = false;
    float *debug_codephase;
    FILE *fp_debug;

    // Ctrl+C will exit loop
    signal(SIGINT, &sigint_handler_gal_task);

    if (debug)
    {
        fp_debug = fopen("./debug_codePhase_data_with_CCP", "wb");
        debug_codephase = (float *)calloc(iq_buff_size, sizeof(float));
    }

    if (use_bits_from_streamer)
    {
        fprintf(stderr, "\nWaiting for navigation message bits ");

        while (1)
        {
            bool flag = false;
            for (int i = 0; i < MAX_CHAN; i++)
            {
                if (!queues[i].empty())
                {
                    flag = true;
                    break;
                }
            }
            if (flag)
            {
                break;
            }
            else
            {
                sleep(1);
                fprintf(stderr, ".");
            }
        }

        fprintf(stderr, "\nBits received - Starting Generator");
        bool reset_index = false;
        fflush(stdout);
    }
    else
    {
        fprintf(stderr, "\nUsing bits from file source\n");
        fflush(stdout);
    }

    long last_nanos;
    long start;

    ////////////////////////////////////////////////////////////
    // Generate baseband signals
    ////////////////////////////////////////////////////////////

    initscr();

    int cp = 0;

    grx = incGalTime(grx, dt);

    fprintf(stderr, "\n Ibit Values: \n");

    for (iumd = 1; iumd < numd; iumd++)
    {
        start = get_nanos();

        // Copy contents of llh received over the locations thread
        memcpy(llh, llhr, 3 * sizeof(double));

        llh[0] = llh[0] / R2D; // convert to RAD
        llh[1] = llh[1] / R2D; // convert to RAD

        llh2xyz(llh, xyz[iumd]);

        for (i = 0; i < MAX_CHAN; i++)
        {
            if (chan[i].prn > 0)
            {

                // Refresh code phase and data bit counters
                range_t rho;
                sv = chan[i].prn - 1;
                eph = eph_vector[sv][current_eph[sv]];

                // Current pseudorange
                computeRange(&rho, eph, &ionoutc, grx, xyz[iumd], chan[i].prn);

                chan[i].azel[0] = rho.azel[0];
                chan[i].azel[1] = rho.azel[1];

                // Update code phase and data bit counters for the first run
                computeCodePhase(&chan[i], rho, dt, grx);

                // Path loss
                path_loss = 30200000.0 / rho.d;

                // Receiver antenna gain
                ibs = (int)((90.0 - rho.azel[1] * R2D) /
                            5.0); // covert elevation to boresight
                ant_gain = ant_pat[ibs];

                // Signal gain
                gain[i] = (path_loss * ant_gain); // scaled by 2^7
            }
        }

        for (isamp = 0; isamp < iq_buff_size; isamp++)
        {
            int i_acc = 0;
            int q_acc = 0;
            double a = get_nanos();
            // fprintf(stderr, "\n");
            for (i = 0; i < MAX_CHAN; i++)
            {
                if (chan[i].prn > 0)
                {
                    if (chan[i].code_phase >= CA_SEQ_LEN_E1)
                    {
                        chan[i].code_phase -= CA_SEQ_LEN_E1;
                        chan[i].ibit++;

                        // 250 symbols = 120 navigation data bits + 10 preamble bits = 1 page
                        if (chan[i].ibit >= N_SYM_PAGE)
                        {
                            chan[i].ibit = 0;
                            chan[i].ipage++;

                            // Generate new page
                            generateNavMsg(grx, &chan[i], advance_fptr);
                        }
                    }

                    int cosPh = cosTable512[((int)(511 * chan[i].carr_phase)) & 511];
                    int sinPh = sinTable512[((int)(511 * chan[i].carr_phase)) & 511];

                    int icode = (int)(chan[i].code_phase * 2);

                    int E1B_chip = chan[i].ca_E1B[icode];
                    int E1C_chip = chan[i].ca_E1C[icode];

                    int databit = chan[i].page[chan[i].ibit] > 0 ? -1 : 1;
                    int secCode = GALILEO_E1_SECONDARY_CODE[chan[i].ibit % 25] > 0 ? -1 : 1;

                    ip = (E1B_chip * databit - E1C_chip * secCode) * cosPh;
                    qp = (E1B_chip * databit - E1C_chip * secCode) * sinPh;

                    // Accumulate for all visible satellites
                    i_acc += ip;
                    q_acc += qp;

                    // Update code phase
                    chan[i].code_phase += chan[i].f_code * delt;

                    // Update carrier phase
                    chan[i].carr_phase += (chan[i].f_carr) * delt;
                    chan[i].carr_phase -= (long)chan[i].carr_phase; // (carr_phase %1)
                }
            }
            // Store I/Q samples into buffer
            iq_buff[isamp * 2] = (short)i_acc * 10.4;
            iq_buff[isamp * 2 + 1] = (short)q_acc * 10.25;
            // advance_fptr = true;
        }

        grx = incGalTime(grx, dt);

        fwrite(iq_buff, sizeof(short), 2 * iq_buff_size, fp);

        // Check and update ephemeris index and satellite allocation every 30 seconds
        if ((int)fmodf(grx.sec - g0.sec, 30) == 0)
        {
            allocateChannel(chan, eph_vector, ionoutc, grx, xyz[iumd], elvmask, &chn_prn_map);
            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                current_eph[sv] = epoch_matcher(grx.sec, eph_vector[sv], current_eph[sv]);
            }
        }

        ////////////////////////////////////////////////////////////
        // Write into FIFO
        ///////////////////////////////////////////////////////////

        if (use_usrp)
        {
            if (!s->galileo_sim.ready)
            {
                // Initialization has been done. Ready to create TX task.
                printf("Galileo signal generator is ready!\n");
                s->galileo_sim.ready = 1;
                pthread_cond_signal(&(s->galileo_sim.initialization_done));
            }

            long int a = get_nanos();
            // Wait utill FIFO write is ready
            pthread_mutex_lock(&(s->galileo_sim.lock));
            while (!is_fifo_write_ready(s))
                pthread_cond_wait(&(s->fifo_write_ready), &(s->galileo_sim.lock));
            pthread_mutex_unlock(&(s->galileo_sim.lock));

            double b = (get_nanos() - a) / 1e6;

            // fprintf(stderr, "\nTime elapsed : %f ms", b);

            // Write into FIFO
            memcpy(&(s->fifo[s->head * 2]), iq_buff,
                   NUM_IQ_SAMPLES * 2 * sizeof(short));

            s->head += (long)NUM_IQ_SAMPLES;
            if (s->head >= FIFO_LENGTH)
                s->head -= FIFO_LENGTH;
            pthread_cond_signal(&(s->fifo_read_ready));
        }

        //
        // Update receiver time
        if (use_bits_from_streamer)
        {

            llh[0] = llh[0] * R2D;
            llh[1] = llh[1] * R2D;

            // Update time counter
            if (!local_fix && tow_fixed)
            {
                double d = (get_nanos() - start) / 1e9;
                grx = incGalTime(grx, dt1);

                printf("\rDT- %f", d);
                printf("\nTOW Fix: %f\n", dt1);

                tow_fixed = true;
                local_fix = true;
            }
        }

        if (verb == TRUE)
        {
            if (use_cursors)
            {
                clear();
                attron(A_REVERSE);
                llh[0] = llh[0] * R2D;
                llh[1] = llh[1] * R2D;

                printw("\n Location: %10f, %10f, %4f - Time: \n", llh[0], llh[1], llh[2]);
                printw("\n Elapsed time %4.1f s\n", subGalTime(grx, g0));
                printw("\n%3s%6s%14s%17s%21s%18s%18s%18s%5s\n", "CH", "PRN", "Azimuth", "Elevation", "Doppler [Hz]", "Code phase", "rx_time", "Pseudorange", "Eph"); //,"Y","M","D","HH","MM","SS");
                attroff(A_REVERSE);
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn > 0 && use_cursors)
                        printw("%3d%6d%14f%17f%21f%18f%18f%18f%5d\n", i, chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].f_carr, chan[i].code_phase, grx.sec, chan[i].rho0.range, current_eph[chan[i].prn - 1]); //, t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec);
                }
                refresh();
            }
        }
        else
            fprintf(stderr, "\rTime into run = %4.1f - %4.1ld", subGalTime(grx, g0), clock() - tstart);

        if (stop_signal_called_gal_task)
            break;

        fflush(stdout);
    }

    tend = clock();
    exit_flag = true;
    fprintf(stderr, "\nDone!\n");
    endwin();

    // Free I/Q buffer
    free(iq_buff);

    // Close file
    fclose(fp);

    if (debug)
        fclose(fp_debug);

    // Process time
    fprintf(stderr, "Process time = %.1f [sec]\n",
            (double)(tend - tstart) / CLOCKS_PER_SEC);

    // return 1;
    return (NULL);
}
