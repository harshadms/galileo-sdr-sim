#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <fstream>
#include "../include/galileo-sdr.h"

using namespace std;

void ConvertD2E(char *str)
{
	while (*str)
	{
		if (*str == 'D')
			*str = 'E';
		str ++;
	}
}

unsigned char GetGalileoUra(double data)
{
	int value = (int)(data * 100);	// convert to cm

	if (value < 0 || value > 6000)
		return 255;
	if (value < 50)
		return (unsigned char)value;
	else if (value < 100)
		return (unsigned char)((value - 50) / 2 + 50);
	else if (value < 200)
		return (unsigned char)((value - 100) / 4 + 75);
	else
		return (unsigned char)((value - 200) / 16 + 100);
}

int ReadContentsData(char *str, double *data, datetime_t *time, bool read_time)
{
    int Second;
    int svid = 0;
	int length = strlen(str);

	ConvertD2E(str);

    if (read_time)
    {    
        sscanf(str+4, "%d %d %d %d %d %d", &(time->y), &(time->m), &(time->d), &(time->hh), &(time->mm), &Second);
	    time->sec = (double)Second;
	    if (str[1] == ' ') svid = 0; else sscanf(str+1, "%2d", &svid);
        if (length > 24 && str[24] != ' ') sscanf(str+23, "%lf", &data[0]); else data[0] = 0.0;
        if (length > 43 && str[43] != ' ') sscanf(str+42, "%lf", &data[1]); else data[1] = 0.0;
        if (length > 62 && str[62] != ' ') sscanf(str+61, "%lf", &data[2]); else data[2] = 0.0;
    }
    else
    {
	    if (length >  5 & str[ 5] != ' ') sscanf(str+ 4, "%lf", &data[0]); else data[0] = 0.0;
        if (length > 24 & str[24] != ' ') sscanf(str+23, "%lf", &data[1]); else data[1] = 0.0;
        if (length > 43 && str[43] != ' ') sscanf(str+42, "%lf", &data[2]); else data[2] = 0.0;
        if (length > 62 && str[62] != ' ') sscanf(str+61, "%lf", &data[3]); else data[3] = 0.0;
    }

    return svid;
}

// https://server.gage.upc.edu/gLAB/HTML/GALILEO_Navigation_Rinex_v3.04.html
int main()
{    
    vector<ephem_t> eph_vector[MAX_SAT];
    ephem_t blank_eph = {0};
    ionoutc_t ionoutc;

    FILE *fp;

    fp = fopen("/home/harshad/Work/galileo-sdr-sim-private/rinex_files/week171.rnx" , "r");
    if(fp == NULL) {
        perror("Error opening file");
        return(-1);
    }

    // Parse header
    // Ionospheric correction
    // UTC

    char str[MAX_CHAR];

    while (1)
    {
        // getline(myfile, line);
        if( fgets (str, MAX_CHAR, fp)==NULL ) 
            break;

        if(strncmp(str + 60, "END OF HEADER", 13) == 0)
            break;

        // Ionospheric corrections ai0, ai1, ai2
        if(strncmp(str + 60, "IONOSPHERIC CORR", 16) == 0)
        {
            ConvertD2E(str);
            sscanf(str + 4, "%lf %lf %lf %lf", &(ionoutc.ai0), &(ionoutc.ai1), &(ionoutc.ai2), &(ionoutc.ai3));
            // char *pch;
            // pch = strtok (str," ");
            // int index=0;

            // while (pch != NULL)
            // { 
            //     if (index > 0)
            //         ConvertD2E(pch);

            //     switch (index)
            //     {
            //         case 1:
            //             ionoutc.ai0 = atof(pch);
            //             break;
            //         case 2:
            //             ionoutc.ai1 = atof(pch);
            //             break;
            //         case 3:
            //             ionoutc.ai2 = atof(pch);
            //             break;
            //         case 4:
            //             ionoutc.ai3 = atof(pch);
            //             break;
            //     }

            //     pch = strtok (NULL," ");
            //     index++;
            // }
        }

        // Time corrections GAUT - GAL to UTC
        if(strncmp(str + 60, "TIME SYSTEM CORR", 16) == 0 && strncmp(str, "GAUT", 4) == 0)
        {
            char ch;
	        int data1, data2;

            ConvertD2E(str);
            ch = str[22]; str[22] = 0; // put string terminator on first data
            sscanf(str + 4, "%lf", &(ionoutc.A0));

            str[22] = ch;	// restore char after first data
            sscanf(str + 22, "%lf %d %d", &(ionoutc.A1), &data1, &data2);
            ionoutc.A2 = 0.0;
            ionoutc.tot = (unsigned char)(data1 >> 12);
            ionoutc.wnt = (short)data2;
            ionoutc.wnlsf = ionoutc.wnt;
            ionoutc.dn = 0;
        }

        // Leap seconds
    }

    // Parse ephemeris data
    while (1)
    {
        if( fgets (str, MAX_CHAR, fp)==NULL ) 
            break;

        // EPH for new satellite - process 8 lines
        if( str[0] == 'E' )
        {
            double data[39];
            datetime_t utctime;
            galtime_t galtime;

            // Extract epoch and SVID
            int svid = ReadContentsData(str, &data[0], &utctime, true);

            // Push back a blank eph
            ephem_t eph;

            // Parse through remaining data lines
            for (int i=0; i < 7; i++)
            {
                fgets(str, MAX_CHAR, fp);
                ReadContentsData(str, &data[i*4+3], &utctime, false);
            }

            //date2gal(&utctime, &galtime);

            eph.svid = svid;
            eph.toc = galtime;
            eph.af0 = data[0];
            eph.af1 = data[1];
            eph.af2 = data[2];
            eph.sqrta = data[10];
            eph.ecc = data[8];
            eph.inc0 = data[15];
            eph.omg0 = data[13];
            eph.aop = data[17];
            eph.m0 = data[6];
            eph.deltan = data[5];
            eph.omgdot = data[18];
            eph.idot = data[19];
            eph.crc = data[16];
            eph.crs = data[4];
            eph.cuc = data[7];
            eph.cus = data[9];
            eph.cic = data[12];
            eph.cis = data[14];
            eph.toe.sec = (int)(data[11] + 0.5); 
            eph.toe.week = (int)data[21];      /* week number */
            eph.iode = (unsigned char)data[3];      /* IODE/AODE */

            eph.svhlth = (unsigned short)data[24];      /* sv health */

			eph.ura = GetGalileoUra(data[23]);
			eph.flag = (unsigned short)data[20];
			eph.bgde5a = data[25];      /* TGD */
			eph.bgde5b = (eph.flag & 0x2) ? data[25] : data[26];      /* TGD for E1/E5b */
			eph.tgd_ext[2] = eph.bgde5a * TGD_GAMME_L5;	/* TGD for E5a */
			eph.tgd_ext[4] = eph.bgde5a * TGD_GAMMA_E5b;	/* TGD for E5b */

            eph.A = eph.sqrta * eph.sqrta;
            eph.n = WGS_SQRT_GM / (eph.sqrta * eph.A) + eph.deltan;
            eph.sq1e2 = sqrt(1.0 - eph.ecc * eph.ecc);
            eph.omg_t = eph.omg0 - OMEGA_EARTH * eph.toe.sec;
            eph.omgkdot = eph.omgdot - OMEGA_EARTH;
            eph.vflg = true;
            eph.PRN = svid;

            eph_vector[svid-1].push_back(eph);
        }
    }
    
    for (int i=0; i<MAX_SAT; i++)
    {
        if (eph_vector[i].size() == 0)
            continue;
        cout << "Loaded " << eph_vector[i].size() << " records for " << i+1 << endl;
    }
}