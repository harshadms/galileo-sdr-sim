#include "../include/galileo-sdr.h"
#include "../include/logging.h"

/* Matches TOW and relevant ephemeris */
int epoch_matcher(galtime_t obsTime, vector<ephem_t> eph, int ind2ex)
{
	// for (int i = index; i < eph.size(); i++)
	// {
	// 	double diff = fabs(obsTime - eph.at(i).toe.sec); 
    //     cout << "EM: " << obsTime << " : " <<  eph.at(i).toe.sec << " - " << eph.at(i).svid + 1 << endl;
	// 	if (diff <= 0)
	// 	{
	// 		index = i - 1; // save index
	// 		break;
	// 	}
	// }
	// return index;
    // Initialize time difference variable using arbitrary large number
	// double diff = 1000000; int index = 0;
	// for (unsigned i = 0; i < eph.size(); i++) {
	// 	double new_diff = fabs(obsTime - eph.at(i).gps_time);
    //     cout << "EM: " << obsTime -  eph.at(i).gps_time << " - " << eph.at(i).svid + 1 << endl;
	// 	if (new_diff < diff) {
	// 		diff = new_diff; // update difference
	// 		index = i; // save index
	// 	}
	// }
    int index = -1;
    double dt;
    for (unsigned i = 0; i < eph.size(); i++) {
        if (eph.at(i).vflg == 1)
        {
            dt = subGalTime(obsTime, eph.at(i).toc);
            
            if (dt>=-SECONDS_IN_HOUR && dt<SECONDS_IN_HOUR)
            {
                //cout << "EM: " << dt << " - " << eph.at(i).svid + 1 << endl;
                index = i;
                break;
            }
        }
    }
	// Index of most appropriate Nav vector
	return index;
}

void convertD2E(char *str)
{
	while (*str)
	{
		if (*str == 'D')
			*str = 'E';
		str ++;
	}
}

unsigned char getGalileoUra(double data)
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

int readContentsData(char *str, double *data, datetime_t *time, bool read_time)
{
    int Second;
    int svid = 0;
	int length = strlen(str);

	convertD2E(str);

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
int readRinexV3(vector<ephem_t> eph_vector[MAX_SAT], ionoutc_t *ionoutc, char *fname)
{
    int eph_count = 0;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        perror("Error opening file");
        return (-1);
    }

    // Default leap seconds (GST-UTC) if not in header
    ionoutc->dtls = 18;

    char str[MAX_CHAR];
    // Parse header
    while (fgets(str, MAX_CHAR, fp) != NULL)
    {
        if (strncmp(str + 60, "END OF HEADER", 13) == 0)
            break;

        // Ionospheric corrections
        if (strncmp(str + 60, "IONOSPHERIC CORR", 16) == 0 && strncmp(str, "GAL", 3) == 0)
        {
            convertD2E(str);
            sscanf(str + 4, "%lf %lf %lf", &(ionoutc->ai0), &(ionoutc->ai1), &(ionoutc->ai2));
        }

        // Time corrections GAUT - GAL to UTC
        if (strncmp(str + 60, "TIME SYSTEM CORR", 16) == 0 && strncmp(str, "GAUT", 4) == 0)
        {
            convertD2E(str);
            char ch = str[22];
            str[22] = 0;
            sscanf(str + 4, "%lf", &(ionoutc->A0));
            str[22] = ch;
            int data1, data2;
            sscanf(str + 22, "%lf %d %d", &(ionoutc->A1), &data1, &data2);
            ionoutc->A2 = 0.0;
            ionoutc->tot = (unsigned char)(data1 >> 12);
            ionoutc->wnt = (short)data2 >> 4;
            ionoutc->wnlsf = (short)data2;
        }

        // Leap seconds
        if (strncmp(str + 60, "LEAP SECONDS", 12) == 0)
        {
            int leap;
            if (sscanf(str, "%d", &leap) == 1)
                ionoutc->dtls = leap;
        }
    }

    // Parse ephemeris data
    while (fgets(str, MAX_CHAR, fp) != NULL)
    {
        if (str[0] != 'E')
            continue; // Only Galileo for now

        double data[39];
        datetime_t utctime;
        galtime_t galtime;

        int svid = readContentsData(str, &data[0], &utctime, true);
        if (svid <= 0 || svid > MAX_SAT)
            continue;

        ephem_t eph;
        bool truncated = false;
        for (int i = 0; i < 7; i++)
        {
            if (fgets(str, MAX_CHAR, fp) == NULL)
            {
                truncated = true;
                break;
            }
            readContentsData(str, &data[i * 4 + 3], &utctime, false);
        }

        if (truncated)
            break;

        date2gal(&utctime, &galtime);

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
        eph.toe.week = (int)data[21];
        
        // GPS Week Rollover Correction (GPSWeek is 10-bit, rolls over every 1024 weeks ~19.7 years)
        // GPS epoch started Jan 6, 1980. As of March 2026, GPS week is ~2410.
        // The 10-bit GPS week in ephemeris rolls: 0-1023 (1980-1999), 0-1023 (1999-2019), 0-1023 (2019-2038)
        // We need to determine which cycle we're in and adjust accordingly.
        int current_gps_week = (int)data[21];  // 10-bit week from ephemeris
        
        // Estimate current GPS week from system time
        time_t now = time(NULL);
        time_t gps_epoch = 315964800;  // Jan 6, 1980 00:00:00 UTC
        long seconds_since_epoch = now - gps_epoch;
        int estimated_gps_week = (int)(seconds_since_epoch / (7 * 86400));
        
        // Determine which 1024-week cycle we should be in
        int week_cycle = (estimated_gps_week / 1024) * 1024;
        
        // Adjust ephemeris week to current cycle
        int adjusted_week = week_cycle + current_gps_week;
        if (abs(estimated_gps_week - adjusted_week) > 512) {
            // Wrong cycle, try adjacent cycle
            adjusted_week = (week_cycle - 1024) + current_gps_week;
            if (abs(estimated_gps_week - adjusted_week) > 512) {
                adjusted_week = (week_cycle + 1024) + current_gps_week;
            }
        }
        
        log_debug("GPS week rollover: raw=%d adjusted=%d (cycle=%d, current=%d)", 
                  current_gps_week, adjusted_week, week_cycle, estimated_gps_week);
        
        eph.toe.week = adjusted_week;  // Use corrected full GPS week
        
        eph.iode = (unsigned char)data[3];
        eph.svhlth = (unsigned short)data[24];
        eph.ura = getGalileoUra(data[23]);
        eph.flag = (unsigned short)data[20];

        // Basic filtering for E1B-I/NAV which is what we simulate
        if (eph.flag != 517 && eph.flag != 257) // 517 = E1B/E5b, 257 = E1B
            continue;

        eph.bgde5a = data[25];
        eph.bgde5b = (eph.flag & 0x2) ? data[25] : data[26];
        eph.tgd_ext[2] = eph.bgde5a * TGD_GAMME_L5;
        eph.tgd_ext[4] = eph.bgde5a * TGD_GAMMA_E5b;

        eph.A = eph.sqrta * eph.sqrta;
        eph.n = WGS_SQRT_GM / (eph.sqrta * eph.A) + eph.deltan;
        eph.sq1e2 = sqrt(1.0 - eph.ecc * eph.ecc);
        eph.omg_t = eph.omg0 - OMEGA_EARTH * eph.toe.sec;
        eph.omgkdot = eph.omgdot - OMEGA_EARTH;
        eph.vflg = 1;
        eph.PRN = svid;
        eph.gps_time = data[27];

        eph_vector[svid - 1].push_back(eph);
        eph_count++;
    }
	// exit(1);

    for (int i=0; i<MAX_SAT; i++)
    {
        if (eph_vector[i].size() == 0)
            continue;
        log_debug("Loaded %zu ephemeris records for satellite %d", eph_vector[i].size(), i+1);
    }

    return eph_count;
}