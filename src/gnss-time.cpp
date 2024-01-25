#include "../include/galileo-sdr.h"

/*! \brief Convert a UTC date into a GAL damte
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

// Source: BOOK called GPS Theory Algorithm & Applications by Guochang Xu (Pg 18-20)
double gps_time(datetime_t *t) {
	// As precaution, check if we have required epoch info

	// Year Month Day of epoch
	double y = t->y; double m = t->m; double d = t->d;
	// Hour Minute Second of epoch
	double hr = t->hh; double min = t->mm; double sec = t->sec;
	// UTC time in hours
	double UTC = hr + min / 60. + sec / 3600.;
	// Taking care of month and year conditioning
	if (m > 2) {
		y = y + 2000;
	}
	else {
		y = y + 2000 - 1;
		m = m + 12;
	}
	// Julian Date
	double jDate = floor(365.25*y) + floor(30.6001*(m + 1)) + d + (UTC / 24) + 1720981.5;
	// GPS Week
	double gpsWeek = floor((jDate - 2444244.5) / 7);
	// GPS Time
	double gpsTime = round(((((jDate - 2444244.5) / 7) - gpsWeek) * 7 * 24 * 3600) * 100) / 100;
	// Return GPS Time in seconds
	return gpsTime;
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

long get_nanos()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}

void set_scenario_start_time(galtime_t *g0, galtime_t gmin, galtime_t gmax, datetime_t *t0, datetime_t *tmin, datetime_t *tmax, bool timeoverwrite, ionoutc_t *ionoutc, int neph, vector<ephem_t> eph1[MAX_SAT])
{
    if (g0->week >= 0) // Scenario start time has been set.
    {
        if (timeoverwrite == TRUE)
        {
            galtime_t gtmp;
            datetime_t ttmp;
            double dsec;

            gtmp.week = g0->week;
            gtmp.sec = (double)(((int)(g0->sec)) / 7200) * 7200.0;

            dsec = subGalTime(gtmp, gmin);

            // Overwrite the UTC reference week number
            ionoutc->wnt = gtmp.week;
            ionoutc->tot = (int)gtmp.sec;

            // Overwrite the TOC and TOE to the scenario start time
            for (int sv = 0; sv < MAX_SAT; sv++)
            {
                for (int i = 0; i < neph; i++)
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
            if (subGalTime(*g0, gmin) < 0.0 || subGalTime(gmax, *g0) < 0.0)
            {
                datetime_t tl;
                gal2date(g0, &tl);
                fprintf(stderr, "Start = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tl.y, tl.m, tl.d, tl.hh, tl.mm, tl.sec, g0->week,
                        g0->sec);

                fprintf(stderr, "ERROR: Invalid start time.\n");
                fprintf(stderr, "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmin->y, tmin->m, tmin->d, tmin->hh, tmin->mm, tmin->sec, gmin.week,
                        gmin.sec);
                fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
                        tmax->y, tmax->m, tmax->d, tmax->hh, tmax->mm, tmax->sec, gmax.week,
                        gmax.sec);
                exit(1);
            }
        }
    }
    else
    {
        g0->sec = gmin.sec;
        g0->week = gmin.week;

        t0 = tmin;
    }

}