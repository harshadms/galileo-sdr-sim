
#include "../include/galileo-sdr.h"


/*! \brief Calculates ionospheric delay (Obliquity model)
 *  \param[out] iono_delay
 *  \param[in] user_llh user LLH coordinate
 *  \param[in] azel satellite Azimuth & Elevation
 */
double ionosphericDelay(double *user_llh, double *azel) {
    double iono_delay = 0.0;
    double E, F;
    E = azel[1] / PI;

    // Obliquity factor
    F = 1.0 + 16.0 * pow((0.53 - E), 3.0);
    iono_delay = F * 5.0e-9 * SPEED_OF_LIGHT;

    return (iono_delay); // No ionospheric delay
}

/*! \brief Calculates ionospheric delay (NequickG model)
 *  \param[out] iono_delay
 *  \param[in] ionoutc IonoUTC model
 *  \param[in] g current time
 *  \param[in] user_llh user LLH coordinate
 *  \param[in] sat_llh satellite LLH coordinate
 *  \param[in] freq frequency
 *  \param[in] azel satellite Azimuth & Elevation
 */
double ionosphericDelay(const ionoutc_t *ionoutc, galtime_t g, double *user_llh, double *sat_llh, double *azel, double freq)
{
    double iono_delay = 0.0;
    if (ionoutc->enable == FALSE) { // if we don't want to take ionosphere delay into account of PR
        iono_delay = 0.0; // No ionospheric delay
        return (iono_delay);
    }
    if (ionoutc->vflg == FALSE) { //if the ai0,ai1,ai2 parameters were not supplied in the rinex file
        iono_delay = ionosphericDelay(user_llh, azel); //calc iono using obliquity model
        return (iono_delay);
    }
    double satLLH[3] = {0};
    double userLLH[3] = {0};
    for (int i=0;i<3;i++) {
        satLLH[i] = sat_llh[i]; userLLH[i] = user_llh[i];
    }
    satLLH[2] = sat_llh[2] / 1000.0;
    userLLH[2] = user_llh[2] / 1000.0;

    return (iono_delay);
}
