#include "../include/galileo-sdr.h"


/*! \brief Converts hex E1 PRN code to binary (-1,1)
 *  \param[in] dest Binary code
 *  \param[in] prn PRN of satellite
 *  \param[in] c1 E1B/E1C selector 
 */
void hex_to_binary_converter(short *dest, bool c1, int prn)
{
    int index = 0;

    // Galileo E1 codes are 4092 chips. 
    // Each hex character represents 4 chips.
    // 1023 * 4 = 4092.
    for (int i = 0; i < 1023; i++)
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

/*! \brief SBOC Modulates input array 
 *  \param[in] in_prn_ca Incoming CA code
 *  \param[in] dest destination
 *  \param[in] len 
 */
void sboc(short *dest, short *in_prn_ca, int len, int m, int n)
{
    // Galileo E1 BOC(1,1) modulation
    // Each chip of the primary code (1.023 Mcps) is mixed with a subcarrier (1.023 MHz)
    // The subcarrier has two 'half-chips' per code chip.
    // Length of 'dest' should be 2 * len.
    
    for (int i = 0; i < len; i++)
    {
        // For BOC(1,1), we have 2 subcarrier chips per code chip
        // Subcarrier sequence for one code chip is [-1, 1] OR [1, -1] depending on implementation
        // According to Galileo ICD, the subcarrier is s_E1_boc = sgn(sin(pi * f_sc * t))
        // This results in two sub-chips per code chip with opposite signs.
        dest[2 * i] = (short)(in_prn_ca[i]);      // First half-chip
        dest[2 * i + 1] = (short)(-in_prn_ca[i]); // Second half-chip (inverted)
    }
}

/*! \brief generate the BOC E1B code sequence for a given Satellite Vehicle PRN
 *  \param[in] prn PRN number of the Satellite Vehicle
 *  \param[out] ca Caller-allocated integer array of 1023 bytes
 */
void codegen_E1B(short *ca, int prn)
{
    // Get the PRN code using the hex to binary converter
    short tmp_ca[CA_SEQ_LEN_E1];
    hex_to_binary_converter(tmp_ca, false, prn - 1);
    // SBOC produces 2 samples per chip, so destination 'ca' must be 2*CA_SEQ_LEN_E1
    sboc(ca, tmp_ca, CA_SEQ_LEN_E1, 1, 1);
}

void codegen_E1C(short *ca, int prn)
{
    // Get the PRN code using the hex to binary converter
    short tmp_ca[CA_SEQ_LEN_E1];
    hex_to_binary_converter(tmp_ca, true, prn - 1);
    // SBOC produces 2 samples per chip, so destination 'ca' must be 2*CA_SEQ_LEN_E1
    sboc(ca, tmp_ca, CA_SEQ_LEN_E1, 1, 1);
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

    // Sagnac effect correction (Earth rotation during signal propagation tau)
    // The ECEF frame rotates CCW with the Earth. To express the satellite's position
    // at the time of transmission in the ECEF frame at the time of reception, we must
    // rotate the satellite's coordinates CW by alpha = omega * tau.
    double alpha = GNSS_OMEGA_EARTH_DOT * tau;
    double sin_a = sin(alpha);
    double cos_a = cos(alpha);
    double x_new = pos[0] * cos_a + pos[1] * sin_a;  // Formal CW rotation
    double y_new = -pos[0] * sin_a + pos[1] * cos_a;
    pos[0] = x_new;
    pos[1] = y_new;

    // New observer to satellite vector and satellite range.
    subVect(los, pos, xyz);
    range = normVect(los);

    rho->d = range;

    // Pseudorange (Initial: geometric + clock bias)
    rho->range = range - SPEED_OF_LIGHT * clk[0];

    // Azimuth and elevation angles.
    double satLLH[3];
    xyz2llh(xyz, llh);     // convert userXYZ to llh
    xyz2llh(pos, satLLH);  // convert satXYZ to llh
    ltcmat(llh, tmat);
    ecef2neu(los, tmat, neu);
    neu2azel(rho->azel, neu);

    // Ionospheric delay
    double frequency = CARR_FREQ;
    rho->iono_delay = ionosphericDelay(ionoutc, g, llh, satLLH, rho->azel, frequency);

    // Tropospheric delay (Saastamoinen model)
    double height = llh[2];
    rho->tropo_delay = troposphericDelay(rho->azel, height);

    // Final Pseudorange: sum all delays (they INCREASE travel time)
    // S04: Re-enable atmospheric corrections for more realistic simulation
    // GNSS-SDR is configured to NOT apply atmospheric corrections when they're already in the pseudorange
    // This provides a middle ground: cleaner simulation than all-geometric, more realistic than totally clean
    rho->range += rho->iono_delay + rho->tropo_delay;
    rho->g = g;

    return;
}

/*! \brief Saastamoinen tropospheric delay model
 *  \param azel Azimuth and Elevation in radians
 *  \param height Receiver height in meters
 *  \return Tropospheric delay in meters
 */
double troposphericDelay(double azel[2], double height)
{
    const double Pressure = 1013.25; // Standard pressure at sea level [hPa]
    const double Temp = 288.15;      // Standard temperature at sea level [K]
    const double Humid = 0.5;       // Relative humidity (typical)

    double el = azel[1]; // Elevation in radians
    if (el < 0.05) el = 0.05; // Elevation mask to prevent singularity

    // Saastamoinen model: ZTD = 0.002277 * (Pressure + (1255/Temp + 0.05) * PartialPressureWaterVapor)
    // Simplified version: 2.31 * exp(-0.000116 * height) / sin(el)
    double ztd = 2.312 * exp(-0.000116 * height);
    return ztd / sin(el);
}

/*! \brief Compute the code phase for a given channel (satellite)
 *  \param chan Channel on which we operate (is updated)
 *  \param[in] rho1 Current range, after \a dt has expired
 *  \param[in dt delta-t (time difference) in seconds
 */
void computeCodePhase(channel_t *chan, range_t rho1, double dt, galtime_t grx)
{
    double rhorate;

    // Pseudorange rate.
    rhorate = (rho1.range - chan->rho0.range) / dt;

    // Carrier and code frequency.
    chan->f_carr = (-rhorate / LAMBDA_E1); 
    chan->f_code = CODE_FREQ_E1 + chan->f_carr * CARR_TO_CODE_E1;

    // Save current pseudorange
    chan->g0 = grx;
    chan->rho0 = rho1;
}
