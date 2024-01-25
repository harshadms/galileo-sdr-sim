#include "../include/galileo-sdr.h"


/*! \brief Converts hex E1 PRN code to binary (-1,1)
 *  \param[in] dest Binary code
 *  \param[in] prn PRN of satellite
 *  \param[in] c1 E1B/E1C selector 
 */
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

/*! \brief SBOC Modulates input array 
 *  \param[in] in_prn_ca Incoming CA code
 *  \param[in] dest destination
 *  \param[in] len 
 */
void sboc(short *dest, short *in_prn_ca, int len, int m, int n)
{
    constexpr uint32_t length_in = CA_SEQ_LEN_E1;
    const auto period = static_cast<uint32_t>(CA_SEQ_LEN_E1 / length_in);

    int i, j, N = 2 * m / n;
    for (i = 0; i < len; i++)
        for (j = 0; j < N; j++)
            dest[N * i + j] = in_prn_ca[i];

    // Mix sub carrier
    for (i = 0; i < N * len / 2; i++)
    {
        dest[2 * i] = -dest[2 * i];
    }
}

/*! \brief generate the BOC E1B code sequence for a given Satellite Vehicle PRN
 *  \param[in] prn PRN number of the Satellite Vehicle
 *  \param[out] ca Caller-allocated integer array of 1023 bytes
 */
void codegen_E1B(short *ca, int prn)
{
    // Get the PRN code using the hex to binary converter
    short *tmp_ca = (short *)malloc(CA_SEQ_LEN_E1 * sizeof(short));
    hex_to_binary_converter(tmp_ca, false, prn - 1);
    sboc(ca, tmp_ca, CA_SEQ_LEN_E1, 1, 1);
}

void codegen_E1C(short *ca, int prn)
{
    // Get the PRN code using the hex to binary converter
    short *tmp_ca = (short *)malloc(CA_SEQ_LEN_E1 * sizeof(short));
    hex_to_binary_converter(tmp_ca, true, prn - 1);
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

    // Add ionospheric delay
	rho->iono_delay = ionosphericDelay(ionoutc, g, llh, rho->azel);

	rho->range += rho->iono_delay;
    rho->g = g;

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
