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

    nequickG_input dataInput;
    bool invalidFlag = false;
    dataInput.ai[0] = ionoutc->ai0;
    dataInput.ai[1] = ionoutc->ai1;
    dataInput.ai[2] = ionoutc->ai2;
    datetime_t currentDate;
    gal2date(&g, &currentDate);
    dataInput.mth = currentDate.m;
    dataInput.UT = currentDate.hh;
    choose_F2_Fm(dataInput.F2_1, dataInput.F2_2, dataInput.Fm3_1, dataInput.Fm3_2, dataInput.mth);
    dataInput.modip = calcMODIP(userLLH);
    dataInput.Az = calcAz(dataInput.ai, dataInput.modip);
    dataInput.Az_R = sqrt(167273 + ((dataInput.Az - 63.7) * 1123.6)) - 408.99;

    double TEC = NeQuickG(userLLH, satLLH, dataInput, invalidFlag);
    double rangeError = TEC * 40.3 / pow(freq, 2);
    iono_delay = rangeError / SPEED_OF_LIGHT_M_S;
    if (invalidFlag) //if the model reached invalid data for calculation
        iono_delay = ionosphericDelay(user_llh, azel); //calc iono using obliquity model
    return (iono_delay);
}

void getCurrentPos(perigee_t ray, double height, double* userLLH, double* resultLLH) {
    double radius_km = sqrt(pow(height, 2) + pow(ray.rp, 2));

    if (radius_km < NEQUICK_G_RE_KM) {
        resultLLH[2] = 0.0;
    }
    resultLLH[2] = (radius_km - NEQUICK_G_RE_KM);

    double tanDelta = height / ray.rp;
    double cosDelta = 1.0 / sqrt(1.0 + pow(tanDelta, 2));
    double sinDelta = tanDelta * cosDelta;

    // latitude
    double sinLat = (sin(userLLH[0]) * cosDelta) + (cos(userLLH[0]) * sinDelta * ray.cossigp);
    double cosLat = sqrt(1.0 - pow(sinLat, 2));
    resultLLH[0] = atan2(sinLat, cosLat);

    // longitude
    double sinDeltaLambda = sinDelta * ray.sinsigp * cos(userLLH[0]);
    double cosDeltaLambda = cosDelta - (sin(userLLH[0]) * sin(resultLLH[0]));
    double deltaLambda = atan2(sinDeltaLambda, cosDeltaLambda);
    resultLLH[1] = deltaLambda + ray.lonp;

    return;
}

void getCurrentPosVertical(perigee_t ray, double* userLLH, double* satLLH, double* resultLLH) {
    // lat & lon
    resultLLH[0] = userLLH[0];
    resultLLH[1] = userLLH[1];

    // height
    resultLLH[2] = satLLH[2];
    return;
}

void cooR2DPath(double s, perigee_t iperigee, double* llh_pointS) {
    double hs = sqrt(pow(s, 2) + pow(iperigee.rp, 2)) - NEQUICK_G_RE_KM;

    double tan_deltaS = s / iperigee.rp;
    double cos_deltaS = 1 / sqrt(1 + pow(tan_deltaS, 2));
    double sin_deltaS = tan_deltaS * cos_deltaS;

    double sin_phiS = (sin(iperigee.latp) * cos_deltaS) + (cos(iperigee.latp) * sin_deltaS * iperigee.cossigp);
    double cos_phiS = sqrt(1 - pow(sin_phiS, 2));
    double phi_s = atan2(sin_phiS, cos_phiS);

    double sinDiff = sin_deltaS * iperigee.sinsigp * cos(iperigee.latp);
    double cosDiff = cos_deltaS - (sin(iperigee.latp) * sin_phiS);
    double lambda_s = atan2(sinDiff, cosDiff) + iperigee.lonp;

    llh_pointS[0] = phi_s;
    llh_pointS[1] = lambda_s;
    llh_pointS[2] = hs;
}

perigee_t calcPerigee(double* llh1, double* llh2, bool invalidFlag) {
    perigee_t perigee;

    perigee.isVertical = (fabs(llh1[0] - llh2[0]) < 1e-5) && (fabs(llh1[1] - llh2[1]) < 1e-5);
    if (perigee.isVertical) {
        perigee.latp = llh1[0];
        perigee.lonp = llh1[1];
        perigee.sinlatp = sin(llh1[0]);
        perigee.coslatp = cos(llh1[0]);
        return perigee;
    }

    double r1 = llh1[2] + NEQUICK_G_RE_KM;
    double r2 = llh2[2] + NEQUICK_G_RE_KM;

    //Eq. 153 -> 156
    double cos_delta = (sin(llh1[0]) * sin(llh2[0])) + (cos(llh1[0]) * cos(llh2[0]) * cos(llh2[1] - llh1[1]));
    double sin_delta = sqrt(1 - pow(cos_delta, 2));
    double zeta = atan2(sin_delta, (cos_delta - (r1 / r2)));
    perigee.rp = r1 * sin(zeta);
    //See F.2.5.1. NeQuick internal function NeqGetRayProperties
    if ((abs(zeta) > 90.0) && (perigee.rp < NEQUICK_G_RE_KM)) {
        invalidFlag = true;
        return perigee;
    }

    if (abs(abs(llh1[0]) - 90) < 1e-10) {
        perigee.latp = (llh1[0] > 0) ? (zeta) : (zeta * (-1));
        if (zeta >= 0)   //Eq. 164
            perigee.lonp = llh2[2] + GNSS_PI;
        else     //(zeta < 0)
            perigee.latp = llh2[2];
    }
    else {
        double sin_sigma, cos_sigma;
        sin_sigma = sin(llh2[1] - llh1[1]) * cos(llh2[0]) / sin_delta;
        cos_sigma = (sin(llh2[0]) - (cos_delta * sin(llh1[0]))) / (sin_delta * cos(llh1[0]));
        double delta_p = (GNSS_PI / 2) - zeta;
        perigee.sinlatp = (sin(llh1[0]) * cos(delta_p)) - (cos(llh1[0]) * sin(delta_p) * cos_sigma);
        perigee.coslatp = sqrt(1 - pow(perigee.sinlatp, 2));
        perigee.latp = atan2(perigee.sinlatp, perigee.coslatp);

        double tempSin, tempCos;
        tempSin = (-1) * sin_sigma * sin(delta_p) / perigee.coslatp;
        tempCos = (cos(delta_p) - (sin(llh1[0]) * perigee.sinlatp)) / (cos(llh1[0]) * perigee.coslatp);
        perigee.lonp = atan2(tempSin, tempCos) + llh1[1];
    }

    double psi;
    if (abs(abs(perigee.latp) - 90) < 1e-10) {
        psi = abs(llh2[0] - perigee.latp); //Eq. 168
        perigee.sinsigp = 0; //Eq. 172
        perigee.cossigp = (perigee.latp > 0) ? -1.0 : 1.0;  //Eq. 173
    }

    else { //Eq. 168 -> 171
        double cos_psi = (perigee.sinlatp * sin(llh2[0])) + (perigee.coslatp * cos(llh2[0]) * cos(llh2[1] - perigee.lonp));
        double sin_psi = sqrt(1 - pow(cos_psi, 2));
        psi = atan2(sin_psi, cos_psi);

        perigee.sinsigp = cos(llh2[0]) * sin(llh2[1] - perigee.lonp) / sin_psi;
        perigee.cossigp = (sin(llh2[0]) - (perigee.sinlatp * cos_psi)) / (perigee.coslatp * sin_psi);
    }

    if (!perigee.isVertical)
        llh1[0] = atan2(perigee.sinlatp, perigee.coslatp);

    return perigee;
}


double Epst(double X, double Y, double Z, double W) {
    double exponent = exp((W - Y) / Z);
    return ((X * exponent) / pow(1 + exponent, 2));
}

double interpolate3D(double z1, double z2, double z3, double z4, double x) {
    double zx;

    if (abs(2 * x) < 1e-10)
        zx = z2;
    else {
        double delta = 2 * x - 1;
        double g1 = z3 + z2;
        double g2 = z3 - z2;
        double g3 = z4 + z1;
        double g4 = (z4 - z1) / 3;
        double a0 = 9 * g1 - g3;
        double a1 = 9 * g2 - g4;
        double a2 = g3 - g1;
        double a3 = g4 - g2;
        zx = (1.0 / 16.0) * (a0 + (a1 * delta) + (a2 * pow(delta, 2)) + (a3 * pow(delta, 3)));
    }
    return zx;
}

double calcMODIP(double *llh) {
    double mu;
    double lat = R2D * llh[0];
    double lon = R2D * llh[1];
    if (lat >= 90) // Extreme Case 1
        mu = 90.0;
    else if (lat <= -90) // Extreme Case 2
        mu = -90.0;
    else { //Main Computations (-90deg < lat < +90deg)
        double lon_index_with_offset = (lon + 180.0) / 10.0;
        int l = int(floor(lon_index_with_offset));
        double y = lon_index_with_offset - double(l);
        if (l < 0)
            l += 36;
        else if (l >= 36)
            l -= 36;

        double a = (lat + 90.0) / 5.0;
        int i = int(floor(a));
        double x = a - i;

        //Subgrid Extraction from MODIP Array
        double z[4][4] = { 0 };
        for (int k = 0; k < 4; k++)
            for (int j = 0; j < 4; j++)
                z[j][k] = modipArr[i + j][l + k];

        double zk[4] = { 0 };
        for (int k = 0; k < 4; k++)
            zk[k] = interpolate3D(z[0][k], z[1][k], z[2][k], z[3][k], x);

        // MODIP[deg]
        mu = interpolate3D(zk[0], zk[1], zk[2], zk[3], y);
    }

    return mu;
}

double calcAz(double *ai, double mu) {
    double Az;
    if ((ai[0] == 0.0) && (ai[1] == 0.0) && (ai[2] == 0.0))
        Az = 63.7;
    else {
        Az = ai[0] + (ai[1] * mu) + (ai[2] * pow(mu, 2));
        if (Az < 0)
            Az = 0;
        else if (Az > 400)
            Az = 400;
    }
    return Az;
}

double calcSolaR2Declination(int mth, int UT) {
    double dy = (30.5 * mth) - 15;
    double t = dy + ((18 - double(UT)) / 24);
    double am = ((0.9856 * t) - 3.289) * D2R;
    double al = am + ((282.634 + (1.916 * sin(am)) + (0.020 * sin(2 * am))) * D2R);
    return (0.39782 * sin(al));
}

double localTime(double lambda, int UT) {
    lambda = lambda * R2D;
    double LT = double(UT) + (lambda / 15);
    if (LT < 0)
        LT += 24.0;
    else if (LT>=24.0)
        LT -= 24.0;
    return LT;
}

double calcSolarZenithAngle(double phi, double LT, double sinDelta, double cosDelta) {
    double cosZenith = (sin(phi) * sinDelta) + (cos(phi) * cosDelta * cos((12 - LT) * GNSS_PI / 12));
    return R2D * atan2(sqrt(1 - pow(cosZenith, 2)), cosZenith);
}

double calcfoE(double phi, double Az, double effectiveZenith, int mth) {
    int seas;
    if (mth == 1 || mth == 2 || mth == 11 || mth == 12)
        seas = -1;
    else if (mth == 3 || mth == 4 || mth == 9 || mth == 10)
        seas = 0;
    else // mth == 5,6,7,8
        seas = -1;

    double ee = exp(0.3 * phi * 180.0 / GNSS_PI);
    double seasp = seas * (ee - 1) / (ee + 1);
    return sqrt((pow(1.112 - (0.019 * seasp), 2) * sqrt(Az) * pow(cos(effectiveZenith * D2R), 0.6)) + 0.49);
}

void interpolatedAzR(double AzR, double *AF2, double *Am3, double *F2_1, double *F2_2, double *Fm3_1, double* Fm3_2) {
    double AzR_divided = AzR / 100.0;
    double tempF21, tempF22;
    for (int j = 0; j < NEQUICK_G_F2_MAX_INDEX_J; j++)
        for (int k = 0; k < NEQUICK_G_F2_MAX_INDEX_K; k++) {
            tempF21 = *(F2_1 + k + j * NEQUICK_G_F2_MAX_INDEX_K);
            tempF22 = *(F2_2 + k + j * NEQUICK_G_F2_MAX_INDEX_K);
            AF2[k + j * NEQUICK_G_F2_MAX_INDEX_K] = (tempF21 * (1 - AzR_divided)) + (tempF22 * AzR_divided);
        }


    for (int j = 0; j < NEQUICK_G_FM3_MAX_INDEX_J; j++)
        for (int k = 0; k < NEQUICK_G_FM3_MAX_INDEX_K; k++)
            Am3[k + j * NEQUICK_G_FM3_MAX_INDEX_K] = (Fm3_1[k + j * NEQUICK_G_FM3_MAX_INDEX_K] * (1 - AzR_divided)) + (Fm3_2[k + j * NEQUICK_G_FM3_MAX_INDEX_K] * AzR_divided);

    return;
}

void foF1_NmF1(double foE, double foF2, double *foF1, double *NmF1) {
    double limit = 1e-6;
    if (foE >= 2.0)
        *foF1 = 1.4 * foE;
    else
        *foF1 = 0;

    if (abs(*foF1 - foF2) < limit)
        *foF1 = *foF1 * 0.85;

    if (*foF1 < 1e-6)
        *foF1 = 0;

    if ((*foF1 <= 0) && (foE > 2))
        *NmF1 = 0.124 * pow(foE + 0.5, 2); //Eq. 38
    else
        *NmF1 = 0.124 * pow(*foF1, 2); //Eq. 39
    return;
}

void calcfoF2(double mu, double llh[3], double* CF2, double* Cm3, double *foF2, double *M3000F2) {
    double foF2_1;
    double m_k[12];
    double p_n[8];
    double s_n[8];
    double c_n[8];
    int n;
    double aux;

    int k = 1;
    m_k[k - 1] = 1;   //Eq. 56
    for (k = 2; k <= 12; k++)
        m_k[k - 1] = pow(sin(mu * D2R), k - 1); //Eq. 57

    for (n = 2; n <= 9; n++) {
        p_n[n - 2] = pow(cos(llh[0]), n - 1); //Eq. 58
        s_n[n - 2] = sin((n - 1) * llh[1]); //Eq. 59
        c_n[n - 2] = cos((n - 1) * llh[1]); //Eq. 60
    }

    foF2_1 = 0;
    for (k = 1; k <= 12; k++)
        foF2_1 += CF2[k - 1] * m_k[k - 1]; //Eq. 61

    int Q[9] = { 12,12,9,5,2,1,1,1,1 };  //Eq. 62-63

    int K[9];                 //Eq. 64
    K[1 - 1] = (-1) * Q[1 - 1];      //Eq. 65
    for (n = 2; n <= 9; n++)
        K[n - 1] = K[n - 2] + (2 * Q[n - 2]); //Eq. 66

    double foF2_n[9];
    foF2_n[0] = foF2_1;
    for (n = 2; n <= 9; n++) {
        aux = 0;
        for (k = 1; k <= Q[n - 1]; k++) //Eq. 67
            aux += ((CF2[K[n - 1] + (2 * k) - 1 - 1] * c_n[n - 2]) + (CF2[K[n - 1] + (2 * k) - 1] * s_n[n - 2])) * m_k[k - 1] * p_n[n - 2];
        foF2_n[n - 1] = aux;
    }

    *foF2 = 0;
    for (n = 1; n <= 9; n++)
        *foF2 += foF2_n[n - 1];  //Eq. 68

    double M3000F2_n[7] = { 0 };
    for (k = 1; k <= 7; k++)
        M3000F2_n[0] += Cm3[k - 1] * m_k[k - 1];  //Eq. 69

    int R[7] = { 7, 8, 6, 3, 2, 1, 1 }; //Eq. 71
    int H[7];                 //Eq. 72
    H[1 - 1] = (-1) * R[1 - 1];      //Eq. 73
    for (n = 2; n <= 7; n++)
        H[n - 1] = H[n - 2] + (2 * R[n - 2]); //Eq. 74

    for (n = 2; n <= 7; n++) {
        aux = 0;
        for (k = 1; k <= R[n - 1]; k++) //Eq. 75
            aux += ((Cm3[H[n - 1] + (2 * k) - 1 - 1] * c_n[n - 2]) + (Cm3[H[n - 1] + (2 * k) - 1] * s_n[n - 2])) * m_k[k - 1] * p_n[n - 2];
        M3000F2_n[n - 1] = aux;
    }


    *M3000F2 = 0;
    for (n = 1; n <= 7; n++)
        *M3000F2 += M3000F2_n[n - 1];  //Eq. 76

    return;
}

void fourierSeriesF2(int UT, double* AF2, double* Am3, double *CF2, double *Cm3) {
    double T = D2R * ((15 * double(UT)) - 180);
    double aux;
    for (int i = 1; i <= 76; i++) {
        aux = 0;
        for (int k = 1; k <= 6; k++)
            aux += AF2[((i - 1) * 13) + (2 * k) - 1] * sin(T * k) + AF2[((i - 1) * 13) + (2 * k)] * cos(T * k);

        CF2[i - 1] = AF2[((i - 1) * 13)] + aux;
    }

    for (int i = 1; i <= 49; i++) {
        aux = 0;
        for (int k = 1; k <= 4; k++)
            aux += Am3[((i - 1) * 9) + (2 * k) - 1] * sin(T * k) + Am3[((i - 1) * 9) + (2 * k)] * cos(T * k);

        Cm3[i - 1] = Am3[((i - 1) * 9)] + aux;
    }

    return;
}

double calc_hmF2(double foE, double foF2, double M3000F2) {
    double foF2_div_foE = foF2 / foE;
    double exponent = exp(20 * (foF2_div_foE - 1.75));
    double p = ((foF2_div_foE * exponent) + 1.75) / (exponent + 1); //Eq. 84

    double deltaM;
    if (foE < 1e-30)
        deltaM = -0.012;
    else
        deltaM = (0.253 / (p - 1.215)) - 0.012;
    double Msquared = pow(M3000F2, 2);

    return ((1490 * M3000F2 * sqrt(((0.0196 * Msquared) + 1) / ((1.2967 * Msquared) - 1))) / (M3000F2 + deltaM)) - 176;
}

void A2_A3(double NmE, double NmF1, double A1, double hmF2, double hmF1, double hmE, double BEtop, double B1bot, double B2bot, double foF1, double *A2, double *A3) {
    double A3a, A2a;

    if (foF1 < 0.5) {
        *A2 = 0;
        *A3 = 4.0 * (NmE - Epst(A1, hmF2, B2bot, hmE));
    }
    else {
        A3a = 4.0 * NmE;
        double exponent;
        for (int k = 1; k <= 5; k++) {
            A2a = 4.0 * (NmF1 - Epst(A1, hmF2, B2bot, hmF1) - Epst(A3a, hmE, BEtop, hmF1));
            exponent = exp(A2a - (0.8 * NmF1));
            A2a = ((A2a * exponent) + (0.8 * NmF1)) / (1 + exponent);
            A3a = 4.0 * (NmE - Epst(A2a, hmF1, B1bot, hmE) - Epst(A1, hmF2, B2bot, hmE));
        }

        *A2 = A2a;
        exponent = exp(60 * (A3a - 0.005));
        *A3 = ((A3a * exponent) + 0.05) / (1 + exponent);
    }

    return;

}

double calcShapeK(int mth, double NmF2, double hmF2, double B2bot, double AzR) {
    double ka, kb;
    double exponent;
    if (mth == 4 || mth == 5 || mth == 6 || mth == 7 || mth == 8 || mth == 9)
        ka = 6.705 - (0.014 * AzR) - (0.008 * hmF2);
    else
        ka = -7.77 + (0.097 * pow(hmF2 / B2bot, 2)) + (0.153 * NmF2);

    exponent = exp(ka - 2);
    kb = ((ka * exponent) + 2) / (1 + exponent);

    exponent = exp(kb - 8);
    return (((8 * exponent) + kb) / (1 + exponent));
}

double calcH0(double B2bot, double k) {
    double Ha = k * B2bot;
    double x = (Ha - 150) / 100;
    double v = (0.041163 * x - 0.183981) * x + 1.424472;
    double H0 = Ha / v;
    return H0;
}

double bottomSideElecDens(double h, double *A, double hmF2, double hmF1, double hmE, double B2bot, double B1top, double B1bot, double BEtop, double BEbot) {
    double N = 0.0;
    double BE, BF1;
    if (h <= hmE)
        BE = BEbot;
    else
        BE = BEtop;

    if (h <= hmF1)
        BF1 = B1bot;
    else
        BF1 = B1top;

    double alpha[3], s[3];
    double exponent;
    if (h < 100.0) {
        double ds[3];
        double fraction;
        exponent = exp(10 / (1 + abs(100.0 - hmF2)));
        alpha[0] = (100.0 - hmF2) / B2bot;
        alpha[1] = ((100.0 - hmF1) / BF1) * exponent;
        alpha[2] = ((100.0 - hmE) / BE) * exponent;
        for (int i = 0; i < 3; i++) {
            if (abs(alpha[i]) > 25) {
                s[i] = 0;
                ds[i] = 0;
            }
            else {
                s[i] = A[i] * exp(alpha[i]) / pow(1 + exp(alpha[i]), 2);
                fraction = (1 - exp(alpha[i])) / (1 + exp(alpha[i]));
                if (i == 0)
                    ds[i] = fraction / B2bot;
                else if (i == 1)
                    ds[i] = fraction / BF1;
                else
                    ds[i] = fraction / BE;
            }
        }
        double sum_ds_Times_s = 0;
        double sum_s = 0;
        for (int i = 0; i < 3; i++) {
            sum_ds_Times_s += ds[i] * s[i];
            sum_s += s[i];
        }


        double BC = 1 - (10 * sum_ds_Times_s / sum_s);
        double z = (h - 100) / 10;
        N = sum_s * exp(1 - (BC * z) - exp(z * (-1))) * 1e11;
    }
    else {
        exponent = exp(10.0 / (1.0 + abs(h - hmF2)));
        alpha[0] = (h - hmF2) / B2bot;
        alpha[1] = ((h - hmF1) / BF1) * exponent;
        alpha[2] = ((h - hmE) / BE) * exponent;
        double sum_s = 0;
        for (int i = 0; i < 3; i++) {
            if (abs(alpha[i]) > 25)
                s[i] = 0;
            else
                s[i] = A[i] * exp(alpha[i]) / pow(1 + exp(alpha[i]), 2);
            sum_s += s[i];
        }
        N = sum_s * 1e11;
    }

    return N;
}

double topSideElecDens(double h, double NmF2, double hmF2, double H0) {
    double N = 0.0;
    double g = 0.125;
    double r = 100.0;
    double delta_h = h - hmF2;
    double z = delta_h / (H0 * (1 + ((r * g * delta_h) / ((r * H0) + (g * delta_h)))));
    double ea = exp(z);
    if (ea > 1e11)
        N = 1e11 * 4 * NmF2 / ea;
    else
        N = 1e11 * 4 * NmF2 * ea / pow(1 + ea, 2);
    return N;
}

double ElecDens(double s, nequickG_input data, double* sLLH) {
    double N = 1.0;
    double mu = calcMODIP(sLLH);

    double AF2[76][13], Am3[49][9];
    interpolatedAzR(data.Az_R, (double*)AF2, (double*)Am3, data.F2_1, data.F2_2, data.Fm3_1, data.Fm3_2);

    double CF2[76], Cm3[49];
    fourierSeriesF2(data.UT, (double*)AF2, (double*)Am3, CF2, Cm3);

    double foF2, M3000F2;
    calcfoF2(mu, sLLH, CF2, Cm3, &foF2, &M3000F2);
    double NmF2 = 0.124 * pow(foF2, 2);

    double sin_deltaSun, cos_deltaSun;
    sin_deltaSun = calcSolaR2Declination(data.mth, data.UT);
    cos_deltaSun = sqrt(1 - pow(sin_deltaSun, 2));

    double LT = localTime(sLLH[1], data.UT);

    double solarZenithAngle = calcSolarZenithAngle(sLLH[0], LT, sin_deltaSun, cos_deltaSun);
    double expChi = exp(12 * (solarZenithAngle - NEQUICK_G_ZENITH0));
    if (expChi > 1e306)
        expChi = 1e306;
    double effectiveSolarZenithAngle = (solarZenithAngle + (90 - (0.24 * exp(20 - (0.2 * solarZenithAngle)))) * expChi) / (1 + expChi);

    double foE = calcfoE(sLLH[0], data.Az, effectiveSolarZenithAngle, data.mth);
    double NmE = 0.124 * pow(foE, 2);

    double hmF2 = calc_hmF2(foE, foF2, M3000F2);

    if (sLLH[2] <= hmF2) {
        double hmE = 120;  //Eq. 78
        double hmF1 = (hmF2 + hmE) / 2;

        double foF1, NmF1;
        foF1_NmF1(foE, foF2, &foF1, &NmF1);

        double B2bot = (0.385 * NmF2) / (0.01 * exp((-3.467) + (0.857 * log(pow(foF2, 2))) + (2.02 * log(M3000F2))));
        double BEbot = 5;
        double B1bot = (hmF1 - hmE) / 2;  //Eq.87
        double BEtop = fmax(B1bot, 7);
        double B1top = 0.3 * (hmF2 - hmF1);

        double Ai[3] = { 0.0 };
        Ai[0] = 4 * NmF2; //Eq.90
        A2_A3(NmE, NmF1, Ai[0], hmF2, hmF1, hmE, BEtop, B1bot, B2bot, foF1, &Ai[1], &Ai[2]);

        N = bottomSideElecDens(sLLH[2], Ai, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot);
    }
    else {
        double B2bot = (0.385 * NmF2) / (0.01 * exp((-3.467) + (0.857 * log(pow(foF2, 2))) + (2.02 * log(M3000F2))));
        double k = calcShapeK(data.mth, NmF2, hmF2, B2bot, data.Az_R);
        double H0 = calcH0(B2bot, k);
        N = topSideElecDens(sLLH[2], NmF2, hmF2, H0);
    }
    return N;
}

double insertToElecDens(double s, perigee_t computedPerigee, nequickG_input data, double* llh) {
    double sLLH[3];
    sLLH[0] = llh[0];
    sLLH[1] = llh[1];
    if (computedPerigee.isVertical) sLLH[2] = s;
    else sLLH[2] = sqrt(pow(s, 2) + pow(computedPerigee.rp, 2)) - NEQUICK_G_RE_KM;
    return ElecDens(s, data, sLLH);
}

void TECGaussKronrod(double heightP1, double heightP2, perigee_t perigee, nequickG_input data, double tol, double *recursion_level, double *result, double* userLLH, double* satLLH, double *currentLLH) {
    *result = 0.0;

    double mid_point = (heightP1 + heightP2) / 2.0;
    double half_diff = (heightP2 - heightP1) / 2.0;

    double K15_integration = 0.0;
    double G7_integration = 0.0;
    int G7_index = 0;

    int i;
    for (i = 0; i < NEQUICK_G_KRONROD_K15_POINT_COUNT; i++) {

        double height_km = mid_point + (half_diff * xi[i]);

        if (!perigee.isVertical)
            getCurrentPos(perigee, height_km, userLLH, currentLLH);
        else
            getCurrentPosVertical(perigee, userLLH, satLLH, currentLLH);

        double TEC = insertToElecDens(height_km, perigee, data, currentLLH);
        K15_integration += (TEC * wi[i]);

        if (i % 2 == 1) {
            G7_integration += (TEC * wig[G7_index]);
            G7_index++;
        }
    }

    K15_integration *= half_diff;
    G7_integration *= half_diff;

    bool errorWithinTolerance = (fabs((K15_integration - G7_integration) / K15_integration) <= tol) ||
        (fabs(K15_integration - G7_integration) <= tol);
    if (errorWithinTolerance || (*recursion_level) == NEQUICK_G_KRONROD_MAX_RECURSION_LEVEL) {
        *result = K15_integration;
        return;
    }

    // Error not acceptable
    // split into two parts to improve accuracy
    // and try again.
    (*recursion_level)++;
    TECGaussKronrod(heightP1, heightP1 + half_diff, perigee, data, tol, recursion_level, result, userLLH, satLLH, currentLLH);
    double secondHalfResult = 0.0;
    TECGaussKronrod(heightP1 + half_diff, heightP2, perigee, data, tol, recursion_level, &secondHalfResult, userLLH, satLLH, currentLLH);
    *result += secondHalfResult;

    (*recursion_level)--;
    return;
}

double NeQuickG(double* userLLH, double* satelliteLLH, nequickG_input modelData, bool &invalidFlag) {
    double TEC = 0.0;

    double currentLLH[3] = { 0 };

    perigee_t rayPerigee;
    rayPerigee = calcPerigee(userLLH, satelliteLLH, invalidFlag);
    if (invalidFlag)
        return 0.0;

    // Radii of objects [km]
    double r1 = userLLH[2] + NEQUICK_G_RE_KM;
    double r2 = satelliteLLH[2] + NEQUICK_G_RE_KM;

    //Intermediate Points [km]
    double s1 = sqrt(pow(r1, 2) - pow(rayPerigee.rp, 2));
    double s2 = sqrt(pow(r2, 2) - pow(rayPerigee.rp, 2));
    double recursion_level = 0.0;

    bool badPos = (satelliteLLH[2] <= 1000.0) || ((satelliteLLH[2] <= 2000.0) && (userLLH[2] >= 1000.0)) || ((satelliteLLH[2] <= 2000.0) && (userLLH[2] < 1000.0));
    if (badPos) {
        invalidFlag = true;
        return 0.0;
    }

    if (userLLH[2] >= 2000.0) {
        if (rayPerigee.isVertical) { s1 = userLLH[2]; s2 = satelliteLLH[2]; }
        TECGaussKronrod(s1, s2, rayPerigee, modelData, 0.01, &recursion_level, &TEC, userLLH, satelliteLLH, currentLLH);
    }
    else if (userLLH[2] >= 1000.0) {
        double sb;
        if (rayPerigee.isVertical) { s1 = userLLH[2]; s2 = satelliteLLH[2]; sb = 2000.0; }
        else sb = sqrt(70076989.44 - pow(rayPerigee.rp, 2));

        double TEC1{}, TEC2{};
        TECGaussKronrod(s1, sb, rayPerigee, modelData, 0.01, &recursion_level, &TEC1, userLLH, satelliteLLH, currentLLH);
        TECGaussKronrod(sb, s2, rayPerigee, modelData, 0.01, &recursion_level, &TEC2, userLLH, satelliteLLH, currentLLH);
        TEC = TEC1 + TEC2;
    }
    else {
        double sa,sb;
        if (rayPerigee.isVertical) { s1 = userLLH[2]; s2 = satelliteLLH[2]; sa = 1000.0; sb = 2000.0; }
        else { sa = sqrt(54334589.44 - pow(rayPerigee.rp, 2)); sb = sqrt(70076989.44 - pow(rayPerigee.rp, 2)); }

        // Integral Segmentation
        double TEC1{}, TEC2{}, TEC3{};
        TECGaussKronrod(s1, sa, rayPerigee, modelData, 0.001, &recursion_level, &TEC1, userLLH, satelliteLLH, currentLLH);
        recursion_level = 0;
        TECGaussKronrod(sa, sb, rayPerigee, modelData, 0.01, &recursion_level, &TEC2, userLLH, satelliteLLH, currentLLH);
        recursion_level = 0;
        TECGaussKronrod(sb, s2, rayPerigee, modelData, 0.01, &recursion_level, &TEC3, userLLH, satelliteLLH, currentLLH);
        TEC = TEC1 + TEC2 + TEC3;
    }

    return (TEC * 1e-13);
}

void choose_F2_Fm(double*& ptrF2_1, double*& ptrF2_2, double*& ptrFm3_1, double*& ptrFm3_2, int month) {
	switch (month) {
	case 1:
		ptrF2_1 = (double*)F2_1_1;
		ptrF2_2 = (double*)F2_2_1;
		ptrFm3_1 = (double*)Fm3_1_1;
		ptrFm3_2 = (double*)Fm3_2_1;
		break;
	case 2:
		ptrF2_1 = (double*)F2_1_2;
		ptrF2_2 = (double*)F2_2_2;
		ptrFm3_1 = (double*)Fm3_1_2;
		ptrFm3_2 = (double*)Fm3_2_2;
		break;
	case 3:
		ptrF2_1 = (double*)F2_1_3;
		ptrF2_2 = (double*)F2_2_3;
		ptrFm3_1 = (double*)Fm3_1_3;
		ptrFm3_2 = (double*)Fm3_2_3;
		break;
	case 4:
		ptrF2_1 = (double*)F2_1_4;
		ptrF2_2 = (double*)F2_2_4;
		ptrFm3_1 = (double*)Fm3_1_4;
		ptrFm3_2 = (double*)Fm3_2_4;
		break;
	case 5:
		ptrF2_1 = (double*)F2_1_5;
		ptrF2_2 = (double*)F2_2_5;
		ptrFm3_1 = (double*)Fm3_1_5;
		ptrFm3_2 = (double*)Fm3_2_5;
		break;
	case 6:
		ptrF2_1 = (double*)F2_1_6;
		ptrF2_2 = (double*)F2_2_6;
		ptrFm3_1 = (double*)Fm3_1_6;
		ptrFm3_2 = (double*)Fm3_2_6;
		break;
	case 7:
		ptrF2_1 = (double*)F2_1_7;
		ptrF2_2 = (double*)F2_2_7;
		ptrFm3_1 = (double*)Fm3_1_7;
		ptrFm3_2 = (double*)Fm3_2_7;
		break;
	case 8:
		ptrF2_1 = (double*)F2_1_8;
		ptrF2_2 = (double*)F2_2_8;
		ptrFm3_1 = (double*)Fm3_1_8;
		ptrFm3_2 = (double*)Fm3_2_8;
		break;
	case 9:
		ptrF2_1 = (double*)F2_1_9;
		ptrF2_2 = (double*)F2_2_9;
		ptrFm3_1 = (double*)Fm3_1_9;
		ptrFm3_2 = (double*)Fm3_2_9;
		break;
	case 10:
		ptrF2_1 = (double*)F2_1_10;
		ptrF2_2 = (double*)F2_2_10;
		ptrFm3_1 = (double*)Fm3_1_10;
		ptrFm3_2 = (double*)Fm3_2_10;
		break;
	case 11:
		ptrF2_1 = (double*)F2_1_11;
		ptrF2_2 = (double*)F2_2_11;
		ptrFm3_1 = (double*)Fm3_1_11;
		ptrFm3_2 = (double*)Fm3_2_11;
		break;
	case 12:
		ptrF2_1 = (double*)F2_1_12;
		ptrF2_2 = (double*)F2_2_12;
		ptrFm3_1 = (double*)Fm3_1_12;
		ptrFm3_2 = (double*)Fm3_2_12;
		break;
	}
	return;
}
