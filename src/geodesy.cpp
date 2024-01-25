#include "../include/galileo-sdr.h"

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

    return;
}

/*! \brief Calculates ionospheric delay
 *  \param[out] iono_delay
 *  \param[in] ionoutc IonoUTC model
 *  \param[in] g current time
 *  \param[in] llh LLH coordinate
 *  \param[in] azel azel
 */
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


/*! \brief Compute dot-product of two vectors
 *  \param[in] x1 First multiplicand
 *  \param[in] x2 Second multiplicand
 *  \returns Dot-product of both multiplicands
 */
double dotProd(const double *x1, const double *x2)
{
    return (x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
}

/*! \brief Checks if a particular satellite is visible
 *  \param[in] eph Current ephemeris
 *  \param[in] g galileo time
 *  \param[in] xyz ecef coords
 *  \param[in] elvMask elevation mask
 *  \param[in] azel azel
 *  \param[in] prn Satellite prn
 *  \returns 0/1 visibility
 */
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

