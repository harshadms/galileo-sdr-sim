#include "galileo-sdr.h"
#include <uhd/usrp/multi_usrp.hpp>
#include <queue>

using namespace std;

struct hash_queue
{
    vector<queue<int>> queue_bits;
    map<int, int> chn_prn_map; // Maps satellites to channels, key is satellite, value is channel
};

typedef struct
{
    pthread_t thread;
    pthread_mutex_t lock;
    // int error;
    uhd::tx_metadata_t md;
    uhd::tx_streamer::sptr stream;
    short *buffer;
    const void **buffer_ptr;
} tx_t;

/*! \brief Structure representing USRP configuration */
typedef struct 
{
    double carr_freq;
    double gain;
    double samp_rate;
    string device_args;
} usrp_conf_t;

typedef struct
{
    pthread_t thread;
    pthread_mutex_t lock;
    // int error;

    int ready;
    pthread_cond_t initialization_done;
} galileo_t;

/*! \brief Structure representing Galileo time */
typedef struct
{
    int week;   /*!< Galileo week number (since August 1999) */
    double sec; /*!< second inside the GPS \a week */
} galtime_t;

/*! \brief Structure representing time */
typedef struct
{                /* time struct */
    time_t time; /* time (s) expressed by standard time_t */
    double sec;  /* fraction of second under 1 s */
} gtime_t;


/*! \brief Structure repreenting UTC time */
typedef struct
{
    int y;      /*!< Calendar year */
    int m;      /*!< Calendar month */
    int d;      /*!< Calendar day */
    int hh;     /*!< Calendar hour */
    int mm;     /*!< Calendar minutes */
    double sec; /*!< Calendar seconds */
} datetime_t;

/*! \brief Structure representing ephemeris of a single satellite */
typedef struct
{
    int vflg; /*!< Valid Flag */
    datetime_t t;
    galtime_t time;
    galtime_t toc; /*!< Time of Clock */
    galtime_t toe; /*!< Time of Ephemeris */
    int iode;      /*!< Isuse of Data, Ephemeris - IODnav*/
    double deltan; /*!< Delta-N (radians/sec) */
    double cuc;    /*!< Cuc (radians) */
    double cus;    /*!< Cus (radians) */
    double cic;    /*!< Correction to inclination cos (radians) */
    double cis;    /*!< Correction to inclination sin (radians) */
    double crc;    /*!< Correction to radius cos (meters) */
    double crs;    /*!< Correction to radius sin (meters) */
    double ecc;    /*!< e Eccentricity */
    double sqrta;  /*!< sqrt(A) (sqrt(m)) */
    double m0;     /*!< Mean anamoly (radians) */
    double omg0;   /*!< Longitude of the ascending node (radians) */
    double inc0;   /*!< Inclination (radians) */
    double aop;
    double omgdot; /*!< Omega dot (radians/s) */
    double idot;   /*!< IDOT (radians/s) */
    double af0;    /*!< Clock offset (seconds) */
    double af1;    /*!< rate (sec/sec) */
    double af2;    /*!< acceleration (sec/sec^2) */
    double bgde5a; /*!< Group delay E5 bias */
    double bgde5b; // = parameters.at(25);	 // BGD_E5b
    double tgd_ext[5];
    double gps_time;
    short flag;
    short ura;

    int svhlth;
    int codeL2;
    // Working variables follow
    double n;       /*!< Mean motion (Average angular velocity) */
    double sq1e2;   /*!< sqrt(1-e^2) */
    double A;       /*!< Semi-major axis */
    double omgkdot; /*!< OmegaDot-OmegaEdot */
    double omg_t;
    int svid;
    int PRN;
} ephem_t;

/*! \brief Structure representing Iono and UTC models */
typedef struct // modified
{
    int enable;
    int vflg;
    double alpha0, alpha1, alpha2, alpha3;
    double beta0, beta1, beta2, beta3;
    double A0, A1, A2;
    int dtls, tot, wnt;
    int dtlsf, dn, wnlsf;
    double ai0, ai1, ai2, ai3; // added
} ionoutc_t;

/*! \brief Structure representing pseudorange */
typedef struct
{
    galtime_t g;
    double range; // pseudorange
    double rate;
    double d; // geometric distance
    double azel[2];
    double iono_delay;
} range_t;

/*! \brief Structure representing a Channel */
typedef struct
{
    int prn;       /*< PRN Number */
    short *ca_E1B; /*< PRN Sequence buffer assign at runtime using malloc*/
    short *ca_E1C; /*< PRN Sequence buffer assign at runtime using malloc*/
    double f_carr; /*< Carrier frequency */
    double f_code; /*< Code frequency */
    int *page;
    double carr_phase;
    double code_phase_p;                 /*< Code phase */
    double code_phase_s;                 /*< Code phase */
    galtime_t g0;                        /*!< GPS time at start */
    int ipage;                           /*!< initial word */
    int ibit;                            /*!< initial bit */
    int icode;                           /*!< initial code */
    int dataBit;                         /*!< current data bit */
    float E1B_chip;                      /*!< current C/A code */
    float E1C_chip;
    double azel[2];
    range_t rho0;
    bool set_code_phase;
    double code_phase;
} channel_t;

typedef struct
{
    char navfile[MAX_CHAR];
    char umfile[MAX_CHAR];
    char outfile[MAX_CHAR];
    char tvfile[MAX_CHAR];
    
    int staticLocationMode;
    int nmeaGGA;
    int iduration;
    int verb;
    galtime_t g0;
    double llh[3];
    int interactive;
    int timeoverwrite;
    int iono_enable;
    bool use_usrp;
    bool use_bit_stream;
} option_t;

typedef struct
{
    option_t opt;

    tx_t tx;
    galileo_t galileo_sim;

    int status;
    short udp_port;
    bool finished;
    short *fifo;
    long head, tail;
    size_t sample_length;

    pthread_cond_t fifo_read_ready;
    pthread_cond_t fifo_write_ready;

    double time;
} sim_t;
