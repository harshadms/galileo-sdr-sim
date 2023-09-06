#include "galileo-sdr.h"
#include <uhd/usrp/multi_usrp.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TX_FREQUENCY 1575420000
#define TX_SAMPLERATE SAMP_RATE
#define TX_BANDWIDTH SAMP_RATE*2

#define NUM_BUFFERS 32
#define SAMPLES_PER_BUFFER (32 * 1024)
#define NUM_TRANSFERS 16
#define TIMEOUT_MS 1000

#define NUM_IQ_SAMPLES (TX_SAMPLERATE / 10)
#define FIFO_LENGTH (NUM_IQ_SAMPLES * 2)

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

typedef struct
{
    pthread_t thread;
    pthread_mutex_t lock;
    // int error;

    int ready;
    pthread_cond_t initialization_done;
} galileo_t;

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

extern void *galileo_task(void *arg);
extern int is_fifo_write_ready(sim_t *s);