#define _CRT_SECURE_NO_WARNINGS

#include "../include/galileo-sdr.h"
#include <math.h>
#include <unistd.h>
#include <vector>

#define EXECUTE_OR_GOTO(label, ...) \
    if (__VA_ARGS__)                \
    {                               \
        return_code = EXIT_FAILURE; \
        goto label;                 \
    }

using namespace std;

bool stop_signal_called = false;

void sigint_handler(int code)
{
    (void)code;
    stop_signal_called = true;
    fprintf(stderr, "Done\n");
}

long int samples_consumed = 0l;

void init_sim(sim_t *s)
{
    printf("\nGalileo SIM initiation started.");

    pthread_mutex_init(&(s->tx.lock), NULL);
    // s->tx.error = 0;

    pthread_mutex_init(&(s->galileo_sim.lock), NULL);
    // s->galileo_sim.error = 0;
    s->galileo_sim.ready = 0;

    pthread_cond_init(&(s->galileo_sim.initialization_done), NULL);

    s->status = 0;
    s->head = 0;
    s->tail = 0;
    s->sample_length = 0;

    pthread_cond_init(&(s->fifo_write_ready), NULL);

    pthread_cond_init(&(s->fifo_read_ready), NULL);

    s->time = 0.0;

    printf("\nGalileo SIM initiated.");
}

void *tx_task(void *arg)
{
    sim_t *s = (sim_t *)arg;
    size_t samples_populated;
    size_t num_samps_sent = 0;
    size_t samples_per_buffer = SAMPLES_PER_BUFFER;
    fprintf(stderr, "\nTx task");
    // sleep(1);
    int k = 0;

    while (1)
    {
        if (stop_signal_called)
            goto out;

        short *tx_buffer_current = s->tx.buffer;
        int buffer_samples_remaining = SAMPLES_PER_BUFFER;
        const void **buffs_ptr = NULL;

        while (buffer_samples_remaining > 0)
        {
            pthread_mutex_lock(&(s->galileo_sim.lock));

            while (get_sample_length(s) == 0)
            {
                pthread_cond_wait(&(s->fifo_read_ready), &(s->galileo_sim.lock));
            }

            samples_populated = fifo_read(tx_buffer_current, buffer_samples_remaining, s);
            pthread_mutex_unlock(&(s->galileo_sim.lock));

            pthread_cond_signal(&(s->fifo_write_ready));

            // Advance the buffer pointer.
            buffer_samples_remaining -= (unsigned int)samples_populated;

            tx_buffer_current += (2 * samples_populated);
        }

        k++;

        buffs_ptr = (const void **)tx_buffer_current;

        // vector<double> tx_buffer_vector; //(tx_buffer_current, tx_buffer_current + SAMPLES_PER_BUFFER);
        //  vector<complex<float>> fc_buffer;
        //  for(int i=0; i < SAMPLES_PER_BUFFER; i=i+2)
        //  {
        //      fc_buffer.push_back(complex<float>(s->tx.buffer[i], s->tx.buffer[i+1]));// * 2], tx_buffer_current[i * 2 + 1]) );
        //      //fprintf(stderr, "\nF - %f/ %d", s->tx.buffer[i], i);
        //  }

        // uhd_tx_streamer_send(s->tx.stream, s->tx.buffer_ptr, SAMPLES_PER_BUFFER, &s->tx.md, 1000, &samples_populated);
        // printf("\nHere %ld\n", tx_buffer_vector.size());

        size_t num_tx_samps = s->tx.stream->send(s->tx.buffer, SAMPLES_PER_BUFFER, s->tx.md, 1000);
        samples_consumed = samples_consumed + num_tx_samps;
        // fprintf(stderr, "\nSent %d", k);
        // tx_buffer_vector.clear();

        if (is_fifo_write_ready(s))
        {
            fprintf(stderr, "\rTime = %4.1f - %d", s->time, k);
            s->time += 0.1;
            fflush(stdout);
        }
        else if (is_finished_generation(s))
        {
            goto out;
        }
    }
out:
    return NULL;
}

int start_tx_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->tx.thread), NULL, tx_task, s);

    return (status);
}

int start_galileo_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->galileo_sim.thread), NULL, galileo_task, s);

    return (status);
}

void usage(char *progname)
{
    printf(
        "Usage: %s [options]\n"
        "Options:\n"
        "  -e <Ephemeris>   RINEX navigation file for Galileo ephemerides (required)\n"
        "  -o <File sink>   File to store IQ samples\n"
        "  -l <location>    Lat,Lon,Hgt (static mode) e.g. 35.274,137.014,100\n"
        "  -t <date,time>   Scenario start time YYYY/MM/DD,hh:mm:ss\n"
        "  -T <date,time>   Overwrite TOC and TOE to scenario start time\n"
        "  -d <duration>    Duration [sec] (max: %.0f)\n"
        "  -a <rf_gain>     Absolute RF gain in [0 ... 60] (default: 30)\n"
        "  -U               Disable USRP (-U 1)\n"
        "  -b               Disable Bit stream (-b 1)\n",
        progname,

        ((double)USER_MOTION_SIZE) / 10.0);

    return;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        usage(argv[0]);
        exit(1);
    }

    // Set default values
    sim_t s;

    s.finished = false;
    s.opt.navfile[0] = 0;
    s.opt.tvfile[0] = 0;
    s.opt.outfile[0] = 0;

    s.opt.g0.week = -1;
    s.opt.g0.sec = 0.0;
    s.opt.iduration = USER_MOTION_SIZE;
    s.opt.verb = FALSE;
    s.opt.llh[0] = 42.3601;
    s.opt.llh[1] = -71.0589;
    s.opt.llh[2] = 2;
    s.opt.interactive = FALSE;
    s.opt.timeoverwrite = FALSE;
    s.opt.iono_enable = TRUE;
    s.udp_port = 5671;
    s.opt.use_usrp = true;
    s.opt.use_bit_stream = true;

    // Options
    int result;
    double duration;
    datetime_t t0;
    double gain = 0;
    std::string device_args = "";
    size_t channel = 1;
    bool verbose = false;
    int return_code = EXIT_SUCCESS;
    char error_string[512];
    double rate = TX_SAMPLERATE;
    double freq = TX_FREQUENCY;
    bool use_usrp = true;
    // Buffer sizes
    size_t samps_per_buff = SAMPLES_PER_BUFFER;
    float *buffer = NULL;
    const void **buffer_ptr = NULL;

    while ((result = getopt(argc, argv, "e:n:o:u:g:l:T:t:d:G:a:p:iI:U:b:v")) != -1)
    {
        switch (result)
        {
        case 'e':
            strcpy(s.opt.navfile, optarg);
            break;

        case 'n':
            strcpy(s.opt.tvfile, optarg);
            break;

        case 'o':
            strcpy(s.opt.outfile, optarg);
            break;

        case 'l':
            sscanf(optarg, "%lf,%lf,%lf", &s.opt.llh[0], &s.opt.llh[1],
                   &s.opt.llh[2]);
            break;

        case 'T':
            s.opt.timeoverwrite = TRUE;
            if (strncmp(optarg, "now", 3) == 0)
            {
                time_t timer;
                struct tm *gmt;

                time(&timer);
                gmt = gmtime(&timer);

                t0.y = gmt->tm_year + 1900;
                t0.m = gmt->tm_mon + 1;
                t0.d = gmt->tm_mday;
                t0.hh = gmt->tm_hour;
                t0.mm = gmt->tm_min;
                t0.sec = (double)gmt->tm_sec;

                date2gal(&t0, &s.opt.g0);

                break;
            }

        case 't':
            sscanf(optarg, "%d/%d/%d,%d:%d:%lf", &t0.y, &t0.m, &t0.d, &t0.hh, &t0.mm,
                   &t0.sec);
            if (t0.y <= 1980 || t0.m < 1 || t0.m > 12 || t0.d < 1 || t0.d > 31 ||
                t0.hh < 0 || t0.hh > 23 || t0.mm < 0 || t0.mm > 59 || t0.sec < 0.0 ||
                t0.sec >= 60.0)
            {
                printf("ERROR: Invalid date and time.\n");
                exit(1);
            }
            t0.sec = floor(t0.sec);
            date2gal(&t0, &s.opt.g0);
            break;

        case 'd':
            duration = atof(optarg);
            s.opt.iduration = (int)(duration * 10.0 + 0.5);
            break;

        case 'G':
            gain = atof(optarg);
            if (gain < 0)
                gain = 0;
            if (gain > 60)
                gain = 60;
            break;

        case 'a':
            device_args = strdup(optarg);
            break;

        case 'p':
            // printf("%s\n", optarg);
            s.udp_port = atoi(optarg);

            break;

        case 'i':
            s.opt.interactive = TRUE;
            break;

        case 'I':
            s.opt.iono_enable = FALSE; // Disable ionospheric correction
            break;

        case 'U':
            s.opt.use_usrp = false;
            use_usrp = false;
            break;

        case 'b':
            s.opt.use_bit_stream = false;
            break;

        case 'v':
            s.opt.verb = true;
            break;

        case ':':

        case '?':
            usage(argv[0]);
            exit(1);

        default:
            break;
        }
    }

    if (s.opt.navfile[0] == 0 && s.opt.tvfile[0] == 0)
    {
        std::cout << s.opt.navfile;
        std::cout << s.opt.tvfile;
        printf("ERROR: Galileo ephemeris/nav_msg file is not specified.\n");
        exit(1);
    }

    if (s.opt.outfile[0] == 0)
    {
        printf("[+] File sink not specified. Using galileosim.bin\n");
        snprintf(s.opt.outfile, sizeof(s.opt.outfile), "galileosim.ishort");
    }

    // Initialize simulator
    init_sim(&s);

    // Allocate FIFOs to hold 0.1 seconds of I/Q samples each.
    s.fifo = (short *)malloc(FIFO_LENGTH * sizeof(short) * 2); // for 32-bit I and Q samples

    if (s.fifo == NULL)
    {
        printf("ERROR: Failed to allocate I/Q sample buffer.\n");
        // goto out;
    }

    if (use_usrp)
    {   
        usrp_conf_t usrp_conf;  
        usrp_conf.carr_freq = GALILEO_E1_FREQ_HZ;
        usrp_conf.gain = gain;
        usrp_conf.device_args = device_args;
        usrp_conf.samp_rate = TX_SAMPLERATE;

        init_usrp(usrp_conf, &s);
    }

    // Start Galileo task.
    s.status = start_galileo_task(&s);
    if (s.status < 0)
    {
        fprintf(stderr, "Failed to start Galileo task.\n");
    }
    else
        printf("\nCreating Galileo task...\n");

    if (use_usrp)
    { // Wait until Galileo task is initialized
        pthread_mutex_lock(&(s.tx.lock));
        while (!s.galileo_sim.ready)
            pthread_cond_wait(&(s.galileo_sim.initialization_done), &(s.tx.lock));

        pthread_mutex_unlock(&(s.tx.lock));

        // Fillfull the FIFO.
        if (is_fifo_write_ready(&s))
            pthread_cond_signal(&(s.fifo_write_ready));

        // // Start TX task
        s.status = start_tx_task(&s);

        if (s.status < 0)
        {
            fprintf(stderr, "Failed to start TX task.\n");
            // goto out;
        }
        else
            printf("Creating TX task...\n");

        // Running...
        printf("Running...\n"
               "Press Ctrl+C to abort.\n");

        // Wainting for TX task to complete.
        pthread_join(s.tx.thread, NULL);
    }
    else
        pthread_join(s.galileo_sim.thread, NULL);

    fprintf(stderr, "\nTotal samples consumed: %ld", samples_consumed);
    printf("\nDone!\n");
}