#include "../include/galileo-sdr.h"

/*
WARNING - HARDCODE CHANNEL TO PRN MAP FOR STREAM TX AS WELL AS STREAM RX - THERE IS NO CHANNEL TO PRN CHECK AS OF NOW
******************** TODO **********************
*/
using namespace std;

#define MAX_CHAN (16)
#define INCOMING_SIZE 9 // Define here how many channels to receive from streamer (+1 for TOW)

int sockinit(short port)
{
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    struct sockaddr_in servaddr;
    memset(&servaddr, 0, sizeof(servaddr));
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(port);
    servaddr.sin_addr.s_addr = inet_addr("127.0.0.1");
    if (connect(sock, (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
    {
        perror("connect");
        exit(1);
    }

    return sock;
}

void sockclose(int s)
{
    close(s);
}

void socksend(int s, void *data, int siz)
{
    send(s, data, siz, 0);
}

long int timem()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec * 1000 + t.tv_usec / 1000;
}

int udpinit(short port)
{
    int sock = socket(AF_INET, SOCK_DGRAM, 0);
    struct sockaddr_in servaddr;
    memset(&servaddr, 0, sizeof(servaddr));
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(port);
    servaddr.sin_addr.s_addr = INADDR_ANY;
    if (::bind(sock, (struct sockaddr *)&servaddr, sizeof(servaddr)) < 0)
    {
        perror("Bind");
        exit(1);
    }
    return sock;
}

int udprecv(int s, void *data, int size)
{
    struct sockaddr from;
    return recvfrom(s, data, size, 0, &from, 0);
}

double raw_bits[INCOMING_SIZE];
double llhr[3] = {42.34937020, -71.08817058, 100};

double dynamic_dt = 0;
double tow_correction;
bool tow_rx;
bool tow_fix = false;
bool use_hash = false;

static long get_nanos2(void)
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}

void *bitstreamer_thread(vector<queue<int>> &queues, map<int, int> &sm, bool *exit_flag)
{
    fprintf(stderr, "\nListening for bit stream\n");

    fprintf(stderr, "Initializing %d channels with 0 bit - %d\n", MAX_CHAN, (int)*exit_flag);

    int s = udpinit(7531);
    tow_rx = false;
    tow_fix = false;
    char chan_str[20];
    char sat_prn_str[20];
    int sat_prn = 0;
    int bit = 0;
    int assigned_chan = 0;
    int content;
    double a;
    while (1)
    {
        // a = get_nanos2();
        udprecv(s, raw_bits, INCOMING_SIZE * sizeof(double));
        // cerr << "\nTime - " << (get_nanos2() - a)/1e6 <<endl;

        for (int i = 0; i < (INCOMING_SIZE - 1); i++)
        {
            content = (int)raw_bits[i];

            // g_queue_push_tail(queues[i], "1");

            sat_prn = content / 10;
            bit = content % 10;

            // sprintf(chan_str,"%d",i);
            snprintf(sat_prn_str, sizeof(sat_prn_str), "%d", sat_prn);

            if (sm.find(sat_prn) != sm.end())
            {
                // fprintf(stderr, "\nReceived data %d - %d - %d", sat_prn, sm->find(sat_prn)->first, sm->find(sat_prn)->second);
                assigned_chan = sm.find(sat_prn)->second;

                // if (sat_prn == 4)
                //     fprintf(stderr, "\nReceived data %d - %d - %d", sat_prn, bit, assigned_chan);

                if (bit == 1)
                    queues.at(assigned_chan).push(1);
                else if (bit == 0)
                    queues.at(assigned_chan).push(-1);
                else
                    queues.at(assigned_chan).push(0);

                // if (sat_prn == 25)
                //     fprintf(stderr, "\nReceived data %d - %d - %d - %d - %d", sat_prn, sm->find(sat_prn)->first, sm->find(sat_prn)->second, assigned_chan, bit);

                // fprintf(stderr,"Chn %d : %d : %d\n", assigned_chan, sat_prn, queues->at(assigned_chan).front());
            }
        }

        if (!tow_rx)
        {
            // udprecv(s, &tow_correction, sizeof(double));
            tow_correction = raw_bits[INCOMING_SIZE - 1] / 1000.0;
            fprintf(stderr, "\nTOW Correction %f: \n", tow_correction);
            tow_fix = true;
            tow_rx = true;
        }
    }
    return NULL;
}

void *dt_thread(void *)
{
    fprintf(stderr, "Listening for run-time range corrections\n");

    int s = udpinit(7532);

    while (1)
    {
        udprecv(s, &dynamic_dt, sizeof(double));
        // printf("DT: %f\n", llhr[0], llhr[1], llhr[2]);
    }
}

void *locations_thread(double llhr_init[3])
{
    fprintf(stderr, "Listening for run-time position updates\n");

    int s = udpinit(7533);

    llhr[0] = llhr_init[0];
    llhr[1] = llhr_init[1];
    llhr[2] = llhr_init[2];

    while (1)
    {
        udprecv(s, llhr, 3 * sizeof(double));
        fprintf(stdout, "Location Update: %f,%f,%f\n", llhr[0], llhr[1], llhr[2]);
    }
}