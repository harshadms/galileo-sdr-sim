#include "../include/galileo-sdr.h"

size_t get_sample_length(sim_t *s)
{
    long length;

    length = s->head - s->tail;
    if (length < 0)
        length += FIFO_LENGTH;

    return ((size_t)length);
}

size_t fifo_read(short *buffer, size_t samples, sim_t *s)
{
    size_t length;
    size_t samples_remaining;
    short *buffer_current = buffer;

    length = get_sample_length(s);

    if (length < samples)
        samples = length;

    length = samples; // return value

    samples_remaining = FIFO_LENGTH - s->tail;

    if (samples > samples_remaining)
    {
        memcpy(buffer_current, &(s->fifo[s->tail * 2]),
               samples_remaining * sizeof(short) * 2);
        s->tail = 0;
        buffer_current += samples_remaining * 2;
        samples -= samples_remaining;
    }

    memcpy(buffer_current, &(s->fifo[s->tail * 2]),
           samples * sizeof(short) * 2);
           
    s->tail += (long)samples;
    if (s->tail >= FIFO_LENGTH)
        s->tail -= FIFO_LENGTH;

    return (length);
}

bool is_finished_generation(sim_t *s) { return s->finished; }

int is_fifo_write_ready(sim_t *s)
{
    int status = 0;

    s->sample_length = get_sample_length(s);
    // fprintf(stderr, "\nFIFO %d - %d", s->sample_length, NUM_IQ_SAMPLES);

    if (s->sample_length < NUM_IQ_SAMPLES)
        status = 1;

    return (status);
}

