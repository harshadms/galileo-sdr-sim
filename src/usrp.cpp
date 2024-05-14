#include "../include/galileo-sdr.h"

void init_usrp(usrp_conf_t usrpconf, sim_t *s)
{
    // Create USRP
    size_t samps_per_buff = SAMPLES_PER_BUFFER;

    std::cout << "Creating the usrp device with: " << usrpconf.device_args << std::endl;

    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(usrpconf.device_args);

    // Lock mboard clocks
    usrp->set_clock_source("internal");
    usrp->set_time_source("internal");

    // Wait for some time to setup and finish
    sleep(1);

    fprintf(stderr, "\nSet Clock Source to %s", usrp->get_clock_source(0).c_str());
    fprintf(stderr, "\nSet Time Source to %s", usrp->get_time_source(0).c_str());

    // set the sample rate
    fprintf(stderr, "\nSetting TX Rate: %f Msps\n", TX_SAMPLERATE / 1e6);
    usrp->set_tx_rate(TX_SAMPLERATE);
    fprintf(stderr, "\nActual TX Rate: %f Msps", usrp->get_tx_rate() / 1e6);

    // set the center frequency
    fprintf(stderr, "\nSetting TX Freq: %f MHz", usrpconf.carr_freq / 1e6);
    usrp->set_tx_freq(usrpconf.carr_freq);
    fprintf(stderr, "\nActual TX Freq: %f MHz", usrp->get_tx_freq() / 1e6);

    // set the rf gain
    fprintf(stderr, "\nSetting TX Gain: %f dB", usrpconf.gain);
    usrp->set_tx_gain(usrpconf.gain);
    fprintf(stderr, "\nActual TX Gain: %f dB", usrp->get_tx_gain());

    std::cout << "\nUsing antenna:" << usrp->get_tx_antenna() << std::endl;

    // setup streamer
    uhd::stream_args_t stream_args("sc16", "sc16"); // complex floats
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    signal(SIGINT, sigint_handler); // stop streaming on STRG+C

    // setup the metadata flags
    uhd::tx_metadata_t tx_md;
    uhd::time_spec_t t = uhd::time_spec_t(0, 100);

    tx_md.start_of_burst = false;
    tx_md.end_of_burst = false;
    tx_md.has_time_spec = false;
    tx_md.time_spec = t;

    s->tx.stream = tx_stream;
    s->tx.md = tx_md;

    // Allocate TX buffer to hold each block of samples to transmit.
    s->tx.buffer = (short *)malloc(SAMPLES_PER_BUFFER * sizeof(short) * 2); // for 16-bit I and Q samples
    s->tx.buffer_ptr = (const void **)&s->tx.buffer;

    if (s->tx.buffer == NULL)
    {
        fprintf(stderr, "ERROR: Failed to allocate TX buffer.\n");
        // goto out;
    }

    size_t i = 0;
    for (i = 0; i < (samps_per_buff * 2); i += 2)
    {
        s->tx.buffer[i] = 0;
        s->tx.buffer[i + 1] = 0;
    }

    // Ctrl+C will exit loop
    // signal(SIGINT, &sigint_handler);
    fprintf(stderr, "Press Ctrl+C to stop streaming...\n");
    fprintf(stderr, "Opening and initializing USRP...\n");
}