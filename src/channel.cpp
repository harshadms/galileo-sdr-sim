#include "../include/galileo-sdr.h"

void init_channel(channel_t *chan, int *allocatedSat)
{
    // Clear all channels
    for (int i = 0; i < MAX_CHAN; i++)
    {
        // Set length of code array at runtime based on samples
        chan[i].set_code_phase = true;
        chan[i].prn = 0;
        chan[i].ca_E1B = (short *)malloc(2 * CA_SEQ_LEN_E1 * sizeof(short));
        chan[i].ca_E1C = (short *)malloc(2 * CA_SEQ_LEN_E1 * sizeof(short));
        chan[i].page = (int *)malloc(PAGE_SIZE * sizeof(int));
    }

    // Clear satellite allocation flag
    for (int sv = 0; sv < MAX_SAT; sv++)
        allocatedSat[sv] = -1;
}

int allocateChannel(channel_t *chan,
                    // map<int, vector<Rinex3Nav::DataGAL>> *navGAL,
                    vector<ephem_t> *eph_vector,
                    ionoutc_t ionoutc, galtime_t grx, double *xyz,
                    double elvMask, map<int, int> *sm,
                    vector<int> current_eph,
                    int *allocatedSat)
{
    int nsat = 0;
    int i, sv;
    double azel[2];

    range_t rho;
    double ref[3] = {0.0};
    double r_ref, r_xyz;
    double phase_ini;

    ephem_t eph;

    for (sv = 0; sv < MAX_SAT; sv++)
    {
        if (!eph_vector[sv].size())
            continue;

        if (!eph_vector[sv][0].vflg)
            continue;
        
        datetime_t t;
        gal2date(&grx, &t);
        double obs_time = gps_time(&t);
        
        current_eph[sv] = epoch_matcher(grx, eph_vector[sv], current_eph[sv]);
        // cout << "here : " << current_eph[sv] << " : " << sv<<  endl;
        if (current_eph[sv] < 0)
            continue;

        eph = eph_vector[sv][current_eph[sv]];
        //cout << "GRX: " << grx.sec << " : " << current_eph[sv] << endl;

        if (checkSatVisibility(eph, grx, xyz, 10, azel, sv + 1) == 1)
        {
            nsat++; // Number of visible satellites

            if (allocatedSat[sv] == -1) // Visible but not allocated
            {
                // Allocated new satellite
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn == 0)
                    {
                        // Initialize channel
                        chan[i].prn = sv + 1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];
                        chan[i].g0 = grx;

                        // Insert latest channel assignment to the map
                        // sm->insert({chan[i].prn, i});

                        // C/A code generation
                        codegen_E1B(chan[i].ca_E1B, chan[i].prn);
                        codegen_E1C(chan[i].ca_E1C, chan[i].prn);

                        // Generate navigation message
                        //generateNavMsg(grx, &chan[i], advance_fptr);
                        generateINavMsg(grx, &chan[i], &eph, &ionoutc);

                        // Initialize pseudorange
                        computeRange(&rho, eph, &ionoutc, grx, xyz, chan[i].prn);
                        chan[i].rho0 = rho;

                        // Initialize carrier phase
                        r_xyz = rho.range;

                        computeRange(&rho, eph, &ionoutc, grx, ref, chan[i].prn);
                        r_ref = rho.range;

                        phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_L1;
                        chan[i].carr_phase = phase_ini - floor(phase_ini);

                        fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.5f\n", chan[i].prn, chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.range, grx.sec);
                        //print_eph2(&eph, sv+1);
                        break;
                    }
                }

                // Set satellite allocation channel
                if (i < MAX_CHAN)
                    allocatedSat[sv] = i;
            }
        }
        else if (allocatedSat[sv] >= 0) // Not visible but allocated
        {
            // Clear channel
            chan[allocatedSat[sv]].prn = 0;

            // Clear satellite allocation flag
            allocatedSat[sv] = -1;
        }
    }
    // advance_fptr = true;
    return (nsat);
}
