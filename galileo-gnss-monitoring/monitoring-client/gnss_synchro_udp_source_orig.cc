#include "gnss_synchro_udp_source.h"
#include "gnss_synchro.pb.h"
#include <sstream>
#include <ncurses.h>
#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <chrono>
using namespace std;
//using boost::asio::ip::udp;

Gnss_Synchro_Udp_Source::Gnss_Synchro_Udp_Source(const unsigned short port) :
    socket{io_service},
    endpoint{boost::asio::ip::udp::v4(), port}
{
    socket.open(endpoint.protocol(), error);  // Open socket.
    socket.bind(endpoint, error);             // Bind the socket to the given local endpoint.
    std::string hostname{"192.168.0.165"};

    sock = ::socket(AF_INET, SOCK_DGRAM, 0);
    destination.sin_family = AF_INET;
    destination.sin_port = htons(7531);
    destination.sin_addr.s_addr = inet_addr(hostname.c_str());

    // Hacky way of removing duplicates - TODO implement a nice buffer either here or at the receiver
    old_tow = 0;
    index = 0;
    sent_tow = false;
    MAX_CHAN = 8;

}

bool Gnss_Synchro_Udp_Source::read_gnss_synchro(gnss_sdr::Observables& stocks)
{
                        std::cerr << "Usage: monitoring-client <port>\n";

    char buff[1500];  // Buffer for storing the received data.

    // This call will block until one or more bytes of data has been received.
    int bytes = socket.receive(boost::asio::buffer(buff));
        
    std::string data(&buff[0], bytes);
    // Deserialize a stock of Gnss_Synchro objects from the binary string.
    return stocks.ParseFromString(data);
}

void Gnss_Synchro_Udp_Source::populate_channels(gnss_sdr::Observables& stocks)
{
    for (std::size_t i = 0; i < stocks.observable_size(); i++)
        {
            gnss_sdr::GnssSynchro ch = stocks.observable(i);
            if (ch.fs() != 0)  // Channel is valid.
                {
                    channels[ch.channel_id()] = ch;
                }
        }
}

bool Gnss_Synchro_Udp_Source::print_table()
{
    
    if (read_gnss_synchro(stocks))
        {
            populate_channels(stocks);

            clear();  // Clear the screen.

            // Print table header.
            attron(A_REVERSE);
            printw("%3s%6s%14s%17s%21s%25s%32s\n", "CH", "PRN", "CN0 [dB-Hz]", "Doppler [Hz]", "prompt_i", "rx_time","TOW_at_current_symbol_ms");
            attroff(A_REVERSE);

            // mycode _ test

            ofstream fout;
            ifstream fin;

            double new_tow;
            double diff;
            int bits[MAX_CHAN];

            for (int i = 0; i < MAX_CHAN; i++)
            {
                bits[i] = 0;
            }

            int bit = 0;
            // Print table contents.
            for (auto const& ch : channels)
                {
                    int channel_id = ch.first;      // Key
                    gnss_sdr::GnssSynchro data = ch.second;  // Value

                    printw("%3d%6d%14f%17f%21f%25f%32lu\n", channel_id, data.prn(), data.cn0_db_hz(), data.carrier_doppler_hz(), data.prompt_i(), data.rx_time(), data.tow_at_current_symbol_ms());

                    if (data.prompt_i() > 0)
                    { 
                        bit = 1;
                    }
                    else
                    {
                        bit = -1;
                    }
                    printw("%3d \n", bit);

                    bits[channel_id] = bit;
                    new_tow = data.tow_at_current_symbol_ms();
                    diff = data.rx_time();
                }
            
            if (old_tow != new_tow)
            {  
                sendto(sock, bits, MAX_CHAN * sizeof(int), 0, reinterpret_cast<sockaddr*>(&destination), sizeof(destination));

                if (!sent_tow)
                {   
                    new_tow = new_tow / 1000;
                    sent_tow = true;
                    sendto(sock, &new_tow, sizeof(double), 0, reinterpret_cast<sockaddr*>(&destination), sizeof(destination));
                }

                old_tow = new_tow;

                
            }
            
            refresh();  // Update the screen.
        }
    else
        {
            return false;
        }

    return true;
}

