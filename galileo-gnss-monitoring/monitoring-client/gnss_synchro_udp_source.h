#ifndef GNSS_SYNCHRO_UDP_SOURCE_H
#define GNSS_SYNCHRO_UDP_SOURCE_H
#include <boost/asio.hpp>
#include <arpa/inet.h> // htons, inet_addr
#include <netinet/in.h> // sockaddr_in
#include <sys/types.h> // uint16_t
#include <sys/socket.h> // socket, sendto
#include <unistd.h> // close
#include "gnss_synchro.pb.h"  // This file is created automatically
                              // by the Protocol Buffers compiler

class Gnss_Synchro_Udp_Source
{
public:
    Gnss_Synchro_Udp_Source(const unsigned short port);
    bool read_gnss_synchro(gnss_sdr::Observables& stocks);
    void populate_channels(gnss_sdr::Observables& stocks);
    bool print_table();

private:
    boost::asio::io_service io_service;
    boost::asio::ip::udp::socket socket;
    boost::system::error_code error;
    boost::asio::ip::udp::endpoint endpoint;
    gnss_sdr::Observables stocks;
    std::map<int, gnss_sdr::GnssSynchro> channels;

    int sock;
    sockaddr_in destination;
    double old_tow;
    int preamble[8] = {1,-1,-1,-1,1,-1,1,1};
    int index;
    bool sent_tow;
    int MAX_CHAN;
    FILE *fp;

};

#endif  // GNSS_SYNCHRO_UDP_SOURCE_H

