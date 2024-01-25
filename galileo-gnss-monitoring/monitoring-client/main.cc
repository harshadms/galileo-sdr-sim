#include "gnss_synchro_udp_source.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <chrono>
#include <thread>
#include <ncurses.h>

int main(int argc, char* argv[])
{
    try
        {
            // Check command line arguments.
            if (argc != 2)
                {
                    // Print help.
                    std::cerr << "Usage: monitoring-client <port>\n";
                    return 1;
                }

            unsigned short port = boost::lexical_cast<unsigned short>(argv[1]);
            Gnss_Synchro_Udp_Source udp_source(port);

            initscr();  // Initialize ncurses.
            printw("Listening on port %d UDP...\n", port);
            refresh();  // Update the screen.
            
            while (true)
                {
                    udp_source.print_table();

                }
        }
    catch (std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }

    return 0;
}

