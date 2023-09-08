# GALILEO-SDR-SIM

GALILEO-SDR-SIM generates Galileo baseband signal data streams, which can be converted to RF using software-defined radio (SDR) platforms.

- Galileo E1B/C signal generation
- Real time navigation message relay
- Real time location update
- USRP TX support 
- File sink

## Requirements

1. g++
2. Cmake
3. Curses
4. UHD (usrp api and dev library)
5. Boost
6. pynput (optional, for real-time location update via keyboard)

```
sudo apt-get install -y libuhd-dev uhd-host gnss-sdr g++ libncurses-dev cmake pkg-config libboost-dev libglib2.0-dev build-essential
pip3 install pynput
```

For offline evaluation, you will also need GNSS-SDR - found here -> https://gnss-sdr.org/

GNSS-SDR installation with package manager on most Linux distros
```
sudo apt install gnss-sdr
```
For more installation options please check their website.
## Installation
```
git clone https://github.com/harshadms/galileo-sdr-sim
cd galileo-sdr-sim
mkdir build && cd build
cmake ../
make
```

## Execution

Usage, help and options
```
Usage: ./usrp_galileo [options]
Options:
  -e <Ephemeris>   RINEX navigation file for Galileo ephemerides (required)
  -n <Navmsg Dir>  Dir containing navigation messages
  -o <File sink>   File to store IQ samples (default: ./galileosim.ishort)
  -u <user_motion> User motion file (dynamic mode)
  -l <location>    Lat,Lon,Hgt (static mode) e.g. 35.274,137.014,100
  -t <date,time>   Scenario start time YYYY/MM/DD,hh:mm:ss
  -T <date,time>   Overwrite TOC and TOE to scenario start time
  -d <duration>    Duration [sec] (default: 300)
  -a <rf_gain>     Absolute RF gain in [0 ... 60] (default: 30)
  -U <use usrp>    Disable USRP (-U 1) (default: true)
  -b <bitstream>   Disable Bit stream (-b 1) (default: true)
  -v <verbose>     Enable verbose output (default: false)
```

Executing the following command will generate a Galileo signal file for the location 42,71,100 starting from 2022/02/20,08:00:01. The generated samples will be stored in **galileo.ishort** as interleaved shorts (I<sub>1</sub>Q<sub>1</sub>, I<sub>2</sub>Q<sub>2</sub>, ... , I<sub>n</sub>Q<sub>n</sub>) @ 2.6e6 samples/sec

```
$ ./usrp_galileo -l 42,71,100 -t 2022/02/20,08:00:01 -n ../tv/20_FEB_2022_GST_08_00_01 -o /tmp/galileosim.bin -U 1 -b 1;
```
This project is still in development. For now, scenario start time is limited to the following because of navigation message availability:

- 11_DEC_2020_GST_08_00_01
- 12_DEC_2020_GST_08_00_01
- 13_DEC_2020_GST_09_00_01
- 15_DEC_2020_GST_08_00_01
- 15_DEC_2020_GST_10_00_01
- 16_DEC_2020_GST_08_00_01
- 17_JAN_2021_GST_08_00_01
- 12_JAN_2021_GST_10_00_01
- 13_JAN_2021_GST_10_00_01
- 25_JAN_2021_GST_08_00_01
- 23_FEB_2021_GST_09_00_31
- 28_FEB_2021_GST_08_00_01
- 20_FEB_2022_GST_08_00_01

Please specify the correct start time with `-t` option and correct navigation message files by specifying `-n`. All navigation messages are stored in `./tv/<date_time>`

## Evaluation

The generated E1B/C signal has been tested with GNSS-SDR and a battery of receivers from u-Blox and Septentrio. Here are some evaluation results.

<object data="./pvt_and_tracking.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="./pvt_and_tracking.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="pvt_and_tracking.pdf">Download PDF</a>.</p>
    </embed>
</object>

eval methodology - [add gnss-sdr instructions and USRP instructions]

## Future work

Here you will find all improvement plans and pending tasks. Please feel free to add suggestions and contribute :)

- [ ] I/NAV message generation (generate_frame.cpp has some PoC)
- [ ] RINEX reader integration
- [ ] Time & Frequency sync with GNSS-SDR monitor
- [ ] Enable multi-SDR support (libhackrf, LimeSuite etc)
- [ ] Fixed point arithmatic to avoid rounding errorss

More ambitious plans:

- [ ] Integrate GPS L1 C/A

#### Acknowledgements

- <a href="https://github.com/osqzss/gps-sdr-sim">GPS-SDR-SIM</a> has been an inspiration to start this project
