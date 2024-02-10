#
# Copyright 2014-2015 Ettus Research LLC
# Copyright 2018 Ettus Research, a National Instruments Company
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

cmake_minimum_required(VERSION 3.5.1)
project(galileo-sdr CXX)

### Configure Compiler ########################################################

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_STANDARD 11)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CURSES_NEED_NCURSES TRUE)

# if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD" AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
#     set(CMAKE_EXE_LINKER_FLAGS "-lthr ${CMAKE_EXE_LINKER_FLAGS}")
#     set(CMAKE_CXX_FLAGS "-stdlib=libc++ ${CMAKE_CXX_FLAGS}")
# endif()


### Set up build environment ##################################################
# Choose a static or shared-library build (shared is default, and static will
# probably need some special care!)
# Set this to ON in order to link a static build of UHD:
option(UHD_USE_STATIC_LIBS OFF)
find_package(PkgConfig)

# To add UHD as a dependency to this project, add a line such as this:
find_package(UHD)
find_package(Boost COMPONENTS thread system program_options REQUIRED)
find_package(Curses REQUIRED)

# The version in  ^^^^^  here is a minimum version.
# To specify an exact version:
#find_package(UHD 3.15.0 EXACT REQUIRED)

# This example also requires Boost.
# Set components here, then include UHDBoost to do the actual finding
set(UHD_BOOST_REQUIRED_COMPONENTS
    program_options
    system
    thread
)

link_directories(${MATH_LIB_DIR} ${CURSES_INCLUDE_DIRS})
    # Name the library both with a full path and as "-lm" to
    # activate the link type switching code for both cases.
    # If the second one links shared then the link will fail.
set(MATH_LIBRARIES ${MATH_LIBRARY} -lm)

# need these include and link directories for the build
include_directories(
    ${UHD_INCLUDE_DIRS} ${CURSES_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${ZeroMQ_LIBRARY} ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src ${CMAKE_BINARY_DIR}/debug_scripts
)

### Make the executable #######################################################
file(GLOB SRC
     "src/*.cpp"
)

file(GLOB INC
     "include/*.h"
)

#add_executable(galileo-sdr galileo-sdr.cpp include/galileo-sdr.h include/FileIO.cpp include/pch.cpp include/TimeUtils.cpp include/StringUtils.cpp include/Rinex3Nav.cpp)
#add_executable(usrp_galileo usrp_galileo.cpp include/galileo-sdr.cpp include/nav_msg_data.cpp include/usrp_galileo.h include/galileo-sdr.h include/constants.h)
#add_executable(generate_frame debug_scripts/generate_frame.cpp)
#add_executable(rinex_reader debug_scripts/rinex_reader.cpp include/galileo-sdr.h)
add_executable(usrp_galileo src/main.cpp  ${SRC} ${INC})

set(LIBM_LIBRARIES)
message(STATUS "******************************************************************************")
message(STATUS "* NOTE: When building your own app, you probably need all kinds of different  ")
message(STATUS "* compiler flags. This is just an example, so it's unlikely these settings    ")
message(STATUS "* exactly match what you require. Make sure to double-check compiler and     ")
message(STATUS "* linker flags to make sure your specific requirements are included.          ")
message(STATUS "******************************************************************************")

# # Shared library case: All we need to do is link against the library, and
# # anything else we need (in this case, some Boost libraries):
# if(NOT UHD_USE_STATIC_LIBS)
#     message(STATUS "Linking against shared UHD library.")
#     target_link_libraries(usrpgps ${UHD_LIBRARIES} ${Boost_LIBRARIES} m Threads::Threads)

# # Shared library case: All we need to do is link against the library, and
# # anything else we need (in this case, some Boost libraries):
# else(NOT UHD_USE_STATIC_LIBS)
#     message(STATUS "Linking against static UHD library.")
#     target_link_libraries(usrpgps
#         # We could use ${UHD_LIBRARIES}, but linking requires some extra flags,
#         # so we use this convenience variable provided to us
#         ${UHD_STATIC_LIB_LINK_FLAG}
#         # Also, when linking statically, we need to pull in all the deps for
#         # UHD as well, because the dependencies don't get resolved automatically
#         ${UHD_STATIC_LIB_DEPS}
#     )
#endif(NOT UHD_USE_STATIC_LIBS)

pkg_check_modules(deps REQUIRED IMPORTED_TARGET glib-2.0)
#target_link_libraries (galileo-sdr ${CMAKE_THREAD_LIBS_INIT})
#target_link_libraries(galileo-sdr PkgConfig::deps)
target_link_libraries (usrp_galileo ${UHD_LIBRARIES} ${Boost_LIBRARIES} ${CURSES_LIBRARIES} ${ZeroMQ_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} PkgConfig::deps)

### Once it's built... ########################################################
# Here, you would have commands to install your program.
# We will skip these in this example.
