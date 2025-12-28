cmake_minimum_required(VERSION 3.10)
project(GPSAcquisition)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_library(FFTW3_LIBRARY NAMES fftw3 fftw3-3 libfftw3-3)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)

if(NOT FFTW3_LIBRARY OR NOT FFTW3_INCLUDE_DIR)
    message(FATAL_ERROR "FFTW3 library not foun")
else()
    include_directories(${FFTW3_INCLUDE_DIR})
endif()

include_directories(include)


add_executable(gps_test 
    src/main.cpp 
    src/Acquisition.cpp 
)


target_link_libraries(gps_test ${FFTW3_LIBRARY} m)