cmake_minimum_required(VERSION 3.13)
project(vaja_VI_1)

set(CMAKE_CXX_STANDARD 14)

add_executable(run main.cpp SpectralHe.cpp SpectralHe.h Helium.cpp Helium.h NumerovHe.cpp NumerovHe.h)
target_link_libraries(run gsl gslcblas)
