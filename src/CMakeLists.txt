cmake_minimum_required(VERSION 3.22)
project(modeling_BLISS_attacks)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgmpxx -lgmp -lntl -lmpfr -lgmp -lcrypto -lboost_program_options")
add_executable(modeling_BLISS_attacks
        ansicolor.h
        Entropy.cpp
        Entropy.h
        Fft.cpp
        Fft.h
        hybrid_demo.cpp
        KeyGen.cpp
        KeyGen.h
        Parameters.h
        rlutil.h
        Sampler.cpp
        Sampler.h
        Setup.cpp
        Setup.h
        Sign.cpp
        Sign.h
        Verify.cpp
        Verify.h)
target_link_libraries(modeling_BLISS_attacks ntl gmpxx gmp mpfr crypto boost_program_options)
