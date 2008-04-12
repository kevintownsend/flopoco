PROJECT(FloPoCo C CXX)
MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is "${CMAKE_C_COMPILER})

INCLUDE(CPack)

SET(CMAKE_BUILD_TYPE debug)

FIND_PATH(MPFR_H
	  mpfr.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of mpfr.h, the include file for GNU MPFR library"
)

FIND_LIBRARY(MPFR_LIB 
             NAMES mpfr 
             PATH /usr/lib/ /usr/local/lib/  
             DOC "Directory of the MPFR library")

IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)


ADD_EXECUTABLE(flopoco
 src/main
 src/utils
 src/Target src/Targets/VirtexIV
 src/Operator
 src/Table 
 src/Wrapper
 src/TestBench
 src/FPMultiplier
 src/IntMultiplier
 src/IntAdder
 src/LZOC
 src/Shifters
 src/LongAcc
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/FPConstMult 
 src/TestCase
)
TARGET_LINK_LIBRARIES(flopoco  mpfr gmp gmpxx)

ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp  mpfr gmp gmpxx)

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp  mpfr gmp gmpxx)
