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

# Detect Sollya's Presence
FIND_LIBRARY(
	SOLLYA_LIB
	NAMES sollya
	PATH /usr/lib/ /usr/local/lib)
IF (SOLLYA_LIB)
	MESSAGE(STATUS "Found Sollya: ${SOLLYA_LIB}")
ELSE (SOLLYA_LIB)
	MESSAGE(STATUS "Sollya not found, HOTBM won't be available")
ENDIF (SOLLYA_LIB)

IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)

# Build HOTBM
IF (SOLLYA_LIB)
	ADD_DEFINITIONS(-DHAVE_HOTBM)
	ADD_LIBRARY(hotbm STATIC
		src/HOTBM.cpp
		src/HOTBM/Function.cc
		src/HOTBM/HOTBMInstance.cc
		src/HOTBM/Param.cc
		src/HOTBM/Exhaustive.cc
		src/HOTBM/Util.cc
		src/HOTBM/PWPolynomial.cc
		src/HOTBM/Polynomial.cc
		src/HOTBM/MPPolynomial.cc
		src/HOTBM/Power.cc
		src/HOTBM/PowerAdHoc.cc
		src/HOTBM/TermROM.cc
		src/HOTBM/Term.cc
		src/HOTBM/Minimax.cc
		src/HOTBM/PowerROM.cc
		src/HOTBM/TermPowMult.cc
	)
	TARGET_LINK_LIBRARIES(hotbm sollya pari xml2 mpfi)
	SET(HOTBM_LIB hotbm)
ENDIF (SOLLYA_LIB)

ADD_EXECUTABLE(flopoco
 src/main
 src/utils
 src/Target src/Targets/VirtexIV
 src/Operator
 src/Table 
 src/Wrapper
 src/TestBench
 src/BigTestBench
 src/FPMultiplier
 src/IntMultiplier
 src/Karatsuba
 src/IntAdder
 src/LZOC
 src/Shifters
 src/LongAcc
 src/FPAdder
 src/DotProduct
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/FPConstMult 
 src/TestCase
 src/FloFP
 src/fpexp/math_lib.cpp
 src/fpexp/signal.cpp
 src/fpexp/fragment.cpp
 src/fpexp/stdfragment.cpp
 src/fpexp/logfragment.cpp
 src/fpexp/gen_table.cpp
 src/fpexp/explore.cpp
 src/FPExp
 src/fplog/FirstInvTable.cpp
 src/fplog/FirstLogTable.cpp
 src/fplog/OtherLogTable.cpp
 src/fplog/SecondInvTable.cpp
 src/fplog/Signal.cpp
 src/FPLog
 src/Mux
)
TARGET_LINK_LIBRARIES(flopoco  mpfr gmp gmpxx ${HOTBM_LIB})

ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp  mpfr gmp gmpxx)

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp  mpfr gmp gmpxx)
