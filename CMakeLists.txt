PROJECT(FloPoCo C CXX)

MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is "${CMAKE_C_COMPILER})

# CMake's way of handling various versions
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


#####################For packaging (source file distrib so far)#######
# Use  make package_source to build a release
# Currently it doesn't work at all

SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

INCLUDE(CPack)


# Uncomment for profiling
#SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg")

# Do you want an executable of 9GB ?
#SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_BUILD_TYPE release)

FIND_PATH(MPFR_H
	  mpfr.h
	  HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/include"
	  DOC "Path of mpfr.h, the include file for GNU MPFR library"
)

FIND_LIBRARY(MPFR_LIB
	NAMES mpfr mpfr.lib
	HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/lib"
	DOC "Directory of the MPFR library")


# Detect Sollya's Presence
FIND_LIBRARY( SOLLYA_LIB sollya )

IF (SOLLYA_LIB)
	ADD_DEFINITIONS(-DHAVE_SOLLYA)
	MESSAGE(STATUS "Found Sollya: ${SOLLYA_LIB}")
ELSE (SOLLYA_LIB)
	MESSAGE(STATUS "Sollya not found, HOTBM and CRFPConstMult won't be available")
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
	TARGET_LINK_LIBRARIES(hotbm  xml2 mpfi ${SOLLYA_LIB})
	SET(HOTBM_LIB hotbm)
ENDIF (SOLLYA_LIB)

IF (HOTBM_LIB)
	ADD_DEFINITIONS(-DHAVE_LNS)
	ADD_LIBRARY(lns STATIC
		src/LNS/LNSAdd
		src/LNS/Cotran
		src/LNS/CotranHybrid
		src/LNS/CotranTables
		src/LNS/LNSAddSub
		src/LNS/LNSMul
		src/LNS/LNSDiv
		src/LNS/LNSSqrt
	)
	TARGET_LINK_LIBRARIES(lns ${HOTBM_LIB})
	SET(LNS_LIB lns)
ENDIF (HOTBM_LIB)


ADD_EXECUTABLE(flopoco
 src/main
 src/utils
 src/Target src/Targets/VirtexIV src/Targets/StratixII
 src/Operator
 src/Table 
 src/Wrapper
 src/TestBench 
 src/BigTestBench
 src/GenericTestBench
 src/IntAdder
 src/IntDualSub
 src/IntMultiplier
 src/Karatsuba
 src/LZOC src/LZOCShifterSticky
 src/Shifters
 src/FPAdder
 src/FPMultiplier
 src/FPDiv.cpp src/FPDiv.hpp
 src/FPDiv/SRT4Step
 src/IntMult2
 src/LongAcc
 src/LongAcc2FP
 src/DotProduct
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/TestCase
 src/TestIOMap
 src/FPNumber
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
 src/fplog/LogRangeRed.cpp
 src/FPLog
)

IF (SOLLYA_LIB)
  TARGET_LINK_LIBRARIES(
    flopoco  
    mpfr gmp gmpxx xml2 mpfi lns hotbm  ${SOLLYA_LIB}
    )
ELSE (SOLLYA_LIB)
  TARGET_LINK_LIBRARIES(
    flopoco  
    mpfr gmp gmpxx)
ENDIF (SOLLYA_LIB)


ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp  mpfr gmp gmpxx)

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp  mpfr gmp gmpxx)
