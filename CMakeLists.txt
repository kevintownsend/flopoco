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

IF(WIN32)
SET(CMAKE_BUILD_TYPE release)
ELSE(WIN32) 
#Compilation flags
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall")
# Uncomment for profiling
#SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Wall")
# Do you want an executable of 12MB ?
SET(CMAKE_BUILD_TYPE debug)
#SET(CMAKE_BUILD_TYPE release)
ENDIF(WIN32)


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

IF(WIN32)
ADD_EXECUTABLE(flopoco
 src/main
 src/utils
 src/Target src/Targets/VirtexIV src/Targets/StratixII
 src/Operator
 src/Signal
 src/Table 
 src/Wrapper
 src/TestBench 
# src/BigTestBench
 src/GenericTestBench
 src/IntAdder
 src/IntNAdder
 src/LongIntAdder
 src/IntDualSub
 src/IntMultiplier
 src/IntMultiplier2
 src/IntMult2
 src/Karatsuba
 src/LZOC src/LZOCShifterSticky
 src/Shifters
 src/FPAdder
 src/FPMultiplier
 src/FPSquarer
 src/FPDiv
 src/division/SRT4Step
 src/FPSqrt
 src/squareroot/PolynomialTable
 src/squareroot/SqrtSRT2Step
 src/LongAcc
 src/LongAcc2FP
 src/DotProduct
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/IntIntKCM src/ConstMult/FirstKCMTable
 src/ConstMult/FPConstMult 
# src/ConstMult/CRFPConstMult 
 src/TestCase
 src/FPNumber
 
 
 src/fpexponential/math_lib
 src/fpexponential/signal
 src/fpexponential/fragment
 src/fpexponential/stdfragment
 src/fpexponential/logfragment
 src/fpexponential/gen_table
 src/fpexponential/explore
 src/FPExp
 
 src/fplogarithm/FirstInvTable
 src/fplogarithm/FirstLogTable
 src/fplogarithm/OtherLogTable
 src/fplogarithm/SecondInvTable
 src/fplogarithm/LogRangeRed
 src/FPLog
 
 src/Collision
 src/IntSquarer
)

ELSE(WIN32)
ADD_EXECUTABLE(flopoco
 src/main
 src/utils
 src/Target src/Targets/VirtexIV src/Targets/StratixII
 src/Operator
 src/Signal
 src/Table 
 src/Wrapper
 src/TestBench 
 src/BigTestBench
 src/GenericTestBench
 src/IntAdder
 src/IntNAdder
 src/LongIntAdder
 src/IntDualSub
 src/IntMultiplier
 src/IntMultiplier2
 src/Karatsuba
 src/LZOC src/LZOCShifterSticky
 src/Shifters
 src/FPAdder
 src/FPMultiplier
 src/FPSquarer
 src/FPDiv
 src/division/SRT4Step
 src/FPSqrt
 src/squareroot/PolynomialTable
 src/squareroot/SqrtSRT2Step
 src/IntMult2
 src/LongAcc
 src/LongAcc2FP
 src/DotProduct
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/IntIntKCM src/ConstMult/FirstKCMTable src/ConstMult/LastKCMTable
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/TestCase
 src/FPNumber
 src/fpexponential/math_lib.cpp
 src/fpexponential/signal.cpp
 src/fpexponential/fragment.cpp
 src/fpexponential/stdfragment.cpp
 src/fpexponential/logfragment.cpp
 src/fpexponential/gen_table.cpp
 src/fpexponential/explore.cpp
 src/FPExp
 src/fplogarithm/FirstInvTable
 src/fplogarithm/FirstLogTable
 src/fplogarithm/OtherLogTable
 src/fplogarithm/SecondInvTable
 src/fplogarithm/LogRangeRed
 src/FPLog
 src/Collision
 src/IntSquarer
)
ENDIF(WIN32) 




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
