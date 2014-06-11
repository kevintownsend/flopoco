PROJECT(FloPoCo C CXX)

MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is " ${CMAKE_C_COMPILER})

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
# 
#Compilation flags
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall")
# Uncomment for profiling
# SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Wall")
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -Wall")
# Do you want an executable of 30MB ?
#SET(CMAKE_BUILD_TYPE DEBUG)
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

FIND_PATH(MPFI_H
	  mpfi.h
	  HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/include"
	  DOC "Path of mpfi.h"
)

FIND_LIBRARY(MPFI_LIB
	NAMES mpfi mpfi.lib
	HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/lib"
	DOC "Directory of the MPFI library")



# Detect Sollya's presence
FIND_LIBRARY( SOLLYA_LIB sollya )

IF (SOLLYA_LIB)
	ADD_DEFINITIONS(-DHAVE_SOLLYA)
	MESSAGE(STATUS "Found Sollya: ${SOLLYA_LIB}")
ELSE (SOLLYA_LIB)
	MESSAGE(STATUS "Sollya not found, many very nice operators won't be available")
ENDIF (SOLLYA_LIB)


IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)

ADD_DEFINITIONS(-DHAVE_LNS)

# Create target for the lexer
 ADD_CUSTOM_TARGET(VHDLLexer echo "Creating VHDLLexer.c VHDLLexer.h")

# Create custom command for flex++/lex (note the outputs)
FIND_PROGRAM(FLEXPP_EXECUTABLE
	NAMES flex
	PATHS /usr/bin
	DOC "Flex++ executable")

IF(FLEXPP_EXECUTABLE)
	MESSAGE(STATUS "Found lex: ${FLEXPP_EXECUTABLE}")
ELSE(FLEXPP_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find LEX.\nOn Debian-based systems install flex++")
ENDIF(FLEXPP_EXECUTABLE)

# necessary to include generated files
INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions)

ADD_CUSTOM_COMMAND(
	SOURCE src/VHDLLexerRule.l
	COMMAND ${FLEXPP_EXECUTABLE} 
	ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/VHDLLexer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/VHDLLexerRule.l
	TARGET VHDLLexer
	OUTPUTS VHDLLexer.cpp VHDLLexer.hpp)

SET_SOURCE_FILES_PROPERTIES(VHDLLexer GENERATED) 


# Create custom command for bison/yacc (note the outputs)
FIND_PROGRAM(BISON_EXECUTABLE
	NAMES bison
	PATHS /usr/bin
	DOC "Bison executable")

IF(BISON_EXECUTABLE)
	MESSAGE(STATUS "Found bison: ${BISON_EXECUTABLE}")
ELSE(BISON_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find Bison.\nOn Debian-based systems install bison")
ENDIF(BISON_EXECUTABLE)

ADD_CUSTOM_COMMAND(
	SOURCE src/FPExpressions/ExpressionParser.y
	COMMAND ${BISON_EXECUTABLE}
	ARGS -y -d --defines=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.h --output=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cpp --name-prefix=FlopocoExpression
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions/ExpressionParser.y
	TARGET ExpressionParser
	OUTPUTS ExpressionParser.cpp)

SET_SOURCE_FILES_PROPERTIES(ExpressionParser GENERATED) 

# Create custom command for flex/lex (note the outputs)
FIND_PROGRAM(FLEX_EXECUTABLE
	NAMES flex
	PATHS /usr/bin
	DOC "Flex executable")

IF(FLEX_EXECUTABLE)
	MESSAGE(STATUS "Found flex: ${FLEX_EXECUTABLE}")
ELSE(FLEX_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find LEX.\nOn Debian-based systems install flex++")
ENDIF(FLEX_EXECUTABLE)

ADD_CUSTOM_COMMAND(
	SOURCE src/FPExpressions/ExpressionScanner.l
	COMMAND ${FLEX_EXECUTABLE} 
	ARGS --prefix=FlopocoExpression --header-file=${CMAKE_CURRENT_BINARY_DIR}/FlopocoScanner.h -o ${CMAKE_CURRENT_BINARY_DIR}/ExpressionScanner.cpp 
    ${CMAKE_SOURCE_DIR}/src/FPExpressions/ExpressionScanner.l
	TARGET ExpressionScanner
	OUTPUTS ExpressionScanner.cpp)

SET_SOURCE_FILES_PROPERTIES(ExpressionScanner GENERATED)


ADD_LIBRARY(FloPoCoLib STATIC

#--------------------------------------------------------------------
#               FRAMEWORK
#--------------------------------------------------------------------

 src/FloPoCo.hpp
 src/Operator
 src/Signal
 src/utils
 src/FlopocoStream
 VHDLLexer
 
 # Resource Estimation related ---------------------------------------
 src/Tools/ResourceEstimationHelper
 
 # Floorplanning related ---------------------------------------------
 src/Tools/FloorplanningHelper

 # Targets -----------------------------------------------------------
 src/Targets/DSP
 src/Target 
 
 src/Targets/Spartan3 
 src/Targets/Virtex4 
 src/Targets/Virtex5 
 src/Targets/Virtex6 
 
 src/Targets/StratixII 
 src/Targets/StratixIII 
 src/Targets/StratixIV
 src/Targets/StratixV
 
 src/Targets/CycloneII 
 src/Targets/CycloneIII 
 src/Targets/CycloneIV 
 src/Targets/CycloneV 

# TestBench-related -------------------------------------------------
 src/TestCase
 src/FPNumber
 src/IEEENumber
 src/Wrapper
 src/TestBench 


# Bit heap  ------------------------------------------------
 src/BitHeap
 src/WeightedBit
 src/Plotter


#--------------------------------------------------------------------
#               OPERATORS          
#--------------------------------------------------------------------

# Shifters/LZOC ----------------------------------------------------
 src/ShiftersEtc/LZOC 
 src/ShiftersEtc/LZOCShifterSticky
 src/ShiftersEtc/Shifters

# pipelined Integer addition -------------------------------------------------
 src/IntAdders/IntAdder  
 src/IntAdders/IntAdderClassical 
 src/IntAdders/IntAdderAlternative 
 src/IntAdders/IntAdderShortLatency
 src/IntAdders/IntComparator
 src/IntAdders/IntDualSub
 src/IntAdders/BasicCompressor
 
# ---------- First generation of fast large adders ------------------------
 src/IntAdders/LongIntAdderAddAddMuxGen1 
 src/IntAdders/LongIntAdderCmpCmpAddGen1 
 src/IntAdders/LongIntAdderCmpAddIncGen1

# ---------- Second Generation of fast large adders -----------------------
 src/IntAdders/IntAdderSpecific 
 src/IntAdders/IntComparatorSpecific 
 src/IntAdders/CarryGenerationCircuit 
 src/IntAdders/LongIntAdderAddAddMuxGen2 
 src/IntAdders/LongIntAdderCmpCmpAddGen2 
 src/IntAdders/LongIntAdderCmpAddIncGen2
 src/IntAdders/LongIntAdderMuxNetwork

# Constant multiplications and divisions
 src/ConstMult/ShiftAddOp 
 src/ConstMult/ShiftAddDag 
 src/ConstMult/IntConstMult 
 src/ConstMult/KCMTable
 src/ConstMult/IntIntKCM 
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/ConstMult/FixRealKCM
 src/ConstMult/FPRealKCM
 src/ConstMult/IntConstDiv 
 src/ConstMult/FPConstDiv 

 
# Multiplication related --------------------------------------------
 src/IntMultipliers/IntMultiplier
 src/IntMultipliers/FixMultAdd
 src/IntMultipliers/FixMultAddBitheap
 src/IntMultipliers/MultiplierBlock
 src/IntMultipliers/IntSquarer


# Floating-point adders and subtracters -----------------------------
 src/FPAddSub/FPAdderDualPath
 src/FPAddSub/FPAdderSinglePath
 src/FPAddSub/FPAdder3Input
 src/FPAddSub/FPAddSub

# Floating-point multipliers and squarers
 src/FPMultSquare/FPMultiplier
#  src/FPMultSquare/FPMultiplierKaratsuba # See TODO
 src/FPMultSquare/FPSquarer



# Generic fixed-point function evaluators ---------------------------
 src/Table 
 src/DualTable 
 
 src/FixFunctions/FixFunction
 src/FixFunctions/FixFunctionByTable
 src/FixFunctions/BasicPolyApprox
 src/FixFunctions/PiecewisePolyApprox
 src/FixFunctions/FixHornerEvaluator
 src/FixFunctions/FixFunctionBySimplePoly
 src/FixConstant

# Trigonometric functions ------------------------------------------
src/Trigs/ConstDiv3ForSinPoly
src/Trigs/FixSinPoly
src/Trigs/FixSinCos

)
IF (SOLLYA_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  mpfr gmp gmpxx xml2 mpfi ${SOLLYA_LIB}
  )
ELSE (SOLLYA_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib 
  mpfr gmp gmpxx xml2 mpfi
  )
ENDIF (SOLLYA_LIB)
SET(FLOPOCO_LIB FloPoCoLib)

ADD_EXECUTABLE(flopoco
 src/main 
)

TARGET_LINK_LIBRARIES(
  flopoco  
	${FLOPOCO_LIB}
  )




ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp  mpfr gmp gmpxx)

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp  mpfr gmp gmpxx)



add_subdirectory(src/random)
target_link_libraries(FloPoCoLib flopoco_random_lib)
