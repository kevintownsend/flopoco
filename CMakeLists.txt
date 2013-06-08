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
  ADD_DEFINITIONS(-DHAVE_LNS)
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
 src/FloPoCo.hpp
 src/Operator
 src/Signal
 src/utils
 src/FlopocoStream
 VHDLLexer
 
 # Resource Estimation related ---------------------------------------
 #--------------------------------------------------------------------
 src/Tools/ResourceEstimationHelper
 #--------------------------------------------------------------------
 
 # Floorplanning related ---------------------------------------------
 #--------------------------------------------------------------------
 src/Tools/FloorplanningHelper
 #--------------------------------------------------------------------

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


 src/Table 
 src/DualTable 
 src/DualTable2
 
 src/FPExpressions/ExpressionParserData
 ExpressionParser 
 ExpressionScanner
 src/FPPipeline

 src/Wrapper
 src/TestBench 
 src/UserDefinedOperator
 
# Addition related --------------------------------------------------
#--------------------------------------------------------------------
 src/IntAdder 
 
 src/IntAddition/IntAdderClassical 
 src/IntAddition/IntAdderAlternative 
 src/IntAddition/IntAdderShortLatency

# ---------- First generation of fast adders ------------------------
 src/IntAddition/LongIntAdderAddAddMuxGen1 
 src/IntAddition/LongIntAdderCmpCmpAddGen1 
 src/IntAddition/LongIntAdderCmpAddIncGen1

# ---------- Second Generation of fast adders -----------------------
 src/IntAddition/IntAdderSpecific 
 src/IntAddition/IntComparatorSpecific 
 src/IntAddition/CarryGenerationCircuit 
 src/IntAddition/LongIntAdderAddAddMuxGen2 
 src/IntAddition/LongIntAdderCmpCmpAddGen2 
 src/IntAddition/LongIntAdderCmpAddIncGen2

 src/IntAddition/LongIntAdderMuxNetwork

# ---------- Multioperand -------------------------------------------
 src/IntMultiAdder src/IntAddition/IntNAdder  src/IntAddition/IntCompressorTree
 src/IntAddition/PopCount src/IntAddition/NewCompressorTree
 src/IntAddition/BasicCompressor
 src/BitHeap
src/WeightedBit
 src/IntComparator
 src/IntDualSub

#---------Plotting-related ---------------------------------------
src/Plotter


 
# Multiplication related --------------------------------------------
#--------------------------------------------------------------------
 src/IntMultiplier
 src/IntMultipliers/IntKaratsuba
 src/IntSquarer
 src/FixMultAdd
 src/IntMultipliers/FormalBinaryProduct
 src/IntMultipliers/GenericBinaryPolynomial
 src/IntMultipliers/IntPower

src/IntMultipliers/MultiplierBlock

src/IntMultipliers/FixSinPoly
src/IntMultipliers/FixXPow3Div6
src/ConstMult/IntConstDiv3

# Shifters/LZOC related ---------------------------------------------
#--------------------------------------------------------------------
 src/LZOC 
 src/LZOCShifterSticky
 src/Shifters

# Fixed-point function evaluation -----------------------------------
#--------------------------------------------------------------------
 src/FixFunctions/Function
 src/FixFunctions/PiecewiseFunction

 src/FixFunctions/FunctionTable

# src/FixFunctions/Multipartite

 src/FixFunctions/PolyCoeffTable
 src/FixFunctions/PolynomialEvaluator
 src/FixFunctions/FunctionEvaluator
 src/FixFunctions/GenericEvaluator

 src/FixFunctions/HOTBM.cpp
 src/FixFunctions/HOTBM/HOTBMInstance.cc
 src/FixFunctions/HOTBM/Param.cc
 src/FixFunctions/HOTBM/Exhaustive.cc
 src/FixFunctions/HOTBM/Util.cc
 src/FixFunctions/HOTBM/PWPolynomial.cc
 src/FixFunctions/HOTBM/Polynomial.cc
 src/FixFunctions/HOTBM/MPPolynomial.cc
 src/FixFunctions/HOTBM/Power.cc
 src/FixFunctions/HOTBM/PowerAdHoc.cc
 src/FixFunctions/HOTBM/TermROM.cc
 src/FixFunctions/HOTBM/Term.cc
 src/FixFunctions/HOTBM/Minimax.cc
 src/FixFunctions/HOTBM/PowerROM.cc
 src/FixFunctions/HOTBM/TermPowMult.cc 

# Fixed-point operators ---------------------------------------------
#--------------------------------------------------------------------
 src/FixSinCos/FixSinCos
 src/FixSinCos/FixSinOrCos
 src/FixSinCos/CordicSinCos

 src/FixFIR
 src/FixDCT

# Floating-Point ----------------------------------------------------
#--------------------------------------------------------------------
 src/FPAdderDualPath
 src/FPAdderSinglePath
 src/FPAdder3Input
 src/FPAddSub

 src/FPMultiplier
 src/FPMultiplierKaratsuba
 src/FPSquarer

 src/FPDiv
 
 src/FPSqrt
 src/FPSqrtPoly
 src/squareroot/PolynomialTable 
 src/squareroot/PolynomialTableCorrectRounded

# src/FP2DNorm
 
 src/DotProduct 
 
 src/ConstMult/IntConstMult src/ConstMult/ShiftAddOp src/ConstMult/ShiftAddDag 
 src/ConstMult/IntIntKCM src/ConstMult/KCMTable
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/ConstMult/FixRealKCM
 src/ConstMult/FPRealKCM

 src/IntConstDiv 
 src/FPConstDiv 

 src/FPLog/FirstInvTable
 src/FPLog/FirstLogTable
 src/FPLog/OtherLogTable
 src/FPLog.cpp

 src/FPExp 
 src/FPPow
 
 src/LongAcc
 src/LongAcc2FP

 src/FPSumOfSquares

# ---------- Complex arithmetic --------------------------------------------- 
src/Complex/FixedComplexAdder
src/Complex/FixedComplexMultiplier

# ---------- Conversion --------------------------------------------- 

 src/FP2Fix
 src/Fix2FP
 src/InputIEEE
 src/OutputIEEE

#---------------------LNS arithmetic--------------------------------
 src/LNS/LNSAdd
 src/LNS/Cotran
 src/LNS/CotranHybrid
 src/LNS/CotranTables
 src/LNS/LNSAddSub
 src/LNS/LNSMul
 src/LNS/LNSDiv
 src/LNS/LNSSqrt
 src/LNS/AtanPow
 src/LNS/LogSinCos


# Applications ------------------------------------------------------
#--------------------------------------------------------------------

# src/apps/CoilInductance/CoilInductance
# src/apps/CoilInductance/CoordinatesTableX
# src/apps/CoilInductance/CoordinatesTableZ
# src/apps/CoilInductance/CoordinatesTableY

 src/apps/Collision
 
 src/apps/FPFMAcc
 src/apps/FPJacobi

# src/apps/TaMaDiCore 
# src/apps/TaMaDiPriorityEncoder 
# src/apps/TaMaDiDecoder 
# src/apps/TaMaDiFIFO 
# src/apps/TaMaDiShiftRegister 
# src/apps/TaMaDiModule 
# src/apps/TaMaDiModuleDummyWrapper
# src/apps/TaMaDiModuleWrapperInterface
# src/apps/TaMaDiDispatcherInterface
# src/apps/TaMaDiDeserializer
# src/apps/TaMaDiSystem

# TestBench-related -------------------------------------------------
#--------------------------------------------------------------------
 src/TestCase
 src/FPNumber
 src/IEEENumber

 src/UtilSollya
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


ADD_EXECUTABLE(fpadder_example
 src/main_minimal 
)
TARGET_LINK_LIBRARIES(
  fpadder_example  
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
