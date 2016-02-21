PROJECT(FloPoCo C CXX)

MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is " ${CMAKE_C_COMPILER})

# CMake's way of handling various versions
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
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
#SET(CMAKE_CXX_FLAGS_DEBUG "-Wall")
# Uncomment for profiling
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Wall")
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -Wall")
# Do you want an executable of 30MB ?
SET(CMAKE_BUILD_TYPE DEBUG)
#SET(CMAKE_BUILD_TYPE release)


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

FIND_LIBRARY(PYTHON_LIB
	NAMES libpython2.7.so
	HINTS "/usr/lib/x86_64-linux-gnu"
	DOC "Directory of the Python library"
)


# Detect Sollya's presence
FIND_LIBRARY( SOLLYA_LIB sollya )

IF (SOLLYA_LIB)
	MESSAGE(STATUS "Found Sollya: ${SOLLYA_LIB}")
ELSE (SOLLYA_LIB)
	MESSAGE(FATAL_ERROR "Sollya not found")
ENDIF (SOLLYA_LIB)

IF ( PYTHON_LIB )
	MESSAGE ( STATUS "Found Python : ${PYTHON_LIB}" )
ELSE ( PYTHON_LIB )
	MESSAGE ( FATAL_ERROR "Could not find PYTHON_LIB" )
ENDIF ( PYTHON_LIB )

IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)


FIND_PATH(WCPG_H
	  wcpg.h
	  DOC "Path of wcpg.h, the include file for GNU MPFR library"
		)

FIND_LIBRARY(WCPG_LIB
	NAMES wcpg
	DOC "Directory of the WCPG library"
	)

IF (WCPG_H AND WCPG_LIB)
    SET(WCPG_FOUND TRUE)
ENDIF (WCPG_H AND WCPG_LIB)

IF ( WCPG_FOUND )
	MESSAGE ( STATUS "Found WCPG : ${WCPG_H} ${WCPG_LIB}" )
  ADD_DEFINITIONS(-DHAVE_WCPG)
ELSE ( WCPG_FOUND )
	MESSAGE ( STATUS "Could not find WCPG, this will affect the FixIIR operator" )
ENDIF ( WCPG_FOUND )


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
	OUTPUT VHDLLexer.cpp VHDLLexer.hpp
	COMMAND ${FLEXPP_EXECUTABLE} 
	ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/VHDLLexer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/VHDLLexerRule.l
	)

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
	OUTPUT ExpressionParser.cpp
	COMMAND ${BISON_EXECUTABLE}
	ARGS -y -d --defines=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.h --output=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cpp --name-prefix=FlopocoExpression
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions/ExpressionParser.y
)

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
	OUTPUT ExpressionScanner.cpp
	COMMAND ${FLEX_EXECUTABLE} 
	ARGS --prefix=FlopocoExpression --header-file=${CMAKE_CURRENT_BINARY_DIR}/FlopocoScanner.h -o ${CMAKE_CURRENT_BINARY_DIR}/ExpressionScanner.cpp 
    ${CMAKE_SOURCE_DIR}/src/FPExpressions/ExpressionScanner.l
	)

SET_SOURCE_FILES_PROPERTIES(ExpressionScanner GENERATED)


ADD_LIBRARY(FloPoCoLib STATIC

#--------------------------------------------------------------------
#               FRAMEWORK
#--------------------------------------------------------------------

# src/FloPoCo.hpp
 src/Operator
 src/UserInterface
 src/Signal
 src/utils
 src/FlopocoStream
 src/TestState
 VHDLLexer
 
 # Resource Estimation related ---------------------------------------
 src/Tools/ResourceEstimationHelper
 
 # Floorplanning related ---------------------------------------------
 src/Tools/FloorplanningHelper

 # Targets -----------------------------------------------------------
 src/Targets/DSP
 src/Target 
 src/TargetModel
 
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
 src/TestBenches/TestCase
 src/TestBenches/FPNumber
 src/TestBenches/IEEENumber
 src/TestBenches/Wrapper
 src/TestBenches/TestBench 


# Bit heap  ------------------------------------------------
 src/BitHeap/BitHeap
 src/BitHeap/WeightedBit
 src/BitHeap/Plotter


#--------------------------------------------------------------------
#               OPERATORS          
#--------------------------------------------------------------------
#src/OperatorPipeline/OperatorPipeline
#src/FPExpressions/ExpressionParserData


# Conversions ----------------------------------------------------
src/Conversions/InputIEEE
src/Conversions/OutputIEEE
src/Conversions/FP2Fix
src/Conversions/Fix2FP

# Shifters/LZOC ----------------------------------------------------
 src/ShiftersEtc/LZOC 
 src/ShiftersEtc/LZOCShifterSticky
 src/ShiftersEtc/Shifters

# FixFilters ----------------------------------------------------------
 src/ShiftReg
 src/FixFilters/FixSOPC
 src/FixFilters/FixFIR
 src/FixFilters/FixHalfSine
 src/FixFilters/FixIIR


# pipelined Integer addition -------------------------------------------------
 src/IntAddSubCmp/IntAdder  
 src/IntAddSubCmp/IntAdderClassical 
 src/IntAddSubCmp/IntAdderAlternative 
 src/IntAddSubCmp/IntAdderShortLatency
 src/IntAddSubCmp/IntComparator
 src/IntAddSubCmp/IntDualSub
 src/IntAddSubCmp/BasicCompressor
 
# ---------- First generation of fast large adders ------------------------
 src/IntAddSubCmp/LongIntAdderAddAddMuxGen1 
 src/IntAddSubCmp/LongIntAdderCmpCmpAddGen1 
 src/IntAddSubCmp/LongIntAdderCmpAddIncGen1

# ---------- Second Generation of fast large adders -----------------------
 src/IntAddSubCmp/IntAdderSpecific 
 src/IntAddSubCmp/IntComparatorSpecific 
 src/IntAddSubCmp/CarryGenerationCircuit 
 src/IntAddSubCmp/LongIntAdderAddAddMuxGen2 
 src/IntAddSubCmp/LongIntAdderCmpCmpAddGen2 
 src/IntAddSubCmp/LongIntAdderCmpAddIncGen2
 src/IntAddSubCmp/LongIntAdderMuxNetwork

# Constant multiplications and divisions
 src/ConstMult/ShiftAddOp 
 src/ConstMult/ShiftAddDag 
 src/ConstMult/IntConstMult
 src/ConstMult/IntConstMCM  
 src/ConstMult/KCMTable
 src/ConstMult/IntIntKCM 
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/ConstMult/FixRealKCM
 src/ConstMult/FPRealKCM
 src/ConstMult/IntConstDiv 
 src/ConstMult/FPConstDiv 

 
# Multiplication related --------------------------------------------
 src/IntMult/IntMultiplier
 src/IntMult/FixMultAdd
 src/IntMult/MultiplierBlock
 src/IntMult/IntSquarer

# Complex numbers ---------------------------------------------------
src/Complex/FixComplexKCM
# TODO for Matei
#src/Complex/FixedComplexAdder

# Floating-point adders and subtracters -----------------------------
 src/FPAddSub/FPAdd
 src/FPAddSub/FPAddDualPath
 src/FPAddSub/FPAddSinglePath
 src/FPAddSub/FPAdd3Input
 src/FPAddSub/FPAddSub

# Floating-point multipliers and squarers
 src/FPMultSquare/FPMult
#  src/FPMultSquare/FPMultKaratsuba # See TODO
 src/FPMultSquare/FPSquare

# Floating-point dividers and square roots
src/FPDivSqrt/SelFunctionTable
src/FPDivSqrt/FPDiv
src/FPDivSqrt/FPSqrt
#src/FPDivSqrt/Tools/NbBitsMin
# src/FPDivSqrt/FPSqrtPoly # Resurrect some day?


# Composite operators
src/FPComposite/FPLargeAcc 
src/FPComposite/LargeAccToFP 
src/FPComposite/FPDotProduct 
src/FPComposite/FP2DNorm


# Generic fixed-point function evaluators ---------------------------
 src/Table 
 src/DualTable 

 src/FixConstant

 src/FixFunctions/FixFunction
 src/FixFunctions/FixFunctionByTable
 src/FixFunctions/BasicPolyApprox
 src/FixFunctions/PiecewisePolyApprox
 src/FixFunctions/FixHornerEvaluator
 src/FixFunctions/FixFunctionBySimplePoly
 src/FixFunctions/FixFunctionByPiecewisePoly

 src/FixFunctions/GenericTable
 src/FixFunctions/BipartiteTable
 src/FixFunctions/Multipartite
 src/FixFunctions/FixFunctionByMultipartiteTable
# src/FixFunctions/QuineMcCluskey

# Trigonometric functions ------------------------------------------
src/Trigs/ConstDiv3ForSinPoly
src/Trigs/FixSinPoly
src/Trigs/FixSinCos
src/Trigs/CordicSinCos
# src/Trigs/CordicSinOrCos
src/Trigs/FixAtan2
src/Trigs/FixAtan2ByCORDIC
src/Trigs/FixAtan2ByRecipMultAtan
src/Trigs/FixAtan2ByBivariateApprox
src/Trigs/Atan2Table

src/Trigs/Fix2DNorm


src/Tools/Plane
src/Tools/Point

# Exp, Log and Power
src/ExpLog/FPExp
src/ExpLog/IterativeLog
src/ExpLog/FPPow

# If you want to add your operator, feel free
src/UserDefinedOperator
)

TARGET_LINK_LIBRARIES(
  FloPoCoLib 
  mpfr gmp gmpxx xml2 mpfi
  )

IF (SOLLYA_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${SOLLYA_LIB}
	)
ENDIF (SOLLYA_LIB)

IF (PYTHON_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${PYTHON_LIB}
  )
ENDIF (PYTHON_LIB)

IF (WCPG_FOUND)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${WCPG_LIB}
  )
ENDIF (WCPG_FOUND)


ADD_EXECUTABLE(flopoco
 src/main 
)
 
TARGET_LINK_LIBRARIES(
  flopoco  
	FloPoCoLib 
  )

# unplugged for now to speed up compilation, to revive some day
#ADD_EXECUTABLE(flopoco_soaktest
# src/main_soaktest
#)
#TARGET_LINK_LIBRARIES(
#  flopoco_soaktest  
#	${FLOPOCO_LIB}
#)



ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp  mpfr gmp gmpxx)

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp  mpfr gmp gmpxx)



add_subdirectory(src/random)
target_link_libraries(FloPoCoLib flopoco_random_lib)
