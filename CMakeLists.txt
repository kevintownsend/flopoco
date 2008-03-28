PROJECT(FloPoCo C CXX)
MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is "${CMAKE_C_COMPILER})

INCLUDE(CPack)


FIND_PATH(MPFR_H
	  mpfr.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of mpfr.h, the include file for GNU MPFR library"
)

FIND_LIBRARY(MPFR_LIB 
             NAMES mpfr 
             PATH /usr/lib/ /usr/local/lib/  
             DOC "Directory of the MPFR library")


ADD_EXECUTABLE(flopoco
 main
 utils
 Target Targets/ClassicalFPGA
 Operator
 Table 
 Wrapper
 Multiplier
 IntMultiplier
 IntAdder
 LZOC
 Shifters
 LongAcc
 ConstMult/IntConstMult ConstMult/ShiftAddOp ConstMult/ShiftAddDag 
 ConstMult/FPConstMult 
)
TARGET_LINK_LIBRARIES(flopoco  mpfr gmp gmpxx)

ADD_EXECUTABLE(fp2bin Tools/fp2bin  utils)
TARGET_LINK_LIBRARIES(fp2bin  mpfr gmp gmpxx)

