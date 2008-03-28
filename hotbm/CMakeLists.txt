PROJECT(HOTBM C CXX)

INCLUDE(CPack)

SET(CMAKE_BUILD_TYPE debug)

SET(IncludePaths /usr/local/include /usr/include)
SET(LibPaths     /usr/local/lib     /usr/lib)

FIND_PATH   (MPFR_H   NAMES mpfr.h PATHS ${IncludePaths} DOC "Path of mpfr.h, the include file for GNU MPFR library")
FIND_LIBRARY(MPFR_LIB NAMES mpfr   PATHS ${LibPaths}     DOC "Directory of the MPFR library")

IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)

ADD_EXECUTABLE(
    hotbm
    src/Exhaustive.cc
    src/Function.cc
    src/HOTBMInstance.cc
    src/Log.cc
    src/Main.cc
    src/Minimax.cc
    src/MPPolynomial.cc
    src/Param.cc
    src/Polynomial.cc
    src/PowerAdHoc.cc
    src/Power.cc
    src/PowerROM.cc
    src/PWPolynomial.cc
    src/Sin.cc
    src/Term.cc
    src/TermPowMult.cc
    src/TermROM.cc
    src/Util.cc
)
TARGET_LINK_LIBRARIES(hotbm mpfr)
