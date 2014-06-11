#ifndef FLOPOCO_HPP
#define FLOPOCO_HPP

// TODO: I guess we should at some point copy here only the public part of each class, 
// to provide a single self-contained include file.

// support the autotools-generated config.h
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Operator.hpp"
#include "FlopocoStream.hpp"


/* resource estimation ---------------------------------------- */
#include "Tools/ResourceEstimationHelper.hpp"
/* resource estimation ---------------------------------------- */


/* floorplanning ---------------------------------------------- */
#include "Tools/ResourceEstimationHelper.hpp"
/* floorplanning ---------------------------------------------- */

/* targets ---------------------------------------------------- */
#include "Target.hpp"
#include "Targets/DSP.hpp"

#include "Targets/Spartan3.hpp"
#include "Targets/Virtex4.hpp"
#include "Targets/Virtex5.hpp"
#include "Targets/Virtex6.hpp"

#include "Targets/StratixII.hpp"
#include "Targets/StratixIII.hpp"
#include "Targets/StratixIV.hpp"
#include "Targets/StratixV.hpp"

#include "Targets/CycloneII.hpp"
#include "Targets/CycloneIII.hpp"
#include "Targets/CycloneIV.hpp"
#include "Targets/CycloneV.hpp"



#include "TestBench.hpp"

/* shifters + lzoc ------------------------------------------- */
#include "ShiftersEtc/Shifters.hpp"
#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/LZOCShifterSticky.hpp"

/* regular pipelined integer adder/ adder+subtracter --------- */
#include "IntAdders/IntAdder.hpp" // includes several other .hpp
#include "IntAdders/IntDualSub.hpp"
#include "IntAdders/IntComparator.hpp"

/* fast large adders ----------------------------------------- */
#include "IntAdders/LongIntAdderAddAddMuxGen1.hpp"
#include "IntAdders/LongIntAdderCmpCmpAddGen1.hpp"
#include "IntAdders/LongIntAdderCmpAddIncGen1.hpp"
#include "IntAdders/IntAdderSpecific.hpp"
#include "IntAdders/LongIntAdderAddAddMuxGen2.hpp"
#include "IntAdders/LongIntAdderCmpCmpAddGen2.hpp"
#include "IntAdders/LongIntAdderCmpAddIncGen2.hpp"
#include "IntAdders/IntComparatorSpecific.hpp"
#include "IntAdders/LongIntAdderMuxNetwork.hpp"

/* Integer and fixed-point multipliers ------------------------ */
#include "IntMultipliers/IntMultiplier.hpp"
#include "IntMultipliers/FixMultAdd.hpp"
// #include "IntMultipliers/IntKaratsuba.hpp"
#include "IntMultipliers/IntSquarer.hpp"
// #include "IntMultipliers/GenericBinaryPolynomial.hpp"
// #include "IntMultipliers/IntPower.hpp"


/* Floating-point adder variants ----------------------------- */
#include "FPAddSub/FPAdderDualPath.hpp"
#include "FPAddSub/FPAdderSinglePath.hpp"
#include "FPAddSub/FPAdder3Input.hpp"
#include "FPAddSub/FPAddSub.hpp"

/* Floating-point multiplier variants-------------------------- */ 
#include "FPMultSquare/FPMultiplier.hpp"
//#include "FPMultiplierKaratsuba.hpp"
#include "FPMultSquare/FPSquarer.hpp"



/* Constant multipliers and dividers ------------------------ */
#include "ConstMult/IntConstMult.hpp"
#include "ConstMult/IntIntKCM.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "ConstMult/FPConstMult.hpp"
#include "ConstMult/CRFPConstMult.hpp"
#include "ConstMult/FPRealKCM.hpp"

#include "ConstMult/IntConstDiv.hpp"
#include "ConstMult/FPConstDiv.hpp"

/* Fixed-point function generators ---------------------*/

#include "FixFunctions/FixFunction.hpp"
#include "FixFunctions/BasicPolyApprox.hpp"
#include "FixFunctions/PiecewisePolyApprox.hpp"
#include "FixFunctions/FixFunctionByTable.hpp"
#include "FixFunctions/FixFunctionBySimplePoly.hpp"

/*  Various trigonometric functions ----------------------------*/
#include "Trigs/FixSinCos.hpp"
// #include "CORDIC/FixedPointSinOrCos.hpp"
// #include "CORDIC/CordicSinCos.hpp"



#if 0



/* multioperand adders --------------------------------------- */
#include "IntMultiAdder.hpp"
#include "IntAdders/IntNAdder.hpp"
#include "IntAdders/IntCompressorTree.hpp"
#include "IntAdders/PopCount.hpp"
#include "IntAdders/BasicCompressor.hpp"
#include "IntAdders/NewCompressorTree.hpp"

/* multiplication-related ------------------------------------ */
#include "IntMultiplier.hpp"
#include "FixMultAdd.hpp"
#include "IntMultipliers/IntKaratsuba.hpp"
#include "IntSquarer.hpp"
#include "IntMultipliers/GenericBinaryPolynomial.hpp"
#include "IntMultipliers/IntPower.hpp"

#include "IntMultipliers/FixSinPoly.hpp"
#include "IntMultipliers/FixXPow3Div6.hpp"
#include "IntMultipliers/MultiplierBlock.hpp"





/* fixed-point function evaluation---------------------------- */
#ifndef _WIN32

#ifdef HAVE_SOLLYA
#endif
#endif
/* fixed-point ----------------------------------------------- */
#ifdef HAVE_SOLLYA

#include "FixedPointFIR.hpp"

#include "FixedPointDCT.hpp"

#endif

/* floating-point -------------------------------------------- */ 
#include "FPMultiplier.hpp"
#include "FPMultiplierKaratsuba.hpp"
#include "FPSquarer.hpp"

#ifndef _WIN32
#ifdef HAVE_SOLLYA
#include "ConstMult/CRFPConstMult.hpp"
#endif
#endif

#include "FPDiv.hpp"
#include "FPExp.hpp" 
#include "FPLog.hpp"
#include "FPPow.hpp"

#include "FPSqrt.hpp"
// #include "FP2DNorm.hpp" // The world is not ready yet 
#include "FPSqrtPoly.hpp"

#include "LongAcc.hpp"
#include "LongAcc2FP.hpp"

#include "DotProduct.hpp"
#include "FPSumOfSquares.hpp"

#include "FPPipeline.hpp"

#include "Fix2FP.hpp"
#include "FP2Fix.hpp"
#include "InputIEEE.hpp"
#include "OutputIEEE.hpp"


/* Complex arithmetic */
#include "Complex/FixedComplexAdder.hpp"
#include "Complex/FixedComplexMultiplier.hpp"


/* test-bench related ---------------------------------------- */
#include "TestBench.hpp"


/* applications ---------------------------------------------- */

/* Coil Inductance application */
//#include "apps/CoilInductance/CoordinatesTableX.hpp"
//#include "apps/CoilInductance/CoordinatesTableZ.hpp"
//#include "apps/CoilInductance/CoordinatesTableY.hpp"
//#include "apps/CoilInductance/CoilInductance.hpp"

/* fast evaluation of the possible intrusion of a point within a 
spheric enclosure -------------------------------------------- */ 
#include "apps/Collision.hpp"

/* a floating-point fused multiply-accumulate operator for the 
use withing matrix-multiplication scenarios ------------------ */ 
#include "apps/FPFMAcc.hpp"

/* a 1D Jacobi computation kernel ---------------------------- */ 
#include "apps/FPJacobi.hpp"

/* logarithmic number system  -------------------------------- */ 
#ifndef _WIN32
#ifdef HAVE_LNS
#include "LNS/LNSAddSub.hpp"
#include "LNS/LNSAdd.hpp"
#include "LNS/CotranTables.hpp"
#include "LNS/Cotran.hpp"
#include "LNS/CotranHybrid.hpp"
#include "LNS/LNSMul.hpp"
#include "LNS/LNSDiv.hpp"
#include "LNS/LNSSqrt.hpp"
#include "LNS/AtanPow.hpp"
#include "LNS/LogSinCos.hpp"
#endif
#endif
#endif /////////////////////////////////////////////////////0

/* misc ------------------------------------------------------ */
#include "Wrapper.hpp"
#include "UserDefinedOperator.hpp"
#include "Plotter.hpp"


#endif //FLOPOCO_HPP
