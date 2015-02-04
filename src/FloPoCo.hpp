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

/* operator pipeline work* ------------------------------------ */
#include "OperatorPipeline/OperatorPipeline.hpp"

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

#include "TestBenches/TestBench.hpp"

/* shifters + lzoc ------------------------------------------- */
#include "ShiftersEtc/Shifters.hpp"
#include "ShiftersEtc/LZOC.hpp"
#include "ShiftersEtc/LZOCShifterSticky.hpp"


/* FixFilters ------------------------------------------------ */
#include "FixFilters/FixSOPC.hpp"
#include "FixFilters/FixFIR.hpp"
#include "FixFilters/FixHalfSine.hpp"
#include "FixFilters/FixIIR.hpp"

#include "ShiftReg.hpp"

/* regular pipelined integer adder/ adder+subtracter --------- */
#include "IntAddSubCmp/IntAdder.hpp" // includes several other .hpp
#include "IntAddSubCmp/IntDualSub.hpp"
#include "IntAddSubCmp/IntComparator.hpp"

/* fast large adders ----------------------------------------- */
#include "IntAddSubCmp/LongIntAdderAddAddMuxGen1.hpp"
#include "IntAddSubCmp/LongIntAdderCmpCmpAddGen1.hpp"
#include "IntAddSubCmp/LongIntAdderCmpAddIncGen1.hpp"
#include "IntAddSubCmp/IntAdderSpecific.hpp"
#include "IntAddSubCmp/LongIntAdderAddAddMuxGen2.hpp"
#include "IntAddSubCmp/LongIntAdderCmpCmpAddGen2.hpp"
#include "IntAddSubCmp/LongIntAdderCmpAddIncGen2.hpp"
#include "IntAddSubCmp/IntComparatorSpecific.hpp"
#include "IntAddSubCmp/LongIntAdderMuxNetwork.hpp"

/* Integer and fixed-point multipliers ------------------------ */
#include "IntMult/IntMultiplier.hpp"
#include "IntMult/FixMultAdd.hpp"
// #include "IntMult/IntKaratsuba.hpp"
#include "IntMult/IntSquarer.hpp"
// #include "IntMult/GenericBinaryPolynomial.hpp"
// #include "IntMult/IntPower.hpp"


/* Floating-point adder variants ----------------------------- */
#include "FPAddSub/FPAddDualPath.hpp"
#include "FPAddSub/FPAddSinglePath.hpp"
#include "FPAddSub/FPAdd3Input.hpp"
#include "FPAddSub/FPAddSub.hpp"

/* Floating-point multiplier variants-------------------------- */ 
#include "FPMultSquare/FPMult.hpp"
//#include "FPMultKaratsuba.hpp" // Resurrect some day?
#include "FPMultSquare/FPSquare.hpp"

#include "FPDivSqrt/FPDiv.hpp"
#include "FPDivSqrt/FPSqrt.hpp"
//#include "FPDivSqrt/FPSqrtPoly.hpp" // Resurrect some day?


/* Constant multipliers and dividers ------------------------ */
#include "ConstMult/IntConstMult.hpp"
#include "ConstMult/IntIntKCM.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "ConstMult/FPConstMult.hpp"
#include "ConstMult/CRFPConstMult.hpp"
#include "ConstMult/FPRealKCM.hpp"

#include "ConstMult/IntConstDiv.hpp"
#include "ConstMult/FPConstDiv.hpp"

/* FP composite operators */ 
#include "FPComposite/FPLargeAcc.hpp"
#include "FPComposite/LargeAccToFP.hpp"
#include "FPComposite/FPDotProduct.hpp"


/* Fixed-point function generators ---------------------*/

#include "FixFunctions/FixFunction.hpp"
#include "FixFunctions/BasicPolyApprox.hpp"
#include "FixFunctions/PiecewisePolyApprox.hpp"
#include "FixFunctions/FixFunctionByTable.hpp"
#include "FixFunctions/FixFunctionBySimplePoly.hpp"
#include "FixFunctions/FixFunctionByPiecewisePoly.hpp"

#include "FixFunctions/BipartiteTable.hpp"
#include "FixFunctions/GenericTable.hpp"

/*  Various elementary functions in fixed or floating point*/
#include "Trigs/FixSinCos.hpp"
#include "Trigs/CordicSinCos.hpp"
#include "Trigs/CordicAtan2.hpp"
#include "Trigs/FixAtan2.hpp"
// #include "Trigs/FixSinOrCos.hpp"  Replug when poly eval fixed
#include "ExpLog/IterativeLog.hpp"
#include "ExpLog/FPExp.hpp"
#include "ExpLog/FPPow.hpp"


#include "Conversions/Fix2FP.hpp"
#include "Conversions/FP2Fix.hpp"
#include "Conversions/InputIEEE.hpp"
#include "Conversions/OutputIEEE.hpp"



/* misc ------------------------------------------------------ */
#include "TestBenches/Wrapper.hpp"
#include "UserDefinedOperator.hpp"




#if 0
// Old stuff removed from older versions, some of which to bring back to life
/* Complex arithmetic */
#include "Complex/FixedComplexAdder.hpp"
#include "Complex/FixedComplexMultiplier.hpp"



/* multioperand adders --------------------------------------- */
#include "IntMultiAdder.hpp"
#include "IntAddSubCmp/IntNAdder.hpp"
#include "IntAddSubCmp/IntCompressorTree.hpp"
#include "IntAddSubCmp/PopCount.hpp"
#include "IntAddSubCmp/BasicCompressor.hpp"
#include "IntAddSubCmp/NewCompressorTree.hpp"
#include "FP2DNorm.hpp"
#include "FPSqrtPoly.hpp"
/* applications ---------------------------------------------- */

/* Coil Inductance application */
//#include "Apps/CoilInductance/CoordinatesTableX.hpp"
//#include "Apps/CoilInductance/CoordinatesTableZ.hpp"
//#include "Apps/CoilInductance/CoordinatesTableY.hpp"
//#include "Apps/CoilInductance/CoilInductance.hpp"

/* fast evaluation of the possible intrusion of a point within a 
spheric enclosure -------------------------------------------- */ 
#include "Apps/Collision.hpp"

/* a floating-point fused multiply-accumulate operator for the 
use withing matrix-multiplication scenarios ------------------ */ 
#include "Apps/FPFMAcc.hpp"

/* a 1D Jacobi computation kernel ---------------------------- */ 
#include "Apps/FPJacobi.hpp"

/* logarithmic number system  -------------------------------- */ 
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
#endif // 0


#endif //FLOPOCO_HPP
