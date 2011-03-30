# Makefile for flopoco/src
SUBDIRS = 
if BUILD_SOLLYA_FILES
SUBDIRS += FPExpressions
endif

ACLOCAL_AMFLAGS = -I ../m4

BUILT_SOURCES = VHDLLexer.cpp
CLEANFILES = $(BUILT_SOURCES) lex.yy.c *.output

noinst_LTLIBRARIES = libflopoco.la

noinst_HEADERS = apps/Collision.hpp apps/FPFMAcc.hpp apps/FPJacobi.hpp\
        ConstMult/IntIntKCM.hpp ConstMult/CRFPConstMult.hpp ConstMult/FixRealKCM.hpp ConstMult/FPRealKCM.hpp ConstMult/FPConstMult.hpp ConstMult/IntConstMult.hpp ConstMult/KCMTable.hpp\
        ConstMult/ShiftAddDag.hpp ConstMult/ShiftAddOp.hpp \
        fplogarithm/FirstInvTable.hpp fplogarithm/FirstLogTable.hpp fplogarithm/OtherLogTable.hpp \
        LNS/AtanPow.hpp LNS/Cotran.hpp LNS/CotranHybrid.hpp LNS/CotranTables.hpp LNS/LNSAdd.hpp LNS/LNSAddSub.hpp LNS/LNSDiv.hpp LNS/LNSMul.hpp LNS/LNSSqrt.hpp \
        squareroot/PolynomialTableCorrectRounded.hpp squareroot/PolynomialTable.hpp \
        Targets/DSP.hpp Targets/Spartan3.hpp Targets/StratixII.hpp Targets/StratixIV.hpp Targets/Virtex4.hpp Targets/Virtex5.hpp Targets/Virtex6.hpp \
        DotProduct.hpp  DualTable.hpp DualTable2.hpp Fix2FP.hpp FP2Fix.hpp FXP_Shift.hpp FlopocoStream.hpp FPAdderDualPath.hpp FPAdderSinglePath.hpp \
        FPDiv.hpp FPExp.hpp FPLog.hpp FPMultiplier.hpp FPNumber.hpp \
	FPSqrt.hpp FPSqrtPoly.hpp FPSquarer.hpp FixedPointFunctions/FunctionEvaluator.hpp FixedPointFunctions/HOTBM.hpp IEEENumber.hpp InputIEEE.hpp IntAdder.hpp IntAdderClassical.hpp IntAdderAlternative.hpp IntAdderShortLatency.hpp IntCompressorTree.hpp IntDualSub.hpp IntKaratsuba.hpp IntMultiplier.hpp \
	IntNAdder.hpp IntSquarer.hpp IntMultipliers/IntTilingMult.hpp IntTruncMultiplier.hpp LongAcc2FP.hpp LongAcc.hpp LongIntAdder.hpp LongIntAdderCmpCmpAdd.hpp LongIntAdderCmpAddAdd.hpp LZOC.hpp LZOCShifterSticky.hpp Operator.hpp OutputIEEE.hpp FixedPointFunctions/PolynomialEvaluator.hpp \
	Shifters.hpp Signal.hpp IntMultipliers/SignedIntMultiplier.hpp IntMultipliers/UnsignedIntMultiplier.hpp IntMultipliers/LogicIntMultiplier.hpp FixedPointFunctions/PolyTableGenerator.hpp FixedPointFunctions/PiecewiseFunction.hpp  FixedPointFunctions/FunctionTable.hpp Table.hpp Target.hpp TestBench.hpp TestCase.hpp UserDefinedOperator.hpp utils.hpp Wrapper.hpp 

libflopoco_la_SOURCES = VHDLLexer.cpp \
			utils.cpp \
			FlopocoStream.cpp \
			Targets/DSP.cpp \
			Target.cpp \
			Targets/Virtex4.cpp \
			Targets/Spartan3.cpp \
			Targets/Virtex6.cpp \
			Targets/StratixII.cpp \
			Targets/Virtex5.cpp \
			Targets/StratixIV.cpp \
			Operator.cpp \
			Signal.cpp \
			Table.cpp \
			DualTable.cpp \
			DualTable2.cpp \
			Wrapper.cpp \
			TestBench.cpp \
			UserDefinedOperator.cpp\
                        IntAdder.cpp IntAdderClassical.cpp IntAdderAlternative.cpp IntAdderShortLatency.cpp \
			IntNAdder.cpp \
			IntCompressorTree.cpp \
                        LongIntAdder.cpp LongIntAdderCmpCmpAdd.cpp LongIntAdderCmpAddAdd.cpp \
                        IntComparator.cpp \
			IntDualSub.cpp \
			IntMultiplier.cpp \
			IntMultipliers/SignedIntMultiplier.cpp \
			IntMultipliers/UnsignedIntMultiplier.cpp \
			IntMultipliers/IntTilingMult.cpp \
			IntMultipliers/LogicIntMultiplier.cpp \
			IntTruncMultiplier.cpp \
			IntKaratsuba.cpp \
			LZOC.cpp \
			LZOCShifterSticky.cpp \
			Shifters.cpp \
			FPAdderDualPath.cpp \
			Fix2FP.cpp \
			FP2Fix.cpp \
			FXP_Shift.cpp \
			FPMultiplier.cpp \
			FPMultiplierKaratsuba.cpp \
			FPMultiplierTiling.cpp \
			FPSquarer.cpp \
			FPDiv.cpp \
			FPSqrt.cpp \
			FPSqrtPoly.cpp \
			squareroot/PolynomialTable.cpp \
			squareroot/PolynomialTableCorrectRounded.cpp \
			LongAcc.cpp \
			LongAcc2FP.cpp \
			DotProduct.cpp \
			ConstMult/IntConstMult.cpp \
			ConstMult/ShiftAddOp.cpp \
			ConstMult/ShiftAddDag.cpp \
			ConstMult/IntIntKCM.cpp \
			ConstMult/KCMTable.cpp \
			ConstMult/FPConstMult.cpp \
			TestCase.cpp \
			FPNumber.cpp \
			IEEENumber.cpp \
			FPExp.cpp \
			fplogarithm/FirstInvTable.cpp \
			fplogarithm/FirstLogTable.cpp \
			fplogarithm/OtherLogTable.cpp \
			FPPow.cpp \
			FPLog.cpp \
			InputIEEE.cpp \
			OutputIEEE.cpp \
			apps/Collision.cpp \
                        apps/FPFMAcc.cpp \
			IntSquarer.cpp

if no_WIN32
libflopoco_la_SOURCES += FPAdderSinglePath.cpp \
			FPAdder3Input.cpp \
                        ConstMult/CRFPConstMult.cpp ConstMult/FixRealKCM.cpp ConstMult/FPRealKCM.cpp \
			FPPipeline.cpp \
			FixedPointFunctions/PolyTableGenerator.cpp \
			FixedPointFunctions/PiecewiseFunction.cpp \
			FixedPointFunctions/FunctionTable.cpp \
			FixedPointFunctions/PolynomialEvaluator.cpp \
			FixedPointFunctions/FunctionEvaluator.cpp \
			apps/FPJacobi.cpp \
			UtilSollya.cc

endif

if BUILD_SOLLYA_FILES
libflopoco_la_SOURCES += VHDLLexer.cpp \
		FixedPointFunctions/HOTBM.cpp \
		FixedPointFunctions/Function.cpp \
		FixedPointFunctions/HOTBM/HOTBMInstance.cc\
		FixedPointFunctions/HOTBM/Param.cc\
		FixedPointFunctions/HOTBM/Exhaustive.cc\
		FixedPointFunctions/HOTBM/Util.cc\
		FixedPointFunctions/HOTBM/PWPolynomial.cc\
		FixedPointFunctions/HOTBM/Polynomial.cc\
		FixedPointFunctions/HOTBM/MPPolynomial.cc\
		FixedPointFunctions/HOTBM/Power.cc\
		FixedPointFunctions/HOTBM/PowerAdHoc.cc\
		FixedPointFunctions/HOTBM/TermROM.cc\
		FixedPointFunctions/HOTBM/Term.cc\
		FixedPointFunctions/HOTBM/Minimax.cc\
		FixedPointFunctions/HOTBM/PowerROM.cc\
		FixedPointFunctions/HOTBM/TermPowMult.cc

libflopoco_la_LIBADD = $(top_builddir)/src/FPExpressions/libFPExpressions.la
libflopoco_la_CPPFLAGS = -I$(top_builddir)/src/FPExpressions/ -I$(top_srcdir)/src/FPExpressions

endif

VHDLLexer.cpp:$(top_srcdir)/src/VHDLLexerRule.l
	${LEX} -o $(top_builddir)/src/VHDLLexer.cpp $(top_srcdir)/src/VHDLLexerRule.l

VHDLLexer.hpp:$(top_builddir)/src/VHDLLexer.cpp
