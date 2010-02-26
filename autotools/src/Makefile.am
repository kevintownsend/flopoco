# Makefile for flopoco/src

ACLOCAL_AMFLAGS = -I ../m4

AM_LFLAGS = -+ -o lex.yy.c --header-file=VHDLLexer.hpp
BUILT_SOURCES = VHDLLexer.hpp VHDLLexer.cpp 
CLEANFILES = $(BUILT_SOURCES) *.output
LEX_OUTPUT_ROOT=VHDLLexer.hpp VHDLLexer.hpp lex.yy

noinst_LTLIBRARIES = libflopoco.la

libflopoco_la_SOURCES = VHDLLexer.lpp \
			utils.cpp \
			FlopocoStream.cpp \
			Targets/DSP.cpp \
			Target.cpp \
			Targets/Virtex4.cpp \
			Targets/Spartan3.cpp \
			Targets/StratixII.cpp \
			Targets/Virtex5.cpp \
			Targets/StratixIV.cpp \
			Operator.cpp \
			Signal.cpp \
			Table.cpp \
			DualTable.cpp \
			TableGenerator.cpp \
			Wrapper.cpp \
			TestBench.cpp \
			IntAdder.cpp \
			IntNAdder.cpp \
			IntCompressorTree.cpp \
			LongIntAdder.cpp \
			IntDualSub.cpp \
			IntMultiplier.cpp \
			SignedIntMultiplier.cpp \
			IntKaratsuba.cpp \
			LZOC.cpp \
			LZOCShifterSticky.cpp \
			Shifters.cpp \
			FPAdder.cpp \
			Fix2FP.cpp \
			FPMultiplier.cpp \
			FPSquarer.cpp \
			FPDiv.cpp \
			FPSqrt.cpp \
			squareroot/PolynomialTable.cpp \
			squareroot/PolynomialTableCorrectRounded.cpp \
			LongAcc.cpp \
			LongAcc2FP.cpp \
			DotProduct.cpp \
			ConstMult/IntConstMult.cpp \
			ConstMult/ShiftAddOp.cpp \
			ConstMult/ShiftAddDag.cpp \
			ConstMult/IntIntKCM.cpp \
			ConstMult/FirstKCMTable.cpp \
			ConstMult/FPConstMult.cpp \
			TestCase.cpp \
			FPNumber.cpp \
			IEEENumber.cpp \
			fpexponential/math_lib.cpp \
			fpexponential/signal.cpp \
			fpexponential/fragment.cpp \
			fpexponential/stdfragment.cpp \
			fpexponential/logfragment.cpp \
			fpexponential/gen_table.cpp \
			fpexponential/explore.cpp \
			FPExp.cpp \
			fplogarithm/FirstInvTable.cpp \
			fplogarithm/FirstLogTable.cpp \
			fplogarithm/OtherLogTable.cpp \
			FPLog.cpp \
			InputIEEE.cpp \
			OutputIEEE.cpp \
			apps/Collision.cpp \
			IntSquarer.cpp

if no_WIN32
libflopoco_la_SOURCES += ConstMult/LastKCMTable.cpp ConstMult/CRFPConstMult.cpp \
			PolynomialEvaluator.cpp \
			FunctionEvaluator.cpp \
			UtilSollya.cc
endif


