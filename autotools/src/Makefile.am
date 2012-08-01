# Makefile for flopoco/src
SUBDIRS = 

ACLOCAL_AMFLAGS = -I ../m4

BUILT_SOURCES = VHDLLexer.cpp \
                ExpressionParser.cpp \
                ExpressionScanner.cpp
CLEANFILES = $(BUILT_SOURCES) lex.yy.c *.output

noinst_LTLIBRARIES = libflopoco.la

noinst_HEADERS = \
          FloPoCo.hpp \
          apps/Collision.hpp \
          apps/FPFMAcc.hpp \
          apps/FPJacobi.hpp \
          ConstMult/IntIntKCM.hpp \
          ConstMult/CRFPConstMult.hpp \
          ConstMult/FixRealKCM.hpp \
          ConstMult/FPRealKCM.hpp \
          ConstMult/FPConstMult.hpp \
          ConstMult/IntConstMult.hpp \
          ConstMult/KCMTable.hpp \
          ConstMult/ShiftAddDag.hpp \
          ConstMult/ShiftAddOp.hpp \
          fplogarithm/FirstInvTable.hpp \
          fplogarithm/FirstLogTable.hpp \
          fplogarithm/OtherLogTable.hpp \
          LNS/AtanPow.hpp \
          LNS/Cotran.hpp \
          LNS/CotranHybrid.hpp \
          LNS/CotranTables.hpp \
          LNS/LNSAdd.hpp \
          LNS/LNSAddSub.hpp \
          LNS/LNSDiv.hpp \
          LNS/LNSMul.hpp \
          LNS/LNSSqrt.hpp \
          squareroot/PolynomialTableCorrectRounded.hpp \
          squareroot/PolynomialTable.hpp \
          Targets/DSP.hpp \
          Targets/Spartan3.hpp \
          Targets/StratixII.hpp \
          Targets/StratixIII.hpp \
          Targets/StratixIV.hpp \
          Targets/Virtex4.hpp \
          Targets/Virtex5.hpp \
          Targets/Virtex6.hpp \
          DotProduct.hpp \
          DualTable.hpp \
          DualTable2.hpp \
          Fix2FP.hpp \
          FP2Fix.hpp \
          FlopocoStream.hpp \
          FPAdderDualPath.hpp \
          FPAdderSinglePath.hpp \
          FPDiv.hpp \
          FPExp.hpp \
          FPLog.hpp \
          FPMultiplier.hpp \
          FPNumber.hpp \
          FPSqrt.hpp \
          FPSqrtPoly.hpp \
          FPSquarer.hpp \
          FixedPointFunctions/FunctionEvaluator.hpp \
          FixedPointFunctions/HOTBM.hpp \
          IEEENumber.hpp \
          InputIEEE.hpp \
          IntAdder.hpp \
          IntAddition/IntAdderClassical.hpp \
          IntAddition/IntAdderAlternative.hpp \
          IntAddition/IntAdderShortLatency.hpp \
          IntAddition/IntCompressorTree.hpp \
          IntMultiAdder.hpp \
          IntDualSub.hpp \
          IntMultipliers/IntKaratsuba.hpp \
          IntMultiplier.hpp \
          IntAddition/IntNAdder.hpp \
          IntSquarer.hpp \
          IntMultipliers/IntTilingMult.hpp \
          IntMultipliers/IntTruncMultiplier.hpp \
          IntMultipliers/FormalBinaryProduct.hpp \
          IntMultipliers/GenericBinaryPolynomial.hpp \
          IntMultipliers/IntPower.hpp \
          IntMultipliers/MultiplierBlock.hpp \
          LongAcc2FP.hpp \
          LongAcc.hpp LZOC.hpp \
          LZOCShifterSticky.hpp \
          Operator.hpp \
          OutputIEEE.hpp \
          FixedPointFunctions/PolynomialEvaluator.hpp \
          Shifters.hpp \
          Signal.hpp \
          IntMultipliers/SignedIntMultiplier.hpp \
          IntMultipliers/UnsignedIntMultiplier.hpp \
          IntMultipliers/LogicIntMultiplier.hpp \
          FixedPointFunctions/PiecewiseFunction.hpp \
          FixedPointFunctions/FunctionTable.hpp \
          FixedPointFunctions/GenericEvaluator.hpp \
					Complex/FixedComplexAdder.hpp \
					Complex/FixedComplexMultiplier.hpp \
          Table.hpp \
          Target.hpp \
          TestBench.hpp \
          TestCase.hpp \
          UserDefinedOperator.hpp \
          utils.hpp \
          Wrapper.hpp \
          FPSumOfSquares.hpp \
          FixSinCos.hpp \
          FPAddSub.hpp

libflopoco_la_SOURCES = \
                        Operator.cpp \
                        Signal.cpp \
                        utils.cpp \
                        FlopocoStream.cpp \
                        VHDLLexer.cpp \
                        Targets/DSP.cpp \
                        Tools/ResourceEstimationHelper.cpp \
                        Tools/FloorplanningHelper.cpp \
                        Target.cpp \
                        Targets/Virtex4.cpp \
                        Targets/Spartan3.cpp \
                        Targets/Virtex6.cpp \
                        Targets/StratixII.cpp \
                        Targets/StratixIII.cpp \
                        Targets/Virtex5.cpp \
                        Targets/StratixIV.cpp \
                        Table.cpp \
                        DualTable.cpp \
                        DualTable2.cpp \
                        FPExpressions/ExpressionParserData.cpp \
                        ExpressionParser.cpp \
                        ExpressionScanner.cpp \
                        FPPipeline.cpp \
                        Wrapper.cpp \
                        TestBench.cpp \
                        UserDefinedOperator.cpp\
                        IntAdder.cpp\
                        IntAddition/IntAdderClassical.cpp\
                        IntAddition/IntAdderAlternative.cpp\
                        IntAddition/IntAdderShortLatency.cpp\
                        IntAddition/LongIntAdderAddAddMuxGen1.cpp \
                        IntAddition/LongIntAdderCmpCmpAddGen1.cpp \
                        IntAddition/LongIntAdderCmpAddIncGen1.cpp \
                        IntAddition/IntAdderSpecific.cpp \
                        IntAddition/IntComparatorSpecific.cpp \
                        IntAddition/CarryGenerationCircuit.cpp \
                        IntAddition/LongIntAdderAddAddMuxGen2.cpp \
                        IntAddition/LongIntAdderCmpCmpAddGen2.cpp \
                        IntAddition/LongIntAdderCmpAddIncGen2.cpp \
                        IntAddition/LongIntAdderMuxNetwork.cpp \
                        IntMultiAdder.cpp \
                        IntAddition/IntNAdder.cpp \
                        IntAddition/IntCompressorTree.cpp \
                        IntAddition/PopCount.cpp \
                        IntAddition/NewCompressorTree.cpp \
                        IntAddition/BasicCompressor.cpp \
                        BitHeap.cpp \
                        IntComparator.cpp \
                        IntDualSub.cpp \
                        IntMultiplier.cpp \
                        IntMultipliers/SignedIntMultiplier.cpp \
                        IntMultipliers/UnsignedIntMultiplier.cpp \
                        IntMultipliers/IntTilingMult.cpp \
                        IntMultipliers/LogicIntMultiplier.cpp \
                        IntMultipliers/IntTruncMultiplier.cpp \
                        IntMultipliers/IntKaratsuba.cpp \
                        IntSquarer.cpp \
                        IntMultipliers/FormalBinaryProduct.cpp \
                        IntMultipliers/GenericBinaryPolynomial.cpp \
                        IntMultipliers/IntPower.cpp \
                        IntMultipliers/MultiplierBlock.cpp \
                        LZOC.cpp \
                        LZOCShifterSticky.cpp \
                        Shifters.cpp \
                        FixedPointFunctions/Function.cpp \
                        FixedPointFunctions/PiecewiseFunction.cpp \
                        FixedPointFunctions/FunctionTable.cpp \
                        FixedPointFunctions/PolyCoeffTable.cpp \
                        FixedPointFunctions/PolynomialEvaluator.cpp \
                        FixedPointFunctions/FunctionEvaluator.cpp \
                        FixedPointFunctions/HOTBM.cpp \
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
                        FixedPointFunctions/HOTBM/TermPowMult.cc \
                        FixSinCos.cpp \
                        CORDIC/FixedPointSinOrCos.cpp \
                        CORDIC/CordicSinCos.cpp \
                        FPAdderDualPath.cpp \
                        FPAdderSinglePath.cpp \
                        FPAdder3Input.cpp \
                        FPAddSub.cpp \
                        FPMultiplier.cpp \
                        FPMultiplierKaratsuba.cpp \
                        FPSquarer.cpp \
                        FPDiv.cpp \
                        FPSqrt.cpp \
                        FPSqrtPoly.cpp \
                        squareroot/PolynomialTable.cpp \
                        squareroot/PolynomialTableCorrectRounded.cpp \
                        DotProduct.cpp \
                        ConstMult/IntConstMult.cpp \
                        ConstMult/ShiftAddOp.cpp \
                        ConstMult/ShiftAddDag.cpp \
                        ConstMult/IntIntKCM.cpp \
                        ConstMult/KCMTable.cpp \
                        ConstMult/FPConstMult.cpp \
                        ConstMult/CRFPConstMult.cpp \
                        ConstMult/FixRealKCM.cpp \
                        ConstMult/FPRealKCM.cpp \
                        IntConstDiv.cpp\
                        FPConstDiv.cpp\
                        fplogarithm/FirstInvTable.cpp \
                        fplogarithm/FirstLogTable.cpp \
                        fplogarithm/OtherLogTable.cpp \
                        FPLog.cpp \
                        FPExp.cpp \
                        FPPow.cpp \
                        LongAcc.cpp \
                        LongAcc2FP.cpp \
                        FPSumOfSquares.cpp \
                        Fix2FP.cpp \
                        InputIEEE.cpp \
                        OutputIEEE.cpp \
                        LNS/LNSAdd.cpp \
                        LNS/Cotran.cpp \
                        LNS/CotranHybrid.cpp \
                        LNS/CotranTables.cpp \
                        LNS/LNSAddSub.cpp \
                        LNS/LNSMul.cpp \
                        LNS/LNSDiv.cpp \
                        LNS/LNSSqrt.cpp \
                        LNS/AtanPow \
                        FixedPointFunctions/GenericEvaluator.cpp \
												Complex/FixedComplexAdder.cpp \
												Complex/FixedComplexMultiplier.cpp \
                        LNS/LogSinCos.cpp \
                        apps/Collision.cpp \
                        apps/FPFMAcc.cpp \
                        apps/FPJacobi.cpp \
                        TestCase.cpp \
                        FPNumber.cpp \
                        IEEENumber.cpp \
                        UtilSollya.cc \
                        FP2Fix.cpp \
                        random/random_main.cpp \
                        random/urng/LutSrRng.cpp

libflopoco_la_CPPFLAGS = -I$(top_srcdir)/src/FPExpressions

VHDLLexer.cpp:$(top_srcdir)/src/VHDLLexerRule.l
	${LEX} -o VHDLLexer.cpp $(top_srcdir)/src/VHDLLexerRule.l

VHDLLexer.hpp:VHDLLexer.cpp

ExpressionParser.cpp: $(top_srcdir)/src/FPExpressions/ExpressionParser.y ExpressionScanner.cpp
	${YACC} -y -d --defines=ExpressionParser.h --output=ExpressionParser.cpp --name-prefix=FlopocoExpression $(top_srcdir)/src/FPExpressions/ExpressionParser.y

ExpressionScanner.cpp:$(top_srcdir)/src/FPExpressions/ExpressionScanner.l
	${LEX} --prefix=FlopocoExpression --header-file=FlopocoScanner.h -o ExpressionScanner.cpp  $(top_srcdir)/src/FPExpressions/ExpressionScanner.l

FlopocoScanner.h:ExpressionScanner.cpp ExpressionParser.cpp

ExpressionParser.h:ExpressionParser.cpp
