
Function: ??
PiecewiseFunction: ???

TODO: manage non-uniform errors; manage non-uniform segmentation 


FixFunction
    could be mostly useless. Currently centralizes Sollya parsing and evaluation to various formats.

FunctionTable: a flopoco operator that is the most basic implementation of a Function


BasicPolyApprox
				Two constructors:
				- the simple one that inputs only target accuracy, and manages all the rest
				- the one that input degree and computes approximation error, to be used by PiecewisePolyApprox

PiecewisePolyApprox: roughly the approximation part of previous FunctionEvaluator

FixPolynomialHornerEvaluator: an operator that inputs an X in [0,1] and a vector of coefficients, and evaluates the corresponding polynomial.

SimplePolyApproximator: a flopoco operator plugging  BasicPolyApprox to FixPolynomialHornerEvaluator

SimpleBitHeapApproximator: an Operator plugging BasicPolyApprox to a bit heap back-end.
													 (with bit heaps you don't do piecewise polynomials, or you do it the HOTBM way)


PiecewisePolyApproximator:
								the flopoco operator corresponding to the previous FunctionEvaluator


				
GenericEvaluator:
	A wrapper class written by Sylvain that should be able to instantiate any of the others.
