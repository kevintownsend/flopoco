#ifndef OperatorPipeline_HPP
#define OperatorPipeline_HPP

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "Operator.hpp"

//this one is kept for the moment to have node class.
#include "FPExpressions/ExpressionParserData.h"


// #include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment
// #include "UtilSollya.hh"

namespace flopoco{

extern vector<Operator*> oplist;

/** The OperatorPipeline class.  */
class OperatorPipeline : public Operator {
public:
    /** Class that assembles floating-point operators starting with
             * an untyped, untimed Python-like description of the computational
             * datapath
             * @param[in] target     The target FPGA object
             * @param[in] filename   The filename containing the datapath
             * @param[in] wE         Exponent width
             * @param[in] wF         Fraction width
            **/
    OperatorPipeline(Target* target, string filename, int wE, int wF);

    /**
            * OperatorPipeline destructor
            */
    ~OperatorPipeline();

    /**
             * Function which generates the VHDL code containing the assembled
             * operators starting from the node containg the output variable
             * @param[in] n    The output variable (one of the output node list)
             * @param[in] top  Boolean describing if this function is called from
             * the statement level, or in some recursion level
            */
    void generateVHDL_c(node* n, bool top);


    /**
         * @brief this function run the optimisation algorithms for our tree
         */
    void optimise_tree();
protected:
    int wE;   /**< Exponent size*/
    int wF;  /**< Significand fraction size */
};

}
#endif
