#ifndef FPPipeline_HPP
#define FPPipeline_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "Operator.hpp"
#include "HOTBM/sollya.h"	// Do NOT use libsollya from user's environment
#include "UtilSollya.hh"

namespace flopoco{

  /** The FPPipeline class.  */
  class FPPipeline : public Operator {
    
  public:
    /* TODO: Doxygen parameters*/ 
    FPPipeline(Target* target, string func, int wE, int wF);

    /**
     * FPPipeline destructor
     */
    ~FPPipeline();
			
  protected:
    int wE;   /**< Exponent size*/ 
    int wF;  /**< Significand fraction size */
    sollya_node_t expTree;  /**< Sollya tree of the expression */
  };

}
#endif
