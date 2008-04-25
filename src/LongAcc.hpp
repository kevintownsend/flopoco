#ifndef LONGACC_HPP
#define LONGACC_HPP
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "Operator.hpp"
#include "Shifters.hpp"

class LongAcc : public Operator
{
public:
  LongAcc(Target* target, int wEX, int wFX, int MaxMSBX, int LSBA, int MSBA);
  ~LongAcc();

  void test_precision(int n);
  void test_precision2();

  int wEX; 
  int wFX; 
  int MaxMSBX;
  int LSBA; 
  int MSBA; 

  // Overloading the virtual functions of Operator
  void output_vhdl(ostream& o, string name);

private:
  Shifter* shifter;
  int sizeAcc;         //  = MSBA-LSBA+1;
  int sizeAccL; // used only for carry-select implementation
  int sizeSummand;     //  = MaxMSBX-LSBA+1;
  int sizeShiftedFrac; //  = sizeSummand + wFX;  to accomodate very small summands
  int maxShift;
  int E0X;
  int sizeShift;
  string summand2cname;
  int c2_chunk_size;
  int c2_pipeline_depth;

  int additionNumberOfChunks;
  int rebalancedAdditionChunkSize;
  int rebalancedAdditionLastChunkSize;

};


#endif
