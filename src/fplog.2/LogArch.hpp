#ifndef LOGARCH_H
#define LOGARCH_H

#include "LogArch.hpp"
#include "FirstInvTable.hpp"
#include "FirstLogTable.hpp"
#include "SecondInvTable.hpp"
#include "OtherLogTable.hpp"

class LZOC;

// An instance of a logarithm architecture, with a bit-true simulator and VHDL output. 

class LogArch
{
 public:

  int wE;

  int wF;

  
  LogArch(int wE, int wF)  ;
  virtual ~LogArch() ;
  void check();
  void describe();
  double simulate(double x);
  void output(std::ostream& o, std::string name, bool pipelined);
  int size_in_LUTs();



  // The input sizes to the successive tables
  int a[42]; 

  // The intermediate precision: at step i, the exp part is bounded by
  //    1 <= m < 1+2^-p[i]
  int p[42]; 

  // The size of the product, and hence of the subword of Z[i] input to the mult
  int psize[42]; 

  // The total size of non-zero bits, will be limited to wF+g
  int s[42]; 

  // The numbers of bits to truncate to limit the size to wF+g
  int t[42];

  // size before truncation, should be equal to s[i]+t[i]
  int sbt[42];



  // The max number of stages
  int stages;


  // The guard bits 
  int gLog;

  // The table objects
  FirstInvTable* it0;
  SecondInvTable* it1;
  FirstLogTable* lt0;
  OtherLogTable* lt[42];


  // This is the size of the virtual mantissa:
  // 1.000(p[i] zeroes)00000xxxxxxx
  // the xxx which will be actually computed have p[i] bits less 
  int sfullZ[42]; 

 private:
  // A global boolean flag disabling the simulation of the fullsize mantissa
  // as soon as it would be more than 64 bits
  int fullSizeSim;

  // The target precision: numbers may be truncated so that their LSB has weight -target_prec 
  int target_prec;  

  // output function
  void showFP(ostream& o, double x);

  // LZOC
  LZOC *lzoc;
};


#endif
