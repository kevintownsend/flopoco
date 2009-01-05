#ifndef FIRSTLOGTABLE_HPP
#define FIRSTLOGTABLE_HPP

#include "FirstInvTable.hpp"

class FirstLogTable : public Table
{
 public:
  FirstLogTable(int p1, int w1,  FirstInvTable* fit);
  
  ~FirstLogTable();

  FirstInvTable* fit;
  
  mpz_class function(int x);

  int    double2input(double x);

  double input2double(int x);
  
  mpz_class double2output(double x);

  double output2double(mpz_class x);

private:

};






#endif //FIRSTLOGTABLE_HPP

