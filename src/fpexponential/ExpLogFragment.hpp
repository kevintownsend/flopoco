#ifndef EXPLOG_FRAGMENT_H
#define EXPLOG_FRAGMENT_H

#include "Fragment.hpp"

class ExpLogFragment : public Fragment
{
public:
  ExpLogFragment(Target* target, int length, Fragment* next_part = 0);
private:
  double area();
  double max_error(double input_error);
  void evalpos(int accuracy, int start, int& overlapping, bool& is_signed);
  void showinfo(int number);
  void generate(std::string prefix);
  int exp_bits, log_bits;
  int product_ibits1, product_ibits2;
};

#endif
