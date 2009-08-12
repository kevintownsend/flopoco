#ifndef EXP_SIMPLE_FRAGMENT_H
#define EXP_SIMPLE_FRAGMENT_H

#include "Fragment.hpp"

class SimpleFragment : public Fragment
{
public:
	SimpleFragment(Target* target, int length, Fragment* next_part = 0);
private:
  double area();
  double max_error(double input_error);
  void evalpos(int accuracy, int start, int& overlapping, bool& is_signed);
  void showinfo(int number);
  void generate(std::string prefix);
	//void write_tbl_declaration(std::string prefix, ostream& o);
	//void write_tbl_arch(std::string prefix, ostream& o);
  int input_bits, output_bits;
  int exp_part1_start, exp_part2_start;
  int product_ibits1, product_ibits2, product_obits;
};

#endif
