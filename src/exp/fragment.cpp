#include <iostream>
#include <fstream>
#include "math_lib.h"
#include "fragment.h"

using namespace std;

Fragment::Fragment(int length, Fragment* next_part) :
length(length), next_part(next_part)
{
}

int Fragment::prepare(double& area, double& max_error, bool with_sign, int accuracy, bool showinfo)
{
  Fragment* fragment;
  int totlength = totallength();
  if (accuracy == -1) accuracy = totlength + 1;

  // �value la position des morceaux
  int start = accuracy - totlength, overlapping = 0, number = 1;
  bool is_signed = with_sign;
  for (fragment = this; fragment != 0; fragment = fragment->next_part) {
    fragment->evalpos(accuracy, start, overlapping, is_signed);
    if (showinfo) fragment->showinfo(number++);
    start = fragment->end - overlapping;
  }

  area = this->area();
  max_error = this->max_error(2.5); // 2.5 � cause de la r�duction d'argument
  return accuracy;
}

void Fragment::generate(ostream& code_file, ostream& table_file)
{
  table_file << "library ieee;\n"
             << "use ieee.std_logic_1164.all;\n"
	     << "use ieee.std_logic_arith.all;\n"
	     << "use ieee.std_logic_unsigned.all;\n\n"
	     << "package pkg_exp_tbl is" << endl;
  write_tbl_declaration(table_file);
  table_file << "end package;\n" << endl;
  write_tbl_arch(table_file);

  code_file << "library ieee;\n"
            << "use ieee.std_logic_1164.all;\n"
	    << "use ieee.std_logic_arith.all;\n"
	    << "use ieee.std_logic_unsigned.all;\n\n"
	    << "package pkg_exp is" << endl;
  write_declaration(code_file);
  code_file << "end package;\n" << endl;
  write_arch(code_file);
}

int Fragment::totallength()
{
  if (next_part == 0)
    return length;
  else
    return length + next_part->totallength();
}

double Fragment::area()
{
  if (next_part == 0)
    return 0;
  else
    return next_part->area();
}

double Fragment::max_error(double input_error)
{
  if (next_part == 0)
    return input_error;
  else
    return next_part->max_error(input_error);
}

void Fragment::showinfo(int number)
{
  cout << "Partie " << number << " : bits " << start + 1 << " a " << end << ", ";
  if (is_signed)
    cout << "signee";
  else
    cout << "non signee";
}

void Fragment::write_declaration(std::ostream& o)
{
  if (next_part != 0) next_part->write_declaration(o);
  int input_size = accuracy - start;
  o << "  component exp_" << input_size << " is\n"
    << "    port (x : in  std_logic_vector(" << input_size << " - 1 downto 0);" << endl
    << "          y : out std_logic_vector(" << input_size + 1 << " - 1 downto 0)";
  if (is_signed)
    o << ';' << endl << "          sign : in std_logic";
  o << ");" << endl << "  end component;" << endl;
}

void Fragment::write_arch(std::ostream& o)
{
  if (next_part != 0) next_part->write_arch(o);
  int input_size = accuracy - start;

  o << "-- exp d'un nombre avec une mantisse de " << accuracy << " bits dont " << accuracy - start << " non nuls" << endl
    << "-- =============================================================" << endl << endl
    << "library ieee;" << endl
    << "use ieee.std_logic_1164.all;" << endl
    << "use ieee.std_logic_arith.all;" << endl
    << "use ieee.std_logic_unsigned.all;" << endl
    << "library work;" << endl
    << "use work.pkg_exp_tbl.all;" << endl;

  if (next_part != 0)
    o << "use work.pkg_exp.exp_" << accuracy - next_part->start << ';' << endl;

  o << endl
    << "entity exp_" << input_size << " is" << endl
    << "  port (x : in  std_logic_vector(" << input_size << " - 1 downto 0);" << endl
    << "        y : out std_logic_vector(" << input_size + 1 << " - 1 downto 0)";
  if (is_signed)
    o << ';' << endl << "        sign : in std_logic";
  o << ");" << endl << "end entity;" << endl << endl
    << "architecture arch of exp_" << input_size << " is" << endl;
}

void Fragment::write_tbl_declaration(std::ostream& o)
{
  if (next_part != 0) next_part->write_tbl_declaration(o);
}

void Fragment::write_tbl_arch(std::ostream& o)
{
  if (next_part != 0) next_part->write_tbl_arch(o);

  o << "-- tables pour le composant exp_" << accuracy - start << endl
    << "-- ===============================" << endl << endl;
}

double table_area(int input_bits, int output_bits)
{
  double base_size;

  if (input_bits <= 7)
    base_size = 1.0 * powOf2(input_bits) / 32.0;
  else if (input_bits == 8)
    base_size = 8.077;
  else if (input_bits == 9)
    base_size = 16.75;
  else
    base_size = 33.538 * powOf2(input_bits) / 1024.0;

  return base_size * static_cast<double>(output_bits);
}

double multiplier_area(int input1, int input2)
{
  return (input1 * input2 * 53) / 100;
}
