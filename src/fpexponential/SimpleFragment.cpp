#include <iostream>
#include <math.h>
#include <mpfr.h>
#include "../utils.hpp"
#include "signal.h"
#include "SimpleFragment.hpp"
#include "SimpleFragmentTable.hpp"

// TODO clean up
#define USE_MPFR 1
using namespace std;
using namespace FloPoCo::import::FPExp;
	extern vector<Operator *> oplist;



void simple_fragment_table_function(int accuracy, double x, bool negative, mpfr_t& result)
  {
    #ifdef USE_MPFR
      mpfr_t value;
      mpfr_init2(value, accuracy + 1);
      mpfr_init2(result, accuracy + (x >= 0));
      mpfr_set_d(value, x, GMP_RNDD);            // exact
      mpfr_exp(result, value, GMP_RNDN);         // arrondi au plus proche
      mpfr_set_d(value, 1 + x, GMP_RNDD);        // exact
      mpfr_sub(result, result, value, GMP_RNDD); // exact
      mpfr_clear(value);
    #else
      result = round(exp(x) - 1 - x, accuracy);
    #endif
  }


SimpleFragment::SimpleFragment(Target* target, int length, Fragment* next_part) :
	Fragment(target, length, next_part)
{
		setCopyrightString("X. Pujol (2007), C. Klein  (2008), F. de Dinechin (2009)");

		// The name is set after exploration is done
}

void SimpleFragment::evalpos(int accuracy, int start, int& overlapping, bool& is_signed)
{
  this->is_signed = is_signed;
  this->start = start;
  this->accuracy = accuracy;
  reallength = length + overlapping;
  end = start + reallength;

  // calcul des paramètres input_bits et output_bits de la table
  // (si output_bits <= 0, aucune table ne sera construite)
  input_bits = min(reallength, accuracy - 2 * start + 1);
  output_bits = accuracy - 2 * start; // prévoit large pour stocker e^x-1-x ~ x^2/2 (x premier morceau)
  exp_part1_start = start - 1;        // prévoit large pour stocker e^x-1 ~ x
  exp_part2_start = end - 1;          // prévoit pour e^morceaux_suivants-1

  // ajuste output_bits et exp_part1_start si on peut gagner un bit
  double max_input = invintpow2(start) - invintpow2(start + input_bits);
  mpfr_t max_output;

  simple_fragment_table_function(accuracy, max_input, false, max_output);
  #ifdef USE_MPFR
    if (mpfr_get_d(max_output, GMP_RNDD) < invintpow2(accuracy - output_bits + 1)) output_bits--;
    mpfr_clear(max_output);
  #else
    if (max_output < invintpow2(accuracy - output_bits + 1)) output_bits--;
  #endif
  if (accuracy - output_bits >= end) exp_part1_start++;

  // paramètres du produit (n'ont un sens que s'il y a un morceau après !)
  // bits significatifs du produit (de valeur >= 2 ^ -accuracy)
  // (si product_obits <= 0, il est inutile de faire le produit)
  product_obits = accuracy - exp_part1_start - exp_part2_start;
  // nombre de bits en entrée du produit
  product_ibits1 = (accuracy - exp_part1_start - end) + 1;
  product_ibits2 = (accuracy - start - exp_part2_start) + 1;

  overlapping = 0;
  // is_signed n'est pas modifié
}

void SimpleFragment::showinfo(int number)
{
  Fragment::showinfo(number);
  cout << ", normal method" << endl;
}

void SimpleFragment::generate(std::string prefix)
{
	if (next_part != 0) {
	  next_part->generate(prefix);
	  oplist.push_back(next_part);
	}

	ostringstream o;
	o << prefix << "_exp" << accuracy - start;
	setName(o.str());

	addInput("x", accuracy - start);
	addOutput("y", accuracy - start+1);


   if (next_part == 0) {
	  vhdl << tab << "-- Here e^x-1 is well approximated by x" <<endl;
	  vhdl << tab << "y <=  '0' & x ;" << endl;
	}
  else {
	  vhdl << tab << declare("part1", end-start) << " <= x" << range(accuracy-start-1, accuracy-end) << ";" << endl;
	  vhdl << tab << "-- using simple tabulation of e^x-1-x" <<endl;
	  SimpleFragmentTable* tbl = new SimpleFragmentTable(Operator::target_, (Fragment*)this);
	  oplist.push_back(tbl);
	  vhdl << tab << "-- expm1 of this part" << endl;
	  inPortMap(tbl, "X", "part1");
	  outPortMap(tbl, "Y", "tbl_out");
	  vhdl << instance(tbl, "component1");

	  int tbl_out_size = getSignalByName("tbl_out") -> width(); // probably equal to output_bits
	  vhdl << tab << declare("expm1_part1",accuracy-start) <<" <= ";
	  if (accuracy - output_bits < end)
		  vhdl << "TODO 1" << ';' << endl;  //;part_1.getPart(start - 1) << " + " << tbl_out.getPart(start - 1) << ';' << endl << endl;
	  else
		  vhdl << "part_1 & " << zeroGenerator((accuracy-start) -(end-start) - tbl_out_size,0) << " & " << use("tbl_out") << ";" << endl << endl;

	  vhdl << tab << declare("part2", accuracy-end) << " <= x" << range(accuracy-end-1, accuracy-accuracy) << ";" << endl;
	  inPortMap(next_part, "x", "part2");
	  outPortMap(next_part, "y", "expm1_part2");
	  vhdl << instance(next_part, "component2");

	  if (product_obits <= 0) {
      /* si on arrive là, il est probable que le découpage ne soit pas optimal
         (la dernière partie est trop petite) */
		  throw string("SimpleFragment.cpp : Not yet implemented, product_obits <= 0"); 
		  vhdl << "  y <= " ; //<< exp_part1.getPart(start - 1) << " + " << exp_part2.getPart(start - 1) << ';' << endl;
	  }
    else {
		 int expm1_part1_size = getSignalByName("expm1_part1") -> width();
		 int expm1_part2_size = getSignalByName("expm1_part2") -> width();
		 vhdl << tab << declare("product", product_ibits1+product_ibits2) << " <= " 
				<< "expm1_part1"<< range(expm1_part1_size-1, expm1_part1_size-product_ibits1) 
				<< "  * expm1_part2" << range(expm1_part2_size-1, expm1_part2_size-product_ibits2)<< ";" << endl;
		 vhdl << tab << "y <= " << "( '0' & expm1_part1"
				<< ") + (" << zeroGenerator(end-start,0) << " & expm1_part2)"
				<< ") + (" << zeroGenerator(accuracy-start+1 - product_obits, 0) << " &  + product"<< range(product_ibits1+product_ibits2-1, product_ibits1+product_ibits2-product_obits) << ");" << endl;
	 }
  }
}




double SimpleFragment::area()
{
  double result = Fragment::area();
  if (output_bits > 0)
    result += table_area(input_bits, output_bits);
  if (next_part != 0)
    result += multiplier_area(product_ibits1, product_ibits2);
  return result;
}

double SimpleFragment::max_error(double input_error)
{
  double result = Fragment::max_error(input_error);
  double fragment_error;

  /* hypothèse : on suppose que le dernier morceau est assez
     long pour que l'influence de l'erreur input_error sur
     le terme quadratique soit négligeable */

  if (output_bits <= 0)
    // approxime e^x par 1+x -> erreur en x^2/2 (termes suivants du DSE négligés)
    fragment_error = invintpow2(2 * start - accuracy + 1);
  else {
    // erreur sur la valeur tabulée de e^x-1-x
    fragment_error = 0.5;
    // éventuellement, erreur quand l'entrée est tronquée
    if (input_bits < reallength)
      /* x tronqué à input_bits bits. on pose e = 2 ^ -(start + input_bits)
	 erreur max en sortie : (x + e) ^ 2 - x ^ 2 = 2ex + e ^ 2
         input_bits choisi pour avoir 2ex <= 2 ^ -accuracy d'où */
      fragment_error += 1.0 + invintpow2((start + input_bits) * 2 - accuracy);
  }

  if (next_part == 0)
    result += fragment_error;
  else if (product_obits <= 0)
    result += fragment_error + 1.0;
  else {
    // majorant de l'exponentielle du premier morceau
    double max_output = exp(invintpow2(start)) - 1;
          /*  ErreurMax = ErreurMax + ErreurMorceau + 3 / 4 + 1 + _
                        ValeurMax * ErreurMax + _
                        (2 ^ -Morceaux(i + 1).Debut) * ErreurMorceau */
    result = result + fragment_error // addition des erreurs précédentes
             + 1.0                   // les entrées du produit sont tronquées
	     + max_output * result   // erreurs précédentes, multipliées
	     + invintpow2(next_part->getStart()) * fragment_error
	     + 1.0;                  // résultat du produit tronqué
  }

  return result;
}
