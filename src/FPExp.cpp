#include <iostream>
#include <fstream>
#include <sstream>

#include "FPExp.hpp"
#include "FloFP.hpp"
#include "utils.hpp"

#include "fpexp/stdfragment.h"
#include "fpexp/explore.h"

using namespace std;

extern char fp_exp_template[];

void copy_file_section(ostream& o, istream& i);

FPExp::FPExp(Target* target, int wE, int wF)
	: wE(wE), wF(wF)
{
	/* Generate unique name */
	{
		std::ostringstream o;
		o << "FPExp_" << wE << "_" << wF;
		unique_name = o.str();
	}

	add_FP_input("x", wE, wF);
	add_FP_output("r", wE, wF);

	int explore_size = wF;
	int exponent_size = wE;

	f = explore(explore_size);
	if (!f) throw std::string("FPExp::FPExp(): No fragment");
}

FPExp::~FPExp()
{
}

// Overloading the virtual functions of Operator
void FPExp::output_vhdl(std::ostream& o, std::string name)
{
	Licence(o, "Cristian KLEIN (2008)");
	stringstream fp_exp, fixp_exp, fixp_exp_tbl;

	int result_length;
	double area, max_error;
	// génère le code pour l'exponentielle en virgule fixe
	result_length = f->prepare(area, max_error);
	f->generate(fixp_exp, fixp_exp_tbl);

	int g = intlog2(max_error) + 2;
	cout
		<< "Estimated area (for fixed-point part): " << area << endl
		<< "Maximum error: " << max_error << endl
		<< "Internal precision: " << result_length << endl
		<< "Precision: " << result_length - g << endl;

	std::string cstInvLog2, cstLog2;
	{
		mpz_class z;

		mpfr_t mp2, mp1, mp;
		mpfr_init2(mp, 2*(wE+wF+g));	// XXX: way too much precision
		mpfr_inits(mp1, mp2, 0);
		mpfr_set_si(mp1, 1, GMP_RNDN);
		mpfr_set_si(mp2, 2, GMP_RNDN);

		mpfr_log(mp, mp2, GMP_RNDN);
		mpfr_mul_2si(mp, mp, wE-1+wF+g, GMP_RNDN);
		mpfr_get_z(z.get_mpz_t(), mp, GMP_RNDN);
		{
			std::ostringstream o;
			o.fill('0');
			o.width(wE-1+wF+g);
			o << z.get_str(2);
			cstLog2 = o.str();
		}

		mpfr_mul_2si(mp, mp, -(wE-1+wF+g), GMP_RNDN);
		mpfr_div(mp, mp1, mp, GMP_RNDN);
		mpfr_mul_2si(mp, mp, wE+1, GMP_RNDN);
		mpfr_get_z(z.get_mpz_t(), mp, GMP_RNDN);
		{
			std::ostringstream o;
			o.fill('0');
			o.width(wE+1);
			o << z.get_str(2);
			cstInvLog2 = o.str();
		}

		mpfr_clears(mp1, mp2, mp, 0);
	}

	// produit l'entité principale à partir du fichier modèle
	stringstream template_file(fp_exp_template);
	copy_file_section(fp_exp, template_file);
	fp_exp
		<< "  component " << unique_name << " is" << endl
		<< "    generic ( wE : positive := " << wE << ";" << endl
		<< "              wF : positive := " << result_length - g << ";" << endl
		<< "              g : positive := " << g << " );" << endl;
	copy_file_section(fp_exp, template_file);
	fp_exp
		<< "entity " << unique_name << " is" << endl
		<< "    generic ( wE : positive := " << wE << ";" << endl
		<< "              wF : positive := " << result_length - g << ";" << endl
		<< "              g : positive := " << g << " );" << endl
		<< "    port ( x : in  std_logic_vector(2+wE+wF downto 0);" << endl
		<< "           r : out std_logic_vector(2+wE+wF downto 0));" << endl
		<< "end entity;" << endl
		<< endl
		<< "architecture arch of " << unique_name << " is" << endl
		<< "  constant cstInvLog2 : std_logic_vector(wE+1 downto 0) := \"" << cstInvLog2 << "\";" << endl
		<< "  constant cstLog2 : std_logic_vector(wE-1+wF+g-1 downto 0) := \"" << cstLog2 << "\";" << endl;
	copy_file_section(fp_exp, template_file);
	fp_exp << "  label2 : exp_" << result_length - 1 << endl;
	copy_file_section(fp_exp, template_file);

	o << fixp_exp_tbl.str() << fixp_exp.str() << fp_exp.str();
}

TestCaseList FPExp::generateStandardTestCases(int n)
{
	// TODO
	return TestCaseList();
}

TestCaseList FPExp::generateRandomTestCases(int n)
{
	Signal& sx = *get_signal_by_name("x");
	Signal& sr = *get_signal_by_name("r");
	Signal  sr_exc = (*get_signal_by_name("r")).getException();
	Signal  sr_sgn = (*get_signal_by_name("r")).getSign();
	Signal  sr_exp = (*get_signal_by_name("r")).getExponent();
	Signal  sr_man = (*get_signal_by_name("r")).getMantissa();

	TestCaseList tcl;	/* XXX: Just like Lyon's Transportion Company. :D */
	FloFP x(wE, wF), r(wE, wF);

	for (int i = 0; i < n; i++)
	{
		x = getLargeRandom(sx.width()-2) + (mpz_class(1) << (wE + wF + 1));
		r = x.exp();

		TestCase tc;
		tc.addInput(sx, x.getSignalValue());
		tc.addExpectedOutput(sr_exc, r.getExceptionSignalValue());
		tc.addExpectedOutput(sr_sgn, r.getSignSignalValue());
		if (r.getExceptionSignalValue() == 1)
		{
			/* Exp only returns faithful rounding */
			tc.addExpectedOutput(sr, r.getRoundedDownSignalValue());
			tc.addExpectedOutput(sr, r.getRoundedUpSignalValue());
		}
		tcl.add(tc);
	}


	return tcl;
}

/* Copie les lignes du flux i au le flux o jusqu'à
   trouver une ligne contenant uniquement un '$'
   ou à arriver à la fin du fichier */
void copy_file_section(ostream& o, istream& i)
{
  const int line_size = 1024;
  int empty_lines = 0;
  char line[line_size];
  while (!(i.eof() || empty_lines > 100)) {
    i.getline(line, line_size);
    if (line[0] == '$' && line[1] == '\0')
      break;
    else if (line[0] == '\0')
      empty_lines++;
    else
      empty_lines = 0;
    if (i.bad() || empty_lines > 100) {
      cerr << "Erreur lors de la lecture du fichier modele" << endl;
      exit(1);
    }
    o << line << endl;
  }
}
