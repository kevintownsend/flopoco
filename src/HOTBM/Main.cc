/* Using Java Coding style & tabs */

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <signal.h>

#include "Function.hh"
#include "Param.hh"
#include "HOTBMInstance.hh"
#include "Exhaustive.hh"

class NoMoreParameters {};
class InterruptedException {};

class Line {
	char *saveptr;
	bool newline;
	char *curline;

public:
	Line() {
		curline = 0;
	}

	void readLine() {
		string buf;
		getline(cin, buf);
		newline = true;
		if (curline == 0)
			free(curline);
		curline = strdup(buf.c_str());
	}
	
	string getString() {
		if (curline == 0)
			throw NoMoreParameters();
			
		char *tok = strtok_r(newline ? curline : 0, " \t\n", &saveptr);
		newline = false;
		if (tok == 0)
			throw NoMoreParameters();

		return string(tok);
	}

	int getInt() {
		return atoi(getString().c_str());
	}
};

Function *f = NULL;
HOTBMInstance *inst = NULL;

void parseBuild(Line line) {
	try {
		if (f) {
			delete f;
			f = NULL;
		}
		if (inst) {
			delete inst;
			inst = NULL;
		}

		string func;	  
		func = line.getString();
		f = new Function(func);

		int wI, wO, n;
		wI = line.getInt();
		if ((wI < 1) || (wI > 60)) throw NoMoreParameters();
		
		wO = line.getInt();
		if ((wO < 1) || (wO > 60)) throw NoMoreParameters();
		
		n = line.getInt();
		if (n < 0) throw NoMoreParameters();

		Param p(wI, wO, n);

		p.alpha = line.getInt();
		if ((p.alpha < 0) || (p.alpha > p.wI)) throw NoMoreParameters();
		p.beta = p.wI - p.alpha;

		for (int i = 0; i <= p.n; i++) {
			func = line.getString();
			if (((func != "rom") && (!i || (func != "powmult")))) throw NoMoreParameters();

			int alpha, beta;

			alpha = line.getInt();
			if ((alpha < 0) || (alpha > p.alpha)) throw NoMoreParameters();

			beta = line.getInt();
			if ((!i && beta) || (i && ((beta < 1) || (beta > p.beta)))) throw NoMoreParameters();

			if (func == "rom")
				p.t[i] = new TermROMParam(alpha, beta);
			else if (func == "powmult") {
				func = line.getString();
				if ((func != "rom") && (func != "adhoc")) throw NoMoreParameters();

				int mu, lambda;

				if (func == "adhoc") {
					mu = line.getInt();
					if ((mu < 0) || (mu > i*beta)) throw NoMoreParameters();
				}

				lambda = line.getInt();
				if ((lambda < 0) || (lambda > (func == "adhoc" ? mu : i*beta))) throw NoMoreParameters();

				PowerParam *pp;
				if (func == "rom")
					pp = new PowerROMParam(beta, lambda);
				else if (func == "adhoc")
					pp = new PowerAdHocParam(beta, mu, lambda);

				int mM, mT;
				
				mM = line.getInt();
				if (mM < 0) throw NoMoreParameters();

				mT = line.getInt();
				if (mT < 0) throw NoMoreParameters();

				int alphas[mM+mT], sigmas[mM+mT];
				for (int j = 0; j < mM+mT; j++) {
					alphas[j] = line.getInt();
					if ((alphas[j] < 0) || (alphas[j] > alpha)) throw NoMoreParameters();

					sigmas[j] = line.getInt();
					int sigmaSum = 0;
					for (int k = 0; k <= j; k++)
						sigmaSum += sigmas[k];
					if ((sigmas[j] < 1) || (sigmaSum > lambda)) throw NoMoreParameters();
				}

				p.t[i] = new TermPowMultParam(alpha, pp, mM, mT, alphas, sigmas);
				delete pp;
			}
		}

		HOTBMInstance::verbose = true;
		inst = new HOTBMInstance(*f, p);
		inst->roundTables();
	}
	catch (NoMoreParameters) {
		cerr << "Invalid build command." << endl << endl;
		cerr << "Usage: build <function> <w_I> <w_O> <n> <alpha> <term T_0> ... <term T_n>" << endl;
		cerr << "with:" << endl;
		cerr << "  <function> = log | sin" << endl << endl;
		cerr << "  <term T_k> =" << endl;
		cerr << "    as TermROM:     rom <alpha_k> <beta_k>" << endl;
		cerr << "    as TermPowMult: powmult <alpha_k> <beta_k> <power k> <mM_k> <mT_k>" << endl;
		cerr << "                            <alpha_{k,1}> <sigma_{k,1}> ..." << endl;
		cerr << "                            <alpha_{k,m_k}> <sigma_{k,m_k}>" << endl << endl;
		cerr << "  <power k> =" << endl;
		cerr << "    as PowerROM:    rom <lambda_k>" << endl;
		cerr << "    as PowerAdHoc:  adhoc <mu_k> <lambda_k>" << endl << endl;
	}
	catch (const char *s) {
		cerr << "Error: " << endl << s << endl;
		cerr << "Aborting." << endl;
	}
	catch (InterruptedException) {
		cerr << "Interrupted" << endl;
	}

}

void parseExplore(Line line) {
	try {
		if (f) {
			delete f;
			f = NULL;
		}
		if (inst) {
			delete inst;
			inst = NULL;
		}

		string func;
		func = line.getString();

		f = new Function(func);

		int wI, wO, n;

		wI = line.getInt();
		if ((wI < 1) || (wI > 60)) throw NoMoreParameters();

		wO = line.getInt();
		if ((wO < 1) || (wO > 60)) throw NoMoreParameters();

		n = line.getInt();
		if (n < 0) throw NoMoreParameters();

		Param p(wI, wO, n);
		Exhaustive ex(*f, p);
		inst = ex.getInstance();
		inst->roundTables();
	}
	catch (NoMoreParameters) {
		cerr << "Invalid explore command." << endl << endl;
		cerr << "Usage: explore <function> <w_I> <w_O> <n>" << endl;
		cerr << "with:" << endl;
		cerr << "  <function> = log | sin" << endl << endl;
	}
    catch (const char *s) {
        cerr << endl << endl << "The following error occured: " << endl;
        cerr << s << endl;
    }
	catch (InterruptedException) {
		cerr << "Interrupted" << endl;
	}
}

void parseDump(Line line) {
	try {
		if (!inst) {
			cerr << "No current HOTBMInstance available." << endl;
			return;
		}

		int n;
		list<int> dumpList;
		string buf;

		n = line.getInt();
		if (n < 1) throw NoMoreParameters();
		for (int i = 0; i < n; i++) {
			buf = line.getString();
			if (buf == "x")
				dumpList.push_back(HOTBM_DUMP_X);
			else if (buf == "function")
				dumpList.push_back(HOTBM_DUMP_FUNCTION);
			else if (buf == "approx")
				dumpList.push_back(HOTBM_DUMP_APPROX);
			else if (buf == "method")
				dumpList.push_back(HOTBM_DUMP_METHOD);
			else if (buf == "round")
				dumpList.push_back(HOTBM_DUMP_ROUND);
			else if (buf == "errpoly")
				dumpList.push_back(HOTBM_DUMP_ERR_POLY);
			else if (buf == "errmethod")
				dumpList.push_back(HOTBM_DUMP_ERR_METHOD);
			else if (buf == "errround")
				dumpList.push_back(HOTBM_DUMP_ERR_ROUND);
			else if (buf == "errmax")
				dumpList.push_back(HOTBM_DUMP_ERR_MAX);
			else if (buf == "inword")
				dumpList.push_back(HOTBM_DUMP_INPUT_WORD);
			else if (buf == "outword")
				dumpList.push_back(HOTBM_DUMP_OUTPUT_WORD);
			else throw NoMoreParameters();
		}

		buf = line.getString();
		ofstream dat(buf.c_str());
		inst->dump(dat, dumpList);
	}
	catch (NoMoreParameters n) {
		cerr << "Invalid dump command." << endl << endl;
		cerr << "Usage: dump <n> <opt_1> ... <opt_n> <file>" << endl;
		cerr << "with:" << endl;
		cerr << "  <opt_i> = x | function | approx | method | round | errpoly | errmethod" << endl;
		cerr << "          | errround | errmax | inword | outword" << endl << endl;
	}
}

volatile bool g_int;
void sighandler_int(int)
{
	signal(SIGINT, sighandler_int);
	g_int = true;
}
void handleSignals()
{
	if (g_int)
	{
		g_int = false;
		throw InterruptedException();
	}
}

int main(int argc, char **argv)
{
	g_int = false;
	signal(SIGINT, sighandler_int);

	try {
		/* Change the working mode of the FPU */
        /* Some tests fail without this */ 
#include <fpu_control.h>
		unsigned short cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
		_FPU_SETCW(cw);

		mpfr_set_default_prec(80);

		Line line;

		while (1) {	/* will use break */
			cerr << "HOTBM> ";
			
			/* Read next input line */
			line.readLine();
			
			string cmd;
			try {
			    /* Parse command */
			    cmd = line.getString();
			} catch (NoMoreParameters) {
			    /* The user has only hit one enter */
			    continue;
			}

			if (cmd == "quit") {
				cerr << "Exiting ..." << endl;
				break;
			}
			else if (cmd == "build") {
				parseBuild(line);
			}
			else if (cmd == "explore") {
				parseExplore(line);
			}
			else if (cmd == "dump") {
				parseDump(line);
			}
			else if (cmd == "tune") {
				try {
					if (!inst) {
						cerr << "No current HOTBMInstance available." << endl;
						continue;
					}
					inst->tune();
				}
				catch (const char *s) {
					cerr << "The following error occured: " << endl;
					cerr << s << endl;
				}
			}
			else if (cmd == "score") {
				if (!inst) {
					cerr << "No current HOTBMInstance available." << endl;
					continue;
				}
				cerr << "Score: ";
				cerr.precision(200); cerr << Exhaustive::score(*inst) << endl;
			}
			else if (cmd == "gen") {
				try {
					if (!inst) {
						cerr << "No current HOTBMInstance available." << endl;
						continue;
					}

					string name = line.getString();
					string file = line.getString();
					ofstream vhdl(file.c_str());
					inst->genVHDL(vhdl, name);
				}
				catch (NoMoreParameters) {
					cerr << "Invalid gen command." << endl << endl;
					cerr << "Usage: gen <name> <file>" << endl << endl;
				}
			}
			else if (cmd == "force") {
				HOTBMInstance::force = !HOTBMInstance::force;
				cerr << "Forcing is now " << (HOTBMInstance::force ? "enabled" : "disabled") << "." << endl;
			}
			else
				cerr << "Unknown command: '" << cmd << "'." << endl;
		}

		if (f)
			delete f;
		if (inst)
			delete inst;

		HOTBMInstance::freeApproxCache();
	}
	catch (NoMoreParameters) {
		cerr << "Internal error: NoMoreParameters not caught" << endl;
		return -1;
	}	
	catch (const char *s) {
		cerr << "Uncaught exception: " << endl << s << endl;
		cerr << "Aborting." << endl;
	}
	catch (InterruptedException) {
		cerr << "Interrupted in top-level." << endl;
	}

	return 0;
}
