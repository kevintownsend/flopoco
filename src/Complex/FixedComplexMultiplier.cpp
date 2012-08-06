#include <fstream>
#include <sstream>
#include "FixedComplexMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator *> oplist;

	//TODO: explore implementation using multiply-accumulate operators
	//FIXME: correct timing of the circuit
	FixedComplexMultiplier::FixedComplexMultiplier(Target* target, int wI_, int wO_, bool signedOperator_, bool threeMultiplications)
		: Operator(target), wI(wI_), wO(wO_), signedOperator(signedOperator_)
	{
		
		ostringstream name;

		setCopyrightString ( "Matei Istoan, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixedComplexMultiplier_" << wI << "_" << wO << "_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "FixedComplexMultiplier_" << wI << "_" << wO << "_uid" << getNewUId();
		setName( name.str() );

		addInput("Xi", 		wI, true);
		addInput("Xr", 		wI, true);
		addInput("Yi", 		wI, true);
		addInput("Yr", 		wI, true);

	

#if 1
		// TODO add an option to exclude -4*-4
		addOutput("Zi",   wO, 2);
		addOutput("Zr",   wO, 2);

		// we compute the two products faithfully on wO bits
		// we add them, we obtain wO+1 bits
		// so after truncating the sum to wO bits the result is faithful
		int g = IntMultiplier::neededGuardBits(wI, wI, wO); 

		BitHeap* bitHeapRe = new BitHeap(this, 1+wO+g);  // will add XrYr - XiYi
		BitHeap* bitHeapIm = new BitHeap(this, 1+wO+g);  // will add XrYi + XiYr
		// Use virtual multipliers that will add their result to the bitHeap
		//IntMultiplier* multXrYr = 
		new IntMultiplier(this, bitHeapRe,
		                  getSignalByName("Xr"),
		                  getSignalByName("Yr"),
		                  wI, wI, wO, 
		                  g, // lsbWeight
		                  false, // negate
		                  signedOperator, 1.0);
		//IntMultiplier* multXiYi = 
		new IntMultiplier(this, bitHeapRe,
		                  getSignalByName("Xi"),
		                  getSignalByName("Yi"),
		                  wI, wI, wO, 
		                  g, // lsbWeight
		                  true, // negate
		                  signedOperator, 1.0);
		// The round bit
		bitHeapRe -> addConstantOneBit(g);
		bitHeapRe -> generateCompressorVHDL();	
		


		//IntMultiplier* multXrYi = 
		new IntMultiplier(this, bitHeapIm,
		                  getSignalByName("Xr"),
		                  getSignalByName("Yi"),
		                  wI, wI, wO, 
		                  g, // lsbWeight
		                  false, // negate
		                  signedOperator, 1.0);
		//IntMultiplier* multXiYr = 
		new IntMultiplier(this, bitHeapIm,
		                  getSignalByName("Xi"),
		                  getSignalByName("Yr"),
		                  wI, wI, wO, 
		                  g, // lsbWeight
		                  true, // negate
		                  signedOperator, 1.0);
		// The round bit
		bitHeapIm -> addConstantOneBit(g);

		//bitHeapIm -> generateCompressorVHDL();			

		vhdl << tab << "Zr <= " << bitHeapRe -> getSumName() << range(wO+g, wO+1) << ";" << endl;
		vhdl << tab << "Zi <= " << bitHeapIm -> getSumName() << range(wO+g, wO+1) << ";" << endl;
		
#else // pre-BitHeap version
		addOutput("Zi",   2*w, 2);
		addOutput("Zr",   2*w, 2);

		if(!threeMultiplications){
			IntMultiplier* multiplyOperator = new IntMultiplier(target, w, w, w, signedOperator, 1.0, inDelayMap("X",getCriticalPath()));
			oplist.push_back(multiplyOperator);
			IntAdder* addOperator =  new IntAdder(target, 2*w, inDelayMap("X",getCriticalPath()));
			oplist.push_back(addOperator);
			
			inPortMap (multiplyOperator, "X", "Xi");
			inPortMap (multiplyOperator, "Y", "Yi");
			outPortMap(multiplyOperator, "R", "XiYi");
			vhdl << instance(multiplyOperator, "MUL_XiYi");
			
			inPortMap (multiplyOperator, "X", "Xr");
			inPortMap (multiplyOperator, "Y", "Yr");
			outPortMap(multiplyOperator, "R", "XrYr");
			vhdl << instance(multiplyOperator, "MUL_XrYr");
			
			inPortMap (multiplyOperator, "X", "Xr");
			inPortMap (multiplyOperator, "Y", "Yi");
			outPortMap(multiplyOperator, "R", "XrYi");
			vhdl << instance(multiplyOperator, "MUL_XrYi");
			
			inPortMap (multiplyOperator, "X", "Xi");
			inPortMap (multiplyOperator, "Y", "Yr");
			outPortMap(multiplyOperator, "R", "XiYr");
			vhdl << instance(multiplyOperator, "MUL_XiYr");
			
			syncCycleFromSignal("XiYr", false);
			
			// invert the sign of XiYi to obtain a subtraction
			vhdl << tab << declare("neg_XiYi", 2*w) << " <= XiYi xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
			
			syncCycleFromSignal("neg_XiYi", false);
			nextCycle();
			
			inPortMap 	(addOperator, "X", 	 "XrYr");
			inPortMap 	(addOperator, "Y", 	 "neg_XiYi");
			inPortMapCst(addOperator, "Cin", "\'1\'");
			outPortMap	(addOperator, "R", 	 "Zr", false);
			vhdl << instance(addOperator, "ADD_XrYrMinXiYi");
			
			inPortMap 	(addOperator, "X", 	 "XrYi");
			inPortMap 	(addOperator, "Y", 	 "XiYr");
			inPortMapCst(addOperator, "Cin", "\'0\'");
			outPortMap	(addOperator, "R", 	 "Zi", false);
			vhdl << instance(addOperator, "ADD_XrYiAddXiYr");
		}
		else{
			try{
				IntMultiplier* multiplyOperator = new IntMultiplier(target, w, w, w, signedOperator, 1.0, inDelayMap("X",getCriticalPath()));
				oplist.push_back(multiplyOperator);
				IntAdder* addOperator =  new IntAdder(target, w, inDelayMap("X",getCriticalPath()));
				oplist.push_back(addOperator);	
				
				vhdl << tab << declare("neg_Yr", w) << " <= Yr xor (" << w-1 << " downto 0 => \'1\');" << endl;	
				
				inPortMap 	(addOperator, "X", 	 "Xr");
				inPortMap 	(addOperator, "Y",   "Xi");
				inPortMapCst(addOperator, "Cin", "\'0\'");
				outPortMap	(addOperator, "R", 	 "XrAddXi");
				vhdl << instance(addOperator, "ADD_XrXi");
			
				inPortMap 	(addOperator, "X", 	 "Yi");
				inPortMap 	(addOperator, "Y",   "neg_Yr");
				inPortMapCst(addOperator, "Cin", "\'1\'");
				outPortMap	(addOperator, "R",   "YiMinYr");
				vhdl << instance(addOperator, "ADD_YiMinYr");
			
				inPortMap 	(addOperator, "X",   "Yi");
				inPortMap 	(addOperator, "Y",   "Yr");
				inPortMapCst(addOperator, "Cin", "\'0\'");
				outPortMap	(addOperator, "R",   "YrAddYi");
				vhdl << instance(addOperator, "ADD_YrAddYi");
			
				syncCycleFromSignal("YrAddYi", false);
				//nextCycle(); 
			
				inPortMap (multiplyOperator, "X", "Yr");
				inPortMap (multiplyOperator, "Y", "XrAddXi");
				outPortMap(multiplyOperator, "R", "K1");
				vhdl << instance(multiplyOperator, "MUL_K1");
			
				inPortMap (multiplyOperator, "X", "Xr");
				inPortMap (multiplyOperator, "Y", "YiMinYr");
				outPortMap(multiplyOperator, "R", "K2");
				vhdl << instance(multiplyOperator, "MUL_K2");
			
				inPortMap (multiplyOperator, "X", "Xi");
				inPortMap (multiplyOperator, "Y", "YrAddYi");
				outPortMap(multiplyOperator, "R", "K3");
				vhdl << instance(multiplyOperator, "MUL_K3");
			
				syncCycleFromSignal("K3", false);
				//nextCycle(); 
			
				vhdl << tab << declare("neg_K3", w) << " <= K3 xor (" << w-1 << " downto 0 => \'1\');" << endl;
			
				syncCycleFromSignal("neg_K3", false);
				//nextCycle();
			
				IntAdder *addOperator2 =  new IntAdder(target, 2*w, inDelayMap("X",getCriticalPath()));
				oplist.push_back(addOperator2);
			
				inPortMap 	(addOperator2, "X",   "K1");
				inPortMap 	(addOperator2, "Y",   "neg_K3");
				inPortMapCst(addOperator2, "Cin", "\'1\'");
				outPortMap	(addOperator2, "R",   "Zr", false);
				vhdl << instance(addOperator2, "ADD_K1MinK3");
			
				inPortMap 	(addOperator2, "X",   "K1");
				inPortMap 	(addOperator2, "Y",   "K2");
				inPortMapCst(addOperator2, "Cin", "\'0\'");
				outPortMap	(addOperator2, "R",   "Zi", false);
				vhdl << instance(addOperator2, "ADD_K1AddK2");
			}catch(std::string str){
				cout << "execution interrupted: " << str << endl;
				exit(1);
			}
		}
#endif
	
	}	


	FixedComplexMultiplier::~FixedComplexMultiplier()
	{
	}
	
	
	void FixedComplexMultiplier::emulate ( TestCase* tc ) {
		mpz_class svXi = tc->getInputValue("Xi");
		mpz_class svYi = tc->getInputValue("Yi");
		mpz_class svXr = tc->getInputValue("Xr");
		mpz_class svYr = tc->getInputValue("Yr");
		
		
		if (! signedOperator){

			// mpz_class svZi = svXr*svYi + svXi*svYr;
			// mpz_class svZr = svXr*svYr - svXi*svYi;
			
			// // Don't allow overflow
			// mpz_clrbit ( svZi.get_mpz_t(), 2*wI );
			// mpz_clrbit ( svZr.get_mpz_t(), 2*w );

			// tc->addExpectedOutput("Zi", svZi);
			// tc->addExpectedOutput("Zr", svZr);
		}
		else{
			mpz_class big1I = (mpz_class(1) << (wI));
			mpz_class big1PI = (mpz_class(1) << (wI-1));
			mpz_class tmpSUB = (mpz_class(1) << (2*wI+1));

			if ( svXi >= big1PI)
				svXi = svXi - big1I;
			if ( svXr >= big1PI)
				svXr = svXi - big1I;

			if ( svYi >= big1PI)
				svYi = svYi - big1I;
			if ( svYr >= big1PI)
				svYr = svYr - big1I;
			
			mpz_class svZr = svXr*svYi + svXi*svYr;
			mpz_class svZi = svXr*svYr - svXi*svYi;
			
			if ( svZr < 0){
				svZr += tmpSUB; 
			}
			if ( svZi < 0){
				svZi += tmpSUB; 
			}
			
			// now truncate to wO bits
			if (wO<2*wI+1){
				svZr = svZr >> (2*wI+1-wO);
				svZi = svZi >> (2*wI+1-wO);
			}

			if (wO>2*wI+1){
				svZr = svZr << (-2*wI+1+wO);
				svZi = svZi << (-2*wI+1+wO);
			}
			tc->addExpectedOutput("Zi", svZi);
			tc->addExpectedOutput("Zr", svZr);
			
			svZr++;
			svZi++;
			mpz_clrbit ( svZr.get_mpz_t(), wO );			// no overflow
			mpz_clrbit ( svZi.get_mpz_t(), wO );			// no overflow
			tc->addExpectedOutput("Zi", svZi);
			tc->addExpectedOutput("Zr", svZr);
			
			
		}
		

	}

}
