#include <fstream>
#include <sstream>
#include "FixedComplexMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator *> oplist;

	//TODO: explore implementation using multiply-accumulate operators
	//FIXME: correct timing of the circuit
	FixedComplexMultiplier::FixedComplexMultiplier(Target* target, int wI_, int wF_, bool signedOperator_, bool hasLessMultiplications)
		: Operator(target), wI(wI_), wF(wF_), signedOperator(signedOperator_)
	{
		signedOperator ? w = 1 + wI + wF : w = wI + wF;
		
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "FixedComplexMultiplier_" << w << "_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "FixedComplexMultiplier_" << w << "_uid" << getNewUId();
		setName( name.str() );

		addInput("Xi", 		w, true);
		addInput("Xr", 		w, true);
		addInput("Yi", 		w, true);
		addInput("Yr", 		w, true);
		addOutput("Zi",   2*w, 2);
		addOutput("Zr",   2*w, 2);

	

#if 1

		int g=IntMultiplier::neededGuardBits(w, w, w); 

		BitHeap* bitHeapRe = new BitHeap(this, w+g);  // will add XrYr - XiYi
		BitHeap* bitHeapIm = new BitHeap(this, w+g);  // will add XrYi + XiYr
		// a virtual multiplier that will add its result to the bitHeap
		IntMultiplier* multXrYr = new IntMultiplier(this, bitHeapRe, 
		                                            getSignalByName("Xr"),
		                                            getSignalByName("Xr"),
		                                            w, w, w, 
		                                            g, // lsbWeight
		                                            false, // subtract
		                                            signedOperator, 1.0);
		IntMultiplier* multXiYi = new IntMultiplier(this, bitHeapRe, 
		                                            getSignalByName("Xr"),
		                                            getSignalByName("Xr"),
		                                            w, w, w, 
		                                            g, // lsbWeight
		                                            true, // subtract
		                                            signedOperator, 1.0);
		
#else // pre-BitHeap version

		if(!hasLessMultiplications){
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

			mpz_class svZi = svXr*svYi + svXi*svYr;
			mpz_class svZr = svXr*svYr - svXi*svYi;
			
			// Don't allow overflow
			mpz_clrbit ( svZi.get_mpz_t(), 2*w );
			mpz_clrbit ( svZr.get_mpz_t(), 2*w );

			tc->addExpectedOutput("Zi", svZi);
			tc->addExpectedOutput("Zr", svZr);
		}else{
			mpz_class big1 = (mpz_class(1) << (w));
			mpz_class big1P = (mpz_class(1) << (w-1));
			mpz_class big2 = (mpz_class(1) << (w));
			mpz_class big2P = (mpz_class(1) << (w-1));

			if ( svXi >= big1P)
				svXi = svXi - big1;
			if ( svXr >= big1P)
				svXr = svXi - big1;

			if ( svYi >= big2P)
				svYi = svYi - big2;
			if ( svYr >= big2P)
				svYr = svYr - big2;
			
			mpz_class svXrYr = svXr*svYr;
			mpz_class svXiYi = svXi*svYi;
			mpz_class svXrYi = svXr*svYi;
			mpz_class svXiYr = svXi*svYr;
			
			if ( svXrYr < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXrYr = tmpSUB + svXrYr; 
			}
			if ( svXiYi < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXiYi = tmpSUB + svXiYi; 
			}
			if ( svXrYi < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXrYi = tmpSUB + svXrYi; 
			}
			if ( svXiYr < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXiYr = tmpSUB + svXiYr; 
			}
			
			mpz_class svZi = svXrYi + svXiYr;
			mpz_class svZr = svXrYr - svXiYi;
			
			//			cout << "Call to emulate() in class FixedComplexMultiplier. Assigning Zi with value: " << svZi << endl;
			//			cout << "Call to emulate() in class FixedComplexMultiplier. Assigning Zr with value: " << svZr << endl;
			 
			// Don't allow overflow
			mpz_clrbit ( svZi.get_mpz_t(), 2*w );
			mpz_clrbit ( svZr.get_mpz_t(), 2*w );
			
			//			cout << "Call to emulate() in class FixedComplexMultiplier. Changing Zi with value : " << svZi << endl;
			tc->addExpectedOutput("Zi", svZi);
			//			cout << "Call to emulate() in class FixedComplexMultiplier. Assigning Zr with value: " << svZr << endl;
			tc->addExpectedOutput("Zr", svZr);
		}
		

	}

}
