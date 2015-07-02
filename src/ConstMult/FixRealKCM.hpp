#ifndef FIXREALKCM_HPP
#define FIXREALKCM_HPP
#include "../Operator.hpp"
#include "../Table.hpp"

#include "../BitHeap/BitHeap.hpp"

namespace flopoco{


	class FixRealKCMTable;

	class FixRealKCM : public Operator
	{
	public:

		/**
		 * @brief Input size will be msbIn-lsbIn+1 if unsigned, msbIn-lsbIn+2 if
		 * signed.  
		 * \todo this is not the standard interface. msbIn should be
		 * the sign bit if signed. Fix this some day
		 */
		FixRealKCM(
				Target* target, 
				bool signedInput, 
				int msbIn, int lsbIn, 
				int lsbOut, 
				string constant, 
				double targetUlpError = 1.0, 
				map<string, double> inputDelays = emptyDelayMap
			);

		FixRealKCM(
				Operator* parentOp, 
				Target* target, 
				Signal* multiplicandX,
				bool signedInput, 
				int msbIn, 
				int lsbIn, 
				int lsbOut, 
				string constant,
				BitHeap* bitheap,
				double targetUlpError = 1.0, 
				map<string, double> inputDelays = emptyDelayMap
			);


		// Overloading the virtual functions of Operator

		void emulate(TestCase* tc);

		static int neededGuardBits(
				Target* target, 
				int wIn, 
				double targetUlpError
			);

		bool signedInput;
		int msbIn;
		int lsbIn;
		int wIn;
		int msbOut;
		int lsbOut;
		int wOut;
		string constant;
		float targetUlpError;
		mpfr_t mpC;
		mpfr_t absC;
		int msbC;
		int g;
		bool negativeConstant;

		/* The heap of weighted bits that will be used to do the additions */
		BitHeap*	bitHeap;    	

		/* The operator which envelops this constant multiplier */
		Operator*	parentOp;
		
		private:
		void init();
		void connectBitHeap(
				FixRealKCMTable** t,
				int* doSize,
				int nbOfTables
			);

		FixRealKCMTable** createTables(
				Target* target,
				int* diSize,
				int nbOfTables,
				int** doSize_target,
				string inputSignalName = "X"
			);
		static int computeTableNumbers(
			Target* target,
			int wIn,
			int** disize_target = nullptr
		);

	};

	class FixRealKCMTable : public Table
	{
		public:
			FixRealKCMTable(
					Target* target, 
					FixRealKCM* mother, 
					int i, 
					int weight, 
					int wIn, 
					int wOut, 
					//if true, all outputs will be prepared for fast sign
					//sign extension (I.E. msb will be inverted)
					bool fastSignExtend, 
					bool last, 
					int logicTable = 1
					);

			mpz_class function(int x);
			FixRealKCM* mother;
			int index;
			//Weight of input lsb
			int weight;
			bool fastSignExtend;
			bool last;
	};

}


#endif
