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
		 * @brief Standalone version of KCM. Input size will be msbIn-lsbIn+1
		 * @param target : target on which we want the KCM to run
		 * @param signedInput : true if input are considered as 2'complemented
		 * 						relative fixed point numbers
		 * 						false if they are considered as positive number
		 * @param msbin : power of two associated with input msb. For unsigned 
		 * 				  input, msb weight will be 2^msb, for signed input, it
		 * 				  will be -2^msb
		 * 	@param lsbIn : power of two of input least significant bit
		 * 	@param lsbOut : desired output precision i.e. output least 
		 * 					significant bit has a weight of 2^lsbOut
		 * 	@param constant : string that describes the constant with sollya
		 * 					  syntax
		 * 	@param targetUlpError : exiged error bound on result. Difference
		 * 							between result and real value should be
		 * 							lesser than targetUlpError * 2^lsbOut.
		 * 							Value has to be in ]0.5 ; 1] (if 0.5 wanted,
		 * 							please consider to create a one bit more
		 * 							precise KCM with a targetUlpError of 1 and
		 * 							truncate the result
		 */
		FixRealKCM(
				Target* target, 
				bool signedInput, 
				int msbIn, 
				int lsbIn, 
				int lsbOut, 
				string constant, 
				double targetUlpError = 1.0, 
				map<string, double> inputDelays = emptyDelayMap
			);

		/**
		 * @brief Incorporated version of KCM. 
		 * @param parentOp : operator frow which the KCM is a subentity
		 * @param target : target on which we want the KCM to run
		 * @param multiplicandX : signal which will be KCM input
		 * @param signedInput : true if input are considered as 2'complemented
		 * 						relative fixed point numbers
		 * 						false if they are considered as positive number
		 * @param msbin : power of two associated with input msb. For unsigned 
		 * 				  input, msb weight will be 2^msb, for signed input, it
		 * 				  will be -2^msb
		 * 	@param lsbIn : power of two of input least significant bit
		 * 	@param lsbOut : desired output precision i.e. output least 
		 * 					significant bit has a weight of 2^lsbOut
		 * 	@param constant : string that describes the constant with sollya
		 * 					  syntax
		 * 	@param bitHeap : bit heap on which the KCM should throw is result
		 * 	@param targetUlpError : exiged error bound on result. Difference
		 * 							between result and real value should be
		 * 							lesser than targetUlpError * 2^lsbOut.
		 * 							Value has to be in ]0.5 ; 1] (if 0.5 wanted,
		 * 							please consider to create a one bit more
		 * 							precise KCM with a targetUlpError of 1 and
		 * 							truncate the result
		 *
		 *  /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
		 * 	WARNING : nothing is added to the bitHeap when constant is set to
		 * 	zero. Your bitHeap must have been fed with another signal despite of
		 * 	what it will not produce an output signal (bitHeap problem)
		 *  /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
		 */
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
				int bitheapLsb,
				double targetUlpError = 1.0, 
				map<string, double> inputDelays = emptyDelayMap
			);


		// Overloading the virtual functions of Operator

		void emulate(TestCase* tc);

		static OperatorPtr parseArguments(
				Target* target,
				vector<string>& args
			);

		static void registerFactory();

		/**
		 * @brief determine how many guards bits will be needed for the KCM to
		 * 			compute result with asked precision
		 * @param wIn : input width
		 * @param targetUlpError : asked precision
		 * @param constant : string describing the constant with sollya syntax
		 * @param lsbIn : weight of input lsb
		 * @param lsbOur : weight of output lsb
		 */
		static int neededGuardBits(
				Target* target, 
				int wIn, 
				double targetUlpError,
				string constant,
				int lsbIn,
				int lsbOut
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
		static int guardBitsFromTableNumber(
					int nbTables,
					double targetUlpError
				);

		/**
		 * @brief handle operator initialisation and constant parsing
		 */
		void init();

		/**
		 *	@brief optimise KCM if constant is or a power of 2
		 *	@return true if a special case was handled (i.e. there is no need to 
		 *	create tables)
		 */
		bool handleSpecialConstant(Operator* op, string inputName = "X");

		/**
		 * @brief branch table output to bitHeap
		 * @param t : table array
		 * @param doSize : array of table output width such that t[i] has an
		 * output width of doSize[i]
		 * @param op : operator from which vhdl stream will be used as output
		 * stream
		 */
		void connectTablesToBitHeap(
				FixRealKCMTable** t,
				int* doSize,
				int nbOfTables,
				Operator* op
			);


		/**
		 * @brief create KCM Tables
		 * @param diSize : precomputed array of table input width
		 * @param nbOfTables : precomputed number of tables
		 * @param doSize_target : where data output width should be located when
		 * 					      function return
		 * @param op : operator from which tables are subentities
		 * @param inputSignalName : name of KCM input signal
		 * @return an array of KCM tables
		 */
		FixRealKCMTable** createTables(
				int* diSize,
				int nbOfTables,
				int** doSize_target,
				Operator* op,
				string inputSignalName = "X"
			);

		/**
		 * @brief compute the number of FixRealKCMTable that will be created by
		 * a KCM
		 * @param target : the target FPGA
		 * @param wIn : Input width
		 * @param msbC : constant most significant bit weight 
		 * @param lsbIn : input lsb weight
		 * @param lsbOut : output lsb weight (result precision)
		 * @param targetUlpError : required precision
		 * @param disize_target : where chunck sizes should be located when
		 * function return
		 */
		static int computeTablesNumber(
			Target* target,
			int wIn,
			int msbC, 
			int lsbIn,
			int lsbOut,
			double targetUlpError,
			int** disize_target = nullptr
		);

		int bitheaplsb;

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
					//if true, table result will be input*C - 1 in order to
					//avoid bit wasting for sign bit
					//e.g. : if cste is -2 and input is 0001 result will be 
					// 11101.
					// For sign extension if output width is 9 we need to add
					// 111100001 and that way we doesn't waste a useless sign bit
					// (as we know that se subproduct sign is the constant sign)
					bool negativeSubproduct, 
					bool last, 
					int logicTable = 1
					);

			mpz_class function(int x);
			FixRealKCM* mother;
			int index;
			//Weight of input lsb
			int weight;
			bool negativeSubproduct;
			bool last;
	};

}


#endif
