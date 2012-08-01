#ifndef IntMultiplierS_HPP
#define IntMultiplierS_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Table.hpp"
#include "BitHeap.hpp"

#include "IntMultipliers/MultiplierBlock.hpp"

namespace flopoco {


		


	/** The IntMultiplier class, getting rid of Bogdan's mess.
	 */
	class IntMultiplier : public Operator {
	
	public:
		/** An elementary LUT-based multiplier, written as a table so that synthesis tools don't infer DSP blocks for it*/
		class SmallMultTable: public Table {
		public:
			int dx, dy, wO;
			bool signedIO;
			SmallMultTable(Target* target, int dx, int dy, int wO, bool signedIO=false );
			mpz_class function(int x);
		};

		/**
		 * The IntMultiplier constructor
		 * @param[in] target           the target device
		 * @param[in] wX             X multiplier size (including sign bit if any)
		 * @param[in] wY             Y multiplier size (including sign bit if any)
		 * @param[in] wOut         wOut size for a truncated multiplier (0 means full multiplier)
		 * @param[in] signedIO     false=unsigned, true=signed
		 * @param[in] ratio            DSP block use ratio
		 **/
		IntMultiplier(Target* target, int wX, int wY, int wOut=0, bool signedIO = false, float ratio = 1.0, map<string, double> inputDelays = emptyDelayMap);


		/**
		 * The virtual IntMultiplier constructor adds all the multiplier bits to some bitHeap, but no more.
		 * @param[in] parentOp      the Operator to which VHDL code will be added
		 * @param[in] bitHeap       the BitHeap to which bits will be added
		 * @param[in] x            a Signal from which the x input should be taken
		 * @param[in] y            a Signal from which the y input should be taken
		 * @param[in] wX             X multiplier size (including sign bit if any)
		 * @param[in] wY             Y multiplier size (including sign bit if any)
		 * @param[in] wOut         wOut size for a truncated multiplier (0 means full multiplier)
		 * @param[in] lsbWeight     the weight, within this BitHeap, corresponding to the LSB of the multiplier output. 
		 *                          Note that there should be enough bits below for guard bits in case of truncation.
		 * @param[in] negate     if true, the multiplier result is subtracted from the bit heap 
		 * @param[in] signedIO     false=unsigned, true=signed
		 * @param[in] ratio            DSP block use ratio
		 **/
		IntMultiplier (Operator* parentOp, BitHeap* bitHeap,  Signal* x, Signal* y, int wX, int wY, int wOut, int lsbWeight, bool negate, bool signedIO, float ratio);

		/** How many guard bits will a truncated multiplier need? Needed to set up the BitHeap of an operator using the virtual constructor */
		static int neededGuardBits(int wX, int wY, int wOut);


		/**
		 *  Destructor
		 */
		~IntMultiplier();

		/**
		 * The emulate function.
		 * @param[in] tc               a test-case
		 */
		void emulate ( TestCase* tc );

		void buildStandardTestCases(TestCaseList* tcl);



	protected:
	
	
		string PP(int i, int j, int uid=-1);
		string PPTbl( int i, int j, int uid=-1);
		string XY(int i, int j, int uid=-1);
		string heap( int i, int j);

		void buildLogicOnly();
		void buildTiling();

		void manageSignBeforeMult();            /**< to be called before either buildHeapLogicOnly or buildHeapTiling **/
		void manageSignAfterMult();            /**< to be called after either buildHeapLogicOnly or buildHeapTiling **/


		/**	builds the multiplier block with logic(smallMultTables) 
		 *@param uid is just a number which helps to form the signal names (for multiple calling of the method
		 )	*/
		// commented by F2D. Is it used?
		//	void buildHeapLogicOnly(MultiplierBlock* mul,int uid);

	
		/**	builds the logic block ( smallMultTables) 
		 *@param topX, topY -top right coordinates 
		 *@param botX, botY -bottom left coordinates 
		 *@param uid is just a number which helps to form the signal names (for multiple calling of the method
		 )	*/
		void buildHeapLogicOnly(int topX, int topY, int botX, int botY, int uid=-1);
		void buildHeapTiling();


		/** splits the inputs in chunks which fits the DSP-s, and send them to the vector of MultiplierBlocks in Bitheap**/
		void splitting(int horDSP, int verDSP, int wxDSP, int wyDSP,int restX,int restY);


		/* Stuff related to SVG output */

		void printConfiguration();
		void printAreaView();
		void printLozengeView();
		void drawLine(int wX, int wY, int wRez, int offsetX, int offsetY, int scalingFactor, bool isRectangle);
		void drawDSP(int i, int xT, int yT, int xB, int yB, int offsetX, int offsetY, int scalingFactor,  bool isRectangle);
		// commented out by F2D. Is it used?
		//	void drawDSPinclined(int i, int xT, int yT, int xB, int yB,  int offsetX, int offsetY, int turnaroundX, double inclinedCoeffB, double inclinedCoeffT);
		void drawTargetFigure(int wX, int wY, int offsetX, int offsetY, int scalingFactor, bool isRectangle);


		int wxDSP, wyDSP;               /**< the width for X/Y in DSP*/
		int wXdecl;                     /**< the width for X as declared*/
		int wYdecl;                     /**< the width for Y  as declared*/
		int wX;                         /**< the width for X after possible swap such that wX>wY */
		int wY;                         /**< the width for Y after possible swap such that wX>wY */
		int wOut;						/**<the size of the output*/
		int wFull;                      /**< size of the full product: wX+wY-1 if signed, wX+wY if unsigned */
		int wTruncated;                 /**< The number of truncated bits, wFull - wOut*/
		int g ;                         /**< the number of guard bits*/
		int maxWeight;                  /**< The max weight for the bit heap of this multiplier, wOut + g*/
		int weightShift;                /**< the shift in weight for a truncated multiplier compared to a full one,  wFull - maxWeight*/
		int signedIO;
		double ratio;
		double maxError;     /**< the max absolute value error of this multiplier, in ulps of the result. Should be 0 for untruncated, 1 or a bit less for truncated.*/  
		double initialCP;     /**< the initial delay, getMaxInputDelays ( inputDelays_ ).*/  
		ofstream fig;	
	private:
		bool useDSP;
		Operator* parentOp;  /**<  For a virtual multiplier, adding bits to some BitHeap, this is a pointer to the Operator that will provide the actual vhdl stream etc. */
		BitHeap* bitHeap;    /**<  The heap of weighted bits that will be used to do the additions */
		// TODO the three following variable pairs seem uglily redundant
		Signal* x;
		Signal* y; 
		string xname;
		string yname;
		string inputName1;
		string inputName2;
		bool negate;  /**< if true this multiplier computes -xy */
		int multiplierUid;
		void initialize();   /**< initialization stuff common to both constructors*/
		vector<DSP*> dsps;
	
	};

}
#endif
