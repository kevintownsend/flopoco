#ifndef FIXREALKCM_HPP
#define FIXREALKCM_HPP
#include "../Operator.hpp"

namespace flopoco{

	extern map<string, double> emptyDelayMap;

	class FixRealKCM : public Operator
	{
	public:
	  FixRealKCM(Target* target, int msbIn, int lsbIn, int lsbOut, string constant, map<string, double> inputDelays = emptyDelayMap);
	  ~FixRealKCM();

	  int msbIn;
	  int lsbIn;
	  int wOut;
	  int lsbOut;
	  string constant;
	};

}
#endif
