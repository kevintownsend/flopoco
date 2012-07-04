#include "Operator.hpp"

#include "utils.hpp"

using namespace flopoco;


class BasicCompressor:public Operator {
	public:
	static string operatorInfo;
	int w;


	public:
	BasicCompressor(Target * target, vector<int> h);

	~BasicCompressor() {
	};

	void emulate(TestCase * tc);
};
 
