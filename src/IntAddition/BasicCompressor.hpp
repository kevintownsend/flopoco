#include "Operator.hpp"

#include "utils.hpp"

using namespace flopoco;


class BasicCompressor:public Operator {
	public:
	vector<int> height;


	public:
	BasicCompressor(Target * target, vector<int> h);

	~BasicCompressor() {
	};

	void emulate(TestCase * tc);
};
 
