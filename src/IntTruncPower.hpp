//#include "Operator.hpp"
//#include "utils.hpp"
#include <vector>
#include <list>
#include <tr1/unordered_map>

//using namespace flopoco;

/* some bit of context here:
 * there's b, a polynomial which is written b_1*1b-1 + b_2*1b-2 + ... +
 * b_n * 1b-n
 * and we try to calculate a power p of b, by evaluating the full
 * formula the simplifying it (see Jérémie Detrey's thesis pp. 40-41)
 */

/* actually MonomialOfBits will be used with a unique size, and using
 * different sizes together won't work
 * if we know this size at compile time it would be better as a template
 * argument */
class MonomialOfBits {
	public:
		MonomialOfBits (size_t n): data(std::vector<bool>(n)) {
		}
		MonomialOfBits (const MonomialOfBits& m) :n(m.n) {
			this->data = m.data;
		}
		MonomialOfBits operator* (const MonomialOfBits& rhs) const;
		MonomialOfBits& operator= (const MonomialOfBits& rhs) {
			this->n = rhs.n;
			this->data = rhs.data;
			return *this;
		}
		bool operator== (const MonomialOfBits& rhs) const {
			return (this->n == rhs.n) && (this->data == rhs.data);
		}
		std::vector<bool> data;
		size_t n;
};

class ProductBit;
/*
struct ProductBitIRElt {
	int times;
	ProductBit data;
};
*/
class ProductBitIR;
class ProductBit {
	public:
		ProductBit (): data (std::list<MonomialOfBits>()) {
		}
		std::list<MonomialOfBits> data;
};

class ProductBitIR {
	public:
		ProductBitIR (const ProductBit& rhs);
		ProductBitIR ()
			:data (std::tr1::unordered_map<MonomialOfBits,int>()) {
		}
		int getTimes (const MonomialOfBits& e) const;
		ProductBitIR& operator+= (const ProductBitIR& rhs);
		ProductBitIR operator* (const ProductBitIR& rhs);
		std::tr1::unordered_map<MonomialOfBits, int> data;
};

class ProductIR {
	public:
		ProductIR (size_t w, int msb) 
			:data (std::vector<ProductBitIR> (w, ProductBitIR()))
			,msb (msb) {
		}
		ProductIR& operator+= (const ProductIR& rhs);
		ProductIR& operator>>= (int n);
		ProductIR& operator*= (const ProductBitIR& rhs);
		ProductIR operator* (const ProductIR& rhs);
		std::vector<ProductBitIR> data;
		int msb;
};

