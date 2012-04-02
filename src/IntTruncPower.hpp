//#include "Operator.hpp"
//#include "utils.hpp"
#include <vector>
#include <list>
#include <map>
#include <iostream>

//using namespace flopoco;

/* some bit of context here:
 * there's b, a polynomial which is written b_1*1b-1 + b_2*1b-2 + ... +
 * b_n * 1b-n
 * and we try to calculate a power p of b, by evaluating the full
 * formula the simplifying it (see Jérémie Detrey's thesis pp. 40-41)
 */

/* actually MonomialOfBits will be used with a unique size, and using
 * different sizes together won't work
 * if we knew this size at compile time it would be better as a template
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
		bool operator< (const MonomialOfBits& rhs) const;
		friend std::ostream& operator<<
			(std::ostream& o, const MonomialOfBits& m);

		std::vector<bool> data;
		size_t n;
};
std::ostream& operator<< (std::ostream& o, const MonomialOfBits& m);

class ProductBit;
/*
struct ProductBitIRElt {
	int times;
	ProductBit data;
};
*/
class ProductBitIR;
class Product;
class ProductBit {
	public:
		ProductBit (): data (std::list<MonomialOfBits>()) {
		}
		ProductBit (const ProductBitIR& rhs);

		std::list<MonomialOfBits> data;
};

class ProductBitIR {
	public:
		ProductBitIR (const ProductBit& rhs);
		ProductBitIR ()
			:data (std::map<MonomialOfBits,int>()) {
		}
		int getTimes (const MonomialOfBits& e) const;
		void addToCoeff (const MonomialOfBits& e, int a);
		ProductBitIR& operator+= (const ProductBitIR& rhs);
		ProductBitIR operator* (const ProductBitIR& rhs);
		friend std::ostream& operator<<
			(std::ostream& o, const ProductBitIR& pbi);

		std::map<MonomialOfBits, int> data;
};
std::ostream& operator<< (std::ostream& o, const ProductBitIR& pbi);

/* data is represented by a *little-endian* vector in the 2 following classes.
   why little ? because within simplify if we have a carry
   on the msb the vector has to grow on the msb side, which
   consequently has to be on the right as std::vectors only
   grow to the right (std::vector<T>::push_back but not push_front) */
class ProductIR {
	public:
		ProductIR (size_t w, int msb) 
			:data (std::vector<ProductBitIR> (w, ProductBitIR()))
			,msb (msb) {
		}
		ProductIR (const Product& rhs);
		ProductIR& operator+= (const ProductIR& rhs);
		ProductIR& operator>>= (int n);
		ProductIR& operator*= (const ProductBitIR& rhs);
		ProductIR operator* (const ProductIR& rhs);
		void simplify (void);
		friend std::ostream& operator<<
			(std::ostream& o, const ProductIR& pi);

		std::vector<ProductBitIR> data;
		int msb;
};
std::ostream& operator<< (std::ostream& o, const ProductIR& pi);

class Product {
	public:
		Product (const ProductIR& rhs);

		std::vector<ProductBit> data;
		int msb;
};

