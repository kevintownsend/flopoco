#ifndef _TAYLOR_HH_
#define _TAYLOR_HH_

#include "Function.hh"
#include "Polynomial.hh"

using namespace std;

class Taylor {
public:
  Taylor(Function &f, int d, double x0);
  ~Taylor();

  Polynomial &getPoly() const;

private:
  Polynomial *poly;
};

#endif // _TAYLOR_HH_
