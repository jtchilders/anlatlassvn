// $Modified: Fri Sep  3 16:03:12 2010 by uwer $
#include <iostream>
#include <iostream>
#include "Reduction.h"  
#include "GGReduction.h"  
#include "ScalarInt.h"

using namespace std;

int main(){
  /* Assume 10 different B-topologies */
  IntType* Bptr[10];
  
  double s12 = 100., m1 = 7., m2 = 8.;

  CoeffCache PVcoeffs;
  PVcoeffs.lookup(Bptr[1],s12,m1,m2);

  cout << Bptr[1][B0] << endl;
  cout << Bptr[1][B00] << endl;
  cout << Bptr[1][B11] << endl;
}
