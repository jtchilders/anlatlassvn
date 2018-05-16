/* $Modified: Wed Feb  3 09:36:24 2010 by uwer $ */
#ifndef GGREDUCTION_H_
#define GGREDUCTION_H_

#define NDEBUG

class Topology;
class data;

#include <vector>
#include <map>
#include <cfloat>
#include <cassert>
#include <string>
#include <stdexcept>
#include "FourMomentum.h"
#include "Matrix.h"

#include "TopologyData.h"
#include "Matrix.h"

#define RED "\033[31m"
#define GREEN "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define RESET "\033[0m"

#define WITHCACHE 
#define WITHEXCEPTIONS

#define ONLYREALPART

#ifdef ONLYREALPART
typedef double IntType;
#else
typedef std::complex<double> IntType;
#endif

typedef std::map<unsigned int,IntType> coeffcache;



class GGException : public std::runtime_error {

 private:
  static int count;

 public:
  GGException(const std::string & msg = "") : std::runtime_error(msg) {
    count++;
  }

  void print();

};



class Topology {

 private:

  data *mydata;
  std::vector<Topology*> subtopos;
  coeffcache cache;
  /*
   *  Make the constructor private so that we can easily keep control over
   *  all created Topologies 
   */
  Topology(){}

 public:
  
  static Topology* factory(const Matrix& Sin);

  ~Topology();

  static void cacheclear();

  Topology* evalSubtopo(int i);
    

  inline double getB() const {
    return mydata->B;
  }

  inline int geti() const {
    return mydata->imin;
  }

  double getZ(int i) const;

  double getRange(int i) const;

  inline double getDetS() const {
    return mydata->detS;
  }

  /*  const std::vector<FourMomentum> & getMomenta() const {
    return(q);
  }

  const  std::vector<double> & getMasses() const {
    return(m);
  }
  */

  double getMass2(int i) const {
    assert( i>0 && i <= mydata->getNpoint() );
    return( - (*(mydata->S))(i,i) / 2.0 );
  }

  
  inline double getSinv(int i, int j) const {
    return ((*(mydata->Sinv))(i,j));
  }

  inline double getS(int i, int j) const {
    return ((*(mydata->S))(i,j));
  }

  Topology* subtopo(int i) const {
    assert(i>0 && i <= mydata->getNpoint() );
    return( subtopos[i-1]);
  }

  void printS(){
    (*(mydata->S)).printMaple();
  }

  friend IntType Int(int d, int n1, int n2, int n3, int n4, int n5,int n6, 
		     Topology* topo);

  friend IntType Int(int d, int n1, int n2, int n3, int n4, int n5, 
		     Topology* topo);

  friend IntType Int(int d, int n1, int n2, int n3, int n4, Topology* topo);

  friend IntType Int(int d, int n1, int n2, int n3, Topology* topo);

};

Topology pentagon(FourMomentum q1, FourMomentum q2, FourMomentum q3,
		  FourMomentum q4, FourMomentum q5, 
		  double m1, double m2, double m3, double m4, double m5);

void pentagon(Topology* & Etopo, const double s12, const double s13, 
	      const double s14, const double s15, const double s23,
	      const double s24, const double s25, const double s34,
	      const double s35, const double s45, const double m1q,
	      const double m2q, const double m3q, const double m4q, 
	      const double m5q);

void hexagon(Topology* & Ftopo, 
	     const double s12, const double s13, const double s14, 
	     const double s15, const double s16, 
	     const double s23, const double s24, const double s25, 
	     const double s26,
	     const double s34, const double s35, const double s36, 
	     const double s45, const double s46,
	     const double s56,
	     const double m1q, const double m2q, const double m3q, 
	     const double m4q, const double m5q, const double m6q);

Topology box(FourMomentum q1, FourMomentum q2, FourMomentum q3,
	     FourMomentum q4, double m1, double m2, double m3, double m4);

data dbox(FourMomentum q1, FourMomentum q2, FourMomentum q3,
	     FourMomentum q4, double m1, double m2, double m3, double m4);

Topology triangle(FourMomentum q1, FourMomentum q2, FourMomentum q3,
		  double m1, double m2, double m3);

Topology bubble(FourMomentum q1, FourMomentum q2,
		  double m1, double m2);

void clearCache();

IntType Int(int d, int n1, int n2, int n3, int n4, int n5, int n6,
	    Topology* topo);

IntType Int(int d, int n1, int n2, int n3, int n4, int n5, Topology* topo);

IntType Int(int d, int n1, int n2, int n3, int n4, Topology* topo);

IntType Int(int d, int n1, int n2, int n3, Topology* topo);

IntType Int(int d, int n1, int n2, double pq, double m1q, double m2q);

IntType Int(int d, int n1, double mq);

void setBCUT(double value);
double setbcut(double value);

#undef NDEBUG
#endif //TOPOLOGY_H_
