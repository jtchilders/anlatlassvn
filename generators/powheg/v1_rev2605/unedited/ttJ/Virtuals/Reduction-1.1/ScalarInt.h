/* $Modified: Fri Jul 21 11:04:41 2006 by puwer $ */
#ifndef _SCALARINT_H_
#define _SCALARINT_H_

#include "Reduction.h"


void initIntegralLibrary();
//added by SA
void finishIntegralLibrary();

IntType ScalarInt(const double mq);
IntType ScalarInt(double pq, double m0q, double m1q);
IntType ScalarInt(const Binput & in);
IntType ScalarInt(const Cinput & in);
IntType ScalarInt(const Dinput & in);

IntType ScalarInt(const std::vector<FourMomentum>& q, 
		  const std::vector<double>& m);

IntType A(double m);
IntType B(double pq,double m0,double m1);
IntType C(double p1q, double p2q, double p5q, 
	  double m0q,double m1q, double m2q);
IntType D(double p1q, double p2q, double p3q, double p4q, 
	  double p5q, double p7q, 
	  double m0q,double m1q, double m2q, double m3q);

IntType diffB0(double pq, double m0, double m1);
#endif // _SCALARINT_H_
