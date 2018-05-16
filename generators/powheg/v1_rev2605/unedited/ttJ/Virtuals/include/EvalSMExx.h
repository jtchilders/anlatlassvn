/* $Modified: Thu Oct  5 11:24:53 2006 by puwer $ */
#ifndef _EVALSME_H_
#define _EVALSME_H_
#include "FourMomentum.h"
#include "Evalf.h"

enum {MINUS = 0, PLUS=1};

#ifdef __cplusplus

#include <complex>

void EvalSMEgg(std::complex<double>S[2][2][2][133],
	     const FourMomentum & p1, const FourMomentum & p2, 
	     const FourMomentum & p3, const FourMomentum & q1,
	     const FourMomentum & q2, const FourMomentum & r1, 
	     const FourMomentum & r2);

void EvalSMEqq(std::complex<double>S[2][2][SMEMAXQQ],
	     const FourMomentum & p1, const FourMomentum & p2, 
	     const FourMomentum & p3, const FourMomentum & q1,
	     const FourMomentum & q2, const FourMomentum & r1, 
	     const FourMomentum & r2);

#define XCOMPLEX std::complex<double>

extern "C" {

#else
#define XCOMPLEX _Complex double
#endif

void EvalSMEgg(XCOMPLEX S[2][2][2][133],
	     const double p1[4],  const double p2[4], 
	      const double p3[4],  const double q1[4],
	      const double q2[4],  const double r1[4], 
	      const double r2[4]);

#ifdef __cplusplus
}
#endif

#undef XCOMPLEX

#endif //_EVALSME_H_
