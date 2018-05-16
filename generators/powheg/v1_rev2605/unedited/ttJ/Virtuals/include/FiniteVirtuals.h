// $Modified: Mon Nov 2 21:34 2010 by S. Alioli $
#include "FourMomentum.h"
#ifndef _FINITEVIRTUALS_H_
#define _FINITEVIRTUALS_H_

double gggtt2virtfin(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qqttg2virtfin(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qgttq2virtfin(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double gqbttqb2virtfin(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

#endif //_FINITEVIRTUALS_H_
