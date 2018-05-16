/* $Modified: Mon Dec 11 18:21:54 2006 by puwer $ */
#include "FourMomentum.h"
#ifndef _SQUAREDMATRIXELEMENTS_H_
#define _SQUAREDMATRIXELEMENTS_H_

#include <complex>

extern "C" {
void gggtthel_(double k1[4], double k2[4], double k3[4], double q1[4], 
	       double q2[4], double r1[4], double r2[4],
	       std::complex<double>& ammm, std::complex<double>& ammp,
	       std::complex<double>& ampm, std::complex<double>& ampp,
	       std::complex<double>& apmm, std::complex<double>& apmp,
	       std::complex<double>& appm, std::complex<double>& appp);

void qqgtthel_(double k1[4], double k2[4], double k3[4], double q1[4], 
	       double q2[4], double r1[4], double r2[4],
	       std::complex<double>& a1mm, std::complex<double>& a1mp,
	       std::complex<double>& a1pm, std::complex<double>& a1pp,
	       std::complex<double>& a2mm, std::complex<double>& a2mp,
	       std::complex<double>& a2pm, std::complex<double>& a2pp,
	       std::complex<double>& a3mm, std::complex<double>& a3mp,
	       std::complex<double>& a3pm, std::complex<double>& a3pp,
	       std::complex<double>& a4mm, std::complex<double>& a4mp,
	       std::complex<double>& a4pm, std::complex<double>& a4pp);
}

double gggtt2(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum & p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qqttg2(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum & p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qgttq2(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum & p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double gqbttqb2(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum & p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double gggtt2virt(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qqttg2virt(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double qgttq2virt(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double gqbttqb2virt(const FourMomentum & p1, const FourMomentum & p2, 
	      const FourMomentum &  p3,
	      const FourMomentum & kq, const FourMomentum & kqb);

double Ioperator_gg(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb);

double Ioperator_qq(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb);

double Ioperator_qg(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb);

double Ioperator_gqb(const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb);

#endif //_SQUAREDMATRIXELEMENTS_H_
