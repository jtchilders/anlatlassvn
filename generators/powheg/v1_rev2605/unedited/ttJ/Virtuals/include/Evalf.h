/* $Modified: Thu Oct  5 11:27:23 2006 by puwer $ */
#ifndef _EVALF_H_
#define _EVALF_H_

#include "FourMomentum.h"
#include "ScalarInt.h"

#define SMEMAXQQ 86

enum PENTAGONS {ALLPENTAGONS=0,
		PENTA1,PENTA11,PENTA12,PENTA13,PENTA14,PENTA15,PENTA16,
		PENTA2,PENTA21,PENTA22,PENTA23,PENTA24,PENTA25,PENTA26,
		PENTA3,PENTA31,PENTA32,PENTA33,PENTA34,PENTA35,PENTA36,
		PENTA4,PENTA41,PENTA42,PENTA43,PENTA44,PENTA45,PENTA46};

void selectPentagon(int i);

void EvalLOfqq(double f[5][SMEMAXQQ], 
	     const FourMomentum & p1, const FourMomentum & p2,
	     const FourMomentum & p3,
	     const FourMomentum & kq, const FourMomentum & kqb);


void EvalLOfgg(double f[7][133], 
	     const FourMomentum & p1, const FourMomentum & p2,
	     const FourMomentum & p3,
	     const FourMomentum & kq, const FourMomentum & kqb);

void masscounterterm_gg(IntType f[7][133], 
			const FourMomentum & p1, const FourMomentum & p2, 
			const FourMomentum & p3,
			const FourMomentum & kq, const FourMomentum & kqb,
			double rterms
			);

void MasslessFermionLoops_gg(IntType f_nfl[12][133], 
			     const FourMomentum & p1, const FourMomentum & p2, 
			     const FourMomentum & p3,
			     const FourMomentum & kq, const FourMomentum & kqb,
			     IntType* Bptr[], IntType* Cptr[],
			     IntType* Dptr[],
			     double rterms
			     );

void TopLoops_gg(IntType f_nfl[12][133], 
		 const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb, 
		 IntType* Bptr[], IntType* Cptr[],
		 IntType* Dptr[],
		 double rterms
		 );

void RemainingDiags_gg(IntType f_nfl[12][133], 
		       const FourMomentum & p1, const FourMomentum & p2, 
		       const FourMomentum & p3,
		       const FourMomentum & kq, const FourMomentum & kqb,
		       IntType* Bptr[], IntType* Cptr[],
		       IntType* Dptr[], 
		       double rterms
		       );

void Ghosts_gg(IntType f_nfl[12][133], 
	       const FourMomentum & p1, const FourMomentum & p2, 
	       const FourMomentum & p3,
	       const FourMomentum & kq, const FourMomentum & kqb, 
	       IntType* Bptr[], IntType* Cptr[],
	       IntType* Dptr[],
	       double rterms
	       );


void pentagon_gg(IntType f_nfl[12][133], 
		 const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb,
		 double rterms
		 );

void masscounterterm_qq(IntType f[5][SMEMAXQQ], 
			const FourMomentum & p1, const FourMomentum & p2, 
			const FourMomentum & p3,
			const FourMomentum & kq, const FourMomentum & kqb,
			double rterms
			);

void MasslessFermionLoops_qq(IntType f_nfl[5][SMEMAXQQ], 
			     const FourMomentum & p1, const FourMomentum & p2, 
			     const FourMomentum & p3,
			     const FourMomentum & kq, const FourMomentum & kqb,
			     IntType* Bptr[], IntType* Cptr[],
			     IntType* Dptr[],
			     double rterms
			     );

void TopLoops_qq(IntType f_nfl[5][SMEMAXQQ], 
		 const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb, 
		 IntType* Bptr[], IntType* Cptr[],
		 IntType* Dptr[],
		 double rterms
		 );

void RemainingDiags_qq(IntType f_nfl[5][SMEMAXQQ], 
		       const FourMomentum & p1, const FourMomentum & p2, 
		       const FourMomentum & p3,
		       const FourMomentum & kq, const FourMomentum & kqb,
		       IntType* Bptr[], IntType* Cptr[],
		       IntType* Dptr[], 
		       double rterms
		       );

void Ghosts_qq(IntType f_nfl[5][SMEMAXQQ], 
	       const FourMomentum & p1, const FourMomentum & p2, 
	       const FourMomentum & p3,
	       const FourMomentum & kq, const FourMomentum & kqb, 
	       IntType* Bptr[], IntType* Cptr[],
	       IntType* Dptr[],
	       double rterms
	       );


void pentagon_qq(IntType f_nfl[5][SMEMAXQQ], 
		 const FourMomentum & p1, const FourMomentum & p2, 
		 const FourMomentum & p3,
		 const FourMomentum & kq, const FourMomentum & kqb,
		 double rterms
		 );


#endif //_EVALF_H_
