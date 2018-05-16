/* $Modified: Mon Sep 25 09:49:31 2006 by puwer $ */
#ifndef _COLOR_H_
#define _COLOR_H_

void LOcolorgg(double CLO[7][7]);

void NLOcolorgg(double CLO[12][12]);

void ggcorrelations(double clo[7][7], double cg1g2[7][7], double cqqb[7][7],
                    double cg1g3[7][7], double cg2g3[7][7], double cg1q[7][7],
                    double cg2q[7][7], double cg3q[7][7], double cg1qb[7][7],
                    double cg2qb[7][7], double cg3qb[7][7]
                    );

void LOcolorqq(double CLO[5][5]);

void qqcorrelations(double clo[5][5], double cqqb[5][5], double cttb[5][5],
		    double cqg3[5][5], double cqbg3[5][5], double cqt[5][5],
		    double cqbt[5][5], double cg3t[5][5], double cqtb[5][5],
		    double cqbtb[5][5], double cg3tb[5][5]
		    );

#endif //_COLOR_H_
