// $Modified: Mon Sep 25 09:48:51 2006 by puwer $
void ggcorrelations(double clo[7][7], double cg1g2[7][7], double cqqb[7][7],
		    double cg1g3[7][7], double cg2g3[7][7], double cg1q[7][7],
		    double cg2q[7][7], double cg3q[7][7], double cg1qb[7][7],
		    double cg2qb[7][7], double cg3qb[7][7]
		    ){
  const double N = 3.0;
#include "colorcorr_gg.dec"
#include "colorcorr_gg.cpp"
}


void LOcolorgg(double CLO[7][7]){
  const double N=3.0;
#include "LO-color-gg.dec"
#include "LO-color-gg.cpp"
}

void NLOcolorgg(double CNLO[12][12]){
  const double N=3.0;
#include "NLO-color-gg.dec"
#include "NLO-color-gg.cpp"
}


void qqcorrelations(double clo[5][5], double cqqb[5][5], double cttb[5][5],
		    double cqg3[5][5], double cqbg3[5][5], double cqt[5][5],
		    double cqbt[5][5], double cg3t[5][5], double cqtb[5][5],
		    double cqbtb[5][5], double cg3tb[5][5]
		    ){
  const double N = 3.0;
#include "colorcorr_qq.dec"
#include "colorcorr_qq.cpp"
}


void LOcolorqq(double CLO[5][5]){
  const double N=3.0;
#include "LO-color-qq.dec"
#include "LO-color-qq.cpp"
}
