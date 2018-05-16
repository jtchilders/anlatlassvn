// F77 interface to C++ class mt2_bisect
//
// mt2_bisect by Cheng and Han can be downloaded from:
// http://particle.physics.ucdavis.edu/hefti/projects/doku.php?id=wimpmass
// (you need the files mt2_bisect.h and mt2_bisect.cpp, see link in section (c):
//  http://particle.physics.ucdavis.edu/hefti/projects/lib/exe/fetch.php?media=mt2-1.01a.tar.gz )
//
// this is just a wrapper function (AvM)
// arguments: [description from mt2_bisect.cpp]
// array pa[0..2], pb[0..2], pmiss[0..2] contains (mass,px,py)
// for the visible particles and the missing momentum. pmiss[0] is not used.


#ifndef MT2_BISECT_WRAP_H
#define MT2_BISECT_WRAP_H

#include "mt2_bisect.h"

extern "C" {

double calc_mt2_(double pa[], double pb[], double pmiss[],
		           const double* mn) {
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta( pa, pb, pmiss );
  mt2_event.set_mn( *mn );
  return mt2_event.get_mt2();
}

} // extern "C"

#endif
