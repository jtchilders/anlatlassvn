/* $Modified: Sat Oct 21 13:52:39 2006 by puwer $ */
#ifndef _FF_H_
#define _FF_H_
#include <cmath>
#include <complex>

extern "C" {
  extern struct {
    double delta;
  } ffcut_;

  void ffini_();
  
  void ffchangeflags_();
  
  void ffexi_();

  void ffxa0_(std::complex<double> & cint, const double & xd0, 
	      const double & xxmuq, const double & mtq, int & ierr);

  void ffxb0_(std::complex<double> & cint,const double & xd0,
	      const double & xxmuq, const double & xxk,
	      const double & xxma, const double & xxmb, int & ierr);

  void ffxc0_(std::complex<double> & cint, double xxpi[6], int & ierr);
 
  void ffxd0_(std::complex<double> & cint, double xxpi[13], int & ierr); 

  /*
   * d/dp^2 B_0(p^2,maq,mbq);
   */
  void ffxdb0_(std::complex<double> & cdb0,std::complex<double> & cdb0p,
	      const double & xp,const double & xma,const double & xmb,
	      int & ier);
}
#endif // _FF_H_
