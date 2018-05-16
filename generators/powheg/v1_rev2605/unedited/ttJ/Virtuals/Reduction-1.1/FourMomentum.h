/* $Modified: Tue Mar 17 16:31:37 2009 by uwer $ */
#ifndef _FOURMOMENTUM_H_
#define _FOURMOMENTUM_H_

#include <string>
#include <iostream>
#include <complex>
#include <vector>
#include <cassert>


//#define RAMBOEVENT

int spaversion();

class FourMomentum {

 private:
  static double tiny;
  double mom[4];
  double mass;
  double mass2;
  std::string name;
  bool ext;

  bool checkmass();

 public:
  FourMomentum(const std::string & s, double m) 
    : name(s), mass(m), mass2(m*m), ext(true) {
    mom[0] = mom[1] = mom[2] = mom[3] = 0.;
  };

  FourMomentum(const std::string & s) : name(s), ext(false) {
    mom[0] = mom[1] = mom[2] = mom[3] = 0.;
  };

  FourMomentum(const std::string & s, const double & p0,
	       const double & p1, const double & p2, 
	       const double & p3) : name(s), ext(false){
    
    mom[0] = p0;
    mom[1] = p1;
    mom[2] = p2;
    mom[3] = p3;
  };
  
  FourMomentum(const double & p0,const double & p1, 
	       const double & p2, const double & p3) : ext(false),
    name("tmp") {
    
    mom[0] = p0;
    mom[1] = p1;
    mom[2] = p2;
    mom[3] = p3;
    
  };
  
  FourMomentum(const double p[4]) : ext(false),
    name("tmp") {
    
    mom[0] = p[0];
    mom[1] = p[1];
    mom[2] = p[2];
    mom[3] = p[3];
    
  };
  
  FourMomentum(double m) : name("tmp"), ext(true), mass(m), mass2(m*m)
    {
      mom[0] = mom[1] = mom[2] = mom[3] = 0.;
    };
  
  FourMomentum() : name("tmp"), ext(false) {
    mom[0] = mom[1] = mom[2] = mom[3] = 0.;
  };
  
  double getmass() const {
    
    if (ext == true){
      return(mass);
    }
    
    return( sqrt( getmass2() ) );
    
  }
  
  double getmass2() const {
    
    if (ext == true){
      return(mass2);
    }
    
    double p = sqrt(fabs(mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3]) );
    
    double m2 = ( mom[0] - p ) * ( mom[0] + p );

    if ( fabs( mom[0] - p  ) < 1.e-15 * p ){
      return(0.0);
    }
    
    return( m2 );
    
  }
  
  
  double getComponent(int i) const {
    return(mom[i]);
  }
  
  double operator()(int i) const {
    return(mom[i]);
  } 
  
  const double* getArray() const {
    return(mom);
  }
  
  static void settiny(double x) {
    tiny = x;
    std::cout << "# FourMomentum::tiny set to " << x << std::endl;
  }
  
  void getF77Array(double p[4]) const {
    p[0] = mom[1];
    p[1] = mom[2];
    p[2] = mom[3];
    p[3] = mom[0];
    return;
  }
  
  double kt(const FourMomentum & ref) const;
  double kl(const FourMomentum & ref) const;
  double eta(const FourMomentum & ref) const;
  double cos(const FourMomentum & ref) const;
  double phi(const FourMomentum & ref1, const FourMomentum ref2) const;
  double y(const FourMomentum & ref) const;
  
  double getBeta() const;
  double norm3() const;
  
  FourMomentum parity() const {
    return( FourMomentum(mom[0],-mom[1],-mom[2],-mom[3]) );
  }
  
  void boost(const FourMomentum & p, double m = 0);
  
  void setFourMomentum(const FourMomentum & p);

  void setFourMomentum(const double p[4]);
  
  void setFourMomentum(const double p0,const double p1, 
		       const double p2, const double p3);
  
  void printMaple() const ;
  void printFortran() const;
  
  

  friend std::ostream& operator<<(std::ostream& os, const FourMomentum & p); 
  
  friend const FourMomentum operator-(const FourMomentum & p1);
  friend const FourMomentum & operator+(const FourMomentum & p1);
  friend const FourMomentum operator-(const FourMomentum & p1, 
				      const FourMomentum & p2);
  friend const FourMomentum operator+(const FourMomentum & p1, 
				      const FourMomentum & p2);
  friend const FourMomentum operator*(const double & s, 
				      const FourMomentum & p);
  friend const FourMomentum operator*(const FourMomentum & p,
				      const double & s );
  friend const FourMomentum operator/(const FourMomentum & p2,
				      const double & s1 );
  
  friend double dot3p(const FourMomentum & p1, const FourMomentum & p2);
  
  friend double dotp(const FourMomentum & p1, const FourMomentum & p2);
  
  friend void cross3p(const FourMomentum& p1, const FourMomentum& p2,
		      FourMomentum& p3);
  
  friend void rotate3(double angle, FourMomentum & n, FourMomentum& x);
  
  friend FourMomentum boost(const FourMomentum & p, const FourMomentum & x);
  
  friend std::complex<double> spa(const FourMomentum & p1, 
				  const FourMomentum & p2);
  
  friend std::complex<double> spb(const FourMomentum & p1, 
				  const FourMomentum & p2);
  
  friend void evalqr(const FourMomentum & kq, const FourMomentum & kqb, 
		     FourMomentum & q1, FourMomentum & q2, 
		     FourMomentum & r1, FourMomentum & r2);
  
  friend void SplitMomenta(const FourMomentum & kq, FourMomentum & q1, 
			   FourMomentum & q2);
  
};



int event(double s, double x[2], FourMomentum & p1, FourMomentum & p2,
	  double & jacobi );

int event(double s, double x[5], FourMomentum & p1, FourMomentum & p2,
	  FourMomentum & p3, double & jacobi );

int event(double s, double x[8], FourMomentum & p1, FourMomentum & p2,
	  FourMomentum & p3,FourMomentum & p4, double & jacobi );

int event(double s, double x[11], FourMomentum & p1, FourMomentum & p2,
	  FourMomentum & p3,FourMomentum & p4,FourMomentum & p5,
	  double & jacobi );

int eventn(const int n, const double s, const double x[], 
           FourMomentum kout[], double & jacobi);

int event(double x[], std::vector<FourMomentum>& kout, double & jacobi);

#ifdef RAMBOEVENT
void RamboEvent(double s, double x[8], FourMomentum & p1, FourMomentum & p2,
		FourMomentum & p3,FourMomentum & p4,
		double & jacobi );
#endif

#endif // _FOURMOMENTUM_H
