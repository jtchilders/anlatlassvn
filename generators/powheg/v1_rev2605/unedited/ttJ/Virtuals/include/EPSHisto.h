// Uncomment to plot histograms of EPS frequency
//#define TESTBADPOINTS

#ifdef TESTBADPOINTS
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <sstream>
#include "Histogram.h"

using namespace std;

class EPSHisto {

 private :

  string name;
  ofstream outfile;
  Histo1d hist;
  int ncalls;
  int nbins;
  double xmin;
  double xmax;
  
 
 public:

  EPSHisto(string name_) {
    name=name_;
    ncalls=0;
    xmin=-10.;
    xmax=10.;
    nbins=200;
    hist.reset(name.c_str(),xmin,xmax,nbins);
    cout<<"# EPS histogram "<<name<<" opened"<<endl;
  }

  ~EPSHisto(){
    hist.rescale(1./ncalls);
    hist.finalize((double) ncalls);
    outfile.open(name.c_str());
    hist.printfile(outfile);
    outfile.close();
    cout<<"# EPS file "<<name<<" closed"<<endl;
  };
 
  void fill(double x_) {
    hist.fill(x_);
    hist.end_correlated_calls();
    ncalls++;
  }
};
  
#endif
