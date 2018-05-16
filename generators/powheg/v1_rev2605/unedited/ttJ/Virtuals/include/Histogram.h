#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#define BINTYPE double

class Histo1d {

 private:
  double xmin,xmax,binsize;
  unsigned long  bincount;
  bool initialized;

  double* current_result;
  double* current_correlated_result;
  double* current_result2; 
  double* accumulated_result;
  double* accumulated_relative_inverse_error2;

  std::string name;
  std::string fname;

  int getBin(double x);
  //const Histo1d & operator=(const Histo1d & h);
  //const Histo1d & operator+=(const Histo1d & h);
  //  friend const Histo1d operator+(const Histo1d & h1, const Histo1d  & h2);
  Histo1d(const Histo1d & ori );

 public:
  Histo1d():initialized(false){};
  Histo1d(std::string name_, double xmin_, double xmax_, int bincount_) 
    : initialized(false) {
    reset(name_, xmin_, xmax_, bincount_);
  };
  ~Histo1d();
  
  void operator++(int);

  void reset(std::string name_, double xmin_, double xmax_, int bincount_);
  void reset();

  void setHist(int i);

  double getBinsize() const;

  void fill(double x);
  void fill(double x, double w);
  void fill2(double x, double w);

  void end_correlated_calls();
  void rescale(double x);
  void finalize(double calls);

  void print() ;
  void printraw();
  void printfile(std::ofstream&);

  double getLeftBound(unsigned long bin) const;
  double getRightBound(unsigned long bin) const ;
  double getCenter(unsigned long bin) const;

};
#endif // HISTOGRAM_H_
