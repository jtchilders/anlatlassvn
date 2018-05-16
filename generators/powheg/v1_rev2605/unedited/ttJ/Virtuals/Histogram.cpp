#include "Histogram.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cassert>

using namespace std;


Histo1d::~Histo1d(){
  delete[] current_result;
  delete[] current_correlated_result;
  delete[] current_result2;
  delete[] accumulated_result;
  delete[] accumulated_relative_inverse_error2;
}

void Histo1d::operator++(int){
}

void Histo1d::reset(std::string name_, double xmin_, double xmax_, 
		    int bincount_){
  name = name_;
  xmin = xmin_;
  xmax = xmax_;
  bincount = bincount_;

  binsize = (xmax-xmin) / bincount;

  if (initialized) {
    delete[] current_result;
    delete[] current_correlated_result;
    delete[] current_result2;
    delete[] accumulated_result;
    delete[] accumulated_relative_inverse_error2;
  }

  current_result = new BINTYPE[bincount];
  current_correlated_result = new BINTYPE[bincount];
  current_result2 = new BINTYPE[bincount];
  accumulated_result = new BINTYPE[bincount];
  accumulated_relative_inverse_error2 = new BINTYPE[bincount];


  ::fill(current_result, current_result + bincount, 0.0);
  ::fill(current_correlated_result, current_correlated_result + bincount, 0.0);
  ::fill(current_result2, current_result2 + bincount, 0.0);
  ::fill(accumulated_result, accumulated_result + bincount, 0.0);
  ::fill(accumulated_relative_inverse_error2, 
	 accumulated_relative_inverse_error2 + bincount, 0.0);

  initialized = true;
}


void Histo1d::reset(){

  ::fill(current_result, current_result + bincount, 0.0);
  ::fill(current_correlated_result, current_correlated_result + bincount, 0.0);
  ::fill(current_result2, current_result2 + bincount, 0.0);
  ::fill(accumulated_result, accumulated_result + bincount, 0.0);
  ::fill(accumulated_relative_inverse_error2, 
	 accumulated_relative_inverse_error2 + bincount, 0.0);

}

void Histo1d::setHist(int i){
  exit(1);
}

double Histo1d::getBinsize() const {
  return(binsize);
}

double Histo1d::getLeftBound(unsigned long bin) const {
  assert(bin < bincount);
  return( xmin + bin*binsize ); 
}

double Histo1d::getRightBound(unsigned long bin) const {
  assert(bin < bincount);
  return( xmin + (bin+1)*binsize ); 
}

double Histo1d::getCenter(unsigned long bin) const {
    assert(bin < bincount);
    return( xmin + binsize/2.0 + bin*binsize ); 
}

int Histo1d::getBin(double x) {
  if ( ( x < xmin ) || (x > xmax) ) {
    return (-1);
  }
  return(static_cast<int>( trunc( ( x - xmin ) / binsize ) ) );
}


void Histo1d::fill(double x){
  fill(x,1.0);
}

void Histo1d::fill(double x, double wgt){
  
  int bin = getBin(x);

  if ( bin < 0 ) {
    return;
  }
  current_correlated_result[ bin ] += wgt;
  
}


void Histo1d::fill2(double x, double wgt){
  
  /*
   * fill2 adds an event wgt to the bin = getBin(x) and 
   * to all the bins to the left
   * 
   */
  int bin = getBin(x);

  if ( bin < 0 ) {
    // outside range, determine wether (x < xmin) or (x>xmax)
    if ( x < xmin ) {
      return;
    } else {
      // x > xmax, add the event to all the bins...
      bin = bincount - 1 ;
    }
  }

  for (int i = 0; i <= bin; i++) {
    current_result[ i ] += wgt;
    current_result2[ i ] += wgt*wgt;
  }
}


void Histo1d::rescale(double x){
  for(int i=0; i < bincount; i++){
    current_result[i] *= x;  
    current_result2[i] *= x*x;
  }
}

void Histo1d::end_correlated_calls(){
  /*
    Ends a sequence of correlated calls, correctly filling the  current_results2 histogram. Then it void the temp histogram that contained the correlated calls.
   */ 
  for(int i=0; i < bincount; i++){
    current_result[i] += current_correlated_result[i];  
    current_result2[i] += current_correlated_result[i]* current_correlated_result[i];
    current_correlated_result[i] = 0.;
  }
}

void Histo1d::finalize(double calls){
  /* 
   * current_result2 is not yet the expectation value f^2/rho^2, 
   * we have to correct for an additional factor 1/N in front
   */
  cout << "# HISTSTART" << endl;
  cout << "# " << name << endl;
  double sum=0, sumacc=0, value,error;
  double current_error2_i;

  for(int i = 0; i < bincount; i++){
    if ( current_result2[i] != 0 ) {
      current_result2[i] *= calls;
      double tmp = sqrt( current_result2[i]);
      current_error2_i = ( tmp - current_result[i] ) 
	* ( tmp + current_result[i] ) / ( calls - 1.) ;
      double current_relative_inverse_error2 = 
	current_result[i]*current_result[i]/current_error2_i;
      accumulated_relative_inverse_error2[i] 
	+= current_relative_inverse_error2;
      accumulated_result[i] += 
	current_relative_inverse_error2 * current_result[i]; 
      if ( accumulated_relative_inverse_error2[i] != 0. ) {
	value = accumulated_result[i] / accumulated_relative_inverse_error2[i]
	  / binsize;
	error = abs(value) / sqrt(accumulated_relative_inverse_error2[i]); 
      } else {
	value = error = 0.;
      }
    } else {
      current_error2_i = 0.;
      accumulated_result[i] = 0.;
      accumulated_relative_inverse_error2[i] = 0.;
    }
    cout << getCenter(i) << " "
	 << current_result[i] / binsize << " " 
	 << sqrt(current_error2_i) / binsize << " " 
	 << value
	 << " " 
	 << error
	 << endl;
    sum += current_result[i];
    sumacc+= value * binsize;

    current_result[i] = current_result2[i] = 0;
  }
  cout << "# Sum this iteration: " << sum << endl;
  cout << "# Accumulated sum   : " << sumacc << endl;
}


void Histo1d::print() {

  
  cout << "# Histogram " << name << endl;
  cout << "# HISTSTART" << endl;
  cout << "# " << name << endl;
  cout << "# DATASTART: centre-of-bin, data, errorestimate " << endl;

  double sum=0,value,error;
  for(int i = 0; i<bincount; i++){
    if ( accumulated_relative_inverse_error2[i] != 0. ) {
      value = accumulated_result[i]/accumulated_relative_inverse_error2[i]
	/ binsize;
      error = abs(value)/ sqrt(accumulated_relative_inverse_error2[i]);
    } else {
      value=error=0.;
    }
    cout << getCenter(i)
	 << " "  
	 << value
	 << " "  
	 << error
	 << endl;
    sum+=value;
      
  }
    
  cout << "# DATAEND:" << endl;
  cout << "# sum =  " << sum*binsize << endl;
  cout << "# HISTEND" << endl;
  cout << endl << endl;


}

void Histo1d::printraw() {
  
  cout << "# Histogram " << name << endl;
  cout << "# HISTSTART" << endl;
  cout << "# " << name << endl;
  cout << "# DATASTART: centre-of-bin, data " << endl;
  for(int i = 0; i<bincount; i++){
    cout << getCenter(i) << " " << current_result[i ] << endl;
  }
  
}

void Histo1d::printfile(ofstream& out) {
  
  out << "(# Histogram " << name << endl;
  out << "(# HISTSTART" << endl;
  out << "(# " << name << endl;
  out << "(# DATASTART: centre-of-bin, data, errorestimate " << endl;
  static  bool firstprint=true;
  if(firstprint) {
    out << " SET BAR Y 0.02 PERMANENT" << endl;
    firstprint=false;
  }
  out << " ("<< name << endl;
  out << " TITLE BOTTOM \" "<<name<<" \""<< endl;
  out << " SET SCALE Y LOG" << endl;   
  out << " SET ORDER X Y DUMMY" << endl;   

  double sum=0,value,error;
  for(int i = 0; i<bincount; i++){
    if ( accumulated_relative_inverse_error2[i] != 0. ) {
      value = accumulated_result[i]/accumulated_relative_inverse_error2[i]
	/ binsize;
      error = abs(value)/ sqrt(accumulated_relative_inverse_error2[i]);
    } else {
      value=error=0.;
    }
    out << getCenter(i)
	 << " "  
	 << value
	 << " "  
	 << error
	 << endl;
    sum+=value;
      
  }
  out << " HIST SOLID" << endl;
  out << " PLOT" << endl;
  out << "(# DATAEND:" << endl;
  out << "(# sum =  " << sum*binsize << endl;
  out << "(# HISTEND" << endl;
  out << endl << endl;


}



//#define TEST
#ifdef TEST

double gauss(){
  double x,test;
  do {
    x = 20.0*drand48()-10.0;
    test = drand48();
  } while ( test > exp( -x*x/10.0 ) );
  return(x);
}
 

int main(){

  Histo1d a("test",0.,1.,100);
  double x;
 
  for (int i = 0; i< 10000; i++){
    x=drand48();
    a.fill(x,x+100.);
    a.fill(x,-100.);
    a.fill(x,20.);
    a.fill(x,-20.);
  };
  
  a.end_correlated_calls();

  a.finalize(10000);
  a.print();


  for (int i = 0; i< 10000; i++){
    x=drand48();
    a.fill(x,x+100.);
    a.fill(x,-100.);
    a.fill(x,20.);
    a.fill(x,-20.);
  };
  
  a.end_correlated_calls();

  a.finalize(10000);
  a.print();


  
  cout << "# Final result ..." << endl;
  cout << "# DATASTART: " << endl;
  a.print();
  cout << "# DATAEND:" << endl;

  return(0);
}

#endif
