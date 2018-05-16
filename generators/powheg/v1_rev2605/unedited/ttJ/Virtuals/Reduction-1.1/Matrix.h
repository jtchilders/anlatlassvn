/* $Modified: Mon Oct  9 13:54:47 2006 by puwer $ */
#ifndef _MATRIX_H_
#define _MATRIX_H_

#define NDEBUG

#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <cfloat>

class Matrix {
 private:
  int ncols,nrows;
  double* store;

  int locate(int m, int n) const {

    assert(m > 0 && m <= nrows && n > 0 && n <= ncols);

    return((m-1)*ncols+n-1);

  }

 public:
  Matrix(const Matrix & M) {
    ncols = M.ncols;
    nrows = M.nrows;
    store = new double[ncols*nrows];
    memcpy(store,M.store,nrows*ncols*sizeof(double));
  };

  //Matrix(int n){ Matrix(n,n); }; 

  Matrix(int rows, int cols){
    ncols = cols;
    nrows = rows;
    store = new double[ncols*nrows];
  };

  ~Matrix(){
    delete[] store;
  }
  
  Matrix& operator=(const Matrix & M){
    if (&M != this){
      assert( (nrows==M.nrows) && (ncols == M.ncols) ) ;
      memcpy(store,M.store,nrows*ncols*sizeof(double));
    } 
    return *this;
  }

  void clean(){
    /*** 
     *** Function to identify possible zeros.
     ***/
    double min=DBL_MAX, max=0.0;
    // Determine minimal and maximal value:
     for (int m=1; m <= nrows; m++)
       for (int n=1; n <= ncols; n++){
	 if (fabs(store[locate(m,n)]) < min)
	   min = fabs(store[locate(m,n)]);
	 if (fabs(store[locate(m,n)]) > max)
	   max = fabs(store[locate(m,n)]);
       }
     /***
      *** If the ratio between an entry and the max value is smaller than
      *** DBL_EPSILON (the difference between 1 and the minimum double 
      *** greater than 1) we consider it the entry as zero.
      ***/
     for (int m=1; m <= nrows; m++)
       for (int n=1; n <= ncols; n++){
	 if (fabs(store[locate(m,n)]) < max * 10.0 * DBL_EPSILON)
	   store[locate(m,n)] = 0.0;
       }
  }

  double operator()(int m, int n) const { 
   return store[locate(m,n)];
  }

  void add(const Matrix & M){
    
    assert((nrows == M.nrows) && (ncols == M.ncols));

    for (int m=1; m <= nrows; m++)
      for (int n=1; n <= ncols; n++)
	store[locate(m,n)] += M.store[M.locate(m,n)] ;

  }

  void sub(const Matrix & M){

    assert((nrows == M.nrows) && (ncols == M.ncols));

    for (int m=1; m <= nrows; m++)
      for (int n=1; n <= ncols; n++)
	store[locate(m,n)] -= M.store[M.locate(m,n)] ;
  }

  void mul(const double & s) {
     for (int m=1; m <= nrows; m++)
	for (int n=1; n <= ncols; n++)
	  store[locate(m,n)] *= s;
  }

  void div(const double & s) {
     for (int m=1; m <= nrows; m++)
	for (int n=1; n <= ncols; n++)
	  store[locate(m,n)] /= s;
  }


  Matrix getSubMatrix(int i) const {

    assert(( i <= nrows ) && (i <= ncols ));

    Matrix tmp = Matrix(nrows-1,ncols-1);
    for (int m=1; m < nrows; m++)
      for (int n=1; n < ncols; n++) {
	tmp.store[tmp.locate(m,n)] = 
	  store[locate(m < i ? m : m + 1, n < i ? n : n + 1)];
      }
    return(tmp);    
  }

  const Matrix transpose() {
    Matrix tmp = Matrix(nrows,ncols);
    for (int m=1; m <= nrows; m++)
      for (int n=1; n <= ncols; n++)
	tmp.store[locate(n,m)] = store[locate(m,n)];
    return tmp;
  }

  const Matrix diaginv() {
    Matrix tmp = Matrix(*this);
    for (int m=1; m <= nrows; m++)
      tmp.store[locate(m,m)] = 1.0/store[locate(m,m)];
    return tmp;
  }
  
  const Matrix diaginv0(int i) {
    Matrix tmp = Matrix(*this);
    for (int m=1; m <= nrows; m++) { 
      if (m!=i) {
	tmp.store[locate(m,m)] = 1.0/store[locate(m,m)];
      } else {
	tmp.store[locate(m,m)] = 0.0;
      }
    }
    return tmp;
  }

  const Matrix mul(const Matrix & M) const {

    assert(nrows == M.ncols);

    double t,y,c,sum;
    Matrix tmp = Matrix(ncols,M.nrows);
    for (int m=1; m <= ncols; m++)
      for (int n=1; n <= M.nrows; n++) {
	sum = 0.0;
	c = 0.0;
	for (int l=1; l <= nrows; l++){
	  y =  store[locate(m,l)]*M.store[locate(l,n)] - c;
	  t = sum + y;
	  c = (t - sum) - y;
	  sum = t;
	  //  tmp.store[locate(m,n)] += store[locate(m,l)]*M.store[locate(l,n)];
	}
	tmp.store[tmp.locate(m,n)] = sum;
      }
    return tmp;
  }

  void negate(){
    for (int m=1; m <= nrows; m++)
      for (int n=1; n <= ncols; n++)
	 store[locate(m,n)] *= -1.0 ; 
  }

  void setValue(int i, int j, double val){

    assert((i<=ncols) && (j<=nrows));
    store[locate(i,j)] = val;
  }


  void setValue(double (*func)(int, int)){
    for (int m=1; m <= nrows; m++)
      for (int n=1; n <= ncols; n++)
	 store[locate(m,n)] = func(m,n); 
  }

  

  int getNColumns() const { return ncols;}
  int getNRows() const { return nrows;}
 
  void printMaple() const ;
  void printC() const ;

  friend const Matrix operator-(const Matrix & M);
  friend const Matrix & operator+(const Matrix & M);
  friend const Matrix operator-(const Matrix & M1, const Matrix & M2);
  friend const Matrix operator+(const Matrix & M1, const Matrix & M2);
  friend const Matrix operator*(const Matrix & M1, const Matrix & M2);
  friend const Matrix operator*(const double & s, const Matrix & M);
  friend const Matrix operator*(const Matrix & M, const double & s );
  friend const Matrix operator/(const Matrix & M, const double & s );

  friend std::ostream& operator<<(std::ostream& os, const Matrix & m); 

};

#undef NDEBUG

#endif //  _MATRIX_H_
