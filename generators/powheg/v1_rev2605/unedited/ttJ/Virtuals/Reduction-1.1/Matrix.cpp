// $Modified: Mon Sep 18 11:35:50 2006 by puwer $
#include <iostream>
#include "Matrix.h"

using namespace std;


void Matrix::printMaple() const {
  cout.precision(15);
  cout << "S:=Matrix([" ;
  for (int i = 1; i<=nrows; i++){
    cout << "[" ;
    for (int j = 1; j<=nrows; j++){
      cout << store[locate(i,j)];
      if (j!=ncols)
	cout <<",";
    }
    cout << "]";
    if (i!=nrows)
	cout <<",";
    cout << endl;
  }
  cout << "]);" << endl;
}

void Matrix::printC() const {
  cout.precision(15);
  for (int i = 1; i<=nrows; i++){
    for (int j = 1; j<=nrows; j++){
      cout << "S.setValue(" << i << "," << j 
	   << store[locate(i,j)] << ");";
    }
  }
}


const Matrix operator-(const Matrix & M){
  Matrix tmp = Matrix(M);
  tmp.negate();
  return(tmp);
}

const Matrix & operator+(const Matrix & M){
  return(M);
}

const Matrix operator+(const Matrix & M1, const Matrix & M2){
  Matrix tmp = Matrix(M1);
  tmp.add(M2);
  return tmp;
}

const Matrix operator-(const Matrix & M1, const Matrix & M2){
  Matrix tmp = Matrix(M1);
  tmp.sub(M2);
  return tmp;
}

const Matrix operator*(const double & s, const Matrix & M){
  Matrix tmp = Matrix(M);
  tmp.mul(s);
  return tmp;
}

const Matrix operator*(const Matrix & M, const double & s ){
  Matrix tmp = Matrix(M);
  tmp.mul(s);
  return tmp;
}

const Matrix operator*(const Matrix & M1, const Matrix & M2){
  return M1.mul(M2);
}


const Matrix operator/(const Matrix & M, const double & s ){
  Matrix tmp = Matrix(M);
  tmp.div(s);
  return tmp;
}


std::ostream& operator<<(std::ostream& os, const Matrix & M){
  for (int m=1; m <= M.nrows; m++){
    os << "[ ";
  for (int n=1; n <= M.ncols; n++) { 
      os << M.store[M.locate(m,n)] << " "; 
    }
    os << "]\n";
  }
  return os;
}
