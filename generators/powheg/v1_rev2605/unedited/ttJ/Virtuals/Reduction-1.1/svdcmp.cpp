/* $Modified: Thu Ja  5 16:47:43 20: Thu Oct  7 22:40:18 2010 by uwer 06 by puwer $ */
#include "Matrix.h"
#include "svdcmp.h"
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <string>


#define SQR(_a) (_a)*(_a)
#define SIGN(_a,_b) ( (_b) >= 0.0 ? fabs(_a) : -fabs(_a))
#define FMAX(_a,_b) ( (_a) > (_b) ? (_a) : (_b) )
#define IMIN(_a,_b) ( (_a) < (_b) ? (_a) : (_b) )

#define NR_END 1
#define FREE_ARG char*

using namespace std;

void nrerror(string error_text)
/* Numerical Recipes standard error handler */
{
  cout << "Numerical Recipes run-time error..." << endl;
  cout << error_text << endl;
  cout << "...now exiting to system..." << endl;
  exit(-1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}



double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmp(const Matrix & A, Matrix & U, Matrix & V, Matrix & W)
{
  //void svdcmp(double **a, int m, int n, double w[], double **v)
  /* 
   * Given a matrix a[1..m][1..n], this routine computes its singular 
   * value decomposition, A = U·W·V^T. The matrix U replaces A on output. 
   * The diagonal matrix of singular values W is output as a vector w[1..n]. 
   * The matrix V (not the transpose V T ) is output as v[1..n][1..n].
   */

  bool printit = false;


  const int m = A.getNColumns();
  const int n = A.getNRows();

  if ((U.getNColumns() != m) || ( U.getNRows()!= n )) {
    cout << "svdcmp: Matrix U not of correct size" << endl;
    exit(-1);
  }

  if ((V.getNColumns() != n) || ( V.getNRows()!= n )) {
    cout << "svdcmp: Matrix V not of correct size" << endl;
    exit(1);
  }

  if ((W.getNColumns() != n) || ( W.getNRows()!= n )) {
    cout << "svdcmp: Matrix W not of correct size" << endl;
    exit(1);
  }

  double a[m+1][n+1];

  for(int jjj=1; jjj <=m ; jjj++)
    for(int kkk=1; kkk <=n ; kkk++){
      a[jjj][kkk] = A(jjj,kkk);
    }

  //  double w[n+1];     
  //double v[n+1][n+1];
  double *w, **v;
  w = dvector(1,n);
  v = dmatrix(1,n,1,n);

  
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  
  double *rv1;
  rv1 = dvector(1,n);
  // double rv1[n+1]; 


  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    //for (its=1;its<=30;its++) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	//if ((double)(fabs(rv1[l])+anorm) == anorm) {
	if ( fabs( rv1[l]/anorm) < DBL_EPSILON){
	  flag=0;
	  break;
	}
      	//if ((double)(fabs(w[nm])+anorm) == anorm) break;
	if ( fabs( w[nm]/anorm) < DBL_EPSILON ) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  //if ((double)(fabs(f)+anorm) == anorm) break;
	  if ( fabs( f/anorm)< DBL_EPSILON ) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (its == 30) {
	printit = true;
	cout << "svdcmp:\n";
	A.printC();
      //	nrerror("no convergence in 30 svdcmp iterations");
      }
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
	}

  for(int jjj=1; jjj <=m ; jjj++)
    for(int kkk=1; kkk <=n ; kkk++)
      U.setValue(jjj,kkk,a[jjj][kkk]);

  for(int jjj=1; jjj <=n ; jjj++)
    for(int kkk=1; kkk <=n ; kkk++)
      V.setValue(jjj,kkk,v[jjj][kkk]);

  for(int jjj=1; jjj <=n ; jjj++)
    for(int kkk=1; kkk <=n ; kkk++)
      (jjj == kkk) ? W.setValue(jjj,jjj,w[jjj]) : W.setValue(jjj,kkk,0.0);


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if (printit) {
    cout << W << endl;
    cout << U << endl;
    cout << V << endl;
  }
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  free_dvector(rv1,1,n);
  free_dvector(w,1,n);
  free_dmatrix(v,1,n,1,n);
}

/* (C) Copr. 1986-92 Numerical Recipes Software Z(2s5(,&2021:N%#(.)+. */
