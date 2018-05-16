// $Modified: Mon Sep 18 11:34:09 2006 by puwer $
#include "TopologyData.h"
#include "GGReduction.h"
#include <stdexcept>

using namespace std;

void data::clear(){

  for ( std::map<uLong,data*>::iterator i=topologies.begin(); 
	i!=topologies.end(); i++ ) {
    data *tmp = (*i).second;
    delete (*tmp).parent;
    delete tmp;
  }
  topologies.clear();
}

data* data::lookup(const Matrix & Sin){
   
  /*
   * Check if the topology is already known:
   */
  int size = Sin.getNRows();
  double Stmp[size][size];
  
  for (int i = 1; i<=size; i++){
    for (int j = 1; j<=size; j++){
      Stmp[i-1][j-1] = Sin(i,j); 
    }
  }
  
  uLong Adler32_id_tmp; 
  uLong CRC32_id_tmp; 
    
  Adler32_id_tmp = adler32(0L, Z_NULL, 0);
  Adler32_id_tmp = adler32(Adler32_id_tmp,(Bytef*)Stmp,sizeof(Stmp));
    
  CRC32_id_tmp = crc32(0L, Z_NULL, 0);
  CRC32_id_tmp = crc32(CRC32_id_tmp,(Bytef*)Stmp,sizeof(Stmp));

  map<uLong,data*>::iterator position = topologies.find(CRC32_id_tmp);
  if ( position != topologies.end() ) {
    if (Adler32_id_tmp != (position->second)->Adler32_id ){
#ifdef WITHEXCEPTIONS
      throw(logic_error("GGreduction: CRC32 checksum not uniq drop event\n"));
#else
      std::cout << "GGreduction: CRC32 checksum not uniq drop event\n";
      exit(1);
#endif
    }
    /*
     * Topology is already known, retrieve it from cache
     */
    return(position->second);
  } else {
    /*
     * Not yet known calculate it...
     */
    
    data* tmp = new data(Sin);
    // Register in the cache:

    topologies[CRC32_id_tmp] = tmp;
    tmp->setAdler(Adler32_id_tmp);
    tmp->setCRC(CRC32_id_tmp);
    
    return(tmp);
  }
}


data::data(const Matrix & Sin){
  parent = 0;
  npoint = Sin.getNColumns();
  
  z = 0;
  range = 0;
  S = new Matrix(npoint,npoint);
  
  *S = Sin;
  /*    for (int i = 1; i<=npoint; i++){
	S->setValue(i,i,-2.0*ms[i-1]*ms[i-1]);
	for (int j = i+1; j<=npoint; j++){
	S->setValue(i,j, dotp(qs[i-1]-qs[j-1],qs[i-1]-qs[j-1])
	-ms[i-1]*ms[i-1]-ms[j-1]*ms[j-1]);
	S->setValue(j,i, (*S)(i,j) );
	}
	}*/
  
  S->clean();
  
  if (npoint > 2) {
    /***
     *** The SVD is only need for npoint > 2 
     ***/
    U = new Matrix(npoint,npoint);
    V = new Matrix(npoint,npoint);
    W = new Matrix(npoint,npoint);
    Sinv = new Matrix(npoint,npoint);
    
    /***
     *** Singular Value Decomposition:
     ***/
    svdcmp(*S,*U,*V,*W);
    
    
    /***
     *** Try to figure out if the matrix is singular
     ***/
    {
      double min = fabs((*W)(npoint,npoint)) , max = min;
      imin = npoint;
      for (int i=1;i<npoint;i++) {
	if (fabs((*W)(i,i)) > max) 
	  max = fabs((*W)(i,i));
	if (fabs((*W)(i,i)) < min) {
	  min = fabs((*W)(i,i));
	  imin = i;
	}
      }
      
      if (min/max < DBL_EPSILON) {
	detS = 0.0;
	z = new double[npoint];
	range = new double[npoint];
	for (int i = 0; i<npoint; i++){
	  /*** 
	   *** Construct the direction of the null space:
	   ***/
	  z[i] = (*V)(i+1,imin) ;
	  /***
	   *** Construct the direction othogonal to the plane which can be
	   *** reached by S. We assume here that only one eigenvalue is zero.
	   *** This guaranteed by the structure of the S matrix.
	   ***/
	  range[i] = (*U)(i+1,imin) ;	
	}
	  
	{
	  double rnorm = 0.0;
	  
	  for (int i = 0; i<npoint; i++) 
	    rnorm += range[i]*range[i];
	  
	  for (int i = 0; i<npoint; i++) 
	    range[i] = range[i] / sqrt(rnorm);
	}

	/*** 
	 *** There is no inverse because det(S) = 0, 
	 *** the equation S*y = b can still have a solution,
	 *** if b is in the range of S. This solution is not uniq 
	 *** but this is not relevant for our purpose here. If 
	 *** b is in the range we get y by 
	 *** y = (V*W2*U^tr) * b 
	 ***/
	(*Sinv) = (*V)*(*W).diaginv0(imin)*(*U).transpose();
	
      } else {
	detS = 1.0;
	for (int i=1;i<=npoint;i++) {
	  detS *= (*W)(i,i); 
	}
	(*Sinv) = (*V)*(*W).diaginv()*(*U).transpose();
	B = 0.0;
	double b[npoint+1];
	for (int i=1;i<=npoint;i++) {
	  b[i] = 0.0;
	  for (int j=1; j<=npoint;j++){
	    b[i] += (*Sinv)(i,j);
	  }
	  B += b[i];
	}
      }
    } 
  } else {
    U = 0;
    V = 0;
    W = 0;
    Sinv = 0;
  }
}
