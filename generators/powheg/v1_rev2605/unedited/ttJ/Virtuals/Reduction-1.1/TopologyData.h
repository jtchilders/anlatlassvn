/* $Modified: Mon Sep 18 11:35:11 2006 by puwer $ */
#ifndef TOPOLOGYDATA_H_
#define TOPOLOGYDATA_H_

class Topology;
class data;

#include "Matrix.h"
#include "zlib.h"
#include "svdcmp.h"
#include <vector>
#include <map>


class data {

 private:
  Topology* parent;
  static std::map<uLong,data*> topologies;
  uLong Adler32_id;    
  uLong CRC32_id; 
  int npoint;

  data(const Matrix& Sin);

 public:  
  int imin;
  double B, detS;

  Matrix *S,*U,*V,*W,*Sinv;
  double *z;
  double *range;

  static data* lookup(const Matrix& Sin);

  void setParent(Topology * parent_tmp){
    parent = parent_tmp;
  }

  Topology* getParent(){
    return(parent);
  }

  static void clear();

  ~data(){
    //std::cout << "data deconstructor called "<<npoint << "-point" << std::endl; 
    delete S;
    delete U;
    delete V;
    delete W;
    delete Sinv;
    //    std::cout << " fix the Topology deconstructor \n"; 
    if (z!=0) delete[] z;
    if (range!=0) delete[] range;
    // std::cout << "done" << std::endl; 
  }
  
  int getNpoint(){ return(npoint); }
  void setAdler(uLong tmp){ Adler32_id = tmp; };
  void setCRC(uLong tmp){ CRC32_id = tmp; };

};

#endif //TOPOLOGYDATA_H_
