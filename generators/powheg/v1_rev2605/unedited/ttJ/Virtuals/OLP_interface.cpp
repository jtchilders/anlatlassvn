#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <sstream>
#include "ttyctrl.h"
#include <fstream>
#include <vector>
#include <string>
#include "FourMomentum.h"
#include "FiniteVirtuals.h"

using namespace std;


void OLP_Start(const string, int &); 
void VirtualInitialize(); 
void SetVirtualScales(int verbose = 0); 

void OLP_EvalSubProcess(const int &, const double *, const double *, const double &, double *);

extern "C" {
  void olp_start_(int & ios, char * filename, int len) {
    OLP_Start(filename,ios); 
  }
  
  void olp_evalsubprocess_(const int & code, const double * p, const double * scales, const double & as, double * virt) {
    OLP_EvalSubProcess(code,p,scales,as,virt);
  }
}


void OLP_Start(const string filename, int & ios) {
  string line,orderline,contractline;
  ifstream in;
  ofstream out;
  vector <string> fileline;
  static bool firstcall=true;
  if(firstcall) { 
    cout<<"  Reading "<<filename.c_str()<<endl;
  } else {
    cout<<"  Virtual OLP initialization "<<endl;
    cout<<endl;
  }
  // First look if the contract file has already been checked and accepted
  in.open(filename.c_str());
  if(in.is_open()) {
    while( !in.eof() ) {
      getline(in,line);
      if (line.find("# contract file") != string::npos) {
	if(firstcall) {
	  cout<<"  Contract file "<<filename.c_str()<<" already accepted"<<endl;
	  ios=1;
	  firstcall=false;
	  return;
	} else {
	  // the contract file has already been accepted and the
	  // OLP codes of subprocesses have been mapped
	  // Now it is possible to initialize virtual routines
	  ios=-1;
	  firstcall=false;
	  VirtualInitialize();
	  return;
	} 
      }
      if (line.find("Error") != string::npos) {
	cout<<"  ERROR: contract file "<<filename.c_str()<<" already checked and discarded"<<endl;
	cout<<"  please remove it and provide a different one"<<endl;
	ios=0;
	return;
      }
    }
  } else {
    cout<<"Error: cannot open "<<filename<<endl;
    ios=0;
    return;
  }
  in.close();
  // Otherwise parse it and add informations
  in.open(filename.c_str());
  if(in.is_open()) {
    while( !in.eof() ) {
      getline(in,line);
      if (line.find("# order file") != string::npos) {
	ios=1;
	orderline=line;
	contractline="# contract file produced by libvirtual";
      }
      if (line.find("MatrixElementSquareType") != string::npos) {
	if (line.find("CHSummed") != string::npos)    {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# CHSummed is supported";
	}
      }
      if (line.find("CorrectionType")!= string::npos) {
	if (line.find("QCD") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# QCD is supported ";
	}
      }
      if (line.find("IRregularisation") != string::npos) {
	if(line.find("CDR") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# CDR is supported";
	}
      }
      if (line.find("MassiveParticleScheme") != string::npos) {
	if(line.find("OnShell") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# OnShell is supported";
	}
      }
      if (line.find("OperationMode") != string::npos) {
	if(line.find("CouplingsStrippedOff") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# CouplingStrippedOff is supported";
	}
      }
      // if (line.find("#ModelFile                ") != string::npos) 
      if (line.find("SubdivideSubprocess") != string::npos) {
	if(line.find("  no") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported flag  \n# no is supported";
	}
      }
      if (line.find("AlphasPower") != string::npos) {
	if(line.find("3") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported value  \n# 3 is supported";
	}
      }
      if (line.find("AlphaPower") != string::npos) {
	if(line.find("0") != string::npos) {  
	  line+=" | OK";
	} else {
	  ios*=0;
	  line+=" | Error: unsupported value  \n# 0 is supported";
	}
      }
      if(line.find("->")!=string::npos ) {
	int flav1,flav2;
	stringstream(line) >> flav1 >> flav2;
	line+= "        | ";
	if(flav1==21 && flav2==21) { 
	  line+="OK  1  # g g channel";
	} else if(flav1 < 0 && flav2 == 21) {
	  line+="OK  2  # qbar g channel";
	} else if(flav1 > 0 && flav2 ==21) {
	  line+="OK  3  # q g channel";
	} else if(flav1==21 && flav2 < 0) {
	  line+="OK  4  # g qbar channel";
	} else if(flav1==21 && flav2 > 0) {
	  line+="OK  5  # g q channel";
	} else if(flav1 < 0 && flav2 > 0) {
	  line+="OK  6  # qbar q channel";
	} else if(flav1 > 0 && flav2 < 0) {
	  line+="OK  7  # q qbar channel";
	} else {
	  ios*=0;
	}
      }
      fileline.push_back(line);
    }
  } else {
    cout<<"Error: cannot open "<<filename<<endl;
    ios=0;
    return;
  }
  in.close();
  out.open(filename.c_str());
  for (int iline=0; iline<fileline.size(); iline++){
    if((fileline[iline].find(orderline) != string::npos) && (ios==1)){
      out<<contractline<<endl;
    } else {
      out<<fileline[iline]<<endl;
    }
  }
  out.close();
  firstcall=false;
  if(ios==1)  cout<<"  Contract file checked and accepted"<<endl;
  return;
}


void OLP_EvalSubProcess(const int & code, const double * p, const double * scales, const double & as, double * virt) 
{
  FourMomentum xp1("p1",p[0],p[1],p[2],p[3]);
  FourMomentum xp2("p2",p[5],p[6],p[7],p[8]);
  FourMomentum xkq("kq",p[10],p[11],p[12],p[13]);
  FourMomentum xkqb("kqb",p[15],p[16],p[17],p[18]);
  FourMomentum xp3("p3",p[20],p[21],p[22],p[23]);

  SetVirtualScales();
  
  switch(code) {
  case 1:  // g g
    virt[2]=gggtt2virtfin(xp1,xp2,xp3,xkq,xkqb)/256.;
    break;

  case 2:  // qbar g
    virt[2]=gqbttqb2virtfin(xp2,xp1,xp3,xkq,xkqb)/96.;
    break;

  case 3:  // q g
    virt[2]=qgttq2virtfin(xp1,xp2,xp3,xkq,xkqb)/96.;
    break;

  case 4:  // g qbar
    virt[2]=gqbttqb2virtfin(xp1,xp2,xp3,xkq,xkqb)/96.;
    break;

  case 5:  // g q
    virt[2]=qgttq2virtfin(xp2,xp1,xp3,xkq,xkqb)/96.;
    break;

  case 6:  // qbar q
    virt[2]=qqttg2virtfin(xp2,xp1,xp3,xkq,xkqb)/36.;
    break;

  case 7:  // q qbar
    virt[2]=qqttg2virtfin(xp1,xp2,xp3,xkq,xkqb)/36.;
    break;
      
  default:
    cout<<"virtual OLP code not allowed"<<endl;
    exit(1);
  }

}
