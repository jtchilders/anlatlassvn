/* $Modified: Tue Dec 16 17:18:27 2008 by uwer $ */
#ifndef DIPOLES_H_
#define DIPOLES_H_

#include <vector>
#include <complex>
#include "FourMomentum.h"
#include "StandardModelParameters.h"

#define MAXPARTICLECOUNT 7

inline int Key2Pol(unsigned int key, int parton){
  return( key & (1<<parton) );
}

inline unsigned int PolKey(const int& pol0, const int& pol1, const int& pol2, 
		    const int& pol3, const int& pol4){
  return( pol0 + (pol1 << 1) + (pol2 << 2) + (pol3 << 3) + (pol4 << 4) );
}

inline unsigned int PolKeyFlip(unsigned int key, int parton){
  return( key ^ (1<<parton) );
}

void printBinary(const unsigned int val);

unsigned int SubprocessKey(std::vector<Particles> & p);


class Dipole {
 public:
  double value;
  double Vdiag;
  /* Alpha value a la Zoltan Nagy, only for massless dipoles: */
  double alpha; 
  std::complex<double> Vpm;
  int empty,i,j,k;
  int emitter;
  int spectator;
  Particles emittertype;
  std::vector<FourMomentum> momenta;
  static std::vector<FourMomentum> auxmomori;
  std::vector<FourMomentum> auxmomcpy;
  std::vector<FourMomentum> * auxmom;
  std::vector<Particles> ptypes; 

  Dipole(int i_, int j_, int k_) : 
  empty(0),auxmom(&auxmomori),i(i_),j(j_),k(k_) {};
  Dipole() : empty(1),auxmom(&auxmomori) {};
  void setaux(std::vector<FourMomentum> & momenta){auxmomori = momenta;};
  void PrintParticles();
  void PrintMomenta();

  friend std::ostream& operator<<(std::ostream& os, const Dipole & d);

};

class SplittingKernels {

 public:
  SplittingKernels(std::vector<FourMomentum> momenta,
		   std::vector<Particles> particles);
  void Kernel( Dipole & );
  void FinalFinal(Dipole &);
  void FinalInitial(Dipole &);
  void InitialFinal(Dipole &);
  void InitialInitial(Dipole &);

  inline bool isquark(int i){
    return StandardModelParameters::isquark(ptypes[i]);
  }

  inline bool isgluon(int i){
    return StandardModelParameters::isgluon(ptypes[i]);
  }

 private:  
  static const double kappa;
  StandardModelParameters& parms;

  double dot[MAXPARTICLECOUNT][MAXPARTICLECOUNT];
  std::vector<FourMomentum> momset;
  std::vector<Particles> ptypes;

};

Dipole Swap(Dipole d);

template<int N> class ColorCorrelations {
public:
  double c[N+1][N+1];
  void set(const ColorCorrelations<N>& src ){
    for ( int i = 1; i <= N; i++){
      for ( int j = 1; j <= N; j++){
	c[i][j] = src.c[i][j];
      }
    }
  }
};

class Process {

 public:
  std::vector<Particles> particles;

  unsigned int key()  {
    return(SubprocessKey(particles));
  }

  void printProcessName(){
    for (std::vector<Particles>::iterator i = particles.begin(); 
	 i!= particles.end(); i++) {
      std::cout << *i ;
    }
  }

  virtual void EvalfDipole(Dipole & d) = 0;  
  virtual double EvalfAmplitudeSquared(const std::vector<FourMomentum> & pset)
    = 0;
};


class Correlator {

private:
  double (*CorrObs)(const std::vector<FourMomentum>& pi);
  std::map<unsigned int, Process* > processes;
  
public:

  Correlator(double (*obs_)(const std::vector<FourMomentum>& pi)) 
    : CorrObs(obs_) {};

  void EvalfDipole(Dipole & d);

  void add(Process* process) {
    processes[process->key()] = process;
    //    process->printProcessName();
  }

};

#endif //DIPOLES_H_ 
