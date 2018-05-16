/* $Modified: Fri Feb  2 13:35:52 2007 by puwer $ */
#ifndef PROCESSES_H_
#define PROCESSES_H_

#include "FourMomentum.h"
#include "Dipoles.h"
#include "Color.h"
#include <vector>
#include <complex>

#define SYMMETRIZE(i_,j_) cij[j_][i_].set(cij[i_][j_])    

class ggttg : public Process  {

private:
  ColorCorrelations<6>  cij[5][5];
  std::complex<double>  ampLO[7][32];

public:
  ggttg();
  const ColorCorrelations<6> (*getColorMatrices())[5][5] {
    return(&cij);
  }
  inline double Color(int i, int j, int c1, int c2) const {
    return(cij[i][j].c[c1][c2]);
  }
  void EvalfAmplitude(const std::vector<FourMomentum> & pset);
  void EvalfAmplitude2(const std::vector<FourMomentum> & pset);
  void EvalfDipole(Dipole & d);
  double EvalfAmplitudeSquared(const std::vector<FourMomentum> & pset);
  std::complex<double> EvalfCorrelation(double KP[7][7]);
};


class qqttg : public Process  {

private:
  ColorCorrelations<4> cij[5][5];
  std::complex<double> ampLO[5][32];

public:
  qqttg();
  const ColorCorrelations<4> (*getColorMatrices())[5][5] {
    return(&cij);
  }
  inline double Color(int i, int j, int c1, int c2) const {
    return(cij[i][j].c[c1][c2]);
  }
  void EvalfAmplitude(const std::vector<FourMomentum> & pset);
  void EvalfDipole(Dipole & d);
  double EvalfAmplitudeSquared(const std::vector<FourMomentum> & pset);
  std::complex<double> EvalfCorrelation(double KP[5][5]);
}; // qqttg


class qgttq : public Process  {

private:
  ColorCorrelations<4> cij[5][5];
  std::complex<double> ampLO[5][32];

public:
  qgttq(); 
  const ColorCorrelations<4> (*getColorMatrices())[5][5] {
    return(&cij);
  }
  inline double Color(int i, int j, int c1, int c2) const {
    return(cij[i][j].c[c1][c2]);
  }
  void EvalfAmplitude(const std::vector<FourMomentum> & pset);
  void EvalfDipole(Dipole & d);
  double EvalfAmplitudeSquared(const std::vector<FourMomentum> & pset);
  std::complex<double> EvalfCorrelation(double KP[5][5]);
  void setFlavour(Particles type) {
    particles[0] = type;
    particles[4] = type;
  }

}; // qgttq



class gqbttqb : public Process  {

private:
  ColorCorrelations<4> cij[5][5];
  std::complex<double> ampLO[5][32];

public:
  gqbttqb();
  const ColorCorrelations<4> (*getColorMatrices())[5][5] {
    return(&cij);
  }
  inline double Color(int i, int j, int c1, int c2) const {
    return(cij[i][j].c[c1][c2]);
  }
  void EvalfAmplitude(const std::vector<FourMomentum> & pset);
  void EvalfDipole(Dipole & d);
  double EvalfAmplitudeSquared(const std::vector<FourMomentum> & pset);
  std::complex<double> EvalfCorrelation(double KP[5][5]);
  void setFlavour(Particles type) {
    particles[1] = type;
    particles[4] = type;
  }

}; // gqbttqb

#endif //PROCESSES_H_
