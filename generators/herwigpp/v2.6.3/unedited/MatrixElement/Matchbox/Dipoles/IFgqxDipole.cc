// -*- C++ -*-
//
// IFgqxDipole.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgqxDipole class.
//

#include "IFgqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFgqxDipole::IFgqxDipole() 
  : SubtractionDipole() {}

IFgqxDipole::~IFgqxDipole() {}

IBPtr IFgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 6 &&
    partons[emission]->mass() == ZERO &&
    partons[spectator]->mass() == ZERO;
}

double IFgqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaS()/
	underlyingBornME()->lastXComb().lastAlphaS(),
	underlyingBornME()->orderInAlphaS());

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaEM()/
	underlyingBornME()->lastXComb().lastAlphaEM(),
	underlyingBornME()->orderInAlphaEW());

  lastME2(res);

  return res;

}

double IFgqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (realEmissionME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaS()/
	underlyingBornME()->lastXComb().lastAlphaS(),
	underlyingBornME()->orderInAlphaS());

  res *=
    pow(realEmissionME()->lastXComb().lastAlphaEM()/
	underlyingBornME()->lastXComb().lastAlphaEM(),
	underlyingBornME()->orderInAlphaEW());

  lastME2(res);

  logME2();

  return res;

}

void IFgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFgqxDipole::Init() {

  static ClassDocumentation<IFgqxDipole> documentation
    ("IFgqxDipole");

  DipoleRepository::registerDipole<IFgqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFgqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFgqxDipole,SubtractionDipole>
describeHerwigIFgqxDipole("Herwig::IFgqxDipole", "HwMatchbox.so");