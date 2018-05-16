// -*- C++ -*-
//
// TwoBodyAllOnCalculator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyAllOnCalculator class.
//

#include "TwoBodyAllOnCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "GenericWidthGenerator.h"

namespace Herwig {
using namespace ThePEG;

Energy TwoBodyAllOnCalculator::partialWidth(Energy2 scale) const {
  return _widthgen->partial2BodyWidth(_mode,sqrt(scale),_mass1,_mass2);
}

}
