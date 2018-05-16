// -*- C++ -*-
//
// MatchboxMElP2lJetJet.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMElP2lJetJet_H
#define HERWIG_MatchboxMElP2lJetJet_H
//
// This is the declaration of the MatchboxMElP2lJetJet class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMEllbarqqbarg.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Luca D'Errico
 *
 * \brief MatchboxMElP2lJetJet implements the matrix element
 * for neutral current, charged lepton 1+1 jet DIS with light quarks.
 *
 * @see \ref MatchboxMElP2lJetJetInterfaces "The interfaces"
 * defined for MatchboxMElP2lJetJet.
 */
class MatchboxMElP2lJetJet: public MatchboxMEBase, public MatchboxMEllbarqqbarg {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMElP2lJetJet();

  /**
   * The destructor.
   */
  virtual ~MatchboxMElP2lJetJet();
  //@}

public:

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const { return phasespace() ? phasespace()->nDim(3) : 5; }

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return 1; }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return 2; }

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector. Derived classes
   * must call this method once internal degrees of freedom are setup
   * and finally return the result of this method.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return false; }

  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  virtual unsigned int nLight() const { return theQuarkFlavours.size(); }

  /**
   * Return the colour correlated matrix element squared with
   * respect to the given two partons as appearing in mePartonData(),
   * suitably scaled by sHat() to give a dimension-less number.
   */
  virtual double colourCorrelatedME2(pair<int,int>) const;

  /**
   * Return the colour and spin correlated matrix element squared for
   * the gluon indexed by the first argument using the given
   * correlation tensor.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   * The lepton flavours to be considered.
   */
  PDVector theLeptonFlavours;

  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

  /**
   * The Z mass squared
   */
  Energy2 theZMass2;

  /**
   * The Z width squared
   */
  Energy2 theZWidth2;

private:

  /**
   * A user defined scale; if zero, the center of mass
   * energy is used.
   */
  Energy theUserScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<MatchboxMElP2lJetJet> initMatchboxMElP2lJetJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMElP2lJetJet & operator=(const MatchboxMElP2lJetJet &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMElP2lJetJet. */
template <>
struct BaseClassTrait<Herwig::MatchboxMElP2lJetJet,1> {
  /** Typedef of the first base class of MatchboxMElP2lJetJet. */
  typedef Herwig::MatchboxMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MatchboxMElP2lJetJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MatchboxMElP2lJetJet>
  : public ClassTraitsBase<Herwig::MatchboxMElP2lJetJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MatchboxMElP2lJetJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MatchboxMElP2lJetJet is implemented. It may also include several, space-separated,
   * libraries if the class MatchboxMElP2lJetJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MatchboxMElP2lJetJet_H */
