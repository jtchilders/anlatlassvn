// -*- C++ -*-
//
// PowhegInclusiveReweight.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegInclusiveReweight_H
#define HERWIG_PowhegInclusiveReweight_H
//
// This is the declaration of the PowhegInclusiveReweight class.
//

#include "Herwig++/MatrixElement/Matchbox/Powheg/ME2byDipoles.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief PowhegInclusiveReweight is used to transform
 * a subtraction dipole into a BBar real emission contribution.
 *
 */
class PowhegInclusiveReweight: public ME2byDipoles {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PowhegInclusiveReweight();

  /**
   * The destructor.
   */
  virtual ~PowhegInclusiveReweight();
  //@}

public:

  /**
   * Evaluate the ratio.
   */
  virtual double evaluate() const;

  /**
   * Return true, if 'Born screening' should be done
   */
  bool bornScreening() const { return theBornScreening; }

  /**
   * Switch on 'Born screening'
   */
  void doBornScreening() { theBornScreening = true; }

  /**
   * Switch off 'Born screening'
   */
  void noBornScreening() { theBornScreening = false; }

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * True, if 'Born screening' should be done
   */
  bool theBornScreening;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegInclusiveReweight & operator=(const PowhegInclusiveReweight &);

};

}

#endif /* HERWIG_PowhegInclusiveReweight_H */
