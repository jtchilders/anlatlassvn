#ifndef AHADIC_Tools_Splitting_Tools_H
#define AHADIC_Tools_Splitting_Tools_H

#include "AHADIC++/Tools/Splitting_Functions.H"
#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Dipole.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace AHADIC {
  class Hadronisation_Parameters;

  struct  DecaySpecs { double weight,massmin,popweight; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;
  
  struct leading {
    enum code {
      none               = 0,
      only_quarks        = 1,
      quarks_and_gluons  = 2,
      quarks_and_gluons2 = 3
    };
  };


  struct PTOrder {
    enum code {
      total       = 3,
      gluon_split = 2,
      gluon_emit  = 1,
      none        = 0,
      flat        = -1
    };
  };
  PTOrder::code DefinePTOrder(const int & ptorder);



  class Splitting_Tools {
  private:
    leading::code         m_leading;
    PTOrder::code         m_ptorder;
    const static ATOOLS::Vec3D s_ex, s_ey, s_ez;
    bool                  m_fourquarks, m_analyse;
    int                   m_masstreatment;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    MODEL::Strong_Coupling * p_as; 
    Splitting_Functions    * p_kernels;
    FlavDecayMap           * p_options;

    Proto_Particle   * p_spect, * p_split, * p_out1, * p_out2;
    double             m_mmin_2, m_pt2max, m_pt2max_factor;

    ATOOLS::Poincare   m_cms;
    ATOOLS::Vec4D      m_mom1,m_mom2,m_mom3,m_mom0,m_nperp,m_lperp;
    bool               m_leadsplit,m_glusplit;

    double             m_kt2max,m_lastpt2;
    double             m_m1,m_m2,m_m3,m_m23,m_m1_2,m_m2_2,m_m3_2,m_m23_2;
    double             m_Q2,m_Q,m_Qt2;
    double             m_phi,m_y,m_z,m_s23,m_kt2,m_kt;
 
    ATOOLS::Flavour    m_flav;
    long int           m_tot,m_d,m_s,m_u,m_reject_y,m_reject_z;

    bool   ProduceKinematics(const bool & first,const bool & vetodiquark);
    bool   SelectFlavour(const bool & vetodiquark);
    bool   FixZRange(double & zmin,double & zmax);
    void   SetInfoTagsForOutgoings() const;
    void   AnalyseKinematics(const ATOOLS::Vec4D & q1,const ATOOLS::Vec4D & q2,
			     const ATOOLS::Vec4D & q3);
    bool   ConstructKinematics();
   public:
    Splitting_Tools(const leading::code &,const PTOrder::code &,
		    const ZForm::code &,
		    MODEL::Strong_Coupling * as,const bool &);
    ~Splitting_Tools();

    void   SetSpectatorAndSplitter(Dipole * dip);
    void   SwapSpectatorAndSplitter(Dipole * dip);
    bool   PrepareKinematics(Dipole * dip1,const bool & first=false,
			     const bool & enforce=false);
    bool   DetermineSplitting(Dipole * dip1,const bool & first,
			      const bool & vetodiquark=true);
    void   AftermathOfSplitting(Dipole * dip1);

    void   SetOptions(FlavDecayMap *);
    void   UpdateSpectator();
    void   GetNewParticles(Proto_Particle *& out1,Proto_Particle *& out2) const;

    const double & PT2() const   { return m_kt2; }
    const double & PT() const    { return m_kt; }
    const double & Z() const     { return m_z; }

    Proto_Particle * Spectator() { return p_spect; }
    Proto_Particle * Splitter()  { return p_split; }

    const leading::code & Leading() const { return m_leading; }
  };

  inline void Splitting_Tools::SetOptions(FlavDecayMap * options) {
    p_options = options;
  }

  inline void Splitting_Tools::UpdateSpectator() {
    p_spect->m_mom = m_mom1;
  }

  inline void Splitting_Tools::GetNewParticles(Proto_Particle *& out1,
					       Proto_Particle *& out2) const {
    out1 = p_out1; out2 = p_out2;
  }


  /*! 
    \file
    \brief Contains the class AHADIC::Splitting_Tools and the struct AHADIC::DecaySpecs.
  */
  /*!
    \class Splitting_Tools
    \brief This class is repsonsible for all initial gluon splittings into constituent (di-)quarks
    and for the test decays of clusters.

    In the splittings, the Splitting_Tools invokes a picture borrowed from perturbation theory:
    It splits a dipole with mass $M$ into two according to a distribution 
    \f$\mbox{\rm d}p_\perp^2/p_\perp^2\,\mbox{\rm d}y N\alpha_s(p_\perp)\f$, 
    where \f$p_\perp\f$ is an analogon of transverse momentum and is given by
    \f[ p^2_\perp = \frac{s_{12}s{23}}{s_{123}}\,. \f]
    Here the \f$s_{ij}\f$ are the invariant masses of the summed corresponding vector after the
    dipole splitting.  They relate to the evolution parameters like
    \f[s_{12,23} = Mp_\perp e^{\mp y}\f].  
    The catch here is to use a non-pertubative strong coupling \f$\alpha_s\f$ and to constrain the
    transverse momenta to be below the perturbative IR cutoff given by the parton shower.

    \todo: Implement strategy to catch mistrials.
    \todo: Set the correct \f$p_\perp^{\rm max}\f$ arguments.  
  */
  /*!
    \var Strong_Coupling * Splitting_Tools::p_as
    \brief A pointer to the non-pertubative strong coupling constant \f$\alpha_S\f$, 
    given by AHADIC::Strong_Coupling.  
  */
  /*!
    \var Constituents * Splitting_Tools::p_constituents
    \brief A pointer to all constituent characteristics, made available through
    AHADIC::Constituents.
  */
  /*!
    \var FlavDecayMap * Splitting_Tools::p_options
    \brief A list of all (kinematically) available options (flavours) in gluon splitting 
    processes, including their popping weights \f$w_q\f$.
  */
  /*!
    \var Dipole * Splitting_Tools::p_dip
    \brief A pointer to the dipole to be split.  

    While for direct dipole decays, this is just a pointer to the dipole, this dipole
    is explicitly constructed in cluster decays.  
  */
  /*!
    \var Proto_Particle * Splitting_Tools::p_spect
    \brief The spectator particle in the dipole decay.

    In constructing the kinematics of the decay it will retain its direction and the recoil 
    just amounts to a compression of its momentum.
  */
  /*!
    \var Proto_Particle * Splitting_Tools::p_split
    \brief The splittnig particle in the dipole decay.  

    In gluon splitting this particle will ultimately be deleted, namely in the corresponding
    method Splitting_Tools::SplitDipole(Dipole * dip,const double pt2max).  In gluon
    emission this particle will essentially take the transverse part of the recoil.  

    \todo:  Update the strategy for the selection of Splitting_Tools::p_split and 
    Splitting_Tools::p_spect in case leading particles are present.  
  */
  /*!
    \var Proto_Particle * Splitting_Tools::p_out1
    \brief A newly constructed outgoing particle, emerging in all dipole splitting processes.
  */
  /*!
    \var Proto_Particle * Splitting_Tools::p_out2
    \brief A newly constructed outgoing particle, emerging in dipole splitting through gluon
    splitting only.
  */
  /*!
    \var double m_M2,m_m1,m_m2,m_m3,m_m12,m_m22,m_m32,m_2p3,m_2p3min
    \brief
  */
  /*!
    \var double m_s12,m_s23,m_pt2veto,m_pt2,m_asmax
    \brief
  */
  /*!
    \var double m_xt2,m_xt2max,m_xt2min,m_y,m_ybound
    \brief
  */
  /*!
    \var double Splitting_Tools::m_pref
    \brief The prefactor \f$N_{\rm split}\f$ in the Sudakov form factor.  

    It is given by either \f$N=C_F\f$ in the case of gluon emission or \f$N=\Sigma_q w_q\f$, 
    the sum of all "popping" weights for the accessible quarks/diquarks \f$q\f$ in gluon
    splitting.
  */
  /*!
    \var ATOOLS::Flavour Splitting_Tools::m_flav
    \brief The flavour of the quark-antiquark (antidiquark-diquark) pair in gluon splitting.
  */
  /*!
    \var ATOOLS::Vec4D Splitting_Tools::m_mom1
    \brief The momenta of the first of the three outgoing particles after the dipole has decayed,
    typically this is the spectator momentum, retaining its original direction.
  */
  /*!
    \var ATOOLS::Vec4D Splitting_Tools::m_mom2
    \brief The momenta of the second of the three outgoing particles after the dipole has decayed.
  */
  /*!
    \var ATOOLS::Vec4D Splitting_Tools::m_mom3
    \brief The momenta of the third of the three outgoing particles after the dipole has decayed.
  */

  /*!
    \fn Splitting_Tools::Splitting_Tools(Strong_Coupling *)
    \brief Explicit constructor. 

    It sets pointers to Splitting_Tools::p_as, Splitting_Tools::p_constituents and it initialises 
    the maximal value of \f$\alpha_S\f$, Splitting_Tools::m_asmax, for the hit-or-miss algorithm, 
    and the minimal mass in any splitting, Splitting_Tools::m_2p3.
  */
  /*!
    \fn bool Splitting_Tools::SplitCluster(Cluster * clu,const double pt2max)
    \brief Method to split a cluster in the sequence of firstly emitting a gluon off the original
    dipole and the splitting one of the two emerging (di-)quark-gluon dipoles.  

    First of all, this method initialises a dipole dip1 from the cluster, which first plays
    the role of the active dipole, Splitting_Tools::p_dip.  For this dipole then the method 
    Splitting_Tools::EmitGluon(const double) is called with the maximally allowed \f$p_\perp\f$ 
    (pt2max) as argument.  Having emitted the gluon, it is set as antitriplet in the
    original dipole dip1 and as triplet in the merging dipole dip2.  To proceed one of those
    two dipoles must be selected for the gluon splitting - the selection proceeds according
    to the mass of the two dipoles.  The winner is then treated by the dipole splitting
    method Splitting_Tools::SplitDipole(Dipole * dip,const double pt2max), where the gluon
    is decayed.  Finally, the emerging two dipoles are translated back into clusters,
    with all connections (left, right, previous) established.

    If the gluon emission or the subsequent dipole splitting fails, the method returns "false",
    triggering a new attempt on the full hadronization process.
  */
  /*!
    \fn bool Splitting_Tools::SplitDipole(Dipole * dip,const double pt2max)
    \brief Method to split any dipole.  

    This method is employed directly in the initial gluon splitting at the beginning of the
    hadronisation phase.  Then it is called by the AHADIC::Gluon_Decayer.  It is also
    needed in the emission of a gluon in cluster splitting and the subsequent splitting of
    the gluon.  Then its actions are steered by Splitting_Tools::EmitGluon(const double) and
    Splitting_Tools::SplitCluster(Cluster * clu,const double pt2max), respectively. 

    The central idea is that first all relevant parameters for the dipole splitting are set.  In
    particular, this involves the squared masses \f$M^2\f$ and \f$m_1^2\f$ of the dipole
    and the spectator, the normalised minimal transverse momentum
    \f[x^2_{\perp,{\rm min}} = \frac{(m_1+m_{\rm min}/2)^2(m_{\rm min})^2}{M^2}\,,\f]
    where \f$m_{\rm min}\f$ is the minimally allowed amss (\f$2m_u\f$ or similar) in
    gluon splitting.  The rapidity bounds are given by
    \f[y_\pm = \mp \frac12\log(x^2_{\perp,{\rm min}})\,,\f]
    and the prefactor \f$N\f$, Splitting_Tools::m_pref, is determined dynamically from the 
    allowed (di-)quark species as \f$N=\Sigma_q w_q\f$.

    Then, the actual splitting is fixed, through the corresponding method
    Splitting_Tools::DetermineSplitting(bool glusplit=true).  The flag there refers to the 
    specific process: gluon emission (false) or splitting (true).   If a viable splitting
    has been found, the kinematics is completely determined, this is done in
    Splitting_Tools::FixKinematics(bool glusplit=true) with the same assignemt for the flag.

    If either the preparation of the dipole for the splitting, the determination of a viable 
    splitting or of its kinematics fails, the method returns false.  In this case the
    full hadronization process should be repeated.
  */
  /*!
    \fn void Splitting_Tools::UpdateSpectator()
    \brief Sets the momentum of the spectator Proto_Particle, Splitting_Tools::p_spect,
    to the internal momentum \f$p_1\f$, Splitting_Tools::m_mom1, calculated and boosted in
    Splitting_Tools::ConstructKinematics() and Splitting_Tools::FixKinematics().
  */
  /*!
    \fn void Splitting_Tools::GetNewParticles(Proto_Particle *& out1,Proto_Particle *& out2)
    \brief Returns the two new particles.

    \todo Maybe better: Return a dipole.
  */
  /*!
    \fn void Splitting_Tools::SetOptions(FlavDecayMap *)
    \brief Sets the options allowed for gluon splitting, as defined in AHADIC::Gluon_Decayer.
  */
  /*!
    \fn Dipole * Splitting_Tools::EmitGluon(const double)
    \brief Method to emit a gluon from a cluster-dipole.

    This method is somewhat the gluon emission version of the full dipole splitting method
    Splitting_Tools::SplitDipole(Dipole * dip,const double pt2max).  Similarly to there,
    first, splitter and spectator are determined, before all relevant invariants, including
    the new flavour, Splitting_Tools::m_flav, are fixed.  They are given by the masses of 
    the dipole and the involved particles, \f$ M^2, m_{1,2,3}, m_{1,2,3}^2\f$, the
    minimal mass of the 23-system, given by \f$m_3\f$, the minimal normalised
    transverse momentum
    \f[x^2_{\perp,{\rm min}} = \frac{(m_1+m_{\rm min}/2)^2(m_3+m_{\rm min}/2)^2}{M^2}\,,\f]
    where the minimal mass for a splitting gluon, \f$m_{\rm min}\f$ is distributed in
    equal proportion over the 12 and the 23 system, in anticipation of the subsequent
    gluon splitting.  The rapidity bounds are given by
    \f[y_\pm = \mp \frac12\log(x^2_{\perp,{\rm min}})\,,\f]
    and the prefactor in the Sudakov form factor, Splitting_Tools::m_pref, is given by 
    \f$C_F = 4/3\f$.

    Then Splitting_Tools::DetermineSplitting(bool glusplit=true) is invoked to emit the
    gluon, before Splitting_Tools::FixKinematics(bool glusplit=true) fixes the kinematics.
    Note that in this case of gluon emission the flag is set to "false" in both methods.  If
    any of both method fails, NULL is returned, otherwise, in the case of success, the 
    emerging 23 dipole is returned.  
  */
  /*!
    \fn void Splitting_Tools::SetSpectatorAndSplitter()
    \brief This method fixes splitter and spectator in dipole decays.
  */
  /*!
    \fn bool Splitting_Tools::DetermineSplitting(bool glusplit=true)
    \brief Determines the actual splitting of a dipole, either a gluon emission (flag is set to
    "false"), or a gluon splitting ("true").

    In an infinite loop Splitting_Tools::SelectPT_Y(bool glusplit=true) is called.  This method 
    will return true as long as a \f$p_\perp > p_{\perp,{\rm min}}\f$ can be found.  If this 
    is not the case, the class will try a minimal splitting with minimal transverse momentum,
    as coded in Splitting_Tools::MinimalDecay(bool glusplit=true).  Only if this also fails,
    the method returns "false", which leads to a repetition of the hadronization step.  
    While the loop delivers viable transverse momenta, the hit-or-miss algorithm corrects
    on physically correct splittings through Splitting_Tools::Verto().
  */
  /*!
    \fn bool Splitting_Tools::SelectPT_Y(bool glusplit=true)
    \brief Determines a naive estimate for a pair of transverse momentum and rapidity.

    The normalised transverse momentum is given by the previous transverse momentum
    and a random number as
    \f[x_\perp^2 = 
    \exp\left\{-\sqrt{\log^2(x^2_{\perp,{\rm old}})-
    \frac{2\log(\#)}{N\alpha_{s, {\rm max}}}}\right\}\f]\,.
    If either argument of the root in the exponent of the above expression is negative or if
    this yields a \f$p_\perp<p_{\perp,{\rm min}}\f$ the method returns false, ultimately
    invoking the forced decay coded in Splitting_Tools::MinimalDecay(bool glusplit=true).
    Otherwise a rapidity \f$y\f$ is selected uniformly in the boundaries and a flavour
    is selected according to its relative "popping" weight.
  */
  /*!
    \fn bool Splitting_Tools::Veto()
    \brief This organises the correction weight in the hit-or-miss algorithm and yields "true"
    for a vetoed emission.

    The vetos are:
    If the rapidity \f$y\f$ is outside the true boundaries given by
    \f[y_\pm^{\rm true} = -log(2\sqrt{x^2_\perp})\cdot\left(1+\sqrt{1-4*x^2_\perp}\right)\,.\f]

    If the kinematics does not work out.  This is determined in 
    Splitting_Tools::ConstructKinematics().

    A veto is performed on the ratio \f$\alpha_s(p_\perp^{\rm phys})/\alpha_s^{\rm max}\f$.

    And emission with a two large physical transverse momentum 
    \f$p_\perp^{\rm phys}>p_\perp^{\rm veto}\f$ are vetoed.
  */
  /*!
    \fn void Splitting_Tools::MinimalDecay(bool glusplit=true)
    \brief The back-up solution of a dipole decay with minimal \f$p_\perp\f$, in case
    the Sudakov form factor fails.  This is some sort of an "enforced" splitting.
    The flag indicates, whether this method is invoked for gluon splitting ("true") or
    gluon emission ("false").  
    
    In the case of gluon splitting, a viable (di-)quark is selected at random 
    (\todo: Maybe fix this to u,d?) and the masses are set.  In the case of gluon 
    emission, the gluon is assumed to have the minimal viable mass.  Then, in both 
    cases, the invariants are given by 
    \f[ s_{23} = (m_2+m_3)^2\;\;\;\mbox{\rm and}\;\;\;
        s_{12} = m_1^2+m_2^2+\frac{m_2}{m_2+m_3}(M^2-m_1^2-s_{23})
               = m_1^2+m_2^2+\frac{1}{2}(M^2-m_1^2-m_2^2-m_3^2-2m_2m_3)\,.\f]
  */
  /*!
    \fn void Splitting_Tools::CalculateInvariants()
    \brief Calculates \f$s_{12}\f$ and \f$s_{23}\f$ from the (evolution) transverse momentum 
    and the rapidity. 

    The relations read 
    \f[s_{12,23} = M^2\sqrt{x^2_\perp}\exp(\mp y)\,.\f]
  */
  /*!
    \fn bool Splitting_Tools::ConstructKinematics()
    \brief Constructs the decay kinematics.

    After some checks whether the \f$s_{ij}\f$ related to the splitting are feasible, the momenta
    of the dipole splitting are constructed using a Sudakov decomposition as
    \f[ p_1     = l^-+a_1l^+ \f]
    \f[ p_{23}  = l^++a_{23}l^- = p_2+p_3 \f]
    \f[ p_{2,3} = z_{2,3}l^++y_{2,3}a_{23}l^-\pm p_\perp\,.\f]

    Obviously, here, \f$p_1\f$ denotes the spectator and \f$p_{23}\f$ the splitter.
    The invariants needed here are given by
    \f[ a_1 = m_1^2/\gamma\;\;\;\mbox{\rm and}\;\;\; a_{23} = s_{23}/\gamma\,,\f]
    where 
    \f[ \gamma = \frac{\tilde M^2+\sqrt{\tilde M^4-4s_23m_1^2}}{2} \f]
    \f[ \tilde M^2 = M^2-s_{23}-m_1^2\,.\f]
    Thus, the \pm vectors are given by
    \f[ l^\pm = \frac{\sqrt{\gamma}}{2}\left(1,0,0,\pm 1\right)^T\,.\f]
    The splitting variables \f$z_{2,3}\f$ satisfy \f$z_2+z_3 = 1\f$ with
    \f[ z_2 = \frac{s_{12}-(m_1^2+m_2^2)-a_1[s_{23}-(m_3^2-m_2^2)]}{\gamma-a_1s_{23}}\,.\f]
    and, similarly, \f$y_2+y_3 = 1\f$ with
    \f[ y_{2,3} = \frac{m_{2,3}^2+p_\perp^2}{z_{2,3}s_{23}}\,.\f]
    From the \f$z\f$, the physical transverse momentum can be deduced:
    \f[ p_\perp^{\rm phys} = z_2z_3s_{23}-z_3m_2^2-z_2m_3^2\,.\f]
    With a random azimuthal angle \f$\phi\f$ this defines the transverse momentum vector.
  */
  /*!
    \fn bool Splitting_Tools::FixKinematics(bool glusplit=true)
    \brief This boosts and rotates the three momenta \f$p_{1,2,3}\f$, Splitting_Tools::m_mom1, 
    Splitting_Tools::m_mom2, and Splitting_Tools::m_mom3, which have previously been constructed 
    in Splitting_Tools::ConstructKinematics() into the correct frame.

    The trick here is that the system of the momenta \f$q_{1,2}\f$ and \f$p_{1,2,3}\f$
    before and after the splitting live in different systems, to overcome this problem,
    both are brought into the same system (c.m. system of the splitting dipole, with
    constituents along the \f$z\f$-axis).  The inverse of the transformations applied on the 
    \f$q_i\f$ are the napplied to the \f$p_i\f$.
  */
}

#endif
