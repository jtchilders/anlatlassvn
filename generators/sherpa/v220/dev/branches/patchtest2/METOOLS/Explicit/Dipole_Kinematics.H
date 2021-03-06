#ifndef METOOLS_Explicit_Dipole_Kinematics_H
#define METOOLS_Explicit_Dipole_Kinematics_H

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Dipole_Info.H"

namespace METOOLS {

  class Dipole_Kinematics {
  protected:

    Current *p_i, *p_j, *p_k, *p_ijt, *p_kt;

    ATOOLS::Vec4D m_pi, m_pj, m_pk, m_Q;

    int m_type, m_swap, m_trig, m_evol;
    double m_kt2c[2];

    double m_z, m_y, m_kt2, m_Q2, m_phase[2], m_res[3];
    double m_mi2, m_mj2, m_mij2, m_mk2, m_ym, m_yp, m_f, m_a;

    Dipole_Info *p_info;
    SizeT_Vector m_pm;

    Current_Vector m_cur;

    ATOOLS::Vec4D_Vector m_p;

    double Lam(const double &s,
	       const double &sb,const double &sc) const;

  public:

    // constructor
    Dipole_Kinematics(Dipole_Info *const info,Current *const i,
		      Current *const j,Current *const k,
		      Current *const ijt,Current *const kt);

    // destructor
    ~Dipole_Kinematics();

    // member functions
    void Evaluate();

    void CheckKT2Min();

    // inline functions
    inline Current *JI() const { return p_i; }
    inline Current *JJ() const { return p_j; }
    inline Current *JK() const { return p_k; }

    inline Current *JIJT() const { return p_ijt; }
    inline Current *JKT() const  { return p_kt;  }

    inline const ATOOLS::Vec4D &PI() const { return m_pi; }
    inline const ATOOLS::Vec4D &PJ() const { return m_pj; }
    inline const ATOOLS::Vec4D &PK() const { return m_pk; }

    inline const ATOOLS::Vec4D &Q() const { return m_Q; }

    inline int Type() const { return m_type; }
    inline int Swap() const { return m_swap; }
    inline int Trig() const { return m_trig; }

    inline double Z() const { return m_z; }
    inline double Y() const { return m_y; }

    inline double YMin() const { return m_ym; }
    inline double YMax() const { return m_yp; }

    inline double KT2() const { return m_kt2; }
    inline double Q2() const  { return m_Q2;  }

    inline double Phase(const int i) const { return m_phase[i]; }
    inline double Res(const int i) const   { return m_res[i];   }

    inline void SetPhase(const double &phase,const int i)
    { m_phase[i]=phase; }
    inline void SetRes(const double &res,const int i)
    { m_res[i]=res; }
    inline void AddRes(const double &res,const int i)
    { m_res[i]+=res; }

    inline void SetCurrents(const Current_Vector &cur)
    { m_cur=cur; m_p.resize(m_cur.size()); }

    inline const ATOOLS::Vec4D_Vector &Momenta() const { return m_p; }

    inline void AddTrig(const int trig) { m_trig&=trig; }
    inline void SetTrig(const int trig) { m_trig=trig;  }

    inline SizeT_Vector &PM() { return m_pm; }

    inline bool Massive() const
    { return m_mi2||m_mj2||m_mk2||m_mij2; }

    inline double F() const { return m_f; }
    inline double A() const { return m_a; }

    inline void SetF(const double &f) { m_f=f; }
    inline void SetA(const double &a) { m_a=a; }

  };// end of class Dipole_Kinematics
  
  std::ostream &operator<<
    (std::ostream &str,const Dipole_Kinematics &k);

}// end of namespace METOOLS

#endif
