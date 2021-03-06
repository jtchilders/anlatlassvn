#include "DIRE/Shower/Alpha_QCD.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "DIRE/Shower/Shower.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#define ZETA3 1.202056903159594

using namespace DIRE;
using namespace ATOOLS;

Alpha_QCD::Alpha_QCD(const Kernel_Key &key):
  Gauge(key)
{
  p_cpl=(MODEL::Running_AlphaS*)
    p_sk->PS()->Model()->GetScalarFunction("alpha_S");
  m_lc=key.p_rd->GetValue<unsigned int>("CSS_CMODE",1);
  m_Nc=key.p_rd->GetValue<unsigned int>("CSS_NCOL",3);
  m_CF=(m_Nc*m_Nc-1.)/(2.0*m_Nc);
  m_CA=m_Nc;
  m_TR=1.0/2.0;
}

void Alpha_QCD::SetLimits()
{
  Shower *ps(p_sk->PS());
  m_fac=(m_type&1)?ps->CplFac(1):ps->CplFac(0);
  double scale=(m_type&1)?ps->TMin(1):ps->TMin(0);
  double scl(CplFac(scale)*scale);
  m_max=(*p_cpl)[Max(p_cpl->ShowerCutQ2(),scl)];
}

double Alpha_QCD::G2(const double &nf) const
{
  return 3.*(67./18.-sqr(M_PI)/6.)-10./9.*nf/2.;
}

double Alpha_QCD::G3(const double &nf) const
{
  return 9.*(245./6.-134./27.*sqr(M_PI)+11./45.*pow(M_PI,4)+22./3.*ZETA3)
    +3.*nf/2.*(-418./27.+40./27.*sqr(M_PI)-56./3.*ZETA3)
    +4./3.*nf/2.*(-55./3.+16.*ZETA3)-16./27.*sqr(nf/2.);
}

double Alpha_QCD::K(const Splitting &s) const
{
  if (!(s.m_kfac&1)) return 0.0;
  double asf=Coupling(s)/(2.0*M_PI), nf=Nf(s);
  if (!(s.m_kfac&4)) return asf*G2(nf);
  return asf*G2(nf)+sqr(asf)*G3(nf)/4.0;
}

double Alpha_QCD::KMax(const Splitting &s) const
{
  if (!(s.m_kfac&1)) return 0.0;
  double asf=CplMax(s)/(2.0*M_PI);
  if (!(s.m_kfac&4)) return asf*G2(3.);
  return asf*G2(3.)+sqr(asf)*G3(3.)/4.0;
}

double Alpha_QCD::Nf(const Splitting &s) const
{
  return p_cpl->Nf(Scale(s));
}

double Alpha_QCD::CplFac(const double &scale) const
{
  return m_fac;
}

double Alpha_QCD::Coupling(const Splitting &s) const
{
  if (s.m_clu&1) return 1.0;
  double scale(Scale(s));
  double scl(CplFac(scale)*scale);
  if (scl<p_cpl->ShowerCutQ2()) return 0.0;
  double cpl=(*p_cpl)(scl);
  if (cpl>m_max) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at \\mu = "
	       <<sqrt(scale)<<" -> q = "<<sqrt(scl)<<"."<<std::endl;
  }
  return cpl;
}

double Alpha_QCD::CplMax(const Splitting &s) const
{
  return m_max;
}

double Alpha_QCD::Solve(const double &as) const
{
  double t0=p_sk->PS()->TMin(m_type&1);
  t0=Max(p_cpl->ShowerCutQ2(),CplFac(t0)*t0);
  double mur2=p_cpl->WDBSolve(as,t0,sqr(rpa->gen.Ecms()));
  msg_Debugging()<<"\\alpha_s("<<sqrt(mur2)<<") = "
		 <<(*p_cpl)(mur2)<<" / "<<as<<"\n";
  return mur2;
}
