#ifndef COMIX_Cluster_Cluster_Algorithm_H
#define COMIX_Cluster_Cluster_Algorithm_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "COMIX/Amplitude/Amplitude.H"
#include "COMIX/Cluster/Color_Setter.H"

namespace PDF { class Cluster_Definitions_Base; }

namespace PHASIC { class Subprocess_Info; }

namespace COMIX {

  struct Cluster_Info {
    Vertex *p_v;
    size_t m_k, m_oew, m_oqcd;
    PDF::CParam m_kt2;
    ATOOLS::Flavour m_mofl;
    Cluster_Info(Vertex *const v,const size_t &k,
		 const PDF::CParam &kt2=PDF::CParam(-1.0,-1.0),
		 const size_t &oew=0,const size_t &oqcd=0,
		 const ATOOLS::Flavour &mofl=ATOOLS::Flavour(kf_none)): 
      p_v(v), m_k(k), m_oew(oew), m_oqcd(oqcd), m_kt2(kt2), m_mofl(mofl) {}
    inline bool operator<(const Cluster_Info &ci) const
    { return p_v<ci.p_v?true:(p_v>ci.p_v?false:m_k<ci.m_k); }
  };// end of struct Cluster_Info

  typedef std::pair<size_t,size_t> Cluster_Key;

  typedef std::pair<Cluster_Key,Cluster_Info>     ClusterInfo_Pair;
  typedef std::multimap<Cluster_Key,Cluster_Info> ClusterInfo_Map;

  typedef std::map<size_t,std::map
		   <size_t,std::map
		    <size_t,std::map<ATOOLS::Flavour,
				     PDF::CParam> > > > Double_Map;

  typedef std::map<size_t,ATOOLS::Vec4D> Vec4D_Map;
  typedef std::map<size_t,size_t>        SizeT_Map;

  typedef std::map<Cluster_Info,double> Vertex_Set;

  typedef std::pair<size_t,double> KT2_Info;
  typedef std::vector<KT2_Info>    KT2Info_Vector;

  class Single_Process;

  class Cluster_Algorithm {
  private:

    Color_Setter m_cs;

    ATOOLS::Mass_Selector *p_ms;

    ATOOLS::Cluster_Amplitude     *p_ampl;
    PDF::Cluster_Definitions_Base *p_clus;

    Amplitude      *p_bg;
    Single_Process *p_xs;

    SizeT_Map m_id;

    bool   m_nosol;
    size_t m_wmode;
    double m_lfrac;

    ATOOLS::ColorID GetColor(Current *const j,
			     Current *const fcur) const;

    int Connected(const Vertex_Vector &dvs,const size_t &idk) const;
    
    PDF::CParam GetMeasure(const size_t &idi,const size_t &idj,
			   const size_t &idk,const ATOOLS::Flavour &mofl,
			   Double_Map &kt2,const SizeT_Map &cid,int cut);

    void CalculateMeasures(const size_t &step,const Vertex_Set &nocl,
			   const Current_Vector &curs,Current *const fcur,
			   ClusterInfo_Map &cinfo,Double_Map &kt2,
			   const SizeT_Map &cid);
    bool CombineWinner(const Cluster_Info &ci,Current_Vector &ccurs,
		       Current *&fcur,ClusterInfo_Map &cinfo);
    bool ClusterStep(const size_t &step,Vertex_Set &nocl,
		     Current_Vector &ccurs,Current *&fcur,
		     ClusterInfo_Map &cinfo,Double_Map &kt2);
    bool Cluster(const size_t &step,const Vertex_Set &nocl,
		 const Current_Vector &ccurs,Current *const fcur,
		 const ClusterInfo_Map &cinfo,KT2Info_Vector &kt2ord,
		 const int complete);

    bool SetRandomColors(Single_Process *const xs);
    bool SetLargeNCColors(Single_Process *const xs);
    bool SetColors(Single_Process *const xs);

    void SetNMax(ATOOLS::Cluster_Amplitude *const ampl,
		 const size_t &id,const size_t &nmax) const;

    KT2Info_Vector UpdateKT2(const KT2Info_Vector &kt2ord,
			     const Cluster_Amplitude *ampl,
			     const int mode=0) const;

  public:

    Cluster_Algorithm(ATOOLS::Mass_Selector *const ms);

    ~Cluster_Algorithm();

    bool Cluster(Single_Process *const xs,const size_t &mode);

    // inline functions
    inline ATOOLS::Cluster_Amplitude *GetAmplitude() { return p_ampl; }

    inline void SetClusterDefinitions
    (PDF::Cluster_Definitions_Base *const cb) { p_clus=cb; }

  };// end of class Cluster_Algorithm

}// end of namespace COMIX

#endif
