#ifndef COMIX_Phasespace_PS_Channel_H
#define COMIX_Phasespace_PS_Channel_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace PHASIC { class Color_Integrator; }

namespace COMIX {

  class Process_Base;
  class PS_Current;

  typedef std::vector<PHASIC::Vegas*>          Vegas_Vector;
  typedef std::map<std::string,PHASIC::Vegas*> Vegas_Map;

  typedef std::map<const Current*,PHASIC::Vegas*> CVegas_Map;
  typedef std::map<size_t,PHASIC::Vegas*>         IVegas_Map;
  typedef std::map<size_t,CVegas_Map>             ICVegas_Map;

  typedef std::map<size_t,std::string> PSId_Map;

  typedef std::map<size_t,std::vector<int> > CId_Map;

  typedef std::map<size_t,const METOOLS::Current*> STCC_Map;

#ifdef USING__Threading
  class PS_Channel;

  struct CDBG_PS_TID {
    pthread_t m_id;
    PS_Channel *p_psc;
    size_t m_s, m_n, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    PSId_Map *p_psid;
    CId_Map *p_cid;
    CDBG_PS_TID(PS_Channel *const psc): 
      p_psc(psc), m_s(2), m_b(0), m_e(0),
      p_psid(new PSId_Map()), p_cid(new CId_Map()) {}
    ~CDBG_PS_TID() { delete p_psid; delete p_cid; }
  };// end of struct CDBG_PS_TID

  typedef std::vector<CDBG_PS_TID*> CDBG_PS_TID_Vector; 
#endif

  class PS_Channel: public PHASIC::Single_Channel {
  protected:

    Process_Base   *p_xs;
    Current_Matrix *p_cur;

    Vertex_Vector  m_vtc;
    Vegas_Map      m_vmap;

    CVegas_Map  m_pcmap;
    IVegas_Map  m_pimap;
    ICVegas_Map m_sicmap, m_ticmap;

    PHASIC::Cut_Data *p_cuts;

    ATOOLS::Vec4D_Vector m_p;

    Double_Vector m_rns, m_wrns;
    Vegas_Vector  m_vgs, m_wvgs;

    size_t m_n, m_nr, m_num, m_lid, m_rid;
    double m_texp, m_sexp, m_thexp, m_mfac;
    double m_aexp, m_srbase;

    int m_bmode, m_omode, m_tmode, m_vmode, m_zmode, m_czmode;
    int m_nvints, m_vsopt, m_nopt;

    SP(PS_Generator) p_gen;

    PSId_Map *p_psid;
    CId_Map  *p_cid;
    
#ifdef USING__Threading 
    CDBG_PS_TID_Vector m_cts;
    pthread_mutex_t    m_vgs_mtx, m_wvgs_mtx;

    static void *TGenerateWeight(void *arg);

    CDBG_PS_TID *GetTId() const;
#endif 

    const std::string      &GetPSId(const size_t &id);
    const std::vector<int> &GetCId(const size_t &id);

    inline size_t CIdCount(const size_t &id) { return GetCId(id).size(); }

    size_t SId(const size_t &id) const;
    double SCut(const size_t &id);

    void FillMoms(const size_t &aid,Int_Vector &cur,size_t n);

    double PropMomenta(const Current *cur,const size_t &id,
		       const double &smin,const double &smax,
		       const double *rn);
    double PropWeight(const Current *cur,const size_t &id,
		      const double &smin,const double &smax,const double &s);

    void TChannelBounds(const size_t &aid,const size_t &lid,
			double &ctmin,double &ctmax,
			const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			const double &s1,const double &s2);
    void SingleTChannelBounds(const size_t &aid,const size_t &lid,
			      double &ctmin,double &ctmax,
			      const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			      const double &s1,const double &s2,const int mode);
    void TChannelMomenta(Current *cur,const size_t &id,const size_t &aid,
			 const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			 ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2,
			 const double &s1,const double &s2,
			 const double *rns);
    double TChannelWeight(Current *cur,const size_t &id,const size_t &aid,
			  const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    void SChannelBounds(const size_t &id,const size_t &lid,
			double &ctmin,double &ctmax);
    void SChannelMomenta(Current *cur,const int type,
			 const ATOOLS::Vec4D &pa,ATOOLS::Vec4D &p1,
			 ATOOLS::Vec4D &p2,const double &s1,
			 const double &s2,const double *rns);
    double SChannelWeight(Current *cur,const int type,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    bool GeneratePoint(Current *const ja,Current *const jb,
		       Current *const jc,Vertex *const v,size_t &nr);
    double GenerateWeight(Current *const ja,Current *const jb,
			  Current *const jc,Vertex *const v,size_t &nr);
    
    bool GenerateWeight(PS_Current *const cur);
    bool GenerateWeight();

    bool GeneratePoint(const size_t &id,size_t &nr,
		       Vertex_Vector &v);
    bool GeneratePoint(Vertex_Vector v);

    bool GenerateChannel(Current *const cur,
			 Vertex_Vector &v);
    bool GenerateChannel(Vertex_Vector &v);
    bool GenerateChannels();

    bool Zero(METOOLS::Vertex *const vtx) const;

    PHASIC::Vegas *GetVegas(const std::string &tag,int ni=-1);

    PHASIC::Vegas *GetPVegas(const Current *cur,const size_t &id);
    PHASIC::Vegas *GetSVegas(const size_t &id,const Current *cur);
    PHASIC::Vegas *GetTVegas(const size_t &id,const Current *cur);

  public :

    // constructor
    PS_Channel(const size_t &nin,const size_t &nout,
	       ATOOLS::Flavour *fl,Process_Base *const ps); 

    // destructor
    ~PS_Channel();

    // member functions
    void GeneratePoint(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts,double *rn);
    void GenerateWeight(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts);

    void ISRInfo(int &type,double &m,double &w);
    void ISRInfo(std::vector<int> &ts,
		 std::vector<double> &ms,std::vector<double> &ws) const;

    int  ChNumber();
    void SetChNumber(int n);

    std::string ChID();

    size_t NChannels() const;

    void AddPoint(double value);

    void MPISync();
    void Optimize();
    void EndOptimize();

    void WriteOut(std::string pid);
    void ReadIn(std::string pid);

  };// end of class PS_Channel

}// end of namespace COMIX

#endif
