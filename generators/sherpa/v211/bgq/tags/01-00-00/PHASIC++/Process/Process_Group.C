#include "PHASIC++/Process/Process_Group.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

Process_Group::~Process_Group()
{
  Clear();
}

size_t Process_Group::Size() const
{
  return m_procs.size();
}

Process_Base *Process_Group::operator[](const size_t &i)
{
  return m_procs[i];
}

Weight_Info *Process_Group::OneEvent(const int wmode,const int mode) 
{
  p_selected=NULL;
  if (p_int->TotalXS()==0.0) {
    p_selected=m_procs[int(ATOOLS::ran->Get()*m_procs.size())];
    return p_selected->OneEvent(mode);
  }
  double disc=p_int->SelectionWeight(wmode)*ATOOLS::ran->Get();
  for (size_t i=0;i<m_procs.size();++i) {
    disc-=dabs(m_procs[i]->Integrator()->SelectionWeight(wmode));
    if (disc<=0.) {
      p_selected=m_procs[i];
      return p_selected->OneEvent(mode);
    }
  }
  msg_Error()<<METHOD<<"(): Cannot select any process. xs = "
	     <<p_int->TotalXS()*rpa->Picobarn()<<" pb."<<std::endl;
  return NULL;
}

#ifdef USING__Threading
void *Process_Group::TDifferential(void *arg)
{
  PG_TID *tid((PG_TID*)arg);
  while (true) {
    // wait for group to signal
    pthread_mutex_lock(&tid->m_s_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_cond_signal(&tid->m_s_cnd);
    if (tid->m_s==0) return NULL;
    // worker routine
    tid->m_bd=tid->m_d=0.0;
    if (tid->m_m&4) {
	for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) {
	  tid->m_d+=tid->p_proc->
	    m_mprocs[tid->m_i]->Differential(*tid->p_p);
	  tid->m_bd+=tid->p_proc->m_mprocs[tid->m_i]->LastB();
	}
    }
    else {
	for (tid->m_i=tid->m_b;tid->m_i<tid->m_e;++tid->m_i) {
	  tid->m_d+=tid->p_proc->
	    m_umprocs[tid->m_i]->Differential(*tid->p_p);
	  tid->m_bd+=tid->p_proc->m_umprocs[tid->m_i]->LastB();
	}
    }
    // signal group to continue
    pthread_cond_wait(&tid->m_t_cnd,&tid->m_t_mtx);
  }
  return NULL;
}
#endif

double Process_Group::Differential(const Vec4D_Vector &p)
{
  m_lastb=m_last=0.0;
#ifdef USING__Threading
  if (m_cts.empty()) {
    for (size_t i=0;i<m_procs.size();i++) {
      m_last+=m_procs[i]->Differential(p);
      m_lastb+=m_procs[i]->LastB();
    }
  }
  else {
    // start calculator threads
    size_t d(m_umprocs.size()/m_cts.size());
    if (m_umprocs.size()%m_cts.size()>0) ++d;
    for (size_t j(0), i(0);j<m_cts.size()&&i<m_umprocs.size();++j) {
      PG_TID *tid(m_cts[j]);
      tid->p_p=&p;
      tid->m_m=1;
      tid->m_b=i;
      tid->m_e=Min(i+=d,m_umprocs.size());
      pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    }
    // suspend calculator threads
    for (size_t j(0), i(0);j<m_cts.size()&&i<m_umprocs.size();++j) {
      i+=d;
      PG_TID *tid(m_cts[j]);
      pthread_mutex_lock(&tid->m_t_mtx);
      pthread_mutex_unlock(&tid->m_t_mtx);
      pthread_cond_signal(&tid->m_t_cnd);
      m_last+=tid->m_d;
      m_lastb+=tid->m_bd;
    }
    // start calculator threads
    d=m_mprocs.size()/m_cts.size();
    if (m_mprocs.size()%m_cts.size()>0) ++d;
    for (size_t j(0), i(0);j<m_cts.size()&&i<m_mprocs.size();++j) {
      PG_TID *tid(m_cts[j]);
      tid->p_p=&p;
      tid->m_m=5;
      tid->m_b=i;
      tid->m_e=Min(i+=d,m_mprocs.size());
      pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    }
    // suspend calculator threads
    for (size_t j(0), i(0);j<m_cts.size()&&i<m_mprocs.size();++j) {
      i+=d;
      PG_TID *tid(m_cts[j]);
      pthread_mutex_lock(&tid->m_t_mtx);
      pthread_mutex_unlock(&tid->m_t_mtx);
      pthread_cond_signal(&tid->m_t_cnd);
      m_last+=tid->m_d;
      m_lastb+=tid->m_bd;
    }
  }
#else
  for (size_t i(0);i<m_procs.size();++i) {
    m_last+=m_procs[i]->Differential(p);
    m_lastb+=m_procs[i]->LastB();
  }
#endif
  if (IsNan(m_last))
    msg_Error()<<METHOD<<"(): "<<om::red
		<<"Cross section is 'nan'."<<om::reset<<std::endl;
  return m_last;
}

void Process_Group::SetScale(const Scale_Setter_Arguments &args)
{
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetScale(args);
}
  
void Process_Group::SetKFactor(const KFactor_Setter_Arguments &args)
{
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetKFactor(args);
}

bool Process_Group::IsGroup() const
{
  return true;
}

void Process_Group::Add(Process_Base *const proc) 
{
  if (proc==NULL) return;
  if (m_procmap.find(proc->Name())!=m_procmap.end())
    THROW(critical_error,"Doubled process '"+proc->Name()+"'");
  m_procmap[proc->Name()]=proc;
  m_oew=Max(m_oew,proc->OrderEW());
  m_oqcd=Max(m_oqcd,proc->OrderQCD());
  if (m_nin>0 && m_nout>0 &&
      (m_nin!=proc->NIn() || m_nout!=proc->NOut())) {
    msg_Error()<<METHOD<<"(): Cannot add process '"
	       <<proc->Name()<<"' to group '"<<m_name<<"'.\n"
	       <<"  Inconsistent number of external legs."<<std::endl; 
    return;
  }  
  m_procs.push_back(proc);
}

bool Process_Group::Remove(Process_Base *const proc) 
{
  for (std::vector<Process_Base*>::iterator xsit=m_procs.begin();
       xsit!=m_procs.end();++xsit) {
    if (*xsit==proc) {
      m_procs.erase(xsit);
      return true;
    }
  }
  return false;
}

bool Process_Group::Delete(Process_Base *const proc) 
{
  if (Remove(proc)) {
    delete proc;
    return true;
  }
  return false;
}

void Process_Group::Clear() 
{
  while (m_procs.size()>0) {
    delete m_procs.back();
    m_procs.pop_back();
  }
}

Process_Base *Process_Group::GetProcess(const std::string &name)
{
  std::map<std::string,Process_Base*>::const_iterator 
    pit(m_procmap.find(name));
  if (pit!=m_procmap.end()) return pit->second;
  if (name==m_name) return this;
  for (size_t i(0);i<m_procs.size();++i)
    if (m_procs[i]->IsGroup()) {
      Process_Base *proc(m_procs[i]->Get<Process_Group>()->GetProcess(name));
      if (proc) return proc;
    }
  return NULL;
}

bool Process_Group::CalculateTotalXSec(const std::string &resultpath,
				       const bool create)
{
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
	  m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
	p_int->ISR()->SetPartonMasses(m_flavs);
      }
    }
    psh->InitCuts();
    for (size_t i=0;i<m_procs.size();++i)
      m_procs[i]->BuildCuts(psh->Cuts());
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  }
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  p_int->SetTotal(0);
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"
	    <<m_name<<"' ("<<(p_gen?p_gen->Name():"")<<")"<<std::endl;
#ifdef USING__Threading
  int helpi;
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(helpi,"PG_THREADS")) helpi=4;
  else msg_Info()<<METHOD<<"(): Set number of threads "<<helpi<<".\n";
  if (m_nout<=2) helpi=0;
  else if (m_nout==3) helpi=Min(helpi,2);
  else if (m_nout==4) helpi=Min(helpi,3);
  if (m_mprocs.size()) m_cts.resize(helpi);
  for (size_t i(0);i<m_cts.size();++i) {
    PG_TID *tid(new PG_TID(this));
    m_cts[i] = tid;
    pthread_cond_init(&tid->m_s_cnd,NULL);
    pthread_cond_init(&tid->m_t_cnd,NULL);
    pthread_mutex_init(&tid->m_s_mtx,NULL);
    pthread_mutex_init(&tid->m_t_mtx,NULL);
    pthread_mutex_lock(&tid->m_s_mtx);
    pthread_mutex_lock(&tid->m_t_mtx);
    tid->m_s=1;
    int tec(0);
    if ((tec=pthread_create(&tid->m_id,NULL,&TDifferential,(void*)tid)))
      THROW(fatal_error,"Cannot create thread "+ToString(i));
  }
#endif
  double totalxs(psh->Integrate()/rpa->Picobarn());
#ifdef USING__Threading
  for (size_t i(0);i<m_cts.size();++i) {
    PG_TID *tid(m_cts[i]);
    tid->m_s=0;
    pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    int tec(0);
    if ((tec=pthread_join(tid->m_id,NULL)))
      THROW(fatal_error,"Cannot join thread"+ToString(i));
    pthread_mutex_unlock(&tid->m_t_mtx);
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_mutex_destroy(&tid->m_t_mtx);
    pthread_mutex_destroy(&tid->m_s_mtx);
    pthread_cond_destroy(&tid->m_t_cnd);
    pthread_cond_destroy(&tid->m_s_cnd);
    delete tid;
  }
  m_cts.clear();
#endif
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->TotalXS()!=0.0) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
  }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}

void Process_Group::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetLookUp(lookup);
}

bool Process_Group::CheckFlavours
(const Subprocess_Info &ii,const Subprocess_Info &fi) const
{
  std::vector<Flavour> cfl;
  for (size_t i(0);i<ii.m_ps.size();++i) cfl.push_back(ii.m_ps[i].m_fl);
  for (size_t i(0);i<fi.m_ps.size();++i) cfl.push_back(fi.m_ps[i].m_fl);
  int charge(0), strong(0);
  size_t quarks(0), nin(ii.m_ps.size());
  for (size_t i(0);i<cfl.size();++i) {
    charge+=i<nin?-cfl[i].IntCharge():cfl[i].IntCharge();
    if (abs(cfl[i].StrongCharge())!=8)
      strong+=i<nin?-cfl[i].StrongCharge():cfl[i].StrongCharge();
    quarks+=cfl[i].IsQuark();
    if (quarks>m_pinfo.m_nmaxq) {
      msg_Debugging()<<METHOD<<"(): '"<<GenerateName(ii,fi)<<"': n_q > "
		     <<m_pinfo.m_nmaxq<<". Skip process.\n";
      return false;
    }
  }
  if (quarks<m_pinfo.m_nminq) {
    msg_Debugging()<<METHOD<<"(): '"<<GenerateName(ii,fi)<<"': n_q < "
		   <<m_pinfo.m_nminq<<". Skip process.\n";
    return false;
  }
  if (charge!=0 || strong!=0) return false;
  bool res(true);
  for (size_t i(0);i<fi.m_ps.size();++i) {
    if (fi.m_ps[i].m_ps.empty()) continue;
    Subprocess_Info cii;
    cii.m_ps.push_back(fi.m_ps[i]);
    if (!CheckFlavours(cii,fi.m_ps[i])) {
      res=false;
      break;
    }
  }
  return res;
}

void Process_Group::SetFlavour(Subprocess_Info &cii,Subprocess_Info &cfi,
			       const ATOOLS::Flavour &fl,const size_t i) const
{
  if (i<m_nin) cii.SetExternal(fl,i);
  else cfi.SetExternal(fl,i-m_nin);
}

bool Process_Group::ConstructProcesses(Process_Info pi,const size_t &ci)
{
  if (ci==m_nin+m_nout) {
    if (!CheckFlavours(pi.m_ii,pi.m_fi)) return false;
    SortFlavours(pi);
    std::string name(GenerateName(pi.m_ii,pi.m_fi));
    if (m_procmap.find(name)!=m_procmap.end()) return false;
    Process_Base *proc(GetProcess(pi));
    if (!proc) return false;
    proc->SetGenerator(Generator());
    proc->Init(pi,p_int->Beam(),p_int->ISR());
    if (!Initialize(proc)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<proc->Name()<<"'\n";
      delete proc;
      m_procmap[name]=NULL;
      return false;
    }
    msg_Debugging()<<METHOD<<"(): Init ok '"
		   <<proc->Name()<<"'\n";
    Add(proc);
    return true;
  }
  bool one(false);
  for (size_t j(0);j<m_flavs[ci].Size();++j) {
    SetFlavour(pi.m_ii,pi.m_fi,m_flavs[ci][j],ci);
    if (ConstructProcesses(pi,ci+1)) one=true;
  }
  return one;
}

bool Process_Group::ConstructProcesses()
{
  Process_Info cpi(m_pinfo);
  bool res(ConstructProcesses(cpi,0));
  return res;
}

void Process_Group::SetUpThreading()
{
#ifdef USING__Threading
  for (size_t i(0);i<m_procs.size();++i) {
    m_procs[i]->SetUpThreading();
    if (m_procs[i]->IsMapped()) m_mprocs.push_back(m_procs[i]);
    else m_umprocs.push_back(m_procs[i]);
  }
#endif
}

void Process_Group::SetGenerator(ME_Generator_Base *const gen) 
{ 
  Process_Base::SetGenerator(gen);
  for (size_t i(0);i<m_procs.size();++i) 
    m_procs[i]->SetGenerator(gen);
}

void Process_Group::SetShower(PDF::Shower_Base *const ps) 
{ 
  Process_Base::SetShower(ps);
  for (size_t i(0);i<m_procs.size();++i) 
    m_procs[i]->SetShower(ps);
}

void Process_Group::SetSelector(const Selector_Key &key)
{
  Process_Base::SetSelector(key);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetSelector(key);
}

void Process_Group::SetFixedScale(const std::vector<double> &s)
{
  Process_Base::SetFixedScale(s);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetFixedScale(s);
}

void Process_Group::SetSelectorOn(const bool on)
{
  Process_Base::SetSelectorOn(on);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetSelectorOn(on);
}

void Process_Group::SetUseBIWeight(bool on)
{
  Process_Base::SetUseBIWeight(on);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetUseBIWeight(on);
}

bool Process_Group::Trigger(const Vec4D_Vector &p)
{
  bool trigger=false;
  for (size_t i(0);i<Size();++i)
    if (m_procs[i]->Trigger(p)) trigger=true;
  return trigger;
}
 
void Process_Group::FillOnshellConditions()
{
  Process_Base::FillOnshellConditions();
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->FillOnshellConditions();
}
