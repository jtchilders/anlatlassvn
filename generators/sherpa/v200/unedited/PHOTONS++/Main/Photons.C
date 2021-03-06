#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/Main/Define_Dipole.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

// define statics
int    PHOTONS::Photons::s_mode          = 2;
bool   PHOTONS::Photons::s_useme         = true;
double PHOTONS::Photons::s_ircutoff      = 1E-3;
int    PHOTONS::Photons::s_ircutoffframe = 0;
double PHOTONS::Photons::s_accu          = 1E-6;
int    PHOTONS::Photons::s_nmax          = std::numeric_limits<double>::max();
double PHOTONS::Photons::s_drcut         = 1000.;
bool   PHOTONS::Photons::s_strict        = false;
double PHOTONS::Photons::s_reducemax     = 1.;

double PHOTONS::Photons::s_alpha                = 0.;
bool   PHOTONS::Photons::s_userunningparameters = false;

// member functions of class Photons

Photons::Photons(Data_Reader* reader) :
  m_name("Photons")
{
  rpa->gen.AddCitation
    (1,"Photons is published under \\cite{Schonherr:2008av}.");
  s_mode          = reader->GetValue<int>("YFS_MODE",2);
  if (s_mode>2) s_mode=2;
  s_useme         = (bool)reader->GetValue<int>("YFS_USE_ME",1);
  s_ircutoff      = reader->GetValue<double>("YFS_IR_CUTOFF",1E-3);
  s_userunningparameters = (bool)reader->GetValue<int>("YFS_USE_RUNNING_PARAMETERS",0);
  std::string irframe
            = reader->GetValue<std::string>("YFS_IR_CUTOFF_FRAME","Multipole_CMS");
  if      (irframe == "Multipole_CMS")      s_ircutoffframe = 0;
  else if (irframe == "Lab")                s_ircutoffframe = 1;
  else if (irframe == "Decayer_Rest_Frame") s_ircutoffframe = 2;
  else {
    s_ircutoffframe = 0;
    msg_Info()<<"value '"<<irframe<<"' for the frame for applying the\n"
              <<"IR cut-off for soft photon radiation unkown ...\n"
              <<"setting it to 'Multipole_CMS' ...\n";
  }
  s_nmax          = reader->GetValue<int>("YFS_MAXEM",-1);
  if (s_nmax<0) s_nmax = std::numeric_limits<double>::max();
  s_drcut         = reader->GetValue<double>("YFS_DRCUT",1000.);
  s_strict        = reader->GetValue<int>("YFS_STRICTNESS",0);
  s_reducemax     = reader->GetValue<double>("YFS_REDUCE_MAXIMUM",1.);
  s_accu          = sqrt(rpa->gen.Accu());
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
		 <<"  Mode: "<<s_mode
		 <<" ,  MEs: "<<(s_mode>1?s_useme:0)
		 <<" ,  nmax: "<<s_nmax
		 <<" ,  strict: "<<s_strict
		 <<" ,  dRcut: "<<s_drcut
		 <<" ,  reducemax: "<<s_reducemax
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<irframe<<" ("<<s_ircutoffframe<<")"
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<"\n}"<<std::endl;
}

Photons::Photons() :
  m_name("Photons")
{
  PRINT_INFO("TODO: check whether running of alphaQED is MSbar");
  PRINT_INFO("TODO: evolve all particle masses in MSbar");
  s_mode          = 2;
  s_useme         = true;
  s_ircutoff      = 1E-3;
  s_ircutoffframe = 0;
  s_nmax          = std::numeric_limits<double>::max();
  s_drcut         = 1000.;
  s_strict        = false;
  s_reducemax     = 1.;
  s_accu          = sqrt(rpa->gen.Accu());
  s_userunningparameters = false;
  m_success       = true;
  m_photonsadded  = false;
  msg_Debugging()<<METHOD<<"(){\n"
		 <<"  Mode: "<<s_mode
		 <<" ,  MEs: "<<(s_mode>1?s_useme:0)
		 <<" ,  nmax: "<<s_nmax
		 <<" ,  strict: "<<s_strict
		 <<" ,  dRcut: "<<s_drcut
		 <<" ,  reducemax: "<<s_reducemax
		 <<" ,  IR cut-off: "<<(s_mode>0?s_ircutoff:0)
		 <<" in frame "<<s_ircutoffframe
		 <<" ,  use running parameters "<<s_userunningparameters
		 <<"\n}"<<std::endl;
}

Photons::~Photons() {
}

bool Photons::AddRadiation(Blob * blob) {

  if (blob->Has(blob_status::needs_extraQED)) {
    ResetAlphaQED();
    Define_Dipole dress(blob);
    dress.AddRadiation();
    m_photonsadded = dress.AddedAnything();
    m_success = dress.DoneSuccessfully();
    if (!blob->MomentumConserved())
      msg_Error()<<"Momentum not conserved:\n"<<*blob<<std::endl;
  }
  return m_photonsadded;
}

