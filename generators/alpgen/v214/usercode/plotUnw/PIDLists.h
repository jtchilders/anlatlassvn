#ifndef PIDLISTS_H
#include <vector>
#include <string>



class PDG_PID{
   enum PID{
      UNDEFINED   = 0,

      // quarks
      d           = 1,
      u           = 2,
      s           = 3,
      c           = 4,
      b           = 5,
      t           = 6,
      b_prime     = 7,
      t_prime     = 8,

      // leptons
      e              = 11,
      nu_e           = 12,
      mu             = 13,
      nu_mu          = 14,
      tau            = 15,
      nu_tau         = 16,
      tau_prime      = 17,
      nu_tau_prime   = 18,

      // gauge bosons
      g              = 21,
      gamma          = 22,
      Z              = 23,
      W              = 24,
      h0             = 25,
      Z_prime        = 32,
      Z_pprime       = 33,
      W_prime        = 34,
      H0             = 35,
      A0             = 36,
      H_plus         = 37,

      // light mesons
      pi0            = 111,
      pi_plus        = 211,

      // light baryons
      p              = 2212,
      n              = 2112

   };

   static const std::map<unsigned int,std::string> pidString;
   static std::map<unsigned int,std::string> SetPIDStrings(){
      std::map<unsigned int,std::string> tmp;

      // leptons
      tmp[PDG_PID::e]            = "e";
      tmp[PDG_PID::nu_e]         = "nu_e";
      tmp[PDG_PID::mu]           = "mu";
      tmp[PDG_PID::nu_mu]        = "nu_mu";
      tmp[PDG_PID::tau]          = "tau";
      tmp[PDG_PID::nu_tau]       = "nu_tau";
      tmp[PDG_PID::tau_prime]    = "tau_prime";

      // gaug bosons
      tmp[PDG_PID::g]            = "g";
      tmp[PDG_PID::gamma]        = "gamma";
      tmp[PDG_PID::Z]            = "Z";
      tmp[PDG_PID::W]            = "W";
      tmp[PDG_PID::h0]           = "h0";
      tmp[PDG_PID::Z_prime]      = "Z_prime";
      tmp[PDG_PID::Z_pprime]     = "Z_pprime";
      tmp[PDG_PID::W_prime]      = "W_prime";
      tmp[PDG_PID::H0]           = "H0";
      tmp[PDG_PID::A0]           = "A0";
      tmp[PDG_PID::H_plus]       = "H_plus";
      
      return tmp;
   }

   static const std::vector<unsigned int> quarks;
   static std::vector<unsigned int> SetQuarks(){
      std::vector<unsigned int> tmp;
      tmp.push_back(1);
      tmp.push_back(2);
      tmp.push_back(3);
      tmp.push_back(4);
      tmp.push_back(5);
      tmp.push_back(6);
      return tmp;
   }

   static const std::vector<unsigned int> quarks

}pid_list;

PIDLists::quarks = PIDLists::SetQuarks();

#endif

