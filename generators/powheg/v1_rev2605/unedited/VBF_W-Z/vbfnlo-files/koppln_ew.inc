c======================================================
c
c  Define common block containing KOPPLN inputs.
c
c======================================================
c
      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf
      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf

c  xmb = mb pole mass,   BMASS_H = mb(m_H) etc.
      double precision xmb, BMASS_H, xmc, CMASS_H, xmtau 
      common /quarkmasses/ xmb, BMASS_H, xmc, CMASS_H, xmtau 
      logical with_anom 	
      common /anom_switch/ with_anom
      logical with_anomHiggs
      common /anomHiggs/ with_anomHiggs
      logical with_kk 	
      common /kk_switch/ with_kk
      logical with_anomwgaj
      common /anomwgaj/ with_anomwgaj


      logical loopSqr_switch, ewcor_switch, FeynH_switch, slha_switch
      logical loopSqr_switchIN
      integer model, HiggsType, ewapprox, sector, ewscheme, floops, MHloops      

      common /cglobalI/ ewapprox, model, HiggsType, ewscheme, floops
      common /cglobalI/ sector, MHloops
      common /cglobalL/ FeynH_switch, slha_switch, ewcor_switch 
      common /cglobalL/ loopSqr_switch, loopSqr_switchIN
