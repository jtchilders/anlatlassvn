      subroutine init_couplings
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'   
      include 'pwhg_st.h'   
      include 'pwhg_flg.h'   
      include 'pwhg_math.h'   
      include 'masses.f'   
      include 'ewcouple.f'   
      include 'qcdcouple.f'   
      include 'alfacut.f'
      include 'includect.f'
      integer ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,ii13
      common/inputprocind/ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,
     1     ii9,ii10,ii11,ii12,ii13
      real * 8 lotopdecaywidth,powheginput,pwhg_alphas
      external lotopdecaywidth,powheginput,pwhg_alphas
c Bornzerodamp on by default
      flg_withdamp=.true.
      flg_bornzerodamp=.true.

c alfafcut parameters (MCFM Nagy)
      aff=1d0
      aif=1d0
      afi=1d0
      aii=1d0
c
      includect=.false.
c
      ph_tmass=powheginput("#tmass")
      if(ph_tmass.lt.0) then
         ph_tmass = 172d0
      endif
c      ph_twidth=1.4381971d0	! Set by call to MCFM routine below
      ph_bmass=powheginput("#bmass")
      if(ph_bmass.lt.0) then
         ph_bmass = 4.75d0
      endif
      write(*,*) ' b mass set to ',ph_bmass
      ph_wmass=80.398d0
      ph_wwidth=2.1054d0
      ph_zmass=91.1876d0
      Gf=1.16639d-5
      gwsq=8d0*ph_wmass**2*Gf/sqrt(2d0)
      gw=sqrt(gwsq)

c sin th Cabibbo, pdg as of June 2012
c No other CKM values needed here
      ph_CKM(1,2)=0.2257d0
      ph_CKM(2,1)=0.2257d0

      
      st_alpha=pwhg_alphas(ph_zmass**2,st_lambda5MSB,-1)
      write(6,*) 'Setting st_alpha at mZ = ',st_alpha

      
      as=st_alpha
      gsq=4d0*pi*as
      ason2pi=as/(2d0*pi)
C     
      mt=ph_tmass
      mb=ph_bmass
      wmass=ph_wmass
      wwidth=ph_wwidth
c Calculate twidth
      ph_twidth=lotopdecaywidth(ph_tmass,ph_bmass,ph_wmass,ph_wwidth)
      twidth=ph_twidth

      write(*,*) ' top width at LO, as computed by MCFM routine, =',
     1     twidth

c POWHEG BOX wants particle masses to be stored here
      kn_masses(ii1)=0
      kn_masses(ii2)=0
      kn_masses(ii3)=ph_tmass
      kn_masses(ii4)=ph_tmass
      kn_masses(ii5)=ph_wmass
      kn_masses(ii6)=ph_wmass
      kn_masses(ii7)=0
      kn_masses(ii8)=0
      kn_masses(ii9)=0
      kn_masses(ii10)=0
      kn_masses(ii11)=ph_bmass
      kn_masses(ii12)=ph_bmass

c kn_minmass must be set to the minimum invariant mass
c of the final state (2 mt in this case, mt+mb in single top, etc.)
      kn_minmass=2*ph_tmass


      end
