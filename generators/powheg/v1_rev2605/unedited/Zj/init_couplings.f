      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      include 'pwhg_physpar.h'
c     to define here rad_bottomthr2
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      real * 8 masswindow_low,masswindow_high
      real *8 powheginput,pwhg_alphas
      external powheginput,pwhg_alphas
      flg_withdamp=.true.
      flg_bornzerodamp=.true.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_Zmass  = 91.188d0     
      ph_Zwidth =  2.486d0

      ph_alphaem = 1d0/128.930d0
      ph_sthw2 = 0.2312d0

c     number of light flavors
      st_nlight = 5
c Masses of light leptons for Z decays:
c     masses for reshuffling procedure of
c     outgoing particles.
c     If no keywords, will use old defaults, so
c     that backward compatibility is preserved
      physpar_ml(1)=powheginput('#lhfm/emass')
      if(physpar_ml(1).lt.0) physpar_ml(1)=0.000511d0
      physpar_ml(2)=powheginput('#lhfm/mumass')
      if(physpar_ml(2).lt.0) physpar_ml(2)=0.1057d0
      physpar_ml(3)=powheginput('#lhfm/taumass')
      if(physpar_ml(3).lt.0) physpar_ml(3)=1.777d0
ccccccccccccccccccccccccccccccc
c     to include Z->bbar
c     old defaults for b and c are 0 GeV, and they are
c     set in init_phys.f. We want now that 
c     - if vdecaymode=5, the bottom mass has always to be set 
c     to a nonzero value (and also the default value has to be
c     nonzero).
c     In this case, we also want to add a correction (1+as/pi).
c     By default this is switched on (ZbbNLO keyword).
c     - if vdecaymode is not 5, we want backward compatibility,
c     (as of revision 2364 this means that there is no reshuffling at all
c     for b and c quarks)
      ZbbNLO=.false.
      if(powheginput("#vdecaymode").eq.5) then
         ZbbNLO=powheginput('#ZbbNLO').ne.0
         if(ZbbNLO) then
            ZbbNLOfac=1+pwhg_alphas(ph_Zmass**2,st_lambda5MSB,5)/pi
            write(*,*) '*****************************************'
            write(*,*) 'Factor (1+as/pi) included for Z->bb decay'
            write(*,*) 'Numerical values is ',ZbbNLOfac
            write(*,*) '*****************************************'
         endif
c$$$  physpar_mq(4)=powheginput('#lhfm/cmass')
c$$$  if(physpar_mq(4).lt.0) physpar_mq(4)=1.5
         physpar_mq(5)=powheginput('#lhfm/bmass')
         if(physpar_mq(5).lt.0) physpar_mq(5)=5.0
c     to improve a bit on splittings that cannot be reshuffled...
         rad_bottomthr2=(2*physpar_mq(5))**2
      endif
ccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2

c     set mass windows around Z-mass peak in unit of ph_Zwidth
c     It is used in the generation of the Born phase space
      masswindow_low = 25
      masswindow_high = 35
      ph_Zmass2low=max(0d0,ph_Zmass-masswindow_low*ph_Zwidth)
      ph_Zmass2low=ph_Zmass2low**2
      ph_Zmass2high=(ph_Zmass+masswindow_high*ph_Zwidth)**2
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_unit_e = sqrt(4*pi*ph_alphaem)

      end



