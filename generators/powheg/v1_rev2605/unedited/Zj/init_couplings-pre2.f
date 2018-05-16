      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      include 'pwhg_physpar.h'
      real * 8 masswindow_low,masswindow_high,zmasslow,zmasshigh
      real * 8 powheginput
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
      physpar_ml(1)=0.511d-3
      physpar_ml(2)=0.1057d0
      physpar_ml(3)=1.777d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2

c     set mass windows around Z-mass peak in unit of ph_Zwidth
c     It is used in the generation of the Born phase space
      zmasslow = powheginput("#min_Z_mass")
      zmasshigh = powheginput("#max_Z_mass")
      if(zmasslow.gt.0) then
         ph_Zmass2low=zmasslow**2
      else
         masswindow_low = 15
         ph_Zmass2low=max(0d0,ph_Zmass-masswindow_low*ph_Zwidth)
         ph_Zmass2low=ph_Zmass2low**2
      endif
      if(zmasshigh.gt.0) then
         ph_Zmass2high=zmasshigh**2
      else
         masswindow_high = 15
         ph_Zmass2high=(ph_Zmass+masswindow_high*ph_Zwidth)**2
      endif

      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_unit_e = sqrt(4*pi*ph_alphaem)

      end



