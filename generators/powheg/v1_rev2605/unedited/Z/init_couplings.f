      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      real * 8 masswindow_low,masswindow_high
      real * 8 mass_low,mass_high
      real * 8 powheginput
      external powheginput
      logical verbose
      parameter(verbose=.true.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/128.89d0
      ph_Zmass = powheginput("#Zmass")
      if (ph_Zmass.le.0d0) ph_Zmass  = 91.1876d0     
      ph_Zwidth = powheginput("#Zwidth")
      if (ph_Zwidth.le.0d0) ph_Zwidth =  2.4952d0
     
      ph_Wmass = powheginput("#Wmass")
      if (ph_Wmass.le.0d0) ph_Wmass  = 80.398d0     
      ph_Wwidth = powheginput("#Wwidth")
      if (ph_Wwidth.le.0d0) ph_Wwidth =  2.141d0

      ph_sthw2 = powheginput("#sthw2")
      if (ph_sthw2.le.0d0) ph_sthw2 = abs(1d0-(ph_Wmass/ph_Zmass)**2)


c     number of light flavors
      st_nlight = 5

c     mass window
      mass_low = powheginput("#mass_low")
      if (mass_low.lt.0d0) mass_low=-1d0
      mass_high = powheginput("#mass_high")
      if (mass_high.lt.0d0) mass_high=-1d0   
      masswindow_low = powheginput("#masswindow_low")
      if (masswindow_low.le.0d0) masswindow_low=30d0
      masswindow_high = powheginput("#masswindow_high")
      if (masswindow_high.le.0d0) masswindow_high=30d0

c     running width
      ph_runwidth = powheginput("#running_width").eq.1d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2

      if(mass_low.ge.0d0) then
         ph_Zmass2low=mass_low**2
      else
         ph_Zmass2low=(max(0d0,ph_Zmass-masswindow_low*ph_Zwidth))**2
      endif
      if (sqrt(ph_Zmass2low).lt.1d0) then
         write(*,*) '*************************************'
         write(*,*) 'WARNING: Z virtuality cutoff at 1 GeV'
         write(*,*) '         to avoid the photon pole    '
         write(*,*) '*************************************'
         ph_Zmass2low=1d0
      endif
     
      if(mass_high.ge.0d0) then
         ph_Zmass2high=mass_high
      else
         ph_Zmass2high=ph_Zmass+masswindow_high*ph_Zwidth
      endif
      ph_Zmass2high=min(kn_sbeams,ph_Zmass2high**2)
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_unit_e = sqrt(4*pi*ph_alphaem)


      if( ph_Zmass2low.ge.ph_Zmass2high ) then
         write(*,*) "Error in init_couplings: mass_low >= mass_high"
         call exit(1)
      endif

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      if(ph_runwidth) write(*,*) 'using running Z width'
      write(*,*) '*************************************'
      write(*,*)
      write(*,*) '*************************************'
      write(*,*) sqrt(ph_Zmass2low),'< M_Z <',sqrt(ph_Zmass2high)
      write(*,*) '*************************************'
      endif
      end

