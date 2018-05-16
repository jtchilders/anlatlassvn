      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'    
      include 'nlegborn.h'      
      include 'pwhg_kn.h'      
      real * 8 masswindow_low,masswindow_high
      logical verbose
      parameter(verbose=.true.)
      integer i,j
      real *8 powheginput
      external powheginput
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_Hmass = powheginput('hmass')
      ph_Hwidth = powheginput('hwidth')


      ph_GF= powheginput('#gfermi') 
      if (ph_GF.le.0d0) ph_GF  = 0.116639D-04     
      ph_topmass = powheginput('#topmass')
      if (ph_topmass.le.0d0) ph_topmass  = 171.3d0
      ph_bmass = powheginput('#bmass')
      if (ph_bmass.le.0d0) ph_bmass  = 4.55d0
      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/137.035999679d0
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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      
      ph_Hmass2 = ph_Hmass**2

c     set mass windows around H-mass peak in unit of ph_Hwidth
c     It is used in the generation of the Born phase space
C     masswindow is an optonal  parameter passed by the user
C     the default vale is 10 
      masswindow_low = powheginput("#masswindow_low")
      if(masswindow_low.lt.0d0) masswindow_low=10d0
      ph_Hmass2low=max(0d0,ph_Hmass-masswindow_low*ph_Hwidth)
      ph_Hmass2low= ph_Hmass2low**2
      masswindow_high = powheginput("#masswindow_high")
      if(masswindow_high.lt.0d0) masswindow_high=10d0
      ph_Hmass2high=ph_Hmass+masswindow_high*ph_Hwidth
      ph_Hmass2high= min(kn_sbeams,ph_Hmass2high**2)
      ph_HmHw = ph_Hmass * ph_Hwidth


      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'H mass = ',ph_Hmass
      write(*,*) 'H width = ',ph_Hwidth
       write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'GF = ',ph_GF
      write(*,*) 'top mass = ',ph_topmass
      write(*,*) '*************************************'
      write(*,*)
      write(*,*) '*************************************'
      write(*,*) sqrt(ph_Hmass2low),' < M_H <',sqrt(ph_Hmass2high)
      write(*,*) '*************************************'
      endif
      end




