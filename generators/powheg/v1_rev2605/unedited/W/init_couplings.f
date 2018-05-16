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
      integer i,j

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_Wmass = powheginput("#Wmass")
      if (ph_Wmass.le.0d0) ph_Wmass  = 80.398d0     
      ph_Wwidth = powheginput("#Wwidth")
      if (ph_Wwidth.le.0d0) ph_Wwidth =  2.141d0

      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/128.89d0
      ph_Zmass = powheginput("#Zmass")
      if (ph_Zmass.le.0d0) ph_Zmass  = 91.1876d0     
      ph_Zwidth = powheginput("#Zwidth")
      if (ph_Zwidth.le.0d0) ph_Zwidth =  2.4952d0
      ph_sthw2 = powheginput("#sthw2")
      if (ph_sthw2.le.0d0) ph_sthw2 = abs(1d0-(ph_Wmass/ph_Zmass)**2)

c     CAVEAT: 
      ph_CKM(1,1) = powheginput("#CKM_Vud")
      if (ph_CKM(1,1).le.0d0) ph_CKM(1,1)=0.975d0 
      ph_CKM(1,2) = powheginput("#CKM_Vus")
      if (ph_CKM(1,2).le.0d0) ph_CKM(1,2)=0.222d0 
      ph_CKM(1,3) = powheginput("#CKM_Vub")
      if (ph_CKM(1,3).le.0d0) ph_CKM(1,3)=1d-10
      ph_CKM(2,1) = powheginput("#CKM_Vcd")
      if (ph_CKM(2,1).le.0d0) ph_CKM(2,1)=0.222d0 
      ph_CKM(2,2) = powheginput("#CKM_Vcs")
      if (ph_CKM(2,2).le.0d0)  ph_CKM(2,2)=0.975d0 
      ph_CKM(2,3) = powheginput("#CKM_Vcb")
      if (ph_CKM(2,3).le.0d0) ph_CKM(2,3)=1d-10
      ph_CKM(3,1) = powheginput("#CKM_Vtd")
      if (ph_CKM(3,1).le.0d0) ph_CKM(3,1)=1d-10
      ph_CKM(3,2) = powheginput("#CKM_Vts")
      if (ph_CKM(3,2).le.0d0) ph_CKM(3,2)=1d-10
      ph_CKM(3,3) = powheginput("#CKM_Vtb")
      if (ph_CKM(3,3).le.0d0) ph_CKM(3,3)=1d0

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
      ph_Wmass2 = ph_Wmass**2


      if(mass_low.ge.0d0) then
         ph_Wmass2low=mass_low**2
      else
         ph_Wmass2low=(max(0d0,ph_Wmass-masswindow_low*ph_Wwidth))**2
      endif
      if(mass_high.gt.0d0) then
         ph_Wmass2high=mass_high
      else
         ph_Wmass2high=ph_Wmass+masswindow_high*ph_Wwidth
      endif
      ph_Wmass2high=min(kn_sbeams,ph_Wmass2high**2)
      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if( ph_Wmass2low.ge.ph_Wmass2high ) then
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
      write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
      write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
      if(ph_runwidth) write(*,*) 'using running W width'
      write(*,*) 'CKM matrix' 
      do i=1,3
         write(*,*) (ph_CKM(i,j),j=1,3)
      enddo
      write(*,*) '*************************************'
      write(*,*)
      write(*,*) '*************************************'
      write(*,*) sqrt(ph_Wmass2low),'< M_W <',sqrt(ph_Wmass2high)
      write(*,*) '*************************************'
      endif
      end

