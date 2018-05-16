      subroutine init_couplings
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      double precision masswindow_low,masswindow_high
      logical verbose
      parameter(verbose=.true.)
      integer i,j
      real * 8 powheginput
      external powheginput

      ph_topmass = powheginput("#topmass")
      if (ph_topmass.le.0d0) ph_topmass = 172.9d0   
     
      ph_topmass2 = ph_topmass**2


c     From now on there are quantities that are relevant for the
c     top decay only.  

      ph_topwidth = powheginput("#topwidth")
      if(ph_topwidth.lt.0d0)  ph_topwidth = 1.31d0

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
      ph_unit_e = sqrt(4*pi*ph_alphaem)

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
      
      
c     masses for reshuffling procedure of
c     outgoing particles
      physpar_ml(1)=powheginput('#lhfm/emass')
      if(physpar_ml(1).lt.0) physpar_ml(1)=0.000511
      physpar_ml(2)=powheginput('#lhfm/mumass')
      if(physpar_ml(2).lt.0) physpar_ml(2)=0.1056
      physpar_ml(3)=powheginput('#lhfm/taumass')
      if(physpar_ml(3).lt.0) physpar_ml(3)=1.777

      physpar_mq(4)=powheginput('#lhfm/cmass')
      if(physpar_mq(4).lt.0) physpar_mq(4)=powheginput('#charmthr')
      if(physpar_mq(4).lt.0) physpar_mq(4)=1.5
      physpar_mq(5)=powheginput('#lhfm/bmass')
      if(physpar_mq(5).lt.0) physpar_mq(5)=powheginput('#bottomthr')
      if(physpar_mq(5).lt.0) physpar_mq(5)=5.0

c     mass window (ONLY ZERO WIDTH APPROX IMPLEMENTED UP TO NOW)
      masswindow_low = powheginput("#masswindow_low")
      if (masswindow_low.le.0d0) masswindow_low=0d0
      masswindow_high = powheginput("#masswindow_high")
      if (masswindow_high.le.0d0) masswindow_high=0d0

c     set mass windows around top-mass peak in unit of ph_topwidth
c     May be used in the generation of the Born phase space in the future

      ph_topmass2low=max(0d0,ph_topmass-masswindow_low*ph_topwidth)
      ph_topmass2low= ph_topmass2low**2
      ph_topmass2high=ph_topmass+masswindow_high*ph_topwidth
      ph_topmass2high= min(kn_sbeams,ph_topmass2high**2)
      ph_topmtopw = ph_topmass * ph_topwidth




      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'top mass = ',ph_topmass
      write(*,*) 'top width = ',ph_topwidth
      write(*,*) 'masses = ',(kn_masses(i),i=1,nlegreal)
      write(*,*) sqrt(ph_topmass2low),'< mt <',sqrt(ph_topmass2high)
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
      write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
      write(*,*) 'CKM matrix' 
      do i=1,3
         write(*,*) (ph_CKM(i,j),j=1,3)
      enddo
      write(*,*) '*************************************'
      write(*,*)
      endif

      end




