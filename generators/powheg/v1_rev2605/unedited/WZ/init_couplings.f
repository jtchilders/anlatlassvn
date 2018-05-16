      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'constants.f'
      include 'zwcouple.f'  
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'anomcoup.f'
      include 'ckm.f'
      include 'cabibbo.f'
      include 'nwz.f'
      include 'cvecbos.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      real * 8 masswindow_low,masswindow_high
      logical verbose
      parameter(verbose=.true.)
      integer j,srint
      real * 8 powheginput
      external powheginput
      real * 8 deltas,asmzopi
      physpar_ml(1)=0.511d-3
      physpar_ml(2)=0.1057d0
      physpar_ml(3)=1.777d0
      physpar_mq(1)=0.33d0     ! up
      physpar_mq(2)=0.33d0     ! down
      physpar_mq(3)=0.50d0     ! strange
      physpar_mq(4)=1.50d0     ! charm
      physpar_mq(5)=4.80d0     ! bottom


      call smcouplings

c     number of light flavors
      st_nlight = 5


      !TM added QCD couplings
      gsq = st_alpha*4d0*pi
      as  = st_alpha
      ason2pi = st_alpha/2d0/pi
      ason4pi = st_alpha/4d0/pi


      
      !TM added z couplings
      Q(-5)=+0.333333333333333d0
      Q(-4)=-0.666666666666667d0
      Q(-3)=+0.333333333333333d0
      Q(-2)=-0.666666666666667d0
      Q(-1)=+0.333333333333333d0
      Q(0)=+0d0
      Q(+1)=-0.333333333333333d0
      Q(+2)=+0.666666666666667d0
      Q(+3)=-0.333333333333333d0
      Q(+4)=+0.666666666666667d0
      Q(+5)=-0.333333333333333d0
      tau=(/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/)
      esq = ph_unit_e**2
      zmass = ph_Zmass
      zwidth = ph_Zwidth
      wmass = ph_Wmass
      wwidth = ph_Wwidth

      call couplz(ph_sthw2)

      xw = ph_sthw2
      gwsq = ph_unit_e**2/ph_sthw2
      gw = dsqrt(gwsq)
      write(*,*)'GW',gw,xw

      ! TM for the different processes the
      ! ew couplings need to be set, as in
      ! chooser.f. For now, ee,mumu,tautau
      ! ---really should depend on idvecdecay
      ! and ideally would be in init_process,
      ! but the above constatns need to be set
      if (((vdecaymodeZ).eq.11).or.((vdecaymodeZ).eq.13).or
     $     .((vdecaymodeZ).eq.15)) then
      q1=-1d0
      l1=le
      r1=re
      elseif(((vdecaymodeZ).eq.12).or.
     . ((vdecaymodeZ).eq.14).or.((vdecaymodeZ).eq.16)) then
      q1=0d0
      l1=ln
      r1=rn
      endif

      ! TM set CKM values
      Vud = 0.974d0
      call setckmmatrix

      !Set anom coulings
      call setanomcoup
      if (anomtgc) then
         write(*,*) 'Anomalous couplings used:'
         write(*,*) 'Delta_g1(Z)', delg1_z
         write(*,*) 'Delta_g1(Gamma)', delg1_g
         write(*,*) 'Lambda(Z)', lambda_z
         write(*,*) 'Lambda(Gamma)', lambda_g
         write(*,*) 'Delta_K(Gamma)', delk_g
         write(*,*) 'Delta_K(Z)', delk_z
      endif


      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'e**2  = ',ph_unit_e**2
      write(*,*) 'Vud,Vus,Vcd,Vcs',Vud,Vus,Vcd,Vcs
      write(*,*) '*************************************'
      endif

      end



