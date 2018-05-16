      subroutine init_couplings
      implicit none
      include "coupl.inc"
      include 'PhysPars.h'
      include 'pwhg_par.h'
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      real * 8 powheginput
      external powheginput
      data called/.false./
      save called
      integer idvecbos,vdecaymode,Vdecmod
      common/cvecbos/idvecbos,vdecaymode
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      data lepmass/0.51099891d-3,0.1056583668d0,1.77684d0/
      

*********************************************************
***********         MADGRAPH                 ************
*********************************************************
c Parameters are read from the MadGraph param_card.dat,
c except the strong coupling constant, which is defined
c somewhere else

      call lh_readin("none")
      call madtophys


c******************************************************
c     Choose the process to be implemented
c******************************************************
c    ID of vector boson produced
      idvecbos=powheginput('idvecbos')
c   decay products of the vector boson
      Vdecmod=powheginput('vdecaymode')
      
      if(idvecbos.eq.24) then
         if (Vdecmod.eq.1) then
            vdecaymode=-11
         elseif (Vdecmod.eq.2) then
            vdecaymode=-13
         elseif (Vdecmod.eq.3) then
            vdecaymode=-15
         else
            write(*,*) 'ERROR: The decay mode you selected ',Vdecmod, 
     $           ' is not allowed '
            call pwhg_exit(1)
         endif
         write(*,*) 
         write(*,*) ' POWHEG: H W+ production and decay ' 
         if (vdecaymode.eq.-11) write(*,*) '  to e+ ve '
         if (vdecaymode.eq.-13) write(*,*) '  to mu+ vmu'
         if (vdecaymode.eq.-15) write(*,*) '  to tau+ vtau'
         write(*,*) 
      elseif(idvecbos.eq.-24) then
         if (Vdecmod.eq.1) then
            vdecaymode=11
         elseif (Vdecmod.eq.2) then
            vdecaymode=13
         elseif (Vdecmod.eq.3) then
            vdecaymode=15
         else
            write(*,*) 'ERROR: The decay mode you selected ',Vdecmod, 
     $           ' is not allowed '
            call pwhg_exit(1)
         endif
         write(*,*) 
         write(*,*) ' POWHEG: H W- production and decay '
         if (vdecaymode.eq.11) write(*,*) '  to e- ve~ '
         if (vdecaymode.eq.13) write(*,*) '  to mu- vmu~'
         if (vdecaymode.eq.15) write(*,*) '  to tau- vtau~'
         write(*,*)    
      else
         write(*,*) 'ERROR: The ID of vector boson you selected ',
     $        idvecbos, ' is not allowed (24: W+ -24: W-)'
         call pwhg_exit(1)
      endif

c     set lepton mass
      decmass=lepmass(Vdecmod)  
           
      if (ph_Wmass2low.lt.decmass**2) then
         write(*,*) 'min_w_mass less than the minimun invariant mass of'
         write(*,*) 'the final-state leptonic system ',decmass
         write(*,*) 'POWHEG aborts'
         call pwhg_exit(-1)
      endif

      end


      subroutine lh_readin(param_name)
c overrides the lh_readin subroutine in MODEL/couplings.f;
c to make it work, rename or delete
c the lh_readin routine in MODEL/couplings.f
      implicit none
      character*(*) param_name
      include 'coupl.inc'
      include 'PhysPars.h'
      double precision  Two, Four, Rt2, Pi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Pi = 3.14159265358979323846d0 )

c     Common to lh_readin and printout
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
      real * 8 powheginput
      external powheginput
c the only parameters relevant for this process are set
c via powheginput. All others are needed for the
c madgraph routines not to blow.

      alpha= 1/132.50698d0
      gfermi = 0.1166390d-4
      alfas = 0.119d0
      zmass = 91.188d0
      tmass = 172.5d0
      lmass = 0d0
      mcMS = 0d0
      mbMS = 0d0
      mtMS = 172.5d0
      mtaMS = 1.777d0
      cmass = 0d0
      bmass = 0d0
      lmass=0d0
      wmass=sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))

      twidth=1.5083d0

      hmass = powheginput('hmass')
      hwidth = powheginput('hwidth')

      ph_Hmass2low=powheginput("min_h_mass")**2
      ph_Hmass2high=powheginput("max_h_mass")**2
      ph_Wmass2low=powheginput("min_w_mass")**2
      ph_Wmass2high=powheginput("max_w_mass")**2    

     
      zwidth=2.441d0
      wwidth=2.0476d0

c     POWHEG CKM matrix
c
c        d     s     b
c    u
c    c
c    t

      Vud=0.97428d0 
      Vus=0.2253d0  
      Vub=0.00347d0 
      Vcd=0.2252d0  
      Vcs=0.97345d0 
      Vcb=0.0410d0  
      Vtd=0.00862d0 
      Vts=0.0403d0  
      Vtb=0.999152d0

c$$$      Vud=1d0
c$$$      Vus=1d-10
c$$$      Vub=1d-10
c$$$      Vcd=1d-10
c$$$      Vcs=1d0
c$$$      Vcb=1d-10
c$$$      Vtd=1d-10
c$$$      Vts=1d-10
c$$$      Vtb=1d0

      end

    

      subroutine madtophys
      implicit none
      include 'coupl.inc'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      real * 8 e_em,g_weak
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb

      
      gal(1) = sqrt(alpha*4*pi)
      e_em=gal(1)      
      ph_unit_e=e_em
      ph_alphaem=e_em**2/(4*pi)
      ph_sthw2=1-(wmass/zmass)**2
      ph_sthw=sqrt(ph_sthw2)
      g_weak=e_em/ph_sthw
      ph_gfermi=sqrt(2d0)*g_weak**2/(8*wmass**2)

      ph_Zmass = zmass
      ph_Wmass = wmass
      ph_Hmass = hmass
      ph_Zwidth = zwidth
      ph_Wwidth = wwidth
      ph_Hwidth = hwidth
      ph_tmass = tmass

      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_Wmass2 = ph_Wmass**2

c     CKM from PDG 2010 (eq. 11.27)
      ph_CKM(1,1)=Vud
      ph_CKM(1,2)=Vus
      ph_CKM(1,3)=Vub
      ph_CKM(2,1)=Vcd
      ph_CKM(2,2)=Vcs
      ph_CKM(2,3)=Vcb
      ph_CKM(3,1)=Vtd
      ph_CKM(3,2)=Vts
      ph_CKM(3,3)=Vtb

      end
