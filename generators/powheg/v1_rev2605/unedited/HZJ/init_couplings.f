      subroutine init_couplings
      implicit none
      include "coupl.inc"
      include 'PhysPars.h'
      real * 8 powheginput
      external powheginput
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      data called/.false./
      save called
      integer idvecbos,vdecaymode,Vdecmod
      common/cvecbos/idvecbos,vdecaymode
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      data lepmass/0.51099891d-3,0.1056583668d0,1.77684d0/
      real * 8 ger,gel
      common/Zlepcoupl/ger,gel
      real * 8 t3lep,qlep,gev,gea

      if(called) then
         return
      else
         called=.true.
      endif

c******************************************************
c     Choose the process to be implemented
c******************************************************

      idvecbos = 23
c     decay products of the vector boson
      Vdecmod=powheginput('vdecaymode')
      
      if (Vdecmod.eq.1) then
         vdecaymode=-11
      elseif (Vdecmod.eq.2) then
         vdecaymode=-13
      elseif (Vdecmod.eq.3) then
         vdecaymode=-15
      elseif (Vdecmod.eq.4) then
         vdecaymode=-12
      elseif (Vdecmod.eq.5) then
         vdecaymode=-14
      elseif (Vdecmod.eq.6) then
         vdecaymode=-16
      else
         write(*,*) 'ERROR: The decay mode you selected ',Vdecmod, 
     $        ' is not allowed '
         call pwhg_exit(1)
      endif
      write(*,*) 
      write(*,*) ' POWHEG: H Z J production and decay ' 
      if (vdecaymode.eq.-11) write(*,*) '         to e+ e- '
      if (vdecaymode.eq.-13) write(*,*) '         to mu+ mu-'
      if (vdecaymode.eq.-15) write(*,*) '         to tau+ tau-'
      if (vdecaymode.eq.-12) write(*,*) '         to antinue nue'
      if (vdecaymode.eq.-14) write(*,*) '         to antinumu numu'
      if (vdecaymode.eq.-16) write(*,*) '         to antinutau nutau'


*********************************************************
***********         MADGRAPH                 ************
*********************************************************
c Parameters are read from the MadGraph param_card.dat,
c except the strong coupling constant, which is defined
c somewhere else
      call setpara("param_card.dat",.true.)
      call madtophys
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set parameter for MadGraph
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Here we assign the coupling of the Zll vertex to the neutrino ones, if 
c     the code has to generate Z -> nu nu
      if (Vdecmod.gt.3) then
         gzl(1)=gzn(1)
         gzl(2)=gzn(2)
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set parameter for GoSam: Z couplings to final-state leptons
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      if (Vdecmod.le.3) then
c     electron         
         t3lep =-1d0/2   
         qlep = -1d0
      else
c     neutrino
         t3lep =+1d0/2   
         qlep = 0d0
      endif
      gev = (t3lep - 2*qlep*ph_sthw**2)/(2*ph_sthw*sqrt(1-ph_sthw**2))
      gea = t3lep/(2*ph_sthw*sqrt(1-ph_sthw**2))
      gel=gev+gea
      ger=gev-gea

      call golem_initialize
 


c     set lepton mass
      if (Vdecmod.gt.3) then
         decmass=0d0
      else
         decmass=lepmass(Vdecmod)   
      endif
         
      if (ph_Zmass2low.lt.4*decmass**2) then
         write(*,*) 'min_z_mass less than the minimun invariant mass of'
         write(*,*) 'the final-state leptonic system ',2*decmass
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

      ph_Zmass2low=powheginput("min_z_mass")**2
      ph_Zmass2high=powheginput("max_z_mass")**2


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

      subroutine set_ebe_couplings
      implicit none
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include "coupl.inc"
c QCD coupling constant
      G=sqrt(st_alpha*4d0*pi)
      GG(1)=-G
      GG(2)=-G

c HEFT coupling
      gh(1) = cmplx( g**2/4d0/PI/(3d0*PI*V), 0d0)
      gh(2) = cmplx( 0d0                   , 0d0)
      ga(1) = cmplx( 0d0                   , 0d0)
      ga(2) = cmplx( g**2/4d0/PI/(2d0*PI*V), 0d0)
      gh4(1) = G*gh(1)
      gh4(2) = G*gh(2)
      ga4(1) = G*ga(1)
      ga4(2) = G*ga(2)

      return
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
      ph_ZmZw = ph_Zmass * ph_Zwidth

      ph_Wmass2 = ph_Wmass**2
      ph_Zmass2 = ph_Zmass**2

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


      
C     initializes all the couplings in the  virtual, code
C     generated by GoSam and sets them equal to the
C     ones defined in the POWHEG BOX.              
      subroutine golem_initialize
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_rnd.h'
      integer ierr
      integer ioerr
      character * 20 param
      character * 20 value
      character * 50 line
      character * 29 path
      real * 8 powheginput
      external powheginput
      integer parallelstage,rndiwhichseed
      common/cpwhg_info/parallelstage,rndiwhichseed
      logical massivetop

      rndiwhichseed=rnd_iwhichseed
      parallelstage=powheginput("#parallelstage")
C     Read from card of top loops should be included
      massivetop = .false.
      if (powheginput("#massivetop").eq.1) massivetop=.true.

C     Parameter definition
      param = 'Nf='
      write(value,'(I1)') st_nlight
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'alpha='
      write(value,'(F20.10)') ph_alphaem
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'GF='
      write(value,'(F20.10)') ph_gfermi
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mZ='
      write(value,'(F20.10)') ph_Zmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mW='
      write(value,'(F20.10)') ph_Wmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mH='
      write(value,'(F20.10)') ph_Hmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mC='
      write(value,'(F20.10)') ph_cmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mB='
      write(value,'(F20.10)') ph_bmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mBMS='
      write(value,'(F20.10)') ph_bmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mT='
      if(massivetop) then
         write(value,'(F20.10)') ph_tmass
      else
         write(value,'(F20.10)') 0d0
      end if
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'me='
      write(value,'(F20.10)') ph_emass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mmu='
      write(value,'(F20.10)') ph_mumass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'mtau='
      write(value,'(F20.10)') ph_taumass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wZ='
      write(value,'(F20.10)') ph_Zwidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wW='
      write(value,'(F20.10)') ph_Wwidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wH='
      write(value,'(F20.10)') ph_hwidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wB='
      write(value,'(F20.10)') ph_bwidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wT='
      write(value,'(F20.10)') ph_twidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'wtau='
      write(value,'(F20.10)') ph_tauwidth
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VUD='
      write(value,'(F20.10)') ph_CKM(1,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVDU='
      write(value,'(F20.10)') ph_CKM(1,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VUS='
      write(value,'(F20.10)') ph_CKM(1,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVSU='
      write(value,'(F20.10)') ph_CKM(1,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VUB='
      write(value,'(F20.10)') ph_CKM(1,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVBU='
      write(value,'(F20.10)') ph_CKM(1,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VCD='
      write(value,'(F20.10)') ph_CKM(2,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVDC='
      write(value,'(F20.10)') ph_CKM(2,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VCS='
      write(value,'(F20.10)') ph_CKM(2,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVSC='
      write(value,'(F20.10)') ph_CKM(2,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VCB='
      write(value,'(F20.10)') ph_CKM(2,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVBC='
      write(value,'(F20.10)') ph_CKM(2,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VTD='
      write(value,'(F20.10)') ph_CKM(3,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVDT='
      write(value,'(F20.10)') ph_CKM(3,1)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VTS='
      write(value,'(F20.10)') ph_CKM(3,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVST='
      write(value,'(F20.10)') ph_CKM(3,2)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'VTB='
      write(value,'(F20.10)') ph_CKM(3,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
      param = 'CVBT='
      write(value,'(F20.10)') ph_CKM(3,3)
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)
      
C     Initialize virtual code
      path = '../GoSam_POWHEG/orderfile.olc'
      
      call OLP_Start(path,ioerr,parallelstage,rndiwhichseed)
      call check_gosam_err('olp_start routine',ierr)
      end
      
      
      subroutine check_gosam_err(param,ierr)
      implicit none
      character *(*) param
      integer ierr
      if (ierr.lt.0) then
         write(*,*)
     $        'GoSam '//param(1:len_trim(param))// ' reports an error.'
         write(*,*) 'The POWHEG BOX aborts'
         call exit(1)
      endif
      end
