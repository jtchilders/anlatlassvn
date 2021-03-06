      subroutine init_couplings
      implicit none
      include "coupl.inc"
      include 'PhysPars.h'
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      real * 8 powheginput
      external powheginput
      data called/.false./
      save called
      if(called) then
         return
      else
         called=.true.
      endif

*********************************************************
***********         MADGRAPH                 ************
*********************************************************
c Parameters are read from the MadGraph param_card.dat,
c except the strong coupling constant, which is defined
c somewhere else
      call setpara("param_card.dat",.true.)

      call madtophys

c      ph_Hmass = powheginput('hmass')
c      ph_Hwidth = powheginput('hwidth')

C      write(*,*) 'hmass', ph_Hmass
C      pause

c Are these needed?
c$$$      physpar_ml(1)=0.511d-3
c$$$      physpar_ml(2)=0.1057d0
c$$$      physpar_ml(3)=1.777d0
c$$$      physpar_mq(1)=0.33d0     ! up
c$$$      physpar_mq(2)=0.33d0     ! down
c$$$      physpar_mq(3)=0.50d0     ! strange
c$$$      physpar_mq(4)=1.50d0     ! charm
c$$$      physpar_mq(5)=4.80d0     ! bottom
      call golem_initialize
      call golem_initialize
      end


      subroutine lh_readin(param_name)
c overrides the lh_readin subroutine in MODEL/couplings.f;
c to make it work, rename or delete
c the lh_readin routine in MODEL/couplings.f
      implicit none
      character*(*) param_name
      include 'coupl.inc'
      double precision  Two, Four, Rt2, Pi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud
      real * 8 powheginput
c the only parameters relevant for this process are set
c via powheginput. All others are needed for the
c madgraph routines not to blow.
      alpha= 7.7585538055706d-03
      gfermi = 0.1166390d-4
      alfas = 0.119d0
      zmass = 91.188d0
      tmass = 174.3d0
      lmass = 0d0
      mcMS = 0d0
      mbMS = 0d0
      mtMS = 174d0
      mtaMS = 1.777d0
      vud = 1d0
      cmass = 0d0
      bmass = 0d0
      lmass=0d0
      hmass = powheginput('hmass')
c      hmass = 120d0
      wmass=sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      twidth=1.5083d0
      hwidth = powheginput('hwidth')
c      hwidth =  5.75308848d-03
      zwidth=2.441d0
      wwidth=2.0476d0
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
      e_em=gal(1)
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


c     CKM from PDG 2010 (eq. 11.27)
      ph_CKM(1,1)=0.97428d0
      ph_CKM(1,2)=0.2253d0
      ph_CKM(1,3)=0.00347d0
      ph_CKM(2,1)=0.2252d0
      ph_CKM(2,2)=0.97345d0
      ph_CKM(2,3)=0.0410d0
      ph_CKM(3,1)=0.00862d0
      ph_CKM(3,2)=0.0403d0
      ph_CKM(3,3)=0.999152d0

      end




C     initializes all the couplings in the  virtual, code
C     generated by GoSam and sets them equal to the
C     ones defined in the POWHEG BOX.
      subroutine golem_initialize
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      integer ierr
      integer ioerr
      character * 20 param
      character * 20 value
      character * 50 line

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

c$$$      param = 'mW='
c$$$      write(value,'(F20.10)') ph_Wmass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

      write(*,*) 'ph_hmass',ph_Hmass

      param = 'mH='
      write(value,'(F20.10)') ph_Hmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)

c$$$      param = 'mC='
c$$$      write(value,'(F20.10)') ph_cmass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)
c$$$
c$$$      param = 'mB='
c$$$      write(value,'(F20.10)') ph_bmass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)
c$$$
c$$$      param = 'mBMS='
c$$$      write(value,'(F20.10)') ph_bmass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

      param = 'mT='
      write(value,'(F20.10)') ph_tmass
      line = trim(param)//trim(adjustl(value))
      call OLP_Option(line,ierr)
      call check_gosam_err(param,ierr)

c$$$      param = 'me='
c$$$      write(value,'(F20.10)') ph_emass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)
c$$$
c$$$      param = 'mmu='
c$$$      write(value,'(F20.10)') ph_mumass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)
c$$$
c$$$      param = 'mtau='
c$$$      write(value,'(F20.10)') ph_taumass
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

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

c$$$      param = 'wB='
c$$$      write(value,'(F20.10)') ph_bwidth
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

c$$$      param = 'wT='
c$$$      write(value,'(F20.10)') ph_twidth
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

c$$$      param = 'wtau='
c$$$      write(value,'(F20.10)') ph_tauwidth
c$$$      line = trim(param)//trim(adjustl(value))
c$$$      call OLP_Option(line,ierr)
c$$$      call check_gosam_err(param,ierr)

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

      call OLP_Start('../GoSam_POWHEG/orderfile.olc',ioerr)
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
      
C     initializes all the couplings in the  virtual, code
C     generated by GoSam and sets them equal to the
C     ones defined in the POWHEG BOX.              
      subroutine golem_initialize
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      integer ierr
      integer ioerr
      character * 20 param
      character * 20 value
      character * 50 line
      
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
      write(value,'(F20.10)') ph_tmass
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
      
      call OLP_Start('../GoSam_POWHEG/orderfile.olc',ioerr)
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
