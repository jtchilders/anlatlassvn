      subroutine setup_PYTHIA_tune
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)
C   
cc      call PYTUNE(320)
      end

      subroutine setup_PYTHIA_parameters
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
c     multiple interactions
      logical mult_inter
      parameter (mult_inter=.false.)
      integer maxev
      common/mcmaxev/maxev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     multiple interactions
c     (MI can increase a lot the execution time)
      if(.not.mult_inter) then
         mstp(81)=20   !No Multiple interactions. Force a call to PYEVNW 
      else
         mstp(81)=21   ! MPI on in the PYEVNW MPI scenario
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     photon radiation off quarks and leptons
c       mstj(41)=12              
c     No photon radiation off quarks and leptons
c      mstj(41)=11              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      mstp(61)=0                !No IS shower
c      mstp(71)=0                !No FS shower
c      mstp(91)=0                !No Primordial kt
c      mstp(131)=0               !No Pile Up
c      mstp(111)=0               !No hadronization

      mstp(64) =3   !use Lambda_MC for IS shower > 6.4.19
c      mstp(64) =1 !use Lambda_MSbar (default)

c     number of warnings printed on the shell
      mstu(26)=20
c     call PYLIST(12) to see the PYTHIA decay table
ccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     tolerate 2% of killed events
      mstu(22)=maxev/50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      end

      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

      subroutine UPINIT
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      double precision parp,pari
      integer mstp,msti
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MDCY,MDME,KFDP
      double precision brat
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer pycomp
      external pycomp
      integer maxev,i
      common/mcmaxev/maxev
      double precision R_jet,ptmin_jet,powheginput
      external powheginput
      common/cjetdefs/R_jet,ptmin_jet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   JET DEFINITIONS MANDATORY TO DEFINE THE PROCESS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      R_jet=powheginput('R_jet')
      ptmin_jet=powheginput('ptmin_jet')

      nevhep=0

c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
c     Make tau stable
      mdcy(pycomp(15),1)=0
c     set stable lighter b-flavoured states: 
      mdcy(pycomp( 521  ) ,1)=0
      mdcy(pycomp( -521 ) ,1)=0
      mdcy(pycomp( 511  ) ,1)=0
      mdcy(pycomp( -511 ) ,1)=0
      mdcy(pycomp( 531  ) ,1)=0
      mdcy(pycomp( -531 ) ,1)=0
      mdcy(pycomp( 5222 ) ,1)=0
      mdcy(pycomp( 5112 ) ,1)=0
      mdcy(pycomp( 5232 ) ,1)=0
      mdcy(pycomp( -5132) ,1)=0
      mdcy(pycomp( 5132 ) ,1)=0
      mdcy(pycomp( 541  ) ,1)=0
      mdcy(pycomp( -541 ) ,1)=0
      mdcy(pycomp( 553  ) ,1)=0
      mdcy(pycomp( -5112) ,1)=0
      mdcy(pycomp( -5222) ,1)=0
      mdcy(pycomp( -5122) ,1)=0
      mdcy(pycomp( 5332 ) ,1)=0
      mdcy(pycomp( -5232) ,1)=0
      mdcy(pycomp( -5332) ,1)=0
      mdcy(pycomp( 5122 ) ,1)=0
c     top decay:
c     force the top to decay always in a W-b
      do i=41,55
         mdme(i,1)=0
      enddo
      mdme(46,1)=1
c     force the W to decay always in a (e,ve)
      do i=190,209
         mdme(i,1)=0
      enddo
      mdme(206,1)=1
c     call PYLIST(12) to see the PYTHIA decay table
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      end

      subroutine UPEVNT
      implicit none
      include 'LesHouches.h'
      real * 8 scalupfac,powheginput
      external powheginput
      logical ini
      data ini/.true./
      save ini
      call lhefreadev(97)
      scalupfac=powheginput('#scalupfac')
      if (scalupfac.lt.0) scalupfac=1d0
      if ((scalupfac.ne.1d0).and.ini) print *,"######  SCALUPFAC ="
     $     ,scalupfac
      ini=.false.
      scalup=scalup*scalupfac
      end


      subroutine upveto
c pythia routine to abort event
      end

      subroutine pyabeg
      call init_hist
      end

      subroutine pyaend
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      open(unit=99,file=pwgprefix(1:lprefix)//'POWHEG+PYTHIA-output.top'
     #     ,status='unknown')
      call pwhgsetout
      call pwhgtopout
      close(99)
      end


      subroutine pyanal
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer mint
      double precision vint
      COMMON/PYINT1/MINT(400),VINT(400)
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      
      if(mint(51).ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call pylist(7)      !hepeup
            call pylist(2)      !all the event
         endif
         return
      endif
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

