      subroutine setup_PYTHIA_tune
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)
C   
ccc      call PYTUNE(320)
      end
c***********************************************************************
      subroutine setup_PYTHIA_parameters
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      real*8 powheginput
      external powheginput
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
      parameter (mult_inter=.true.)
      integer maxev
      common/mcmaxev/maxev
c     multiple interactions
c     (MI can increase a lot the execution time)
      if(.not.mult_inter) then
      mstp(81)=20   !No Multiple interactions. Force a call to PYEVNW 
      else
      mstp(81)=21   ! MPI on in the PYEVNW MPI scenario
      endif
c     photon radiation off quarks and leptons
c       mstj(41)=12              
c     No photon radiation off quarks and leptons
      mstj(41)=11              

c     if using wgrad2, there must be no photon radiation
      if(powheginput('wgrad2').eq.1)then
      mstj(41)=11    !time-like QCD branching of quarks & gluons 
                     !pg 317 of hep-ph/0603175
      mstp(61)=1     !space-like QCD radiation if IS is hadronic
                     !pg 208
      mstp(71)=1     !master switch for FS radiation
                     !pg 208
      endif
c      mstp(61)=0                !No IS shower
c      mstp(71)=0                !No FS shower
c      mstp(91)=0                !No Primordial kt
c      mstp(131)=0               !No Pile Up
c      mstp(111)=0               !No hadronization
      mstp(64) =3 !use Lambda_MC for IS shower > 6.4.19
c      mstp(64) =1 !use Lambda_MSbar (default)

c     number of warnings printed on the shell
      mstu(26)=20
c     call PYLIST(12) to see the PYTHIA decay table

c     tolerate 2% of killed events
      mstu(22)=maxev/50
      end
c***********************************************************************
      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end
c***********************************************************************
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
      integer maxev
      common/mcmaxev/maxev
      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
      if (lprup(1).eq.10015)  mdcy(pycomp(15),1)=0
      if (lprup(1).eq.9985)  mdcy(pycomp(15),1)=0
      end
c***********************************************************************
      subroutine UPEVNT
      implicit none
      call lhefreadev(97)
      end
c***********************************************************************
      subroutine upveto
c pythia routine to abort event
      end
c***********************************************************************
      subroutine pyabeg
      call init_hist
      end
c***********************************************************************
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
c***********************************************************************
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
c WZGRAD edit begin
      logical flg_saverand
      common/pwhg_flg_saverand/flg_saverand
c WZGRAD edit end
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
c WZGRAD edit begin
      flg_saverand=.true.
c WZGRAD edit end
      call analysis(xwgtup)
      call pwhgaccumup 
      end

