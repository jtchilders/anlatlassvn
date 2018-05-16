      subroutine setup_PYTHIA_tune
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)
C   
      call PYTUNE(320)
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
      parameter (mult_inter=.true.)
      integer maxev
      common/mcmaxev/maxev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     multiple interactions
c     (MI can increase a lot the execution time)
      if(.not.mult_inter) mstp(81)=0   !No Multiple interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     photon radiation off leptons
c      mstj(41)=3                !No photon radiations off leptons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      mstp(61)=0                !No IS shower
c      mstp(71)=0                !No FS shower
c      mstp(91)=0                !No Primordial kt
c      mstp(131)=0               !No Pile Up
c      mstp(111)=0               !No hadronization

c     Change lambda in alfa running.
c     Keep all commented to let PYTHIA do its default
c     mstp(64) =3 !use Lambda_MC for IS shower
c     mstp(64) =1               !use Lambda_MSbar (default)

c     number of warnings printed on the shell
      mstu(26)=20

c     print the values of some PYTHIA parameters
c      call PYSTAT(5)

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
      integer maxev
      common/mcmaxev/maxev

      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0

      if(lprup(1).eq.1004) then
c make D+ and D- stable (using PDG codes)
         mdcy(pycomp(411),1)=0
         mdcy(pycomp(-411),1)=0
c make D0 and D0bar stable
         mdcy(pycomp(421),1)=0
         mdcy(pycomp(-421),1)=0
      elseif(lprup(1).eq.1005) then
c make B+ and B- stable
         mdcy(pycomp(521),1)=0
         mdcy(pycomp(-521),1)=0
c make B0 and B0bar stable
         mdcy(pycomp(511),1)=0
         mdcy(pycomp(-511),1)=0
      endif
      end

      subroutine UPEVNT
      call lhefreadev(97)
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
      include 'LesHouches.h'
      include 'hepevt.h'
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end
