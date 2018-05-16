      subroutine setup_PYTHIA_tune
C --------------------------------------------------------------- C
C - N.B. PYTUNE(ITUNE) must be called before the call to PYINIT - C
C --------------------------------------------------------------- C
C -  100         A :  Rick Field's CDF Tune A     (Oct 2002)
C -  103        DW :  Rick Field's CDF Tune DW    (Apr 2006)
C -  320 Perugia 0 :  Perugia update of S0-Pro    (Feb 2009)
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
      integer maxev
      common/mcmaxev/maxev

c      mstj(41) =3 ! Photon radiation off leptons. Not recommended to touch
c                  ! this -- causes interference with the pT ordered shower.
c      mstp(61) =0 ! No IS shower
c      mstp(71) =0 ! No FS shower
c      mstp(81) =0 ! No Multiple interactions (MI increases execution time).
c      mstp(91) =0 ! No Primordial kt
c      mstp(131)=0 ! No Pile Up
c      mstp(111)=0 ! No hadronization
      
c     N.B.
c     ====
c     For the case of jet production the following parameter setting
c     limits the transverse momentum of secondary scatterings, due
c     to multiple parton interactions, to be less than that of the
c     primary interaction (see POWHEG Dijet paper arXiv:1012.3380
c     [hep-ph] sec. 4.1 and also the PYTHIA Manual).
      mstp(86)=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     maximum number of errors before pythia aborts (def=10)
      mstu(22)=10
c     number of warnings printed on the shell
      mstu(26)=20
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
c     Make tau stable
      mdcy(pycomp(15),1)=0 ! tau stable. Done by Tauola too
      end

      subroutine UPEVNT
      implicit none
      call lhefreadev(97)
c      call lhefinitemasses      
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
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end
