      subroutine setup_PYTHIA_tune
C   100       A : Rick Field's CDF Tune A                     (Oct 2002)
C   103      DW : Rick Field's CDF Tune DW                    (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro                (Feb 2009)
C   
c     call PYTUNE(320)
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

      mstp(81)=20      !No Multiple interactions. Force a call to PYEVNW 
c         mstp(81)=21   ! MPI on in the PYEVNW MPI scenario


      PARJ(90)= 2*10**4         ! tauola to prevent pythia to radiate photon 
c     off leptons. The parameter, divided by 2,  represents the threshold 
c     in GeV below which leptons do not radiate

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0  ! to reduce number of photons
c     Make tau stable
      mdcy(pycomp(15),1)=0 ! tau stable. Done by Tauola too

c     set stable lighter b-flavoured states: needed to analize single-top
c     events in a reasonable simple way (see the analize subroutine), but,
c     strictly speaking, not necessary for the program.
      mdcy(pycomp(  521 ) ,1)=0
      mdcy(pycomp( -521 ) ,1)=0
      mdcy(pycomp(  511 ) ,1)=0
      mdcy(pycomp( -511 ) ,1)=0
      mdcy(pycomp(  531 ) ,1)=0
      mdcy(pycomp( -531 ) ,1)=0
      mdcy(pycomp(  541 ) ,1)=0
      mdcy(pycomp( -541 ) ,1)=0
      mdcy(pycomp(  553 ) ,1)=0
      mdcy(pycomp(  5212) ,1)=0
      mdcy(pycomp(  5222) ,1)=0
      mdcy(pycomp( -5222) ,1)=0
      mdcy(pycomp(  5112) ,1)=0
      mdcy(pycomp( -5112) ,1)=0
      mdcy(pycomp(  5232) ,1)=0
      mdcy(pycomp( -5232) ,1)=0
      mdcy(pycomp(  5132) ,1)=0
      mdcy(pycomp( -5132) ,1)=0
      mdcy(pycomp(  5122) ,1)=0
      mdcy(pycomp( -5122) ,1)=0
      mdcy(pycomp(  5332) ,1)=0
      mdcy(pycomp( -5332) ,1)=0


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
      integer maxev
      common/mcmaxev/maxev
      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c  SISTEMARE *****************************************
c     if (lprup(1).eq.10015)  mdcy(pycomp(15),1)=0
      end

      subroutine UPEVNT
      implicit none
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

