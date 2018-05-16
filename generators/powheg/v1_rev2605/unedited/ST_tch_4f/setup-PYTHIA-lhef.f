      subroutine setup_PYTHIA_tune
C --------------------------------------------------------------- C
C - N.B. PYTUNE(ITUNE) must be called before the call to PYINIT - C
C --------------------------------------------------------------- C
C   100         A : Rick Field's CDF Tune A      (Oct 2002)
C   103        DW : Rick Field's CDF Tune DW     (Apr 2006)
C   320 Perugia 0 : "Perugia" update of S0-Pro   (Feb 2009)
c$$$      call PYTUNE(320)
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
      double precision dummy,powheginput
      external powheginput
      
c      dummy=powheginput("#dummy")


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Useful settings to interface POWHEG with PYTHIA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     to use external PDF
c     ...................
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     multiple interactions
c     (MI can increase a lot the execution time)
      if(.not.mult_inter) then
         mstp(81)=20            ! No MPI, but FORCE a call to PYEVNW (pt-ordering !!)
      else
         mstp(81)=21            ! MPI switched on, in the PYEVNW MPI scenario
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     photon radiation off quarks and leptons
c      mstj(41)=12
c     No photon radiation off quarks and leptons
       mstj(41)=11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      mstp(61)=0                !No IS shower
c      mstp(71)=0                !No FS shower
c      mstp(91)=0                !No Primordial kt
c      mstp(131)=0               !No Pile Up
c      mstp(111)=0               !No hadronization

c       mstp(64) =3 !use Lambda_MC for IS shower
c       mstp(64) =1 !use Lambda_MSbar (default)

c$$$c     number of events printed on the shell
c$$$c      maxpr=2

c     number of warnings printed on the shell
      mstu(26)=20
ccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     tolerate 2% of killed events
c      mstu(22)=maxev/50
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
      integer KCHG
      double precision PMAS,PARF,VCKM
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      integer pycomp
      external pycomp
      integer maxev
      common/mcmaxev/maxev
      integer i
      nevhep=0
c read the header first, so lprup is set
      call lhefreadhdr(97)

c     Make PI0 stable as in herwig default
      mdcy(pycomp(111),1)=0
c     Make taus stable
      mdcy(pycomp(  15 ), 1)=0
      mdcy(pycomp( -15 ), 1)=0
c     set stable lighter b-flavoured states: needed to analize single-top
c     events in a reasonable simple way (see the analize subroutine), but,
c     strictly speaking, not necessary for the program.
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

c     PYTHIA masses for heavy quarks
      pmas(pycomp(5),1)=4.75
      pmas(pycomp(-5),1)=4.75
      pmas(pycomp(6),1)=172.5
      pmas(pycomp(-6),1)=172.5

      pmas(pycomp(24),1)=80.398
      pmas(pycomp(-24),1)=80.398
c$$$      pmas(pycomp(24),2)=0.
c$$$      pmas(pycomp(-24),2)=0.

c$$$c     PYTHIA ckm matrix: for consistency it should be the                       
c$$$c     same as the one used in POWHEG                                            
c$$$      VCKM(1,1)=powheginput('CKM_ud')
c$$$      VCKM(1,2)=powheginput('CKM_us')
c$$$      VCKM(1,3)=powheginput('CKM_ub')
c$$$      VCKM(2,1)=powheginput('CKM_cd')
c$$$      VCKM(2,2)=powheginput('CKM_cs')
c$$$      VCKM(2,3)=powheginput('CKM_cb')
c$$$      VCKM(3,1)=powheginput('CKM_td')
c$$$      VCKM(3,2)=powheginput('CKM_ts')
c$$$      VCKM(3,3)=powheginput('CKM_tb')

c     top decay:
c     relevant only when POWHEG is run with spin correlations switched off;
c     in that case, force t->e ve b decay.
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

