      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.false.)
      logical vetotrunc
      parameter (vetotrunc=.false.)
ccccccccccccccccccccccccccccccccccccccccc
      logical asmcnlo
      parameter (asmcnlo=.true.)

      INTEGER N,NSTEP,I,JPR0,JPR
C QQIN IS THE EVENT FILE
      CHARACTER*50 QQIN
      COMMON/VVJIN/QQIN
      REAL*8 TMPLAM,GAMT0
      INTEGER IPDF
      CHARACTER * 70 LHAPDF
      LOGICAL LHACRTL,OLDFORM
      PARAMETER (LHACRTL=.TRUE.)
      LOGICAL ENDOFRUN
      COMMON/CENDOFRUN/ENDOFRUN
cccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Useful settings to interface POWHEG with HERWIG

      if(vetotrunc) then
         TRUNSH = .true.
         PTVETO = .true.
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     to use external PDF
c$$$      MODPDF(1)=10050
c$$$      MODPDF(2)=10050
c$$$      AUTPDF(1)='HWLHAPDF'
c$$$      AUTPDF(2)='HWLHAPDF'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     do not use soft me correction
      SOFTME=.FALSE.
      HARDME=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     tolerate 2% of killed events (default is 1%)
c      MAXER=MAXEV/50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     underlying event
      if(.not.uevent) PRSOF=0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!
c     intrinsic pt
      ptrms=2.5d0
      write(*,*) '**************************'
      write(*,*) 'Initial pt-spreading=',ptrms,' GeV'
      write(*,*) '**************************'
!!!!!!!!!!!!!!!!!!!!!!!!!!

c     number of events printed on the shell
      maxpr=2

c     do not print vertexes
      prvtx=.false.

cccccccccccccccccccccccccccccccccccccccccccc
c     to do exactly the same as in mc@nlo.
c     Here I will overwrite some parameters already set above, 
c     but I want to re-set all as in mc@nlo.
      if(asmcnlo) then

         PRESPL=.FALSE.

         PTRMS=2.5D0
         
c     UE already switched off (see above)

C Select W/Z boson decay modes
         MODBOS(1)=5
         MODBOS(2)=5

      WRITE(*,*)'Enter Lambda_QCD, <0 for Herwig default'
      TMPLAM=-1
      IF(TMPLAM.GE.0.D0)QCDLAM=TMPLAM
C
      WRITE(*,*)'Enter Z mass, width'
      RMASS(200)=91.17
      GAMZ=2.495
      WRITE(*,*)'Enter W mass, width'
      RMASS(198)=80.398
      GAMW=0.0
      RMASS(199)=RMASS(198)
      WRITE(*,*)'Enter top mass, width'
      RMASS(6)=172.5
      GAMT0=0.0  
      WRITE(*,*)'Enter Higgs (SM) boson mass, width'
      RMASS(201)=120.0
      GAMH=0.0049
      WRITE(*,*)'Enter quark (d,u,s,c,b) and gluon masses'
      RMASS(1)= 0.32
      RMASS(2)= 0.32
      RMASS(3)= 0.5
      RMASS(4)= 1.55
      RMASS(5)= 4.75
      RMASS(13)= 0.75

      DO I=1,5
         RMASS(I+6)=RMASS(I)
      ENDDO
C Set electron and muon masses equal to zero to avoid rounding problems
      RMASS(121)=0.D0
      RMASS(123)=0.D0
      RMASS(127)=0.D0
      RMASS(129)=0.D0
C NO SOFT AND HARD ME CORRECTIONS (ALREADY INCLUDED IN MC@NLO)
      SOFTME=.FALSE.
      HARDME=.FALSE.
      ZMXISR=0                  ! No photon radiation from ISR
      NOWGT=.FALSE.
C     NEGATIVE WEIGHTS ALLOWED
      NEGWTS=.TRUE.
      MAXPR=2
      MAXER=MAXEV/100
      LRSUD=0
      LWSUD=77
C     IN THE CASE HERWIG PDFS ARE USED, ADOPT MRST
      NSTRU=8
      PRVTX=.FALSE.
      PTMIN=0.5
      NRN(1)=1973774260         !ER: random number seed
      NRN(2)=1099242306         !ER: random number seed
C     THE FOLLOWING SHOULD BE USED ONLY IN WEIGHTED MODE
      IF(.NOT.NOWGT)THEN
         WGTMAX=1.000001D0
         AVABW=1.000001D0
      ENDIF
C FOR TOP PRODUCTION (HARMLESS ELSEWHERE)
      RLTIM(6)=1.D-23 
      RLTIM(12)=1.D-23
      endif
cccccccccccccccccccccccccccccccccccccccccccc

      end

      subroutine getmaxev(maxev)
      integer maxev
C---  Opens input file and counts number of events, setting MAXEV
      call opencount(maxev)
      end
      
      subroutine UPINIT
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      CALL HWUSTA('MU+     ')
      CALL HWUSTA('MU-     ')
      CALL HWUSTA('TAU-    ')
      CALL HWUSTA('TAU+    ')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     the following HWUSTA calls are exactly as in MC@NLO.
c     They are needed to analize single-top events in a reasonable simple way
c     (see the analize subroutine), but, strictly speaking, they are not necessary
c     for the program.
      CALL HWUSTA('B+      ')!
      CALL HWUSTA('B-      ')!
      CALL HWUSTA('B_D0    ')!
      CALL HWUSTA('B_DBAR0 ')!
      CALL HWUSTA('B_S0    ')!
      CALL HWUSTA('B_SBAR0 ')!
      CALL HWUSTA('SIGMA_B+')!
      CALL HWUSTA('SIGMA_B-')!
      CALL HWUSTA('XI_B0   ')!
      CALL HWUSTA('XI_B+   ')!
      CALL HWUSTA('XI_B-   ')!
      CALL HWUSTA('B_C+    ')!
      CALL HWUSTA('B_C-    ')!
      CALL HWUSTA('UPSLON1S')!
      CALL HWUSTA('SGM_BBR+')!
      CALL HWUSTA('SGM_BBR-')!
      CALL HWUSTA('LMD_BBR0')!
      CALL HWUSTA('OMEGA_B-')!
      CALL HWUSTA('XI_BBAR0')!
      CALL HWUSTA('OMG_BBR+')!
      CALL HWUSTA('LMBDA_B0')!
c     call HWIODK to see the HERWIG decay table

c     HERWIG masses for heavy quarks
c     (antiparticles with i+6)
      RMASS(5) = 4.75
      RMASS(11)= 4.75
      RMASS(6) = 172.5
      RMASS(12)= 172.5

      RMASS(198)=80.398
      RMASS(199)=RMASS(198)
      GAMW=0d0

ccccccccccccccccccccccccccccccccccccccccc
c$$$c     HERWIG ckm matrix: for consistency it should be the                       
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
cccccccccccccccccccccccccccccccccccccccccccc
      end
      

      subroutine UPEVNT
      implicit none
c     Needed to force t->e ve b decay.
      CALL HWMODK(6,1d0,10005,12,-11,5,0,0)
      CALL HWMODK(-6,1d0,10005,-12,11,-5,0,0)
      call lhefreadev(97)
c      call HWIODK(66,1,0)
c      stop
      end


      subroutine hwabeg
      call init_hist
      end

      subroutine hwaend
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      open(unit=99,file=pwgprefix(1:lprefix)//'POWHEG+HERWIG-output.top'
     #     ,status='unknown')
      call pwhgsetout
      call pwhgtopout
      close(99)
      end
      

      subroutine hwanal
      INCLUDE 'HERWIG65.INC'
      include 'LesHouches.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)

      if (ierror.ne.0) then
         if(verbose) then
            write(*,*) 'Killed event'
            write(*,*) 'Scalup= ',scalup
            call HWUPUP         !hepeup
            call hwuepr         !all the event
         endif
         return
      endif
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

