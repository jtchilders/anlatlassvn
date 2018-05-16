      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.true.)
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c      PTMIN=100.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Useful settings to interface POWHEG with HERWIG

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     tolerate 2% of killed events (default is 1%)
      MAXER=MAXEV/50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     underlying event
      if(.not.uevent) PRSOF=0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      CALL HWUSTA('TAU-    ')
      CALL HWUSTA('TAU+    ')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     the following HWUSTA calls are exactly as in MC@NLO.
c     They are needed to analize single-top events in a reasonable simple way
c     (see the analize subroutine), but, strictly speaking, they are not necessary
c     for the program.
      CALL HWUSTA('B+      ')
      CALL HWUSTA('B-      ')
      CALL HWUSTA('B_D0    ')
      CALL HWUSTA('B_DBAR0 ')
      CALL HWUSTA('B_S0    ')
      CALL HWUSTA('B_SBAR0 ')
      CALL HWUSTA('SIGMA_B+')
      CALL HWUSTA('SIGMA_B-')
      CALL HWUSTA('XI_B0   ')
      CALL HWUSTA('XI_B+   ')
      CALL HWUSTA('XI_B-   ')
      CALL HWUSTA('B_C+    ')
      CALL HWUSTA('B_C-    ')
      CALL HWUSTA('UPSLON1S')
      CALL HWUSTA('SGM_BBR+')
      CALL HWUSTA('SGM_BBR-')
      CALL HWUSTA('LMD_BBR0')
      CALL HWUSTA('OMEGA_B-')
      CALL HWUSTA('XI_BBAR0')
      CALL HWUSTA('OMG_BBR+')
      CALL HWUSTA('LMBDA_B0')

c     Needed to force t->e ve b decay.
      CALL HWMODK(6,ONE,100,12,-11,5,0,0)
      CALL HWMODK(-6,ONE,100,-12,11,-5,0,0)
c     call HWIODK to see the HERWIG decay table
      end
      

      subroutine UPEVNT
      implicit none
      call lhefreadev(97)
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

