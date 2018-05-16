      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
      MAXPR=2
c     do not use soft me correction
      SOFTME=.FALSE.     

c     tolerate 2% of killed events (default is 1%)
      MAXER=MAXEV/50

c     no underlying event
      PRSOF=0d0

c     intrinsic pt
      ptrms=2.5d0
      write(*,*) '**************************'
      write(*,*) 'Initial pt-spreading=',ptrms,' GeV'
      write(*,*) '**************************'

c     do not print vertexes
      prvtx=.false.

c     set stable lighter b-flavoured states: see the following calls to HWUSTA

c     top decay:
c     relevant only when POWHEG is run with spin correlations switched off;
c     in that case, force t->e ve b decay.
c     see the following calls to HWMODK
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
c      CALL HWUINC

C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')

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
      if (lprup(1).eq.10015) then
         CALL HWUSTA('TAU-    ')
         CALL HWUSTA('TAU+    ')
      endif
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

