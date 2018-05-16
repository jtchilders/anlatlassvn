      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.true.)
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c      PTMIN=100.
      MAXPR=2

      ptrms=2.5d0
      write(*,*)
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*) ' INITIAL p_T SPREADING OF ',ptrms,' GEV    '
      write(*,*) '*******************************************'
      write(*,*) '*******************************************'
      write(*,*)
      if(.not.uevent) then
         PRSOF=0d0
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' NO UNDERLYING EVENT  WILL BE GENERATED    '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      else
         write(*,*)
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*) ' UNDERLYING EVENT  WILL BE GENERATED     '
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
         write(*,*)   
      endif
c----DO NOT USE SOFT ME CORRECTION     
      SOFTME=.FALSE.
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

