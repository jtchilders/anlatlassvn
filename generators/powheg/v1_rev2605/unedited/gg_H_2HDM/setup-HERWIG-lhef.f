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
      integer hdecaymode

      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)

      hdecaymode=lprup(1)-10000
      if ((hdecaymode.lt.-1).or.(hdecaymode.gt.12)) then
         write(*,*) "Higgs decay channel not allowed"
         stop
      endif
      if (hdecaymode.eq.-1) then
         iproc=-1600
      else
         iproc=-1600-hdecaymode ! Les Houches interface
      endif

C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      CALL HWUSTA('TAU-    ')
      CALL HWUSTA('TAU+    ')
      END

      subroutine UPEVNT
      implicit none
      integer hdecaymode
      include 'LesHouches.h'
      call lhefreadev(97)
      hdecaymode=lprup(1)-10000
      if (hdecaymode.eq.-1) then
         istup(3)=2    ! needed in order not to decay the Higss boson 
      endif
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
      integer hdecaymode
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
      hdecaymode=lprup(1)-10000
      if (hdecaymode.gt.0) then
      xwgtup=xwgtup*brhig(hdecaymode)
      endif
      call analysis(xwgtup)
      call pwhgaccumup 
      end

