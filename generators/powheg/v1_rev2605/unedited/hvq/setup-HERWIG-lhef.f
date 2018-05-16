      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.true.)
      if(.not.uevent) PRSOF=0d0
c----DO NOT USE SOFT ME CORRECTION     
      SOFTME=.FALSE.
      end


      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

      subroutine UPINIT
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)
      if(lprup(1).eq.1004) then
         call hwusta('D+      ')
         call hwusta('D-      ')
         call hwusta('D0      ')
         call hwusta('DBAR0   ')
      elseif(lprup(1).eq.1005) then
         call hwusta('B_D0    ')
         call hwusta('B_DBAR0 ')
         call hwusta('B-      ')
         call hwusta('B+      ')
      endif
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      end

      subroutine UPEVNT
      implicit none
      include 'LesHouches.h'
      logical ini
      save ini
      data ini/.true./
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
      if (ierror.ne.0) then
         return
      endif
      nevhep=nevhep+1
      if(idwtup.eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end
      
