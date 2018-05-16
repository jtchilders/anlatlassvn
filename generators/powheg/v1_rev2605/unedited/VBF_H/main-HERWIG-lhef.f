      PROGRAM HWIGPR
C---COMMON BLOCKS ARE INCLUDED AS FILE herwig6510.h
      INCLUDE '../herwig6510.inc'
      include 'LesHouches.h'
      include 'PhysPars.h'
      integer n
      logical uevent 
      parameter (uevent=.false.)
c     we need to tell to the analysis file which program is running it
      character * 6 WHCPRG
      integer iun
      common/cWHCPRG/WHCPRG
      WHCPRG='HERWIG'
C---PROCESS; set to negative for user supplied me
      iproc=-1                  ! Les Houches interface
C--- Opens input file and counts number of events, setting MAXEV;
c    MAXEV must be set before HWIGIN call.
      call opencount(maxev)
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c     set Higgs boson mass and width
      call init_couplings
c      PTMIN=100.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ptrms=2.5d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      maxpr=2
      if(.not.uevent) PRSOF=0d0
c----DO NOT USE SOFT ME CORRECTION     
      SOFTME=.FALSE.
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      CALL HWUSTA('HIGGS   ')
C---USER'S INITIAL CALCULATIONS
      CALL HWABEG
C---INITIALISE ELEMENTARY PROCESS
      CALL HWEINI
C---LOOP OVER EVENTS
      DO N=1,maxev
C---INITIALISE EVENT
      CALL HWUINE
C---GENERATE HARD SUBPROCESS
      CALL HWEPRO
      if(nup.eq.0) goto 111
C---GENERATE PARTON CASCADES
      CALL HWBGEN
C---DO HEAVY OBJECT DECAYS
      CALL HWDHOB
C---DO CLUSTER FORMATION
      CALL HWCFOR
C---DO CLUSTER DECAYS
      CALL HWCDEC
C---DO UNSTABLE PARTICLE DECAYS
      CALL HWDHAD
C---DO HEAVY FLAVOUR HADRON DECAYS
      CALL HWDHVY
C---ADD SOFT UNDERLYING EVENT IF NEEDED
      CALL HWMEVT
C---FINISH EVENT
      CALL HWUFNE    
C---USER'S EVENT ANALYSIS
      CALL HWANAL
      if (mod(nevhep,20000).eq.0) then
         write(*,*) "# of events processed =",nevhep
         call hwaend
      endif
      ENDDO
 111  continue
C---TERMINATE ELEMENTARY PROCESS
      CALL HWEFIN
C---USER'S TERMINAL CALCULATIONS
      CALL HWAEND

c      call newunit(iun)
c      open(unit=iun,file='HERWIGcounters.dat',status='unknown')
c      call printcnt(iun)
c      close(iun)

      END

      subroutine UPINIT
      implicit none
      call lhefreadhdr(97)
      end

      subroutine UPEVNT
      implicit none
      include 'LesHouches.h'
      logical ini
      save ini
      data ini/.true./
      call lhefreadev(97)

      if (ini) then
         ini = .false.
         write(*,*)
         write(*,*)
         write(*,*) '**************************************'
         write(*,*) '**************************************'
         write(*,*) '***         NO HIGGS DECAY         ***'
         write(*,*) '**************************************'
         write(*,*) '**************************************'
         write(*,*)
         write(*,*)
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCCCCCCC        HIGGS BOSON ANALYSIS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      istup(3)=2                ! needed in order not to decay the Higss boson 
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
      INCLUDE '../herwig6510.inc'
      include 'LesHouches.h'
      if (ierror.ne.0) then
         return
      endif
      if(idwtup.eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end
