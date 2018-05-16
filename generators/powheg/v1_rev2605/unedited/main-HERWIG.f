      PROGRAM HWIGPR
C---COMMON BLOCKS ARE INCLUDED AS FILE herwig6510.h
      INCLUDE 'HERWIG65.INC'
      include 'LesHouches.h'
      integer n
c     we need to tell to the analysis file which program is running it
      character * 6 WHCPRG
      integer iun
      common/cWHCPRG/WHCPRG
      WHCPRG='HERWIG'
C---PROCESS; set to negative for user supplied me
      iproc=-1                  ! Les Houches interface
C--- Opens input file and counts number of events, setting MAXEV;
c    MAXEV must be set before HWIGIN call.
      call getmaxev(maxev)
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---SETUP INITIAL PARAMETERS
      call setup_HERWIG_parameters
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---USER'S INITIAL CALCULATIONS
      CALL HWABEG
C---INITIALISE ELEMENTARY PROCESS
      CALL HWEINI
C---LOOP OVER EVENTS
      DO N=1,maxev
C---  INITIALISE EVENT
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
         if ((nevhep.gt.0).and.(mod(nevhep,20000).eq.0)) then
            write(*,*) "# of events processed =",nevhep
            call hwaend
         endif
      ENDDO
 111  continue
C---  TERMINATE ELEMENTARY PROCESS
      CALL HWEFIN
      write(*,*) 'At the end NEVHEP is ',nevhep
C---  USER'S TERMINAL CALCULATIONS
      CALL HWAEND

      END

