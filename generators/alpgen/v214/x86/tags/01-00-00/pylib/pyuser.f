*--   Author :   Andrea Messina
*--   modified by Stephen Mrenna
C-----------------------------------------------------------------------------
      PROGRAM PYUSER
C-----------------------------------------------------------------------------
      
c     INCLUDE 'pythia6300.inc'
      
C...EXTRA COMMONBLOCK TO TRANSFER RUN INFO.
      INTEGER NEV
      DATA NEV /0/
      COMMON/PRIV/MODE,NLIM  
      COMMON /PROC/ IPROC
C...USER PROCESS EVENT COMMON BLOCK.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &    ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &    VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE/HEPEUP/
      
C...SET NUMBER OF EVENTS
C      NEV=2500
C...SWITCH PROCESS MODE; AGREES WITH IDWTUP CODE (+-1,+-2,+-3,+-4).
      EXTERNAL PYDATA


      MODE=3

C INPUT HERE YOUR PARAMETER DEFAULTS:
C.....(1) it is always safer to use PYGIVE, since you don't
C.........have to include common blocks directly
C.....(2) using the string "VAR=" prints out the current value of VAR
C.....(3) PYGIVE gives a record in the log file of what was changed
C.....
C.....probably DON'T want to change these from defaults
C.....ISR on by default, set to 2 for QED and QCD, 1 for just QCD, 0 for none
      call pygive('MSTP(61)=')                ! INITIAL STATE RADIATION
C.....FSR on by default, set to 0 for none 
      call pygive('MSTP(71)=')                ! FINAL STATE RADIATION
C.....Multiple interactions on (NEW MODEL)
C.......Is this consistent with old parton shower?
      call pygive('MSTP(81)=')                ! PILE UP
C.....NOOOOOOO!!  This uses a toy model of underlying event
C.....      call pygive('MSTP(82)=0')                ! GAUSSIAN PILE UP
C.....on by default, set to 0 to turn off
      call pygive('MSTP(82)=')
      call pygive('MSTP(111)=')               ! HADRONIZATION
C
C NEW FOR ALPGEN VERSION V2.0 AND PYTHIA VERSION>=6.226 OR >=6.320
C CRUCIAL FOR JET-PARTON MATCHING:
      call pygive('MSTP(143)=1')       ! CALL UPVETO, ALLOW JET-PARTON MATCHING
C      CALL PYSTAT(2)            ! DECAY TABLES AND PARTICLES INFORMATIONS

C...SET PI0 STABLE TO TRIM EVENT LISTINGS.
C.....NOT GOOD FOR FULL EVENT SIMULATIONS!!!!!!!!!
      call pygive('MDCY(C111,1)=0')
      
      
C...INITIALIZE WITH EXTERNAL PROCESS.
      CALL PYINIT('USER',' ',' ',0D0)
C INITIALIZE BOOKING OF USER HTSOGRAMS, ANALYSIS
c      CALL INITHBOOK

C   SELECT NEV>0 TO SET AN UPPER LIMIT TO THE NUMBER OF GENERATED EVENTS
      NEV=0     ! RUN THROUGH ALL EVENTS OF THE ALPGEN FILE
C      NEV=100   ! STOP AFTER 100 EVENTS

C...EVENT LOOP; GENERATE EVENT; CHECK IT WAS OBTAINED OR QUIT.
      IEV=0
 1    IEV=IEV+1

      CALL PYEVNT
      CALL PYHEPC(1)
      IF(NEV.NE.0.AND.IEV.EQ.NEV) GOTO 140
C NUP=0: REACHED END OF UNWEIGHTED EVENT FILE
      IF(NUP.EQ.0) GOTO 140
      IF(IEV.LE.2) THEN
C DEBUGGING:
c         CALL PYLIST(7)
c         CALL PYLIST(2)
c         CALL PYLIST(5)
      ENDIF
C
C CALL HERE YOUR ANALYSIS ROUTINES:
c      IF(IEV.le.1000) CALL FILLNTUPLE(IPROC)
c      IF (MOD(IEV,1000).EQ.0)   PRINT *, IEV
 130  GOTO 1
      
C...  STATISTICS AND HISTOGRAMS.
 140  CALL PYSTAT(1)
C
C TERMINATE YOUR ANALYSIS:
c      CALL CLOSENTUPLE
      
      END
   

