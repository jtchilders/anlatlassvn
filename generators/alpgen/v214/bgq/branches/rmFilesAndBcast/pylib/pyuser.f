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


C     Set Perugia 2011 Tune
      call pygive('MSTP(5)=356')

C     Some custom masses
      call pygive('PMAS(6,1)=172.5')
      call pygive('PMAS(24,1)=80.399')
      call pygive('PMAS(23,1)=91.1876')
      call pygive('PARU(102)=0.23133')
      call pygive('PARJ(90)=20000.')
      call pygive('MDCY(15,1)=0')
      call pygive('MSTP(143)=1')


      
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
      print *,"Event # ",IEV

      CALL PYEVNT
      print *,"Done with gen event"
      CALL PYHEPC(1)
      CALL PYLIST(5)
      IF(NEV.NE.0.AND.IEV.EQ.NEV) GOTO 140
C NUP=0: REACHED END OF UNWEIGHTED EVENT FILE
      IF(NUP.EQ.0) GOTO 140
      IF(IEV.LE.2) THEN
C DEBUGGING:
c         CALL PYLIST(7)
c         CALL PYLIST(2)
c         CALL PYLIST(5)
      ENDIF
      IF(IEV.EQ.100) GOTO 140
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
   

