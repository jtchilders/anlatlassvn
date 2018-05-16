      SUBROUTINE ALSHCD(CSHO)
      CHARACTER*3 CSHO
      CSHO='HER'
      END

      SUBROUTINE ALSHEN
      INCLUDE 'HERWIG65.INC'
      WGTSUM=WGTSUM-evwgt
      WSQSUM=WSQSUM-evwgt**2
      END

      SUBROUTINE ALSHER(I)
      INCLUDE 'HERWIG65.INC'
      INTEGER I
      I=0
      MAXER=MAX(10,NEVHEP/100)
      IF (IERROR.NE.0.OR.EVWGT.EQ.0) THEN
        I=1
      ENDIF 
      RETURN
      END
CDECK  ID>, UPINIT.
*CMZ :-        -13/02/02  07.20.46  by  Peter Richardson
*-- Author :    Michelangelo Mangano
C----------------------------------------------------------------------
      SUBROUTINE ALSHIN(I)
C----------------------------------------------------------------------
C     subroutine to initialise the events
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'alpsho.inc'
      INTEGER I
C COMMUNICATION CODE
      I=1
C FILL IN HERWIG SPECIFIC INPUTS
      IF(MC.NE.0) RMASS(4)=MC
      IF(MB.NE.0) RMASS(5)=MB
      IF(MT.NE.0) RMASS(6)=MT
      IF(MW.NE.0) THEN
        RMASS(198)=MW
        RMASS(199)=MW
      ENDIF
      IF(MZ.NE.0) RMASS(200)=MZ
      IF(MH.NE.0) RMASS(201)=MH
cc
      END

C-----------------------------------------------------------------------
      SUBROUTINE ALSFIN
C-----------------------------------------------------------------------
C     HERWIG END OF FILE TREATMENT
C-----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      include 'alpsho.inc'
C--   GUP Run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &     IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &     XMAXUP(MAXPUP),LPRUP(MAXPUP)
      DOUBLE PRECISION XSEC,EFF
      XSEC=XSECUP(1)
      WGTMAX=XSEC*1.0D-3
      WBIGST=XSEC*1.0D-3
      WGTSUM=WGTSUM*REAL(NEVHEP-1)/REAL(NEVHEP)
      WSQSUM=WSQSUM*REAL(NEVHEP-1)/REAL(NEVHEP)
      WGTSUM=WGTSUM*XSEC*1.0D-3
      WSQSUM=WSQSUM*XSEC**2*1.0D-6
C FILL HERWIG RUN LOG FILE
      WRITE(NUNITOUT,*) " "
      WRITE(NUNITOUT,*) "END OF RUN LOG: "
      WRITE(NUNITOUT,*) "NUMBER OF INPUT PROCESSED:",NWGTS
      WRITE(NUNITOUT,*)"NUMBER OF EVENTS GENERATED:",NEVHEP-1
      EFF=REAL(NEVHEP-1)/REAL(NWGTS)
      WRITE(NUNITOUT,*)"MATCHING EFFICIENCY:",EFF
      WRITE(NUNITOUT,*)"FINAL CROSS SECTION (PB):",XSEC*EFF
      WRITE(NUNITOUT,*)"INTEGRATED LUMINOSITY (PB-1):",REAL(NEVHEP-1)
     $     /(XSEC*EFF)
      CALL HWUGUP
      END
