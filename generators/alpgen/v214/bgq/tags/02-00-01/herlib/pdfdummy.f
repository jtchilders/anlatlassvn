*CMZ :-        -26/04/91  11.11.54  by  Bryan Webber
*-- Author :    Bryan Webber
C----------------------------------------------------------------------
      SUBROUTINE PDFSET(PARM,VAL)
C----------------------------------------------------------------------
C     DUMMY SUBROUTINE: DELETE AND SET MODPDF(I)
C     IN MAIN PROGRAM IF YOU USE PDFLIB CERN-LIBRARY
C     PACKAGE FOR NUCLEON STRUCTURE FUNCTIONS
C----------------------------------------------------------------------
      DOUBLE PRECISION VAL(20)
      CHARACTER*20 PARM(20)
      WRITE (6,10)
   10 FORMAT(/10X,'PDFSET CALLED BUT NOT LINKED')
      STOP
      END
*CMZ :S        E26/04/91  11.11.54  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE STRUCTM(X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
C-----------------------------------------------------------------------
C     DUMMY SUBROUTINE: DELETE IF YOU USE PDFLIB CERN-LIBRARY
C     PACKAGE FOR NUCLEON STRUCTURE FUNCTIONS
C-----------------------------------------------------------------------
      DOUBLE PRECISION X,QSCA,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU
      WRITE (6,10)
  10  FORMAT(/10X,'STRUCTM CALLED BUT NOT LINKED')
      STOP
      END
