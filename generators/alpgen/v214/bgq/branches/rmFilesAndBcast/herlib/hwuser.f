      program hwigpr
      include 'HERWIG65.INC'
      EXTERNAL HWUDAT
      integer icount,ipveto
c---process - don't change
      iproc=-100
c---initialise other common blocks
      call hwigin
c---user can reset parameters at
c   this point, otherwise values
c   set in hwigin will be used.
      maxpr=1
      prvtx = .false.    ! to exclude the  production vtx in printout
c remove
c      call hwigup
c---compute parameter-dependent constants
      call hwuinc
c---call hwusta to make any particle stable
      call hwusta('PI0     ')
c---user's initial calculations       
      call hwabeg
c---initialise elementary process
      call hweini
c---loop over events
      icount=0
c---initialise event
 1    CALL hwuine
c      icount=icount+1
c      if(mod(icount,1000).eq.0)write(*,*) 'processed',icount,' events'
c---generate hard subprocess
      call hwepro
c matching-driven veto
      CALL UPVETO(IPVETO)
      IF(IPVETO.NE.0) THEN
        CALL HWWARN('UPVETO',-1)
        GOTO 100
      ENDIF
c---generate parton cascades
      call hwbgen
c---do heavy objects deycays
      call hwdhob
c---do cluster hadronization
c      call hwcfor
c---do cluster decay
c      call hwcdec
c---do unstable particle decays
      call hwdhad
c---do heavy flavour decays
      call hwdhvy
c---  add soft underlying event if needed 
c      call hwmevt
c-- event generation completed, wrap up event .... 
 100  CALL HWUFNE
c-- and carry out user analysis
      CALL HWANAL
c      if(icount.gt.10000) call hwugup
      GOTO 1
  999 CONTINUE
      END

C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
  999 END                                         
C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C user analysis routine
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      IF(IERROR.NE.0) RETURN
  999 END  
C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C  user routine to output analysis results
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
 999  END

      SUBROUTINE TIMEL(TRES)
C---DUMMY SUBROUTINE FOR CPU TIME REMAINING.  DELETE AND
C   LINK TO CERN PROGRAM LIBRARY TO GET CORRECT VALUE
      TRES=1000
      END
cC-----------------------------------------------------------------------
cCDECK  ID>,  TIMEL.
c*CMZ :-        -28/06/01  16.55.32  by  Bryan Webber
c*-- Author :    Bryan Webber
cC-----------------------------------------------------------------------
c      SUBROUTINE TIMEL(TRES)
cC-----------------------------------------------------------------------
cC     DUMMY TIME SUBROUTINE: DELETE AND REPLACE BY SYSTEM
cC     ROUTINE GIVING TRES = CPU TIME REMAINING (SECONDS)
cC-----------------------------------------------------------------------
c      REAL TRES
c      LOGICAL FIRST
c      DATA FIRST/.TRUE./
c      SAVE FIRST
c      IF (FIRST) THEN
c      WRITE (6,10)
c   10 FORMAT(/10X,'SUBROUTINE TIMEL CALLED BUT NOT LINKED.'/
c     &        10X,'DUMMY TIMEL WILL BE USED. DELETE DUMMY'/
c     &        10X,'AND LINK CERNLIB FOR CPU TIME REMAINING.')
c      FIRST=.FALSE.
c      ENDIF
c      TRES=1E10
c      END
