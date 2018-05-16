      subroutine setup_HERWIG_parameters
      include 'HERWIG65.INC'
      include 'LesHouches.h'
      logical uevent 
      parameter (uevent=.false.)
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
c      PTMIN=100.
      ptrms=2.5d0
      print *
      print *,'* UPINIT: Intrinsic pT spread set to = ',ptrms,' GeV *'
      print *
      maxpr=2
      if(.not.uevent) PRSOF=0d0
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
      logical uevent 
      parameter (uevent=.false.)
      nevhep=0
c first call lhefreadhdr; this sets lprup;
      call lhefreadhdr(97)

C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')
      CALL HWUSTA('TAU-    ')
      CALL HWUSTA('TAU-    ')
      end

      subroutine UPEVNT

      INCLUDE 'HERWIG65.INC'

      include 'LesHouches.h'

      LOGICAL MASSIVE_FS_PARTONS,DEBUGGING
      INTEGER IXX,JXX,NFS
      REAL*8  PARTON_MOM(4,5)
      REAL*8  PTOT(4),PIN(5),POUT(5),MPOUT(5)
      REAL*8  MODP2(5),MASS2(5),SHAT
      REAL*8  SCALING,THRESHOLD,TMP1,TMP2,NUM,DENOM,DIFF

      logical ini
      save ini
      data ini/.true./
      call lhefreadev(97)

      if (ini) then
         ini = .false.
         write(*,*)
         write(*,*) '*******************************'
         write(*,*) '*******************************'
         write(*,*) '*                             *'
         write(*,*) '*  setup-HERWIG-lhef: UPEVNT  *'
         write(*,*) '*                             *'
         write(*,*) '*   Initialization finished   *'
         write(*,*) '*                             *'
         write(*,*) '*******************************'
         write(*,*) '*******************************'
         write(*,*)
      endif


C - If there are more than 4 partons i.e. 2 or more final state
C - partons and if we are going to be doing rescaling (option
C - MASSIVE_FS_PARTONS=.TRUE.) then we
C - start here:

C - 09/03/11 - KH - This old reshuffling code putting quarks on-shell
C - has been disabled here by setting MASSIVE_FS_PARTONS=.FALSE. .
C - Its job is now being done by the call to the "blessed" lhefinitemasses
C - routine in Born.f. The code is kept here for historical value only
C - (you never know ...) and can be removed.
C -      MASSIVE_FS_PARTONS=.TRUE.
      MASSIVE_FS_PARTONS=.FALSE.

      DEBUGGING=.FALSE.           ! Suppress debugging output
      IF(NUP.GE.4.AND.MASSIVE_FS_PARTONS.EQV..TRUE.) THEN

C - Initialise incoming, outgoing and total momentum vectors
         DO IXX=1,4
            PTOT(IXX)=0d0
            PIN(IXX) =0d0
            POUT(IXX)=0d0
         ENDDO

C - Calculate initial incoming, outgoing and total momentum in the lab
         DO IXX=1,NUP
            IF(ISTUP(IXX).EQ.-1) THEN
               DO JXX=1,4
                  PTOT(JXX)=PTOT(JXX)+PUP(JXX,IXX)
                  PIN(JXX) =PIN(JXX) +PUP(JXX,IXX)
               ENDDO
            ELSE
               DO JXX=1,4
                  PTOT(JXX)=PTOT(JXX)-PUP(JXX,IXX)
                  POUT(JXX)=POUT(JXX)+PUP(JXX,IXX)
               ENDDO
            ENDIF
         ENDDO

C - Calculate the partonic COM energy
         SHAT=POUT(4)*POUT(4)-POUT(1)*POUT(1)
     $       -POUT(2)*POUT(2)-POUT(3)*POUT(3)
         PIN(5) =SQRT(SHAT)
         POUT(5)=SQRT(SHAT)

C - Show debugging output
         IF(DEBUGGING.OR.(ABS(PTOT(1)).GT.1D-5).OR.
     $                   (ABS(PTOT(2)).GT.1D-5).OR.
     $                   (ABS(PTOT(3)).GT.1D-5).OR.
     $                   (ABS(PTOT(4)).GT.1D-5)) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'UPEVNT: pre rescaling debugging output'
            WRITE(6,*) 'PTOT = ',PTOT(1),PTOT(2),PTOT(3),PTOT(4)
            WRITE(6,*) 'PIN  = ',PIN(1) ,PIN(2) ,PIN(3) ,PIN(4) ,PIN(5)
            WRITE(6,*) 'POUT = ',POUT(1),POUT(2),POUT(3),POUT(4),POUT(5)
            DO IXX=1,NUP
               WRITE(6,*) 
     $           IXX,IDUP(IXX),
     $           PUP(1,IXX),PUP(2,IXX),PUP(3,IXX),PUP(4,IXX),PUP(5,IXX),
     $           RMASS(ABS(IDUP(IXX)))
            ENDDO
         ENDIF

C - Initialise the vectors of |p|, masses and the counters for the
C - threshold energy and number of final state particles 
         NFS=0
         THRESHOLD=0
         DO IXX=1,NUP
            MODP2(IXX)=0
            MASS2(IXX)=0
         ENDDO

C - Work out the momenta and their magnitudes^2, & the HW masses^2 in
C - the partonic COM frame: PARTON_MOM(JXX,IXX), MODP2(IXX),MASS2(IXX)
         DO IXX=1,NUP
            CALL HWULF4(POUT,PUP(1,IXX),PARTON_MOM(1,IXX))
            IF(ISTUP(IXX).NE.-1) THEN
               NFS=NFS+1
               THRESHOLD=THRESHOLD+RMASS(ABS(IDUP(IXX)))
            ENDIF
            DO JXX=1,3
               MODP2(IXX)=MODP2(IXX)+PARTON_MOM(JXX,IXX)**2
            ENDDO
            MASS2(IXX)=RMASS(ABS(IDUP(IXX)))**2
         ENDDO

C - If the threshold condition is met then we attempt to solve for the
C - final state momentum rescaling factor using the Newton Raphson method
         SCALING=1D0               ! Initialisation of rescaling factor
         IF(SQRT(SHAT).GE.THRESHOLD) THEN
            TMP1=1D0 ! Initial guess solution (no rescaling: rescaling = 1)
            DIFF=1D0 ! Initialisation of the convergence measure
            DO WHILE(ABS(DIFF/TMP1).GT.1D-10)
               NUM=-SQRT(SHAT)     ! Initialisation of the root equation
               DENOM=0D0           ! Initialisation of its derivative 
               DO IXX=1,NUP
                  IF(ISTUP(IXX).NE.-1) THEN
                     NUM=NUM+SQRT(TMP1**2*MODP2(IXX)+MASS2(IXX))
                     DENOM=DENOM+TMP1*MODP2(IXX)
     $                          /SQRT(TMP1**2*MODP2(IXX)+MASS2(IXX))
                  ENDIF
               ENDDO
               TMP2=TMP1-NUM/DENOM ! New best guess
               DIFF=TMP2-TMP1      ! To test convergence   
               TMP1=TMP2           ! Next trial value equals new best guess
               IF(TMP1.LE.0) THEN
                  WRITE(6,*) 'main-HERWIG-lhef UPEVNT warning:'
                  WRITE(6,*) 'Negative value encountered in solving'
                  WRITE(6,*) 'for the momentum rescaling factor.'
                  WRITE(6,*) 'Momenta will not be rescaled to the HW'
                  WRITE(6,*) 'mass shell values'
                  GOTO 300
               ENDIF
            ENDDO
            SCALING=TMP1
         ENDIF
C - Rescale the momenta
         DO IXX=1,NUP
            IF(ISTUP(IXX).NE.-1) THEN
               DO JXX=1,3
                  PARTON_MOM(JXX,IXX)=SCALING*PARTON_MOM(JXX,IXX)
               ENDDO
               PARTON_MOM(4,IXX)=
     $              SQRT(PARTON_MOM(1,IXX)*PARTON_MOM(1,IXX)
     $                  +PARTON_MOM(2,IXX)*PARTON_MOM(2,IXX)
     $                  +PARTON_MOM(3,IXX)*PARTON_MOM(3,IXX)
     $                  +MASS2(IXX))
            ENDIF
         ENDDO
C - Boost everything back to the lab frame
 300     DO JXX=1,3
            MPOUT(JXX)=-POUT(JXX)
         ENDDO
         MPOUT(4)=POUT(4)
         MPOUT(5)=POUT(5)
         DO IXX=1,NUP
            CALL HWULF4(MPOUT,PARTON_MOM(1,IXX),PUP(1,IXX))
C - Calculate the mass components as E^2-|p|^2
            PUP(5,IXX)=(PUP(4,IXX)+PUP(3,IXX))*(PUP(4,IXX)-PUP(3,IXX))
     $                - PUP(1,IXX)*PUP(1,IXX) - PUP(2,IXX)*PUP(2,IXX)
C - Sometimes the near massless particles can experience rounding
C - errors in the calculation of their mass components, almost always
C - the initial state ones. This is dealt with here:
            IF(PUP(5,IXX).LT.0D0) THEN
               IF(DEBUGGING.OR.PUP(5,IXX).LT.-1D-6) THEN
               WRITE(6,*) 'main-HERWIG-lhef: UPEVNT warning'
               WRITE(6,*) 'PUP(JXX,',IXX,') has E^2-|p|^2 = ',PUP(5,IXX)
               WRITE(6,*) 'Setting PUP(JXX,IXX) = -SQRT(',PUP(5,IXX),')'
               ENDIF
               PUP(5,IXX)=SQRT(-PUP(5,IXX))
            ELSE
               PUP(5,IXX)=SQRT( PUP(5,IXX))
            ENDIF
         ENDDO

C - Initialise incoming, outgoing and total momentum vectors
         DO IXX=1,4
            PTOT(IXX)=0d0
            PIN(IXX) =0d0
            POUT(IXX)=0d0
         ENDDO

C - Calculate initial incoming, outgoing and total momentum in the lab
         DO IXX=1,NUP
            IF(ISTUP(IXX).EQ.-1) THEN
               DO JXX=1,4
                  PTOT(JXX)=PTOT(JXX)+PUP(JXX,IXX)
                  PIN(JXX) =PIN(JXX) +PUP(JXX,IXX)
               ENDDO
            ELSE
               DO JXX=1,4
                  PTOT(JXX)=PTOT(JXX)-PUP(JXX,IXX)
                  POUT(JXX)=POUT(JXX)+PUP(JXX,IXX)
               ENDDO
            ENDIF
         ENDDO

C - Recalculate the partonic COM energy
         SHAT=POUT(4)*POUT(4)-POUT(1)*POUT(1)
     $       -POUT(2)*POUT(2)-POUT(3)*POUT(3)
         PIN(5) =SQRT(SHAT)
         POUT(5)=SQRT(SHAT)

C - Show debugging output
         IF(DEBUGGING.OR.(ABS(PTOT(1)).GT.1D-5).OR.
     $                   (ABS(PTOT(2)).GT.1D-5).OR.
     $                   (ABS(PTOT(3)).GT.1D-5).OR.
     $                   (ABS(PTOT(4)).GT.1D-5)) THEN
            WRITE(6,*) 'UPEVNT: post rescaling debugging output'
            WRITE(6,*) 'PTOT = ',PTOT(1),PTOT(2),PTOT(3),PTOT(4)
            WRITE(6,*) 'PIN  = ',PIN(1) ,PIN(2) ,PIN(3) ,PIN(4) ,PIN(5)
            WRITE(6,*) 'POUT = ',POUT(1),POUT(2),POUT(3),POUT(4),POUT(5)
            DO IXX=1,NUP
               WRITE(6,*) 
     $           IXX,IDUP(IXX),
     $           PUP(1,IXX),PUP(2,IXX),PUP(3,IXX),PUP(4,IXX),PUP(5,IXX),
     $           RMASS(ABS(IDUP(IXX)))
            ENDDO
         ENDIF

      ENDIF

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
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end
