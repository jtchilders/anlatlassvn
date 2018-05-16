C**** DECEMBER 1990  :
C     RENAME TO MLMGEN **** WHOLE NEW SET OF PDF: DO,EHLQ,DFLM,HMRS,MT
C****  JULY 1990         -- ADD MULTIHISTOGRAMMING
C****  JAN 29 , 1990 -- USE WEBBER ROUTINE FOR STRUCTURE FUNCTIONS
C****  JUNE 12   , 1989  -- ADD SCATTER PLOT ROUTINES, AND INCLUDE SOME
C                           UTILITIES (VSUM, DOT)
C****  DECEMBER 2, 1988  -- ADD EHLQ1,GHR
C****  NOVEMBER 9, 1988
C
C
C****  MLMGEN.FOR ********
C--- GENERAL UTILITY SUBROUTINES:
C--- STRUCTURE FUNCTIONS, VEGAS, HISTOGRAMMING&PLOTTING


C-----------------------------------------------------------------------
C------- START STRUCTURE FUNCTION SECTION -------------------------------
C--------------------------------------------------------------------------


c-------------------------------------------------------------------
      subroutine pdfconv(nin,nout,type)
c-------------------------------------------------------------------
c converts ALPHA convention for PDF namings to hvqpdf conventions
c
c THIS ROUTINE SHOULD APPEAR AS WELL IN ALPSHO.F, NAMED AS PDFCONVH
C
      implicit none
      integer nin,nout,i
      character*25 type
      character*25 pdftyp(50,4)
      character*3 nn
      data pdftyp/
c cteq sets
     $     'CTEQ4M ','CTEQ4L ','CTEQ4HJ',
     $     'CTEQ5M ','CTEQ5L ','CTEQ5HJ',
     $     'CTEQ6M ','CTEQ6L ','CTEQ6L1',41*'CTEQ6xx ',
C MRST SETS
     $     'MRST99 ',
     $     'MRST01; as=0.119','MRST01; as=0.117','MRST01; as=0.121'
     $     ,'MRST01J; as=0.121','MRST02LO',44*' ',
C MSTW SETS
     $     'MSTW2008lo','MSTW2008nlo',40*'MSTW2008lo68clxx','MRSTLO*'
     $     ,'MRSTLO**',6*' ',
C CTQ6.6 SETS
     $     'CTQ6.6',49*' '/
      integer pdfmap(50,4)
      data pdfmap/
     $   81,83,88,   101,103, 104,   131,133,134, 41*0,
     $  111,  185,186,187,188,189, 44*0,
     $  201, 202, 48*0,
     $  301, 49*0/
c
      do i=10,50
         pdfmap(i,1)=135+(i-10)
         write(nn,'(I3)') 100+i-10
         pdftyp(i,1)='CTEQ61.'//nn(2:3)
      enddo
      do i=3,42
         pdfmap(i,3)=200+i
         write(nn,'(I2)') i-2
         pdftyp(i,3)='MSTW2008lo68cl'//nn(1:2)
      enddo
      pdfmap(43,3)=243
      pdftyp(43,3)='MRST LO*'
      pdfmap(44,3)=244
      pdftyp(44,3)='MRST LO**'
      do i=2,45
         pdfmap(i,4)=300+i
         write(nn,'(I3)') 100+i-2
         pdftyp(i,4)='CTQ66.'//nn(2:3)
      enddo
C CTEQ LO MODIFIED SETS (these are included as data files ctq66.45/46
      pdfmap(46,4)=346
      pdfmap(47,4)=347
      pdftyp(46,4)='CT09MC1'
      pdftyp(47,4)='CT09MC2'
C
      nout=pdfmap(mod(nin ,100),1+nin /100)
      type=pdftyp(mod(nin ,100),1+nin /100)
      
      end

C-------------------------------------------------------------------------
      SUBROUTINE PRNTSF(iunit)
C     prints details of the structure function sets
C-------------------------------------------------------------------------
      integer iunit
      WRITE(iunit,100)
     $  'NDNS     Set         Lambda_4    Lambda_5_2loop   Scheme'
      WRITE(iunit,100)
     $  '   1     CTEQ4M       .298        .202              MS  '
     $ ,'   2     CTEQ4L       .298        .202              MS  '
     $ ,'   3     CTEQ4HJ      .298        .202              MS  '
      WRITE(iunit,100)
     $  '   4     CTEQ5M       .326        .226  (as=0.118)  MS  ' 
     $ ,'   5     CTEQ5L *     .192        .144  (asLO=0.127)  MS  ' 
     $ ,'   6     CTEQ5HJ      .326        .226  (as=0.118)  MS  ' 
      WRITE(iunit,100)
     $  '   7     CTEQ6M       .326        .226  (as=0.118)  MS  ' 
     $ ,'   8     CTEQ6L       .326        .226  (as=0.118)  MS  ' 
     $ ,'   9     CTEQ6L1      .215        .167  (asLO=0.130)  MS  ' 
     $ ,'10-50    CTEQ6xx      .326        .226  (as=0.118)  MS  ' 
      WRITE(iunit,100)
     $  ' 101     MRST99 COR01 .321        .220              MS  '
      WRITE(iunit,100)
     $  ' 102     MRST2001     .342        .239  (as=0.119 ) MS  '
     $ ,' 103     MRST2001     .310        .214  (as=0.117 ) MS  '
     $ ,' 104     MRST2001     .378        .267  (as=0.121 ) MS  '
     $ ,' 105     MRST2001J    .378        .267  (as=0.121)  MS  '
     $ ,' 106     MRST2002LO * .215         .167  (asLO=0.13) MS  '
      WRITE(iunit,100)
     $  ' 201     MSTW2008LO *  .322*      .255* (asLO=0.139 ) MS  '
     $ ,' 202     MSTW2008NLO   .365       .255  (as=0.120 ) MS  '
c     $ ,' 202     MSTW2008NLO   .0996      .127  (as=0.120 ) MS  '
     $ ,' 203-42  MSTW2008LO68CL * .322    .255* (asLO=0.139 ) MS  '
     $ ,' 243     MRST-LOstar *  .365      .255  (as=0.120 ) MS  '
     $ ,' 244     MRST-LOstarstar *  .280  .190  (as=0.115 ) MS  '
      WRITE(iunit,100)
     $  ' 301     CTQ66        .326        .226 (as=0.118 ) MS  '
     $ ,' 302-45  CTQ66        .326        .226 (as=0.118 ) MS  '
     $ ,' 346     CT09MC1 *    .215*       .167*(asLO=0.130 ) MS  '
     $ ,' 347     CT09MC2 *    .326        .226 (as=0.118 ) MS  '
C ---------------------------------------------------------------------------
      WRITE(iunit,100)                             
     $  '  PDF sets followed by * are obtained from a 1-loop analysis,'
     $ ,'  and the relative values of Lambda  refer to 1-loop. '
     $ ,'  The MSbar scheme is used by default with 1-loop'
     $ ,'  structure functions.'
     $ ,'  In all cases the values of Lambda and loop order are set'
     $ ,'  automatically by the code, The user only needs to input ndns'
c      WRITE(iunit,100)
c     $ ,'  Sets 81-89 are the CTEQ4 fits, H.L. Lai et al.,'
c     $ ,'  CTEQ-604, hep-ph/9606399, (81=default, 82=DIS scheme,'
c     $ ,'  83=leading order, 84-87=variable Lambda, 88=High-et jet fit,'
c     $ ,'  89=low momentum evolution)'
 100  FORMAT(1X,A,100(/,1X,A))
      END

      SUBROUTINE PDFPAR(J,IH,XLAM,SCHE,NL,IRET)
      PARAMETER (NPDF=347)
C LAMBDA VALUES (lAMBDA_5FLAVOUR_2LOOP) FOR DIFFERENT PARTON DENSITIES
      IMPLICIT REAL * 8 (A-H,O-Z)
      CHARACTER * 2 SCHE,SCH(NPDF)
      INTEGER NL,NLOOP(NPDF)
      DIMENSION XLA(NPDF)
      data NLOOP/347*2/
      DATA SCH/4*'MS',3*'DI',2*'  ',
     $ 11*'MS',
     $  4*'DI','MS','DI',2*'MS',2*'  ',
     $  3*'MS',6*' ',
     $  3*'MS','DG','MS',5*'  ',
     $  5*'MS','DG',2*'MS',3*'  ',
     $  3*'MS','DI','MS',
c CTEQ3M
     $  2*'MS','DI', 2*'  ',
c MRSAp, MRSG, MRSalpha
     $  8*'MS',2*'  ',
c CTEQ4
     $  'MS','DI',7*'MS','  ',
c MRSR and MRST
     $  9*'MS','  ',
c CTEQ5
     $  'MS','DI',6*'MS',2*'  ',
c MRST99
     $  12*'MS',8*'  ',
C CTEQ6
     $  'MS', 'DI', 2*'MS', 41*'MS',5*'  ',
C MRS2002NNLO
     $  4*'MS',
C MRS2001 NLO, 2002LO
     $  5*'MS',11*'  ',
C MSTW2008, MRST LO*, LO**
     $  44*'MS',56*'  ',
C CTQ66
     $  47*'MS'/
      DATA XLA/         
c 1 DO
     $ .34D0,.68D0,.34D0,.49D0,
     $ .101D0,.173D0,.250D0,2*0.D0,
c 10     MRSA mod
     $ .151d0,.122D0,.122D0,.083D0,.101D0,.130D0,.155D0,3*.140d0,.151d0,
c 21 MT S1
     $ .138D0,.125D0,.123D0,.097D0,.138d0,2*.156d0,.245d0,2*0.D0,
     $  3*.122D0,6*0.D0,
     $  0.68D0,4*.130D0,5*0.D0,
     $  0.68D0,7*.130D0,3*0.D0,
c 61 CTEQ1M
     $  2*0.152D0,0.220D0,0.164D0,0.125D0,
c 66 CTEQ3M
     $  .158d0,.132d0,.164d0,2*0.D0,
     $  0.152D0,0.170D0,
c MRSA-alpha dependent
     $  0.09936d0,0.1396d0,0.1903d0,0.2526d0,0.3276d0,0.4162d0,2*0.D0,
c The values given above for the MRSXXX sets are consistent with the 
c alfas(Mz) given by MRS. The values    
c     $  .094d0,0.130d0,0.178d0,0.237d0,0.309d0,0.396d0/
c are on the other hand consistent with the Lambda_4 given by MRS
c CTEQ4
     $  3*0.2018d0,0.1396d0,0.1687d0,
     $  0.2392d0,0.2811d0,0.2018d0,0.1793d0,0d0,
c 91-94 MRSR
     $  .159d0,.237d0,.159d0,.237d0,
c 95-99 MRST                       
     $  3*.220d0,.164d0,.288d0,0d0,
c 101-108 cteq5                    
     $  2*.226d0,.1437,2*.226d0,2*1.d-8,.226d0,2*0.d0,
c 111-122 MRST99
     $  3*.220d0,.164d0,.288d0,.224d0,.215d0,5*.220d0,8*0d0,
C 131-174 CTEQ6: 6M, 6DI, 6L, 6L1, 40xx
     $  3*.226D0,0.166748113d0,41*.226D0,5*0d0,
C 181-184, MRS2002NNLO
     $  3*.196d0,.226d0,
C 185-189, MRS2001nlo + MRST2002LO
     $  .239,.214,.267,.267,.1668, 11*0d0,
C 201-242 (first is lambda_5(LO), second lambda_5(NLO), equal by accident
C 243-244: MRST LO*, LO**
     $ .255,.255,40*.255,0.255D0,0.190D0,56*0D0,
C 301-345 CTQ66
     $ 45*0.226D0,0.167d0,0.226D0/
      IRET=0                             
      IF(J.LT.1.OR.J.GT.NPDF) THEN
        WRITE(*,*) ' PDF SET ',J,' NOT EXISTING'
        IRET=1
        RETURN
      ENDIF
      XLAM = XLA(J)
C SCHEME
      SCHE = SCH(J)
      IF(XLAM.EQ.0.OR.SCHE.EQ.'  ') THEN
        WRITE(*,*) ' PDF SET ',J,' NOT EXISTING'
        IRET=1
        RETURN
      ENDIF
C NLOOPS
C DEFAULT IS NLOOP=2, WITH FOLLOWING EXCEPTIONS:
C CTEQ5L
      NLOOP(103)=1
C CTEQ6L1
      NLOOP(134)=1
C MRST2002LO
      NLOOP(189)=1
C MSTW2008
      NLOOP(201)=1
      do i=203,242
         NLOOP(i)=1
      enddo
C CT09MC1
      NLOOP(346)=1
C
      NL=NLOOP(J)
      END

C--------------------------------------------------
C- STRUCTURE FUNCTION MAIN PROGRAM
C--------------------------------------------------
      SUBROUTINE MLMPDF(NDNS,IH,Q2,X,FX,NF)
      REAL FX(-NF:NF),DISF(13)
      INTEGER IPAR(-6:6)
      DATA IPAR/12,11,10,9,7,8,13,2,1,3,4,5,6/
C Fix to prevent undefined math operations for x=1.
C Assumes that all structure functions vanish for x=1.
      IF(1-X.EQ.0) THEN
         DO J=-NF,NF
            FX(J) = 0
         ENDDO
         RETURN
      ENDIF
C
      IH0=IH
      IF(IH.EQ.0) IH0=1
C
      IF(NDNS.LE.89) THEN
C-- CTEQ4 FITS          
        ISET=NDNS-80
        CALL CTEQ4(ISET,IH0,Q2,X,FX,NF)
      ELSEIF(NDNS.LE.107) THEN
C-- CTEQ5 FITS          
        ISET=NDNS-100
        CALL CTEQ5(ISET,IH0,Q2,X,FX,NF)
      ELSEIF(NDNS.LE.122) THEN
C-- MRST99 sets        
        ISET=NDNS-70
#ifndef DISABLE_MRST_PDFS
        CALL HMRS(ISET,IH0,Q2,X,FX,NF)
#endif
      ELSEIF(NDNS.LE.175) THEN
C-- CTEQ6 FITS          
        ISET=NDNS-130
        if(iset.ge.5) iset=iset-5+200
        CALL CTEQ6(ISET,IH0,Q2,X,FX,NF)
C--MRST2002 FITS
      ELSEIF(NDNS.LE.189) THEN
        ISET=NDNS-(184-56)
#ifndef DISABLE_MRST_PDFS
        CALL HMRS(ISET,IH0,Q2,X,FX,NF)
#endif
C--EMPTY SLOTS
      ELSEIF(NDNS.LE.200) THEN
        WRITE(*,*) ' STRUCTURE FUNCTION SET NOT DEFINED , STOP'
        STOP
C--MSTW2008 FITS
      ELSEIF(NDNS.LE.244) THEN
        ISET=NDNS
#ifndef DISABLE_MRST_PDFS
        CALL HMRS(ISET,IH0,Q2,X,FX,NF)        
#endif       
C--EMPTY SLOTS
      ELSEIF(NDNS.LE.300) THEN
        WRITE(*,*) ' STRUCTURE FUNCTION SET NOT DEFINED , STOP'
        STOP
C--CTEQ6.6 FITS
      ELSEIF(NDNS.LE.347) THEN
        ISET=NDNS-300
        iset=iset+400-1
        CALL CTEQ6(ISET,IH0,Q2,X,FX,NF)
      ENDIF
      IF(IH.EQ.0) THEN
        FX(1)  = 0.5 * ( FX(1)+FX(2) )
        FX(-1) = 0.5 * ( FX(-1)+FX(-2) )
        FX(2)  = FX(1)
        FX(-2) = FX(-1)
      ENDIF
      END

C-------------------------------------------------------------
      subroutine qmaxlim
C------------------------------------------------------------
      implicit none
      integer icount,ilmax
      data icount/0/,ilmax/1/
      icount=icount+1
      if(log10(real(icount)).gt.ilmax) then
         write(*,*)' Q > Qmax in str. functions more than 10**',
     $        ilmax,' times'                          
          ilmax=ilmax+1
      endif
      end
C-------------------------------------------------------------
      subroutine qminlim
C------------------------------------------------------------
      implicit none
      integer icount,ilmax
      data icount/0/,ilmax/1/
      icount=icount+1
      if(log10(real(icount)).gt.ilmax) then
         write(*,*)' Q < Qmin in str. functions more than 10**',
     $        ilmax,' times'                          
          ilmax=ilmax+1
      endif
      end
C-------------------------------------------------------------
      subroutine xminlim
C------------------------------------------------------------
      implicit none
      integer icount,ilmax
      data icount/0/,ilmax/1/
      icount=icount+1
      if(log10(real(icount)).gt.ilmax) then
         write(*,*)' x < xmin in str. functions more than 10**',
     $        ilmax,' times'                          
          ilmax=ilmax+1
      endif
      end
C-------------------------------------------------------------
      subroutine xmaxlim
C------------------------------------------------------------
      implicit none
      integer icount,ilmax
      data icount/0/,ilmax/1/
      icount=icount+1
      if(log10(real(icount)).gt.ilmax) then
         write(*,*)' x > xmax in str. functions more than 10**',
     +        ilmax,' times'                          
          ilmax=ilmax+1
      endif
      end


#ifndef DISABLE_MRST_PDFS
C
C------------------------------------------------------------
C----- START HMRS ------------------------------
      SUBROUTINE  HMRS(MODE,IH,Q2,X,FX,NF)
      REAL FX(-NF:NF)
      DOUBLE PRECISION FXMSTW(-6:6)
      REAL*8 DX,DQ,UPV,DOV,SEA,USEA,DSEA,STR,CHR,BOT,GLU,XPHOT
      REAL*8 XMIN,XMAX,QSQMIN,QSQMAX,QSQ
      REAL*8 IXMIN,IXMAX,IQSQMIN,IQSQMAX
      CHARACTER prefix*50
      INTEGER IMODE,IH0
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA INI/0/
      IF(INI.GT.0) GO TO 1
        IF(MODE.EQ.10)QSQMIN=0.625D0
        IF(MODE.GT.30.AND.MODE.LE.52) QSQMIN=1.25D0
        IF(MODE.GT.30.AND.MODE.LE.52) QSQMAX=1.D7
        IF(MODE.GT.52.AND.MODE.LE.56) QSQMIN=1.25D0
        IF(MODE.GT.52.AND.MODE.LE.56) QSQMAX=1.D7
        IF(MODE.GT.200) THEN
           XMIN = 1.D-6
           XMAX = 1.D0
           QSQMIN = 1.D0
           QSQMAX = 1.D9
        ENDIF
        ILXMIN=0                    
        ILXMAX=0
        ILQSQMIN=0
        ILQSQMAX=0
1     CONTINUE
      IF(ABS(IH).GE.3) CALL NOSETP       
      IH0=IH
      IF(ABS(IH).EQ.2) IH0=ISIGN(1,IH)
      Q=SQRT(Q2)
      DQ=DBLE(Q)
      DX=DBLE(X)
      QSQ=DQ**2
      IF(DX.LT.XMIN) call xminlim
      IF(DX.GT.XMAX) call xmaxlim
      IF(QSQ.LT.QSQMIN) call qminlim
      IF(QSQ.GT.QSQMAX) call qmaxlim
      IF(MODE.LE.61) THEN
        INI=1
        IF(MODE.GT.60)THEN
          CALL mrstlo
     $              (DX,DQ,MODE-60,UPV,DOV,USEA,DSEA,STR,CHR,BOT,GLU)
        ELSEIF(MODE.GT.56)THEN
          CALL mrst2001
     $              (DX,DQ,MODE-56,UPV,DOV,USEA,DSEA,STR,CHR,BOT,GLU)
c        ELSEIF(MODE.GT.52)THEN
c          CALL mrst0201127
c     $              (DX,DQ,MODE-52,UPV,DOV,USEA,DSEA,STR,CHR,BOT,GLU)
        ELSEIF(MODE.GT.39)THEN
          CALL MRS99(DX,DQ,MODE-40,UPV,DOV,USEA,DSEA,STR,CHR,BOT,GLU)
        ENDIF
        FX(0)=SNGL(GLU)
        FX(-IH0)=SNGL(USEA)
        FX(-2*IH0)=SNGL(DSEA)
        FX(IH0)  =SNGL(UPV+USEA)
        FX(2*IH0)=SNGL(DOV+DSEA)
        IF(NF.GE.3) FX(3)=SNGL(STR)
        IF(NF.GE.4) FX(4)=SNGL(CHR)
        IF(NF.GE.5) FX(5)=SNGL(BOT)
        IF(NF.eq.6) FX(6)=0
        DO I=3,NF
           FX(-I)=FX(I)
        ENDDO
      ENDIF
      IF(MODE.GT.200)THEN
         IF(INI.EQ.0) THEN
           IF(MODE.GE.203.AND.MODE.LE.242)THEN
             prefix='mstw2008lo.68cl'
             IMODE=mode-202
           ELSE
             IF(MODE.eq.201) prefix='mstw2008lo'
             IF(MODE.eq.202) prefix='mstw2008nlo'
             IF(MODE.eq.243) prefix='mrstlostar'
             IF(MODE.eq.244) prefix='mrstlostarstar'
             IMODE=0
           ENDIF
           INI=1
         ENDIF
         CALL GetAllPDFsAlt(prefix,IMODE,DX,DQ,FXMSTW,XPHOT)
         fx(0)=REAL(FXMSTW(0))
         fx(ih0)=REAL(FXMSTW(2))
         fx(2*ih0)=REAL(FXMSTW(1))
         fx(-ih0)=REAL(FXMSTW(-2))
         fx(-2*ih0)=REAL(FXMSTW(-1))
         DO I=3,NF
            FX(I*ih0)=REAL(FXMSTW(I))
            FX(-I*ih0)=REAL(FXMSTW(-I))
         ENDDO
      ENDIF
      DO I=-NF,NF
       FX(I)=FX(I)/X
       if(fx(i).lt.0) fx(i)=0
      ENDDO
C...TRANSFORM PROTON INTO NEUTRON
      IF(ABS(IH).EQ.2) THEN
        T=FX(1)
        FX(1)=FX(2)
        FX(2)=T
        T=FX(-1)
        FX(-1)=FX(-2)
        FX(-2)=T
      ENDIF
      END
C
C
C
C
      SUBROUTINE MRSCHECK(VAL,MODE)
      REAL * 8 VAL
      CHARACTER * 10 NAME(42)
      DATA NAME/'HMRSB','KMRSB','HMRSB135','HMRSB160',
     $          'HMRSB200','HMRSB235','MRSS0','MRSD0','MRSD-',
     $          'MRSA','MRSAP','MRSG','MRS105','MRS110','MRS115',
     $          'MRS120','MRS125','MRS130',2*'EMPTY',
     $          'MRSR1','MRSR2','MRSR3','MRSR4',
     $          'MRST','MRSTH','MRSTL','MRSTM','MRSTP','EMPTY',
     $          'MRST991','MRST992','MRST993','MRST994',
     $          'MRST995','MRST996','MRST997','MRST998',
     $          'MRST999','MRST9910','MRST9911','MRST9912'/
      IF(ABS(VAL-0.00232D0).LT.0.000001) THEN                      
         IMODE = 10
      ELSEIF(ABS(VAL-0.03058D0).LT.0.000001) THEN                      
         IMODE = 1
      ELSEIF(ABS(VAL-0.01727D0).LT.0.000001) THEN
         IMODE = 2
      ELSEIF(ABS(VAL-0.01683D0).LT.0.000001) THEN
         IMODE = 3
      ELSEIF(ABS(VAL-0.01663D0).LT.0.000001) THEN
         IMODE = 4
      ELSEIF(ABS(VAL-0.01601D0).LT.0.000001) THEN
         IMODE = 5
      ELSEIF(ABS(VAL-0.01571D0).LT.0.000001) THEN
         IMODE = 6
      ELSEIF(ABS(VAL-0.01356D0).LT.0.000001) THEN
         IMODE = 7
      ELSEIF(ABS(VAL-0.00527D0).LT.0.000001) THEN
         IMODE = 8
      ELSEIF(ABS(VAL-0.00474D0).LT.0.000001) THEN
         IMODE = 9
      ELSEIF(ABS(VAL-0.00383D0).LT.0.000001) THEN
         IMODE = 10         
      ELSEIF(ABS(VAL-0.00341D0).LT.0.000001) THEN                      
         IMODE = 11         
      ELSEIF(ABS(VAL-0.00269D0).LT.0.000001) THEN                      
         IMODE = 12         
      ELSEIF(ABS(VAL-0.00429D0).LT.0.000001) THEN                      
         IMODE = 13         
      ELSEIF(ABS(VAL-0.00350D0).LT.0.000001) THEN                      
         IMODE = 14         
      ELSEIF(ABS(VAL-0.00294D0).LT.0.000001) THEN                      
         IMODE = 15         
      ELSEIF(ABS(VAL-0.00273D0).LT.0.000001) THEN                      
         IMODE = 16         
      ELSEIF(ABS(VAL-0.00195D0).LT.0.000001) THEN                      
         IMODE = 17         
      ELSEIF(ABS(VAL-0.00145D0).LT.0.000001) THEN                      
         IMODE = 18         
c
      ELSEIF(ABS(VAL-0.00150D0).LT.0.000001) THEN                      
         IMODE = 21
      ELSEIF(ABS(VAL-0.00125D0).LT.0.000001) THEN                      
         IMODE = 22
      ELSEIF(ABS(VAL-0.00181D0).LT.0.000001) THEN                      
         IMODE = 23
      ELSEIF(ABS(VAL-0.00085D0).LT.0.000001) THEN                      
         IMODE = 24         
      ELSEIF(ABS(VAL-0.00561D0).LT.0.000001) THEN                      
         IMODE = 25
      ELSEIF(ABS(VAL-0.00510D0).LT.0.000001) THEN                      
         IMODE = 26
      ELSEIF(ABS(VAL-0.00408D0).LT.0.000001) THEN                      
         IMODE = 27
      ELSEIF(ABS(VAL-0.00586D0).LT.0.000001) THEN                      
         IMODE = 28
      ELSEIF(ABS(VAL-0.00410D0).LT.0.000001) THEN                      
         IMODE = 29          
c
      ELSEIF(ABS(VAL-0.00524D0).LT.0.000001) THEN                      
         IMODE = 31
      ELSEIF(ABS(VAL-0.00497D0).LT.0.000001) THEN                      
         IMODE = 32
      ELSEIF(ABS(VAL-0.00398D0).LT.0.000001) THEN                      
         IMODE = 33
      ELSEIF(ABS(VAL-0.00585D0).LT.0.000001) THEN                      
         IMODE = 34
      ELSEIF(ABS(VAL-0.00384D0).LT.0.000001) THEN                      
         IMODE = 35
      ELSEIF(ABS(VAL-0.00177D0).LT.0.000001) THEN                      
         IMODE = 36
      ELSEIF(ABS(VAL-0.00593D0).LT.0.000001) THEN                      
         IMODE = 37
      ELSEIF(ABS(VAL-0.00541D0).LT.0.000001) THEN                      
         IMODE = 38
      ELSEIF(ABS(VAL-0.91673D0).LT.0.000001) THEN                      
         IMODE = 39
      ELSEIF(ABS(VAL-0.00525D0).LT.0.000001) THEN                      
         IMODE = 40
      ELSEIF(ABS(VAL-0.89447D0).LT.0.000001) THEN                      
         IMODE = 41
      ELSEIF(ABS(VAL-0.00515D0).LT.0.000001) THEN                      
         IMODE = 42
      ELSE
         WRITE(*,*) ' MRSCHECK: ERROR,'
         WRITE(*,*) ' NO TABLE MATCHING THE ENTRY HAS BEEN FOUND'
         STOP
      ENDIF        
      IF(IMODE.NE.MODE) THEN
         WRITE(*,*) ' MRSCHECK: ERROR,'
         WRITE(*,*) ' MRSCHECK: MODE CORRESPONDS TO ',NAME(MODE)
         WRITE(*,*) ' MRSCHECK: TABLES ARE ',NAME(IMODE)
         STOP
      ENDIF
      WRITE(*,*)' MRSCHECK: MODE ',NAME(MODE)
      END

C

      subroutine mrs99(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C****************************************************************C
C								 C
C     This is a package for the new **corrected** MRST parton    C
C     distributions. The format is similar to the previous       C
C     (1998) MRST series.                                        C
C								 C
C     NOTE: 7 new sets are added here, corresponding to shifting C
C     the small x HERA data up and down by 2.5%, and by varying  C
C     the charm and strange distributions, and by forcing a      C
C     larger d/u ratio at large x.                               C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV, MSbar factorization is assumed,     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C								 C
C     NAMING SCHEME:                                             C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
C  2     COR02  higher gluon          300      0.1175   0.00497  C
C  3     COR03  lower gluon           300      0.1175   0.00398  C
C  4     COR04  lower a_s             229      0.1125   0.00585  C
C  5     COR05  higher a_s            383      0.1225   0.00384  C
C  6     COR06  quarks up             303.3    0.1178   0.00497  C
C  7     COR07  quarks down           290.3    0.1171   0.00593  C
C  8     COR08  strange up            300      0.1175   0.00524  C
C  9     COR09  strange down          300      0.1175   0.00524  C
C  10    C0R10  charm up              300      0.1175   0.00525  C
C  11    COR11  charm down            300      0.1175   0.00524  C
C  12    COR12  larger d/u            300      0.1175   0.00515  C
C						                 C
C      The corresponding grid files are called cor01.dat etc.    C
C							  	 C
C      The reference is:                                         C
C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
C      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
C                                                                C
C      Comments to : W.J.Stirling@durham.ac.uk                   C
C                                                                C
C								 C
C****************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin) call qminlim
      if(q2.gt.qsqmax) call qmaxlim
      if(x.lt.xmin) call xminlim
      if(x.gt.xmax) call xmaxlim
      if(mode.eq.1) then
        call mrs991(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
      return
      end

      subroutine mrs991(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=1,file='cor01',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     $		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
	do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      call mrscheck(f(1,1,1),31)
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0
      init=1
      close(1)
  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     $	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
	  fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
c end of MRS99      
      subroutine mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 NLO parton           C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0110215                                  C
C                                                               C
C  There are 4 pdf sets corresponding to mode = 1, 2, 3, 4      C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.323 C
C  corresponding to alpha_s(M_Z) of 0.119                       C
C  This set reads a grid whose first number is 0.00927          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.290         C
C  corresponding to alpha_s(M_Z) of 0.117                       C
C  This set reads a grid whose first number is 0.00953          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.362         C
C  corresponding to alpha_s(M_Z) of 0.121                       C
C  This set reads a grid whose first number is 0.00889          C
C                                                               C
C  Mode=4 gives the set MRST2001J which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.353(alpha_s(M_Z) = 0.121 C 
C  This set reads a grid whose first number is 0.00826          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin) call qminlim
      if(q2.gt.qsqmax) call qmaxlim
      if(x.lt.xmin) call xminlim
      if(x.gt.xmax) call xmaxlim
      if(mode.eq.1) then
        call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.2) then
        call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.3) then
        call mrst3(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      elseif(mode.eq.4) then
        call mrst4(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
      endif 
      return
      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     $f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     $cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     $ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf119.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     $		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end
 
      subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     $f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     $cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     $ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf117.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     $		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst3(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     $f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     $cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     $ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='alf121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     $		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrst4(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     $f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     $cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     $ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='j121.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     $		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      subroutine mrstlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C***************************************************************C
C								C
C  This is a package for the new MRST 2001 LO parton            C
C  distributions.                                               C     
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0201xxx                                  C
C                                                               C
C  There is 1 pdf set corresponding to mode = 1                 C
C                                                               C
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.220 C
C  corresponding to alpha_s(M_Z) of 0.130                       C
C  This set reads a grid whose first number is 0.02868          C
C                                                               C
C   This subroutine uses an improved interpolation procedure    C 
C   for extracting values of the pdf's from the grid            C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin) call qminlim
      if(q2.gt.qsqmax) call qmaxlim
      if(x.lt.xmin) call xminlim
      if(x.gt.xmax) call xmaxlim
      if(mode.eq.1) then
        call mrst1lo(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu) 
      endif 
      return
      end

      subroutine mrst1lo(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     $f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     $cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     $ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     $	      1d-4,2d-4,4d-4,6d-4,8d-4,
     $	      1d-3,2d-3,4d-3,6d-3,8d-3,
     $	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     $	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     $	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     $	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     $	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     $        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     $        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     $        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     $        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        open(unit=33,file='lo2002.dat',status='old')
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(33,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     $		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call jeppe1(nx,nq,xxl,qql,f1,cc1)
      call jeppe1(nx,nq,xxl,qql,f2,cc2)
      call jeppe1(nx,nq,xxl,qql,f3,cc3)
      call jeppe1(nx,nq,xxl,qql,f4,cc4)
      call jeppe1(nx,nq,xxl,qql,f6,cc6)
      call jeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue
      
      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then 
      call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then 
      call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

#endif
c DISABLE_MRST_PDFS

      subroutine jeppe1(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      PARAMETER(NNX=49,MMY=37)
      dimension xx(nx),yy(my),ff(nx,my),ff1(NNX,MMY),ff2(NNX,MMY),
     $ff12(NNX,MMY),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     $cl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     $		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     $		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     $		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     $		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     $		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     $		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     $		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     $		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     $		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     $		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     $		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     $		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     $		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     $		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     $ff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     $ff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     $ff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine jeppe2(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=locx(xx,nx,x)
      m=locx(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     $       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

      integer function locx(xx,nx,x)
      implicit real*8(a-h,o-z)
      dimension xx(nx)
      if(x.le.xx(1)) then
      locx=1
      return
      endif
      if(x.ge.xx(nx)) then 
      locx=nx-1  
      return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
      jl=jm
      else
      ju=jm
      endif
      go to 1
    2 locx=jl
      return
      end


      real*8 function  polderiv(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     $(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end
c end of mrs 2002
c

C----- START CTEQ4 FITS ------------------------------
      SUBROUTINE  CTEQ4(ISET,IH,Q2,X,FX,NF)
      REAL FX(-NF:NF) 
      REAL*8 DX,DQ,CTQ4FN
C
      IF(ABS(IH).GE.3) CALL NOSETP
      IH0=IH
      IF(ABS(IH).EQ.2) IH0=ISIGN(1,IH)
      Q=SQRT(Q2)
      DQ=DBLE(Q)
      DX=DBLE(X)
C The set CTEQ4A3 (iset=6 in the CTEQ convention) is identical to
C the set CTEQ4M, and was not inserted in our package
      IF(ISET.GE.6)ISET=ISET+1
C The function CTQ4FN return the parton distribution inside the proton.
C The division by the factor DX is NOT needed
      FX(0)=SNGL(CTQ4FN(ISET,0,DX,DQ))
      FX(IH0)=SNGL(CTQ4FN(ISET,1,DX,DQ))
      FX(2*IH0)=SNGL(CTQ4FN(ISET,2,DX,DQ))
      FX(-IH0)=SNGL(CTQ4FN(ISET,-1,DX,DQ))
      FX(-2*IH0)=SNGL(CTQ4FN(ISET,-2,DX,DQ))
      DO I=3,NF
        FX(I)=SNGL(CTQ4FN(ISET,I,DX,DQ))
      ENDDO
      DO I=-NF,-3
        FX(I)=SNGL(CTQ4FN(ISET,I,DX,DQ))
      ENDDO
C...TRANSFORM PROTON INTO NEUTRON
      IF(ABS(IH).EQ.2) THEN
        T=FX(1)
        FX(1)=FX(2)
        FX(2)=T
        T=FX(-1)
        FX(-1)=FX(-2)
        FX(-2)=T
      ENDIF
      END

C============================================================================
C                CTEQ Parton Distribution Functions: Version 4
C                          June 21, 1996
C
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C   Ref: MSUHEP-60426, CTEQ-604, e-Print Archive: hep-ph/9606399
C
C   This package contains 9 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C   Iset   PDF      Description             Alpha_s(Mz)  Q0(GeV)  Table_File
C ---------------------------------------------------------------------------
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ for which Q starts at a lower value of 0.7 GeV.  
C   The Table_Files are assumed to be in the working directory.
C   
C   The function Ctq4Fn (Iset, Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar)
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.

C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq4Fn (Iset, Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Character Flnm(10)*11
      Common
     $ / CtqPar2 / Nx, Nt, NfMx
     $ / QCDtable /  Alambda, Nfl, Iorder
      Data (Flnm(I), I=1,10)
     $ / 'cteq4m', 'cteq4d', 'cteq4l'
     $ , 'cteq4a1', 'cteq4a2', 'cteq4m', 'cteq4a4'
     $ , 'cteq4a5', 'cteq4hj', 'cteq4lq' /
      Data Isetold, Isetmin, Isetmax / -987, 1, 10 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         If (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in Ctq4Fn :', Iset
	    Stop
	 Endif
	 IU= NextUt()
         Open(IU, File=Flnm(Iset), Status='OLD', Err=100)
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq4Fn: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq4Fn: ', Q
	Stop
      Endif
      If (Iparton .lt. -NfMx .or. Iparton .gt. NfMx) Then
	Print *, 'Iparton out of range in Ctq4Fn: ', Iparton
	Stop
      Endif

      Ctq4Fn = PartonX (Iparton, X, Q)
      if(Ctq4Fn.lt.0.D0)  Ctq4Fn = 0.D0

      Return

 100  Print *, ' Data file ', Flnm(Iset), ' cannot be opened '
     >//'in Ctq4Fn!!'
      Stop
C                             ********************
      End

      Function NextUt()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 50, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUt = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     $ / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     $ / CtqPar2 / Nx, Nt, NfMx
     $ / XQrange / Qini, Qmax, Xmin
     $ / QCDtable /  Alambda, Nfl, Iorder
     $ / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Common 
     $ / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     $ / CtqPar2 / Nx, Nt, NfMx
     $ / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin) Then
         ixrange=ixrange+1
         if(ixrange.eq.1) call xminlim
         If (Jx .LT. 0) Jx = 0  
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  call qminlim
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  call qmaxlim
      Endif                     

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint_dd (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint_dd (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END

      SUBROUTINE POLINT_DD (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
C          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C--- END CTEQ4 FITS -----------------------------
C    
C----- START CTEQ5 FITS ------------------------------
C Cteq5m1 (fitted form, the grid is not necessary) added on mar-23-2001 by SF
c This set seemingly supersedes Cteq5m, which was affected (?) by a bug
c in the evolution code
      SUBROUTINE  CTEQ5(ISET,IH,Q2,X,FX0,NF)
      REAL FX0(2*NF+1)  
      REAL FX(-100:100)  
      REAL*8 DX,DQ,CTQ5PDF,CTQ5PD
      double precision PDFS(-100:100)
      DATA INIT/0/ 
C
      Q=SQRT(Q2)
      DQ=DBLE(Q)
      DX=DBLE(X)
      IF(ISET.NE.1.and.ISET.NE.3)THEN
        IF(INIT.EQ.0) THEN
          CALL SETCTQ5(ISET)
          INIT=1
        ENDIF
        DO I=-NF,NF
          PDFS(I)=CTQ5PDF(I,DX,DQ)
        ENDDO
      ELSE
        DO I=-NF,NF
          PDFS(I)=CTQ5PD(ISET,I,DX,DQ,IRET) 
        ENDDO
      ENDIF
C                         
      IF(ABS(IH).GE.3) CALL NOSETP
      IH0=IH
      IF(ABS(IH).EQ.2) IH0=ISIGN(1,IH)
C The function CTQ5PDF return the parton distribution inside the proton.
C The division by the factor DX is NOT needed
      FX(0)=SNGL(PDFS(0))
      FX(IH0)=SNGL(PDFS(1))
      FX(2*IH0)=SNGL(PDFS(2))
      FX(-IH0)=SNGL(PDFS(-1))
      FX(-2*IH0)=SNGL(PDFS(-2))
      DO I=3,NF              
        FX(I)=SNGL(PDFS(I))
      ENDDO          
      DO I=-NF,-3
        FX(I)=SNGL(PDFS(I))
      ENDDO          
C...TRANSFORM PROTON INTO NEUTRON
      IF(ABS(IH).EQ.2) THEN
        T=FX(1)
        FX(1)=FX(2)
        FX(2)=T
        T=FX(-1)
        FX(-1)=FX(-2)
        FX(-2)=T
      ENDIF
c    
      do i=-nf,nf
        fx0(i+nf+1)=fx(i)
      enddo
c
      END

C============================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             March 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.0 < Q < 10,000 (GeV). 
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / K719CtqPar2 / Nx, Nt, NfMx
     > / K719QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = Parton5X (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION Parton5X (IPRTN, X, Q)
C                     
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / K719CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / K719CtqPar2 / Nx, Nt, NfMx
     > / K719XQrange / Qini, Qmax, Xmin
      Dimension Fq(M1), Df(M1)
      data ixrange/0/
      data iqmnrng/0/
      data iqmxrng/0/
      save ixrange,iqmnrng,iqmxrng
C

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      if(x.lt.xmin) call xminlim
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X < Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      if(Q.lt.Qini) call Qminlim
      if(Q.gt.Qmax) call Qmaxlim
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  then
           iqmnrng=iqmnrng+1
           if(iqmnrng.eq.1) Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q < Qini, extrapolation used; Q, Qini =', Q, Qini
         endif
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax) then
           iqmxrng=iqmxrng+1
           if(iqmxrng.eq.1) Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
         endif
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint5 (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint5 (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)
                  
      Parton5X = Ftmp
C
      RETURN
C                        ****************************
      END

      SUBROUTINE POLINT5 (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
C          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

      Subroutine SetCtq5 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=7)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5m', 'cteq5d', 'cteq5l', 'cteq5hj'
     > , 'cteq5hq', 'cteq5f3', 'cteq5f4' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl5 (IU)
         Close (IU)   
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
 101  Print*, Tablefile, ' cannot be opened '
      Stop
c      Print*, 'Please input the .tbl file:'
c      Read (*,'(A)') Tablefile
c      Goto 21
C                             ********************
      End

      Subroutine ReadTbl5 (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / K719CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / K719CtqPar2 / Nx, Nt, NfMx
     > / K719XQrange / Qini, Qmax, Xmin
     > / K719QCDtable /  Alambda, Nfl, Iorder
     > / K719Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C
c end cteq5 (grids)
c begin cteq5 fitted
C   CTEQ5M1 and CTEQ5L Parton Distribution Functions in Parametrized Form
C                             
C               September 15, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C   hep-ph/9903282
C
C   The CTEQ5M1 set given here is an updated version of the original CTEQ5M
C     set posted, in the table version, on the Web page of CTEQ.
C     The differences between CTEQ5M and CTEQ5M1 are insignificant for almost
C     all applications. 
C   The improvement is in the QCD evolution which is now more accurate, and
C   which agrees completely with the benchmark work of the HERA 96/97 Workshop.

C   The differences between the parametrized and the corresponding table ver-
C sions (on which it is based) are of similar order as between the two version.
C    
C!! Because accurate parametrizations over a wide range of (x,Q) is hard to
C   obtain, only the most widely used sets CTEQ5M and CTEQ5L are available 
C   in parametrized form for now. 

C   These parametrizations were obtained by Jon Pumplin.
C
C                    ******************************
C  Iset   PDF        Description                 Alpha_s(Mz)  Lam4  Lam5
C ---------------------------------------------------------------------------
C   1    CTEQ5M1  Standard NLO MSbar scheme         0.118     326   226
C   3    CTEQ5L   Leading Order                     0.127     192   146
C ---------------------------------------------------------------------------
C   Note the Qcd-lambda values given for CTEQ5L is for the leading order
C     form of Alpha_s!!  Alpha_s(Mz) gives the absolute calibration.

C  The two Iset value are adopted to agree with the standard table versions.

C   The following user-callable routines are provided:
C 
C     FUNCTION Ctq5Pd (Iset, Iprtn, X, Q, Irt) 
C         returns the PROBABILITY density for a GIVEN flavor;
C
C     FUNCTION Ctq5df (Iset, Iprtn, X, Q, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.
C
C     SUBROUTINE Ctq5Pds(Iset, Pdf, X, Q, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C
C   The arguments of these routines are as follows: 
C
C   Iset is the set number:  1 for CTEQ5M1 or 3 for CTEQ5L  
C
C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark and 1 as U-quark, 
C               which might be different from your labels.
C
C   X, Q are the usual x, Q; 
C
C   Irt is an error code: 0 if there was no error; 1 or more if (x,q) was 
C   outside the range of validity of the parametrization.
C       
C  Range of validity:
C  
C     The range of (x, Q) covered by this parametrization of the QCD evolved
C     parton distributions is 1E-6 < x < 1 ; 1.1 GeV < Q < 10 TeV.  Of course,
C     the PDF's are constrained by data only in a subset of that region; and 
C     the assumed DGLAP evolution is unlikely to be valid for all of it either.
C
C     The range of (x, Q) used in the CTEQ5 round of global analysis is 
C     approximately 0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for 
C     fixed target experiments; 0.0001 < x < 0.3 from HERA data; and   
C     Q^2 up to 40,000 GeV^2 from Tevatron inclusive Jet data.
C
C   DOUBLE PRECISION is used throughout in these routines, but conversion to 
C   SINGLE PRECISION is possible by removing the Implicit Double Precision statements. 
C
C **************************************************************************

C ********************************************************
      FUNCTION CTQ5PD(ISET, IPARTON, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

c if called at a point (x,q) that is outside the region that was 
c actually parametrized, return a value of 0, and set the error code IRT=1.  
c The user can remove the following IF statement to receive instead an 
c extrapolated value, which may be wildly unphysical.
      if((x .lt. 1.e-6). or. (x .gt. 1.) 
     &	 .or. (q .lt. .99) .or. (q .gt. 10000.)) then
         ctq5pd = 0.d0
         irt = 1
         return
      endif

      irt = 0
      if(iset .eq. 3) then
         ctq5pd = ctq5L(iparton,x,q)
      elseif(iset .eq. 1) then
         ctq5pd = ctq5Mi(iparton,x,q)
      else
         print *,'iset=',iset,' has not been parametrized.' 
	   print '(/A)', 'Use the interpolation-table version instead.'
         stop
      endif

      return
      end

C ********************************************************
      FUNCTION CTQ5DF(ISET, IFL, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CTQ5DF = X * CTQ5PD(ISET, IPARTON, X, Q, IRT)
        
      RETURN
      END

C ********************************************************
      SUBROUTINE CTQ5PDS(ISET, PDF, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION PDF (-6:6)

      IRT = 0

      DO IFL= -6,2
         PDF(IFL) = CTQ5PD(ISET,IFL,X,Q,IRT1)
         IRT = IRT + IRT1

         IF (IFL .LE. -3) THEN
            PDF(-IFL) = PDF(IFL)
         ENDIF

      ENDDO

      RETURN
      END

c --------------------------------------------------------------------------
	double precision function ctq5MI(ifl,x,q)
c Parametrization of cteq5MI parton distribution functions (J. Pumplin 9/99).
c ifl: 1=u,2=d,3=s,4=c,5=b;0=gluon;-1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar.
c --------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	integer ifl

	ii = ifl
	if(ii .gt. 2) then
	   ii = -ii
	endif

	if(ii .eq. -1) then
	   sum = faux5MI(-1,x,q)
	   ratio = faux5MI(-2,x,q)
	   ctq5MI = sum/(1.d0 + ratio)

	elseif(ii .eq. -2) then
	   sum = faux5MI(-1,x,q)
	   ratio = faux5MI(-2,x,q)
	   ctq5MI = sum*ratio/(1.d0 + ratio)

	elseif(ii .ge. -5) then
	   ctq5MI = faux5MI(ii,x,q)

	else
	   ctq5MI = 0.d0 

	endif

	return
	end

c ---------------------------------------------------------------------
      double precision function faux5MI(ifl,x,q)
c auxiliary function for parametrization of CTEQ5MI (J. Pumplin 9/99).
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ifl

      parameter (nex=8, nlf=2)
      dimension am(0:nex,0:nlf,-5:2)
      dimension alfvec(-5:2), qmavec(-5:2)
      dimension mexvec(-5:2), mlfvec(-5:2)
      dimension ut1vec(-5:2), ut2vec(-5:2)
      dimension af(0:nex)

      data mexvec( 2) / 8 /
      data mlfvec( 2) / 2 /
      data ut1vec( 2) /  0.5141718E+01 /
      data ut2vec( 2) / -0.1346944E+01 /
      data alfvec( 2) /  0.5260555E+00 /
      data qmavec( 2) /  0.0000000E+00 /
      data (am( 0,k, 2),k=0, 2)
     & /  0.4289071E+01, -0.2536870E+01, -0.1259948E+01 /
      data (am( 1,k, 2),k=0, 2)
     & /  0.9839410E+00,  0.4168426E-01, -0.5018952E-01 /
      data (am( 2,k, 2),k=0, 2)
     & / -0.1651961E+02,  0.9246261E+01,  0.5996400E+01 /
      data (am( 3,k, 2),k=0, 2)
     & / -0.2077936E+02,  0.9786469E+01,  0.7656465E+01 /
      data (am( 4,k, 2),k=0, 2)
     & /  0.3054926E+02,  0.1889536E+01,  0.1380541E+01 /
      data (am( 5,k, 2),k=0, 2)
     & /  0.3084695E+02, -0.1212303E+02, -0.1053551E+02 /
      data (am( 6,k, 2),k=0, 2)
     & / -0.1426778E+02,  0.6239537E+01,  0.5254819E+01 /
      data (am( 7,k, 2),k=0, 2)
     & / -0.1909811E+02,  0.3695678E+01,  0.5495729E+01 /
      data (am( 8,k, 2),k=0, 2)
     & /  0.1889751E-01,  0.5027193E-02,  0.6624896E-03 /

      data mexvec( 1) / 8 /
      data mlfvec( 1) / 2 /
      data ut1vec( 1) /  0.4138426E+01 /
      data ut2vec( 1) / -0.3221374E+01 /
      data alfvec( 1) /  0.4960962E+00 /
      data qmavec( 1) /  0.0000000E+00 /
      data (am( 0,k, 1),k=0, 2)
     & /  0.1332497E+01, -0.3703718E+00,  0.1288638E+00 /
      data (am( 1,k, 1),k=0, 2)
     & /  0.7544687E+00,  0.3255075E-01, -0.4706680E-01 /
      data (am( 2,k, 1),k=0, 2)
     & / -0.7638814E+00,  0.5008313E+00, -0.9237374E-01 /
      data (am( 3,k, 1),k=0, 2)
     & / -0.3689889E+00, -0.1055098E+01, -0.4645065E+00 /
      data (am( 4,k, 1),k=0, 2)
     & /  0.3991610E+02,  0.1979881E+01,  0.1775814E+01 /
      data (am( 5,k, 1),k=0, 2)
     & /  0.6201080E+01,  0.2046288E+01,  0.3804571E+00 /
      data (am( 6,k, 1),k=0, 2)
     & / -0.8027900E+00, -0.7011688E+00, -0.8049612E+00 /
      data (am( 7,k, 1),k=0, 2)
     & / -0.8631305E+01, -0.3981200E+01,  0.6970153E+00 /
      data (am( 8,k, 1),k=0, 2)
     & /  0.2371230E-01,  0.5372683E-02,  0.1118701E-02 /

      data mexvec( 0) / 8 /
      data mlfvec( 0) / 2 /
      data ut1vec( 0) / -0.1026789E+01 /
      data ut2vec( 0) / -0.9051707E+01 /
      data alfvec( 0) /  0.9462977E+00 /
      data qmavec( 0) /  0.0000000E+00 /
      data (am( 0,k, 0),k=0, 2)
     & /  0.1191990E+03, -0.8548739E+00, -0.1963040E+01 /
      data (am( 1,k, 0),k=0, 2)
     & / -0.9449972E+02,  0.1074771E+01,  0.2056055E+01 /
      data (am( 2,k, 0),k=0, 2)
     & /  0.3701064E+01, -0.1167947E-02,  0.1933573E+00 /
      data (am( 3,k, 0),k=0, 2)
     & /  0.1171345E+03, -0.1064540E+01, -0.1875312E+01 /
      data (am( 4,k, 0),k=0, 2)
     & / -0.1014453E+03, -0.5707427E+00,  0.4511242E-01 /
      data (am( 5,k, 0),k=0, 2)
     & /  0.6365168E+01,  0.1275354E+01, -0.4964081E+00 /
      data (am( 6,k, 0),k=0, 2)
     & / -0.3370693E+01, -0.1122020E+01,  0.5947751E-01 /
      data (am( 7,k, 0),k=0, 2)
     & / -0.5327270E+01, -0.9293556E+00,  0.6629940E+00 /
      data (am( 8,k, 0),k=0, 2)
     & /  0.2437513E-01,  0.1600939E-02,  0.6855336E-03 /

      data mexvec(-1) / 8 /
      data mlfvec(-1) / 2 /
      data ut1vec(-1) /  0.5243571E+01 /
      data ut2vec(-1) / -0.2870513E+01 /
      data alfvec(-1) /  0.6701448E+00 /
      data qmavec(-1) /  0.0000000E+00 /
      data (am( 0,k,-1),k=0, 2)
     & /  0.2428863E+02,  0.1907035E+01, -0.4606457E+00 /
      data (am( 1,k,-1),k=0, 2)
     & /  0.2006810E+01, -0.1265915E+00,  0.7153556E-02 /
      data (am( 2,k,-1),k=0, 2)
     & / -0.1884546E+02, -0.2339471E+01,  0.5740679E+01 /
      data (am( 3,k,-1),k=0, 2)
     & / -0.2527892E+02, -0.2044124E+01,  0.1280470E+02 /
      data (am( 4,k,-1),k=0, 2)
     & / -0.1013824E+03, -0.1594199E+01,  0.2216401E+00 /
      data (am( 5,k,-1),k=0, 2)
     & /  0.8070930E+02,  0.1792072E+01, -0.2164364E+02 /
      data (am( 6,k,-1),k=0, 2)
     & / -0.4641050E+02,  0.1977338E+00,  0.1273014E+02 /
      data (am( 7,k,-1),k=0, 2)
     & / -0.3910568E+02,  0.1719632E+01,  0.1086525E+02 /
      data (am( 8,k,-1),k=0, 2)
     & / -0.1185496E+01, -0.1905847E+00, -0.8744118E-03 /

      data mexvec(-2) / 7 /
      data mlfvec(-2) / 2 /
      data ut1vec(-2) /  0.4782210E+01 /
      data ut2vec(-2) / -0.1976856E+02 /
      data alfvec(-2) /  0.7558374E+00 /
      data qmavec(-2) /  0.0000000E+00 /
      data (am( 0,k,-2),k=0, 2)
     & / -0.6216935E+00,  0.2369963E+00, -0.7909949E-02 /
      data (am( 1,k,-2),k=0, 2)
     & /  0.1245440E+01, -0.1031510E+00,  0.4916523E-02 /
      data (am( 2,k,-2),k=0, 2)
     & / -0.7060824E+01, -0.3875283E-01,  0.1784981E+00 /
      data (am( 3,k,-2),k=0, 2)
     & / -0.7430595E+01,  0.1964572E+00, -0.1284999E+00 /
      data (am( 4,k,-2),k=0, 2)
     & / -0.6897810E+01,  0.2620543E+01,  0.8012553E-02 /
      data (am( 5,k,-2),k=0, 2)
     & /  0.1507713E+02,  0.2340307E-01,  0.2482535E+01 /
      data (am( 6,k,-2),k=0, 2)
     & / -0.1815341E+01, -0.1538698E+01, -0.2014208E+01 /
      data (am( 7,k,-2),k=0, 2)
     & / -0.2571932E+02,  0.2903941E+00, -0.2848206E+01 /

      data mexvec(-3) / 7 /
      data mlfvec(-3) / 2 /
      data ut1vec(-3) /  0.4518239E+01 /
      data ut2vec(-3) / -0.2690590E+01 /
      data alfvec(-3) /  0.6124079E+00 /
      data qmavec(-3) /  0.0000000E+00 /
      data (am( 0,k,-3),k=0, 2)
     & / -0.2734458E+01, -0.7245673E+00, -0.6351374E+00 /
      data (am( 1,k,-3),k=0, 2)
     & /  0.2927174E+01,  0.4822709E+00, -0.1088787E-01 /
      data (am( 2,k,-3),k=0, 2)
     & / -0.1771017E+02, -0.1416635E+01,  0.8467622E+01 /
      data (am( 3,k,-3),k=0, 2)
     & / -0.4972782E+02, -0.3348547E+01,  0.1767061E+02 /
      data (am( 4,k,-3),k=0, 2)
     & / -0.7102770E+01, -0.3205337E+01,  0.4101704E+00 /
      data (am( 5,k,-3),k=0, 2)
     & /  0.7169698E+02, -0.2205985E+01, -0.2463931E+02 /
      data (am( 6,k,-3),k=0, 2)
     & / -0.4090347E+02,  0.2103486E+01,  0.1416507E+02 /
      data (am( 7,k,-3),k=0, 2)
     & / -0.2952639E+02,  0.5376136E+01,  0.7825585E+01 /

      data mexvec(-4) / 7 /
      data mlfvec(-4) / 2 /
      data ut1vec(-4) /  0.2783230E+01 /
      data ut2vec(-4) / -0.1746328E+01 /
      data alfvec(-4) /  0.1115653E+01 /
      data qmavec(-4) /  0.1300000E+01 /
      data (am( 0,k,-4),k=0, 2)
     & / -0.1743872E+01, -0.1128921E+01, -0.2841969E+00 /
      data (am( 1,k,-4),k=0, 2)
     & /  0.3345755E+01,  0.3187765E+00,  0.1378124E+00 /
      data (am( 2,k,-4),k=0, 2)
     & / -0.2037615E+02,  0.4121687E+01,  0.2236520E+00 /
      data (am( 3,k,-4),k=0, 2)
     & / -0.4703104E+02,  0.5353087E+01, -0.1455347E+01 /
      data (am( 4,k,-4),k=0, 2)
     & / -0.1060230E+02, -0.1551122E+01, -0.1078863E+01 /
      data (am( 5,k,-4),k=0, 2)
     & /  0.5088892E+02, -0.8197304E+01,  0.8083451E+01 /
      data (am( 6,k,-4),k=0, 2)
     & / -0.2819070E+02,  0.4554086E+01, -0.5890995E+01 /
      data (am( 7,k,-4),k=0, 2)
     & / -0.1098238E+02,  0.2590096E+01, -0.8062879E+01 /

      data mexvec(-5) / 6 /
      data mlfvec(-5) / 2 /
      data ut1vec(-5) /  0.1619654E+02 /
      data ut2vec(-5) / -0.3367346E+01 /
      data alfvec(-5) /  0.5109891E-02 /
      data qmavec(-5) /  0.4500000E+01 /
      data (am( 0,k,-5),k=0, 2)
     & / -0.6800138E+01,  0.2493627E+01, -0.1075724E+01 /
      data (am( 1,k,-5),k=0, 2)
     & /  0.3036555E+01,  0.3324733E+00,  0.2008298E+00 /
      data (am( 2,k,-5),k=0, 2)
     & / -0.5203879E+01, -0.8493476E+01, -0.4523208E+01 /
      data (am( 3,k,-5),k=0, 2)
     & / -0.1524239E+01, -0.3411912E+01, -0.1771867E+02 /
      data (am( 4,k,-5),k=0, 2)
     & / -0.1099444E+02,  0.1320930E+01, -0.2353831E+01 /
      data (am( 5,k,-5),k=0, 2)
     & /  0.1699299E+02, -0.3565802E+02,  0.3566872E+02 /
      data (am( 6,k,-5),k=0, 2)
     & / -0.1465793E+02,  0.2703365E+02, -0.2176372E+02 /

      if(q .le. qmavec(ifl)) then
         faux5MI = 0.d0
         return
      endif

      if(x .ge. 1.d0) then
         faux5MI = 0.d0
         return
      endif

      tmp = log(q/alfvec(ifl))
      if(tmp .le. 0.d0) then
         faux5MI = 0.d0
         return
      endif

      sb = log(tmp)
      sb1 = sb - 1.2d0
      sb2 = sb1*sb1

      do i = 0, nex
         af(i) = 0.d0
         sbx = 1.d0
         do k = 0, mlfvec(ifl)
            af(i) = af(i) + sbx*am(i,k,ifl)
            sbx = sb1*sbx
         enddo
      enddo

      y = -log(x)
      u = log(x/0.00001d0)

      part1 = af(1)*y**(1.d0+0.01d0*af(4))*(1.d0+ af(8)*u)
      part2 = af(0)*(1.d0 - x) + af(3)*x 
      part3 = x*(1.d0-x)*(af(5)+af(6)*(1.d0-x)+af(7)*x*(1.d0-x))
      part4 = ut1vec(ifl)*log(1.d0-x) + 
     &	      AF(2)*log(1.d0+exp(ut2vec(ifl))-x)

      faux5MI = exp(log(x) + part1 + part2 + part3 + part4)

c include threshold factor...
      faux5MI = faux5MI * (1.d0 - qmavec(ifl)/q)

      return
      end
c --------------------------------------------------------------------------
	double precision function ctq5L(ifl,x,q)
c Parametrization of cteq5L parton distribution functions (J. Pumplin 9/99).
c ifl: 1=u,2=d,3=s,4=c,5=b;0=gluon;-1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar.
c --------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	integer ifl

	ii = ifl
	if(ii .gt. 2) then
	   ii = -ii
	endif

	if(ii .eq. -1) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5L = sum/(1.d0 + ratio)

	elseif(ii .eq. -2) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5L = sum*ratio/(1.d0 + ratio)

	elseif(ii .ge. -5) then
	   ctq5L = faux5L(ii,x,q)

	else
	   ctq5L = 0.d0 

	endif

	return
	end

c ---------------------------------------------------------------------
      double precision function faux5L(ifl,x,q)
c auxiliary function for parametrization of CTEQ5L (J. Pumplin 9/99).
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ifl

      parameter (nex=8, nlf=2)
      dimension am(0:nex,0:nlf,-5:2)
      dimension alfvec(-5:2), qmavec(-5:2)
      dimension mexvec(-5:2), mlfvec(-5:2)
      dimension ut1vec(-5:2), ut2vec(-5:2)
      dimension af(0:nex)

      data mexvec( 2) / 8 /
      data mlfvec( 2) / 2 /
      data ut1vec( 2) /  0.4971265E+01 /
      data ut2vec( 2) / -0.1105128E+01 /
      data alfvec( 2) /  0.2987216E+00 /
      data qmavec( 2) /  0.0000000E+00 /
      data (am( 0,k, 2),k=0, 2)
     & /  0.5292616E+01, -0.2751910E+01, -0.2488990E+01 /
      data (am( 1,k, 2),k=0, 2)
     & /  0.9714424E+00,  0.1011827E-01, -0.1023660E-01 /
      data (am( 2,k, 2),k=0, 2)
     & / -0.1651006E+02,  0.7959721E+01,  0.8810563E+01 /
      data (am( 3,k, 2),k=0, 2)
     & / -0.1643394E+02,  0.5892854E+01,  0.9348874E+01 /
      data (am( 4,k, 2),k=0, 2)
     & /  0.3067422E+02,  0.4235796E+01, -0.5112136E+00 /
      data (am( 5,k, 2),k=0, 2)
     & /  0.2352526E+02, -0.5305168E+01, -0.1169174E+02 /
      data (am( 6,k, 2),k=0, 2)
     & / -0.1095451E+02,  0.3006577E+01,  0.5638136E+01 /
      data (am( 7,k, 2),k=0, 2)
     & / -0.1172251E+02, -0.2183624E+01,  0.4955794E+01 /
      data (am( 8,k, 2),k=0, 2)
     & /  0.1662533E-01,  0.7622870E-02, -0.4895887E-03 /

      data mexvec( 1) / 8 /
      data mlfvec( 1) / 2 /
      data ut1vec( 1) /  0.2612618E+01 /
      data ut2vec( 1) / -0.1258304E+06 /
      data alfvec( 1) /  0.3407552E+00 /
      data qmavec( 1) /  0.0000000E+00 /
      data (am( 0,k, 1),k=0, 2)
     & /  0.9905300E+00, -0.4502235E+00,  0.1624441E+00 /
      data (am( 1,k, 1),k=0, 2)
     & /  0.8867534E+00,  0.1630829E-01, -0.4049085E-01 /
      data (am( 2,k, 1),k=0, 2)
     & /  0.8547974E+00,  0.3336301E+00,  0.1371388E+00 /
      data (am( 3,k, 1),k=0, 2)
     & /  0.2941113E+00, -0.1527905E+01,  0.2331879E+00 /
      data (am( 4,k, 1),k=0, 2)
     & /  0.3384235E+02,  0.3715315E+01,  0.8276930E+00 /
      data (am( 5,k, 1),k=0, 2)
     & /  0.6230115E+01,  0.3134639E+01, -0.1729099E+01 /
      data (am( 6,k, 1),k=0, 2)
     & / -0.1186928E+01, -0.3282460E+00,  0.1052020E+00 /
      data (am( 7,k, 1),k=0, 2)
     & / -0.8545702E+01, -0.6247947E+01,  0.3692561E+01 /
      data (am( 8,k, 1),k=0, 2)
     & /  0.1724598E-01,  0.7120465E-02,  0.4003646E-04 /

      data mexvec( 0) / 8 /
      data mlfvec( 0) / 2 /
      data ut1vec( 0) / -0.4656819E+00 /
      data ut2vec( 0) / -0.2742390E+03 /
      data alfvec( 0) /  0.4491863E+00 /
      data qmavec( 0) /  0.0000000E+00 /
      data (am( 0,k, 0),k=0, 2)
     & /  0.1193572E+03, -0.3886845E+01, -0.1133965E+01 /
      data (am( 1,k, 0),k=0, 2)
     & / -0.9421449E+02,  0.3995885E+01,  0.1607363E+01 /
      data (am( 2,k, 0),k=0, 2)
     & /  0.4206383E+01,  0.2485954E+00,  0.2497468E+00 /
      data (am( 3,k, 0),k=0, 2)
     & /  0.1210557E+03, -0.3015765E+01, -0.1423651E+01 /
      data (am( 4,k, 0),k=0, 2)
     & / -0.1013897E+03, -0.7113478E+00,  0.2621865E+00 /
      data (am( 5,k, 0),k=0, 2)
     & / -0.1312404E+01, -0.9297691E+00, -0.1562531E+00 /
      data (am( 6,k, 0),k=0, 2)
     & /  0.1627137E+01,  0.4954111E+00, -0.6387009E+00 /
      data (am( 7,k, 0),k=0, 2)
     & /  0.1537698E+00, -0.2487878E+00,  0.8305947E+00 /
      data (am( 8,k, 0),k=0, 2)
     & /  0.2496448E-01,  0.2457823E-02,  0.8234276E-03 /

      data mexvec(-1) / 8 /
      data mlfvec(-1) / 2 /
      data ut1vec(-1) /  0.3862583E+01 /
      data ut2vec(-1) / -0.1265969E+01 /
      data alfvec(-1) /  0.2457668E+00 /
      data qmavec(-1) /  0.0000000E+00 /
      data (am( 0,k,-1),k=0, 2)
     & /  0.2647441E+02,  0.1059277E+02, -0.9176654E+00 /
      data (am( 1,k,-1),k=0, 2)
     & /  0.1990636E+01,  0.8558918E-01,  0.4248667E-01 /
      data (am( 2,k,-1),k=0, 2)
     & / -0.1476095E+02, -0.3276255E+02,  0.1558110E+01 /
      data (am( 3,k,-1),k=0, 2)
     & / -0.2966889E+01, -0.3649037E+02,  0.1195914E+01 /
      data (am( 4,k,-1),k=0, 2)
     & / -0.1000519E+03, -0.2464635E+01,  0.1964849E+00 /
      data (am( 5,k,-1),k=0, 2)
     & /  0.3718331E+02,  0.4700389E+02, -0.2772142E+01 /
      data (am( 6,k,-1),k=0, 2)
     & / -0.1872722E+02, -0.2291189E+02,  0.1089052E+01 /
      data (am( 7,k,-1),k=0, 2)
     & / -0.1628146E+02, -0.1823993E+02,  0.2537369E+01 /
      data (am( 8,k,-1),k=0, 2)
     & / -0.1156300E+01, -0.1280495E+00,  0.5153245E-01 /

      data mexvec(-2) / 7 /
      data mlfvec(-2) / 2 /
      data ut1vec(-2) /  0.1895615E+00 /
      data ut2vec(-2) / -0.3069097E+01 /
      data alfvec(-2) /  0.5293999E+00 /
      data qmavec(-2) /  0.0000000E+00 /
      data (am( 0,k,-2),k=0, 2)
     & / -0.6556775E+00,  0.2490190E+00,  0.3966485E-01 /
      data (am( 1,k,-2),k=0, 2)
     & /  0.1305102E+01, -0.1188925E+00, -0.4600870E-02 /
      data (am( 2,k,-2),k=0, 2)
     & / -0.2371436E+01,  0.3566814E+00, -0.2834683E+00 /
      data (am( 3,k,-2),k=0, 2)
     & / -0.6152826E+01,  0.8339877E+00, -0.7233230E+00 /
      data (am( 4,k,-2),k=0, 2)
     & / -0.8346558E+01,  0.2892168E+01,  0.2137099E+00 /
      data (am( 5,k,-2),k=0, 2)
     & /  0.1279530E+02,  0.1021114E+00,  0.5787439E+00 /
      data (am( 6,k,-2),k=0, 2)
     & /  0.5858816E+00, -0.1940375E+01, -0.4029269E+00 /
      data (am( 7,k,-2),k=0, 2)
     & / -0.2795725E+02, -0.5263392E+00,  0.1290229E+01 /

      data mexvec(-3) / 7 /
      data mlfvec(-3) / 2 /
      data ut1vec(-3) /  0.3753257E+01 /
      data ut2vec(-3) / -0.1113085E+01 /
      data alfvec(-3) /  0.3713141E+00 /
      data qmavec(-3) /  0.0000000E+00 /
      data (am( 0,k,-3),k=0, 2)
     & /  0.1580931E+01, -0.2273826E+01, -0.1822245E+01 /
      data (am( 1,k,-3),k=0, 2)
     & /  0.2702644E+01,  0.6763243E+00,  0.7231586E-02 /
      data (am( 2,k,-3),k=0, 2)
     & / -0.1857924E+02,  0.3907500E+01,  0.5850109E+01 /
      data (am( 3,k,-3),k=0, 2)
     & / -0.3044793E+02,  0.2639332E+01,  0.5566644E+01 /
      data (am( 4,k,-3),k=0, 2)
     & / -0.4258011E+01, -0.5429244E+01,  0.4418946E+00 /
      data (am( 5,k,-3),k=0, 2)
     & /  0.3465259E+02, -0.5532604E+01, -0.4904153E+01 /
      data (am( 6,k,-3),k=0, 2)
     & / -0.1658858E+02,  0.2923275E+01,  0.2266286E+01 /
      data (am( 7,k,-3),k=0, 2)
     & / -0.1149263E+02,  0.2877475E+01, -0.7999105E+00 /

      data mexvec(-4) / 7 /
      data mlfvec(-4) / 2 /
      data ut1vec(-4) /  0.4400772E+01 /
      data ut2vec(-4) / -0.1356116E+01 /
      data alfvec(-4) /  0.3712017E-01 /
      data qmavec(-4) /  0.1300000E+01 /
      data (am( 0,k,-4),k=0, 2)
     & / -0.8293661E+00, -0.3982375E+01, -0.6494283E-01 /
      data (am( 1,k,-4),k=0, 2)
     & /  0.2754618E+01,  0.8338636E+00, -0.6885160E-01 /
      data (am( 2,k,-4),k=0, 2)
     & / -0.1657987E+02,  0.1439143E+02, -0.6887240E+00 /
      data (am( 3,k,-4),k=0, 2)
     & / -0.2800703E+02,  0.1535966E+02, -0.7377693E+00 /
      data (am( 4,k,-4),k=0, 2)
     & / -0.6460216E+01, -0.4783019E+01,  0.4913297E+00 /
      data (am( 5,k,-4),k=0, 2)
     & /  0.3141830E+02, -0.3178031E+02,  0.7136013E+01 /
      data (am( 6,k,-4),k=0, 2)
     & / -0.1802509E+02,  0.1862163E+02, -0.4632843E+01 /
      data (am( 7,k,-4),k=0, 2)
     & / -0.1240412E+02,  0.2565386E+02, -0.1066570E+02 /

      data mexvec(-5) / 6 /
      data mlfvec(-5) / 2 /
      data ut1vec(-5) /  0.5562568E+01 /
      data ut2vec(-5) / -0.1801317E+01 /
      data alfvec(-5) /  0.4952010E-02 /
      data qmavec(-5) /  0.4500000E+01 /
      data (am( 0,k,-5),k=0, 2)
     & / -0.6031237E+01,  0.1992727E+01, -0.1076331E+01 /
      data (am( 1,k,-5),k=0, 2)
     & /  0.2933912E+01,  0.5839674E+00,  0.7509435E-01 /
      data (am( 2,k,-5),k=0, 2)
     & / -0.8284919E+01,  0.1488593E+01, -0.8251678E+00 /
      data (am( 3,k,-5),k=0, 2)
     & / -0.1925986E+02,  0.2805753E+01, -0.3015446E+01 /
      data (am( 4,k,-5),k=0, 2)
     & / -0.9480483E+01, -0.9767837E+00, -0.1165544E+01 /
      data (am( 5,k,-5),k=0, 2)
     & /  0.2193195E+02, -0.1788518E+02,  0.9460908E+01 /
      data (am( 6,k,-5),k=0, 2)
     & / -0.1327377E+02,  0.1201754E+02, -0.6277844E+01 /

      if(q .le. qmavec(ifl)) then
         faux5L = 0.d0
         return
      endif

      if(x .ge. 1.d0) then
         faux5L = 0.d0
         return
      endif

      tmp = log(q/alfvec(ifl))
      if(tmp .le. 0.d0) then
         faux5L = 0.d0
         return
      endif

      sb = log(tmp)
      sb1 = sb - 1.2d0
      sb2 = sb1*sb1

      do i = 0, nex
         af(i) = 0.d0
         sbx = 1.d0
         do k = 0, mlfvec(ifl)
            af(i) = af(i) + sbx*am(i,k,ifl)
            sbx = sb1*sbx
         enddo
      enddo

      y = -log(x)
      u = log(x/0.00001d0)

      part1 = af(1)*y**(1.d0+0.01d0*af(4))*(1.d0+ af(8)*u)
      part2 = af(0)*(1.d0 - x) + af(3)*x 
      part3 = x*(1.d0-x)*(af(5)+af(6)*(1.d0-x)+af(7)*x*(1.d0-x))
      part4 = ut1vec(ifl)*log(1.d0-x) + 
     &	      AF(2)*log(1.d0+exp(ut2vec(ifl))-x)

      faux5L = exp(log(x) + part1 + part2 + part3 + part4)

c include threshold factor...
      faux5L = faux5L * (1.d0 - qmavec(ifl)/q)

      return
      end

c end cteq5 fitted
C    
C----- START CTEQ6 FITS ------------------------------
C Cteq6, added by P. Nason on 4-2-2002
C Cteq61, update by mlm on july 10 2003
      SUBROUTINE  CTEQ6(ISET,IH,Q2,X,FX0,NF)
      REAL FX0(2*NF+1)  
      REAL FX(-100:100)  
      REAL*8 DX,DQ,CTQ6PDF,CTQ6PD,PDFS(-100:100)
      DATA INIT/0/ 
C
      Q=SQRT(Q2)
      DQ=DBLE(Q)
      DX=DBLE(X)
      IF(INIT.EQ.0) THEN
         CALL SETCTQ6(ISET)
         INIT=1
      ENDIF
      DO I=-NF,NF
         PDFS(I)=CTQ6PDF(I,DX,DQ)
      ENDDO
C                         
      IF(ABS(IH).GE.3) CALL NOSETP
      IH0=IH
      IF(ABS(IH).EQ.2) IH0=ISIGN(1,IH)
C The function CTQ6PDF return the parton distribution inside the proton.
C The division by the factor DX is NOT needed
      FX(0)=SNGL(PDFS(0))
      FX(IH0)=SNGL(PDFS(1))
      FX(2*IH0)=SNGL(PDFS(2))
      FX(-IH0)=SNGL(PDFS(-1))
      FX(-2*IH0)=SNGL(PDFS(-2))
      DO I=3,NF              
        FX(I)=SNGL(PDFS(I))
      ENDDO          
      DO I=-NF,-3
        FX(I)=SNGL(PDFS(I))
      ENDDO          
C...TRANSFORM PROTON INTO NEUTRON
      IF(ABS(IH).EQ.2) THEN
        T=FX(1)
        FX(1)=FX(2)
        FX(2)=T
        T=FX(-1)
        FX(-1)=FX(-2)
        FX(-2)=T
      ENDIF
!
      do i=-nf,nf
        fx0(i+nf+1)=fx(i)
      enddo
!
      END

C============================================================================
C                CTEQ Parton Distribution Functions: version 6.0-6.6
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C                             August 6, 2003, v6.11
C                             December 12, 2004, v6.12
C                             December 4, 2006, v6.5 (CTEQ6.5M series added)
C                             March 23, 2007, v6.51 (CTEQ6.5S/C series added)
C                             April 24, 2007, v6.52 (minor improvement)
C                             March 30, 2008, v6.6 
C                             March 23, 2010, v6.6as (AS seies added)
C                             March 29, 2010, v6.6as (QCD_Lambda added for AS)
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
C       JHEP 0310:046(2003), hep-ph/0303013
C
C   Ref[3]: "Neutrino dimuon Production and Strangeness Asymmetry of the Nucleon"
C       By: F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky, H.L. Lai, S. Kretzer, J.F. Owens, W.K. Tung
C       Eur. Phys. J. C40:145(2005), hep-ph/0312323
C
C   Ref[4]: "CTEQ6 Parton Distributions with Heavy Quark Mass Effects"
C       By: S. Kretzer, H.L. Lai, F. Olness, W.K. Tung
C       Phys. Rev. D69:114005(2004), hep-ph/0307022
C
C   Ref[5]: "Heavy Quark Mass Effects in Deep Inelastic Scattering and Global QCD Analysis"
C       By : W.K. Tung, H.L. Lai, A. Belyaev, J. Pumplin, D. Stump, C.-P. Yuan
C       JHEP 0702:053(2007), hep-ph/0611254
C
C   Ref[6]: "The Strange Parton Distribution of Nucleon: Global Analysis and Applications"
C       By : H.L. Lai, P. Nadolsky, J. Pumplin, D. Stump, W.K. Tung, C.-P. Yuan
C       JHEP 0704:089,2007, hep-ph/0702268
C
C   Ref[7]: "The Charm Content of the Nucleon"
C       By : J. Pumplin, H.L. Lai, W.K. Tung
C       Phys.Rev.D75:054029,2007, hep-ph/0701220

C   Ref[8]: "Implications of CTEQ global analysis for collider observables"
C       By : P. M. Nadolsky, H.-L. Lai, Q.-H. Cao, J. Huston, J. Pumplin, D. R. Stump, W.-K. Tung, C.-P. Yuan
C       arXiv:0802.0007 [hep-ph], submitted to Phys. Rev. D. 
C

C   Ref[9]: TBA

C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
C   (4) 5 special sets for strangeness study from Ref[3].
C   (5) 1 special set for heavy quark study from Ref[4].
C   (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref[5].
C   (7) 8 sets of PDFs resulting from the strangeness study, Ref[6].
C   (8) 7 sets of PDFs resulting from the charm study, Ref[7].
C   (9) CTEQ6.6M and its 44 up/down eigenvector sets from Ref[8].
C  (10) Fits with nonperturbative charm from the study in  Ref[8].
C  (11) Fits with alternative values of the strong coupling strength from the study in Ref[9].


C  Details about the calling convention are:
C --------------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File   Ref
C ================================================================================
C Standard, "best-fit", sets:                 
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl    [1]
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl    [1]
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl    [1]
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl   [1]
C 200    CTEQ6.1M: updated CTEQ6M (see below, under "uncertainty" section)     [2]
C 400    CTEQ6.6M; the 2008 set (see below, under "uncertainty" section)       [8]
C
C --------------------------
C  Special sets with nonperturbative charm at Q_0=1.3 GeV from Ref [8]
C --------------------------
C 450    CTEQ6.6C1   BHPS model for IC    0.118     326   226    ctq66.c1.pds
C 451    CTEQ6.6C2   BHPS model for IC    0.118     326   226    ctq66.c2.pds
C 452    CTEQ6.6C3   Sea-like model       0.118     326   226    ctq66.c3.pds
C 453    CTEQ6.6C4   Sea-like model       0.118     326   226    ctq66.c4.pds
C     Momentum Fraction carried by c+cbar=2c at Q0=1.3 GeV:
C    Iset:     451  452   453   454 
C Mom. frac:  0.01 0.035  0.01  0.035


C --------------------------
C  Special CTEQ6.6AS sets with alternative values of strong coupling strength [9]
C --------------------------
C 470    CTEQ6.6AS1                       0.116           202    ctq66.as1.pds
C 471    CTEQ6.6AS2                       0.117           214    ctq66.as2.pds
C 472    CTEQ6.6AS3                       0.119           239    ctq66.as3.pds
C 473    CTEQ6.6AS4                       0.120           251    ctq66.as4.pds

C --------------------------
C Special sets for strangeness study:  Ref.[3]
C --------------------------
C  11    CTEQ6A   Class A                 0.118     326   226    cteq6sa.pds
C  12    CTEQ6B   Class B                 0.118     326   226    cteq6sb.pds
C  13    CTEQ6C   Class C                 0.118     326   226    cteq6sc.pds
C  14    CTEQ6B+  Large [S-]              0.118     326   226    cteq6sb+.pds
C  15    CTEQ6B-  Negative [S-]           0.118     326   226    cteq6sb-.pds
C --------------------------
C Special set for Heavy Quark study:   Ref.[4]
C --------------------------
C  21    CTEQ6HQ                          0.118     326   226    cteq6hq.pds
C --------------------------
C Released sets for strangeness study:  Ref.[6]
C -------------------------- s=sbr
C  30    CTEQ6.5S0   Best-fit             0.118     326   226    ctq65.s+0.pds
C  31    CTEQ6.5S1   Low s+               0.118     326   226    ctq65.s+1.pds
C  32    CTEQ6.5S2   High s+              0.118     326   226    ctq65.s+2.pds
C  33    CTEQ6.5S3   Alt Low s+           0.118     326   226    ctq65.s+3.pds
C  34    CTEQ6.5S4   Alt High s+          0.118     326   226    ctq65.s+4.pds
C -------------------------- s!=sbr
C          strangeness asymmetry <x>_s-
C  35    CTEQ6.5S-0  Best-fit    0.0014    0.118     326   226    ctq65.s-0.pds
C  36    CTEQ6.5S-1  Low        -0.0010    0.118     326   226    ctq65.s-1.pds
C  37    CTEQ6.5S-2  High        0.0050    0.118     326   226    ctq65.s-2.pds
C --------------------------
C Released sets for charm study:  Ref.[7]
C --------------------------
C  40    CTEQ6.5C0   no intrinsic charm   0.118     326   226    ctq65.c0.pds
C  41    CTEQ6.5C1   BHPS model for IC    0.118     326   226    ctq65.c1.pds
C  42    CTEQ6.5C2   BHPS model for IC    0.118     326   226    ctq65.c2.pds
C  43    CTEQ6.5C3   Meson cloud model    0.118     326   226    ctq65.c3.pds
C  44    CTEQ6.5C4   Meson cloud model    0.118     326   226    ctq65.c4.pds
C  45    CTEQ6.5C5   Sea-like model       0.118     326   226    ctq65.c5.pds
C  46    CTEQ6.5C6   Sea-like model       0.118     326   226    ctq65.c6.pds
C
C     Momentum Fraction carried by c,cbar at Q0=1.3 GeV:
C    Iset:charm  ,cbar     | Iset:charm  ,cbar     | Iset:charm  ,cbar
C    41: 0.002857,0.002857 | 43: 0.003755,0.004817 | 45: 0.005714,0.005714
C    42: 0.010000,0.010000 | 44: 0.007259,0.009312 | 46: 0.012285,0.012285
C
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects, Ref[5]:  central fit: CTEQ6.5M (=CTEQ65.00)
C                              -----------------------
C  3xx  CTEQ65.xx  +/- sets               0.118     326   226    ctq65.xx.pds
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 300      is CTEQ65.00 (=CTEQ6.5M),
C             301/302 are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects and free strangeness, Ref[8]:  
C                central fit: CTEQ6.6M (=CTEQ66.00)
C                              -----------------------
C  4xx  CTEQ66.xx  +/- sets               0.118     326   226    ctq66.xx.pds
C        where xx = 01-44: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 400      is CTEQ66.00 (=CTEQ6.6M),
C             401/402 are CTEQ66.01/02, +/- sets of 1st eigenvector, ... etc.

C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for 
C    *  10^-8 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.6 series;
C    *  10^-7 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.5S/C series;
C    *  10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV) for CTEQ6, CTEQ6.1 series;
C
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   nadolsky@physics.smu.edu, pumplin@pa.msu.edu or hllai@tmue.edu.tw.
C
C===========================================================================

      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / Ctq6Par2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0d0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Ctq6Pdf = 0D0
        Return
      Endif

      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif

      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if (Ctq6Pdf.lt.0.D0) Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=8)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.', 'ctq66.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
c ct09MC: 444 -> 446
      Data Isetmin4,Isetmax4 /400,446/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetmnC08,IsetmxC08 /450,453/
      Data IsetmnAS08,IsetmxAS08 /470,473/

      Data IsetHQ /21/
      Common /Setchange/ Isetch
      Common /Ctq6Jset/ Jset
      save
#ifdef USE_MPI
      include 'mpif.h'
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / Ctq6Par1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / Ctq6Par2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
#endif

      Jset=Iset
C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
        fmtpds=.true.

        If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                  Iset = 1,2,3 for 6m, 6d, 6l
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'.tbl'
        Elseif (Iset.eq.4) Then
C                                                             4  (2nd LO fit)
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'1.tbl'
        Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(1)//nn//'.tbl'
        Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(5)//nn(2:3)//'.tbl'
        Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15
          If(Iset.eq.11) then
            Tablefile=Flnm(6)//'a.pds'
          Elseif(Iset.eq.12) then
            Tablefile=Flnm(6)//'b.pds'
          Elseif(Iset.eq.13) then
            Tablefile=Flnm(6)//'c.pds'
          Elseif(Iset.eq.14) then
            Tablefile=Flnm(6)//'b+.pds'
          Elseif(Iset.eq.15) then
            Tablefile=Flnm(6)//'b-.pds'
          Endif
        Elseif (Iset.eq.IsetHQ) Then
C                                                               21
          TableFile='cteq6hq.pds'
        Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'s+'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37
          Is = Iset - 5
          write(nn,'(I2)') Is
          Tablefile=Flnm(7)//'s-'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'c'//nn(2:2)//'.pds'
        Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340
          write(nn,'(I3)') Iset
          Tablefile=Flnm(7)//nn(2:3)//'.pds'
        Elseif (Iset.ge.Isetmin4 .and. Iset.le.Isetmax4) Then
C                                                    (Cteq6.6)  400 - 444   
          write(nn,'(I3)') Iset
          Tablefile=Flnm(8)//nn(2:3)//'.pds'
        Elseif (Iset.ge.IsetmnC08 .and. Iset.le.IsetmxC08) Then
C                                                   (Cteq6.6C)  450 - 453
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.ge.IsetmnAS08 .and. Iset.le.IsetmxAS08) Then
C                                                   (Cteq6.6AS)  470 - 473
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'as'//nn(3:3)//'.pds'
        Else
          Print *, 'Invalid Iset number in SetCtq6 :', Iset
          Stop
        Endif
        IU= NextUn()
#ifdef USE_MPI
        if(mpirank.eq.0) then
         Open(IU, File=Tablefile, Status='OLD', Err=100)
 21      Call Readpds (IU,fmtpds)
         Close (IU)
        endif
        call MPI_BCAST(Al,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Al'
         call exit(-1)
        endif
        call MPI_BCAST(XV,MXX+1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting XV'
         call exit(-1)
        endif
        call MPI_BCAST(TV,MXQ+1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting TV'
         call exit(-1)
        endif
        call MPI_BCAST(UPD,MXPQX,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting UPD'
         call exit(-1)
        endif
        call MPI_BCAST(Nx,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Nx'
         call exit(-1)
        endif
        call MPI_BCAST(Nt,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Nt'
         call exit(-1)
        endif
        call MPI_BCAST(NfMx,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting NfMx'
         call exit(-1)
        endif
        call MPI_BCAST(MxVal,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting MxVal'
         call exit(-1)
        endif
        call MPI_BCAST(Qini,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Qini'
         call exit(-1)
        endif
        call MPI_BCAST(Qmax,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Qmax'
         call exit(-1)
        endif
        call MPI_BCAST(Xmin,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Xmin'
         call exit(-1)
        endif
        call MPI_BCAST(Alambda,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Alambda'
         call exit(-1)
        endif
        call MPI_BCAST(Nfl,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Nfl'
         call exit(-1)
        endif
        call MPI_BCAST(Iorder,1,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Iorder'
         call exit(-1)
        endif
        call MPI_BCAST(Amass,6,MPI_DOUBLE_PRECISION,0,
     $                 MPI_COMM_WORLD,retval)
        if(retval.ne.0) then
         write(6,*) ' Error broadcasting Amass'
         call exit(-1)
        endif


#else
        Open(IU, File=Tablefile, Status='OLD', Err=100)
 21     Call Readpds (IU,fmtpds)
        Close (IU)
#endif
        Isetold=Iset
        Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >  //'in SetCtq6!!'
#ifdef USE_MPI
      call MPI_ABORT(-1)
#else
      Stop
#endif
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / Ctq6Par1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / Ctq6Par2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      Common /Ctq6Jset/ Jset
      Data IsetmnAS08,IsetmxAS08 /470,473/

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      If((Jset.ge.IsetmnAS08) .and.(Jset.le.IsetmxAS08)) then
c    for CTEQ6.6AS series, Al is not Lambda_QCD, but Qbase so the Qqrids would be the same as standard CTEQ6.6 series. 
c    Hardwired Lambda_QCD for those sets
         If(Jset.eq.470) then
	    Alambda=.2018d0
         elseIf(Jset.eq.471) then
	    Alambda=.2138d0
         elseIf(Jset.eq.472) then
	    Alambda=.2392d0
         elseIf(Jset.eq.473) then
	    Alambda=.2526d0
	 endif
      else
          Alambda = Al
      endif

      Read  (Nu, '(A)') Line
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

      Common
     > / Ctq6Par1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / Ctq6Par2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END

C
#ifndef DISABLE_MRST_PDFS

C START MSTW
C----------------------------------------------------------------------
C--   Fortran interpolation code for MSTW PDFs, building on existing
C--   MRST Fortran code and Jeppe Andersen's C++ code.
C--   Three user interfaces:
C--    call GetAllPDFs(prefix,ih,x,q,upv,dnv,usea,dsea,
C--                    str,sbar,chm,cbar,bot,bbar,glu,phot)
C--    call GetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)
C--    xf = GetOnePDF(prefix,ih,x,q,f)
C--   See enclosed example.f for usage.
C--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>.
C----------------------------------------------------------------------
C--   Alternative LHAPDF-like interface: return PDFs in an array.
      subroutine GetAllPDFsAlt(prefix,ih,x,q,xpdf,xphoton)
      implicit none
      integer ih,f
      double precision x,q,xpdf(-6:6),xphoton,xvalence,GetOnePDF
      character*(*) prefix

      do f = 1, 6
         xpdf(f) = GetOnePDF(prefix,ih,x,q,f) ! quarks
         xvalence = GetOnePDF(prefix,ih,x,q,f+6) ! valence quarks
         xpdf(-f) = xpdf(f) - xvalence ! antiquarks
      end do
      xpdf(0) = GetOnePDF(prefix,ih,x,q,0) ! gluon
      xphoton = GetOnePDF(prefix,ih,x,q,13) ! photon
      
      return
      end

C----------------------------------------------------------------------

C--   Get only one parton flavour 'f', using PDG notation:
C--    f =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
C--      = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
C--   Can also get valence quarks directly:
C--    f =  7, 8, 9,10,11,12.
C--      = dv,uv,sv,cv,bv,tv.
C--   Photon: f = 13.
      double precision function GetOnePDF(prefix,ih,x,q,f)
      implicit none
      logical warn,fatal
      parameter(warn=.false.,fatal=.true.)
C--   Set warn=.true. to turn on warnings when extrapolating.
C--   Set fatal=.false. to return zero instead of terminating when
C--    invalid input values of x and q are used.
      integer ih,f,nhess,nx,nq,np,nqc0,nqb0,nqc,nqb,n,m,ip,io,
     &     alphaSorder,alphaSnfmax,nExtraFlavours
      double precision x,q,xmin,xmax,qsqmin,qsqmax,mc2,mb2,eps,
     &     dummy,qsq,xlog,qsqlog,res,res1,anom,ExtrapolatePDF,
     &     InterpolatePDF,distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      parameter(nx=64,nq=48,np=12,nqc0=4,nqb0=14,
     &     nqc=nq-nqc0,nqb=nq-nqb0)
      parameter(xmin=1d-6,xmax=1d0,qsqmin=1d0,qsqmax=1d9,eps=1d-6)
      parameter(nhess=2*20)
      character set*2,prefix*(*),filename*60,oldprefix(0:nhess)*50
      character dummyChar,dummyWord*50
      double precision ff(np,nx,nq)
      double precision qq(nq),xx(nx),cc(np,0:nhess,nx,nq,4,4)
      double precision xxl(nx),qql(nq)
C--   Store distance along each eigenvector, tolerance,
C--   heavy quark masses and alphaS parameters in COMMON block.
      common/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      save
      data xx/1d-6,2d-6,4d-6,6d-6,8d-6,
     &     1d-5,2d-5,4d-5,6d-5,8d-5,
     &     1d-4,2d-4,4d-4,6d-4,8d-4,
     &     1d-3,2d-3,4d-3,6d-3,8d-3,
     &     1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     &	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     &	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     &	   .5d0,.525d0,.55d0,.575d0,.6d0,.625d0,.65d0,.675d0,
     &     .7d0,.725d0,.75d0,.775d0,.8d0,.825d0,.85d0,.875d0,
     &     .9d0,.925d0,.95d0,.975d0,1d0/
      data qq/1.d0,
     &     1.25d0,1.5d0,0.d0,0.d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,
     &     1d1,1.2d1,0.d0,0.d0,2.6d1,4d1,6.4d1,1d2,
     &     1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     &     1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     &     1.8d6,3.2d6,5.6d6,1d7,1.8d7,3.2d7,5.6d7,1d8,
     &     1.8d8,3.2d8,5.6d8,1d9/

      if (f.lt.-6.or.f.gt.13) then
         print *,"Error: invalid parton flavour = ",f
         stop
      end if

      if (ih.lt.0.or.ih.gt.nhess) then
         print *,"Error: invalid eigenvector number = ",ih
         stop
      end if

C--   Check if the requested parton set is already in memory.
      if (oldprefix(ih).ne.prefix) then

C--   Start of initialisation for eigenvector set "i" ...
C--   Do this only the first time the set "i" is called,
C--   OR if the prefix has changed from the last time.

C--   Check that the character arrays "oldprefix" and "filename"
C--   are large enough.
         if (len_trim(prefix).gt.len(oldprefix(ih))) then
            print *,"Error in GetOnePDF: increase size of oldprefix"
            stop
         else if (len_trim(prefix)+7.gt.len(filename)) then
            print *,"Error in GetOnePDF: increase size of filename"
            stop
         end if

         write(set,'(I2.2)') ih  ! convert integer to string
C--   Remove trailing blanks from prefix before assigning filename.
         filename = prefix(1:len_trim(prefix))//'.'//set//'.dat'
C--   Line below can be commented out if you don't want this message.
         write(6,*) "Reading PDF grid from ",
     $              filename(1:len_trim(filename))
         open(unit=33,file=filename,iostat=io,status='old')
         if (io.ne.0) then
            print *,"Error in GetOnePDF: can't open ",
     &           filename(1:len_trim(filename))
            stop
         end if

C--   Read header containing heavy quark masses and alphaS values.
         read(33,*) 
         read(33,*)
         read(33,*) dummyChar,dummyWord,dummyWord,dummyChar,
     &        distance,tolerance
         read(33,*) dummyChar,dummyWord,dummyChar,mCharm
         read(33,*) dummyChar,dummyWord,dummyChar,mBottom
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSQ0
         read(33,*) dummyChar,dummyWord,dummyChar,alphaSMZ
         read(33,*) dummyChar,dummyWord,dummyWord,dummyChar,
     &        alphaSorder,alphaSnfmax
         read(33,*) dummyChar,dummyWord,dummyChar,nExtraFlavours
         read(33,*)
         read(33,*)
         mc2=mCharm**2
         mb2=mBottom**2
         qq(4)=mc2
         qq(5)=mc2+eps
         qq(14)=mb2
         qq(15)=mb2+eps

C--   Check that the heavy quark masses are sensible.
         if (mc2.lt.qq(3).or.mc2.gt.qq(6)) then
            print *,"Error in GetOnePDF: invalid mCharm = ",mCharm
            stop
         end if
         if (mb2.lt.qq(13).or.mb2.gt.qq(16)) then
            print *,"Error in GetOnePDF: invalid mBottom = ",mBottom
            stop
         end if
         
C--   The nExtraFlavours variable is provided to aid compatibility
C--   with future grids where, for example, a photon distribution
C--   might be provided (cf. the MRST2004QED PDFs).
         if (nExtraFlavours.lt.0.or.nExtraFlavours.gt.1) then
            print *,"Error in GetOnePDF: invalid nExtraFlavours = ",
     &           nExtraFlavours
            stop
         end if

C--   Now read in the grids from the grid file.
         do n=1,nx-1
            do m=1,nq
               if (nExtraFlavours.gt.0) then
                  if (alphaSorder.eq.2) then ! NNLO
                     read(33,'(12(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,12)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     read(33,'(10(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9),ff(12,n,m)
                  end if
               else             ! nExtraFlavours = 0
                  if (alphaSorder.eq.2) then ! NNLO
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(11(1pe12.4))',iostat=io)
     &                 (ff(ip,n,m),ip=1,11)
                  else          ! LO or NLO
                     ff(10,n,m) = 0.d0 ! = chm-cbar
                     ff(11,n,m) = 0.d0 ! = bot-bbar
                     ff(12,n,m) = 0.d0 ! = photon
                     read(33,'(9(1pe12.4))',iostat=io)
     &                    (ff(ip,n,m),ip=1,9)
                  end if
               end if
               if (io.ne.0) then
                  print *,"Error in GetOnePDF reading ",filename
                  stop
               end if
            enddo
         enddo

C--   Check that ALL the file contents have been read in.
         read(33,*,iostat=io) dummy
         if (io.eq.0) then
            print *,"Error in GetOnePDF: not at end of ",filename
            stop
         end if
         close(unit=33)

C--   PDFs are identically zero at x = 1.
         do m=1,nq
            do ip=1,np
               ff(ip,nx,m)=0d0
            enddo
         enddo

         do n=1,nx
            xxl(n)=log10(xx(n))
         enddo
         do m=1,nq
            qql(m)=log10(qq(m))
         enddo

C--   Initialise all parton flavours.
         do ip=1,np
            call InitialisePDF(ip,np,ih,nhess,nx,nq,nqc0,nqb0,
     &           xxl,qql,ff,cc)
         enddo

         oldprefix(ih) = prefix

C--   ... End of initialisation for eigenvector set "ih".

      end if                    ! oldprefix(ih).ne.prefix

C----------------------------------------------------------------------

      qsq=q*q
C--   If mc2 < qsq < mc2+eps, then qsq = mc2+eps.
      if (qsq.gt.qq(nqc0).and.qsq.lt.qq(nqc0+1)) qsq = qq(nqc0+1)
C--   If mb2 < qsq < mb2+eps, then qsq = mb2+eps.
      if (qsq.gt.qq(nqb0).and.qsq.lt.qq(nqb0+1)) qsq = qq(nqb0+1)
      
      xlog=log10(x)
      qsqlog=log10(qsq)

      res = 0.d0

      if (f.eq.0) then          ! gluon
         ip = 1
      else if (f.ge.1.and.f.le.5) then ! quarks
         ip = f+1
      else if (f.le.-1.and.f.ge.-5) then ! antiquarks
         ip = -f+1
      else if (f.ge.7.and.f.le.11) then ! valence quarks
         ip = f
      else if (f.eq.13) then    ! photon
         ip = 12
      else if (abs(f).ne.6.and.f.ne.12) then
         if (warn.or.fatal) print *,"Error in GetOnePDF: f = ",f
         if (fatal) stop
      end if

      if (x.le.0.d0.or.x.gt.xmax.or.q.le.0.d0) then

         if (warn.or.fatal) print *,"Error in GetOnePDF: x,qsq = ",
     &        x,qsq,q,xmax
         if (fatal) stop

      else if (abs(f).eq.6.or.f.eq.12) then ! set top quarks to zero
         
         res = 0.d0

      else if (qsq.lt.qsqmin) then ! extrapolate to low Q^2

         if (warn) then
            print *, "Warning in GetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         if (x.lt.xmin) then    ! extrapolate to low x

            res = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         else                   ! do usual interpolation
            
            res = InterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(qsqmin),nx,nq,xxl,qql,cc)
            res1 = InterpolatePDF(ip,np,ih,nhess,xlog,
     &           log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
               res = res - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(qsqmin),nx,nq,xxl,qql,cc)
               res1 = res1 - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &              log10(1.01D0*qsqmin),nx,nq,xxl,qql,cc)
            end if
            
         end if

C--   Calculate the anomalous dimension, dlog(xf)/dlog(qsq),
C--   evaluated at qsqmin.  Then extrapolate the PDFs to low
C--   qsq < qsqmin by interpolating the anomalous dimenion between
C--   the value at qsqmin and a value of 1 for qsq << qsqmin.
C--   If value of PDF at qsqmin is very small, just set
C--   anomalous dimension to 1 to prevent rounding errors.
         if (abs(res).ge.1.D-5) then
            anom = (res1-res)/res/0.01D0
         else
            anom = 1.D0
         end if
         res = res*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin)

      else if (x.lt.xmin.or.qsq.gt.qsqmax) then ! extrapolate

         if (warn) then
            print *, "Warning in GetOnePDF, extrapolating: f = ",f,
     &           ", x = ",x,", q = ",q
         end if

         res = ExtrapolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)
         
         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - ExtrapolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if

      else                      ! do usual interpolation
         
         res = InterpolatePDF(ip,np,ih,nhess,xlog,
     &        qsqlog,nx,nq,xxl,qql,cc)

         if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
            res = res - InterpolatePDF(ip+5,np,ih,nhess,xlog,
     &           qsqlog,nx,nq,xxl,qql,cc)
         end if
            
      end if
      
      GetOnePDF = res

      return
      end

C----------------------------------------------------------------------

      subroutine InitialisePDF(ip,np,ih,nhess,nx,my,myc0,myb0,
     &     xx,yy,ff,cc)
      implicit none
      integer nhess,ih,nx,my,myc0,myb0,j,k,l,m,n,ip,np
      double precision xx(nx),yy(my),ff(np,nx,my),
     &     ff1(nx,my),ff2(nx,my),ff12(nx,my),ff21(nx,my),
     &     yy0(4),yy1(4),yy2(4),yy12(4),z(16),
     &     cl(16),cc(np,0:nhess,nx,my,4,4),iwt(16,16),
     &     polderiv1,polderiv2,polderiv3,d1,d2,d1d2,xxd

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     &     -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     &     2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     &     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     &     0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     &     0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     &     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     &     9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     &     -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     &     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     &     -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     &     4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/

      do m=1,my
         ff1(1,m)=polderiv1(xx(1),xx(2),xx(3),
     &        ff(ip,1,m),ff(ip,2,m),ff(ip,3,m))
         ff1(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff(ip,nx-2,m),ff(ip,nx-1,m),ff(ip,nx,m))
         do n=2,nx-1
            ff1(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff(ip,n-1,m),ff(ip,n,m),ff(ip,n+1,m))
         enddo
      enddo

C--   Calculate the derivatives at qsq=mc2,mc2+eps,mb2,mb2+eps
C--   in a similar way as at the endpoints qsqmin and qsqmax.
      do n=1,nx
         do m=1,my
            if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff2(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff(ip,n,m),ff(ip,n,m+1),ff(ip,n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff2(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff(ip,n,m-2),ff(ip,n,m-1),ff(ip,n,m))
            else
               ff2(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff(ip,n,m-1),ff(ip,n,m),ff(ip,n,m+1))
            end if
         end do
      end do

C--   Calculate the cross derivatives (d/dx)(d/dy).
      do m=1,my
         ff12(1,m)=polderiv1(xx(1),xx(2),xx(3),
     &        ff2(1,m),ff2(2,m),ff2(3,m))
         ff12(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx),
     &        ff2(nx-2,m),ff2(nx-1,m),ff2(nx,m))
         do n=2,nx-1
            ff12(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1),
     &           ff2(n-1,m),ff2(n,m),ff2(n+1,m))
         enddo
      enddo

C--   Calculate the cross derivatives (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff21(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2),
     &              ff1(n,m),ff1(n,m+1),ff1(n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff21(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m),
     &              ff1(n,m-2),ff1(n,m-1),ff1(n,m))
            else
               ff21(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1),
     &              ff1(n,m-1),ff1(n,m),ff1(n,m+1))
            end if
         end do
      end do

C--   Take the average of (d/dx)(d/dy) and (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            ff12(n,m)=0.5*(ff12(n,m)+ff21(n,m))
         end do
      end do

      do n=1,nx-1
         do m=1,my-1
            d1=xx(n+1)-xx(n)
            d2=yy(m+1)-yy(m)
            d1d2=d1*d2
            
            yy0(1)=ff(ip,n,m)
            yy0(2)=ff(ip,n+1,m)
            yy0(3)=ff(ip,n+1,m+1)
            yy0(4)=ff(ip,n,m+1)
            
            yy1(1)=ff1(n,m)
            yy1(2)=ff1(n+1,m)
            yy1(3)=ff1(n+1,m+1)
            yy1(4)=ff1(n,m+1)
            
            yy2(1)=ff2(n,m)
            yy2(2)=ff2(n+1,m)
            yy2(3)=ff2(n+1,m+1)
            yy2(4)=ff2(n,m+1)
            
            yy12(1)=ff12(n,m)
            yy12(2)=ff12(n+1,m)
            yy12(3)=ff12(n+1,m+1)
            yy12(4)=ff12(n,m+1)
            
            do k=1,4
               z(k)=yy0(k)
               z(k+4)=yy1(k)*d1
               z(k+8)=yy2(k)*d2
               z(k+12)=yy12(k)*d1d2
            enddo
            
            do l=1,16
               xxd=0.d0
               do k=1,16
                  xxd=xxd+iwt(k,l)*z(k)
               enddo
               cl(l)=xxd
            enddo
            l=0
            do k=1,4
               do j=1,4
                  l=l+1
                  cc(ip,ih,n,m,k,j)=cl(l)
               enddo
            enddo
         enddo
      enddo
      return
      end
#endif
c DISABLE_MRST_PDFS

C----------------------------------------------------------------------

      double precision function InterpolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,locxMSTW,l,m,n,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,t,u

      n=locxMSTW(xx,nx,x)
      m=locxMSTW(yy,my,y)
      
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(ip,ih,n,m,l,4)*u+cc(ip,ih,n,m,l,3))*u
     .        +cc(ip,ih,n,m,l,2))*u+cc(ip,ih,n,m,l,1)
      enddo

      InterpolatePDF = z

      return
      end

C----------------------------------------------------------------------

      double precision function ExtrapolatePDF(ip,np,ih,nhess,x,y,
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,locxMSTW,n,m,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4),
     &     x,y,z,f0,f1,z0,z1,InterpolatePDF
      
      n=locxMSTW(xx,nx,x)           ! 0: below xmin, nx: above xmax
      m=locxMSTW(yy,my,y)           ! 0: below qsqmin, my: above qsqmax
      
C--   If extrapolation in small x only:
      if (n.eq.0.and.m.gt.0.and.m.lt.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(1),y,nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(2),y,nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
C--   If extrapolation into large q only:
      else if (n.gt.0.and.m.eq.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,x,yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,x,yy(my-1),nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
C--   If extrapolation into large q AND small x:
      else if (n.eq.0.and.m.eq.my) then
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(1),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(1),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         f0 = InterpolatePDF(ip,np,ih,nhess,xx(2),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolatePDF(ip,np,ih,nhess,xx(2),yy(my-1),nx,my,xx,yy,
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.1.d-3.and.z1.gt.1.d-3) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if
      else
         print *,"Error in ExtrapolatePDF"
         stop
      end if

      ExtrapolatePDF = z      

      return
      end

C----------------------------------------------------------------------

      integer function locxMSTW(xx,nx,x)
C--   returns an integer j such that x lies inbetween xx(j) and xx(j+1).
C--   nx is the length of the array with xx(nx) the highest element.
      implicit none
      integer nx,jl,ju,jm
      double precision x,xx(nx)
      if(x.eq.xx(1)) then
         locxMSTW=1
         return
      endif
      if(x.eq.xx(nx)) then
         locxMSTW=nx-1  
         return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
         jl=jm
      else
         ju=jm
      endif
      go to 1
    2 locxMSTW=jl
      return
      end

C----------------------------------------------------------------------

      double precision function polderiv1(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x1 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv1=(x3*x3*(y1-y2)+2.d0*x1*(x3*(-y1+y2)+x2*(y1-y3))
     &     +x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv2(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x2 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv2=(x3*x3*(y1-y2)-2.d0*x2*(x3*(y1-y2)+x1*(y2-y3))
     &     +x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv3(x1,x2,x3,y1,y2,y3)
C--   returns the estimate of the derivative at x3 obtained by a
C--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv3=(x3*x3*(-y1+y2)+2.d0*x2*x3*(y1-y3)+x1*x1*(y2-y3)
     &     +x2*x2*(-y1+y3)+2.d0*x1*x3*(-y2+y3))/
     &     ((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

C----------------------------------------------------------------------

C END MSTW


C-------------------------------------------------------------------
      SUBROUTINE NOSETP
      WRITE(*,*) ' SET OF STRUCTURE FUNCTIONS NOT IMPLEMENTED'
      WRITE(*,*) ' OR'
      WRITE(*,*) ' HADRON TYPE NOT DESCRIBED BY THE REQUESTED SET:'
      WRITE(*,*) ' IH =    1    2    3    -1    -2    -3    4    5'
      WRITE(*,*) ' HAD=    P    N    PI+  PBAR  NBAR  PI-  PH   EL'
      STOP
      END
C
C----------------------------------------------------------------------------
C-------------------------------------------------------------------
C------- ALPHA QCD -------------------------------------
c Program to calculate alfa strong with nf flavours,
c as a function of lambda with 5 flavors.
c The value of alfa is matched at the thresholds q = mq.
c When invoked with nf < 0 it chooses nf as the number of
c flavors with mass less then q.
c
      function alfas(q2,xlam,nloop,inf)
      implicit real * 8 (a-h,o-z)
      data olam/0.d0/,pi/3.14159d0/
      data xmb/4.5d0/,xmc/1.5d0/
      integer nloop
      if(xlam.ne.olam) then
        olam = xlam
        b5  = (33-2*5)/pi/12
        b4  = (33-2*4)/pi/12
        b3  = (33-2*3)/pi/12
        if(nloop.eq.1) then
           bp5=0
           bp4=0
           bp3=0
        elseif(nloop.eq.2) then
           bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
           bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
           bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        endif 
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     $        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     $        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
        if( q .gt. xmb ) then
          nf = 5
        elseif( q .gt. xmc ) then
          nf = 4
        else
          nf = 3
        endif
      endif
      if    ( nf .eq. 5 ) then
        alfas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        alfas = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
      elseif( nf .eq. 3 ) then
        alfas = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        stop
      endif
      return
      end

c     same a alfas, but using different xlam and nloop. to avoid
c     multiple reinitializations
      function alfas_clu(q2,xlam,nloop,inf)
      implicit real * 8 (a-h,o-z)
      data olam/0.d0/,pi/3.14159d0/
      data xmb/4.5d0/,xmc/1.5d0/
      integer nloop
      if(xlam.ne.olam) then
        olam = xlam
        b5  = (33-2*5)/pi/12
        b4  = (33-2*4)/pi/12
        b3  = (33-2*3)/pi/12
        if(nloop.eq.1) then
           bp5=0
           bp4=0
           bp3=0
        elseif(nloop.eq.2) then
           bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
           bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
           bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        endif 
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     $        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     $        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
        if( q .gt. xmb ) then
          nf = 5
        elseif( q .gt. xmc ) then
          nf = 4
        else
          nf = 3
        endif
      endif
      if    ( nf .eq. 5 ) then
        alfas_clu = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        alfas_clu = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) +
     $       c45 )
      elseif( nf .eq. 3 ) then
        alfas_clu = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) +
     $       c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        stop
      endif
      return
      end

c-------------------------------------------
c Program to calculate as with nf flavours
c as a function of lambda with nf flavours
c
      function alfa(q,xlam,nloop,nf)
      data pi/3.14159/
      xlp = float(nloop-1)
      b  = (33-2*nf)/pi/12
      bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf)  * xlp
      t = 2 * log( q/xlam )
      xlt = log( t )
      alfa = 1/(b * t) -  bp/(b * t)**2 * xlt
      return
      end

c----------------------------------------------------------
c Program to get lambda_nf from as_nf at the scale q
c
      function xlambd(as,q,nloop,nf)
      implicit real * 8 (a-h,o-z)
      data pi/3.14159d0/
      xlp = float(nloop-1)
      b  = (33-2*nf)/pi/12
      bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf) * xlp
      t  = 1/b/as
c-----------------------------------------------------------
c Solve the equation
c
    1 xlt = log(t)
      ot = t
c-----------------------------------------------------------
c Solve the equation
c Value and Derivative of alfa with respect to t
c
      as0  = 1/b/t - bp*xlt/(b*t)**2
      as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3
      t  = (as-as0)/as1 + t
      if(abs(ot-t)/ot.gt..00001)goto 1
      xlambd = q/exp(t/2)
      return
      end


      SUBROUTINE MWARN(ROUT)
      CHARACTER*(*) ROUT
      WRITE(*,*) '***********************************************'
      WRITE(*,*) '***** WARNING CALLED FROM ROUTINE ',ROUT,':'
      END

