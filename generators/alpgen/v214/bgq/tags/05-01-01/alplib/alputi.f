C**********************************************************************
C    SIMPLE HISTOGRAMMING PACKAGE --  SIMPLIFIED VERSION OF HBOOK
C    BY Michelangelo Mangano    NOVEMBER 1988
C    LAST REVISED NOVEMBER 9, 1988
C    LAST REVISED JUNE 12, 1989  (ADD SCATTER PLOTS)
C    LAST REVISED oct 1990 (Add multi-plots on one page, routines MULTITOP,
C         			MTFILL,...)
C**********************************************************************
C
C Fills up to 100 histograms with up to 100 bins. 
C Gives a data file (to be specified in the calling program by assigning 
C a file name to unit 98) and a topdrawer file (to be specified in the 
C calling program by assigning a file name to unit 99).
C
C INITIALIZATION:
C Call once INIHIST; this just resets a few counters and logicals
C Call MBOOK(N,'TITLE',DEL,XMIN,XMAX) for each histogram to be booked.
C N (an integer) is the label of the histogram;
C 'TITLE' is the name of the histogram (no more then 100 characters);
C DEL (real*4) is the bin size;
C XMIN (real*4) is the lower limit of the first bin;
C XMAX (real*4)is the upper limit of the last  bin
C Example:
C      call mbook(2,'pt distribution',1.,10,70)
C This call initializes histogram number 2, called 'pt distribution';
C The bin size will be 1. (possibly GeV, if that's what you want), the
C first bin being  10<x<11. and the last one being 69.<x<70
C
C FILLING:
C When it's time, call MFILL(N,X,Y); this will add Y (real*4) to the bin 
C in which X (real*4) happens to be, within histogram N. 
C
C PLAYING AROUND:
C At the end of the day you may want to sum, divide, cancel, etc.etc.
C various histograms (bin by bin). Then you call MOPERA(I,'O',J,K,X,Y). 
C The 1-character string O can take the following values:
C +  : sums       X*(hist I) with Y*(hist J) and puts the result in hist K;
C -  : subtracts  X*(hist I) with Y*(hist J) and puts the result in hist K;
C *  : multiplies X*(hist I) with Y*(hist J) and puts the result in hist K;
C /  : divides    X*(hist I) with Y*(hist J) and puts the result in hist K;
C F  : multiplies hist I by the factor X, and puts the result in hist K;
C R  : takes the square root of  hist  I, and puts the result in hist K;if
C      the value at a given bin is less than or equal to 0, puts 0 in K
C S  : takes the square      of  hist  I, and puts the result in hist K;
C L  : takes the log_10 of  hist  I, and puts the result in hist K; if the
C      value at a given bin is less than or equal to 0, puts 0 in K
C M  : statistical analysis; if I contains the weights (let's say WGT),
C      J contains variable times weight (F*WGT) and K contains the
C      variable squared times the weight (F**2*WGT), then, after using 'M',
C      J will contain the average value of the variable <F> and K will 
C      contain the sigma of the average: sigma=sqrt(<F**2>-<F>**2).
C      If WGT=1. for all the entries, then it is enough to put I=J, and
C      it is not necessary to book a hist with the weights.
C V  : estimates errors for vegas evaluation of differential distributions.
C      Fill I with the values of
C      the functions do integrate times the Vegas weight (fun*wgt); fill
C      J with fun**2*wgt; then K will contain an estimate of the error
C      of the integration. Putting X=1/(#of iterations) performs the 
C      avegare over the iterations, and gives the right normalization to 
C      the differential distribution, I, and to the errors, K. J stays the same.
C
C FINAL ACCOUNTING:
C Now we can finalize our histograms; MFINAL(N) will calculate the integral
C of the histogram N, the mean value of the X variable and its RMS.
C If we now want to renormalize the hist's, we can call MNORM(N,X), which
C will normalize the integral to X  -- CAUTION: do not call MNORM before
C MFINAL, it will blow up.
C
C OUTPUT:
C To get a .dat file containing the values of the histograms, together with
C some information (like integral, mean values, etc.etc.) call MPRINT(N),
C for each hist N that you want in the .dat file. Before the call to MPRINT
C you want to open unit 98 and give it a name:                       
C     OPEN(UNIT=98,NAME='NAME.DAT',STATUS='NEW')
C If you want a topdrawer file with a plot of the hist values, call 
C MTOP(N,M,'X','Y','SCALE'). The points of the plot will be taken from histogram
C N, the error bars from histogram M. 'SCALE', character*(*), determines
C the scale for y, logarithmic or linear (SCALE=LOG,LIN). 
C If you do not want error bars, keep
C a histogram of zeros, or just call a hist that had not been booked.
C X will appear as a 'bottom title', and Y will appear as a 'left title'.
C The top title is by default the name of the histogram itself.
C A little box below the plot will contain some information on the plot
C itself. Before calling MTOP,
C     OPEN(UNIT=99,NAME='NAME.TOP',STATUS='NEW')
C--------------------------------------------------------------------------
C
C  COMMON/HISTO/  Histogram N
C                           
C   BOOK(N),      Three-letter character-string: 'NO' if histogram was not 
C		  Booked, 'YES' otherwise
C   TITLE(N),     Title of the histogram
C
C   HMIN(N),      Min value of x range
C   HMAX(N),      Max value of x range
C   HDEL(N),      Bin width
C   NBIN(N),      Total number of bins
C   USCORE(N),    Total integral of underscores with x < HMIN(N)
C   OSCORE(N),    Total integral of onderscores with x > HMAX(N)
C   IUSCORE(N),   Number of entries with x < HMIN(N)
C   IOSCORE(N),   Number of entries with x > HMAX(N)
C   IENT(N),      Total number of entries within x range HMIN(N)<x<HMAX(N)
C   HINT(N),      Integral of the histogram within HMIN(N)<x<HMAX(N)
C   HAVG(N),      Average value of x, weighted over the x range of the histo
C   HSIG(N),      Quadratic dispersion of x around the average
C   HIST(N,L),    Value of bin L-th
C   XHIS(N,L),    Central x value of bin L-th
C   IHIS(N,L),    Number of entries within bin L-th
C   NHIST         Total number of booked histograms
C


      SUBROUTINE INIHIST
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST
      COMMON/HISTO2/HIST2(10,100,100),XHIS2(10,100),YHIS2(10,100)
     &,HDEL2(2,10),XPROJ(10,100),YPROJ(10,100),HMIN2(2,10)
     &,HMAX2(2,10),NBIN2(2,10),IHIS2(10,100,100),IOSCORE2(10)
     &,IENT2(10),NHIST2,HAVG2(2,10),HINT2(10),HSIG2(2,10),BOOK2(10),
     & TITLE2(10)
      CHARACTER TITLE*25,BOOK*3
      CHARACTER TITLE2*25,BOOK2*3
      double precision hist
      DO 1, I=1,200             
   1  BOOK(I)=' NO'
      DO 2 I=1,10
   2  BOOK2(I)=' NO'
      END  
    
      SUBROUTINE MBOOK(N,TIT,DEL,XMIN,XMAX)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST
      CHARACTER TITLE*25,BOOK*3
      CHARACTER*(*) TIT
      double precision hist
      NHIST = MAX(N,NHIST)
      IF(BOOK(N)(1:1).EQ.'Y') THEN
	 CALL MBWARN('MBOOK')
         WRITE(*,*) 'Histogram',N,TITLE(N),' already in use. '
         WRITE(*,*) 'superseded by ',TIT
      ENDIF
      BOOK(N) = 'YES'
      TITLE(N) = ' '//TIT
1     HDEL(N) = DEL
      NBIN(N) = INT((XMAX-XMIN)/DEL)
      IF(NBIN(N).GT.100) THEN
	WRITE(*,*) 'TOO MANY BINS (',NBIN(N),') REQUIRED IN HIST ',N
	WRITE(*,*) 'RE-ENTER BIN SIZE DELTA (OLD BIN = ',DEL,' ):'
	READ(*,*) DEL
	GO TO 1
      ENDIF
      HMIN(N) = XMIN
      HMAX(N) = NBIN(N)*DEL+XMIN
c      IF(HMAX(N).NE.XMAX) THEN
c	 CALL MBWARN('MBOOK')
c         WRITE(*,*)
c     #'Histogram ', TIT, ' Change of upper limit:',xmax,'-->',HMAX(N)
c      ENDIF
      IENT(N) = 0
      IUSCORE(N) = 0
      IOSCORE(N) = 0
      USCORE(N) = 0
      OSCORE(N) = 0
      HAVG(N) = 0
      HINT(N) = 0
      HSIG(N) = 0
      DO I=1,NBIN(N)
         XHIS(N,I)=HMIN(N)+HDEL(N)*(FLOAT(I)-0.5)
         IHIS(N,I)=0
         HIST(N,I)=0
      ENDDO
      END

      SUBROUTINE MBOOK_BIN(N,TIT,XEDGE)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST
      COMMON/HISTO_BIN/XLOW(200,100),XHIGH(200,100),XDEL(200,100)
      CHARACTER TITLE*25,BOOK*3
      REAL*4 XEDGE(100)
      CHARACTER*(*) TIT
      double precision hist
      NHIST = MAX(N,NHIST)
      IF(BOOK(N)(1:1).EQ.'Y') THEN
	 CALL MBWARN('MBOOK')
         WRITE(*,*) 'Histogram',N,TITLE(N),' already in use. '
         WRITE(*,*) 'superseded by ',TIT
      ENDIF
      BOOK(N) = 'YES'
      TITLE(N) = ' '//TIT
      HDEL(N) = -1
      HMIN(N) = XEDGE(1)
      IENT(N) = 0
      IUSCORE(N) = 0
      IOSCORE(N) = 0
      USCORE(N) = 0
      OSCORE(N) = 0
      HAVG(N) = 0
      HINT(N) = 0
      HSIG(N) = 0
      DO I=1,99
        DEL=XEDGE(I+1)-XEDGE(I)
        IF(DEL.GT.0) THEN
          XDEL(N,I)=DEL
          XHIS(N,I)=XEDGE(I)+0.5*DEL
          XLOW(N,I)=XEDGE(I)
          XHIGH(N,I)=XEDGE(I+1)
          HMAX(N) = XEDGE(I+1)
          NBIN(N) = I
          IHIS(N,I)=0
          HIST(N,I)=0
        ENDIF
      ENDDO
      END


      SUBROUTINE MFILL(N,X,Y)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      XI=((X-HMIN(N))/HDEL(N))+1
      I=INT(XI)
      IF(I.LT.1) THEN
         USCORE(N) = USCORE(N) + Y
         IUSCORE(N) = IUSCORE(N) + 1
      ELSEIF(I.GT.NBIN(N)) THEN
         OSCORE(N) = OSCORE(N) + Y
         IOSCORE(N) = IOSCORE(N) + 1
      ELSE
         IENT(N)=IENT(N)+1
         IHIS(N,I)=IHIS(N,I)+1
         HIST(N,I)=HIST(N,I)+dble(Y)
      ENDIF
      END

      SUBROUTINE MFILL_BIN(N,X,Y)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      COMMON/HISTO_BIN/XLOW(200,100),XHIGH(200,100),XDEL(200,100)
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      DO I=1,NBIN(N)
        IF(X.GE.XLOW(N,I).AND.X.LT.XHIGH(N,I)) THEN 
          IF(I.LT.1) THEN
            USCORE(N) = USCORE(N) + Y
            IUSCORE(N) = IUSCORE(N) + 1
          ELSEIF(I.GT.NBIN(N)) THEN
            OSCORE(N) = OSCORE(N) + Y
            IOSCORE(N) = IOSCORE(N) + 1
          ELSE
            IENT(N)=IENT(N)+1
            IHIS(N,I)=IHIS(N,I)+1
            HIST(N,I)=HIST(N,I)+dble(Y)
          ENDIF
          GOTO 999
        ENDIF
      ENDDO
 999  CONTINUE
      END


      SUBROUTINE MINTEG(NIN,NOUT,IDIR,IPOW)
C If IPOW=1 performs the integral of the distribution contained in histogram
C NIN up to the value specified by the abscissa (if IDIR=1) or from this
C value on (if IDIR=-1). The resulting integral distribution is put into 
C NOUT, which is automatically booked if NOUT.ne.NIN .  Choosing IPOW=2
C the routine will return the square root of the integral of the squares,
C as is required, for example, for the propagation of the mean quadratic error
C of a given distribution. Overscores and underscores are included.
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      CHARACTER*14  C
      double precision hist
      DIMENSION C(2) 
      DATA C/' INTEG BELOW X',' INTEG ABOVE X'/
      IF(BOOK(NIN)(1:1).NE.'Y') RETURN
      M = NBIN(NIN)                                           
      I = (IDIR + 3)/2
      IF(NOUT.NE.NIN.AND.BOOK(NOUT)(1:1).NE.'Y') THEN
      	CALL MBOOK(NOUT,TITLE(NIN)//C(I), 
     &                HDEL(NIN),HMIN(NIN),HMAX(NIN))
      ENDIF
      IF(IDIR.EQ.1) THEN
        HIST(NOUT,1) = dble(SUMPOW(real(hist(nin,1)),USCORE(NIN),IPOW))
         IHIS(NOUT,1) = IHIS(NIN,1) + IUSCORE(NIN)
         XHIS(NOUT,1) = XHIS(NIN,1) + HDEL(NIN)/2
         DO L=2,M                      
           HIST(NOUT,L) = dble(SUMPOW(real(HIST(NIN,L)),real(HIST(NOUT,L
     $          -1)),IPOW))
            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L-1) 
            XHIS(NOUT,L) = XHIS(NIN,L) + HDEL(NIN)/2
         ENDDO
         OSCORE(NOUT) = dble(SUMPOW(OSCORE(NIN),real(HIST(NIN,M)),IPOW))
         IOSCORE(NOUT) = IOSCORE(NIN) + IHIS(NIN,M)
      ELSEIF(IDIR.EQ.-1) THEN
         HIST(NOUT,M) = dble(SUMPOW(real(HIST(NIN,M)),OSCORE(NIN),IPOW))
         IHIS(NOUT,M) = IHIS(NIN,M) + IOSCORE(NIN)
         XHIS(NOUT,M) = XHIS(NIN,M) - HDEL(NIN)/2
         DO L=M-1,1,-1                        
           HIST(NOUT,L) = dble(SUMPOW(real(HIST(NIN,L)),real(HIST(NOUT,L
     $          +1)),IPOW))
            IHIS(NOUT,L) = IHIS(NIN,L) + IHIS(NOUT,L+1)
            XHIS(NOUT,L) = XHIS(NIN,L) - HDEL(NIN)/2
         ENDDO
         USCORE(NOUT) = dble(SUMPOW(USCORE(NIN),real(HIST(NIN,1)),IPOW))
         IUSCORE(NOUT) = IUSCORE(NIN)+IHIS(NIN,1)
      ELSE                                 
         CALL MBWARN('MINTEG')
         WRITE(*,*) 'Wrong idir in minteg: OPERATION NOT PERFORMED'
         RETURN
      ENDIF
      END

      FUNCTION SUMPOW(X,Y,IPOW)
      IF(IPOW.EQ.1) THEN
         SUMPOW = X + Y
      ELSEIF(IPOW.EQ.2) THEN
         SUMPOW = SQRT(X**2+Y**2)
      ELSEIF(IPOW.EQ.0) THEN
         CALL MBWARN('SUMPOW')
         WRITE(*,*)'Error: IPOW=0 not allowed in SUMPOW'
      ELSE
         SUMPOW = (X**IPOW+Y**IPOW)**(1./IPOW)
      ENDIF
      END

      SUBROUTINE MOPERA(I,OPER,J,K,X,Y)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      CHARACTER OPER*1
      double precision hist
      IF(NBIN(I).NE.NBIN(J).AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &    '*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M'.OR.OPER.EQ.'A')) THEN
	  CALL MBWARN('MOPERA')
          WRITE(*,*) I,J                               
  20      FORMAT(' ****** INCOMPATIBLE OPERATION HIST ',I2,' &',I2,
     &    '*******'/)
          RETURN
      ENDIF
      IF(BOOK(I)(1:1).NE.'Y'.OR.BOOK(J)(1:1).NE.'Y') RETURN
      IF(BOOK(K)(1:1).NE.'Y') 
     &  CALL MBOOK(K,TITLE(I),HDEL(I),HMIN(I),HMAX(I))
      IF(OPER.EQ.'E') THEN
c If I contains the accumulated weights, J the accumulated squares of the
c weights and IHIS(J,1) the number of accumulated entries, 'E' will add
c the average value of I to K and will put in J the quadratic dispersion.
         IF(IHIS(J,1).NE.0) THEN
            XXX = 1./IHIS(J,1)
         ELSE
            XXX = 0
         ENDIF
         DO L=1,NBIN(I)
            XSUM   = HIST(I,L)
            XSUMSQ = HIST(J,L)
            HIST(K,L)=HIST(K,L) + XXX*HIST(I,L)
            HIST(J,L)=XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
         ENDDO
         IENT(K)=IENT(K)+IENT(I)
         XSUM = USCORE(I)
         XSUMSQ = USCORE(J)
         USCORE(K) = USCORE(K)+XXX*XSUM
         USCORE(J) = XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
         XSUM = OSCORE(I)
         XSUMSQ = OSCORE(J)
         OSCORE(K) = OSCORE(K)+XXX*XSUM
         OSCORE(J) = XXX*SQRT(ABS(XSUMSQ-XSUM**2*XXX))
      ELSEIF(OPER.EQ.'Q') THEN
         DO L=1,NBIN(I)
            HIST(K,L) = SQRT(HIST(J,L)**2+HIST(I,L)**2)
         ENDDO
         USCORE(K) = SQRT(USCORE(J)**2+USCORE(I)**2)
         OSCORE(K) = SQRT(OSCORE(J)**2+OSCORE(I)**2)
      ELSEIF(OPER.EQ.'A') THEN
         DO L=1,NBIN(I)
            HIST(J,L) = HIST(J,L) + HIST(I,L)
            IHIS(J,L) = IHIS(J,L) + IHIS(I,L)
            HIST(K,L) = HIST(K,L) + HIST(I,L)**2
            IHIS(K,L) = IHIS(K,L) + 1
            IENT(K) = IENT(K)+1
            HIST(I,L) = 0
            IHIS(I,L) = 0
         ENDDO
         IENT(J) = IENT(J)+IENT(I)
         IUSCORE(J) = IUSCORE(J) + IUSCORE(I)
         USCORE(J) = USCORE(J) + USCORE(I)
         IOSCORE(J) = IOSCORE(J) + IOSCORE(I)
         OSCORE(J) = OSCORE(J) + OSCORE(I)
         IUSCORE(K) = IUSCORE(K) + 1
         USCORE(K) = USCORE(K) + USCORE(I)**2
         IOSCORE(K) = IOSCORE(K) + 1
         OSCORE(K) = OSCORE(K) + OSCORE(I)**2
         IENT(I) = 0
         IUSCORE(I) = 0
         IOSCORE(I) = 0
         USCORE(I) = 0
         OSCORE(I) = 0
      ELSE
        DO L=1,NBIN(I)
      	IF(OPER.EQ.'+') THEN
       	  HIST(K,L)=X*HIST(I,L) + Y*HIST(J,L)
      	ELSEIF(OPER.EQ.'-') THEN
      	  HIST(K,L)=X*HIST(I,L) - Y*HIST(J,L)
      	ELSEIF(OPER.EQ.'*') THEN
      	  HIST(K,L)=X*HIST(I,L) * Y*HIST(J,L)
      	ELSEIF(OPER.EQ.'/') THEN
          IF(Y.EQ.0..OR.HIST(J,L).EQ.0.) THEN
            HIST(K,L)=0.
          ELSE
            HIST(K,L)=X*HIST(I,L) / (Y*HIST(J,L))
          ENDIF
       	ELSEIF(OPER.EQ.'F') THEN
      	  HIST(K,L)=X*HIST(I,L)
      	ELSEIF(OPER.EQ.'R') THEN
          IF(HIST(I,L).GT.0.) THEN
            HIST(K,L)=X*SQRT(HIST(I,L))
          ELSE                           
            HIST(K,L)=0.
          ENDIF
      	ELSEIF(OPER.EQ.'S') THEN
          HIST(K,L)=X*HIST(I,L)**2
      	ELSEIF(OPER.EQ.'L') THEN  
          IF(HIST(I,L).EQ.0..OR.J.EQ.0.) THEN
             HIST(K,L)=0.
           ELSE
             HIST(K,L)=X*LOG10(Y*HIST(I,L))
           ENDIF
      	ELSEIF(OPER.EQ.'M') THEN
           IF(I.NE.J) XNORM=HIST(I,L)
           IF(I.EQ.J) XNORM=FLOAT(IHIS(J,L))
           IF(XNORM.NE.0.) THEN
             XAVG=HIST(J,L)/XNORM
             HIST(K,L)=
     &       SQRT(ABS(-XAVG**2+HIST(K,L)/XNORM)/FLOAT(IHIS(I,L)))
             HIST(J,L)=XAVG 
           ELSE 
             HIST(K,L)=0.
             HIST(J,L)=0.
           ENDIF
      	ELSEIF(OPER.EQ.'V') THEN                 
           XAVG=HIST(I,L)*X
           XSQAVG=HIST(J,L)*X
           XNORM=FLOAT(IHIS(I,L))*X
           IF(XNORM.NE.0.) THEN
              HIST(K,L)=SQRT(ABS(XSQAVG-XAVG**2)/XNORM)
              HIST(I,L)=XAVG
           ELSE  
              HIST(K,L)=0.
           ENDIF 
      	ELSE 
	 CALL MBWARN('MOPERA')
         WRITE(*,*) OPER
   5     FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
         RETURN
        ENDIF
        ENDDO
      	IF(OPER.EQ.'+') THEN
       	  USCORE(K)=X*USCORE(I) + Y*USCORE(J)  
       	  OSCORE(K)=X*OSCORE(I) + Y*OSCORE(J)  
      	ELSEIF(OPER.EQ.'-') THEN     
      	  USCORE(K)=X*USCORE(I) - Y*USCORE(J)
      	  OSCORE(K)=X*OSCORE(I) - Y*OSCORE(J)
      	ELSEIF(OPER.EQ.'*') THEN     
      	  USCORE(K)=X*USCORE(I) * Y*USCORE(J)
      	  OSCORE(K)=X*OSCORE(I) * Y*OSCORE(J)
      	ELSEIF(OPER.EQ.'/') THEN     
          IF(Y.EQ.0..OR.USCORE(J).EQ.0.) THEN
            USCORE(K)=0.
          ELSE
            USCORE(K)=X*USCORE(I) / (Y*USCORE(J))
          ENDIF
          IF(Y.EQ.0..OR.OSCORE(J).EQ.0.) THEN
            OSCORE(K)=0.
          ELSE
            OSCORE(K)=X*OSCORE(I) / (Y*OSCORE(J))
          ENDIF
       	ELSEIF(OPER.EQ.'F') THEN
      	  USCORE(K)=X*USCORE(I)
      	  OSCORE(K)=X*OSCORE(I)
      	ELSEIF(OPER.EQ.'R') THEN
          IF(USCORE(I).GT.0.) THEN
            USCORE(K)=X*SQRT(USCORE(I))
          ELSE                           
            USCORE(K)=0.
          ENDIF     
          IF(OSCORE(I).GT.0.) THEN
            OSCORE(K)=X*SQRT(OSCORE(I))
          ELSE                           
            OSCORE(K)=0.
          ENDIF     
      	ELSEIF(OPER.EQ.'S') THEN
          USCORE(K)=X*USCORE(I)**2
          OSCORE(K)=X*OSCORE(I)**2
      	ELSEIF(OPER.EQ.'L') THEN  
          IF(USCORE(I).EQ.0..OR.J.EQ.0.) THEN
             USCORE(K)=0.
           ELSE
             USCORE(K)=X*LOG10(Y*USCORE(I))
           ENDIF                         
          IF(OSCORE(I).EQ.0..OR.J.EQ.0.) THEN
             OSCORE(K)=0.
           ELSE
             OSCORE(K)=X*LOG10(Y*OSCORE(I))
           ENDIF                         
        ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE MCOPY(NIN,NOUT)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
C
      IF(BOOK(NIN)(1:1).NE.'Y') RETURN
      IF(BOOK(NOUT)(1:1).NE.'Y') 
     &  CALL MBOOK(NOUT,TITLE(NIN),HDEL(NIN),HMIN(NIN),HMAX(NIN))
      HDEL(NOUT)=HDEL(NIN)
      HMIN(NOUT)=HMIN(NIN)
      HMAX(NOUT)=HMAX(NIN)
      NBIN(NOUT)=NBIN(NIN)
      IENT(NOUT)=IENT(NIN)
      IUSCORE(NOUT)=IUSCORE(NIN)
      IOSCORE(NOUT)=IOSCORE(NIN)
      USCORE(NOUT)=USCORE(NIN)
      OSCORE(NOUT)=OSCORE(NIN)
      HAVG(NOUT)=HAVG(NIN)
      HINT(NOUT)=HINT(NIN)
      HSIG(NOUT)=HSIG(NIN)
      DO I=1,NBIN(NOUT)
        HIST(NOUT,I)=HIST(NIN,I)
        XHIS(NOUT,I)=XHIS(NIN,I)
        IHIS(NOUT,I)=IHIS(NIN,I)
      ENDDO
      END

      SUBROUTINE MZERO(N)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      BOOK(N)='RES'
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0.
      HINT(N)=0.
      DO 1 I=1,NBIN(N)
   1  HIST(N,I)=0.
      END                   

      SUBROUTINE MSAVE
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST
      DIMENSION HISTS(200,100),XHISS(200,100),HDELS(200),HMINS(200)
     $     ,HMAXS(200),USCORES(200),OSCORES(200),NBINS(200),IHISS(200
     $     ,100),IUSCORES(200),IOSCORES(200),IENTS(200),HAVGS(200)
     $     ,HINTS(200),HSIGS(200),BOOKS(200),TITLES(200)
      CHARACTER TITLE*25,BOOK*3,TITLEL*25,BOOKL*3
      double precision hist
      SAVE 
      DO N=1,200
         DO I=1,100
            HISTS(N,I)=HIST(N,I)
            XHISS(N,I)=XHIS(N,I)
         ENDDO
         HDELS(N)=HDEL(N)
         HMINS(N)=HMIN(N)
         HMAXS(N)=HMAX(N)
         USCORES(N)=USCORE(N)
         OSCORES(N)=OSCORE(N)
         NBINS(N)=NBIN(N)
         HSIGS(N)=HSIG(N)
         IENTS(N)=IENT(N)     
         IUSCORES(N)=IUSCORE(N)
         IOSCORES(N)=IOSCORE(N)      
         HAVGS(N)=HAVG(N)         
         HINTS(N)=HINT(N)   
         NHISTS=NHIST
      ENDDO
      RETURN
      ENTRY MRESTORE
      DO N=1,200
         DO I=1,100
            HIST(N,I)=HISTS(N,I)
            XHIS(N,I)=XHISS(N,I)
         ENDDO
         HDEL(N)=HDELS(N)
         HMIN(N)=HMINS(N)
         HMAX(N)=HMAXS(N)
         USCORE(N)=USCORES(N)
         OSCORE(N)=OSCORES(N)
         NBIN(N)=NBINS(N)
         HSIG(N)=HSIGS(N)
         IENT(N)=IENTS(N)     
         IUSCORE(N)=IUSCORES(N)
         IOSCORE(N)=IOSCORES(N)      
         HAVG(N)=HAVGS(N)         
         HINT(N)=HINTS(N)   
         NHISTS=NHISTS
      ENDDO
      END

      SUBROUTINE MDUMP(N)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      DIMENSION HISTL(100),XHISL(100),IHISL(100)
      CHARACTER TITLE*25,BOOK*3,TITLEL*25,BOOKL*3
      double precision hist
      SAVE HISTL,XHISL,HDELL,HMINL,HMAXL,USCOREL,OSCOREL,
     &NBINL,IHISL,IUSCOREL,IOSCOREL,IENTL,HAVGL,HINTL,HSIGL,
     &BOOKL,TITLEL
      BOOKL=BOOK(N)                                                           
      TITLEL=TITLE(N)                   
      HDELL=HDEL(N)  
      HMINL=HMIN(N)
      HMAXL=HMAX(N)
      USCOREL=USCORE(N)
      OSCOREL=OSCORE(N)
      NBINL=NBIN(N)
      HSIGL=HSIG(N)
      IENTL=IENT(N)     
      IUSCOREL=IUSCORE(N)
      IOSCOREL=IOSCORE(N)      
      HAVGL=HAVG(N)         
      HINTL=HINT(N)   
      DO 1 I=1,NBIN(N)
      XHISL(I)=XHIS(N,I)
      IHISL(I)=IHIS(N,I)
   1  HISTL(I)=HIST(N,I)
      RETURN
      ENTRY MFETCH(N)
      BOOK(N)=BOOKL                   
      TITLE(N)=TITLEL                   
      HDEL(N)=HDELL  
      HMIN(N)=HMINL
      HMAX(N)=HMAXL
      USCORE(N)=USCOREL
      OSCORE(N)=OSCOREL
      NBIN(N)=NBINL
      HSIG(N)=HSIGL
      IENT(N)=IENTL     
      IUSCORE(N)=IUSCOREL
      IOSCORE(N)=IOSCOREL      
      HAVG(N)=HAVGL         
      HINT(N)=HINTL   
      DO 2 I=1,NBINL
      XHIS(N,I)=XHISL(I)
      IHIS(N,I)=IHISL(I)
   2  HIST(N,I)=HISTL(I)
      END               

      SUBROUTINE MRESET(N)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      BOOK(N)='RES'
      END

      SUBROUTINE PUTTAG(J,NAME)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
c Per marcare un istogramma
      CHARACTER * (*) NAME, TAG
      double precision hist
      BOOK(J) = NAME
      RETURN
      ENTRY GETTAG(J,TAG)
      TAG = BOOK(J)
      END

      SUBROUTINE MFINAL(N)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      AVG=0
      XIN=0                                  
      SIG=0
      IF=0
      DO J=1,NBIN(N)
         X=HIST(N,J)
 	 AVG=AVG+X*XHIS(N,J)
         XIN=XIN+X
	 IF(X.NE.0) IF=1
      ENDDO             
      IF(XIN.EQ.0) GO TO 10
      AVG = AVG/XIN
      DO J=1,NBIN(N)
         SIG=HIST(N,J)*(XHIS(N,J)-AVG)**2+SIG
      ENDDO
      SIG=SQRT(ABS(SIG/XIN))
 10   CONTINUE
      HINT(N) = XIN
      HAVG(N) = AVG
      HSIG(N) = SIG
      IF(IF.EQ.0) BOOK(N)='RES'
      END               

      SUBROUTINE MNORM(N,X,IOPT)
C NORMALIZES PLOT TO INTEGRAL=X.
C IOPT=0: OVER(UNDER)FLOWS NOT INCLUDED
C IOPT=1: OVERFLOWS  INCLUDED
C IOPT=2: UNDERFLOWS INCLUDED
C IOPT=3:  ALLFLOWS  INCLUDED
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3
      double precision hist
      IF(BOOK(N)(1:1).NE.'Y')RETURN
      IF(HINT(N).EQ.0.) THEN
	CALL MBWARN('MNORM')
	WRITE(*,*)' INTEGRAL HIST ',N,' IS ZERO: CANNOT RENORMALIZE'
	RETURN               
      ELSE
	DEN=HINT(N)
        IF(IOPT.EQ.1) THEN
          DEN=DEN+OSCORE(N)
        ELSEIF(IOPT.EQ.2) THEN
          DEN=DEN+USCORE(N)
        ELSEIF(IOPT.EQ.3) THEN
          DEN=DEN+USCORE(N)+OSCORE(N)
        ENDIF                
        Y=X/DEN
      ENDIF
      DO 1, I=1,NBIN(N)
    1 HIST(N,I)=HIST(N,I)*Y
      HINT(N)=X            
      OSCORE(N)=OSCORE(N)*Y
      USCORE(N)=USCORE(N)*Y
      END                  

      SUBROUTINE MTOP(N,M,BTIT,LTIT,SCALE)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3,NOW*24
      CHARACTER*(*) LTIT,BTIT,SCALE
      double precision hist
      DATA INI/0/
      IF(INI.EQ.0) THEN
         call USRTIME(NOW)
      INI=1
      ENDIF
      IF(BOOK(N)(:1).NE.'Y') RETURN
      WRITE(99,100) now,TITLE(N),BTIT,LTIT,SCALE,HMIN(N),HMAX(N),N
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 9.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET FONT DUPLEX '/1X, 
     &' title 10.3 9 angle -90 ','"ALPGEN, ',A,'"'/1X, 
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A,'"',/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ',/1X,
     &' ( BEGIN HIST=',I3)
      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      WRITE(99,'(3X,F10.4,2(2X,E15.4))')  
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
    1 CONTINUE
      WRITE(99,200) N
  200 FORMAT(' ( END HIST=',I3,/,'   PLOT')
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &     ,IOSCORE(N),USCORE(N),OSCORE(N)
  300 FORMAT( /1x,                               
     &' BOX 6.25 0.9 SIZE 7.5 1.2'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' SET FONT DUPLEX '/1X,
     &' TITLE 2.8 1.2 "INT =',1PE10.3,'   AVG =',1PE10.3,
     &             '   RMS =',1PE10.3,'"',/1X,
     &' TITLE 2.8 0.9 "Entries =',I8,2x,'Undersc =',I6,2X
     &                                 ,'Oversc =',I6,'"',/1X,
     &' TITLE 2.8 0.6 "Intgr ufloat=',1PE10.3,'  Intgr ofloat=',
     &      1PE10.3,'"',/1X
     $     ,' SET TITLE SIZE -2')                            
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END
 
            
      SUBROUTINE MULTITOP(NH,NE,N,M,BTIT,LTIT,SCA)
      COMMON/HISTO/HIST(200,100),XHIS(200,100),HDEL(200),HMIN(200)
     &,HMAX(200),USCORE(200),OSCORE(200)
     &,NBIN(200),IHIS(200,100),IUSCORE(200),IOSCORE(200)
     &,IENT(200),HAVG(200),HINT(200),HSIG(200),BOOK(200),TITLE(200)
     &,NHIST                                                      
      CHARACTER TITLE*25,BOOK*3,SCALE*3,NOW*24
      CHARACTER*(*) LTIT,BTIT,SCA
      double precision hist
      CHARACTER*7  PLOT(4)
      DATA PLOT/'SOLID','DASHES','DOTS','DOTDASH'/
C  PLOT SIZE, CORNERS
      DATA WIDTH,HEIGHT/11.5,8.5/,XCORN,YCORN/1.5,1./
C  PLOT VERSUS TEXT FRACTION                  
      DATA XPFRAC,YPFRAC/0.75,0.75/,XTFRAC,YTFRAC/0.25,0.25/
C  DEFAULT SIZES                                           
      REAL LAB0,LABS               
      DATA TIT0,LAB0,TIC0/-3,-3,0.06/
      DATA INI/0/                                          
      IF(INI.EQ.0) THEN
         CALL USRTIME(NOW)
         IFRAME=0        
         INI=1         
      ENDIF
      IF(SCA.EQ.'REF') THEN
	IFRAME=0
	RETURN
      ENDIF
      IF(BOOK(NH)(:1).NE.'Y') RETURN
      IFRMAX=N*M         
      IFRAME=IFRAME+1
      IF(IFRAME.GT.IFRMAX.OR.N.NE.NOLD.OR.M.NE.MOLD) THEN
      	IFRAME=1
        WRITE(99,202)   
         WRITE(99,71) NOW
 71      FORMAT(' title 12 9 angle -90 SIZE 1.5','"ALPGEN, ',A,'"'/1X)
  1     FORMAT(' SET FONT DUPLEX',/,'  SET TITLE SIZE 2',/,
     +       ' TITLE 12.8 9 ANGLE -90 ','" MLM   ',I2,'-',I2,1X,I2,':'
     $       ,I2,'"')
      ENDIF                                
      IF(IFRAME.EQ.1) THEN
    	I=1
	J=1
      ELSEIF(IFRAME.LE.IFRMAX) THEN
	IF(I.LE.N) I=I+1
        IF(I.GT.N) THEN
		I=1
		J=J+1
	ENDIF
      ENDIF
      IF(N.EQ.NOLD) GO TO 10
      NS=N-1
      XD=WIDTH/FLOAT(N)
      SRED=SQRT(FLOAT(N*M))
      TITS=TIT0/SRED          
      LABS=LAB0/SRED
      TICS=TIC0/SRED
      XTIT0=0.55*XPFRAC*XD
      NOLD=N            
10    IF(M.EQ.MOLD) GO TO 20
      YD=HEIGHT/FLOAT(M)
      YTIT0=0.06*YD
      MOLD=M        
20    CONTINUE
      XL=(I-1)*XD + XCORN
      YL=(M-J)*YD + YCORN
      XU=XL+XD*XPFRAC
      YU=YL+YD*YPFRAC        
      IP=0
      FMN=MAX(HINT(NH)*NBIN(NH),1.E12)
      FMX=-FMN                        
      XMX=0.
      DO IBIN=1,NBIN(NH)
	X=HIST(NH,IBIN)
      	IF(X.NE.0.) FMX=MAX(FMX,X)
      	IF(X.NE.0.) FMN=MIN(FMN,X)
	XMX=MAX(XMX,ABS(X))
      ENDDO                
      IF(XMX.EQ.0.) GO TO 203
      SCALE=SCA
50    IF(SCALE.EQ.'LIN') THEN
	IF(FMN.GE.0.)	FMIN=0.
	IF(FMN.LT.0.)	FMIN=FMN*1.3
	IF(FMX.GT.0.)	FMAX=FMX*1.3
	IF(FMX.LT.0.)	FMAX=0.
      ELSEIF(SCALE.EQ.'LOG') THEN 
        IF(FMN.LE.0.) THEN
	     	SCALE='LIN'
		GO TO 50
	ENDIF
	FMIN=LOG10(FMN)
	FMAX=LOG10(FMX)
       	FEXP=AINT(FMAX+(FMAX-FMIN)*0.2+1)
	FMIN=10.**AINT(LOG10(FMN))
	IF(FMIN.GT.FMX) FMIN=FMIN/10.
	FMAX=10.**(FEXP)                 
      ENDIF                         
      WRITE(99,100) TITS,LABS,TICS,XL,XU,YL,YU
100   FORMAT(2X,'SET FONT DUPLEX',/,                           
     *       2X,'SET TITLE SIZE ',F8.4,/,
     *       2X,'SET LABEL SIZE ',F8.4,/,
     *       2X,'SET TICKS TOP OFF SIZE ',F8.4,/,
     *       2X,'SET WINDOW X ',F8.4,' TO ',F8.4,/,
     *       2X,'SET WINDOW Y ',F8.4,' TO ',F8.4)
      XTIT=XL+XTIT0
      YTIT=YU+YTIT0
      WRITE(99,101) XL,YTIT,TITLE(NH)(1:25)
101   FORMAT('  TITLE ',2(F8.4,1X),'"',A,'"')                  
      YTIT=YTIT-2.*YTIT0
      WRITE(99,102) XTIT,YTIT,HINT(NH)
102   FORMAT('  TITLE ',2(F8.4,1X),'" INT=',1PE10.3,'"')                  
      YTIT=YTIT-YTIT0
      WRITE(99,103) XTIT,YTIT,IENT(NH)
103   FORMAT('  TITLE ',2(F8.4,1X),'" ENT=',I9,'"')                  
      YTIT=YTIT-YTIT0                         
      IF(USCORE(NH).NE.0.) THEN
        WRITE(99,104) XTIT,YTIT,USCORE(NH)
104     FORMAT('  TITLE ',2(F8.4,1X),'" UFL=',1PE10.3,'"')                  
        YTIT=YTIT-YTIT0                      
      ENDIF
      IF(OSCORE(NH).NE.0.) THEN
        WRITE(99,105) XTIT,YTIT,OSCORE(NH)
105     FORMAT('  TITLE ',2(F8.4,1X),'" OFL=',1PE10.3,'"')                  
        YTIT=YTIT-YTIT0                      
      ENDIF      
      WRITE(99,106) XTIT,YTIT,XU,YTIT,XTIT,YTIT,XTIT,YU
106   FORMAT(2X,'SET ORD X Y ',/,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT',/,
     *       2X,2(F8.4,1X),/,2(F8.4,1X),/,
     *       2X,'JOIN TEXT')                                    
      WRITE(99,108) TITS*1.5
108   FORMAT(2X,'SET TITLE SIZE ',F8.4)
      WRITE(99,107) BTIT,XL-0.75*XD*XTFRAC,YL+(YU-YL)/3.,LTIT,SCALE,
     * HMIN(NH),HMAX(NH),FMIN,FMAX
107   FORMAT(                                           
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE ',f10.5,f10.5,' ANGLE 90 ','"',A,'"',/1X,
     &' SET SCALE Y ',A,/1X,
     &' SET TICKS TOP OFF   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET LIMITS Y ',1PE10.3,' ',1PE10.3,/1X,
     &' SET ORDER X Y DY')
C                       
C  END HEADER , FILL TOPDRAWER WITH DATA
C
      ENTRY MTFILL(NH,NE,N,M,BTIT,LTIT,SCA)
      IP=IP+1                             
      IF(IP.GT.4) IP=1
      WRITE(99,110) TITLE(NH),HINT(NH),IENT(NH),USCORE(NH),OSCORE(NH),NH
 110  FORMAT(' ( ',A,/,' ( INT=',1PE10.3,'  ENTRIES=',I12,' UFLOW='
     $     ,1PE10.3,' OFLOW=',1PE10.3,/1X,' ( BEGIN HIST=',I3)
      DO 200 IBIN=1,NBIN(NH)           
      IF(HIST(NH,IBIN).EQ.0..AND.HIST(NE,IBIN).EQ.0.) GO TO 200
      WRITE(99,'(3X,F10.4,2(2X,E15.4))')  
     &          XHIS(NH,IBIN),HIST(NH,IBIN),HIST(NE,IBIN)
200   CONTINUE
      WRITE(99,201)  NH,PLOT(IP)
      IF(BOOK(NE).NE.'NO')   WRITE(99,*)  '  PLOT'
201   FORMAT(' ( END HIST=',I3,/2X,'HIST ',A)           
202   FORMAT('   NEW PLOT',/,/)
203   RETURN     
      END                     

      SUBROUTINE NEWPLOT
      WRITE(99,202) 
202   FORMAT('   NEW PLOT',/,/)
      CALL MULTITOP(1,1,1,1,' ',' ','REF')
      END                         

C*******************************************************************
C     END OF THE HISTOGRAMMING PACKAGE
C*******************************************************************

C***** GENERIC UTILITIES *********************************
C
C   VSUM(P,Q,K)    returns K(I)=P(I)+Q(I) , I=1,4
C   DOT(P,Q)       Lorentz scalar product (+---)
C   RNDINT(N,INDX) returns a random ordering of 1,2,..,N
C   UTSORT(A,N,K,IOPT):
C     Sort A(N) into ascending order
C     IOPT = 1 : return sorted A and index array K
C     IOPT = 2 : return index array K only

      SUBROUTINE VSUM(P1,P2,Q)
      REAL P1(4),P2(4),Q(4)
      DO I=1,4
      Q(I)=P1(I)+P2(I)
      END DO
      END

      SUBROUTINE RNDINT(N,INDX)
C   returns a random ordering of the string of integers (1,2,..,N)
      real*8 rn,rangen2
      DIMENSION X(500),INDX(N)
      REAL*8 X
      INTEGER INDX
      IF(N.GT.500) then
        WRITE(6,*) 'WARNING**RNDINT - more than 500 entries'
        stop
      endif
      DO I=1,N
      X(I)=rangen2(1)
      ENDDO
      CALL alusor(X,N,INDX,2)
      END

*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALUSOR(A,N,K,IOPT)
C-----------------------------------------------------------------------
C     Sort A(N) into ascending order
C     IOPT = 1 : return sorted A and index array K
C     IOPT = 2 : return index array K only
C-----------------------------------------------------------------------
      DOUBLE PRECISION A(N),B(500)
      INTEGER N,I,J,IOPT,K(N),IL(500),IR(500)
      IF (N.GT.500) then
         write(*,*) 'array too large to sort in ALUSOR'
         stop
      endif 
      IL(1)=0
      IR(1)=0
      DO 10 I=2,N
      IL(I)=0
      IR(I)=0
      J=1
   2  IF(A(I).GT.A(J)) GOTO 5
   3  IF(IL(J).EQ.0) GOTO 4
      J=IL(J)
      GOTO 2
   4  IR(I)=-J
      IL(J)=I
      GOTO 10
   5  IF(IR(J).LE.0) GOTO 6
      J=IR(J)
      GOTO 2
   6  IR(I)=IR(J)
      IR(J)=I
  10  CONTINUE
      I=1
      J=1
      GOTO 8
  20  J=IL(J)
   8  IF(IL(J).GT.0) GOTO 20
   9  K(I)=J
      B(I)=A(J)
      I=I+1
      IF(IR(J)) 12,30,13
  13  J=IR(J)
      GOTO 8
  12  J=-IR(J)
      GOTO 9
  30  IF(IOPT.EQ.2) RETURN
      DO 31 I=1,N
  31  A(I)=B(I)
 999  END


      SUBROUTINE MBWARN(ROUT)
      CHARACTER*(*) ROUT
      WRITE(*,*) '***********************************************'
      WRITE(*,*) '***** WARNING CALLED FROM ROUTINE ',ROUT,':'
      END
