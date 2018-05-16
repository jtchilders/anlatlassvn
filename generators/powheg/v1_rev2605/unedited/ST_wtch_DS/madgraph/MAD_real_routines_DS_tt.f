      subroutine choose_real_process_tt(p,flav,amp2)
      real * 8 p(0:3,1: 5)
      integer flav( 5)
      real * 8 amp2

      if ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttbbx_twmbx(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttbbx_twmdx(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttbbx_twmsx(p,amp2)

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttbxb_twmbx(p,amp2)

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttbxb_twmdx(p,amp2)

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttbxb_twmsx(p,amp2)

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttgg_twmbx(p,amp2)

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttgg_twmdx(p,amp2)

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttgg_twmsx(p,amp2)

      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuux_twmbx(p,amp2)

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuux_twmbx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuux_twmbx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuux_twmbx(p,amp2)







      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuux_twmdx(p,amp2)

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuux_twmdx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuux_twmdx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuux_twmdx(p,amp2)







      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuux_twmsx(p,amp2)

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuux_twmsx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuux_twmsx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuux_twmsx(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sttuxu_twmbx(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuxu_twmdx(p,amp2)

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuxu_twmdx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuxu_twmdx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sttuxu_twmdx(p,amp2)










      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuxu_twmsx(p,amp2)

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuxu_twmsx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuxu_twmsx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sttuxu_twmsx(p,amp2)




      else
         amp2=0d0  !not a double resonant subprocess



      endif
      end

      SUBROUTINE Sttbbx_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- b~  
C  
C Crossing   1 is b b~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbbx_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmbx(maxamps), jamp2ttbbx_twmbx(0:maxamps)
      common/to_ampsttbbx_twmbx/  amp2ttbbx_twmbx,       jamp2ttbbx_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbbx_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbbx_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbbx_twmbx(ihel)=0d0
              jamp2ttbbx_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbbx_twmbx(0))
              jamp2ttbbx_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbbx_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbbx_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbbx_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmbx(maxamps), jamp2ttbbx_twmbx(0:maxamps)
      common/to_ampsttbbx_twmbx/  amp2ttbbx_twmbx,       jamp2ttbbx_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttbbx_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbbx_twmbx =ttbbx_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbbx_twmbx(i)=amp2ttbbx_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbbx_twmbx(i)=Jamp2ttbbx_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttbbx_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- d~  
C  
C Crossing   1 is b b~ -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbbx_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmdx(maxamps), jamp2ttbbx_twmdx(0:maxamps)
      common/to_ampsttbbx_twmdx/  amp2ttbbx_twmdx,       jamp2ttbbx_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbbx_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbbx_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbbx_twmdx(ihel)=0d0
              jamp2ttbbx_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbbx_twmdx(0))
              jamp2ttbbx_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbbx_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbbx_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbbx_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbbx_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbbx_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmdx(maxamps), jamp2ttbbx_twmdx(0:maxamps)
      common/to_ampsttbbx_twmdx/  amp2ttbbx_twmdx,       jamp2ttbbx_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttbbx_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbbx_twmdx =ttbbx_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbbx_twmdx(i)=amp2ttbbx_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbbx_twmdx(i)=Jamp2ttbbx_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttbbx_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- s~  
C  
C Crossing   1 is b b~ -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbbx_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmsx(maxamps), jamp2ttbbx_twmsx(0:maxamps)
      common/to_ampsttbbx_twmsx/  amp2ttbbx_twmsx,       jamp2ttbbx_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbbx_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbbx_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbbx_twmsx(ihel)=0d0
              jamp2ttbbx_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbbx_twmsx(0))
              jamp2ttbbx_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbbx_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbbx_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbbx_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbbx_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbbx_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbbx_twmsx(maxamps), jamp2ttbbx_twmsx(0:maxamps)
      common/to_ampsttbbx_twmsx/  amp2ttbbx_twmsx,       jamp2ttbbx_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttbbx_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbbx_twmsx =ttbbx_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbbx_twmsx(i)=amp2ttbbx_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbbx_twmsx(i)=Jamp2ttbbx_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttbxb_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- b~  
C  
C Crossing   1 is b~ b -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbxb_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmbx(maxamps), jamp2ttbxb_twmbx(0:maxamps)
      common/to_ampsttbxb_twmbx/  amp2ttbxb_twmbx,       jamp2ttbxb_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbxb_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbxb_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbxb_twmbx(ihel)=0d0
              jamp2ttbxb_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbxb_twmbx(0))
              jamp2ttbxb_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbxb_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbxb_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbxb_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbxb_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbxb_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmbx(maxamps), jamp2ttbxb_twmbx(0:maxamps)
      common/to_ampsttbxb_twmbx/  amp2ttbxb_twmbx,       jamp2ttbxb_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttbxb_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbxb_twmbx =ttbxb_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbxb_twmbx(i)=amp2ttbxb_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbxb_twmbx(i)=Jamp2ttbxb_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttbxb_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- d~  
C  
C Crossing   1 is b~ b -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbxb_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmdx(maxamps), jamp2ttbxb_twmdx(0:maxamps)
      common/to_ampsttbxb_twmdx/  amp2ttbxb_twmdx,       jamp2ttbxb_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbxb_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbxb_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbxb_twmdx(ihel)=0d0
              jamp2ttbxb_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbxb_twmdx(0))
              jamp2ttbxb_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbxb_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbxb_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbxb_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbxb_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbxb_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmdx(maxamps), jamp2ttbxb_twmdx(0:maxamps)
      common/to_ampsttbxb_twmdx/  amp2ttbxb_twmdx,       jamp2ttbxb_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttbxb_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbxb_twmdx =ttbxb_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbxb_twmdx(i)=amp2ttbxb_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbxb_twmdx(i)=Jamp2ttbxb_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttbxb_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- s~  
C  
C Crossing   1 is b~ b -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttbxb_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmsx(maxamps), jamp2ttbxb_twmsx(0:maxamps)
      common/to_ampsttbxb_twmsx/  amp2ttbxb_twmsx,       jamp2ttbxb_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttbxb_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttbxb_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttbxb_twmsx(ihel)=0d0
              jamp2ttbxb_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttbxb_twmsx(0))
              jamp2ttbxb_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttbxb_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttbxb_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttbxb_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttbxb_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttbxb_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttbxb_twmsx(maxamps), jamp2ttbxb_twmsx(0:maxamps)
      common/to_ampsttbxb_twmsx/  amp2ttbxb_twmsx,       jamp2ttbxb_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttbxb_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttbxb_twmsx =ttbxb_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttbxb_twmsx(i)=amp2ttbxb_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttbxb_twmsx(i)=Jamp2ttbxb_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttgg_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- b~  
C  
C Crossing   1 is g g -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttgg_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmbx(maxamps), jamp2ttgg_twmbx(0:maxamps)
      common/to_ampsttgg_twmbx/  amp2ttgg_twmbx,       jamp2ttgg_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttgg_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    3/          
      DATA jamp2ttgg_twmbx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttgg_twmbx(ihel)=0d0
              jamp2ttgg_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttgg_twmbx(0))
              jamp2ttgg_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttgg_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttgg_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttgg_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttgg_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttgg_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   3,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  11, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmbx(maxamps), jamp2ttgg_twmbx(0:maxamps)
      common/to_ampsttgg_twmbx/  amp2ttgg_twmbx,       jamp2ttgg_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,4   ),GWFTB ,AMP(2   ))          
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,10  ))      
      CALL FVOXXX(W(1,3   ),W(1,10  ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTB ,AMP(3   ))          
      JAMP(   1) = -AMP(   1)+AMP(   3)
      JAMP(   2) = -AMP(   2)-AMP(   3)
      ttgg_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttgg_twmbx =ttgg_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttgg_twmbx(i)=amp2ttgg_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttgg_twmbx(i)=Jamp2ttgg_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttgg_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- d~  
C  
C Crossing   1 is g g -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttgg_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmdx(maxamps), jamp2ttgg_twmdx(0:maxamps)
      common/to_ampsttgg_twmdx/  amp2ttgg_twmdx,       jamp2ttgg_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttgg_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    3/          
      DATA jamp2ttgg_twmdx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttgg_twmdx(ihel)=0d0
              jamp2ttgg_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttgg_twmdx(0))
              jamp2ttgg_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttgg_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttgg_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttgg_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttgg_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttgg_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   3,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  11, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmdx(maxamps), jamp2ttgg_twmdx(0:maxamps)
      common/to_ampsttgg_twmdx/  amp2ttgg_twmdx,       jamp2ttgg_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,4   ),GWFTD ,AMP(2   ))          
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,10  ))      
      CALL FVOXXX(W(1,3   ),W(1,10  ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTD ,AMP(3   ))          
      JAMP(   1) = -AMP(   1)+AMP(   3)
      JAMP(   2) = -AMP(   2)-AMP(   3)
      ttgg_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttgg_twmdx =ttgg_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttgg_twmdx(i)=amp2ttgg_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttgg_twmdx(i)=Jamp2ttgg_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttgg_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- s~  
C  
C Crossing   1 is g g -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttgg_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmsx(maxamps), jamp2ttgg_twmsx(0:maxamps)
      common/to_ampsttgg_twmsx/  amp2ttgg_twmsx,       jamp2ttgg_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttgg_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    3/          
      DATA jamp2ttgg_twmsx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttgg_twmsx(ihel)=0d0
              jamp2ttgg_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttgg_twmsx(0))
              jamp2ttgg_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttgg_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttgg_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttgg_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttgg_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttgg_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   3,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  11, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttgg_twmsx(maxamps), jamp2ttgg_twmsx(0:maxamps)
      common/to_ampsttgg_twmsx/  amp2ttgg_twmsx,       jamp2ttgg_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL IOVXXX(W(1,5   ),W(1,9   ),W(1,4   ),GWFTS ,AMP(2   ))          
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,10  ))      
      CALL FVOXXX(W(1,3   ),W(1,10  ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTS ,AMP(3   ))          
      JAMP(   1) = -AMP(   1)+AMP(   3)
      JAMP(   2) = -AMP(   2)-AMP(   3)
      ttgg_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttgg_twmsx =ttgg_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttgg_twmsx(i)=amp2ttgg_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttgg_twmsx(i)=Jamp2ttgg_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuux_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- b~  
C  
C Crossing   1 is u u~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuux_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmbx(maxamps), jamp2ttuux_twmbx(0:maxamps)
      common/to_ampsttuux_twmbx/  amp2ttuux_twmbx,       jamp2ttuux_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuux_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuux_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuux_twmbx(ihel)=0d0
              jamp2ttuux_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuux_twmbx(0))
              jamp2ttuux_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuux_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuux_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuux_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuux_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuux_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmbx(maxamps), jamp2ttuux_twmbx(0:maxamps)
      common/to_ampsttuux_twmbx/  amp2ttuux_twmbx,       jamp2ttuux_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttuux_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuux_twmbx =ttuux_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuux_twmbx(i)=amp2ttuux_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuux_twmbx(i)=Jamp2ttuux_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuux_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- d~  
C  
C Crossing   1 is u u~ -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuux_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmdx(maxamps), jamp2ttuux_twmdx(0:maxamps)
      common/to_ampsttuux_twmdx/  amp2ttuux_twmdx,       jamp2ttuux_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuux_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuux_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuux_twmdx(ihel)=0d0
              jamp2ttuux_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuux_twmdx(0))
              jamp2ttuux_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuux_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuux_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuux_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuux_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuux_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmdx(maxamps), jamp2ttuux_twmdx(0:maxamps)
      common/to_ampsttuux_twmdx/  amp2ttuux_twmdx,       jamp2ttuux_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttuux_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuux_twmdx =ttuux_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuux_twmdx(i)=amp2ttuux_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuux_twmdx(i)=Jamp2ttuux_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuux_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- s~  
C  
C Crossing   1 is u u~ -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuux_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmsx(maxamps), jamp2ttuux_twmsx(0:maxamps)
      common/to_ampsttuux_twmsx/  amp2ttuux_twmsx,       jamp2ttuux_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuux_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuux_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuux_twmsx(ihel)=0d0
              jamp2ttuux_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuux_twmsx(0))
              jamp2ttuux_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuux_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuux_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuux_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuux_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuux_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuux_twmsx(maxamps), jamp2ttuux_twmsx(0:maxamps)
      common/to_ampsttuux_twmsx/  amp2ttuux_twmsx,       jamp2ttuux_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      JAMP(   1) = -AMP(   1)
      ttuux_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuux_twmsx =ttuux_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuux_twmsx(i)=amp2ttuux_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuux_twmsx(i)=Jamp2ttuux_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuxu_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- b~  
C  
C Crossing   1 is u~ u -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuxu_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmbx(maxamps), jamp2ttuxu_twmbx(0:maxamps)
      common/to_ampsttuxu_twmbx/  amp2ttuxu_twmbx,       jamp2ttuxu_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuxu_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuxu_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuxu_twmbx(ihel)=0d0
              jamp2ttuxu_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuxu_twmbx(0))
              jamp2ttuxu_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuxu_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuxu_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuxu_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuxu_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuxu_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmbx(maxamps), jamp2ttuxu_twmbx(0:maxamps)
      common/to_ampsttuxu_twmbx/  amp2ttuxu_twmbx,       jamp2ttuxu_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttuxu_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuxu_twmbx =ttuxu_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuxu_twmbx(i)=amp2ttuxu_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuxu_twmbx(i)=Jamp2ttuxu_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuxu_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- d~  
C  
C Crossing   1 is u~ u -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuxu_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmdx(maxamps), jamp2ttuxu_twmdx(0:maxamps)
      common/to_ampsttuxu_twmdx/  amp2ttuxu_twmdx,       jamp2ttuxu_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuxu_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuxu_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuxu_twmdx(ihel)=0d0
              jamp2ttuxu_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuxu_twmdx(0))
              jamp2ttuxu_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuxu_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuxu_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuxu_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuxu_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuxu_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmdx(maxamps), jamp2ttuxu_twmdx(0:maxamps)
      common/to_ampsttuxu_twmdx/  amp2ttuxu_twmdx,       jamp2ttuxu_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttuxu_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuxu_twmdx =ttuxu_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuxu_twmdx(i)=amp2ttuxu_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuxu_twmdx(i)=Jamp2ttuxu_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sttuxu_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- s~  
C  
C Crossing   1 is u~ u -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ttuxu_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmsx(maxamps), jamp2ttuxu_twmsx(0:maxamps)
      common/to_ampsttuxu_twmsx/  amp2ttuxu_twmsx,       jamp2ttuxu_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixttuxu_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    1/          
      DATA jamp2ttuxu_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ttuxu_twmsx(ihel)=0d0
              jamp2ttuxu_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ttuxu_twmsx(0))
              jamp2ttuxu_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ttuxu_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ttuxu_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ttuxu_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ttuxu_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ttuxu_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   1,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   7, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ttuxu_twmsx(maxamps), jamp2ttuxu_twmsx(0:maxamps)
      common/to_ampsttuxu_twmsx/  amp2ttuxu_twmsx,       jamp2ttuxu_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      JAMP(   1) = +AMP(   1)
      ttuxu_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ttuxu_twmsx =ttuxu_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ttuxu_twmsx(i)=amp2ttuxu_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ttuxu_twmsx(i)=Jamp2ttuxu_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
