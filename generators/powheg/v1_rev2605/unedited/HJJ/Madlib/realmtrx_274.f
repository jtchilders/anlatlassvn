      SUBROUTINE SREALMTRX_274(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> h s s~ g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "nexternal.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
      INTEGER NGRAPHS
      PARAMETER (NGRAPHS=  84)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 REALMTRX_274
      REAL*8 ZERO
      PARAMETER(ZERO=0d0)
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I,L,K
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      INTEGER NGOOD,igood(ncomb),jhel
      data ngood /0/
      save igood,jhel
      REAL*8 hwgt
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_Ramps_274/  amp2,       jamp2

      integer j,jj
      integer max_bhel
      parameter ( max_bhel =          32 )

      INTEGER NCOLOR
      DATA NCOLOR   /   6/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 6) / 1, 1,-1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 6) / 1, 2, 3, 4, 5, 6/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      DO IHEL=1,NGRAPHS
          amp2(ihel)=0d0
      ENDDO
      jamp2(0)=dble(NCOLOR)
      DO IHEL=1,int(jamp2(0))
          jamp2(ihel)=0d0
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=REALMTRX_274(P1,NHEL(1,IHEL),IHEL,JC(1))              
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION REALMTRX_274(P,NHEL,HELL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> h s s~ g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  84,NEIGEN=  6) 
      include "nexternal.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  97, NCOLOR=   6) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL), HELL
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
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_Ramps_274/  amp2,       jamp2
      integer max_bhel
      parameter ( max_bhel =          32 )
      include "coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            9/                                       
      DATA (CF(i,1  ),i=1  ,6  ) /    64,   -8,    1,    1,   -8,   10/    
C               T[ 4, 5, 2, 1, 6]                                          
      DATA Denom(2  )/            9/                                       
      DATA (CF(i,2  ),i=1  ,6  ) /    -8,   64,   -8,   10,    1,    1/    
C               T[ 4, 5, 2, 6, 1]                                          
      DATA Denom(3  )/            9/                                       
      DATA (CF(i,3  ),i=1  ,6  ) /     1,   -8,   64,    1,   10,   -8/    
C               T[ 4, 5, 6, 2, 1]                                          
      DATA Denom(4  )/            9/                                       
      DATA (CF(i,4  ),i=1  ,6  ) /     1,   10,    1,   64,   -8,   -8/    
C               T[ 4, 5, 1, 6, 2]                                          
      DATA Denom(5  )/            9/                                       
      DATA (CF(i,5  ),i=1  ,6  ) /    -8,    1,   10,   -8,   64,    1/    
C               T[ 4, 5, 1, 2, 6]                                          
      DATA Denom(6  )/            9/                                       
      DATA (CF(i,6  ),i=1  ,6  ) /    10,    1,   -8,   -8,    1,   64/    
C               T[ 4, 5, 6, 1, 2]                                          
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL SXXXXX(P(0,3   ),+1*IC(3   ),W(1,3   ))                         
      CALL OXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL VXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL JVSHXX(W(1,2   ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL FVOXXX(W(1,4   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL JIOXXX(W(1,5   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,9   ),G ,AMP(1   ))              
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL JIOXXX(W(1,10  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,11  ))     
      CALL JVSHXX(W(1,11  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,12  ),G ,AMP(2   ))              
      CALL JVVXXX(W(1,11  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,13  ))      
      CALL VVSHXX(W(1,6   ),W(1,13  ),W(1,3   ),GH ,AMP(3   ))             
      CALL FVIXXX(W(1,10  ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,4   ),W(1,6   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,8   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,15  ))     
      CALL IOVXXX(W(1,5   ),W(1,15  ),W(1,6   ),GG ,AMP(5   ))             
      CALL VVVSXX(W(1,11  ),W(1,2   ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(6   ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,7   ),W(1,11  ),G ,AMP(7   ))              
      CALL IOVXXX(W(1,10  ),W(1,8   ),W(1,6   ),GG ,AMP(8   ))             
      CALL FVOXXX(W(1,4   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,16  ))     
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,17  ))     
      CALL JIOXXX(W(1,17  ),W(1,16  ),GG ,ZERO    ,ZERO    ,W(1,18  ))     
      CALL VVSHXX(W(1,6   ),W(1,18  ),W(1,3   ),GH ,AMP(9   ))             
      CALL FVOXXX(W(1,16  ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,19  ))     
      CALL JIOXXX(W(1,5   ),W(1,19  ),GG ,ZERO    ,ZERO    ,W(1,20  ))     
      CALL VVSHXX(W(1,6   ),W(1,20  ),W(1,3   ),GH ,AMP(10  ))             
      CALL JIOXXX(W(1,5   ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,21  ))     
      CALL JVVXXX(W(1,21  ),W(1,7   ),G ,ZERO    ,ZERO    ,W(1,22  ))      
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,22  ),G ,AMP(11  ))              
      CALL UVVAXX(W(1,21  ),W(1,7   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,23  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,1   ),W(1,23  ),G ,ZERO    ,AMP(12  ))     
      CALL JVVXXX(W(1,21  ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,24  ))      
      CALL JVSHXX(W(1,24  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,25  ))     
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,25  ),G ,AMP(13  ))              
      CALL UVVAXX(W(1,21  ),W(1,1   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,26  ))                                                      
      CALL UTSAXX(W(1,26  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,27  ))     
      CALL VVTAXX(W(1,6   ),W(1,2   ),W(1,27  ),G ,ZERO    ,AMP(14  ))     
      CALL JVVXXX(W(1,24  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,28  ))      
      CALL VVSHXX(W(1,6   ),W(1,28  ),W(1,3   ),GH ,AMP(15  ))             
      CALL JVTAXX(W(1,2   ),W(1,26  ),G ,ZERO    ,ZERO    ,W(1,29  ))      
      CALL VVSHXX(W(1,6   ),W(1,29  ),W(1,3   ),GH ,AMP(16  ))             
      CALL JVVXXX(W(1,7   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,30  ))      
      CALL FVIXXX(W(1,5   ),W(1,30  ),GG ,ZERO    ,ZERO    ,W(1,31  ))     
      CALL IOVXXX(W(1,31  ),W(1,4   ),W(1,6   ),GG ,AMP(17  ))             
      CALL FVOXXX(W(1,4   ),W(1,30  ),GG ,ZERO    ,ZERO    ,W(1,32  ))     
      CALL IOVXXX(W(1,5   ),W(1,32  ),W(1,6   ),GG ,AMP(18  ))             
      CALL VVVSXX(W(1,24  ),W(1,2   ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(19  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,7   ),W(1,24  ),G ,AMP(20  ))              
      CALL VVTAXX(W(1,6   ),W(1,7   ),W(1,26  ),G ,ZERO    ,AMP(21  ))     
      CALL VVVXXX(W(1,6   ),W(1,30  ),W(1,21  ),G ,AMP(22  ))              
      CALL UVVAXX(W(1,7   ),W(1,1   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,33  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,21  ),W(1,33  ),G ,ZERO    ,AMP(23  ))     
      CALL FVIXXX(W(1,5   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,34  ))     
      CALL JIOXXX(W(1,34  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,35  ))     
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,35  ),G ,AMP(24  ))              
      CALL JIOXXX(W(1,5   ),W(1,16  ),GG ,ZERO    ,ZERO    ,W(1,36  ))     
      CALL JVSHXX(W(1,36  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,37  ))     
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,37  ),G ,AMP(25  ))              
      CALL JVVXXX(W(1,36  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,38  ))      
      CALL VVSHXX(W(1,6   ),W(1,38  ),W(1,3   ),GH ,AMP(26  ))             
      CALL FVIXXX(W(1,34  ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,39  ))     
      CALL IOVXXX(W(1,39  ),W(1,4   ),W(1,6   ),GG ,AMP(27  ))             
      CALL FVOXXX(W(1,16  ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,40  ))     
      CALL IOVXXX(W(1,5   ),W(1,40  ),W(1,6   ),GG ,AMP(28  ))             
      CALL VVVSXX(W(1,36  ),W(1,2   ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(29  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,7   ),W(1,36  ),G ,AMP(30  ))              
      CALL IOVXXX(W(1,34  ),W(1,16  ),W(1,6   ),GG ,AMP(31  ))             
      CALL FVOXXX(W(1,4   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,41  ))     
      CALL JIOXXX(W(1,10  ),W(1,41  ),GG ,ZERO    ,ZERO    ,W(1,42  ))     
      CALL VVSHXX(W(1,6   ),W(1,42  ),W(1,3   ),GH ,AMP(32  ))             
      CALL JIOXXX(W(1,17  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,43  ))     
      CALL JVSHXX(W(1,43  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,44  ))     
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,44  ),G ,AMP(33  ))              
      CALL JVSHXX(W(1,1   ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,45  ))     
      CALL FVOXXX(W(1,4   ),W(1,45  ),GG ,ZERO    ,ZERO    ,W(1,46  ))     
      CALL JIOXXX(W(1,5   ),W(1,46  ),GG ,ZERO    ,ZERO    ,W(1,47  ))     
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,47  ),G ,AMP(34  ))              
      CALL JVVXXX(W(1,1   ),W(1,43  ),G ,ZERO    ,ZERO    ,W(1,48  ))      
      CALL VVSHXX(W(1,6   ),W(1,48  ),W(1,3   ),GH ,AMP(35  ))             
      CALL FVIXXX(W(1,17  ),W(1,45  ),GG ,ZERO    ,ZERO    ,W(1,49  ))     
      CALL IOVXXX(W(1,49  ),W(1,4   ),W(1,6   ),GG ,AMP(36  ))             
      CALL FVOXXX(W(1,46  ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,50  ))     
      CALL IOVXXX(W(1,5   ),W(1,50  ),W(1,6   ),GG ,AMP(37  ))             
      CALL VVVSXX(W(1,1   ),W(1,43  ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(38  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,45  ),W(1,43  ),G ,AMP(39  ))              
      CALL IOVXXX(W(1,17  ),W(1,46  ),W(1,6   ),GG ,AMP(40  ))             
      CALL FVOXXX(W(1,41  ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,51  ))     
      CALL JIOXXX(W(1,5   ),W(1,51  ),GG ,ZERO    ,ZERO    ,W(1,52  ))     
      CALL VVSHXX(W(1,6   ),W(1,52  ),W(1,3   ),GH ,AMP(41  ))             
      CALL JVVXXX(W(1,21  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,53  ))      
      CALL JVSHXX(W(1,53  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,54  ))     
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,54  ),G ,AMP(42  ))              
      CALL UVVAXX(W(1,21  ),W(1,2   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,55  ))                                                      
      CALL UTSAXX(W(1,55  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,56  ))     
      CALL VVTAXX(W(1,6   ),W(1,1   ),W(1,56  ),G ,ZERO    ,AMP(43  ))     
      CALL JVVXXX(W(1,21  ),W(1,45  ),G ,ZERO    ,ZERO    ,W(1,57  ))      
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,57  ),G ,AMP(44  ))              
      CALL UVVAXX(W(1,21  ),W(1,45  ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,58  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,2   ),W(1,58  ),G ,ZERO    ,AMP(45  ))     
      CALL JVVXXX(W(1,1   ),W(1,53  ),G ,ZERO    ,ZERO    ,W(1,59  ))      
      CALL VVSHXX(W(1,6   ),W(1,59  ),W(1,3   ),GH ,AMP(46  ))             
      CALL JVTAXX(W(1,1   ),W(1,55  ),G ,ZERO    ,ZERO    ,W(1,60  ))      
      CALL VVSHXX(W(1,6   ),W(1,60  ),W(1,3   ),GH ,AMP(47  ))             
      CALL JVVXXX(W(1,45  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,61  ))      
      CALL FVIXXX(W(1,5   ),W(1,61  ),GG ,ZERO    ,ZERO    ,W(1,62  ))     
      CALL IOVXXX(W(1,62  ),W(1,4   ),W(1,6   ),GG ,AMP(48  ))             
      CALL FVOXXX(W(1,4   ),W(1,61  ),GG ,ZERO    ,ZERO    ,W(1,63  ))     
      CALL IOVXXX(W(1,5   ),W(1,63  ),W(1,6   ),GG ,AMP(49  ))             
      CALL VVVSXX(W(1,1   ),W(1,53  ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(50  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,45  ),W(1,53  ),G ,AMP(51  ))              
      CALL VVTAXX(W(1,6   ),W(1,45  ),W(1,55  ),G ,ZERO    ,AMP(52  ))     
      CALL VVVXXX(W(1,6   ),W(1,61  ),W(1,21  ),G ,AMP(53  ))              
      CALL UVVAXX(W(1,45  ),W(1,2   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,64  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,21  ),W(1,64  ),G ,ZERO    ,AMP(54  ))     
      CALL JIOXXX(W(1,5   ),W(1,41  ),GG ,ZERO    ,ZERO    ,W(1,65  ))     
      CALL JVSHXX(W(1,65  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,66  ))     
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,66  ),G ,AMP(55  ))              
      CALL FVIXXX(W(1,5   ),W(1,45  ),GG ,ZERO    ,ZERO    ,W(1,67  ))     
      CALL JIOXXX(W(1,67  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,68  ))     
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,68  ),G ,AMP(56  ))              
      CALL JVVXXX(W(1,1   ),W(1,65  ),G ,ZERO    ,ZERO    ,W(1,69  ))      
      CALL VVSHXX(W(1,6   ),W(1,69  ),W(1,3   ),GH ,AMP(57  ))             
      CALL FVIXXX(W(1,67  ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,70  ))     
      CALL IOVXXX(W(1,70  ),W(1,4   ),W(1,6   ),GG ,AMP(58  ))             
      CALL FVOXXX(W(1,41  ),W(1,45  ),GG ,ZERO    ,ZERO    ,W(1,71  ))     
      CALL IOVXXX(W(1,5   ),W(1,71  ),W(1,6   ),GG ,AMP(59  ))             
      CALL VVVSXX(W(1,1   ),W(1,65  ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(60  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,45  ),W(1,65  ),G ,AMP(61  ))              
      CALL IOVXXX(W(1,67  ),W(1,41  ),W(1,6   ),GG ,AMP(62  ))             
      CALL FVIXXX(W(1,10  ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,72  ))     
      CALL JIOXXX(W(1,72  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,73  ))     
      CALL VVSHXX(W(1,6   ),W(1,73  ),W(1,3   ),GH ,AMP(63  ))             
      CALL FVIXXX(W(1,17  ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,74  ))     
      CALL JIOXXX(W(1,74  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,75  ))     
      CALL VVSHXX(W(1,6   ),W(1,75  ),W(1,3   ),GH ,AMP(64  ))             
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,76  ))      
      CALL FVOXXX(W(1,4   ),W(1,76  ),GG ,ZERO    ,ZERO    ,W(1,77  ))     
      CALL JIOXXX(W(1,5   ),W(1,77  ),GG ,ZERO    ,ZERO    ,W(1,78  ))     
      CALL VVSHXX(W(1,6   ),W(1,78  ),W(1,3   ),GH ,AMP(65  ))             
      CALL JVSHXX(W(1,21  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,79  ))     
      CALL JVVXXX(W(1,79  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,80  ))      
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,80  ),G ,AMP(66  ))              
      CALL UVVAXX(W(1,79  ),W(1,2   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,81  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,1   ),W(1,81  ),G ,ZERO    ,AMP(67  ))     
      CALL JVVXXX(W(1,1   ),W(1,79  ),G ,ZERO    ,ZERO    ,W(1,82  ))      
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,82  ),G ,AMP(68  ))              
      CALL UVVAXX(W(1,1   ),W(1,79  ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,83  ))                                                      
      CALL VVTAXX(W(1,6   ),W(1,2   ),W(1,83  ),G ,ZERO    ,AMP(69  ))     
      CALL JVVXXX(W(1,21  ),W(1,76  ),G ,ZERO    ,ZERO    ,W(1,84  ))      
      CALL VVSHXX(W(1,6   ),W(1,84  ),W(1,3   ),GH ,AMP(70  ))             
      CALL UVVAXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,ZERO    ,       
     &     W(1,85  ))                                                      
      CALL JVTAXX(W(1,21  ),W(1,85  ),G ,ZERO    ,ZERO    ,W(1,86  ))      
      CALL VVSHXX(W(1,6   ),W(1,86  ),W(1,3   ),GH ,AMP(71  ))             
      CALL JVSHXX(W(1,76  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,87  ))     
      CALL FVIXXX(W(1,5   ),W(1,87  ),GG ,ZERO    ,ZERO    ,W(1,88  ))     
      CALL IOVXXX(W(1,88  ),W(1,4   ),W(1,6   ),GG ,AMP(72  ))             
      CALL FVOXXX(W(1,4   ),W(1,87  ),GG ,ZERO    ,ZERO    ,W(1,89  ))     
      CALL IOVXXX(W(1,5   ),W(1,89  ),W(1,6   ),GG ,AMP(73  ))             
      CALL VVVSXX(W(1,21  ),W(1,76  ),W(1,6   ),W(1,3   ),DUM1 ,GH4 ,      
     &     AMP(74  ))                                                      
      CALL VVVXXX(W(1,6   ),W(1,76  ),W(1,79  ),G ,AMP(75  ))              
      CALL VVTAXX(W(1,6   ),W(1,79  ),W(1,85  ),G ,ZERO    ,AMP(76  ))     
      CALL VVVXXX(W(1,6   ),W(1,87  ),W(1,21  ),G ,AMP(77  ))              
      CALL UTSAXX(W(1,85  ),W(1,3   ),GH ,ZERO    ,ZERO    ,W(1,90  ))     
      CALL VVTAXX(W(1,6   ),W(1,21  ),W(1,90  ),G ,ZERO    ,AMP(78  ))     
      CALL FVIXXX(W(1,5   ),W(1,76  ),GG ,ZERO    ,ZERO    ,W(1,91  ))     
      CALL JIOXXX(W(1,91  ),W(1,4   ),GG ,ZERO    ,ZERO    ,W(1,92  ))     
      CALL VVSHXX(W(1,6   ),W(1,92  ),W(1,3   ),GH ,AMP(79  ))             
      CALL JVVSXX(W(1,21  ),W(1,2   ),W(1,3   ),DUM1 ,GH4 ,ZERO    ,       
     &     ZERO    ,W(1,93  ))                                             
      CALL VVVXXX(W(1,6   ),W(1,1   ),W(1,93  ),G ,AMP(80  ))              
      CALL JVVSXX(W(1,1   ),W(1,21  ),W(1,3   ),DUM1 ,GH4 ,ZERO    ,       
     &     ZERO    ,W(1,94  ))                                             
      CALL VVVXXX(W(1,6   ),W(1,2   ),W(1,94  ),G ,AMP(81  ))              
      CALL JVVSXX(W(1,1   ),W(1,2   ),W(1,3   ),DUM1 ,GH4 ,ZERO    ,       
     &     ZERO    ,W(1,95  ))                                             
      CALL FVIXXX(W(1,5   ),W(1,95  ),GG ,ZERO    ,ZERO    ,W(1,96  ))     
      CALL IOVXXX(W(1,96  ),W(1,4   ),W(1,6   ),GG ,AMP(82  ))             
      CALL FVOXXX(W(1,4   ),W(1,95  ),GG ,ZERO    ,ZERO    ,W(1,97  ))     
      CALL IOVXXX(W(1,5   ),W(1,97  ),W(1,6   ),GG ,AMP(83  ))             
      CALL VVVXXX(W(1,6   ),W(1,95  ),W(1,21  ),G ,AMP(84  ))              
      JAMP(   1) = +AMP(   1)-AMP(   5)+AMP(  11)+AMP(  12)-AMP(  18)
     &             +AMP(  22)-AMP(  23)-AMP(  41)+AMP(  42)+AMP(  43)
     &             +AMP(  46)+AMP(  47)+AMP(  49)+AMP(  50)+AMP(  51)
     &             +AMP(  52)-AMP(  53)+AMP(  54)+AMP(  55)+AMP(  57)
     &             -AMP(  59)+AMP(  60)+AMP(  61)+AMP(  65)+AMP(  66)
     &             +AMP(  67)+AMP(  70)+AMP(  71)+AMP(  73)+AMP(  74)
     &             -AMP(  75)+AMP(  76)-AMP(  77)+AMP(  78)+AMP(  80)
     &             +AMP(  83)-AMP(  84)
      JAMP(   2) = -AMP(   1)+AMP(   2)-AMP(   3)-AMP(   6)+AMP(   7)
     &             -AMP(   8)-AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)
     &             +AMP(  15)-AMP(  16)+AMP(  19)-AMP(  20)-AMP(  21)
     &             -AMP(  32)-AMP(  42)-AMP(  43)-AMP(  44)-AMP(  45)
     &             -AMP(  46)-AMP(  47)-AMP(  50)-AMP(  51)-AMP(  52)
     &             -AMP(  55)+AMP(  56)-AMP(  57)-AMP(  60)-AMP(  61)
     &             -AMP(  62)-AMP(  66)-AMP(  67)+AMP(  68)+AMP(  69)
     &             -AMP(  80)+AMP(  81)
      JAMP(   3) = -AMP(   2)+AMP(   3)-AMP(   4)+AMP(   6)-AMP(   7)
     &             +AMP(  13)+AMP(  14)-AMP(  15)+AMP(  16)-AMP(  17)
     &             -AMP(  19)+AMP(  20)+AMP(  21)-AMP(  22)+AMP(  23)
     &             +AMP(  44)+AMP(  45)+AMP(  48)+AMP(  53)-AMP(  54)
     &             -AMP(  56)-AMP(  58)-AMP(  63)-AMP(  68)-AMP(  69)
     &             -AMP(  70)-AMP(  71)+AMP(  72)-AMP(  74)+AMP(  75)
     &             -AMP(  76)+AMP(  77)-AMP(  78)+AMP(  79)-AMP(  81)
     &             +AMP(  82)+AMP(  84)
      JAMP(   4) = -AMP(   9)-AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)
     &             +AMP(  15)-AMP(  16)+AMP(  19)-AMP(  20)-AMP(  21)
     &             +AMP(  24)-AMP(  25)+AMP(  26)+AMP(  29)-AMP(  30)
     &             -AMP(  31)+AMP(  33)-AMP(  34)+AMP(  35)+AMP(  38)
     &             +AMP(  39)-AMP(  40)-AMP(  42)-AMP(  43)-AMP(  44)
     &             -AMP(  45)-AMP(  46)-AMP(  47)-AMP(  50)-AMP(  51)
     &             -AMP(  52)-AMP(  66)-AMP(  67)+AMP(  68)+AMP(  69)
     &             -AMP(  80)+AMP(  81)
      JAMP(   5) = -AMP(  10)+AMP(  13)+AMP(  14)-AMP(  15)+AMP(  16)
     &             +AMP(  18)-AMP(  19)+AMP(  20)+AMP(  21)-AMP(  22)
     &             +AMP(  23)+AMP(  25)-AMP(  26)-AMP(  28)-AMP(  29)
     &             +AMP(  30)+AMP(  34)-AMP(  37)+AMP(  44)+AMP(  45)
     &             -AMP(  49)+AMP(  53)-AMP(  54)-AMP(  65)-AMP(  68)
     &             -AMP(  69)-AMP(  70)-AMP(  71)-AMP(  73)-AMP(  74)
     &             +AMP(  75)-AMP(  76)+AMP(  77)-AMP(  78)-AMP(  81)
     &             -AMP(  83)+AMP(  84)
      JAMP(   6) = +AMP(  11)+AMP(  12)+AMP(  17)+AMP(  22)-AMP(  23)
     &             -AMP(  24)-AMP(  27)-AMP(  33)-AMP(  35)-AMP(  36)
     &             -AMP(  38)-AMP(  39)+AMP(  42)+AMP(  43)+AMP(  46)
     &             +AMP(  47)-AMP(  48)+AMP(  50)+AMP(  51)+AMP(  52)
     &             -AMP(  53)+AMP(  54)-AMP(  64)+AMP(  66)+AMP(  67)
     &             +AMP(  70)+AMP(  71)-AMP(  72)+AMP(  74)-AMP(  75)
     &             +AMP(  76)-AMP(  77)+AMP(  78)-AMP(  79)+AMP(  80)
     &             -AMP(  82)-AMP(  84)
      REALMTRX_274 = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          REALMTRX_274 =REALMTRX_274+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
      END
       
       