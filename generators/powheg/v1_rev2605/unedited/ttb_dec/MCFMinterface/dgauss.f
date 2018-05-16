      DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)
      DOUBLE PRECISION W(12),X(12),A,B,EPS,DELTA,CONST,AA,BB,Y,C1,C2,S8,
     1                 S16,U,F
      DATA CONST /1.0D-25/
      DATA W
     1       / 0.10122 85362 90376 25915 25313 543D0,
     2         0.22238 10344 53374 47054 43559 944D0,
     3         0.31370 66458 77887 28733 79622 020D0,
     4         0.36268 37833 78361 98296 51504 493D0,
     5         0.02715 24594 11754 09485 17805 725D0,
     6         0.06225 35239 38647 89286 28438 370D0,
     7         0.09515 85116 82492 78480 99251 076D0,
     8         0.12462 89712 55533 87205 24762 822D0,
     9         0.14959 59888 16576 73208 15017 305D0,
     A         0.16915 65193 95002 53818 93120 790D0,
     B         0.18260 34150 44923 58886 67636 680D0,
     C         0.18945 06104 55068 49628 53967 232D0 /
      DATA X
     1       / 0.96028 98564 97536 23168 35608 686D0,
     2         0.79666 64774 13626 73959 15539 365D0,
     3         0.52553 24099 16328 98581 77390 492D0,
     4         0.18343 46424 95649 80493 94761 424D0,
     5         0.98940 09349 91649 93259 61541 735D0,
     6         0.94457 50230 73232 57607 79884 155D0,
     7         0.86563 12023 87831 74388 04678 977D0,
     8         0.75540 44083 55003 03389 51011 948D0,
     9         0.61787 62444 02643 74844 66717 640D0,
     A         0.45801 67776 57227 38634 24194 430D0,
     B         0.28160 35507 79258 91323 04605 015D0,
     C         0.09501 25098 37637 44018 53193 354D0 /
      DELTA=CONST*DABS(A-B)
      DGAUSS=0.
      AA=A
    5 Y=B-AA
      IF(DABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5D0*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 1 I = 1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I = 5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(DABS(S16-S8) .GT. EPS*(1.0D0+DABS(S16))) GO TO 4
      DGAUSS=DGAUSS+S16
      AA=BB
      GO TO 5
    4 Y=0.5D0*Y
      IF(DABS(Y) .GT. DELTA) GO TO 2
      PRINT 7
c      DGAUSS=0.
      RETURN
    7 FORMAT(1X,37HDGAUSS ... TOO HIGH ACCURACY REQUIRED)
      END
