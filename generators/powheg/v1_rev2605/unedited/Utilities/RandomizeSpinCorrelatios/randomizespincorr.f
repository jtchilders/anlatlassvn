      program randomizespincorr
      implicit none
      character * 200 filein,fileout
      character * 500 inline
      integer maxres
      parameter (maxres=100)
      integer reslist(maxres),i,j,k,in,out,ios
      include 'LesHouches.h'
      write(*,*) ' enter name of input LHE file'
      read(*,*) filein
      write(*,*) ' enter name of output LHE file'
      read(*,*) fileout
      write(*,*) ' which resonances to randomize;'
      write(*,*) ' enter sequence of integers labeling the resonance'
      write(*,*) ' position in the LHE record. Enter a negative number'
      write(*,*) ' to terminate'
      do j=1,maxres
         read(*,*) reslist(j)
         if(reslist(j).lt.0) exit
      enddo
      if(j.gt.maxres) then
         write(*,*)' increase maxres'
         call exit(-1)
      endif
      in=11
      open(unit=in,file=filein,status='old')
      out=12
      open(unit=out,file=fileout,status='unknown')
      do k=1,1000000000
 1       read(in,'(a)',iostat=ios) inline
         if(ios.eq.-1) goto 998
         if(inline.ne.'<event>') then
            write(out,'(a)')trim(inline)
         else
            write(out,'(a)') '<event>'
            read(in,*,iostat=ios)
     1           nup,idprup,xwgtup,scalup,aqedup,aqcdup
            if(ios.eq.-1) goto 998
            if(ios.ne.0) goto 1
            do i=1,nup   
               read(in,*,iostat=ios) idup(i),istup(i),mothup(1,i),
     &             mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &              vtimup(i),spinup(i)
               if(ios.eq.-1) goto 998
               if(ios.ne.0) goto 1
            enddo
            do j=1,maxres
               if(reslist(j).gt.0) then
                  call randomrotate(reslist(j))
               else
                  exit
               endif
            enddo
c write out LH record
            write(out,210)nup,idprup,xwgtup,scalup,aqedup,aqcdup
            do i=1,nup
               write(out,220) idup(i),istup(i),mothup(1,i),
     &             mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &              vtimup(i),spinup(i)
            enddo
         endif
      enddo
 998  continue
      close(in)
      close(out)
      return
 210  format(1p,2(1x,i6),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end
            

      subroutine randomrotate(ind)
      implicit none
      integer ind
      include 'LesHouches.h'
      real * 8 pres(5),vec(3),beta,r(3,3)
      logical sonof
      integer j
      beta=sqrt(pup(1,ind)**2+pup(2,ind)**2+pup(3,ind)**2)/pup(4,ind)
      vec(1)=pup(1,ind)/(beta*pup(4,ind))
      vec(2)=pup(2,ind)/(beta*pup(4,ind))
      vec(3)=pup(3,ind)/(beta*pup(4,ind))
      call uniformrot(r)
      do j=3,nup
         if(sonof(ind,j)) then
            call mboost5(1,vec,-beta,pup(:,j),pup(:,j))
            call matrixmultvec(r,pup(1:3,j))
            call mboost5(1,vec,beta,pup(:,j),pup(:,j))
         endif
      enddo
      end

      function sonof(m,k)
      implicit none
      logical sonof
      integer m,k
      include  'LesHouches.h'
      integer j,kcurr
      integer ngenerations
      parameter (ngenerations=4)
      kcurr=mothup(1,k)
      do j=1,ngenerations
         if(kcurr.eq.m) then
            sonof = .true.
            return
         endif
         kcurr = mothup(1,kcurr)
         if(kcurr.eq.0) then
            sonof = .false.
            return
         endif
      enddo
      sonof=.false.
      end


      subroutine mboost5(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(5,m),vout(5,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
         vout(5,ipart)=vin(5,ipart)
      enddo
      end


      subroutine matrixmultvec(r,v)
      implicit none
      real * 8 r(3,3),v(3),res(3)
      integer j
      do j=1,3
         res(j)=r(j,1)*v(1)+r(j,2)*v(2)+r(j,3)*v(3)
      enddo
      v = res
      end


      subroutine uniformrot(R)
      implicit none
c     Generate a uniformly distributed rotation
      real * 8 r(3,3)
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 costh,sinth,phi,gamma,sing,cosg,norm
      real * 8 random
      external random
      costh=2*random()-1
      sinth=sqrt(abs(1-costh**2))
      phi=2*pi*random()
c First axis in random direction
      r(1,1)=costh
      r(2,1)=sinth*sin(phi)
      r(3,1)=sinth*cos(phi)
c now pick a vector orthogonal to the first axis
      if(costh.gt.0.5d0) then
         norm=sqrt(r(1,1)**2+r(2,1)**2)
         r(1,2)=r(2,1)/norm
         r(2,2)=-r(1,1)/norm
         r(3,2)=0
      else
         norm=sqrt(r(2,1)**2+r(3,1)**2)
         r(1,2)=0
         r(2,2)=r(3,1)/norm
         r(3,2)=-r(2,1)/norm
      endif
c Now totate r(:,2) around r(:,1) of an arbitrary angle
      gamma = 2*pi * random()
      sing = sin(gamma)
      cosg = cos(gamma)
      call mrotate(r(:,1),sing,cosg,r(:,2))
c Last axis is cross product of 1 and 2
      r(1,3)=r(2,1)*r(3,2)-r(3,1)*r(2,2)
      r(2,3)=r(3,1)*r(1,2)-r(1,1)*r(3,2)
      r(3,3)=r(1,1)*r(2,2)-r(2,1)*r(1,2)
      end



      subroutine mrotate(dir,sinphi,cosphi,vec)
c Rotates vector vec counterclockwise around the direction
c dir (|dir|=1) with angle phi, given sin phi and cos phi.
      implicit none
      real * 8 sinphi,cosphi,dir(3),vec(3)
      real * 8 dircrossvec(3),dirdotvec
      integer i
      dircrossvec(1)=dir(2)*vec(3)-dir(3)*vec(2)
      dircrossvec(2)=dir(3)*vec(1)-dir(1)*vec(3)
      dircrossvec(3)=dir(1)*vec(2)-dir(2)*vec(1)
      dirdotvec=dir(1)*vec(1)+dir(2)*vec(2)+dir(3)*vec(3)
      do i=1,3
         vec(i)=vec(i)+sinphi*dircrossvec(i)
     #        -(1-cosphi)*(vec(i)-dir(i)*dirdotvec)
      enddo
      end

      function random()
      real * 8 random,r(1)
      call rm48(r,1)
      random=r(1)
      end


      SUBROUTINE RM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
C   Input and output entry points: RM48IN, RM48UT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++ 
C!!!                    output by RM48UT)                            ++ 
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++  
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C for 32-bit machines, use IMPLICIT DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL,TWOM49, ONE, ZERO
      save /R48ST1/
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,I10,2X,2I10)') ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.
      HALF = 0.5
      ZERO = 0.
      DO 2 II= 1, 97
      S = 0.
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM49 = T
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
      WRITE(6,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
c Bug! fixed by P. Nason, 13/12/2012
c          DO 40 IDUM = 1, NTOT
          DO 40 IDUM = 1, NOW
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
C             Replace exact zeros by 2**-49
         IF (UNI .EQ. ZERO)  THEN
            RVEC(IVEC) = TWOM49
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END
