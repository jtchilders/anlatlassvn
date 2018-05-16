cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c compute tensor coefficients needed for diagrams with massless internal legs
c
c	in:  q1(0:3),q2(0:3)	external momenta
c	     musq		regul. sccale
c
c	out: Cij(2,4)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c notation for scalar triangles is as:
c
c			   			1
c      C[k1,k2,m1s,m2s,m3s] ~ --------------------------------------------
c				[k^2-m1s][(k+k1)^2-m2s][(k+k1+k2)^2-m3s]    
c
c	we abbreviate C0[l1,l2]   = C[l1,l2,0,0,0]
c
c
c	only the case m1s = m2s = m3s = 0 is considered here
c
c
c	for the C2j the finite terms are determined according to the 
c	prescription of [OZ] 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tens_tri(q1,q2,musq,Cij)
      
      implicit none
      
      double precision q12,q13,q23,musq
      double precision q1s,q2s,q3s
      double precision q1(0:3),q2(0:3),q3(0:3)
      
      double precision zero
      parameter (zero = 0.d0)
      
      double complex Cij(0:2,4)
      
      double precision dotrr
      external dotrr
      
      double complex C0fin,B0fin
      external C0fin,B0fin
      
      double complex B012,B112,B023,B123,B013,B113,C00
      double precision f1,f2
      double complex r1,r2,r3,r4,r5,r6
      double complex C11,C12,C24,C21,C22,C23
      
      logical debug
      parameter (debug = .false.) 
      
c      double precision small
c      parameter (small = 1d-9) ! define as 'zero'
      include 'cint_param.inc'
      
      integer mu,i,j
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	initialize:
	do i = 0,2	  
	   do j = 1,4
	      cij(i,j) = 0d0
	   enddo
	enddo      	
c
c compute q4 from momentum conservation:

	do mu = 0,3
	   q3(mu) = -q1(mu)-q2(mu)
	enddo
c
c compute invariants from momenta:

      q12 = dotrr(q1(0),q2(0))
      q1s = dotrr(q1(0),q1(0))
      q2s = dotrr(q2(0),q2(0))
      q3s = dotrr(q3(0),q3(0))
	
      if (abs(q1s).lt.small) q1s = 0d0        
      if (abs(q2s).lt.small) q2s = 0d0        
      if (abs(q3s).lt.small) q3s = 0d0        
      
       
      ! B012 = B0[q1s,mus]:
      if (q1s.eq.0d0) then
         B012 = 0d0
	 B112 = 0d0
      else 	 
         B012 = B0fin(q1s,musq)
	 B112 = -0.5d0*B012
      endif	
       
      ! B023 = B0[q2s,mus]
      if (q2s.eq.0d0) then
         B023 = 0d0
	 B123 = 0d0
      else 	 
         B023 = B0fin(q2s,musq)
	 B123 = -0.5d0*B023
      endif	 
 
       ! B013 = B0[q1s+2*q12+q2s,mus] = B0[q3s,mus]
      if (q3s.eq.0d0) then
         B013 = 0d0
	 B113 = 0d0
      else 	 
         B013 = B0fin(q3s,musq)
	 B113 = -0.5d0*B013
      endif	 
     
      !C00 = C0[q1s,q2s,q12,mus] (note difference in notation of args.) 
      C00 = C0fin(q1s,q2s,q3s,musq) 	
      
      Cij(0,1) = C00

      f1 = -q1s
      f2 = -q2s-2d0*q12

      r1 =  0.5d0*(B013-B023+f1*C00) 
      r2 =  0.5d0*(B012-B013+f2*C00)

      C11 = (-(q2s*r1) + q12*r2)/(q12**2 - q1s*q2s)   
      C12 = (q12*r1 - q1s*r2)/(q12**2 - q1s*q2s)
      
      Cij(1,1) = C11
      Cij(1,2) = C12

      r3 = 0.5d0*(f1*C11+B113+B023)
      r4 = 0.5d0*(f2*C11+B112-B113) 
      r5 = 0.5d0*(f1*C12+B113-B123) 
      r6 = 0.5d0*(f2*C12-B113)

      C24 = -0.25d0*(f1*C11+f2*C12-B023-1d0) !finite terms of [OZ]

c      C24 = -0.25d0*(f1*C11+f2*C12-B023+1d0)   ! pole terms of [OZ]
c      C24 = -0.25d0*(f1*C11+f2*C12-B023)	! original [PV]

      C21 = (c24*q2s - q2s*r3 + q12*r4)/(q12**2 - q1s*q2s)
      C22 = (c24*q1s + q12*r5 - q1s*r6)/(q12**2 - q1s*q2s)
      C23 =-((c24*q12 - q12*r3 + q1s*r4)/(q12**2 - q1s*q2s))
	
      Cij(2,1) = C21
      Cij(2,2) = C22
      Cij(2,3) = C23
      Cij(2,4) = C24
      
      end
