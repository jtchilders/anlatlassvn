cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c compute tensor coefficients needed for diagrams with massless internal legs
c
c	in:  q1(0:3),q2(0:3),q3(0:3)	external momenta
c	     musq			regul. sccale
c
c	out: Dij(3,13)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c notation for massless scalar boxes is as:
c
c					    1
c  D[k1,k2,k3] ~ ----------------------------------------------
c		  [k^2][(k+k1)^2][(k+k1+k2)^2][(k+k1+k2+k3)^2]  
cc
c	only the case of massless internal propagators is considered here
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tens_box(q1,q2,q3,musq,Dij)
      
      implicit none
      
      double precision q12,q13,q23,musq
      double precision q1s,q2s,q3s,q4s
      double precision q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      double precision q1pq2(0:3),q2pq3(0:3)
      
      double precision zero
      parameter (zero = 0.d0)
      
      double complex Cij(0:2,4),Dij(0:3,13)
      
      double precision dotrr
      external dotrr
      
      double complex C0fin,B0fin,D0fin
      external C0fin,B0fin,D0fin
      
      double precision f1,f2,f3
      double complex r20,r21,r22
      double complex C0123,C0124,C0134,C0234,D00
      double complex D11,D12,D13
      
      double complex Cij_123(0:2,4),Cij_124(0:2,4),
     &               Cij_134(0:2,4),Cij_234(0:2,4)    
      double complex r30,r31,r32,r33,r34,r35,r36,r37,r38
      double complex D21,D22,D23,D24,D25,D26,D27

      double complex r40,r41,r42,r43,r44,r45,r46,r47,r48,r49,
     &      	     r50,r51,r52,r53,r54
      double complex D31,D32,D33,D34,D35,D36,D37,D38,D39,
     &      	     D310,D311,D312,D313
      
c      double precision small
c      parameter (small = 1d-6) ! define as 'zero'

      include 'dint_param.inc'
         
      logical debug
      parameter (debug = .false.) 
      
      integer mu,i,j
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	initialize:
	do i = 0,3	  
	   do j = 1,13
	      dij(i,j) = 0d0
	   enddo
	enddo      	
	do i = 0,2	  
	   do j = 1,4
	      cij_123(i,j) = 0d0
	      cij_124(i,j) = 0d0
	      cij_134(i,j) = 0d0
	      cij_234(i,j) = 0d0
	   enddo
	enddo      	
c
c compute q4 from momentum conservation:

	do mu = 0,3
	   q4(mu) = -q1(mu)-q2(mu)-q3(mu)
	enddo
c
c compute invariants from momenta:

      q12 = dotrr(q1(0),q2(0))
      q13 = dotrr(q1(0),q3(0))
      q23 = dotrr(q2(0),q3(0))
      
      q1s = dotrr(q1(0),q1(0))
      q2s = dotrr(q2(0),q2(0))
      q3s = dotrr(q3(0),q3(0))
      q4s = dotrr(q4(0),q4(0))

c for massless quarks in q(q1)->VVq(q4):
      q1s = 0d0
      q4s = 0d0

      if (abs(q1s).lt.small) q1s = 0d0        
      if (abs(q2s).lt.small) q2s = 0d0        
      if (abs(q3s).lt.small) q3s = 0d0        
      if (abs(q4s).lt.small) q4s = 0d0        

      do mu = 0,3
         q1pq2(mu) = q1(mu)+q2(mu)
         q2pq3(mu) = q2(mu)+q3(mu)
      enddo

      call tens_tri(q1,q2,musq,Cij_123)
      call tens_tri(q1,q2pq3,musq,Cij_124)
      call tens_tri(q1pq2,q3,musq,Cij_134)
      call tens_tri(q2,q3,musq,Cij_234)      

	!C0123 = C0[q1s,q2s,q12,mus];
c      C0123 = C0fin(q1s,q2s,2d0*q12+q1s+q2s,musq)    
      C0123 = Cij_123(0,1)
      
	!C0124 = C0[q1s,q2s+q3s+2*q23,q12+q13,mus]
c      C0124 = C0fin(q1s,q2s+q3s+2d0*q23,q4s,musq)
      C0124 = Cij_124(0,1)

	!C0134 = C0[q1s+q2s+2*q12,q3s,q13+q23,mus]
c      C0134 = C0fin(q1s+q2s+2d0*q12,q3s,q4s,musq)
      C0134 = Cij_134(0,1)

	!C0234 = C0[q2s,q3s,q23,mus]
c      C0234 = C0fin(q2s,q3s,2d0*q23+q2s+q3s,musq)
      C0234 = Cij_234(0,1)
      
	!D00 = D0[q1s, q2s, q3s, q4s, q12, q23, mus]
      D00 = D0fin(q1s+q2s+2d0*q12,q2s+q3s+2d0*q23,q1s,q2s,q3s,q4s,musq)	
      
      Dij(0,1) = D00
      
      f1 = -q1s
      f2 = -q2s-2.d0*q12
      f3 = -q4s+q1s+q2s+2d0*q12
 
      r20 = 0.5d0*(f1*D00+C0134-C0234)
      r21 = 0.5d0*(f2*D00+C0124-C0134)
      r22 = 0.5d0*(f3*D00+C0123-C0124)

      D11 = (q23**2*r20 - q2s*q3s*r20 + q12*q3s*r21 + q13*q2s*r22 -
     &       q23*(q13*r21 + q12*r22))/
     &(-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D12 = (q12*q3s*r20 + q13**2*r21 - q1s*q3s*r21 + q1s*q23*r22 -
     &       q13*(q23*r20 + q12*r22))/
     &(-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D13 = (q13*q2s*r20 + q1s*q23*r21 - q12*(q23*r20 + q13*r21) + 
     &       q12**2*r22 - q1s*q2s*r22)/
     &(-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))
               
      Dij(1,1) = D11
      Dij(1,2) = D12
      Dij(1,3) = D13
      
ccccccccccc
    
      r30 = 0.5d0*(f1*D11+Cij_134(1,1)+C0234)
      r31 = 0.5d0*(f2*D11+Cij_124(1,1)-Cij_134(1,1))
      r32 = 0.5d0*(f3*D11+Cij_123(1,1)-Cij_124(1,1))
      r33 = 0.5d0*(f1*D12+Cij_134(1,1)-Cij_234(1,1))
      r34 = 0.5d0*(f2*D12+Cij_124(1,2)-Cij_134(1,1))
      r35 = 0.5d0*(f3*D12+Cij_123(1,2)-Cij_124(1,2))
      r36 = 0.5d0*(f1*D13+Cij_134(1,2)-Cij_234(1,2))
      r37 = 0.5d0*(f2*D13+Cij_124(1,2)-Cij_134(1,2))
      r38 = 0.5d0*(f3*D13-Cij_124(1,2))
      
      D27 = -0.5d0*(f1*D11+f2*D12+f3*D13-C0234)
      
      D21=(d27*(-q23**2 + q2s*q3s) + q23**2*r30 - q2s*q3s*r30 + 
     &    q12*q3s*r31 +q13*q2s*r32 - q23*(q13*r31 + q12*r32))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D22 = (d27*(-q13**2 + q1s*q3s) + q12*q3s*r33 + q13**2*r34 - 
     &    q1s*q3s*r34 +q1s*q23*r35 - q13*(q23*r33 + q12*r35))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D23 = (d27*(-q12**2 + q1s*q2s) + q13*q2s*r36 + q1s*q23*r37 -
     &    q12*(q23*r36 + q13*r37) + q12**2*r38 - q1s*q2s*r38)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D24 = (d27*(q13*q23 - q12*q3s) + q12*q3s*r30 + q13**2*r31 - 
     &    q1s*q3s*r31 +q1s*q23*r32 - q13*(q23*r30 + q12*r32))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D25 = (d27*(q12*q23 - q13*q2s) + q13*q2s*r30 + q1s*q23*r31 -
     &    q12*(q23*r30 + q13*r31) + q12**2*r32 - q1s*q2s*r32)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D26=(d27*(q12*q13 - q1s*q23) + q13*q2s*r33 + q1s*q23*r34 -
     &    q12*(q23*r33 + q13*r34) + q12**2*r35 - q1s*q2s*r35)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      Dij(2,1) = D21
      Dij(2,2) = D22
      Dij(2,3) = D23
      Dij(2,4) = D24
      Dij(2,5) = D25
      Dij(2,6) = D26
      Dij(2,7) = D27
         
ccccccccccc

      r40 = 0.5d0*(f1*D27+Cij_134(2,4)-Cij_234(2,4))
      r41 = 0.5d0*(f2*D27+Cij_124(2,4)-Cij_134(2,4))
      r42 = 0.5d0*(f3*D27+Cij_123(2,4)-Cij_124(2,4))
      r43 = 0.5d0*(f1*D21+Cij_134(2,1)-C0234)
      r44 = 0.5d0*(f2*D21+Cij_124(2,1)-Cij_134(2,1))
      r45 = 0.5d0*(f3*D21+Cij_123(2,1)-Cij_124(2,1))
      r46 = 0.5d0*(f1*D22+Cij_134(2,1)-Cij_234(2,1))
      r47 = 0.5d0*(f2*D22+Cij_124(2,2)-Cij_134(2,1))
      r48 = 0.5d0*(f3*D22+Cij_123(2,2)-Cij_124(2,2))
      r49 = 0.5d0*(f1*D23+Cij_134(2,2)-Cij_234(2,2))
      r50 = 0.5d0*(f2*D23+Cij_124(2,2)-Cij_134(2,2))
      r51 = 0.5d0*(f3*D23-Cij_124(2,2))
      r52 = 0.5d0*(f1*D24+Cij_134(2,1)+Cij_234(1,1))
      r53 = 0.5d0*(f2*D24+Cij_124(2,3)-Cij_134(2,1))
      r54 = 0.5d0*(f3*D24+Cij_123(2,3)-Cij_124(2,3))

c compute D311,D312,D313 first, 
c because they are needed as input for other D3j:
      D311 = (q23**2*r40 - q2s*q3s*r40 + q12*q3s*r41 + q13*q2s*r42 -
     &    q23*(q13*r41 + q12*r42))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D312 = (q12*q3s*r40 + q13**2*r41 - q1s*q3s*r41 + q1s*q23*r42 -
     &    q13*(q23*r40 + q12*r42))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D313 = (q13*q2s*r40 + q1s*q23*r41 - q12*(q23*r40 + q13*r41) + 
     &    q12**2*r42 - q1s*q2s*r42)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))
     
c compute D31 to D310 next:      
      D31= (-2*d311*(q23**2 - q2s*q3s) + q23**2*r43 - q2s*q3s*r43 +
     &    q12*q3s*r44 + q13*q2s*r45 - q23*(q13*r44 + q12*r45))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))
     
      D32 = (-2*d312*(q13**2 - q1s*q3s) + q12*q3s*r46 + q13**2*r47 -
     &    q1s*q3s*r47 + q1s*q23*r48 - q13*(q23*r46 + q12*r48))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D33 = (-2*d313*(q12**2 - q1s*q2s) + q13*q2s*r49 + q1s*q23*r50 -
     &    q12*(q23*r49 + q13*r50) + q12**2*r51 - q1s*q2s*r51)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D34 = (2*d311*(q13*q23 - q12*q3s) + q12*q3s*r43 + q13**2*r44 -
     &    q1s*q3s*r44 + q1s*q23*r45 - q13*(q23*r43 + q12*r45))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D35 =(2*d311*(q12*q23 - q13*q2s) + q13*q2s*r43 + q1s*q23*r44 -
     &    q12*(q23*r43 + q13*r44) + q12**2*r45 - q1s*q2s*r45)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D36 = (2*d312*(q13*q23 - q12*q3s) + q23**2*r46 - q2s*q3s*r46 +
     &    q12*q3s*r47 + q13*q2s*r48 - q23*(q13*r47 + q12*r48))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D37 = (2*d313*(q12*q23 - q13*q2s) + q23**2*r49 - q2s*q3s*r49 +
     &    q12*q3s*r50 + q13*q2s*r51 - q23*(q13*r50 + q12*r51))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D38 = (2*d312*(q12*q13 - q1s*q23) + q13*q2s*r46 + q1s*q23*r47 -
     &    q12*(q23*r46 + q13*r47) + q12**2*r48 - q1s*q2s*r48)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D39 = (2*d313*(q12*q13 - q1s*q23) + q12*q3s*r49 + q13**2*r50 -
     &    q1s*q3s*r50 + q1s*q23*r51 - q13*(q23*r49 + q12*r51))/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      D310 = (d311*q12*q13 + d312*q12*q23 - d311*q1s*q23 - d312*q13*q2s -
     &    q12*q23*r52 + q13*q2s*r52 - q12*q13*r53 + q1s*q23*r53 +
     &    q12**2*r54 - q1s*q2s*r54)/
     &  (-2*q12*q13*q23 + q13**2*q2s + q12**2*q3s + q1s*(q23**2 - q2s*q3s))

      Dij(3,1)  = D31
      Dij(3,2)  = D32
      Dij(3,3)  = D33
      Dij(3,4)  = D34
      Dij(3,5)  = D35
      Dij(3,6)  = D36
      Dij(3,7)  = D37
      Dij(3,8)  = D38
      Dij(3,9)  = D39
      Dij(3,10) = D310
      Dij(3,11) = D311
      Dij(3,12) = D312
      Dij(3,13) = D313
  
 
  
      end
