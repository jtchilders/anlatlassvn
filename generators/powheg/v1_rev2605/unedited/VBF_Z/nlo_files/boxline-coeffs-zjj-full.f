*******************************************************************
*
*     this subroutine gives coefficients of the various Lorentz structures
*     contributing to selfenergy+vertex+box corrections to 
*     quark line which emits two massive/virtual vector bosons
*
*     q(p1)->V1(q1) V2(q2) q(p2)
c
c
*******************************************************************

      subroutine boxcorr(p1,q1,q2,p2,musq,c1,c2,c3,c4)
      
c       IN:
c	p1,q1,q2,p2: "external" momenta (p1->q1+q2+p2)
c	musq: regularization scale squared
c
c	OUT:
c	c1(si),c2(rho),c3(si,rho,om),c4(om): 
c		coefficients of the various Lorentz structres	
c
c	ATTENTION: Born-type contributions are subtracted;
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      implicit none
      
      double precision p1(0:3),q1(0:3),q2(0:3),p2(0:3),musq
      double precision mq1(0:3),mp2(0:3),mq2(0:3),mq1mq2(0:3),p1mq1(0:3)
      integer mu,nu,si,rho,om,i,j
      double precision gmnr(0:3,0:3)

      double precision p12,p1q2,p2q2,q2s,q1s,q12
      double precision p12sq,q2ssq,p2q2sq,p1q2sq
      double precision p12tr,q2str,p2q2tr,p1q2tr
      double precision mp12
       
      double complex Cij(0:2,4),Dij(0:3,13)
      
      double precision dotrr
      external dotrr
      
      double complex C0fin,B0fin,D0fin
      external C0fin,B0fin,D0fin
      
      double complex B0q1,B0q2,B0pq 
      
      double complex C0123,C0124,C0134,C0234,D00       
      double complex Cij_123(0:2,4),Cij_124(0:2,4),
     &               Cij_134(0:2,4),Cij_234(0:2,4)    
      double complex C0cr,Cijcr(0:2,4)    
      
      double complex Cij12311,Cij12312,Cij12321,Cij12322,Cij12323,Cij12324
      double complex Cij12411,Cij12412,Cij12421,Cij12422,Cij12423,Cij12424
      double complex Cij13411,Cij13412,Cij13421,Cij13422,Cij13423,Cij13424
      double complex Cij23411,Cij23412,Cij23421,Cij23422,Cij23423,Cij23424
      double complex Cijcr11,Cijcr12,Cijcr21,Cijcr22,Cijcr23,Cijcr24

      double complex Dij11,Dij12,Dij13
      double complex Dij21,Dij22,Dij23,Dij24,Dij25,Dij26,Dij27
      double complex Dij31,Dij32,Dij33,Dij34,Dij35,Dij36,Dij37,
     % 		     Dij38,Dij39,Dij310,Dij311,Dij312,Dij313
      
      double complex c(4,10)

      double complex c1(0:3),c2(0:3),c3(0:3,0:3,0:3),c4(0:3)
      
      double precision small
      parameter (small = 1d-9) ! define as 'zero'
      
      logical box_type,vtx_type
      parameter (box_type = .true.,vtx_type=.true.)
      
      ! take from global.inc later:
      double precision pi,pi2o3
      parameter (pi = 3.141592653589793d0, pi2o3= pi**2/3d0)
c local parameter (don't mix with cvirt!)
      double precision bvirt 
      parameter (bvirt = -7d0+pi2o3) 
      double precision fpi
      parameter (fpi=4d0*pi)

*
*******************************************************************

c need to check momentum conservation here (test sign of momenta)
	
	p12  = dotrr(p1(0),p2(0))
	p1q2 = dotrr(p1(0),q2(0))
	p2q2 = dotrr(p2(0),q2(0))
	q2s  = dotrr(q2(0),q2(0))
	q1s  = dotrr(q1(0),q1(0))
	q12  = dotrr(q1(0),q2(0))
	
	do mu = 0,3
	   mq1(mu) = -q1(mu)
	   mq2(mu) = -q2(mu)
	   mp2(mu) = -p2(mu)
	   mq1mq2(mu) = -q1(mu)-q2(mu)
	   p1mq1(mu) = p1(mu)-q1(mu)
	enddo 
	mp12 = -p12
	
	
	if (abs(q1s).lt.small) q1s = 0d0  	
	if (abs(q2s).lt.small) q2s = 0d0  	
	 
	p12sq  = p12**2 
	q2ssq  = q2s**2
	p2q2sq = p2q2**2 
	p1q2sq = p1q2**2  
	
	p12tr  = p12*p12sq 
	q2str  = q2s*q2ssq
	p2q2tr = p2q2*p2q2sq 
	p1q2tr = p1q2*p1q2sq 
	
c adjust scale to comply with normalization of integrals:
	musq = musq/fpi	
	
	if (box_type) then 
	 ! D00 = D0((l1+l2)^2,(l2+l3)^2,l1s,l2s,l3s,l4s,musq)
	 ! for p1 -> q1+q2+p2
c        D00 = D0fin(2d0*p2q2+q2s,q1s+q2s+2d0*q12,0d0,q1s,q2s,0d0,musq)	
	
	call tens_box(p1,mq1,mq2,musq,Dij)
		
	D00   = Dij(0,1)	
	Dij11 = Dij(1,1)
	Dij12 = Dij(1,2)
	Dij13 = Dij(1,3)
	Dij21 = Dij(2,1)
	Dij22 = Dij(2,2)
	Dij23 = Dij(2,3)
	Dij24 = Dij(2,4)
	Dij25 = Dij(2,5)
	Dij26 = Dij(2,6)
	Dij27 = Dij(2,7)
		
	Dij31 = Dij(3,1)
	Dij32 = Dij(3,2)
	Dij33 = Dij(3,3)
	Dij34 = Dij(3,4)
	Dij35 = Dij(3,5)
	Dij36 = Dij(3,6)
	Dij37 = Dij(3,7)
	Dij38 = Dij(3,8)
	Dij39 = Dij(3,9)
	Dij310 = Dij(3,10)
	Dij311 = Dij(3,11)
	Dij312 = Dij(3,12)
	Dij313 = Dij(3,13)
	 
	
	endif !box_type
     	
	call tens_tri(p1,mq1,musq,Cij_123)		!l1 = p1,l2 = -q1
c	C0123 = C0fin(0d0,q1s,2d0*p2q2+q2s,musq)

	C0123    = Cij_123(0,1)
	Cij12311 = Cij_123(1,1)
	Cij12312 = Cij_123(1,2)
	Cij12321 = Cij_123(2,1)
	Cij12322 = Cij_123(2,2)
	Cij12323 = Cij_123(2,3)
	Cij12324 = Cij_123(2,4)	
	
	
	call tens_tri(p1,mq1mq2,musq,Cij_124)		!l1 = p1,l2 = -q1-q2
c	C0124 = C0fin(0d0,q1s+q2s+2d0*q12,0d0,musq)
	
	C0124    = Cij_124(0,1)
	Cij12411 = Cij_124(1,1)
	Cij12412 = Cij_124(1,2)
	Cij12421 = Cij_124(2,1)
	Cij12422 = Cij_124(2,2)
	Cij12423 = Cij_124(2,3)
	Cij12424 = Cij_124(2,4)
	

	call tens_tri(p1mq1,mq2,musq,Cij_134)		!l1 = p1-q1,l2 = -q2
c	C0134 = C0fin(2d0*p2q2+q2s,q2s,0d0,musq)
	
	C0134    = Cij_134(0,1)
	Cij13411 = Cij_134(1,1)
	Cij13412 = Cij_134(1,2)
	Cij13421 = Cij_134(2,1)
	Cij13422 = Cij_134(2,2)
	Cij13423 = Cij_134(2,3)
	Cij13424 = Cij_134(2,4)

	call tens_tri(mq1,mq2,musq,Cij_234)		!l1 = -q1,l2 = -q2
c	C0234 = C0fin(q1s,q2s,2d0*mp12,musq)
	
	C0234    = Cij_234(0,1)
	Cij23411 = Cij_234(1,1)
	Cij23412 = Cij_234(1,2)
	Cij23421 = Cij_234(2,1)
	Cij23422 = Cij_234(2,2)
	Cij23423 = Cij_234(2,3)
	Cij23424 = Cij_234(2,4)
	
	
	call tens_tri(mp2,mq2,musq,Cijcr)		!l1 = -p2,l2 = -q2
c	C0cr = C0fin(0d0,q2s,2d0*p2q2+q2s,musq)
	
	C0cr    = Cijcr(0,1)
	Cijcr11 = Cijcr(1,1)
	Cijcr12 = Cijcr(1,2)
	Cijcr21 = Cijcr(2,1)
	Cijcr22 = Cijcr(2,2)
	Cijcr23 = Cijcr(2,3)
	Cijcr24 = Cijcr(2,4)
		
	B0q1 = B0fin(q1s,musq)
        B0q2 = B0fin(q2s,musq) 
	B0pq = B0fin(2*p2q2 + q2s,musq)

c initialize:
	do i = 1,4
	do j = 1,10
	    c(i,j) = 0d0
	enddo
	enddo	

c Lorentz structure L1 = ubar(p2).gamma_rho.u(p1) 
c comes with coefficient 
c	C1.{si} = p1.{si}*c(1,1) + p2.{si}*c(1,2) + q2.{si}*c(1,3)	
c
      if (vtx_type) then
      c(1,1) = 0d0

      c(1,2) = ((-2*(4*Cij12324 + 4*Cijcr24 - 4*C0123*p12 - 
     &  4*Cij12311*p12 - 4*C0123*p1q2 - 4*Cij12311*p1q2 -
     &  4*Cij12312*p2q2 + 4*Cijcr23*p2q2 + 2*C0cr*q2s - 
     &  2*Cij12312*q2s + 2*Cijcr11*q2s + 2*Cijcr23*q2s - 2*B0q1 - 
     &  2*B0q2 + B0pq))/(2*p2q2 + q2s)+
     &  (-2*bvirt)/(2*p2q2 + q2s))

      c(1,3) = -4*(Cijcr12 + Cijcr22)
     
      endif !vtx-type
      
      if (box_type) then
c"p1-coeff:"
	c(1,1) = c(1,1) + (
     &   -4*(C0234 - 2*Dij27 + 2*Dij311 - 2*Dij312 + 2*Dij22*p12 - 
     &   2*Dij24*p12 + 2*Dij22*p1q2 - 
     &   2*Dij24*p1q2 + 2*Dij25*p1q2 - 2*Dij26*p1q2 + 2*Dij12*p2q2 - 
     &   2*Dij13*p2q2 - 2*Dij22*p2q2 + 2*Dij24*p2q2 + Dij12*q2s - 
     &   Dij13*q2s - Dij22*q2s + Dij24*q2s - Dij25*q2s + Dij26*q2s))
    
c"p2-coeff:"
	c(1,2) = c(1,2) + (
     &  -4*(-C0234 - Cij23411 + 2*Dij312 - 2*D00*p12 - 2*Dij11*p12 - 
     &  2*Dij12*p12 - 2*Dij24*p12 -  
     &  2*Dij12*p1q2 + 2*Dij13*p1q2 - 2*Dij24*p1q2 + 
     &  2*Dij26*p1q2 - 2*Dij12*p2q2 + 2*Dij13*p2q2 - 2*Dij22*p2q2 + 
     &  2*Dij26*p2q2 - Dij12*q2s + Dij13*q2s - Dij22*q2s + Dij26*q2s))
    
c"q2-coeff: "
	c(1,3) = c(1,3) + (
     &  -4*(-C0234 + 4*Dij27 + 2*Dij312 - 2*Dij313 - 2*Dij12*p12 + 
     &  2*Dij13*p12 - 2*Dij22*p12 + 2*Dij26*p12 - 
     &  2*Dij12*p1q2 + 2*Dij13*p1q2 - 2*Dij22*p1q2 - 2*Dij23*p1q2 + 
     &  4*Dij26*p1q2 + 2*Dij22*p2q2 - 2*Dij26*p2q2 + Dij22*q2s + 
     &  Dij23*q2s - 2*Dij26*q2s))
      
      
      endif!box-type

c	C1.{si} = p1.{si}*c11 + p2.{si}*c12 + q2.{si}*c13	
	do si = 0,3
	   c1(si) = p1(si)*c(1,1) + p2(si)*c(1,2) + q2(si)*c(1,3)
	enddo
	
c -----------------	
c	
c Lorentz structure L2 = ubar(p2).gamma_si.u(p1) 
c comes with coefficient 
c	C2.{rho} = p1.{rho}*c21 + p2.{rho}*c22 + q2.{rho}*c23	
c
      if (vtx_type) then
      c(2,1) = -4*(C0123 + Cij12311 - Cij12322 + Cij12323)
  
      c(2,2) = -4*(Cij12312 + Cij12322)

      c(2,3) = -4*(Cij12312 + Cij12322)
      endif ! vtx-type
      
      if (box_type) then
c  "p1-coeff:"
	c(2,1) = c(2,1) + (
     &   -4*(Cij23411 + 2*Dij27 + 2*Dij311 - 2*Dij312 - 
     &   2*Dij11*p12 + 2*Dij12*p12 - 2*Dij21*p12 - 2*Dij22*p12 + 
     &   4*Dij24*p12 - 2*Dij12*p2q2 + 2*Dij13*p2q2 + 4*Dij22*p2q2 - 
     &   4*Dij24*p2q2 + 2*Dij25*p2q2 - 2*Dij26*p2q2 - Dij12*q2s + 
     &   Dij13*q2s + Dij22*q2s - Dij24*q2s + Dij25*q2s - Dij26*q2s))
    
c"p2-coeff:"
	c(2,2) = c(2,2) + (
     &   -4*(2*Dij27 + 2*Dij312 - 2*Dij22*p1q2 + 2*Dij24*p1q2 + 
     &   Dij12*q2s - Dij13*q2s + Dij22*q2s - Dij26*q2s))
    
c"q2-coeff: "
	c(2,3) = c(2,3) + (
     &   -4*(-C0234 + 4*Dij27 + 2*Dij312 - 2*Dij313 - 2*Dij22*p1q2 + 
     &   2*Dij24*p1q2 - 2*Dij25*p1q2 + 2*Dij26*p1q2 - 
     &   2*Dij12*p2q2 + 2*Dij13*p2q2 + Dij22*q2s + 
     &   Dij23*q2s - 2*Dij26*q2s))
    
      endif!box-type
	
c	C2.{rho} = p1.{rho}*c21 + p2.{rho}*c22 + q2.{rho}*c23	
	do rho = 0,3
	   c2(rho) = p1(rho)*c(2,1) + p2(rho)*c(2,2) + q2(rho)*c(2,3)
	enddo
c
c -----------------	
c	
c Lorentz structure L3 = ubar(p2).gamma_om.u(p1) 
c comes with coefficient 
c C3.{si,rho,om} = (g_{si}.{rho}*c31     + q2.{si}*q2.{rho}*c32 +
c		   q2.{si}*p2.{rho}*c33 + q2.{si}*p1.{rho}*c34 +
c		   p2.{si}*q2.{rho}*c35 + p1.{si}*q2.{rho}*c36 +
c		   p2.{si}*p2.{rho}*c37 + p1.{si}*p1.{rho}*c38 +
c		   p1.{si}*p2.{rho}*c39 + p2.{si}*p1.{rho}*c310)*q2.{omm} 

c there are no vtx-type contributions 

	if (box_type) then
	
c"g-coeff:"
	c(3,1) = -4*(C0234 - Cij23411 + Cij23412 - 4*Dij27 + 
     &   2*Dij312 - 2*Dij313 + 2*Dij22*p12 - 
     &   2*Dij24*p12 + 2*Dij25*p12 - 2*Dij26*p12 + 2*Dij22*p1q2 - 
     &   2*Dij24*p1q2 + 2*Dij25*p1q2 - 2*Dij26*p1q2 - 4*Dij22*p2q2 - 
     &   2*Dij23*p2q2 + 6*Dij26*p2q2 - 2*Dij22*q2s - 2*Dij23*q2s + 
     &	 4*Dij26*q2s)

c"q*q-coeff:"
	c(3,2) = -8*(Dij22 + Dij23 - 2*Dij26 + Dij32 - 
     &		  Dij33 - 3*Dij38 + 3*Dij39)

c"q2.{si}*p2.{rho}:"
	c(3,3) = -8*(Dij22 + Dij23 - 2*Dij26 + Dij32 - 2*Dij38 + Dij39)

c"q2.{si}*p1.{rho}:"
	c(3,4) = -8*(Dij12 - Dij13 + Dij23 + Dij24 - Dij25 - Dij26 - 
     &		  2*Dij310 - Dij32 + Dij36 + Dij37 + 2*Dij38 - Dij39)

c"p2.{si}*q2.{rho}:"
	c(3,5) = -8*(Dij12 - Dij13 + 2*Dij22 + Dij23 - 
     &		  3*Dij26 + Dij32 - 2*Dij38 + Dij39)

c"p1.{si}*q2.{rho}:"
	c(3,6) = 8*(2*Dij310 + Dij32 - Dij36 - Dij37 - 2*Dij38 + Dij39)

c"p2.{si}*p2.{rho}:"
	c(3,7) = -8*(Dij12 - Dij13 + 2*Dij22 - 2*Dij26 + Dij32 - Dij38)

c"p1.{si}*p1.{rho}:"
	c(3,8) = 8*(Dij22 - Dij24 + Dij25 - Dij26 - 2*Dij310 - 
     &		 Dij32 - Dij34 + Dij35 +  2*Dij36 + Dij38)

c"p1.{si}*p2.{rho}:"
	c(3,9) = 8*(Dij25 - Dij26 + Dij310 + Dij32 - Dij36 - Dij38)

c"p2.{si}*p1.{rho}:"
	c(3,10) = -8*(Dij12 - Dij13 - Dij22 + 2*Dij24 - Dij25 - 
     &		   Dij310 - Dij32 + Dij36 + Dij38)

	endif !box-type

c C3.{si,rho,om} = (g_{si}.{rho}*c31     + q2.{si}*q2.{rho}*c32 +
c		   q2.{si}*p2.{rho}*c33 + q2.{si}*p1.{rho}*c34 +
c		   p2.{si}*q2.{rho}*c35 + p1.{si}*q2.{rho}*c36 +
c		   p2.{si}*p2.{rho}*c37 + p1.{si}*p1.{rho}*c38 +
c		   p1.{si}*p2.{rho}*c39 + p2.{si}*p1.{rho}*c310)*q2.{omm}  
        
	do mu = 0,3
	   do nu = 0,3
	      gmnr(mu,nu) = 0d0
	   enddo
	enddo           
	gmnr(0,0) =  1d0
        gmnr(1,1) = -1d0
        gmnr(2,2) = -1d0
        gmnr(3,3) = -1d0

	do si = 0,3
	   do rho = 0,3
	      do om = 0,3
	 	 c3(si,rho,om) = (
     &		    gmnr(si,rho)*c(3,1) + q2(si)*q2(rho)*c(3,2) +
     &		  q2(si)*p2(rho)*c(3,3) + q2(si)*p1(rho)*c(3,4) +
     &		  p2(si)*q2(rho)*c(3,5) + p1(si)*q2(rho)*c(3,6) +
     &		  p2(si)*p2(rho)*c(3,7) + p1(si)*p1(rho)*c(3,8) +
     &		  p1(si)*p2(rho)*c(3,9) + p2(si)*p1(rho)*c(3,10))*q2(om)	     
	      enddo
	   enddo
	enddo
c
c
c -----------------	
c	
c Lorentz structure L4 = ubar(p2).gamma_si.gamma_om.gamma_rho.u(p1) 
c comes with coefficient 
c	C4.{om} =  q2.{om}*c4	

	if (vtx_type) then
	
       c(4,1) = ((-4*Cij12324 - 4*Cijcr24 + 4*C0123*p12 + 
     &   4*Cij12311*p12 + 4*C0123*p1q2 + 4*Cij12311*p1q2 + 
     &   4*C0cr*p2q2 + 4*Cij12312*p2q2 + 4*Cijcr11*p2q2 + 
     &   4*Cijcr12*p2q2 + 2*Cij12312*q2s + 2*Cijcr12*q2s + 2*B0q1 + 
     &   2*B0q2 - B0pq)/(2*p2q2 + q2s) +
     &	 (-bvirt/(2*p2q2 + q2s)))
    
     	endif !(vtx_type)
	
	if(box_type) then
	
	c(4,1) = c(4,1) + (
     &   2*(-C0234 + Cij23411 - Cij23412 + 6*Dij27 + 2*D00*p12 + 
     &	 2*Dij11*p12 - 2*Dij22*p12 + 2*Dij24*p12 - 2*Dij25*p12 + 
     &   2*Dij26*p12 - 2*Dij22*p1q2 + 2*Dij24*p1q2 - 2*Dij25*p1q2 + 
     &   2*Dij26*p1q2 + 4*Dij22*p2q2 + 2*Dij23*p2q2 - 6*Dij26*p2q2 + 
     &   2*Dij22*q2s + 2*Dij23*q2s - 4*Dij26*q2s))
    
     	endif !(box_type)
	
     
        do om = 0,3
	   c4(om) = q2(om)*c(4,1)
	enddo
 	    
        return
	end
	
