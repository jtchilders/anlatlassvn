c this file contains the following scalar triangle functions:
c
c	short:			Binoth's notation:
c
c	C0i1e(s,mus)	 		I3d(s,0,0;0,0,0)
c	C0i2e(s1,s2,mus) 		I3d(s1,s2,0;0,0,0)
c	C0i3e(s1,s2,s3,mus)		I3d(s1,s2,s3;0,0,0)     
c
c	C1i1e(s,mis,mus)		I3d(s,0,0;0,mis,0)
c	C1i2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,mis,0
c	C1d2e(s1,s2,mis,mus)		I3d(s1,s2,0;0,0,mis)
c	C2i2e(s1,s2,m2s,m3s,mus)	I3d(s1,s2,0;0,mis,mis)	

c	C2i3e(s1,s2,s3,m2s,m3s,mus)	I3d(s1,s2,s3;0,mis,mis) 
c      
c      IR divergent integrals are taken from Dittmaier hep-ph/0308246;
c      other integrals taken from Binoth et al., 0709.3513 [hep-ph]    
c      
c      
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i1e(s,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: one external mass scale s
c	      no  internal mass	
c
c		g_1e = 1/s
c
c		A_1e = 1
c		B_1e = log(4*pi*mus/(-s))
c
c		C_1e = 	1/2*log(4*pi*mus/(-s))**2 - pi^2/6
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s,0,0;0,0,0)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,mus
      double precision fpm
      double complex sbar
      double complex c1e
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
          
      sbar = s+ieps
      fpm = fpi*mus
c
c for comparison:
c	fpm = mus

      	c1e = 0.5d0*(log((-sbar)/fpm))**2-pi2o6
     
	C0i1e = c1e/s
		            
      return
      end      
      
      
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i2e(s1,s2,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      no  internal mass	
c
c		g_2e = 1/(s1-s2)
c
c		A_2e = 0
c		B_2e = log((-s2)/(-s1))
c
c		C_2e =	 1/2*log((-s1)/4*pi*mus)**2 
c			-1/2*log((-s2)/4*pi*mus)**2 
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s1,s2,0;0,0,0)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mus
      double precision fpm
      double complex s1bar,s2bar
      double complex c2e
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
          
      s1bar = s1+ieps
      s2bar = s2+ieps
      fpm = fpi*mus
c
c for comparison:
c	fpm = mus

      	c2e = (log((-s1bar)/fpm))**2-(log((-s2bar)/fpm))**2
     
	C0i2e = 0.5d0*c2e/(s1-s2)
		            
      return
      end      
      
            
c---------  scalar triangle function: massless internal lines  --------------
c
      complex*16 function C0i3e(s1,s2,s3,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: three external mass scales s1,s2,s3

c		A = 0
c		B = 0
c	     
c	taken from [Binoth et al.], see (A.7)
c
c	Binoth's notation: I3d(s1,s2,s3;0,0,0)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,s3,mus
      double complex xp,xm,yp,ym
      double complex gki,c1
      
      double precision zero
      parameter (zero=0d0)
      double precision lam
      complex*16 ieps,li2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc,lam,eta
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o3=pi**2/3.d0)
          
      xp = (s1+s3-s2-sqrt(lam(s1,s2,s3)-ieps*s1))/(2.d0*s1)
      xm = (s1+s3-s2+sqrt(lam(s1,s2,s3)-ieps*s1))/(2.d0*s1)
      
      yp = 1.-xm
      ym = 1.-xp	

      c1 =  2.d0*li2(-xm/yp)+2.d0*li2(-ym/xp)+pi2o3
     &	   +0.5d0*(log(xm/yp))**2+0.5d0*(log(ym/xp))**2  
     &	   +0.5d0*(log(xp/yp))**2-0.5d0*(log(xm/ym))**2  
            
      gki = sqrt(lam(s1,s2,s3)-ieps*s1)	
	
      C0i3e = -c1/gki
		            
      return
      end      
     
       
     
         
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1i1e(s,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi0s,mi1s,mi2s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi0s][(q+q1)^2-mi1s][(q+q1+q2)^2-mi2s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: one external mass scale s attached to
c	      one internal mass	scale mis
c
c		g_1i1e = 1/s
c
c		A_1i1e = 0
c		B_1i1e = log(mis/(mis-s))
c
c		C_1i1e = log(4*pi*mus/(mis))*log(mis/(mis-s))
c			+log((mis-s)/mis)**2
c			+li2((s)/mis)
c
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(sq,0,0;0,mis,0)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,mis,mus
      double precision fpm
      double complex sbar
      double complex c1
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi)
          
      sbar = s+ieps
      fpm = fpi*mus
c
c for comparison:
c	fpm = mus

      	c1 = log(fpm/(mis))*log(mis/(mis-sbar))
     &	    +log((mis-sbar)/mis)**2
     &	    +li2(sbar/mis)
     
	C1i1e = c1/s
		            
      return
      end      
      
       
         
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1i2e(s1,s2,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      one internal mass	scale mis in between
c
c		g_1i2e = 1/(s1-s2)
c
c		A_1i2e = 0
c		B_1i2e = log((mis-s2)/(mis-s1))
c
c		C_1i2e = log(4*pi*mus/(mis))*log((mis-s2)/(mis-s1))
c			+log((mis-s1)/mis)**2
c			-log((mis-s2)/mis)**2
c			+li2((s1)/mis)-li2((s2)/mis)
c
c
c	     
c	(see my notes master-integrals.tex)
c
c	Binoth's notation: I3d(s1,s2,0;0,mis,0)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mis,mus
      double precision fpm
      double complex s1bar,s2bar
      double complex c1
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi)
          
      s1bar = s1+ieps
      s2bar = s2+ieps
      fpm = fpi*mus
c
c for comparison:
c	fpm = mus

      	c1 = log(fpm/(mis))*log((mis-s2bar)/(mis-s1bar))
     &	    +(log((mis-s1bar)/mis))**2-(log((mis-s2bar)/mis))**2
     &	    +li2(s1bar/mis)-li2(s2bar/mis)
     
	C1i2e = c1/(s1-s2)
		            
      return
      end      
     
   
         
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C1d2e(s1,s2,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      one internal mass	scale mis on opposite side
c
c	attention: not symmetric in s1,s2!!
c		
c		g_1d2e = 1/(s1-s2)
c
c		A_1d2e = 0
c		B_1d2e = 0
c
c		C_1d2e =
c
c
c	     
c	taken from [Binoth et al.]
c
c	Binoth's notation: I3d(s1,s2,0;0,0,mis)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,mis,mus
      double precision x0,fpm
      double complex x1bar,s1bar
      double complex c1
      
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2,rfunc
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc
      double precision pi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6.d0)
          
      s1bar = s1+ieps
      
      x0 = s1/(s1-s2)
      x1bar = s1bar/(s1-s2+mis)

      c1 = rfunc(dcmplx(x0),x1bar)-pi2o6+li2(1.-1/(x0-ieps))
	
      C1d2e = c1/(s2-s1)
		            
      return
      end      
         
           
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C2i2e(s1,s2,m2s,m3s,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: two external mass scales s1,s2
c	      two equal internal mass scales m3s=m2s 

c		A = 0
c		B = 0
c	
c	(does not depend on mus up to order eps^0)
c
c	     
c	taken from [Binoth et al.], (A.11)
c
c	Binoth's notation: I3d(s1,s2,0;0,mis,mis)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,m2s,m3s,mis,mus
      double complex s1bar,s2bar
      double complex xp,xm,eta0
      double complex c1,clog
      double precision gki,x0
      
      double precision zero
      parameter (zero=0d0)
      double complex czero
      parameter (czero=(0d0,0d0))
      complex*16 ieps,li2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc,eta
          
      if (m3s.ne.m2s) then
      	print*,'C2i2e is implemented for two equal internal mases only!'
	print*,'set C2i2e = 0'
	C2i2e = 0d0
	return
      endif	
      
      mis = m2s
	
      s2bar = s2+ieps
      s1bar = s1+ieps
      
      x0 = 1.d0-s1/s2
      
      xp = 0.5d0*(1+sqrt(1d0-4*mis/s2bar)) 	
      xm = 0.5d0*(1-sqrt(1d0-4*mis/s2bar)) 	
	
      eta0 = eta(1d0-(s1bar/mis)*x0,mis/(mis-s2bar*x0*(1d0-x0))) 
      
      if (eta0.eq.czero) then
        clog = czero
      else
        clog = -eta0*log((1d0-x0)/(-x0))
      endif		
 	
 
      c1 = li2(s1bar/mis)-li2(1d0/xp)-li2(1d0/xm)
     &    +rfunc(dcmplx(x0),mis/s1bar)+rfunc(dcmplx(1d0-x0),xm)
     &    -rfunc(dcmplx(x0),xm)
      
      gki = s2-s1
	
      C2i2e = (c1+clog)/gki
      		            
      return
      end      
     
           
c---------  scalar triangle function: massive internal lines  --------------
c
      complex*16 function C2i3e(s1,s2,s3,m2s,m3s,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c
c	here: three external mass scales s1,s2,s3
c	      two equal internal mass scales m2s=m3s 

c		A = 0
c		B = 0
c	     
c	taken from [Binoth et al.]
c
c	Binoth's notation: I3d(s1,s2,s3;0,mis,mis)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s1,s2,s3,m3s,m2s,mis,mus
      double complex x0bar,x1bar,s1bar,s2bar
      double complex xp,xm,etap,etam
      double complex gki,c1
      
      double precision zero
      parameter (zero=0d0)
      double precision lam
      complex*16 ieps,li2,rfunc,eta
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc,lam,eta
      double precision pi,fpi,pi2o6
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi,pi2o6=pi**2/6.d0)
          
      if (m3s.ne.m2s) then
      	print*,'C2i3e is implemented for two equal internal mases only!'
	print*,'set C2i3e = 0'
	C2i3e = 0d0
	return
      endif	
	
      s2bar = s2+ieps
      s1bar = s1+ieps
     	
      xp = (s1+s2-s3-sqrt(lam(s1,s2,s3)-ieps*s2))/(2.d0*s2)	
      xm = (s1+s2-s3+sqrt(lam(s1,s2,s3)-ieps*s2))/(2.d0*s2)	
      
      x0bar = 0.5d0*(1.d0-sqrt(1.-4*m3s/s2bar))
      x1bar = (s1bar-m3s)/(s1-s3)
      
      etap = eta(-s3*xp-s1*(1.-xp)+m3s-ieps,1/(m3s-xp*(1.-xp)*s2bar))
      etam = eta(-s3*xm-s1*(1.-xm)+m3s-ieps,1/(m3s-xm*(1.-xm)*s2bar))
      

      c1 = rfunc(xm,x1bar)-rfunc(xp,x1bar)+rfunc(1.-xm,x0bar)
     &	  -rfunc(1.-xp,x0bar)-rfunc(xm,x0bar)+rfunc(xp,x0bar)
     &    -etam*log((1.-xm)/(-xm))+etap*log((1.-xp)/(-xp))      
      
      gki = sqrt(lam(s1,s2,s3)-ieps*s2)	
	
      C2i3e = c1/gki
		            
      return
      end      
     
      
c ------------------------------------------------------------------------         
c---------  auxiliary functions for computation of loop integrals  --------
c
      complex*16 function rfunc(y0,z) 
c	     
c	taken from [Binoth et al.], see Eq.(A.4)
c
c ------------------------------------------------------------------------
      implicit none

      double complex y0,z
      
      double complex z1,z2,eta1,eta2
      
      double precision zero
      parameter (zero=0d0)
      complex*16 li2,eta
      external li2,eta     
 
      z1 = y0/(y0-z)
      z2 = (y0-1.d0)/(y0-z)
            
      eta1 = eta(dcmplx(-z),dcmplx(1/(y0-z)))
      eta2 = eta(dcmplx(1.-z),dcmplx(1/(y0-z)))
      
      rfunc = li2(dcmplx(z1))-li2(dcmplx(z2))+eta1*log(z1)-eta2*log(z2)

      
		            
      return
      end      
         
c-------------------------------------------------------------------
c
      double precision function lam(x,y,z) 
c	     
c	taken from [Binoth et al.], see Eq.(A.2)
c
c -------------------------------------------------------------------
      implicit none
      
      double precision x,y,z
 
      lam = x**2+y**2+z**2-2*x*y-2*y*z-2*x*z
		            
      return
      end  
                
c-------------------------------------------------------------------------
***********************************************************************
        FUNCTION ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
c        IMPLICIT   LOGICAL(A-Z)
	
	implicit none
	                                        
        COMPLEX*16 ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
 
	if (((IM1.eq.0d0).and.(DREAL(C1).lt.0d0)).or.
     &	    ((IM2.eq.0d0).and.(DREAL(C2).lt.0d0)).or.
     &	    ((IM12.eq.0d0).and.(DREAL(C1*C2).lt.0d0))) then
	  write(*,*) 'eta function on cut !!!'
	  write(*,*) 'C1    = ',C1
	  write(*,*) 'C2    = ',C2
	  write(*,*) 'C1*C2 = ',C1*C2
	  stop
	endif
                                                                      
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            
