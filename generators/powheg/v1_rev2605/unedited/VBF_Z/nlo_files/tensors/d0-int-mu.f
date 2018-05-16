c this file contains the following scalar box functions:
c
c	short:			   Binoth's notation:
c
c    D02mh(s,t,q3s,q4s,mus)	   I4d(q3s,q4s,0,0;s,t;0,0,0,0) changed order!
c    D02me(s,t,q2s,q4s,mus)	   I4d(q2s,0,q4s,0;t,s;0,0,0,0)
c    D01m(s,t,q4s,mus)  	   I4d(q4s,0,0,0;t,s;0,0,0,0)
c
c    D1i1e(s,t,q3s,mis,mus)	   I4d(q3s,0,0,0;s,t;0,mis,0,0)
c    D2i2e(s,t,q2s,q3s,mis,mus)    I4d(q2s,q3s,0,0;t,s;0,mis,mis,0)
c    D2i2d(s,t,q1s,q3s,mis,mus)    I4d(q1s,0,q3s,0;s,t;0,0,mis,mis)
c
c  boxes with massless propagators are taken from Duplancic, Nizic;
c  massive boxes are taken from Binoth et al., 0709.3513 [hep-ph] ;
c
c	typo in (A.19) for D2i2e is corrected
c
c      
c---------  scalar box function: massless internal lines  --------------
c
      complex*16 function D02mh(s,t,q3s,q4s,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c	     
c	     mi1s=mi2s=mi3s=mi4s=0 and q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	case a) two adjacent massive external legs ("2mh")
c		-> q1s=0,q2s=0,q3s.ne.0,q4s.ne.0
c
c		 g_2mh = 1/s/t
c		 f_2mh = (s+t-q4s)/(s*t)
c
c		 A_2mh = 1
c		 B_2mh =  ln[(-q4s-ieps)/(-s-ieps)]
c			 -ln[(-t-ieps)/(4*pi*mus)]
c	 		 +ln[(-q3s-ieps)/(-t-ieps)]
c
c		 C1_2mh = 2*Li2[1-q4s/s]-2*Li2[1-q3s/t] 
c			 +2*Li2[1-s*f_2mh]+2*Li2[1-t*f_2mh] 
c			 -2*Li2[1-q4s*f_2mh] 
c
c		 C2_2mh = -pi^2/6
c			  +1/2*ln^2[-s/(4*pi*mus)]
c			      +ln^2[-t/(4*pi*mus)]
c	 		  -1/2*ln^2[(-q4s)/(4*pi*mus)]
c			  -1/2*ln^2[(-q3s)/(4*pi*mus)]
c			  +1/2*ln^2[(-q4s)/(-s)]
c			  -1/2*ln^2[(-q3s)/(-q4s)]
c			  +1/2*ln^2[(-q3s)/(-s)]
c
c		(disregarded ieps here for simplicity)	
c
c	taken from Duplancic, Nizic;
c	(see my notes master-integrals.tex)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q3s,q4s,mus
      double precision f2mh,fpm
      double complex sbar, tbar, q3sbar, q4sbar
      double complex c12mh,c22mh
      double complex D0fin
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
     
     
      sbar = s+ieps
      tbar = t+ieps
      q3sbar = q3s+ieps
      q4sbar = q4s+ieps
      
      f2mh = (s+t-q4s)/(s*t)
      
      fpm = fpi*mus
c
c for comparison:
c	fpm = mus!-sbar 
                
      	c12mh = 2*li2(1.-q4sbar/sbar)-2*li2(1.-q3sbar/tbar) 
     &	       +2*li2(1.-sbar*f2mh)+2*li2(1.-tbar*f2mh) 
     &	       -2*li2(1.-q4sbar*f2mh)
     
     	c22mh = -pi2o6+0.5*(log((-sbar)/fpm))**2
     &		     +(log((-tbar)/fpm))**2
     &		 -0.5*(log((-q4sbar)/fpm))**2
     &		 -0.5*(log((-q3sbar)/fpm))**2
     &		 +0.5*(log((-q4sbar)/(-sbar)))**2
     &		 -0.5*(log((-q3sbar)/(-q4sbar)))**2
     &		 +0.5*(log((-q3sbar)/(-sbar)))**2

	D02mh = (c12mh+c22mh)/s/t
		            
      return
      end      
      
      
      
c---------  scalar box function: massless internal lines  --------------
c
      complex*16 function D02me(s,t,q2s,q4s,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c	     
c	     mi1s=mi2s=mi3s=mi4s=0 and q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	case b) two opposite massive external legs ("2me")
c		-> q1s=0,q3s=0,q2s.ne.0,q4s.ne.0
c
c		 g_2me = 1/(s*t-q2s*q4s)
c		 f_2me = (s+t-q2s-q4s)/(s*t-q2s*q4s)
c
c		 A_2me = 0
c		 B_2me =  2*ln[(-q2s-ieps)/(-t-ieps)]
c			 +2*ln[(-q4s-ieps)/(-s-ieps)]
c
c		 C1_2me = 2*Li2[1-s*f_2me]+2*Li2[1-t*f_2me] 
c		         -2*Li2[1-q2s*f_2me]-2*Li2[1-q4s*f_2me] 
c
c		 C2_2me = +ln^2[-s/(4*pi*mus)]
c			  +ln^2[-t/(4*pi*mus)]
c	 		  -ln^2[(-q4s)/(4*pi*mus)]
c			  -ln^2[(-q2s)/(4*pi*mus)]
c			  
c
c		(disregarded ieps here for simplicity)	
c
c	taken from Duplancic, Nizic;
c	(see my notes master-integrals.tex)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q2s,q4s,mus
      double precision g2me,f2me,fpm
      double complex sbar, tbar, q2sbar, q4sbar
      double complex c12me,c22me
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
          
      sbar = s+ieps
      tbar = t+ieps
      q2sbar = q2s+ieps
      q4sbar = q4s+ieps
      
      g2me = 1/(s*t-q2s*q4s)
      f2me = (s+t-q2s-q4s)/(s*t-q2s*q4s)
      
      fpm = fpi*mus

                
      	c12me = 2*li2(1.-sbar*f2me)+2*li2(1.-tbar*f2me) 
     &	       -2*li2(1.-q2sbar*f2me)-2*li2(1.-q4sbar*f2me)
     
     	c22me =   (log((-sbar)/fpm))**2+(log((-tbar)/fpm))**2
     &		 -(log((-q4sbar)/fpm))**2-(log((-q2sbar)/fpm))**2
     
	D02me = (c12me+c22me)*g2me
	
		            
      return
      end      
          
      
c---------  scalar box function: massless internal lines  --------------
c
      complex*16 function D01m(s,t,q4s,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c	     
c	     mi1s=mi2s=mi3s=mi4s=0 and q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	case c) one massive external legs ("1m")
c		-> q1s=0,q2s=0,q3s=0,q4s.ne.0
c
c		 g_1m = 1/(s*t)
c		 f_1m = (s+t-q4s)/(s*t)
c
c		 A_1m = 2
c		 B_1m = -2*ln[(-s-ieps)/(4*pi*mus)]
c			-2*ln[(-t-ieps)/(4*pi*mus)]
c			+2*ln[(-q4s-ieps)/(4*pi*mus)]
c
c		 C1_1m = 2*Li2[1-s*f_1m]+2*Li2[1-t*f_1m] 
c		        -2*Li2[1-q4s*f_1m] 
c
c		 C2_1m = -2*pi^2/3
c			 +ln^2[-s/(4*pi*mus)]
c			 +ln^2[-t/(4*pi*mus)]
c	 		 -ln^2[(-q4s)/(4*pi*mus)]
c			  
c
c		(disregarded ieps here for simplicity)	
c
c	taken from Duplancic, Nizic;
c	(see my notes master-integrals.tex)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q4s,mus
      double precision g1m,f1m,fpm
      double complex sbar, tbar, q4sbar
      double complex c11m,c21m
      double complex D0fin
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2
      parameter (ieps=(0d0,1d-16))
      external li2
      double precision pi,pi2o3,fpi,pia
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o3=pi**2/3d0,fpi=4d0*pi)
          
      sbar = s+ieps
      tbar = t+ieps
      q4sbar = q4s+ieps
      
      f1m = (s+t-q4s)/(s*t)
      
      fpm = fpi*mus

      	c11m = 2*li2(1.-sbar*f1m)+2*li2(1.-tbar*f1m) 
     &	      -2*li2(1.-q4sbar*f1m)
     
     	c21m =   +(log((-sbar)/fpm))**2+(log((-tbar)/fpm))**2
     &		 -(log((-q4sbar)/fpm))**2-2*pi2o3
     
	D01m = (c11m+c21m)/s/t
		            
      return
      end      
      
      
      
      
c -----------------------------------------------------------------------            
c ----------  scalar box function: massive internal lines  --------------
c
      complex*16 function D1i1e(s,t,q3s,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and
c	     
c	     mi1s=mi2s=mi3s=mi4s=0 and q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	result taken from Binoth et al
c	
c	Binoth's notation: I4d(s1,0,0,0;0,mis,0,0)
c	
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q3s,mis,mus
      double precision s1
      double complex sbar
      double precision x0,x1
      double complex x0bar,x1bar,x2bar
      double precision b1,b2,b3,b4,bb
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2,rfunc
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
      double complex c0i1e,c1i2e,c1d2e,c1i1e,D1i1e6
      external c0i1e,c1i2e,c1d2e,c1i1e
      
      ! change to notation of Binoth:
      s1 = q3s
      
      sbar = s+ieps
      
      x0 = s/(s+t-s1)
      x0bar = x0+ieps/(t-mis)
      
      x1 = s/(s-s1)
      x1bar = x1-ieps
      
      x2bar = sbar/(s-s1+mis)
     
      D1i1e6 =   (-t+mis)/t/(s1-s-t)
     & 		*(rfunc(dcmplx(x0),x2bar)+li2(1.-1/x0bar)-pi2o6)
     &		- mis/t/(s1-s)
     &		*(rfunc(dcmplx(x1),x2bar)+li2(1.-1/x1bar)-pi2o6)

      b1 = -1/(mis-t)      
      b2 = (s1-t)/s/(mis-t)      
      b3 = (t*(s-s1)+mis*(s+2.d0*t-s1))/s/(mis-t)**2
      b4 = -t/s/(mis-t)      
      bb = 2.d0*t*(s+t-s1)/s/(mis-t)**2
      
      D1i1e = b1*c0i1e(s,mus) + b2*c1i2e(s1,t,mis,mus)
     &	     +b3*c1d2e(s,s1,mis,mus) + b4*c1i1e(t,mis,mus)
     &	     +bb*D1i1e6 

      return
      end      
      
      
c ----------  scalar box function: massive internal lines  --------------
c
      complex*16 function D2i2e(s,t,q2s,q3s,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps, two equal internal masses mis and
c	     
c	     q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	A = 0
c	B = -log((-tb+mis)/mis)+log((-s1+mis)/(-sb+mis))
c
c
c	result taken from Binoth et al, TYPO CORRECTED IN (A.19)
c	
c	Binoth's notation: I4d(s1,s2,0,0;sb,tb;0,mis,mis,0)
c			=  I4d(q2s,q3s,0,0;t,s;0,mis,mis,0)	
c	
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q2s,q3s,mis,mus
      double precision s1,s2,sb,tb
      double complex s1bar,s2bar,sbar,tbar
      double precision x0,x1,x2
      double complex xp,xm,x1bar,x2bar
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2,rfunc
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc
      double precision pi,pi2o6,fpi
      parameter (pi=3.14159 26535 89793d0)
      parameter (pi2o6=pi**2/6d0,fpi=4d0*pi)
      
      double precision fpm,sqrm
      double complex ilog2,ilogp,irfun,ifin
      double precision gki
      
      ! change to notation of Binoth:
      s1 = q2s
      s2 = q3s
      sb = t
      tb = s
      
      s1bar = s1+ieps
      s2bar = s2+ieps
      tbar = tb+ieps
      sbar = sb+ieps
      
      xp = 0.5d0*(1d0+sqrt(1.-4.d0*mis/s2bar))
      xm = 0.5d0*(1d0-sqrt(1.-4.d0*mis/s2bar))
      
      x0 = sb/(sb+tb-s1)
      
      x1 = mis/tb
      x1bar = (mis-ieps)/tb
      
      x2 = (sb-mis)/(sb-s1)
      x2bar = (sbar-mis)/(sb-s1)
          
      fpm = fpi*mus
      sqrm = sqrt(fpm)
     
      
      ilog2 = (log((-tbar+mis)/sqrm))**2-(log(mis/sqrm))**2
     & 	     -(log((-s1bar+mis)/sqrm))**2+(log((-sbar+mis)/sqrm))**2
     
      irfun = -2.d0*rfunc(dcmplx(x0),x1bar)+2.d0*rfunc(dcmplx(x0),x2bar)
     &        +rfunc(dcmplx(1.-x1),xm)-rfunc(dcmplx(x1),xm)
     &        -rfunc(dcmplx(1.-x2),xm)+rfunc(dcmplx(x2),xm)
     
      ilogp = -(log(-s2bar)+log(xp-x1)+log(x1-xm))*log((1.-x1bar)/(-x1bar))
     &	      +(log(-s2bar)+log(xp-x2)+log(x2-xm))*log((1.-x2bar)/(-x2bar))
     
      ifin = ilog2+irfun+ilogp
      
      gki = sb*tb-(sb+tb-s1)*mis
      
      
      D2i2e = ifin/gki      
     

      return
      end      
   
      
c ----------  scalar box function: massive internal lines  --------------
c
      complex*16 function D2i2d(s,t,q1s,q3s,mis,mus) 
c
c notation for scalar integrals is as:
c	
c	D0(q1,q2,q3,mi1s,mi2s,mi3s,mi4s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			       	     1
c      -------------------------------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s][(q+q1+q2+q3)^2-mi4s]	
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and two equal internal masses mis
c		
c	 and q4 = -q1-q2-q3
c
c	     s = (q1+q2)^2,   t = (q2+q3)^2
c
c
c	A = 0
c	B = 0
c
c
c	result taken from Binoth et al
c	
c	Binoth's notation: I4d(s1=q1s,0,s3=q3s,0;s,t;0,0,mis,mis)
c	
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q1s,q3s,mis,mus
      double precision s1,s3
      double complex s1bar,s3bar,sbar,tbar
      double complex xp,xm,x0p,x0m
      double complex detS, J
        
      double precision zero
      parameter (zero=0d0)
      complex*16 ieps,li2,rfunc
      parameter (ieps=(0d0,1d-16))
      external li2,rfunc
       
      double complex ip,im,gki
      
      
      ! change to notation of Binoth:
      s1 = q1s
      s3 = q3s
      
      s1bar = s1+ieps
      s3bar = s3+ieps
      tbar = t+ieps
      sbar = s+ieps

      
      detS = (s*t-s1*s3+mis*(s-t))**2+4.d0*mis*(s*t-s1*s3)*(-s+s1+mis)      
      J = 2*(s3-t)*(s*t-s1*s3+mis*(s-t))+4.d0*(s*t-s1*s3)*(-s+s1+mis)
      
      xp = ((s*t-s1*s3)+mis*(s-t)+sqrt(detS-J*ieps))/2.d0/(s*t-s1*s3)
      xm = ((s*t-s1*s3)+mis*(s-t)-sqrt(detS-J*ieps))/2.d0/(s*t-s1*s3)
      
      x0p = 0.5d0*(1.+sqrt(1.-4.d0*mis/s3bar))
      x0m = 0.5d0*(1.-sqrt(1.-4.d0*mis/s3bar))
          
      ip = rfunc(xp,mis/tbar)-rfunc(1.-xp,mis/sbar)
     &    +rfunc(1.-xp,x0m)-rfunc(xp,x0m)
     &    +log((1.-xp)/(-xp))*
     &	        (log(sbar/s1bar)+log(tbar/s3bar)
     &		+log(xp-mis/tbar)
     &		-log(x0p-xp)-log(xp-x0m)+log(1.-xp-mis/sbar))
     
      im = rfunc(xm,mis/tbar)-rfunc(1.-xm,mis/sbar)
     &    +rfunc(1.-xm,x0m)-rfunc(xm,x0m)
     &    +log((1.-xm)/(-xm))*
     &	        (log(sbar/s1bar)+log(tbar/s3bar)
     &		+log(xm-mis/tbar)
     &		-log(x0p-xm)-log(xm-x0m)+log(1.-xm-mis/sbar))
      
      gki = sqrt(detS-J*ieps)
      
      D2i2d = (im-ip)/gki
      
     

      return
      end      
   
