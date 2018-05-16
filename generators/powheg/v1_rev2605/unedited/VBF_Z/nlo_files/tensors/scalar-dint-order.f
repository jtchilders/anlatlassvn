c	in this routine suitable box and triangle functions are selected,
c	depending on internal and external masses
c
c	functions implemented are:
c       complex*16 function D0fin(s,t,q1s,q2s,q3s,q4s,mus) 
c       complex*16 function C0fin(q1s,q2s,q3s,mus): arbitrary ordering of args. 
c       complex*16 function C0mfin(q1s,q2s,q3s,m1s,m2s,m3s,mus) 
c       complex*16 function B0fin(qs,mus) 
c
c
c      
c---------  scalar box function: massless internal lines  --------------
c
      complex*16 function D0fin(s,t,q1s,q2s,q3s,q4s,mus) 
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
c	      =  Gamma(1+eps) (  ----- - ----- + D0fin ) +O(eps)
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
c
c	case b) one massive external legs ("1m")
c
c
c	here only selected configurations are considered:
c		q1s = 0, q4s = 0
c		q2s,q3s can be 0 or finite
c		at least one mass must be non-zero
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision s,t,q1s,q2s,q3s,q4s,mus
      double precision zero
      parameter (zero=0d0)
      
      double complex D02mh,D02me,D01m
      external D02mh,D02me,D01m
      
      logical debug
      parameter (debug = .false.)
c      parameter (debug = .true.)

      integer icount

c      double precision small
c      parameter (small = 1d-9) ! define as 'zero'
       include 'dint_param.inc'
     
c ------------------------------------------------------------------------
	
	if (abs(q1s).lt.small) q1s = 0d0  	
	if (abs(q2s).lt.small) q2s = 0d0  	
	if (abs(q3s).lt.small) q3s = 0d0   	
	if (abs(q4s).lt.small) q4s = 0d0 

      icount = 0
     
      if (mus.le.0d0) then
         print*,' WARNING: non-positive mus in D0fin: ',mus
         print*,' Set D0fin = 0 '
         D0fin= 0
         return 
      endif
      
      ! only the cases with one or two non-zero external masses are considered;
      ! q1s and q4s must be zero:      
      if ((q1s).ne.zero) then ! require q1s = 0
	  print*,'q1s must be zero'
          icount = icount+1
      elseif ((q2s*q3s*q4s).ne.zero) then ! no D04m or D03m
	  print*,'D03m is not implemented'
          icount = icount+1
      elseif ((abs(q2s)+abs(q3s)+abs(q4s)).eq.zero) then ! no D00m
	  print*,'D00m is not implemented'
          icount = icount+1
      else
         if (q4s.eq.zero) then
	    if ((q2s*q3s).ne.zero) then !D2mh
	       D0fin = D02mh(t,s,q2s,q3s,mus)
	       if (debug) print*,'d02mh called (a) '
	    elseif ((q2s).eq.zero) then !D1mh
	       D0fin = D01m(t,s,q3s,mus)
	       if (debug) print*,'d01m called with q3s'
	    elseif ((q3s).eq.zero) then  !D!mh
	       D0fin = D01m(s,t,q2s,mus)
	       if (debug) print*,'d01m called with q2s'
	    else
               icount = icount+1
	    endif 
         endif  
      endif  
        
      if (icount.ne.0) then
	    print*,'D0 is not defined for this configuration:'
	    print*,'q1s = ',q1s, ' q2s = ',q2s
	    print*,'q3s = ',q3s, ' q4s = ',q4s
            print*,' -> set D0fin = 0 '
	
	    D0fin = 0
	    return
      endif	 
       
      return
      end            
c      
c------ scalar triangle function: massless propagators -------------------------
c
      complex*16 function C0fin(q1s,q2s,q3s,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			   1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
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
c	     mi1s=mi2s=mi3s=0 and q3 = -q1-q2
c
c
c	(see my notes master-integrals.tex)
c
c ------------------------------------------------------------------------
      implicit none
  
      double precision q1s,q2s,q3s,mus
      double precision zero
      parameter (zero=0d0)
      
      double complex C0i1e,C0i2e,C0i3e
      external C0i1e,C0i2e,C0i3e
      
      logical debug
      parameter (debug = .false.)
c      parameter (debug = .true.)
c      
c      double precision small
c      parameter (small = 1d-9) ! define as 'zero'
      include 'cint_param.inc'
     
c ------------------------------------------------------------------------
	
      if (abs(q1s).lt.small) q1s = 0d0        
      if (abs(q2s).lt.small) q2s = 0d0        
      if (abs(q3s).lt.small) q3s = 0d0        
     
      if (mus.le.0d0) then
         print*,' WARNING: non-positive mus in C0fin_exact: ',mus
         print*,' Set C0fin_exact = 0 '
         C0fin= 0
         return 
      endif
      if ((q1s.eq.0d0).and.(q2s.eq.0d0).and.(q3s.eq.0d0)) then
         print*,' WARNING: at least one argument of C0 must be non-zero '
         print*,' Set C0fin_exact = 0 '
         C0fin= 0
         return 
      endif      
      
      if ((q1s*q2s*q3s).ne.zero) then
	    C0fin = C0i3e(q1s,q2s,q3s,mus)
	    if (debug) print*,'c0i3e called'
      elseif (q3s.eq.zero) then
         if (q2s.eq.zero) then
	    C0fin = C0i1e(q1s,mus)
	    if (debug) print*,'c0i1e called'	
         elseif (q1s.eq.zero) then
	    C0fin = C0i1e(q2s,mus)
	    if (debug) print*,'c0i1e called'	
	 else
	    C0fin = C0i2e(q1s,q2s,mus)	
	    if (debug) print*,'c0i2e called (a)'	
	 endif        	 		
      elseif (q2s.eq.zero) then
         if (q1s.eq.zero) then
	    C0fin = C0i1e(q3s,mus)
	    if (debug) print*,'c0i1e called'	
	 else
	    C0fin = C0i2e(q3s,q1s,mus)
	    if (debug) print*,'c0i2e called (b)'	
	 endif        	 	
      else
	 C0fin = C0i2e(q2s,q3s,mus) 
	 if (debug) print*,'c0i2e called (c)'     	 
      endif 
      
      
      return
      end      
      
c     
c--------- scalar triangle function with internal masses ------------------
c
      complex*16 function C0mfin(q1s,q2s,q3s,m1s,m2s,m3s,mus) 
c
c notation for scalar integrals is as:
c	
c	C0(q1,q2,mi1s,mi2s,mi3s) = (2*pi*mu)^(4-d)/(i*pi^2)
c
c	Int d^dq
c			   1
c      ----------------------------------------------
c	[q^2-mi1s][(q+q1)^2-mi2s][(q+q1+q2)^2-mi3s]
c
c
c
c	  	         	  A_k     B_k
c	= Gamma(1+eps) * g_k* (  ----- - ----- + C1_k + C2_k ) +O(eps)
c	          		 eps^2    eps
c
c
c	with d = 4-2 eps and q3 = -(q1+q2)
c
c
c	(see my notes master-integrals.tex)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision q1s,q2s,q3s,m1s,m2s,m3s,mus
      double precision zero
      parameter (zero=0d0)
 	      
          
      double complex C1i1e,C1i2e,C1d2e,C2i2e,C2i3e
      external C1i1e,C1i2e,C1d2e,C2i2e,C2i3e
      
      logical debug
      parameter (debug = .false.)
      
     
      if (mus.le.0d0) then
         print*,' WARNING: non-positive mus in C0m: ',mus
         print*,' Set C0m = 0 '
         C0mfin= 0.d0
         return 
      endif
      if (q1s.eq.0d0) then
         print*,' WARNING: first argument of C0m must be non-zero '
         print*,' Set C0m = 0 '
         C0mfin= 0.d0
         return 
      endif
      if ((m1s+m2s+m3s).eq.0d0) then
         print*,' WARNING: at least one mass argument of C0m must be non-zero '
	 print*,' call C0 (function for massless propagators) instead'
         print*,' Set C0m = 0 '
         C0mfin= 0.d0
         return 
      endif
      if ((m1s+m2s+m3s).ne.(abs(m1s)+abs(m2s)+abs(m3s))) then
         print*,' WARNING: negative mass in C0m is not allowed '
	 print*,' Set C0m = 0 '
         C0mfin= 0.d0
         return 
      endif
      
      if (q3s.eq.zero) then
         if (q2s.eq.zero) then
	    if ((m1s+m3s).eq.zero) then 
	       C0mfin = C1i1e(q1s,m2s,mus)! m1s attached to massive external leg
	       if (debug) print*,'c1i1e called'	
	    else
	       goto 123
	    endif
	 else !q2s.ne.zero
	    if ((m1s+m3s).eq.zero) then
               C0mfin = C1i2e(q1s,q2s,m2s,mus) !m1s in between q1s and q2s
	       if (debug) print*,'c1i2e called'	
	    elseif ((m1s+m2s).eq.zero) then  
               C0mfin = C1d2e(q1s,q2s,m3s,mus) !m3s not between q1s and q2s
	       if (debug) print*,'c1d2e called'	
	    elseif (((m3s*m2s).ne.zero).and.(m1s.eq.zero)) then  
               C0mfin = C2i2e(q1s,q2s,m2s,m3s,mus)
	       if (debug) print*,'c2i2e called'	
	    else
	       goto 123
	    endif	       
	 endif	
      elseif ((q2s*q1s).ne.zero) then
         if((m1s.eq.zero).and.((m2s*m3s).ne.zero)) then     	     
            C0mfin = C2i3e(q1s,q2s,q3s,m2s,m3s,mus)
	    if (debug) print*,'c2i3e called'	
	 else
	    goto 123
	 endif  	    
      else 	
 123	 print*,'C0m is not defined for this configuration:'
	 print*,'q1s = ',q1s, ' q2s = ',q2s,'q3s = ',q3s
  	 print*,'m1s = ',m1s, ' m2s = ',m2s,'m3s = ',m3s
         print*,' -> set C0m = 0 '    
      	 C0mfin = 0d0	    
      endif 
     	    
      

      
      return
      end      
      
c------ scalar bubble function: massless propagators --------------------
c
      complex*16 function B0fin(qs,mus) 
c
c notation for scalar integrals is as:
c	
c	B0(q1,mi1s,mi2s) = (2*pi*mu)^(4-d)/(i*pi^2)
c	
c				    1
c   	     Int d^dq   ---------------------------
c			[q^2-mi1s][(q+q1)^2-mi2s]
c
c
c
c	  	         	  B_k
c	= Gamma(1+eps) * g_k* (  ----- + C1_k  ) +O(eps)
c	          		  eps
c
c
c	with d = 4-2 eps and
c	     
c	     mi1s=mi2s=0, write q1s = qs
c
c
c	(see my notes master-integrals.tex )
c
c	result for B0 is taken from Beenakker, PhD thesis (G.5)
c
c ------------------------------------------------------------------------
      implicit none
      
      double precision qs,mus
      double complex qsbar
      
      complex*16 ieps
      parameter (ieps=(0d0,1d-16))
      double precision pi,fpi,fpm
      parameter (pi=3.14159 26535 89793d0)
      parameter (fpi=4d0*pi)
     
      logical debug
      parameter (debug = .false.)
      
c      double precision small
c      parameter (small = 1d-9) ! define as 'zero'
      include 'int_param.inc'
c
c ------------------------------------------------------------------------
	
      if (abs(qs).lt.small) qs = 0d0  	

      if (mus.le.0d0) then
         print*,' WARNING: non-positive mus in B0fin: ',mus
         print*,' Set B0fin = 0 '
         B0fin= 0
         return 
      endif
      if (qs.eq.0d0) then
c         print*,' WARNING: momentum argument of B0 must be non-zero '
c         print*,' Set B0fin = 0 '
         B0fin= 0
         return 
      endif
      
      qsbar = qs+ieps
      fpm = fpi*mus
     
      B0fin = log(fpm/(-qsbar))+2.d0
       
      return
      end      
      
      
