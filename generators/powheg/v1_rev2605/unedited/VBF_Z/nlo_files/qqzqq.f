c
c returns Born (nlo=0) or 
c finite part of virtual contributions (NOT sum of Born + virtuals!) 
c
c---------------------------------------------------------------------------
c      
      subroutine qqzqq(pbar,sign, nlo, L,k,ans)
      implicit none
c
c	Last modified for POWHEG by Barbara Jaeger: 2012 June
C
C  QQZQQ calculates the matrix elements**2 for electroweak
c  lepton pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 e+ e-
C
C  Crossing related processes can be computed as well. 
c  Pauli interference terms for identical fermions are neglected. 
c  In particular, only the t-channel exchange of elctroweak bosons 
c  is considered. s-channl production of 2 weak bosons is NOT implemented.
c
c  summation of lepton helicities is performed in this file
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return resk = 2Re(M_Born^* M_virt)
c          NLO = 0       return resk = |M_born|^2  
c   for steering the calculation of finite box contributions
c   the following additional options are implemented
c          NLO = +4      set all finite box contributions to 0
c                        but include cvirt*M_Born contribution
c          NLO = -4      finite Born-type virtual contributions only
C
      include '../../include/pwhg_math.h'
      include '../../include/pwhg_st.h'

      include 'global.inc'
      include 'tensorl-hel.inc'

c electroweak couplings:
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c alfas, scales etc
      include 'scales.inc'
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+nv), musq
      double precision res(6),resv(6)
      real*8 ans
      double precision  p(0:3,4+nv), p21(0:4), p43(0:4)
      integer  sign(4+nv), nlo, mu, i, j,  k, 
     1         isig, isig1, isig3,h
      integer  ifl(4,6), js1, js3, L, is1, is3, is,kl
      
      double complex matv(6,-1:1,-1:1,-1:1,3)
      double complex mat(6,-1:1,-1:1,-1:1,3)
      double complex mm(6,-1:1,-1:1,-1:1)
      double complex mv12(6,-1:1,-1:1,-1:1), mv34(6,-1:1,-1:1,-1:1)

      double complex prop21(4), prop43(4)    
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2)
      
      double complex ma5,ma6
      double complex maa,maz,mza,mzz
      	
      double complex propt1(-1:1,-1:1,4,2), propt2(-1:1,-1:1,4,2),
     1		     propt(-1:1,-1:1,5:6,2)
     
      double complex contract_Tjj, dotcc, dotrc, dotcr, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotcr,dotqj, s1c
      logical linit
      
      data linit /.true./

      logical v_only, b_only
      parameter (v_only=.false.,b_only=.false.)
c      parameter (v_only=.true.,b_only=.false.)

      double complex m1,m2,m3,m4,maux
      double complex my1,my2,my3,my4
            
      save ifl,  linit
      double complex  zero
      parameter (zero = (0d0,0d0) )
      integer ii,ll
c
c for SMEs:
       integer nu,ka
c   
       double complex gmnc(0:3,0:3)
       double precision gmnr(0:3,0:3)
	
       double complex j21_a(0:3,-1:1),j43_a(0:3,-1:1)     
       double complex j21_abc(0:3,0:3,0:3,-1:1),j43_abc(0:3,0:3,0:3,-1:1)     
 
       double precision dotrr
       double complex s1r, sc3, contract_Tjjj
       external s1r, sc3, contract_Tjjj,dotrr
     
       double complex p2ze(-1:1),p2ae(-1:1)
       double complex p4ze(-1:1),p4ae(-1:1)
       double complex j21a_jqq(-1:1,-1:1)
       double complex j21a_ae(-1:1,-1:1),j21a_ze(-1:1,-1:1)
       double complex j21abc_zqj(-1:1,-1:1,-1:1),j21abc_aqj(-1:1,-1:1,-1:1)
       double complex j21abc_jqz(-1:1,-1:1,-1:1),j21abc_jqa(-1:1,-1:1,-1:1)
       double complex j43a_jqq(-1:1,-1:1)
       double complex j43a_ae(-1:1,-1:1),j43a_ze(-1:1,-1:1)
       double complex j43abc_zqj(-1:1,-1:1,-1:1),j43abc_aqj(-1:1,-1:1,-1:1)
       double complex j43abc_jqz(-1:1,-1:1,-1:1),j43abc_jqa(-1:1,-1:1,-1:1)
       double complex mae1(-1:1,-1:1,-1:1),mze1(-1:1,-1:1,-1:1)
       double complex mae2(-1:1,-1:1,-1:1),mze2(-1:1,-1:1,-1:1)
       double complex mae3(-1:1,-1:1,-1:1),mze3(-1:1,-1:1,-1:1)
       double complex mae4(-1:1,-1:1,-1:1),mze4(-1:1,-1:1,-1:1)
       double complex me1(6,-1:1,-1:1,-1:1)
       double complex me2(6,-1:1,-1:1,-1:1)
       double complex me3(6,-1:1,-1:1,-1:1)
       double complex me4(6,-1:1,-1:1,-1:1)
       double precision p2e(0:3),p2esq
       double precision p1e(0:3),p1esq
       double precision p4e(0:3),p4esq
       double precision p3e(0:3),p3esq
       double complex p2jqq(-1:1),p4jqq(-1:1)
       double complex qec(0:3),p43c(0:3),p21c(0:3)

      integer hmin,hmax,hstep
      common /hval/hmin,hmax,hstep

c stuff for virtual corrections

      double complex cc1(0:3),cc2(0:3),cc3(0:3,0:3,0:3),cc4(0:3)
      
      double complex c1ze(-1:1),c1ae(-1:1),c2ze(-1:1),c2ae(-1:1)
      double complex c2jqq(-1:1),c1jqq(-1:1)
      double complex con1z(-1:1,-1:1,-1:1),con1a(-1:1,-1:1,-1:1)
      double complex con2z(-1:1,-1:1,-1:1),con2a(-1:1,-1:1,-1:1)
      double complex con3z(-1:1,-1:1,-1:1),con3a(-1:1,-1:1,-1:1)
      double complex con4z(-1:1,-1:1,-1:1),con4a(-1:1,-1:1,-1:1)
      double complex mev2(6,-1:1,-1:1,-1:1),mev1(6,-1:1,-1:1,-1:1)
      double complex mev4(6,-1:1,-1:1,-1:1),mev3(6,-1:1,-1:1,-1:1)
c
c variables for powheg:
      logical nc_type,cc_type
      double precision q2_up,q2_lo,rup,rlo,lrup,lrlo
      double precision cvirtl    

c identify "bad" points (low qsq):
      logical qdamp      
c
      double precision c2,c2o4pi 
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4  ! compute box-corr. if T

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	
c
c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         print*,'qsqAmin in Born/virt=',qsqamin
         print*,'damping factor of 1d-20 below qsqAmin '
         if (v_only) print*,' graphs with vertex topology only'
         if (b_only) print*,' graphs with box topology only'
         linit = .false.
         kl = 1                  ! uucc
         ifl(1,kl) = 3
         ifl(2,kl) = 3
         ifl(3,kl) = 3
         ifl(4,kl) = 3
         kl = 2                  ! uuss
         ifl(1,kl) = 3
         ifl(2,kl) = 3
         ifl(3,kl) = 4
         ifl(4,kl) = 4
         kl = 3                  ! ddcc
         ifl(1,kl) = 4
         ifl(2,kl) = 4
         ifl(3,kl) = 3
         ifl(4,kl) = 3
         kl = 4                  ! ddss
         ifl(1,kl) = 4
         ifl(2,kl) = 4
         ifl(3,kl) = 4
         ifl(4,kl) = 4
         kl = 5                  ! udsc
         ifl(1,kl) = 3
         ifl(2,kl) = 4
         ifl(3,kl) = 4
         ifl(4,kl) = 3
         kl =6                   ! ducs
         ifl(1,kl) = 4
         ifl(2,kl) = 3
         ifl(3,kl) = 3
         ifl(4,kl) = 4
      endif

      if (k.le.4) then
         nc_type = .true.
         cc_type = .false.
      else
         cc_type = .true.
         nc_type = .false.
      endif

      do kl = 1,6 ! flavor combination
         do isig1 = -1,1,2
            do isig3 = -1,1,2
            do h = -1,1,2
               do i = 1,3 ! topology
                  mat(kl,isig1,isig3,h,i) = 0d0
                  matv(kl,isig1,isig3,h,i) = 0d0
               enddo
	    enddo   
            enddo
         enddo
      enddo
c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo	 
         p21(mu) =   p(mu,2) - p(mu,1)
         p43(mu) =   p(mu,4) - p(mu,3)	 
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2    

      qdamp = .false.
cc eliminate contributions from too low qsq in t-channel:
      if ( (abs(p21(4)).lt.qsqAmin).or.(abs(p43(4)).lt.qsqAmin)) then
      	qdamp = .true.
      endif

c get the vector boson propagator factors in complex mass scheme
c
      prop21(1) = 1/p21(4)
      prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
      prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
      prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      prop43(4) = prop43(3)

c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            do kl = 1,4
               propt1(isig1,isig3,kl,2) = 
     1            clr(ifl(1,kl),1,isig1)*clr(ifl(3,kl),1,isig3)*prop43(1) 
     2          + clr(ifl(1,kl),2,isig1)*clr(ifl(3,kl),2,isig3)*prop43(2)
               propt2(isig1,isig3,kl,2) = 
     1            clr(ifl(2,kl),1,isig1)*clr(ifl(3,kl),1,isig3)*prop43(1) 
     2          + clr(ifl(2,kl),2,isig1)*clr(ifl(3,kl),2,isig3)*prop43(2)	               
               propt1(isig1,isig3,kl,1) = 
     1            clr(ifl(1,kl),1,isig1)*clr(ifl(3,kl),1,isig3)*prop21(1)
     2          + clr(ifl(1,kl),2,isig1)*clr(ifl(3,kl),2,isig3)*prop21(2)
               propt2(isig1,isig3,kl,1) = 
     1            clr(ifl(1,kl),1,isig1)*clr(ifl(4,kl),1,isig3)*prop21(1)
     2          + clr(ifl(1,kl),2,isig1)*clr(ifl(4,kl),2,isig3)*prop21(2)
            enddo
         enddo
      enddo
c
      do kl = 5,6
     	 propt(-1,-1,kl,2) = clr(3,3,-1)**2*prop43(3)		     
     	 propt(-1,-1,kl,1) = clr(3,3,-1)**2*prop21(3)
      enddo   
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    
c compute MEs via "SMEs" :
	    
c construct Lorentz strings for SMEs:
c 	attention: some routines require complex/real arguments!!
c
c initialize:
	do isig = -1,1
	do mu = 0,3
	   j21_a(mu,isig) = 0d0
	   j43_a(mu,isig) = 0d0
	   do nu = 0,3
	   do ka = 0,3
	   	j21_abc(mu,nu,ka,isig) = 0d0
	   	j43_abc(mu,nu,ka,isig) = 0d0
	   enddo !ka
	   enddo !nu
	enddo !mu
	enddo !isig
	
        do mu = 0,3
	   do nu = 0,3
	      gmnc(mu,nu) = 0d0
	      gmnr(mu,nu) = 0d0
	   enddo
	enddo   
        gmnc(0,0) =  1d0 
        gmnc(1,1) = -1d0 
        gmnc(2,2) = -1d0
        gmnc(3,3) = -1d0
        
	gmnr(0,0) =  1d0
        gmnr(1,1) = -1d0
        gmnr(2,2) = -1d0
        gmnr(3,3) = -1d0
	
	do isig = -1,1,2
     
        ! ubar.gamma_mu.u
        do mu = 0,3
           j21_a(mu,isig) = 
     &  	s1r(psi(1,isig,2),gmnr(0,mu),(mu.eq.0),isig,psi(1,isig,1))
     
     	   j43_a(mu,isig) = 
     &  	s1r(psi(1,isig,4),gmnr(0,mu),(mu.eq.0),isig,psi(1,isig,3))     
        enddo 
	
	! rename:
	do mu = 0,3
	   jqq(mu,isig,1) =  j21_a(mu,isig)
	   jqq(mu,isig,2) =  j43_a(mu,isig)
        enddo
	
        ! ubar.gamma_mu.gamma_nu.gamma_ka.u
	do mu = 0,3
	   do nu = 0,3
	      do ka = 0,3
	   	 j21_abc(mu,nu,ka,isig) = sc3(psi(1,isig,2),
     &	   		gmnc(0,mu),gmnr(0,nu),gmnc(0,ka),psi(1,isig,1),isig)
     
	   	 j43_abc(mu,nu,ka,isig) = sc3(psi(1,isig,4),
     &	   		gmnc(0,mu),gmnr(0,nu),gmnc(0,ka),psi(1,isig,3),isig)
	      enddo !ka
	   enddo !nu
	enddo !mu	

	enddo !isig

        if (b_only) goto 444

c---------------------------------------------------------------
c
c for (V) topology:
c
       do h = hmin,hmax,hstep
c
C NC contributions:
       if (nc_type) then 
       do isig1 = -1,1,2
         do isig3 = -1,1,2
            maa = contract_Tjj(aaee(0,0,L,h),j21_a(0,isig1),j43_a(0,isig3))
            maz = contract_Tjj(azee(0,0,L,h),j21_a(0,isig1),j43_a(0,isig3))
            mza = contract_Tjj(zaee(0,0,L,h),j21_a(0,isig1),j43_a(0,isig3))
            mzz = contract_Tjj(zzee(0,0,L,h),j21_a(0,isig1),j43_a(0,isig3))

c            do k = 1,4
               mat(k,isig1,isig3,h,1) = 
     1              maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2            + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3            + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4            + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
 		
c            enddo         
      	enddo !isig3
      enddo !isig1
	    	
C CC contributions:
      else ! cc_type
c
      isig1 = -1
      isig3 = -1
        if (k.eq.5) then    
         ma5 = contract_Tjj(CCee(0,0,L,h),j21_a(0,isig1),j43_a(0,isig3))
         mat(5,isig1,isig3,h,1) = ma5*clr(3,3,-1)**2
        elseif (k.eq.6) then   			       
         ma6 = contract_Tjj(CCee6(0,0,L,h),j43_a(0,isig3),j21_a(0,isig1))
         mat(6,isig1,isig3,h,1) = ma6*clr(3,3,-1)**2
        endif !k
      
      endif !nc/cc

      enddo !h

 444  continue
      if (v_only) goto 999  
	
c---------------------------------------------------------------
c
c for (B) topology:

	do mu = 0,3
	   p2e(mu) =  p(mu,2)+qe(mu)
	   p1e(mu) =  p(mu,1)-qe(mu)
	   p4e(mu) =  p(mu,4)+qe(mu)
	   p3e(mu) =  p(mu,3)-qe(mu)
	enddo
	
	p2esq = dotrr(p2e(0),p2e(0)) 
	p1esq = dotrr(p1e(0),p1e(0)) 
	p4esq = dotrr(p4e(0),p4e(0)) 
	p3esq = dotrr(p3e(0),p3e(0)) 
	
	do h = hmin,hmax,hstep
	   p2ze(h) = dotrc(p(0,2),ze(1,h))
	   p2ae(h) = dotrc(p(0,2),ae(1,h))
	   p4ze(h) = dotrc(p(0,4),ze(1,h))
	   p4ae(h) = dotrc(p(0,4),ae(1,h))
	enddo

	do isig = -1,1,2
	   p2jqq(isig) = dotrc(p(0,2),jqq(0,isig,2))
	   p4jqq(isig) = dotrc(p(0,4),jqq(0,isig,1))
	enddo

        do mu = 0,3
           qec(mu)  = dcmplx(qe(mu))
           p43c(mu) = dcmplx(p43(mu))
           p21c(mu) = dcmplx(p21(mu))
        enddo !mu

	do isig1 = -1,1,2
	do isig3 = -1,1,2
	
cccccccccccccc

	! evaluate q1->V1 V2 q2 with V1 -> J43, V2 -> Z:
	
	j21a_jqq(isig1,isig3) = dotcc(j21_a(0,isig1),jqq(0,isig3,2))

	do h = hmin,hmax,hstep
	j21abc_zqj(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     %		ze(1,h),qec(0),jqq(0,isig3,2))		
	j21abc_aqj(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     %		ae(1,h),qec(0),jqq(0,isig3,2))	
          	
	
	mae2(isig1,isig3,h) = (2d0*j21a_jqq(isig1,isig3)*p2ae(h) + 
     %			      j21abc_aqj(isig1,isig3,h))
	
	mze2(isig1,isig3,h) = (2d0*j21a_jqq(isig1,isig3)*p2ze(h) + 
     %			      j21abc_zqj(isig1,isig3,h))
        enddo !h

ccccc
	
	! evaluate q3->V1 V2 q4 with V1 -> J21, V2 -> Z:
	
	j43a_jqq(isig3,isig1) = dotcc(j43_a(0,isig3),jqq(0,isig1,1))

	do h = hmin,hmax,hstep
	j43abc_zqj(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     %	        ze(1,h),qec(0),jqq(0,isig1,1))	      
	j43abc_aqj(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     %		ae(1,h),qec(0),jqq(0,isig1,1))	
          	
	
	mae4(isig3,isig1,h) = (2d0*j43a_jqq(isig3,isig1)*p4ae(h) + 
     %			      j43abc_aqj(isig3,isig1,h))
	
	mze4(isig3,isig1,h) = (2d0*j43a_jqq(isig3,isig1)*p4ze(h) + 
     %			      j43abc_zqj(isig3,isig1,h))
        enddo !h


        enddo !isig1
	enddo !isig3
	
cccccccccccccc

	! evaluate q1->V2 V1 q2 with V1 -> J43, V2 -> Z:
	
	do h = hmin,hmax,hstep
	do isig1 = -1,1,2
	
	j21a_ae(isig1,h) = dotcc(j21_a(0,isig1),ae(1,h))
	j21a_ze(isig1,h) = dotcc(j21_a(0,isig1),ze(1,h))	
	
	do isig3 = -1,1,2

	j21abc_jqz(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     %		jqq(0,isig3,2),p43c(0),ze(1,h))
	j21abc_jqa(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     %		jqq(0,isig3,2),p43c(0),ae(1,h))
          	
	
	mae1(isig1,isig3,h) = (2d0*j21a_ae(isig1,h)*p2jqq(isig3) + 
     %			      j21abc_jqa(isig1,isig3,h))
	
	mze1(isig1,isig3,h) = (2d0*j21a_ze(isig1,h)*p2jqq(isig3) + 
     %			      j21abc_jqz(isig1,isig3,h))

        enddo !isig3
	enddo !isig1

	
cccccccccccccc

	! evaluate q1->V2 V1 q2 with V1 -> J21, V2 -> Z:
	
	do isig3 = -1,1,2
	
	j43a_ae(isig3,h) = dotcc(j43_a(0,isig3),ae(1,h))
	j43a_ze(isig3,h) = dotcc(j43_a(0,isig3),ze(1,h))    
	
	do isig1 = -1,1,2

	j43abc_jqz(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     %	      jqq(0,isig1,1),p21c(0),ze(1,h))
	j43abc_jqa(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     %	      jqq(0,isig1,1),p21c(0),ae(1,h))
          	
	
	mae3(isig3,isig1,h) = (2d0*j43a_ae(isig3,h)*p4jqq(isig1) + 
     %			      j43abc_jqa(isig3,isig1,h))
	
	mze3(isig3,isig1,h) = (2d0*j43a_ze(isig3,h)*p4jqq(isig1) + 
     %			      j43abc_jqz(isig3,isig1,h))

        enddo !isig3
	enddo !isig1
	

        if (nc_type) then
c	do k = 1,4
	do isig1 = -1,1,2
 	do isig3 = -1,1,2
	
	me2(k,isig1,isig3,h) = (mae2(isig1,isig3,h)*clr(ifl(2,k),1,isig1)+    
     $			        mze2(isig1,isig3,h)*clr(ifl(2,k),2,isig1))*
     %			     propt1(isig1,isig3,k,2)/p2esq
	
	me1(k,isig1,isig3,h) = (mae1(isig1,isig3,h)*clr(ifl(1,k),1,isig1)+
     $		                mze1(isig1,isig3,h)*clr(ifl(1,k),2,isig1))*
     %			     propt2(isig1,isig3,k,2)/p1esq
	
	me4(k,isig1,isig3,h) = (mae4(isig3,isig1,h)*clr(ifl(4,k),1,isig3)+
     $		                mze4(isig3,isig1,h)*clr(ifl(4,k),2,isig3))*
     %			     propt1(isig1,isig3,k,1)/p4esq
	
	me3(k,isig1,isig3,h) = (mae3(isig3,isig1,h)*clr(ifl(3,k),1,isig3)+
     $		                mze3(isig3,isig1,h)*clr(ifl(3,k),2,isig3))*
     %			     propt2(isig1,isig3,k,1)/p3esq
     
     	mat(k,isig1,isig3,h,2) = me1(k,isig1,isig3,h)+me2(k,isig1,isig3,h)     
     	mat(k,isig1,isig3,h,3) = me3(k,isig1,isig3,h)+me4(k,isig1,isig3,h)

        enddo !isig3
	enddo !isig1
     
c     	enddo !k

        else ! cc_type

c	do k = 5,6
	isig1 = -1
 	isig3 = -1
	
	me2(k,isig1,isig3,h) =  (mae2(isig1,isig3,h)*clr(ifl(2,k),1,isig1)+    
     $			         mze2(isig1,isig3,h)*clr(ifl(2,k),2,isig1))*
     %			     propt(isig1,isig3,k,2)/p2esq
	
	me1(k,isig1,isig3,h) =(mae1(isig1,isig3,h)*clr(ifl(1,k),1,isig1)+
     $		               mze1(isig1,isig3,h)*clr(ifl(1,k),2,isig1))*
     %			     propt(isig1,isig3,k,2)/p1esq
	
	me4(k,isig1,isig3,h) = (mae4(isig3,isig1,h)*clr(ifl(4,k),1,isig3)+
     $		                mze4(isig3,isig1,h)*clr(ifl(4,k),2,isig3))*
     %			     propt(isig3,isig1,k,1)/p4esq
	
	me3(k,isig1,isig3,h) = (mae3(isig3,isig1,h)*clr(ifl(3,k),1,isig3)+
     $		                mze3(isig3,isig1,h)*clr(ifl(3,k),2,isig3))*
     %			     propt(isig3,isig1,k,1)/p3esq
     
     	mat(k,isig1,isig3,h,2) = me1(k,isig1,isig3,h)+me2(k,isig1,isig3,h)     
     	mat(k,isig1,isig3,h,3) = me3(k,isig1,isig3,h)+me4(k,isig1,isig3,h)

c       enddo !k

       endif !nc/cc
       
       enddo !h
     
c----------------------      

	if (lbox.and.(.not.qdamp)) then

c BOX CORRECTION TO UPPER FERMION LINE		
c
	! evaluate q1->V1 V2 q2 with V1 -> J43, V2 -> Z:

	    musq = -p21(4)		
	    call boxcorr(p(0,1),p43(0),qe(0),p(0,2),musq,cc1,cc2,cc3,cc4)
	    
	    do h = hmin,hmax,hstep
	    
	    ! for L1:
	    ! (L1.J43) = j21a_jqq(isig1,isig3) 
	    ! (c1.eps2):	    
	    c1ze(h) = dotcc(cc1(0),ze(1,h))
	    c1ae(h) = dotcc(cc1(0),ae(1,h))	
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con1z(isig1,isig3,h) = j21a_jqq(isig1,isig3)*c1ze(h)
	    	con1a(isig1,isig3,h) = j21a_jqq(isig1,isig3)*c1ae(h)
	    enddo
	    enddo
	   
	    ! for L2:
	    ! (L2.eps) = j21a_ae(isig1),j21a_ze(isig1)
	    ! (c2.J43):
	    do isig = -1,1,2
	   	c2jqq(isig) = dotcc(cc2(0),jqq(0,isig,2))
	    enddo 	   
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con2z(isig1,isig3,h) = j21a_ze(isig1,h)*c2jqq(isig3)
	    	con2a(isig1,isig3,h) = j21a_ae(isig1,h)*c2jqq(isig3)
	    enddo
	    enddo
	    
	    ! for L3:
	    ! need c3_si,rho,om . eps2^si * J43^rho * L3^om
	    ! L3^om = j21_a(om,isig1)
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con3z(isig1,isig3,h) = contract_Tjjj(cc3(0,0,0),
     &			   ze(1,h),jqq(0,isig3,2),j21_a(0,isig1))		
		con3a(isig1,isig3,h) = contract_Tjjj(cc3(0,0,0),
     &			   ae(1,h),jqq(0,isig3,2),j21_a(0,isig1))		
	    enddo
	    enddo
	    

	    ! for L4:
	    ! L4^si,om,rho = j21_abc(si,om,rho,isig1)
	    ! need L4_si,om,rho . eps2^si * c4^om * J43^rho
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con4z(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     &						ze(1,h),cc4(0),jqq(0,isig3,2))
     
     		con4a(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     &						ae(1,h),cc4(0),jqq(0,isig3,2))
	    enddo
	    enddo

	if (nc_type) then    
c	do k = 1,4
	do isig1 = -1,1,2
 	do isig3 = -1,1,2

	mev2(k,isig1,isig3,h) = (clr(ifl(2,k),1,isig1)*
     &		(con1a(isig1,isig3,h)+con2a(isig1,isig3,h)+
     &		 con3a(isig1,isig3,h)+con4a(isig1,isig3,h))+
     &		 	       clr(ifl(2,k),2,isig1)*
     $		(con1z(isig1,isig3,h)+con2z(isig1,isig3,h)+
     &		 con3z(isig1,isig3,h)+con4z(isig1,isig3,h)))*
     %			     propt1(isig1,isig3,k,2)
   
 
        enddo !isig3
	enddo !isig1
     
c     	enddo !k

        else !cc

c	do k = 5,6
	isig1 = -1
 	isig3 = -1
    
	mev2(k,isig1,isig3,h) = (clr(ifl(2,k),1,isig1)*
     &		(con1a(isig1,isig3,h)+con2a(isig1,isig3,h)+
     &		 con3a(isig1,isig3,h)+con4a(isig1,isig3,h))+
     &		 	       clr(ifl(2,k),2,isig1)*
     $		(con1z(isig1,isig3,h)+con2z(isig1,isig3,h)+
     &		 con3z(isig1,isig3,h)+con4z(isig1,isig3,h)))*
     %			     propt(isig1,isig3,k,2)
        
c	enddo !k

        endif !nc/cc
c	
	enddo !h

c----------------------      
c----------------------      

c now the same with crossed boxline routine:
	
	! evaluate q1->V2 V1 q2 with V1 -> J43, V2 -> Z:

	    musq = -p21(4)		
	    call boxcorr_cross(p(0,1),p43(0),qe(0),p(0,2),musq,cc1,cc2,cc3,cc4)
	    
	    do h = hmin,hmax,hstep
	    
	    
	    ! for L1:
	    ! (L1.J43) = j21a_jqq(isig1,isig3) 
	    ! (c1.eps2):	    
	    c1ze(h) = dotcc(cc1(0),ze(1,h))
	    c1ae(h) = dotcc(cc1(0),ae(1,h))	
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con1z(isig1,isig3,h) = j21a_jqq(isig1,isig3)*c1ze(h)
	    	con1a(isig1,isig3,h) = j21a_jqq(isig1,isig3)*c1ae(h)
	    enddo
	    enddo
	   
	    ! for L2:
	    ! (L2.eps) = j21a_ae(isig1),j21a_ze(isig1)
	    ! (c2.J43):
	    do isig = -1,1,2
	   	c2jqq(isig) = dotcc(cc2(0),jqq(0,isig,2))
	    enddo 	   
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con2z(isig1,isig3,h) = j21a_ze(isig1,h)*c2jqq(isig3)
	    	con2a(isig1,isig3,h) = j21a_ae(isig1,h)*c2jqq(isig3)
	    enddo
	    enddo
	    
	    ! for L3:
	    ! need c3_si,rho,om . eps2^si * J43^rho * L3^om
	    ! L3^om = j21_a(om,isig1)
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con3z(isig1,isig3,h) = contract_Tjjj(cc3(0,0,0),
     &			   ze(1,h),jqq(0,isig3,2),j21_a(0,isig1))		
		con3a(isig1,isig3,h) = contract_Tjjj(cc3(0,0,0),
     &			   ae(1,h),jqq(0,isig3,2),j21_a(0,isig1))		
	    enddo
	    enddo
	    

	    ! for L4:
	    ! L4^si,om,rho = j21_abc(si,om,rho,isig1)
	    ! need L4_si,om,rho . eps2^si * c4^om * J43^rho
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con4z(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     &						ze(1,h),cc4(0),jqq(0,isig3,2))
     
     		con4a(isig1,isig3,h) = contract_Tjjj(j21_abc(0,0,0,isig1),
     &						ae(1,h),cc4(0),jqq(0,isig3,2))
	    enddo
	    enddo	    


	if (nc_type) then    
c	do k = 1,4
	do isig1 = -1,1,2
 	do isig3 = -1,1,2

	mev1(k,isig1,isig3,h) = (clr(ifl(1,k),1,isig1)*
     &		(con1a(isig1,isig3,h)+con2a(isig1,isig3,h)+
     &		 con3a(isig1,isig3,h)+con4a(isig1,isig3,h))+
     &		 	       clr(ifl(1,k),2,isig1)*
     $		(con1z(isig1,isig3,h)+con2z(isig1,isig3,h)+
     &		 con3z(isig1,isig3,h)+con4z(isig1,isig3,h)))*
     %			     propt2(isig1,isig3,k,2)
   
     	matv(k,isig1,isig3,h,2) = mev1(k,isig1,isig3,h)+
     &                            mev2(k,isig1,isig3,h)    

        enddo !isig3
	enddo !isig1
     
c     	enddo !k

        else !cc

c	do k = 5,6
	isig1 = -1
 	isig3 = -1
    
	mev1(k,isig1,isig3,h) = (clr(ifl(1,k),1,isig1)*
     &		(con1a(isig1,isig3,h)+con2a(isig1,isig3,h)+
     &		 con3a(isig1,isig3,h)+con4a(isig1,isig3,h))+
     &		 	       clr(ifl(1,k),2,isig1)*
     $		(con1z(isig1,isig3,h)+con2z(isig1,isig3,h)+
     &		 con3z(isig1,isig3,h)+con4z(isig1,isig3,h)))*
     %			     propt(isig1,isig3,k,2)

     	matv(k,isig1,isig3,h,2) = mev1(k,isig1,isig3,h)+
     &                            mev2(k,isig1,isig3,h)  

c	enddo !k

        endif !nc/cc
	
	enddo !h

c----------------------      
c	
c BOX CORRECTION TO LOWER FERMION LINE		

	! evaluate q3->V1 V2 q4 with V1 -> J21, V2 -> Z:

	    musq = -p43(4)		
	    call boxcorr(p(0,3),p21(0),qe(0),p(0,4),musq,cc1,cc2,cc3,cc4)
	    
	    do h = hmin,hmax,hstep
	    
	    ! for L1:
	    ! (L1.J21) = j43a_jqq(isig3,isig1) 
	    ! (c1.eps2):	    
	    c1ze(h) = dotcc(cc1(0),ze(1,h))
	    c1ae(h) = dotcc(cc1(0),ae(1,h))	
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con1z(isig3,isig1,h) = j43a_jqq(isig3,isig1)*c1ze(h)
	    	con1a(isig3,isig1,h) = j43a_jqq(isig3,isig1)*c1ae(h)
	    enddo
	    enddo
	    
	    ! for L2:
	    ! (L2.eps) = j43a_ae(isig1),j43a_ze(isig1)
	    ! (c2.J21):
	    do isig = -1,1,2
	   	c2jqq(isig) = dotcc(cc2(0),jqq(0,isig,1))
	    enddo 	   
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con2z(isig3,isig1,h) = j43a_ze(isig3,h)*c2jqq(isig1)
	    	con2a(isig3,isig1,h) = j43a_ae(isig3,h)*c2jqq(isig1)
	    enddo
	    enddo
   
	    ! for L3:
	    ! need c3_si,rho,om . eps2^si * J21^rho * L3^om
	    ! L3^om = j43_a(om,isig1)
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con3z(isig3,isig1,h) = contract_Tjjj(cc3(0,0,0),
     &			   ze(1,h),jqq(0,isig1,1),j43_a(0,isig3))		
		con3a(isig3,isig1,h) = contract_Tjjj(cc3(0,0,0),
     &			   ae(1,h),jqq(0,isig1,1),j43_a(0,isig3))		
	    enddo
	    enddo

	    ! for L4:
	    ! L4^si,om,rho = j43_abc(si,om,rho,isig3)
	    ! need L4_si,om,rho . eps2^si * c4^om * J21^rho
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con4z(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     &						ze(1,h),cc4(0),jqq(0,isig1,1))
     
     		con4a(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     &						ae(1,h),cc4(0),jqq(0,isig1,1))
	    enddo
	    enddo

        if (nc_type) then    
c	do k = 1,4
	do isig1 = -1,1,2
 	do isig3 = -1,1,2

	mev4(k,isig1,isig3,h) = (clr(ifl(4,k),1,isig3)*
     &		(con1a(isig3,isig1,h)+con2a(isig3,isig1,h)+
     &		 con3a(isig3,isig1,h)+con4a(isig3,isig1,h))+
     &		 	       clr(ifl(4,k),2,isig3)*
     $		(con1z(isig3,isig1,h)+con2z(isig3,isig1,h)+
     &		 con3z(isig3,isig1,h)+con4z(isig3,isig1,h)))*
     %			     propt1(isig1,isig3,k,1)
          
 
        enddo !isig3
	enddo !isig1
     
c     	enddo !k

        else !cc

c	do k = 5,6
	isig1 = -1
 	isig3 = -1
    
	mev4(k,isig1,isig3,h) = (clr(ifl(4,k),1,isig3)*
     &		(con1a(isig3,isig1,h)+con2a(isig3,isig1,h)+
     &		 con3a(isig3,isig1,h)+con4a(isig3,isig1,h))+
     &		 	       clr(ifl(4,k),2,isig3)*
     $		(con1z(isig3,isig1,h)+con2z(isig3,isig1,h)+
     &		 con3z(isig3,isig1,h)+con4z(isig3,isig1,h)))*
     %			     propt(isig1,isig3,k,1)
        
  
c	enddo !k

        endif !nc/cc
	
	enddo !h

c=================================

c BOX CORRECTION TO LOWER FERMION LINE		

	! evaluate q3->V1 V2 q4 with V1 -> J21, V2 -> Z:

	    musq = -p43(4)		
	    call boxcorr_cross(p(0,3),p21(0),qe(0),p(0,4),musq,cc1,cc2,cc3,cc4)
	    
	    do h = hmin,hmax,hstep
	    
	    ! for L1:
	    ! (L1.J21) = j43a_jqq(isig3,isig1) 
	    ! (c1.eps2):	    
	    c1ze(h) = dotcc(cc1(0),ze(1,h))
	    c1ae(h) = dotcc(cc1(0),ae(1,h))	
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con1z(isig3,isig1,h) = j43a_jqq(isig3,isig1)*c1ze(h)
	    	con1a(isig3,isig1,h) = j43a_jqq(isig3,isig1)*c1ae(h)
	    enddo
	    enddo
	    
	    ! for L2:
	    ! (L2.eps) = j43a_ae(isig1),j43a_ze(isig1)
	    ! (c2.J21):
	    do isig = -1,1,2
	   	c2jqq(isig) = dotcc(cc2(0),jqq(0,isig,1))
	    enddo 	   
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
	    	con2z(isig3,isig1,h) = j43a_ze(isig3,h)*c2jqq(isig1)
	    	con2a(isig3,isig1,h) = j43a_ae(isig3,h)*c2jqq(isig1)
	    enddo
	    enddo
   
	    ! for L3:
	    ! need c3_si,rho,om . eps2^si * J21^rho * L3^om
	    ! L3^om = j43_a(om,isig1)
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con3z(isig3,isig1,h) = contract_Tjjj(cc3(0,0,0),
     &			   ze(1,h),jqq(0,isig1,1),j43_a(0,isig3))		
		con3a(isig3,isig1,h) = contract_Tjjj(cc3(0,0,0),
     &			   ae(1,h),jqq(0,isig1,1),j43_a(0,isig3))		
	    enddo
	    enddo

	    ! for L4:
	    ! L4^si,om,rho = j43_abc(si,om,rho,isig3)
	    ! need L4_si,om,rho . eps2^si * c4^om * J21^rho
	    do isig1 = -1,1,2
	    do isig3 = -1,1,2
		con4z(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     &						ze(1,h),cc4(0),jqq(0,isig1,1))
     
     		con4a(isig3,isig1,h) = contract_Tjjj(j43_abc(0,0,0,isig3),
     &						ae(1,h),cc4(0),jqq(0,isig1,1))
	    enddo
	    enddo
	    
        if (nc_type) then    
c	do k = 1,4
	do isig1 = -1,1,2
 	do isig3 = -1,1,2

	mev3(k,isig1,isig3,h) = (clr(ifl(3,k),1,isig3)*
     &		(con1a(isig3,isig1,h)+con2a(isig3,isig1,h)+
     &		 con3a(isig3,isig1,h)+con4a(isig3,isig1,h))+
     &		 	       clr(ifl(3,k),2,isig3)*
     $		(con1z(isig3,isig1,h)+con2z(isig3,isig1,h)+
     &		 con3z(isig3,isig1,h)+con4z(isig3,isig1,h)))*
     %			     propt2(isig1,isig3,k,1)
   
     	matv(k,isig1,isig3,h,3) = mev3(k,isig1,isig3,h)+
     &                            mev4(k,isig1,isig3,h)
    
 
        enddo !isig3
	enddo !isig1
     
c     	enddo !k

        else !cc

c	do k = 5,6
	isig1 = -1
 	isig3 = -1
    
	mev3(k,isig1,isig3,h) = (clr(ifl(3,k),1,isig3)*
     &		(con1a(isig3,isig1,h)+con2a(isig3,isig1,h)+
     &		 con3a(isig3,isig1,h)+con4a(isig3,isig1,h))+
     &		 	       clr(ifl(3,k),2,isig3)*
     $		(con1z(isig3,isig1,h)+con2z(isig3,isig1,h)+
     &		 con3z(isig3,isig1,h)+con4z(isig3,isig1,h)))*
     %			     propt(isig1,isig3,k,1)

     	matv(k,isig1,isig3,h,3) =   mev3(k,isig1,isig3,h)+
     &                              mev4(k,isig1,isig3,h) 

c	enddo !k

        endif !nc/cc
	
	
	enddo !h
	
c=================================
c
	endif !lbox
c
c --------------------------------------------------------------
     
c sum the graphs, square them and map them onto res

c i = 1		V
c i = 2,3       B

 999   continue
 
 	if (b_only) then
	  do kl = 1,6
	  do isig1 = -1,1,2
	  do isig3 = -1,1,2
	    do h = hmin,hmax,hstep
	      mat(kl,isig1,isig3,h,1) = 0d0
	      matv(kl,isig1,isig3,h,1) = 0d0
	    enddo !h
	  enddo !isig3
          enddo !isig1
	  enddo !kl   
	endif
	if (v_only) then
	  do kl = 1,6
	  do isig1 = -1,1,2
	  do isig3 = -1,1,2
	    do h = hmin,hmax,hstep
	      mat(kl,isig1,isig3,h,2:3) = 0d0
	      matv(kl,isig1,isig3,h,2:3) = 0d0
	    enddo !h
	  enddo !isig3
          enddo !isig1
	  enddo !kl   
	endif

c      do k = 1,6
	res(k) = 0
        resv(k) = 0
        do h = hmin,hmax,hstep
	do isig1 = -1,1,2
	    do isig3 = -1,1,2
	       mm(k,isig1,isig3,h) = 0
	       do i = 1,3
		  mm(k,isig1,isig3,h) = mm(k,isig1,isig3,h) + 
     1  				mat(k,isig1,isig3,h,i)
	       enddo !i
	      
	       res(k) = res(k) + dreal(mm(k,isig1,isig3,h))**2
     &  		       + dimag(mm(k,isig1,isig3,h))**2
		mv12(k,isig1,isig3,h) = 0d0
		mv34(k,isig1,isig3,h) = 0d0
ccc

		if (lnlo.and.(.not.qdamp)) then
		  mv12(k,isig1,isig3,h) = matv(k,isig1,isig3,h,2) 
		  mv34(k,isig1,isig3,h) = matv(k,isig1,isig3,h,3) 

c  add Born type term and multiply by F_q = alphas*C_2/4pi;
c  the factor pi^2/3+9/2 for the born term is 
c  after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
c
c  comply with POWHEG normalization:
c  vertex corrections in dreg:
c
c   V^nu = (4*Pi)^ep * Gamma(1+ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7)*B^nu
c
c         = (4*Pi)^ep / Gamma(1-ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7-Pi^2/3)*B^nu
c
c     The factor  (4*Pi)^ep/Gamma(1-ep) IS NOT RETURNED by this subroutine
c     and it's thought as factorized in front of the real counterterms too.
c
c     if (4*Pi)^ep / Gamma(1-ep) is collected in front then cvirt:
c     set parameter in dred (cvirtl = -8d0)
c
c     squared momentum of the weak boson connected with the upper line
          q2_up = p21(4)
c     squared momentum of the weak boson connected with the lower line
          q2_lo = p43(4)

          rup = st_muren2/(-q2_up)
          if (rup.lt.0d0.and..not.qdamp) then
              write(*,*) 'Error in setvirtual: q2_up should be < 0!!'
              write(*,*) 'q2_up = ', q2_up
              write(*,*) 'st_muren2 = ', st_muren2
              ! set to dummy value:
              	rup = 1d0
c              stop
          endif      
          rlo = st_muren2/(-q2_lo)
          if (rlo.lt.0d0.and..not.qdamp) then
              write(*,*) 'Error in setvirtual: q2_lo should be < 0!!'
              write(*,*) 'q2_lo = ', q2_lo
              write(*,*) 'st_muren2 = ', st_muren2
              ! set to dummy value:
              	rlo = 1d0
c              stop
          endif      
          lrup = log(rup)
          lrlo = log(rlo) 

          cvirtl = -8d0

                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3,h) = als(1,1)*c2o4pi*
     1                ( mv12(k,isig1,isig3,h) + 
     2                  mm(k,isig1,isig3,h)*(-lrup**2-3*lrup+cvirtl) )
                     mv34(k,isig1,isig3,h) = als(2,1)*c2o4pi*
     1                ( mv34(k,isig1,isig3,h) + 
     2                  mm(k,isig1,isig3,h)*(-lrlo**2-3*lrlo+cvirtl) )
                  else
                     mv12(k,isig1,isig3,h) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3,h)
                     mv34(k,isig1,isig3,h) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3,h)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1  	       mm(k,isig1,isig3,h)   *
     1  	    conjg( mv12(k,isig1,isig3,h)+mv34(k,isig1,isig3,h) )  )
	       endif
	       
	   enddo ! isig3
	 enddo  !isig1 
	enddo !h
	
         if (nlo.eq.0) then
            res(k) = res(k)*9d0
         elseif (nlo.gt.0) then
            if (qdamp) then
               !fakevirt:
                res(k) = res(k)*9d0   ! Born-type
                res(k) = res(k)*0.2d0 ! 'fakevirt' factor
            else !no damping  
c sum of Born+virt:
c            	res(k) = (res(k)+resv(k))*9d0! 9 is the color sum factor
c virt only:
                res(k) = (resv(k))*9d0      ! 9 is the color sum factor
	    endif
         else
	    ! evaluate only parts of virtuals:
	    print*,'you called virt. amplitudes with nlo =',lnlo
	    print*,'this option is not supported by this code version'
	    stop
c	    res(k) = resv(k)*9d0 ! 9 is the color sum factor
         endif !lnlo
     
c      enddo !k
c          
c eliminate contributions from too low qsq in t-channel:
c      if ( (abs(p21(4)).lt.qsqAmin).or.(abs(p43(4)).lt.qsqAmin)) then
      if (qdamp) then 
          res(k) = res(k)*1d-20
      endif

      ans = res(k)

      return
      end


