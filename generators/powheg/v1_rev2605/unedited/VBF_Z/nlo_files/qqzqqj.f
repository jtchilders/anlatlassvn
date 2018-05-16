c      
      subroutine qqzqqj(pbar,sign,qbar,gsign,k,ans)
      implicit none
c
c	Last modified for POWHEG by Barbara Jaeger:  June 2012
C
C  QQZQQJ calculates the matrix elements**2 for electroweak
c  weak boson production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 g e+ e-
C
C  Crossing related processes are also computed. Pauli interference terms for
c  identical fermions are neglected. In particular, only the
c  t-channel exchange of elctroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
c
c
C  This code includes only real emission contributions, i.e.
c
c      return uucc = |M_real|^2   etc.
c
c	fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c		upper (12) line (j=2) or lower (34) line (j=3)
c	l is the gluon polarization in the kartesian basis (l=1,2)
c		l=0 stands for building blocks without gluon emission
c	k is the process ID (1:uucc,2:uuss,3:ddcc;4:ddss;5:udsc;6:ducs)
c	isig1/isig3 are the helicities of parton1/3 
c
c	lepton helicities are summed over in this routine 
c
c---------------------------------------------------------------------
c
c
      include '../../include/pwhg_math.h'
      include '../../include/pwhg_st.h'
      include 'global.inc'
      include 'tensorl-hel.inc'
c
c electroweak couplings are taken from KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c
c alfas, scales etc
      include 'scales.inc'
c
c variables for the main part of the program
c
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,4+nv),qbar(0:4),musq
      double precision ans(3)
      double precision res(6,2:3),resv(6)
      double precision  p(0:3,4+nv),q(0:4), p21(0:4,2:3), p43(0:4,2:3),
     1                  pq(0:4,4)
      integer  sign(4+nv),gsign, mu, i, j,jj, k, kk, id,kl,
     1         isig, isig1, isig3,is,ka,kb,kf,il,ib,ll,bos
      integer  ifl(4,6), js1, js3, is1, is3
      integer  l   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)
      integer jmin, jmax
      logical jlog2,jlog3
      double complex prop21(4,2:3), prop43(4,2:3)!,
c     1		     prop_ee(4,2:3),prop_uu(4,2:3)
      double complex mat(6,-1:1,-1:1,2:3,0:2,-1:1,3)
      integer id1,id2,id3
      double complex mm(6,-1:1,-1:1,2:3,2,-1:1)
      double complex maa, maz, mza, mzz, mzz5, mzz6
      double complex  m1,m2 ! for checks only
      double complex  matot,ratot,raz,raa! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4),jqq(0:5,-1:1,2,-1:1,0:2), 
     1 		     braketg(2,-1:1,4,2), jh1(0:5,-1:1), jh2(0:5,-1:1)
      double complex psize(2,-1:1,-1:1,4),psiae(2,-1:1,-1:1,4),
     1               jez(0:5,-1:1,-1:1,4,0:2,-1:1),jea(0:5,-1:1,-1:1,4,0:2,-1:1)
      double precision fqze(0:4,4), fqae(0:4,4),
     1 		       dummy(0:4)
      double complex braketgze(2,-1:1,4,0:2,-1:1),braketzeg(2,-1:1,4,0:2,-1:1),
     3		     braketgae(2,-1:1,4,0:2,-1:1),braketaeg(2,-1:1,4,0:2,-1:1)  
      double precision  pgze(0:4,4),pzeg(0:4,4),
     3			pgae(0:4,4),paeg(0:4,4)
      double precision  pga(0:4,4),pgz(0:4,4),
     3			pag(0:4,4),pzg(0:4,4)
      double complex psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4)
      double complex jaegi(0:5,-1:1),jgaei(0:5,-1:1),jaeg0i(0:5,-1:1),
     1		     jaegii(0:5,-1:1),jgaeii(0:5,-1:1),jaeg0ii(0:5,-1:1),
     1		     jzegi(0:5,-1:1),jgzei(0:5,-1:1),jzeg0i(0:5,-1:1),
     1		     jzegii(0:5,-1:1),jgzeii(0:5,-1:1),jzeg0ii(0:5,-1:1)
      double complex zee(4:5,2:3)
      double complex m1aae,m1aze,m1zae,m1zze,z1aze,z1zze,
     1 		     m2aae,m2aze,m2zae,m2zze,m1e,m2e,
     1		     mz1aae,mz1aze,mz1zae,mz1zze,
     1		     mz2aae,mz2aze,mz2zae,mz2zze,
     1		     m3aae,m3aze,m3zae,m3zze,z3aze,z3zze,
     3 		     m4aae,m4aze,m4zae,m4zze,m3e,m4e,
     3		     mz3aae,mz3aze,mz3zae,mz3zze,
     3		     mz4aae,mz4aze,mz4zae,mz4zze
      double complex mezz(4,-1:1,-1:1,2:3,2,-1:1),meza(4,-1:1,-1:1,2:3,2,-1:1),
     #		     meaz(4,-1:1,-1:1,2:3,2,-1:1),meaa(4,-1:1,-1:1,2:3,2,-1:1),
     #		     zezz(-1:1,-1:1,2:3,2:3,2),zeza(-1:1,-1:1,2:3,2:3,2),
     #		     zeaz(-1:1,-1:1,2:3,2:3,2),zeaa(-1:1,-1:1,2:3,2:3,2)
      double complex  ma(4,2:3,0:2), mz(4,2:3,0:2)
      double complex propt(-1:1,-1:1,6,2:3,2),
     1		     propbbe(2:3), propbbu(2:3),
     1		     propbbea(6,-1:1,-1:1,2:3),
     1		     propbbez(6,-1:1,-1:1,2:3),
     1		     propbbua(6,-1:1,-1:1,2:3),
     1		     propbbuz(6,-1:1,-1:1,2:3)
      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1 		     dotcc, dotrc, dotqj, s1c
      external contract_Tjj, !contract_T1j,contract_T2j,
     1 	       dotcc, dotrc, dotqj, s1c
      integer lh
      double complex m1hc(-1:1),m1kc(2),m1hh(-1:1),m1kh(2)
      double complex m1hj(-1:1,2:3),m1kj(2,2:3)
      double complex m1haj(-1:1,2:3),m1hzj(-1:1,2:3),
     1 		     m1hla(-1:1,2:3),
     1 		     m1hlz(-1:1,2:3),
     1 		     m1kaj(2,2:3),m1kzj(2,2:3),
     1 		     m1kua(2,2:3),m1kla(2,2:3),
     1 		     m1kuz(2,2:3),m1klz(2,2:3)
      double complex fac
      double complex im
      parameter (im = (0d0,1d0))
      logical ldebug, linit,lerror,lgauge,ldebugm
      integer n,m

      logical nc_type,cc_type
      
      integer h,hmin,hmax,hstep
      common /hval/hmin,hmax,hstep
     
      save ifl, linit, lerror,ldebug,ldebugm,lgauge
	
      data linit /.true./, lerror /.false./, ldebug /.false./, 
     $	   ldebugm /.false./, lgauge /.false./
     
      logical vdebug,bdebug
      parameter (vdebug = .false.,bdebug = .false.)

      logical v_only,b_only
      parameter (v_only = .false.,b_only=.false.)
c      
c  ---------------------------------------------------------------------
c  ---------------------------------------------------------------------
c 
c initialize & precompute stuff needed below:
c
c fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)

c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         print*,'qsqAmin in real emission=',qsqamin
         print*,'damping factor of 1d-20 below qsqAmin '
         if (v_only) print*,'real graphs with vertex topology only'
         if (b_only) print*,'real graphs with box topology only'
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
      
      if (gsign.eq.1) then	!final state gluon
        jlog2 = .true.		! can couple to upper/lower line
	jlog3 = .true. 
      	jmin = 2
	jmax = 3
      else		       !initial state gluon -> only:
	 if (sign(1).ne.sign(2)) then  !gluon from upper line
	       jlog2 = .true.
	       jlog3 = .false. 
	       jmin = 2
	       jmax = 2
	 else			       !gluon from lower line
	       jlog2 = .false.
	       jlog3 = .true. 
	       jmin = 3
	       jmax = 3
       endif
      endif	
                  
      do kl = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
	       do j = 2,3
	          do l = 0,2
                  do h = -1,1,2
                     do i = 1,3
                  	mat(kl,isig1,isig3,j,l,h,i)  = 0
     		     enddo
		  enddo !h   
		  enddo	
               enddo
            enddo
         enddo
      enddo
c
c identify fermion line sign factors (for 1 3 -> 2 4 etc.)
c
c      is1 = sign(1)
c      is3 = sign(3)
c      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
c      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1

c fix is1 such that is1 = +1 for   q1 ->  q2 g 
c		    is1 = -1 for  q1b -> q2b g 
c		    is1 =  0 for    g -> q1b q2
c
c      for is3:     is3 = +1 for   q3 ->  q4 g 
c		    is3 = -1 for  q3b -> q4b g 
c		    is3 =  0 for    g -> q3b q4

	is1 = (sign(1)+sign(2))/2  
	is3 = (sign(3)+sign(4))/2
		
c (is1,is3 are fixed here and don't change throughout this run of the program)


c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
	 q(mu) = qbar(mu)*gsign
	 
         p21(mu,3) = p(mu,2) - p(mu,1)	! g from 34 line
         p21(mu,2) = p21(mu,3) + q(mu) 	! g from 12 line
         p43(mu,2) = p(mu,4) - p(mu,3)
         p43(mu,3) = p43(mu,2) + q(mu)	 
	 
c	 pez(mu) = qe(mu)
      enddo
      
      do j= 2,3
        p21(4,j) = p21(0,j)**2 - p21(1,j)**2 - p21(2,j)**2 - p21(3,j)**2
        p43(4,j) = p43(0,j)**2 - p43(1,j)**2 - p43(2,j)**2 - p43(3,j)**2     
      enddo 

c  ---------------------------------------------------------------------
c
c get the vector boson propagator factors
c
c depending on value of j, gluon is attached to respective quark line or not;
c no V is attached here 

      do j = 2,3
         prop21(1,j) = 1/p21(4,j)	      
         prop21(2,j) = 1/dcmplx(p21(4,j)-xm2(2),xmg(2))
         prop21(3,j) = 1/dcmplx(p21(4,j)-xm2(3),xmg(3))
         prop21(4,j) = prop21(3,j)

         prop43(1,j) = 1/p43(4,j)
         prop43(2,j) = 1/dcmplx(p43(4,j)-xm2(2),xmg(2))
         prop43(3,j) = 1/dcmplx(p43(4,j)-xm2(3),xmg(3))
         prop43(4,j) = prop43(3,j)
      enddo
c
c  ---------------------------------------------------------------------

c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents (with no gluon radiation) 
c	J21^mu=jqq(mu,isig1,1,is1,0), J43^mu=jqq(mu,isig3,2,is3,0) 
c
        call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),
     #						jqq(0,-1,1,is1,0))      
        call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),
     #						jqq(0,-1,2,is3,0))

c
c  Get the gluon polarization vector and the gluon emission spinors
      do l = 1,2	! 2 gluon polarizations
         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors
	 
c for gauge check:
        if (lgauge) then ! contract amplitude with q rather than eps(q)
	 do mu = 0,3
	 	eps(mu,l) = qbar(mu)
	 enddo		 
	endif
	 	 
         do isig = -1,1,2	! fermion helicity 
 
c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
            
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig1
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2
            call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,eps(0,l),
     $           braketg(1,isig,3,l),pq(0,3))      ! |q,3>_l,isig3
            call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,eps(0,l),
     $           braketg(1,isig,4,l),pq(0,4))      ! <4,q|_l,isig4

         enddo !isig
	          
c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)
 
      enddo !l            
       
c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:
c
c	gluon from upper line:
      do l = 1, 2 ! gluon polarizations
         call curr6(1,psi(1,-1,2),p(0,2),braketg(1,-1,1,l),pq(0,1),jh1)	
c                                            =  <2|vertex|q,1>_l,isig1
         call curr6(1,braketg(1,-1,2,l),pq(0,2),psi(1,-1,1),p(0,1),jh2)	
c                                            =  <2,q|vertex|1>_l,isig1
         do isig = -1,1,2 ! fermion helicity
            do mu = 0,5
	       jqq(mu,isig,1,is1,l) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig1
            enddo !mu
         enddo !isig
         

c	gluon from lower line:
         call curr6(1,psi(1,-1,4),p(0,4),braketg(1,-1,3,l),pq(0,3),jh1)	
c                                           =   <4|gam.mu|q,3>_l,isig3
         call curr6(1,braketg(1,-1,4,l),pq(0,4),psi(1,-1,3),p(0,3),jh2)	
c                                           =   <4,q|gam.mu|3>_l,isig3
         do isig = -1,1,2
            do mu = 0,5
               jqq(mu,isig,2,is3,l) = jh1(mu,isig) + jh2(mu,isig)
c                            = (<4|gam.mu|q,3>+<4,q|gam.mu|3>)_l,isig3
            enddo !mu
         enddo !isig
      enddo !l

      if (b_only) goto 444 

c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c for (V) topology:
c
c neutral current:
c
      if (nc_type) then

      do h = hmin,hmax,hstep
      do l = 1,2	! gluon polarization
        do isig1 = -1,1,2  ! fermion1 helicity
          do isig3 = -1,1,2  ! fermion3 helicity
	 
           if (jlog2) then
	   j = 2 ! g from upper line
	   
	   maa = contract_Tjj(aaee(0,0,j,h),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           maz = contract_Tjj(azee(0,0,j,h),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           mza = contract_Tjj(zaee(0,0,j,h),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))
           mzz = contract_Tjj(zzee(0,0,j,h),jqq(0,isig1,1,is1,l),
     #	   					jqq(0,isig3,2,is3,0))

c           do k = 1,4
              mat(k,isig1,isig3,j,l,h,1) = 
     1    	   maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2    	 + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3    	 + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4    	 + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
c           enddo !k
	   endif !jlog2
	 
           if (jlog3) then
	   j = 3 ! g from lower line
	   
	   maa = contract_Tjj(aaee(0,0,j,h),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           maz = contract_Tjj(azee(0,0,j,h),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           mza = contract_Tjj(zaee(0,0,j,h),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))
           mzz = contract_Tjj(zzee(0,0,j,h),jqq(0,isig1,1,is1,0),
     #	   					jqq(0,isig3,2,is3,l))

c           do k = 1,4
              mat(k,isig1,isig3,j,l,h,1) = 
     1    	   maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2    	 + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3    	 + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4    	 + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
c           enddo !k
	   endif !jlog3
c
          enddo !isig3
        enddo !isig1
      enddo !l
      
      enddo !h 
      endif !nc
     
c 
c -----------------------
c
c charged current (k=5,6):
c      
      if (cc_type) then

      do h = hmin,hmax,hstep
      do l = 1,2	! gluon polarization (kart. basis)
      
      if (jlog2) then
      j =  2 ! g from upper line

      if (k.eq.5) then
         mzz5 = contract_Tjj(CCee(0,0,j,h),jqq(0,-1,1,is1,l),
     #						jqq(0,-1,2,is3,0))
         mat(5,-1,-1,j,l,h,1) = mzz5*clr(3,3,-1)**2
      elseif (k.eq.6) then
         mzz6 = contract_Tjj(CCee6(0,0,j,h),jqq(0,-1,2,is3,0),
     #      					jqq(0,-1,1,is1,l))
         mat(6,-1,-1,j,l,h,1) = mzz6*clr(3,3,-1)**2  
      endif !k   

      endif !jlog2
  
      if (jlog3) then
      j =  3 ! g from lower line
            
      if (k.eq.5) then
         mzz5 = contract_Tjj(CCee(0,0,j,h),jqq(0,-1,1,is1,0),
     #						jqq(0,-1,2,is3,l))
         mat(5,-1,-1,j,l,h,1) = mzz5*clr(3,3,-1)**2

      elseif (k.eq.6) then
         mzz6 = contract_Tjj(CCee6(0,0,j,h),jqq(0,-1,2,is3,l),
     #						jqq(0,-1,1,is1,0))
         mat(6,-1,-1,j,l,h,1) = mzz6*clr(3,3,-1)**2    
      endif!k

      endif !jlog3
      enddo ! l-loop
      enddo !h
      endif !cc

 444  continue
     
      if (v_only) goto 999 ! skip box topology
c 
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c  prepare box diagrams: attach V to external spinors
c 
c      isig = +-1   : left- and righthanded spinors coupling to A/Z's
c
c  Notation for virtual 2-component spinors and momenta
c
c  Z -> e+e- attached to quark number i: psize(*,isig,h,i) with momentum fqze(mu,i)
c  A -> e+e- attached to quark number i: psiae(*,isig,h,i) with momentum fqae(mu,i)
c  
c  the fermion current corresponding to a quark line with the real emitted W+
c  attached next to quark number i is stored in jwp(mu,isig,is,i,l). Similarly
c  jwm(mu,isig,i,l) is the corresponding current for a W- attached next to 
c  quark i; for l = 0 jwp/jwm doesn't contain a gluon; for l = 1,2 jwp/jwm 
c  includes sum of possible gluon couplings on quark leg i
c 
      l = 0    ! no gluon 
 
 	do i = 1,3,2  ! fermion ID
		
		if (i.eq.1) then
		   is = is1
		else
		   is = is3   
		endif   
        
	 do h = hmin,hmax,hstep
	 do isig = -1,1,2
	
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,ze(1,h),
     1              psize(1,isig,h,i),fqze(0,i))
         call ket2c(psi(1,isig,i),.true.,p(0,i),isig,qe,ae(1,h),
     1              psiae(1,isig,h,i),fqae(0,i))

        
	 call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,ze(1,h),
     1              psize(1,isig,h,i+1),fqze(0,i+1))
         call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,qe,ae(1,h),
     1              psiae(1,isig,h,i+1),fqae(0,i+1))


	enddo !isig     
	  
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psize(1,-1,h,i),fqze(0,i), jez(0,-1,is,i,l,h)   )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1                 psiae(1,-1,h,i),fqae(0,i), jea(0,-1,is,i,l,h)   )

         call curr6(1,psize(1,-1,h,i+1),fqze(0,i+1),
     1                 psi(1,-1,i),p(0,i), jez(0,-1,is,i+1,l,h) )
         call curr6(1,psiae(1,-1,h,i+1),fqae(0,i+1),
     1                 psi(1,-1,i),p(0,i), jea(0,-1,is,i+1,l,h) )
     
      
	enddo !h
 	enddo ! i loop

c
c -----------------------------------------------------------------------
c
c keep structure of LO code, but replace jez etc. with 
c jez(mu,isig,is,i,l=0:2,h) etc. to take gluon radiation into account


c attach A/Z to f1 or f2 / f3 or f4:

 	do i = 1,3,2  ! fermion ID (isigi=-1 or 1)  

		if (i.eq.1) then
		   is = is1
		else
		   is = is3   
		endif  
 
 ! gluon radiation from fermion i / i+1
 	 do l = 1,2
	 do h = hmin,hmax,hstep

	 do isig = -1,1,2
            	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                     isig,qe,ze(1,h),braketgze(1,isig,i,l,h),pgze(0,i))
            	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,qe,ze(1,h),braketgze(1,isig,i+1,l,h),pgze(0,i+1))
		
		call ket2r(psize(1,isig,h,i),.false.,fqze(0,i),isig,
     $	    		    q,eps(0,l),braketzeg(1,isig,i,l,h),pzeg(0,i))      
            	call bra2r(psize(1,isig,h,i+1),.false.,fqze(0,i+1),isig,
     $	    		q,eps(0,l),braketzeg(1,isig,i+1,l,h),pzeg(0,i+1))      

            	call ket2c(braketg(1,isig,i,l),.false.,pq(0,i),
     $                     isig,qe,ae(1,h),braketgae(1,isig,i,l,h),pgae(0,i))
            	call bra2c(braketg(1,isig,i+1,l),.false.,pq(0,i+1),
     $                 isig,qe,ae(1,h),braketgae(1,isig,i+1,l,h),pgae(0,i+1))
		
		call ket2r(psiae(1,isig,h,i),.false.,fqae(0,i),isig,
     $	    		    q,eps(0,l),braketaeg(1,isig,i,l,h),paeg(0,i))      
            	call bra2r(psiae(1,isig,h,i+1),.false.,fqae(0,i+1),isig,
     $	    		q,eps(0,l),braketaeg(1,isig,i+1,l,h),paeg(0,i+1))      
	 enddo !isig


   	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    			braketgze(1,-1,i,l,h),pgze(0,i),jzegi) 			
      	    call curr6(1,braketgze(1,-1,i+1,l,h),pgze(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jzegii)		       
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			        braketzeg(1,-1,i,l,h),pzeg(0,i),jgzei)
      	    call curr6(1,braketzeg(1,-1,i+1,l,h),pzeg(0,i+1),	
     $			        psi(1,-1,i),p(0,i),jgzeii)
   
    	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $	    			braketgae(1,-1,i,l,h),pgae(0,i),jaegi) 			
      	    call curr6(1,braketgae(1,-1,i+1,l,h),pgae(0,i+1),
     $	    			psi(1,-1,i),p(0,i),jaegii)		       
      	    call curr6(1,psi(1,-1,i+1),p(0,i+1),
     $			        braketaeg(1,-1,i,l,h),paeg(0,i),jgaei)
      	    call curr6(1,braketaeg(1,-1,i+1,l,h),paeg(0,i+1),	
     $			        psi(1,-1,i),p(0,i),jgaeii)

 
  ! gluon radiation from fermion i+1 / i
     	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    			 psize(1,-1,h,i),fqze(0,i),jzeg0i)
     	    call curr6(1,psize(1,-1,h,i+1),fqze(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jzeg0ii)
     	    
	    call curr6(1,braketg(1,-1,i+1,l),pq(0,i+1),
     $	    			 psiae(1,-1,h,i),fqae(0,i),jaeg0i)
     	    call curr6(1,psiae(1,-1,h,i+1),fqae(0,i+1),
     $	    			braketg(1,-1,i,l),pq(0,i),jaeg0ii)

   	 
	    do mu = 0,5
	    do isig = -1,1,2
	    	   
 		   jez(mu,isig,is,i,l,h) = jzegi(mu,isig)+
     $			jgzei(mu,isig)+jzeg0i(mu,isig)  ! Ze & g emission from i/i+1 line 
  		   	
		   jez(mu,isig,is,i+1,l,h) = jzegii(mu,isig)+
     $			jgzeii(mu,isig)+jzeg0ii(mu,isig)   
	    	    		   
 		   jea(mu,isig,is,i,l,h) = jaegi(mu,isig)+
     $			jgaei(mu,isig)+jaeg0i(mu,isig)    ! Ae & g emission from i/i+1 lin
		   
		   jea(mu,isig,is,i+1,l,h) = jaegii(mu,isig)+
     $			jgaeii(mu,isig)+jaeg0ii(mu,isig)  	    	   
    	    	   
	    enddo !isig
	    enddo  !mu       

	enddo !h       
       enddo ! l = 1,2
c
       enddo ! i loop     
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c
c for (B) topology:
c
c
	do isig1 = -1,1,2
	   do isig3 = -1,1,2
c
	   
	   do il = 1,2 ! two possible g pols.
	   
	   do j = jmin,jmax  ! g from upper or lower line
	   
	   if (j.eq.2) then
	       l = il  ! loop over pol. of g from upper line
	       ll = 0  ! no g from lower line
	   else !j=3
	       ll = il  ! loop over pol. of g from lower line
	       l = 0    ! no g from upper line
	   endif
c
c e+e- from upper line:
        do h = hmin,hmax,hstep
c	
	mezz(1,isig1,isig3,j,il,h) = 
     &          dotcc(jez(0,isig1,is1,1,l,h),jqq(0,isig3,2,is3,ll))
	meaa(1,isig1,isig3,j,il,h) = 
     &          dotcc(jea(0,isig1,is1,1,l,h),jqq(0,isig3,2,is3,ll))
c	
	mezz(2,isig1,isig3,j,il,h) =
     &           dotcc(jez(0,isig1,is1,2,l,h),jqq(0,isig3,2,is3,ll))
	meaa(2,isig1,isig3,j,il,h) = 
     &          dotcc(jea(0,isig1,is1,2,l,h),jqq(0,isig3,2,is3,ll))
c	
	mezz(3,isig1,isig3,j,il,h) = 
     &          dotcc(jqq(0,isig1,1,is1,l),jez(0,isig3,is3,3,ll,h))
	meaa(3,isig1,isig3,j,il,h) = 
     &          dotcc(jqq(0,isig1,1,is1,l),jea(0,isig3,is3,3,ll,h))
c	
	mezz(4,isig1,isig3,j,il,h) = 
     &          dotcc(jqq(0,isig1,1,is1,l),jez(0,isig3,is3,4,ll,h))
	meaa(4,isig1,isig3,j,il,h) = 
     &          dotcc(jqq(0,isig1,1,is1,l),jea(0,isig3,is3,4,ll,h))
	
	enddo !h	
c
c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators:
c  not needed as resulting contributions cancel in sum of all diagrams	
	
        enddo !j
	enddo !il

c -------------------------------

      if (nc_type) then

	do j = jmin,jmax
  
c        do k = 1,4 
                  
      	      propbbez(k,isig1,isig3,j) = 
     1 	     	prop43(2,j)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
     
      	      propbbea(k,isig1,isig3,j) = 
     1 	     	prop43(1,j)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)    
                  
      	      propbbuz(k,isig1,isig3,j) = 
     1 	     	prop21(2,j)*clr(ifl(2,k),2,isig1)*clr(ifl(3,k),2,isig3)
     
      	      propbbua(k,isig1,isig3,j) = 
     1 	     	prop21(1,j)*clr(ifl(2,k),1,isig1)*clr(ifl(3,k),1,isig3)    

	
	do l = 1,2
        do h = hmin,hmax,hstep
     
c for NC make use of: ifl(1,k) = ifl (2,k) and ifl(3,k) = ifl (4,k) 
c 	gauge terms cancel in sum of NC contributions from leg1&2/3&4 
	     
      mat(k,isig1,isig3,j,l,h,2) = 
     1 ((mezz(1,isig1,isig3,j,l,h)+mezz(2,isig1,isig3,j,l,h))*
     1		clr(ifl(1,k),2,isig1) + 
     1	(meaa(1,isig1,isig3,j,l,h)+meaa(2,isig1,isig3,j,l,h))*
     1		clr(ifl(1,k),1,isig1))*
     1 (propbbea(k,isig1,isig3,j)+propbbez(k,isig1,isig3,j))
      
      
       mat(k,isig1,isig3,j,l,h,3) = 
     3 ((mezz(3,isig1,isig3,j,l,h)+mezz(4,isig1,isig3,j,l,h))*
     3	       clr(ifl(4,k),2,isig3) + 
     3  (meaa(3,isig1,isig3,j,l,h)+meaa(4,isig1,isig3,j,l,h))*
     3	       clr(ifl(4,k),1,isig3))*
     3	(propbbua(k,isig1,isig3,j)+propbbuz(k,isig1,isig3,j))
     
	enddo !h
	enddo !l
c	enddo !k
	enddo !j

      endif !nc  
	
      enddo !isig3
      enddo !isig1	

c -----------------

      if (cc_type) then

      do j = jmin,jmax 

c      do k = 5,6 ! charged current
      isig1 = -1
      isig3 = -1
            	    
      propbbe(j) =  prop43(3,j)*clr(3,3,-1)**2
      propbbu(j) =  prop21(3,j)*clr(3,3,-1)**2
	    
      do l = 1,2 
        do h = hmin,hmax,hstep
       
      mat(k,isig1,isig3,j,l,h,2) = 
     1 ((mezz(1,isig1,isig3,j,l,h))*clr(ifl(1,k),2,isig1) + 
     1	(meaa(1,isig1,isig3,j,l,h))*clr(ifl(1,k),1,isig1))*
     1		propbbe(j)    
     2 +
     2 ((mezz(2,isig1,isig3,j,l,h))*clr(ifl(2,k),2,isig1) + 
     2	(meaa(2,isig1,isig3,j,l,h))*clr(ifl(2,k),1,isig1))*
     2		propbbe(j)  

      mat(k,isig1,isig3,j,l,h,3) =
     3 ((mezz(3,isig1,isig3,j,l,h))*clr(ifl(3,k),2,isig3) + 
     3  (meaa(3,isig1,isig3,j,l,h))*clr(ifl(3,k),1,isig3))*
     3		propbbu(j)      
     4 +
     4 ((mezz(4,isig1,isig3,j,l,h))*clr(ifl(4,k),2,isig3) + 
     4  (meaa(4,isig1,isig3,j,l,h))*clr(ifl(4,k),1,isig3))*
     4		propbbu(j)    

	  enddo !h
	enddo !l  
	
c      enddo !k

	enddo !j

      endif !nc/cc  
c
c
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------

c sum the graphs, square them and map them onto uucc, uuss etc.

c i = 1		V
c i = 2,3       B
 
 999	continue
 	if (b_only) then
	do kl = 1,6
	  do isig1 = -1,1,2
	  do isig3 = -1,1,2
            do j = 2,3
	      do h = hmin,hmax,hstep
		do l = 1,2
		   mat(kl,isig1,isig3,j,l,h,1) = 0d0
		enddo !l
	      enddo !h
	    enddo !j
	  enddo !isig3
          enddo !isig1
	enddo !kl     
	endif

      if (nc_type) then  
c      do k = 1,4
 	 do j = 2,3     
 	    res(k,j) = 0
            do h = hmin,hmax,hstep
	    do isig1 = -1,1,2
	       do isig3 = -1,1,2
 	    	  do l = 1,2
              	     mm(k,isig1,isig3,j,l,h) = 0
               	     do i = 1,3
                        mm(k,isig1,isig3,j,l,h) = 
     1                 		 mm(k,isig1,isig3,j,l,h) + 
     1		    	         mat(k,isig1,isig3,j,l,h,i)
      		     enddo !i
              	     res(k,j) = res(k,j) 
     &		       	       + dreal(mm(k,isig1,isig3,j,l,h))**2
     &                         + dimag(mm(k,isig1,isig3,j,l,h))**2
	          enddo !l
	      enddo !isig3		     
           enddo !isig1
	   enddo !h
           res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
	 enddo !j
c      enddo !k      
 
      else !cc          
     
c      do k = 5,6
 	 do j = 2,3     
 	    res(k,j) = 0
            do h = hmin,hmax,hstep
 	    do l = 1,2
               mm(k,-1,-1,j,l,h) = 0
               do i = 1,3
            	  mm(k,-1,-1,j,l,h) = 
     1      		   mm(k,-1,-1,j,l,h) + 
     1	    		   mat(k,-1,-1,j,l,h,i)
     	       enddo !i
               res(k,j) = res(k,j) 
     &	    		 + dreal(mm(k,-1,-1,j,l,h))**2
     &      		 + dimag(mm(k,-1,-1,j,l,h))**2
	    enddo !l
	    enddo !h
            res(k,j) = res(k,j)*12d0*fpials(j)   ! C_2*9 is the color factor
 	 enddo !j
c      enddo !k     

      endif !nc/cc     
            
      if (jmin.eq.3) then
         res(k,2) = 0d0
      elseif (jmax.eq.2) then
         res(k,3) = 0d0
      endif	

      if ((abs(p21(4,2)).lt.qsqAmin).or.(abs(p21(4,3)).lt.qsqAmin).or.
     &    (abs(p43(4,2)).lt.qsqAmin).or.(abs(p43(4,3)).lt.qsqAmin)) then       
         res(k,3) = res(k,3)*0.5d-20
         res(k,2) = res(k,2)*0.5d-20
      endif

		    
      do j = 2,3
         ans(j) = res(k,j)
      enddo
      ans(1) = res(k,3)+res(k,2)  
 
      return
      end

