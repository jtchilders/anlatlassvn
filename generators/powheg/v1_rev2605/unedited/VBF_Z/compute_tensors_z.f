ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute leptonic tensors for VBF pp->Zjj 
c
      subroutine compute_tensors_z(pin)
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
      include 'nlegborn.h'
      include 'pwhg_math.h'
c
      integer nlegs
      parameter (nlegs=nlegborn)
      real*8 pin(0:3,nlegs)
c
c vbfnlo stuff:
      include 'global.inc'
      include 'coupl.inc'
      include 'tensorl-hel.inc'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3),q34(0:4,3),qww(0:3),qq

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu

      integer h,hmin,hmax,hstep
      common /hval/hmin,hmax,hstep
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      hmin = -1 !helicity range for decay leptons
      hmax =  1
      hstep = 2
c
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,2
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,5)
         p(mu,4) = pin(mu,6) 
         p(mu,5) = 0d0   
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset leptonic tensors:
      call vtoll_reset
c compute leptonic tensors:

         do mu = 0,3
            qe(mu) = v(mu,1)+v(mu,2)
         enddo
         qe(4) = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2

	do h = hmin,hmax,hstep ! sum over lepton helicities
	
         CALL IXXXXX(v(0,1),ZERO ,+h,-1,lep(1,h))	      !e+	
         CALL OXXXXX(v(0,2),ZERO ,-h, 1,lem(1,h))	      !e- 

	 CALL JIOXXX(lep(1,h),lem(1,h),GZL ,ZMASS,ZWIDTH,ze(1,h))     !Zl
	 CALL JIOXXX(lep(1,h),lem(1,h),GAL ,ZERO ,ZERO  ,ae(1,h))     !Al	
	 
c leptonic tensors
         j = 1

C V1V2 --> e+ e- (Vi=A,Z) tensors for NC process (k=1...4)
	 call VVtoll(1,h,q12(0,j),v,aaee(0,0,j,h),azee(0,0,j,h),
     #		 zaee(0,0,j,h),zzee(0,0,j,h))    	

C W+W- --> e+ e- for CC (k=5)
         call WWtoll(1,h,q12(0,j),v,CCee(0,0,j,h))   
	 	 
C W-W+ --> e+ e- for CC (k=6)	
         call WWtoll(1,h,q34(0,j),v,CCee6(0,0,j,h))  

	enddo !h

      return
      end
