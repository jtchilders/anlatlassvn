ccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute leptonic tensors 
c
      subroutine compute_tensors(pin)
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
      include 'tensorww.inc'
      real*8 p(0:3,np), v(0:3,nv)
      real*8 q12(0:4,3), q34(0:4,3),qq

      complex*16 zero
      parameter (zero=(0d0,0d0))

      double complex dotrc
      external dotrc
c
c declare local variables
      integer i,j,mu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do mu = 0,3
         p(mu,1) = pin(mu,1)
         p(mu,2) = pin(mu,2)   
         do i = 1,4
            v(mu,i) = pin(mu,i+2)
         enddo !i     
         p(mu,3) = pin(mu,7)
         p(mu,4) = pin(mu,8) 
         p(mu,5) = 0d0   
      enddo ! mu

c compute invariants:
      call calc_invariants(p,v,q12,q34,1)
c
c reset leptonic tensors:
      call vtowpm_reset
c compute leptonic tensors:
         CALL IXXXXX(v(0,1),ZERO , 1,-1,wel) !e+ 
         CALL OXXXXX(v(0,2),ZERO ,-1, 1,wve) !ve 
         CALL IXXXXX(v(0,3),ZERO , 1,-1,wmu) !mu+ 
         CALL OXXXXX(v(0,4),ZERO ,-1, 1,wvu) !vm 
         CALL JIOXXX(wel,wve,GWF ,WMASS,WWIDTH,we)    !W+->e+ve
         CALL JIOXXX(wmu,wvu,GWF ,WMASS,WWIDTH,wu)    !W+->mu+vm

         do mu = 0,3
            qe(mu)  = v(mu,1)+v(mu,2)
            qu(mu)  = v(mu,3)+v(mu,4)
            qww(mu) = qe(mu) + qu(mu)
         enddo
         qe(4)  = qe(0)**2-qe(1)**2-qe(2)**2-qe(3)**2
         qu(4)  = qu(0)**2-qu(1)**2-qu(2)**2-qu(3)**2
         qww(4) = qww(0)**2-qww(1)**2-qww(2)**2-qww(3)**2

	 ! determine "reduced" polarization vectors:
            qq = qe(0)*qww(0)-qe(1)*qww(1)-qe(2)*qww(2)-qe(3)*qww(3)
            xxew = dotrc(qww,we)/qq
            
	    qq = qu(0)*qww(0)-qu(1)*qww(1)-qu(2)*qww(2)-qu(3)*qww(3)
            xxuw = dotrc(qww,wu)/qq  
	    
            do mu = 0,3
               wet(mu+1) = we(mu+1) - xxew*qe(mu)
               wut(mu+1) = wu(mu+1) - xxuw*qu(mu)
           enddo
           do mu = 5,6
               wet(mu) = we(mu)
               wut(mu) = wu(mu)
           enddo
c	 
c leptonic tensors
C for W+W+->e+ve mu+vm
         j = 1

         call wpptowpp(q12(0,j),q34(0,j),v,wppwpp(0,0,j))  	 
 
c t-channel processes:	 
C for V1 V2->e+ve
         call WPVtoWP(1,1,q34(0,j),v,awen(0,0,1,j),zwen(0,0,1,j)) !Wu on upper
         call WPVtoWP(1,2,q12(0,j),v,awen(0,0,2,j),zwen(0,0,2,j)) !Wu on lower
c        
         call WPVtoWP(1,1,q12(0,j),v,waen(0,0,2,j),wzen(0,0,2,j)) !Wu on lower
         call WPVtoWP(1,2,q34(0,j),v,waen(0,0,1,j),wzen(0,0,1,j)) !Wu on upper

c for V1 V2->mu+vm
         call WPVtoWP(2,1,q34(0,j),v,awmu(0,0,1,j),zwmu(0,0,1,j)) !We on upper
         call WPVtoWP(2,2,q12(0,j),v,awmu(0,0,2,j),zwmu(0,0,2,j)) !We on lower
c        
         call WPVtoWP(2,1,q12(0,j),v,wamu(0,0,2,j),wzmu(0,0,2,j)) !We on lower
         call WPVtoWP(2,2,q34(0,j),v,wamu(0,0,1,j),wzmu(0,0,1,j)) !We on upper

      return
      end
