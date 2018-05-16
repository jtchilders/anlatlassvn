      subroutine spinormu(N,p,zamub)
c---Calculate spinor products dotted in with a vector k
c---extended to deal with negative energies ie with all momenta outgoing
C   zamub=<i-|mu|j-> 
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),rt(mxpart),k(4),kp,km,flip(mxpart),
     & mu1(4),mu2(4),mu3(4),mu4(4),metric(4)
      double complex pr(mxpart),pl(mxpart),f(mxpart),kr,kl,
     & zamub(mxpart,4,mxpart)
      integer i,j,N,mu
      data mu1/1d0,0d0,0d0,0d0/      
      data mu2/0d0,1d0,0d0,0d0/      
      data mu3/0d0,0d0,1d0,0d0/      
      data mu4/0d0,0d0,0d0,1d0/      
      data metric/-1d0,-1d0,-1d0,+1d0/      
      save mu1,mu2,mu3,mu4,metric
C--setup components for vector which is contracted in
      do mu=1,4
      if (mu == 1) k(:)=mu1(:)
      if (mu == 2) k(:)=mu2(:)
      if (mu == 3) k(:)=mu3(:)
      if (mu == 4) k(:)=mu4(:)
      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=dcmplx(+k(3),-k(2))
      kl=dcmplx(+k(3),+k(2))

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
            zamub(j,mu,j)=2d0*p(j,mu)
            
C-----positive energy case
            if (p(j,4) .gt. 0d0) then
            flip(j)=1d0
            f(j)=cone
            else
            flip(j)=-1d0
            f(j)=im
            endif
            rt(j)=dsqrt(flip(j)*(p(j,4)+p(j,1)))
            pr(j)=dcmplx(flip(j)*p(j,3),-flip(j)*p(j,2))
            pl(j)=Dconjg(pr(j))
      enddo
      do i=1,N
         do j=1,i
         zamub(i,mu,j)=metric(mu)*f(i)*f(j)
     & *(pr(i)*pl(j)*dcmplx(kp/(rt(i)*rt(j)))
     &    -pr(i)*kl*dcmplx(rt(j)/rt(i))
     &    -dcmplx(rt(i)/rt(j))*kr*pl(j)+dcmplx(rt(i)*rt(j)*km))
         zamub(j,mu,i)=flip(i)*flip(j)*Dconjg(zamub(i,mu,j))

      enddo
      enddo
      enddo

      return
      end
