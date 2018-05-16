      subroutine fillBmunu(p1,p2,p3,p4,p5,p6,p,Bmunu)
      implicit none
      include 'constants.f'
      integer mu,nu,p1,p2,p3,p4,p5,p6
      double precision Bmunu(4,4,6),B5(4,4),B6(4,4),p(mxpart,4)
      double complex fab(4,2),fba(4,2),fqed(4,2)
      Bmunu(:,:,:)=0d0
      call W2jetampmu(p1,p2,p3,p4,p5,p6,fab,fba,fqed)
      do mu=1,4
      do nu=1,4
      B5(mu,nu)=
     & +dble(fab(mu,1)*dconjg(fab(nu,1))
     &      +fab(mu,2)*dconjg(fab(nu,2)))
     & +dble(fba(mu,1)*dconjg(fba(nu,1))
     &      +fba(mu,2)*dconjg(fba(nu,2)))
     & -dble(fqed(mu,1)*dconjg(fqed(nu,1))
     &      +fqed(mu,2)*dconjg(fqed(nu,2)))/xnsq
      enddo
      enddo
      call W2jetampmu(p1,p2,p3,p4,p6,p5,fab,fba,fqed)
      do mu=1,4
      do nu=1,4
      B6(mu,nu)=
     & +dble(fab(mu,1)*dconjg(fab(nu,1))
     &      +fab(mu,2)*dconjg(fab(nu,2)))
     & +dble(fba(mu,1)*dconjg(fba(nu,1))
     &      +fba(mu,2)*dconjg(fba(nu,2)))
     & -dble(fqed(mu,1)*dconjg(fqed(nu,1))
     &      +fqed(mu,2)*dconjg(fqed(nu,2)))/xnsq
      enddo
      enddo

      Bmunu(:,:,p5)=B5(:,:)
      Bmunu(:,:,p6)=B6(:,:)


C--optional gauge choice A4=0
      do mu=1,4
      do nu=1,4
      if ((mu == 4) .or. (nu == 4)) then
      Bmunu(mu,nu,p5)=0d0
      Bmunu(mu,nu,p6)=0d0
      else
      Bmunu(mu,nu,p5)=B5(mu,nu)
     & -p(p5,mu)*B5(4,nu)/p(p5,4)
     & -p(p5,nu)*B5(mu,4)/p(p5,4)
     & +p(p5,mu)*p(p5,nu)*B5(4,4)/p(p5,4)**2
      Bmunu(mu,nu,p6)=B6(mu,nu)
     & -p(p6,mu)*B6(4,nu)/p(p6,4)
     & -p(p6,nu)*B6(mu,4)/p(p6,4)
     & +p(p6,mu)*p(p6,nu)*B6(4,4)/p(p6,4)**2
      endif
      enddo
      enddo

      return
      end 
