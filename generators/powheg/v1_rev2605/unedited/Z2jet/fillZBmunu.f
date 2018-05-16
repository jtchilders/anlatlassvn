      subroutine fillZBmunu(p1,p2,p3,p4,p5,p6,p,Bmunu)
      implicit none
      include 'constants.f'
      integer mu,nu,p1,p2,p3,p4,p5,p6,pq,pl
      double precision Bmunu(4,4,6,2,2),B5(4,4,2,2),B6(4,4,2,2),
     & p(mxpart,4)
      double complex fab(4,2,2,2),fba(4,2,2,2),fqed(4,2,2,2)
      Bmunu(:,:,:,:,:)=0d0
      call Z2jetampmu(p1,p2,p3,p4,p5,p6,fab,fba,fqed)
      do mu=1,4
      do nu=1,4
      do pq=1,2
      do pl=1,2
      B5(mu,nu,pq,pl)=
     & +dble(fab(mu,1,pq,pl)*dconjg(fab(nu,1,pq,pl))
     &      +fab(mu,2,pq,pl)*dconjg(fab(nu,2,pq,pl)))
     & +dble(fba(mu,1,pq,pl)*dconjg(fba(nu,1,pq,pl))
     &      +fba(mu,2,pq,pl)*dconjg(fba(nu,2,pq,pl)))
     & -dble(fqed(mu,1,pq,pl)*dconjg(fqed(nu,1,pq,pl))
     &      +fqed(mu,2,pq,pl)*dconjg(fqed(nu,2,pq,pl)))/xnsq
      enddo
      enddo
      enddo
      enddo
      call Z2jetampmu(p1,p2,p3,p4,p6,p5,fab,fba,fqed)
      do mu=1,4
      do nu=1,4
      do pq=1,2
      do pl=1,2
      B6(mu,nu,pq,pl)=
     & +dble(fab(mu,1,pq,pl)*dconjg(fab(nu,1,pq,pl))
     &      +fab(mu,2,pq,pl)*dconjg(fab(nu,2,pq,pl)))
     & +dble(fba(mu,1,pq,pl)*dconjg(fba(nu,1,pq,pl))
     &      +fba(mu,2,pq,pl)*dconjg(fba(nu,2,pq,pl)))
     & -dble(fqed(mu,1,pq,pl)*dconjg(fqed(nu,1,pq,pl))
     &      +fqed(mu,2,pq,pl)*dconjg(fqed(nu,2,pq,pl)))/xnsq
      enddo
      enddo
      enddo
      enddo

      do mu=1,4
      do nu=1,4
      do pq=1,2
      do pl=1,2
      Bmunu(mu,nu,p5,pq,pl)=B5(mu,nu,pq,pl)
      Bmunu(mu,nu,p6,pq,pl)=B6(mu,nu,pq,pl)
      enddo
      enddo
      enddo
      enddo


C--optional gauge choice A4=0
      do mu=1,4
      do nu=1,4
      do pq=1,2
      do pl=1,2
      if ((mu == 4) .or. (nu == 4)) then
       Bmunu(mu,nu,p5,pq,pl)=0d0
       Bmunu(mu,nu,p6,pq,pl)=0d0
      else
       Bmunu(mu,nu,p5,pq,pl)=B5(mu,nu,pq,pl)
     & -p(p5,mu)*B5(4,nu,pq,pl)/p(p5,4)
     & -p(p5,nu)*B5(mu,4,pq,pl)/p(p5,4)
     & +p(p5,mu)*p(p5,nu)*B5(4,4,pq,pl)/p(p5,4)**2
      Bmunu(mu,nu,p6,pq,pl)=B6(mu,nu,pq,pl)
     & -p(p6,mu)*B6(4,nu,pq,pl)/p(p6,4)
     & -p(p6,nu)*B6(mu,4,pq,pl)/p(p6,4)
     & +p(p6,mu)*p(p6,nu)*B6(4,4,pq,pl)/p(p6,4)**2
      endif
      enddo
      enddo
      enddo
      enddo
C--optional gauge choice A4=0

      Bmunu(:,:,:,:,:)=8d0*Bmunu(:,:,:,:,:)
      return
      end 
