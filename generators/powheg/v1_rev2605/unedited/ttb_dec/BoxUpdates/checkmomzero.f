      subroutine checkmass(n,p)
      implicit none
      integer n
      real * 8 p(0:3,n)
      real * 8 s(0:3)
      integer mu,j
      do mu=0,3
         s(mu)=0
      enddo
      do j=1,2
         do mu=0,3
            s(mu)=s(mu)+p(mu,j)
         enddo
      enddo
      do j=3,n
         do mu=0,3
            s(mu)=s(mu)-p(mu,j)
         enddo
      enddo
      write(*,*) 'mass',sqrt(s(0)**2-s(1)**2-s(2)**2-s(3)**2)
      end

 
      subroutine checkmomzero(n,p)
      implicit none
      integer n
      real * 8 p(0:3,n)
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 s(0:3),r(0:3)
      real * 8 ep
      parameter (ep=1d-10)
      integer j,k,ires
      do j=2,n
         if(j.eq.2.or.flst_isres(j)) then
            if(j.eq.2) then
               s=p(:,1)+p(:,2)
               r=s
               ires=0
            else
               ires=j
               s=p(:,j)
               r=p(:,j)
            endif
            do k=3,nlegreal-1
               if(flst_bornres(k,1).eq.ires) then
                  s=s-p(:,k)
                  r=r+p(:,k)
               endif
            enddo
            if(flst_bornres(kn_emitter,1).eq.ires) then
               s=s-p(:,nlegreal)
               r=r+p(:,nlegreal)
            endif
            if(s(0)**2+s(1)**2+s(2)**2+s(3)**2.ne.0d0
     1           .and. (s(0)**2+s(1)**2+s(2)**2+s(3)**2)
     2           /(r(0)**2+r(1)**2+r(2)**2+r(3)**2).gt.ep) then
               write(*,*) ' momentum check not working',s
               write(*,*) ' for ',j,' decay products'
            endif
         endif
      enddo
      end

        
