      subroutine setreal(pin,rflav,amp2real)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'qcdcouple.f'
      include 'constants.f'
      real * 8 pin(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2real
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      integer i

c --- set scale dependent QCD coupling 
      gsq = st_alpha*fourpi
      ason2pi = st_alpha/twopi

      do i=1,nlegreal
         p(i,4) = pin(0,i)
         p(i,1:3) = pin(1:3,i)
      enddo

      p(1,:)=-p(1,:)
      p(2,:)=-p(2,:)

      call qqb_ww_g(p,msq)

      amp2real = msq(rflav(1),rflav(2))
      amp2real = amp2real/ason2pi

      end
