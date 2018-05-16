      subroutine setreal(pin,rflav,amp2real)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'qcdcouple.f'
c     vector boson id and decay
      include 'cvecbos.h'
      real * 8 pin(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2real

      integer i,mxpart
      parameter (mxpart=12)
      double precision p(mxpart,4),msq(-5:5,-5:5)
      real * 8 suppfact4e
      external suppfact4e
      ason2pi = st_alpha/2d0/pi

      do i=1,nlegreal
         p(i,4) = pin(0,i)
         p(i,1:3) = pin(1:3,i)
      enddo

      p(1,:)=-p(1,:)
      p(2,:)=-p(2,:)


      call qqb_zz_g(p,msq)


      amp2real = msq(rflav(1),rflav(2))
      amp2real = amp2real/ason2pi

c phase space suppression of (36)(45) singularities
      amp2real = suppfact4e(pin,rflav) * amp2real

      end
