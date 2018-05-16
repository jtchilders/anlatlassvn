      logical function isodd(k)
      implicit none
      integer k
      isodd=((k/2)*2.ne.k)
      end
