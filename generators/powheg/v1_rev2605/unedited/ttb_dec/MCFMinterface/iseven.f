      logical function iseven(k)
      implicit none
      integer k
      iseven=((k/2)*2.eq.k)
      end
