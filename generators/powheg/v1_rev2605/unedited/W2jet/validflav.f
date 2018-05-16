c this should be a process dependent function
      function validflav(lflav)
      implicit none
      include 'nlegborn.h'
      logical validflav
      integer lflav(nlegborn)
c Up to Z+1J, only gg initial states are invalid.
c Since we arrive here always after one clustering, this
c is sufficient up to Z+2j
      if(lflav(1).eq.0.and.lflav(2).eq.0) then
         validflav = .false.
      else
         validflav = .true.
      endif
      end

