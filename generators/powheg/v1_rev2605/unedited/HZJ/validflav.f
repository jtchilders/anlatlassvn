
c this should be a process dependent function
      function validflav(lflav)
      implicit none
      include 'nlegborn.h'
      logical validflav
      integer lflav(nlegborn)
c      if(lflav(1).eq.0.and.lflav(2).eq.0) then
c         validflav = .false.
c      else
      validflav = .true.
c      endif
      end
