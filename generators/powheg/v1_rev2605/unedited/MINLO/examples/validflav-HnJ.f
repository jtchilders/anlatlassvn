
c this should be a process dependent function
      function validflav(lflav)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      logical validflav
      integer lflav(nlegborn)
      integer onem
      parameter (onem=1000000)
      integer j
c the only invalid flavours is q qbar ->H with no other partons
      if(lflav(1).ne.0.or.lflav(2).ne.0) then
         do j=3,nlegborn
            if(abs(lflav(j)).le.st_nlight) then
               validflav = .true.
               return
            endif
         enddo
         validflav = .false.
         return
      endif
      validflav = .true.
      end

