

      subroutine matchcolour(n,flav,colour)
      implicit none
      integer n,flav(4),colour(2,4)
      integer j,k,is,iwhere
c Check that color assignment is compatible with flavour.
c If not, conjugate the colour in the whole colour chain.
c This can be done consistently if the two ends of the colour chain
c have opposite flavour.
c We assumes that the colour array is a self consistent colour assignment.

c first conjugate initial flavours and colours, in order to handle
c all particles as outgoing.
      flav(1)=-flav(1)
      flav(2)=-flav(2)
      call colour_conj(colour(:,1))
      call colour_conj(colour(:,2))

      do k=1,n
         if(flav(k).eq.0) then
            if(colour(1,k).eq.0.or.colour(2,k).eq.0) goto 999
         elseif(abs(flav(k)).le.6) then
            if(flav(k).gt.0) then
               is = 1
            else
               is = 2
            endif
c The following is normal assignment
            if(colour(is,k).gt.0.and.colour(3-is,k).eq.0) cycle
c the following is impossible to handle:
            if(colour(1,k).eq.0.and.colour(2,k).eq.0) goto 999
            if(colour(1,k).ne.0.and.colour(2,k).ne.0) goto 999
c the remaining possibility is colour(is,k)=0 and colour(3-is,k)!=0;
c Correct by exchanging colour and anticolour in the whole colour connected
c chain, provided the last element has opposite sign of flavour
            iwhere = k
 10         continue
            do j=1,n
               if(j.ne.iwhere) then
                  if(colour(is,j).eq.colour(3-is,iwhere)) then
                     if(flav(j).eq.0) then
                        if(colour(1,j).ne.0.and.colour(2,j).ne.0) then
                           call colour_conj(colour(:,iwhere))
                           iwhere = j
                           goto 10
                        else
c this cannot be handled
                           goto 999
                        endif
                     elseif(flav(j)*flav(k).lt.0) then
c found a connected opposite flavour;
c The following can't be, colour(is,j)!=0;
                        if(colour(3-is,j).ne.0) goto 999
                        call colour_conj(colour(:,iwhere))
                        call colour_conj(colour(:,j))
                        exit
                     else
                        goto 999
                     endif
                  endif
               endif
            enddo
         endif
      enddo

c Straighten up colour and flavour
      flav(1)=-flav(1)
      flav(2)=-flav(2)
      call colour_conj(colour(:,1))
      call colour_conj(colour(:,2))

      return

 999  continue
      write(*,*)
     1     ' matchcolour: incompatible colour-flavour configuration',
     2     ' better to switch off smartsig in sigborn.f ...'
      call pwhg_exit(-1)
      end
