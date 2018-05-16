      subroutine setupckmallowed(nwz,diagonal)
      implicit none
      include 'ckmallowed.f'
      integer i,k,nwz
      logical diagonal
c--- initialize array for possible CKM combinations
      do i=-5,5
      do k=-5,5
      ckmallowed(i,k)=.false.
      enddo
      enddo
      
c--- diagonal only
      if (nwz .eq. +1) then
         ckmallowed(2,-1)=.true.
         ckmallowed(-1,2)=.true.
         ckmallowed(-3,4)=.true.
         ckmallowed(4,-3)=.true.
         if (diagonal) return
         ckmallowed(2,-3)=.true.
         ckmallowed(4,-1)=.true.
         ckmallowed(-3,2)=.true.
         ckmallowed(-1,4)=.true.
      elseif (nwz .eq. -1) then
         ckmallowed(-2,1)=.true.
         ckmallowed(-4,3)=.true.
         ckmallowed(1,-2)=.true.
         ckmallowed(3,-4)=.true.
         if (diagonal) return
         ckmallowed(-2,3)=.true.
         ckmallowed(-4,1)=.true.
         ckmallowed(3,-2)=.true.
         ckmallowed(1,-4)=.true.
      else
         write(6,*) 'setupckmallowed:Unallowed value of nwz= ',nwz
         call exit(-1)
      endif
      return
      end
