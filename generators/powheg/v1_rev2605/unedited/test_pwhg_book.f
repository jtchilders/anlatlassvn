      implicit none
      real * 8 x,random
      integer j
      call pwhginihist
      call pwhgbookup(2,'test','LIN',1d-2,0d0,1d0)
      do j=1,10000
         x=random()
         call pwhgfill(2,x,x+100)
         call pwhgfill(2,x,-100d0)
         call pwhgfill(2,x,20d0)
         call pwhgfill(2,x,-20d0)
         call pwhgaccumup
      enddo
      call pwhgaddout
      open(unit=99,file='1.top',status='unknown')
      call pwhgtopout
      close(99)
      do j=1,10000
         x=random()
         call pwhgfill(2,x,x+100)
         call pwhgfill(2,x,-100d0)
         call pwhgfill(2,x,20d0)
         call pwhgfill(2,x,-20d0)
         call pwhgaccumup
      enddo
      call pwhgaddout
      open(unit=99,file='2.top',status='unknown')
      call pwhgtopout
      close(99)
      end
