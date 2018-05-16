c---- scalar product for the two spinors
      subroutine psp(sp1,sp2,r)
      implicit none
      integer i
      double complex sp1(4),sp2(4)
      double complex r 
      r=dcmplx(0d0,0d0)
      do i=1,4
      r = r + sp1(i)*sp2(i)
      enddo
      return
      end
