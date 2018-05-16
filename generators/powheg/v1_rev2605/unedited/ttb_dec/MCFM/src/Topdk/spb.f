c-----multiplication of a bar spinor with vslash from the right
      subroutine spb(sp,v,f) 
      implicit none
      double complex sp(4),v(4),f(4),im
      parameter(im=(0d0,1d0))
      f(1)=sp(3)*(v(1)-v(4))-sp(4)*(v(2)+im*v(3))
      f(2)=sp(4)*(v(4)+v(1))-sp(3)*(v(2)-im*v(3))
      f(3)=sp(1)*(v(4)+v(1))+sp(2)*(v(2)+im*v(3))
      f(4)=sp(2)*(v(1)-v(4))+sp(1)*(v(2)-im*v(3))
      return
      end
