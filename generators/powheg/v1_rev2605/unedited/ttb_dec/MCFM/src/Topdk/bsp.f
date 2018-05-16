c-----multiplication of a spinor with vslash from the left
      subroutine bsp(v,sp,f) 
      implicit none
      double complex sp(4),v(4),f(4),im
      parameter(im=(0d0,1d0))
      f(1)=(v(1)+v(4))*sp(3)+(v(2)-im*v(3))*sp(4)
      f(2)=(v(1)-v(4))*sp(4)+(v(2)+im*v(3))*sp(3)
      f(3)=(v(1)-v(4))*sp(1)-(v(2)-im*v(3))*sp(2)
      f(4)=(v(1)+v(4))*sp(2)-(v(2)+im*v(3))*sp(1)
      return
      end
