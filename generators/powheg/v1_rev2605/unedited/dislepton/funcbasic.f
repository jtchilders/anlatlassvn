c     some simple auxiliary functions for matrix elements
c     2012-04 AvM


c     Kronecker delta
      real *8 function kronecker(i,j)
      implicit none
      integer i,j
      if (i.eq.j) then
        kronecker = 1d0
      else
        kronecker = 0d0
      endif
      return
      end



c     exchanges two integers
      subroutine swapi(i1,i2)
      implicit none
      integer i1,i2,itmp
      itmp=i1
      i1=i2
      i2=itmp
      end



c     scalar product of two four vectors (Bjorken-Drell metric)
      real *8 function sprod(p,k)
      implicit none
      real *8 p(0:3),k(0:3)
      sprod = (p(0)*k(0) - p(3)*k(3)) - p(1)*k(1) - p(2)*k(2)
      end



c     norm squared of a four vector
      real *8 function momnormsq(p)
      implicit none
      real *8 p(0:3), sprod
      momnormsq = sprod(p,p)
      end



c     norm squared of sum of 2 four vectors
      real *8 function momsum2normsq(p1,p2)
      implicit none
      real *8 p1(0:3),p2(0:3),psum(0:3)
      real *8 sprod
      integer j
      do j=0,3
        psum(j) = p1(j) + p2(j)
      enddo
      momsum2normsq = sprod(psum,psum)
      end



c     norm squared of sum of 3 four vectors
      real *8 function momsum3normsq(p1,p2,p3)
      implicit none
      real *8 p1(0:3),p2(0:3),p3(0:3),psum(0:3)
      real *8 sprod
      integer j
      do j=0,3
        psum(j) = p1(j) + p2(j) + p3(j)
      enddo
      momsum3normsq = sprod(psum,psum)
      end
