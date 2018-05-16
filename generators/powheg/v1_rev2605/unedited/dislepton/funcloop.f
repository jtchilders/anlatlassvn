c     wrappers for loop integrals
c     2012-04 AvM


      complex *16 function I1fin(m1sq,musq)
      implicit none
      complex *16 qlI1
      real *8 m1sq, musq
      I1fin = qlI1(m1sq,musq,0)
      return
      end



      complex *16 function I2fin(p1sq,m1sq,m2sq,musq)
      implicit none
      complex *16 qlI2
      real *8 p1sq, m1sq, m2sq, musq
      I2fin = qlI2(p1sq,m1sq,m2sq,musq,0)
      return
      end



      complex *16 function I3fin(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq)
      implicit none
      complex *16 qlI3
      real *8 p1sq, p2sq, p3sq, m1sq, m2sq, m3sq, musq
      I3fin = qlI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,musq,0)
      return
      end
