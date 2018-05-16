      subroutine w2jetsqlc(p,msq1,msq2)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex qcd1(-1:1,-1:1),qcd2(-1:1,-1:1)
      double precision p(mxpart,4),msq1,msq2
      integer i1,i2,i3,i4,i5,i6
      call spinoru(6,p,za,zb)
      call subqcd(1,2,3,4,5,6,za,zb,qcd1)
      call subqcd(1,2,3,4,6,5,za,zb,qcd2)

      msq1= abs(qcd1(+1,+1))**2+abs(qcd1(+1,-1))**2
     .     +abs(qcd1(-1,+1))**2+abs(qcd1(-1,-1))**2

      msq2= abs(qcd2(+1,+1))**2+abs(qcd2(+1,-1))**2
     .     +abs(qcd2(-1,+1))**2+abs(qcd2(-1,-1))**2

      return
      end
