      double complex function A0Hggggpmmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'MCFM_Include/constants.f'
      include 'MCFM_Include/zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiggggpmmm
      A0Hggggpmmm=A0phiggggpmmm(j1,j2,j3,j4,za,zb)
c     .           +A0phiggggmppp(j1,j2,j3,j4,zb,za) ! This term is zero
      return
      end
