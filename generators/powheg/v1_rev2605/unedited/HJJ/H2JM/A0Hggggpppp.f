      double complex function A0Hggggpppp(j1,j2,j3,j4,za,zb)
      include 'MCFM_Include/constants.f'
      include 'MCFM_Include/zprods_decl.f'
      integer j1,j2,j3,j4
      double complex A0phiggggpppp,A0phiggggmmmm
      A0Hggggpppp=A0phiggggpppp(j1,j2,j3,j4,za,zb)
     .           +A0phiggggmmmm(j1,j2,j3,j4,zb,za)
      return
      end
