c      -*- fortran -*-
ccccccccccccccccccccccccccccccccccc
c     POWHEG COMMON BLOCK (PHYSICAL PARAMETERS)
      common/cewparam/alphaem_pow,gfermi_pow,sthw2_pow,topmass_pow,
     #     mtau_pow,mhiggs_pow,mcharm_pow,mbottom_pow,
     #     zmass_pow,wmass_pow 
      real *8 alphaem_pow,gfermi_pow,sthw2_pow,topmass_pow,
     #     mtau_pow,mhiggs_pow,mcharm_pow,mbottom_pow,
     #     zmass_pow,wmass_pow

      common/cwidth/topwidth_pow,wwidth_pow
      real *8 topwidth_pow,wwidth_pow

      common/cckm/CKM_pow,CKM
      real *8 CKM_pow(3,3),CKM(1:6,1:6)

c     POWHEG COMMON BLOCK (PHASE SPACE GENERATION)
      integer psgen
      real *8 numgamma(2)
      common/cphasespace/numgamma,psgen

c     MADGRAPH COMMON BLOCK
      include 'coupl.inc'

c     T/TBAR
      integer ttype
      common/cttype/ttype
cccccccccccccccccccccccccccccccc  
