      double complex function Vpole(sij)
      implicit none
c---  DKS Eq. 2.12
!      include 'epinv.f'
!      include 'epinv2.f'
!      include 'scale.f'
      include '../include/pwhg_st.h'
      !TM added these
      double precision musq,epinv,epinv2
      double precision sij
      double complex Lnrat,xl12
        
      !TM set variables
      musq = st_muren2
      epinv = 0d0
      epinv2 = 0d0

      xl12=Lnrat(-sij,musq)

      Vpole=!-epinv*epinv2+epinv*(-1.5d0+xl12)
     .   -0.5d0*xl12**2+1.5d0*xl12-3.5d0

      return
      end
