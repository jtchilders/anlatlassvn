      double complex function Vpole(sij)
      implicit none
c---  DKS Eq. 2.12

      include '../include/pwhg_st.h'
      double precision musq,epinv,epinv2
      double precision sij
      double complex Lnrat,xl12
        

      musq = st_muren2
      epinv=0d0
      epinv2=0d0
      
      epinv=1d0
      epinv2=1d0
      xl12=Lnrat(-sij,musq)

c      Vpole=-epinv*epinv2+epinv*(-1.5d0+xl12)
c     .   -0.5d0*xl12**2+1.5d0*xl12-3.5d0
      
      Vpole = -0.5d0*xl12**2+1.5d0*xl12-3.5d0


      return
      end
