
c-------------------------------------------------------------------

      subroutine calc_invariants!(p,v,q12,q34,Lmax)
      implicit none
c
c determine the momentum transfers of the scattering weak bosons 
c and define values for the factorization scales and alphas in vbfnlo format
c
c INPUT  p,v         external particle momenta in standard notation
c        Lmax        number of momentum configurations, L=1 for direct term
c                    L=2,...,Lmax for ptilde momenta for subtraction terms
c                    Lmax not in [1,3] prints initialization information
c
c OUTPUT q12(0:4,L)  momentum of weak boson attached to upper fermion line
c        q34(0:4,L)  momentum of weak boson attached to lower fermion line
c                    fourth component is abs(qij^2)
c
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'scales.inc'
      include 'vbfnlo-files/global.inc'
c
      real*8 p(0:3,np,3), v(0:3,nv), q12(0:4,3), q34(0:4,3), qsq
      integer Lmax, L, mu
c
cccccccccccccccccccccccccccccccccccccccccccccc
c
      qsq = st_muren2
      do L = 1,3
         als(1,L) =  st_alpha
         als(2,L) =  st_alpha       
         mursq(1,L) = qsq
	 mursq(2,L) = qsq
      enddo         


      end
