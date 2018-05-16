      subroutine bornzerodamp(alr,r0,rc,rs,dampfac)
c given the R_alpha region (i.e. the alr) and the associated
c real contribution r (without pdf factor),
c returns in dampfac the damping factor to be applied to
c the real contribution to implement Born zero suppression
      implicit none
      include 'pwhg_flg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      integer alr
      real * 8 r0,rc,rs,dampfac
      real * 8 m2eg,k2ga
      real * 8 dotp
      external dotp
c In case of photon emission, we need to damp the region where the
c electron is collinear to the photon.
      if( flg_withdamp .and. flst_alr(5,alr).eq.0 ) then
         m2eg = 2*dotp(kn_preal(0,3),kn_preal(0,5))+kn_masses(3)**2
         k2ga = kn_preal(1,5)**2+kn_preal(2,5)**2
         dampfac=m2eg/(k2ga+m2eg)
      else
         dampfac=1
      endif
      end
