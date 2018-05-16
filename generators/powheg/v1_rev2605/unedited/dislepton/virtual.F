c     slepton pair production
c     virtual corrections at NLO SQCD:
c     parton parton -> slepton slepton
c     (effectively just Drell-Yan type corrections to q qbar -> gamma/Z)
c     2012-04 Andreas v. Manteuffel


c     pick SUSY masses relevant for specific initial state
      subroutine set_squark_params(quarkid, sqlmass, sqrmass)
      implicit none
#include "PhysPars.h"
      integer quarkid
      real *8 sqlmass, sqrmass

      select case(abs(quarkid))
      case (1) ! d
        sqlmass = real(Mass_MSf(1,4,1))
        sqrmass = real(Mass_MSf(2,4,1))
      case (2) ! u
        sqlmass = real(Mass_MSf(1,3,1))
        sqrmass = real(Mass_MSf(2,3,1))
      case (3) ! s
        sqlmass = real(Mass_MSf(1,4,2))
        sqrmass = real(Mass_MSf(2,4,2))
      case (4) ! c
        sqlmass = real(Mass_MSf(1,3,2))
        sqrmass = real(Mass_MSf(2,3,2))
      case (5) ! b
        sqlmass = real(Mass_MSf(1,4,3))
        sqrmass = real(Mass_MSf(2,4,3))
      case default
        write(*,*) 'encountered unhandled incoming quark ID ', quarkid
        stop
      end select
      end



c     returns 2 Re(M_B * M_V)/(as/(2pi)),
c     where M_B is the Born amplitude and
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
c     A factor (4 Pi)^eps/Gamma[1-eps] with eps=(4-d)/2 is implied.
      subroutine setvirtual(p,vflav,twicevirtual)
      implicit none
      include 'nlegborn.h'
      real *8 p(0:3,nlegborn), twicevirtual
      integer vflav(nlegborn)
      complex *16 cv

      call compvirtual(p,vflav,cv)
      twicevirtual = 2d0 * dreal(cv)

      end



c     returns (M_B * M_V)/(as/(2pi)) as complex number,
c     where M_B is the Born amplitude and
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
c     A factor (4 Pi)^ep/Gamma(1-ep) with ep=(4-d)/2 is implied.
      subroutine compvirtual(p,vflav,virtual)
      implicit none
#include "PhysPars.h"
      include 'pwhg_st.h'    ! renormalization scale
      include 'pwhg_math.h'  ! pi
      include 'nlegborn.h'
      real *8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      complex *16 virtual
      real *8 msql, msqr, mgl
      real *8 m2sql, m2sqr, m2gl
      real *8 lmu,s,sprod,bornlr,bornrl
      real *8 rgl, rsq1l, rsq1r, musq
      real *8 rglMrsq1l, rglMrsq1r
      complex *16 ii, vlr, vrl, pre
      complex *16 I1fin, I2fin, I3fin
      real *8 O, z2
      logical verbose
      parameter (O=0d0,z2=pi**2/6d0)
      integer hSM,hSUSY,hWF,hAltNorm
      ! flags to select different contributions:
      ! hSUSY = 0 or 1: omit/include SUSY contributions
      ! hWF   = 0 or 1: omit/include wave function contrib
      ! (setting all to 1 is default and includes everything)
      ! debug flag to switch to alternative normalization factor:
      ! hAltNorm = 0 or 1: 0 for default as in subroutine description
      !                    1 for (4 pi)^ep * exp(-gamma ep) prefactor
      parameter (hSM=1,hSUSY=1,hWF=1,hAltNorm=0,verbose=.true.)
      logical firsttime
      data firsttime/.true./
      save firsttime

      ! initialise QCDLoop library for scalar loop integrals
      if (firsttime) then
        call qlinit
        firsttime = .false.
      endif

      ! squarks for incoming massless flavors assumed non-mixing
      ! => Born factorises from hel. amp.
      call compbornhel(p,vflav,bornlr,bornrl)

      if (vflav(1).ne.(-vflav(2)))  stop "invalid initial state"
      call set_squark_params(vflav(1), msql, msqr)
      mgl = real(Mass_MGl)
      m2gl = mgl**2
      m2sql = msql**2
      m2sqr = msqr**2

      ii = dcmplx(0d0, 1d0)
      s = 2d0 * sprod(p(0,1),p(0,2))
      lmu = log(st_muren2/s)
      musq  = st_muren2
      rgl = m2gl / s
      rsq1l = m2sql / s
      rsq1r = m2sqr / s
      rglMrsq1l = rgl - rsq1l
      rglMrsq1r = rgl - rsq1r

      vlr =  + hSUSY * ( 1.D0/2.D0 - 1.D0/2.D0*hWF*rglMrsq1l**(-1)*rgl
     &     + 1.D0/4.D0*hWF )
      vlr = vlr + hSM * (  - 4.D0 - 3.D0/2.D0*ii*pi - 3.D0/2.D0*Lmu - 
     &    Lmu*ii*pi - 1.D0/2.D0*Lmu**2 + 3.D0*z2 + 1.D0/2.D0*hAltNorm*
     &    z2 )
      vlr = vlr + I1fin(m2gl,musq)*hSUSY*s**(-1) * (
     &     - 1.D0 - 1.D0/2.D0*hWF*rglMrsq1l**(-2)*rgl )
      vlr = vlr + I1fin(m2sql,musq)*hSUSY*s**(-1) * ( 1.D0
     &    + 1.D0/2.D0*hWF*rglMrsq1l**(-2)*rgl
     &    + 1.D0/2.D0*hWF*rglMrsq1l**(-1) )
      vlr = vlr + I2fin(s,m2sql,m2sql,musq)*hSUSY * (
     &    1.D0/2.D0 - rsq1l + rgl )
      vlr = vlr + I3fin(O,O,s,m2sql,m2gl,m2sql,musq)*hSUSY*
     &    s * ( rsq1l**2 + rgl - 2.D0*rgl*rsq1l + rgl**2 )

      vrl =  + hSUSY * ( 1.D0/2.D0 - 1.D0/2.D0*hWF*rglMrsq1r**(-1)*rgl
     &     + 1.D0/4.D0*hWF )
      vrl = vrl + hSM * (  - 4.D0 - 3.D0/2.D0*ii*pi - 3.D0/2.D0*Lmu - 
     &    Lmu*ii*pi - 1.D0/2.D0*Lmu**2 + 3.D0*z2 + 1.D0/2.D0*hAltNorm*
     &    z2 )
      vrl = vrl + I1fin(m2gl,musq)*hSUSY*s**(-1) * (  - 1.D0
     &    - 1.D0/2.D0*hWF*rglMrsq1r**(-2)*rgl )
      vrl = vrl + I1fin(m2sqr,musq)*hSUSY*s**(-1) * ( 1.D0
     &    + 1.D0/2.D0*hWF*rglMrsq1r**(-2)*rgl 
     &    + 1.D0/2.D0*hWF*rglMrsq1r**(-1) )
      vrl = vrl + I2fin(s,m2sqr,m2sqr,musq)*hSUSY * (
     &    1.D0/2.D0 - rsq1r + rgl )
      vrl = vrl + I3fin(O,O,s,m2sqr,m2gl,m2sqr,musq)*hSUSY*
     &    s * ( rsq1r**2 + rgl - 2.D0*rgl*rsq1r + rgl**2 )

      pre =  + CF

      virtual = pre *(1d0/2)**2 *(bornlr*vlr + bornrl*vrl)

      end
