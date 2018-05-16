c     slepton pair production
c     real emission contributions at NLO SQCD:
c     parton parton -> slepton slepton + parton
c     2012-04 Andreas v. Manteuffel

c      -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6
c      t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t

c     neutral current slepton pairs + jet
c     q(1) qbar(2) -> Z/gamma -> slepton^-(3) slepton^+(4) gluon(5)
c     and crossed versions
c     average over incoming helicites
      subroutine setreal(p,rflav,amp2real)
      implicit none
#include "PhysPars.h"
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      real *8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real *8 amp2real
      real *8 pin(0:3,5), pd(0:3)
      integer i1, i2, i5, k, NA, xsign
      real *8 s125, r15, r25, rd, r1d, r2d, r5d
      real *8 gs2, delta_ij, amp2, pre
      real *8 e_q, t3_q, gl_q, gr_q
      real *8 e_l, t3_l, gl_l, gr_l
      complex *16 cgl_sl, cgr_sl
      real *8 aapre,aalr,aarl,zapre,zalr,zarl,zzpre,zzlr,zzrl
      real *8 propfac_a, averagefac
      complex *16 propfac_z
      real *8 kronecker
      real *8 momnormsq, momsum2normsq, momsum3normsq
      real *8 alemeff, e4
      logical verbose
c      parameter (verbose=.true.)
      logical init_real
      data init_real/.true./

      if (init_real) then
         verbose = .true.
         init_real = .false.
      else
         verbose = .false.
      endif

      ! prepare momenta for crossings
      do k=0,3
        pin(k,1) =  p(k,1)
        pin(k,2) =  p(k,2)
        pin(k,5) = -p(k,5)
        pd(k) = (p(k,3) - p(k,4))
      enddo
      ! set up some constants
      NA = NC**2 - 1
      gs2 = 4*pi*st_alpha

      ! spin averaged square is pquark <-> pantiquark symmetric,
      ! need only to identify gluon momentum
      ! default numbering is quark(1) + antiquark(2) -> gluon(5) + X
      i1 = 1
      i2 = 2
      i5 = 5
      if (rflav(5).eq.0) then ! q q -> sl sl g
        xsign = 1
        averagefac =  (1d0/2)**2 * (1d0/NC)**2
      else if (rflav(1).eq.0) then ! g q -> sl sl q
        call swapi(i5, i1)
        xsign = -1
        averagefac =  (1d0/2)**2 * (1d0/NC) * (1d0/NA)
      else if (rflav(2).eq.0) then ! q g -> sl sl q
        call swapi(i5, i2)
        xsign = -1
        averagefac =  (1d0/2)**2 * (1d0/NC) * (1d0/NA)
      else
        write(*,*) "found real flavours ", rflav
        write(*,*) "this is probably an error, check init_processes"
        stop
      endif

      s125 = momsum3normsq(pin(0,i1),pin(0,i2),pin(0,i5))
      if (verbose) write(*,*) "s125=", s125
      rd = momnormsq(pd)/s125
      r15 = momsum2normsq(pin(0,i1),pin(0,i5))/s125
      r25 = momsum2normsq(pin(0,i2),pin(0,i5))/s125
      r1d = momsum2normsq(pin(0,i1),pd)/s125
      r2d = momsum2normsq(pin(0,i2),pd)/s125
      r5d = momsum2normsq(pin(0,i5),pd)/s125

      if ((abs(rflav(i1)).ge.1).and.(abs(rflav(i1)).le.6)
     &    .and.(rflav(i1).eq.(-rflav(i2)*xsign))) then
        ! photon and Z propagator factors
        propfac_a  = 1d0
        propfac_z=dcmplx(s125)/dcmplx(s125-ph_Zmass2,ph_Zmass*ph_Zwidth)
        ! quark couplings
        if (mod(abs(rflav(i1)),2).eq.0) then
          ! up-type quark
          if (verbose) write(*,*) "up-type quark selected"
          e_q = 2d0 / 3
          t3_q = 1d0 / 2
        else
          ! down-type quark
          if (verbose) write(*,*) "down-type quark selected"
          e_q = -1d0 / 3
          t3_q = -1d0 / 2
        endif
        gl_q = -(t3_q - e_q*ph_sthw2)/(ph_sthw*ph_cthw)
        gr_q = -(0d0  - e_q*ph_sthw2)/(ph_sthw*ph_cthw)
        ! lepton couplings (charged lepton l-)
        e_l  = -1d0
        t3_l = -1d0 / 2
        gl_l = -(t3_l - e_l*ph_sthw2)/(ph_sthw*ph_cthw)
        gr_l = -(0d0  - e_l*ph_sthw2)/(ph_sthw*ph_cthw)
        ! slepton couplings
        delta_ij = kronecker(abs(rflav(3)),abs(rflav(4)))
        cgl_sl = gl_l*ph_c_U3l*conjg(ph_c_U4l)
        cgr_sl = gr_l*ph_c_U3r*conjg(ph_c_U4r)
        ! |gamma|^2       (for specific q hel)
        aapre = delta_ij *e_q**2 *e_l**2 *propfac_a**2
        aalr = aapre
        aarl = aapre
        ! 2 Re(Z gamma^*) (for specific q hel)
        zapre = 2d0*dreal(delta_ij*e_q*e_l*propfac_a*propfac_z
     &          *(cgl_sl+cgr_sl))
        zalr = zapre*gl_q
        zarl = zapre*gr_q
        ! |Z|^2           (for specific q hel)
        zzpre = abs(propfac_z)**2 *(abs(cgl_sl+cgr_sl))**2
        zzlr = zzpre*gl_q**2
        zzrl = zzpre*gr_q**2

        ! electromagnetic unit charge to power 4
        e4 = 16d0 * pi**2 * alemeff(s125)**2

        pre = (NA/2d0)
     &    * (1d0/s125) * ( 2* gs2 * e4 / r15 / r25 )
     &    * (
     &       - 2*rd + 4*rd**2 - 2*r5d*rd - 3*r2d*rd + r2d*r5d
     &       - 3*r1d*rd + r1d*r5d + 2*r1d*r2d + 2*r25*rd - 3*r25*rd**2
     &       + r25*r5d*rd + r25*r2d*rd + 4*r25*r1d*rd - r25*r1d*r5d
     &       - r25*r1d*r2d - r25*r1d**2 - r25**2*rd + 2*r15*rd
     &       - 3*r15*rd**2 + r15*r5d*rd + 4*r15*r2d*rd - r15*r2d*r5d
     &       - r15*r2d**2 + r15*r1d*rd - r15*r1d*r2d - r15**2*rd
     &    )

        ! averaged |M|^2
        amp2 = xsign*averagefac*pre*(aalr + zalr + zzlr
     &                             + aarl + zarl + zzrl)
      else
        amp2 = 0d0
      endif

      amp2real = amp2/(st_alpha/(2*pi))

      end
