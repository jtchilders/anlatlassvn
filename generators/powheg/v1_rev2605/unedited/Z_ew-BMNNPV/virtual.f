c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual
      real * 8 virt_st
      complex*16 virt_ew
      real * 8 born,dummy(0:3,0:3)
      real *8 s,dotp
      external dotp
*
      if(vflav(1).eq.22.and.vflav(2).eq.22) then
          virtual = 0d0
          return
      endif
*
      s = 2d0*dotp(p(0,1),p(0,2))
      call compborn(p,vflav,born,dummy)
      virt_st = pi2 - 8 - 3*log(st_muren2/s) - log(st_muren2/s)**2
      virt_st = virt_st * cf * born
*
      call virtual_ew(p,vflav,virt_ew)
*
      virtual = virt_st + 2d0*dble(virt_ew) / 4d0 / nc
     +                                      /(st_alpha/(2d0*pi))
*
      end


c$$$c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c$$$c     where M_B is the Born amplitude and 
c$$$c     M_V is the finite part of the virtual amplitude
c$$$c     The as/(2pi) factor is attached at a later point
c$$$
c$$$c     Like setvirtual, but it calculates only the independent
c$$$c     virtual contributions and fills directly the array 'virt_arr'.
c$$$c     
c$$$c     At the moment (revision 12), no calls to setvirtual_fast are present
c$$$c     in the main code. Therefore, this routine is left dummy.
c$$$
c$$$      subroutine setvirtual_fast(virt_arr)
c$$$      implicit none
c$$$      include 'nlegborn.h'
c$$$      include 'pwhg_flst.h'
c$$$      real * 8 virt_arr(flst_nborn)
c$$$
c$$$      write(*,*) 'Error: setvirtual_fast is not implemented yet'
c$$$      stop
c$$$
c$$$      end

*
** subroutines for the calculation of virtual corrections to dy 
** production at hadron colliders. 
**
** (p1)u/d \            / l+ (p3)
**          \----------/
**          |\________/    (etc)
**          |/   Z/g  \
**          /          \
** (p2)u/d~/            \ l- (p4)
*

      subroutine virtual_ew(p,flav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'mathx.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real* 8 p(0:3,nlegborn)
      integer flav(nlegborn)
      complex*16 virtual
*
      real*8 dotp
      external dotp

      complex*16 flo
      external flo

      real*8 chargeofparticle
      external chargeofparticle
*
      real*8 mlep2
      common/leptmass/mlep2

      complex*16 bornamp(0:1,0:1)
      common/bornamps/bornamp
*
      real*8 s,t,u
      integer sig,tau
      complex*16 self,vert,box 

      real*8 pq(0:3),pqb(0:3),pl(0:3),plb(0:3)
      real*8 ptmpu(0:3),ptmpt(0:3),ptmps(0:3)
      real*8 mqbig2
      integer nu

      qq = chargeofparticle(flav(1))

      if (qq.gt.0d0) then
          gq(1) = gu(1)
          gq(0) = gu(0)
          i3q   = 0.5d0*cone
      elseif (qq.lt.0d0) then
          gq(1) = gd(1)
          gq(0) = gd(0)
          i3q   = -0.5d0*cone
      endif

      if(flav(1).gt.0) then
          do nu=0,3
              pq(nu)  = p(nu,1)
              pqb(nu) = p(nu,2)
          enddo
      else
          do nu=0,3
              pqb(nu) = p(nu,1)
              pq(nu)  = p(nu,2)
          enddo
      endif

      do nu=0,3
          pl(nu)  = p(nu,3)
          plb(nu) = p(nu,4)
      enddo

      if (abs(flav(1)).eq.5) then
          if (iftopinloop) then
              mqbig2 = mt2
          else
              mqbig2 = 0d0
          endif
      else
          mqbig2 = 0d0
      endif
*
** mandelstam invariants
*
      do nu=0,3
          ptmps(nu) = pq(nu) + pqb(nu)
          ptmpt(nu) = pq(nu) - pl(nu)
          ptmpu(nu) = pq(nu) - plb(nu)
      enddo
      s = dotp(ptmps,ptmps)
      t = dotp(ptmpt,ptmpt)
      u = dotp(ptmpu,ptmpu)

      a(1,1) = 2d0*u*cone
      a(0,0) = a(1,1)
      a(0,1) = 2d0*t*cone
      a(1,0) = a(0,1)
      chiz = s/(s*cone - ph_Zmass2 + ii*ph_ZmZw)
*
** calculate born amplitudes
*
      do sig=0,1
          do tau=0,1
              bornamp(sig,tau)=flo(sig,tau,s)
          enddo
      enddo
*
      call deltaself(s,self)

      call deltavert(s,mqbig2,vert)

      call deltabox (s,t,u,mqbig2,box)
*
* eq. (2.10) of Dittmaier-Kraemer PRD65 073007
*
      virtual = self + vert + box

      end subroutine virtual_ew
*
** Born amplitude
*
      complex*16 function flo(sig,tau,s)
      implicit none
      include'pwhg_physpar.h'
*
      integer sig,tau
      real*8 s
*
      flo = - el2_scheme/s * ( qq*ql + gq(sig)*gl(tau)*chiz )
*
      end function flo

*
** Subroutine for the calculation of the 
** contribution to delvirt due to Z/g self-energy
**
**  u~/d \            / l+ 
**        \    f     /
**         \___/\___/    (etc)
**         /   \/   \
**        /    f~    \
**  u/d~ /            \ l-
**
**        s = (p1+p2)**2
*
      subroutine deltaself(s,self)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
      include'PhysPars.h'
*
      real*8 s
      complex*16 self
*
      complex*16 bornamp(0:1,0:1)
      common/bornamps/bornamp
*
      complex*16 szztmz2,szztpmz2,sazts,sazt0
      complex*16 saats,szzts,saztmz2,saatp0

      complex*16 deltamz2,delaa,delzz,delaz,delza
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save deltamz2,delaa,delzz,delaz,delza
*
      integer sig,tau

      if (ifirst.eq.0) then
          ifirst = 1

          call sigmazzt(dble(mz2)*cone,szztmz2)
          call sigmazztp(dble(mz2)*cone,szztpmz2)
          call sigmaazt(zero,sazt0)
          call sigmaazt(dble(mz2)*cone,saztmz2)
          call sigmaaatp(zero,saatp0,.true.)
*
* eq. (3.19) of ArXiv:0709.1075  (Denner Fortschritte) 
* and eq. (4.9) of ArXiv:0505042
*
          if (complexmasses) then
              deltamz2  = szztmz2 + ii*dimag(mz2)*szztpmz2
          else
              deltamz2  = (szztmz2)
          endif
*
* eq. (3.28) of Dittmaier-Huber 0911.2329 and (4.9/10) of hep-ph:0505042
*
          delaa = - saatp0

          if (complexmasses) then
              delzz = - szztpmz2
          else
              delzz = - (szztpmz2)
          endif

          delza =   2d0/mz2*sazt0

          if (complexmasses) then
              delaz = - 2d0/dble(mz2)*saztmz2+(mz2/dble(mz2)-cone)*delza
          else
              delaz = dble(- 2d0/dble(mz2)*saztmz2)
          endif

      endif

      call sigmaaat(s*cone,saats,.true.)
      call sigmazzt(s*cone,szzts)
      call sigmaazt(s*cone,sazts)
*
* eq. (3.27) of Dittmaier-Huber 0911.2329
*
      self = zero
*
      do sig=0,1
          do tau=0,1
          if (complexmasses) then
              self = self 
     +         + el2_scheme * ( qq*ql/s**2 * (saats + delaa*s)
     +                  + gq(sig)*gl(tau)/(s*cone-mz2)**2
     +                    * (szzts-deltamz2+delzz*(s*cone-mz2))
     +                  - (ql*gq(sig)+qq*gl(tau))/(s*(s*cone-mz2))
     +                    * (sazts + 0.5d0*delaz*s
     +                                  + 0.5d0*delza*(s*cone-mz2)) )
     +          *a(sig,tau)*conjg(a(sig,tau))*conjg(bornamp(sig,tau))
          else
              self = self
     +         + el2_scheme * ( qq*ql/s**2 * (saats + delaa*s)
     +                  + gq(sig)*gl(tau)
     +                        /(s*cone-mz2)/(s*cone-mz2+ii*ph_ZmZw)
     +                    * (szzts-deltamz2+delzz*(s*cone-mz2))
     +                  - (ql*gq(sig)+qq*gl(tau))
     +                      /(s*(s*cone-mz2+ii*ph_ZmZw))
     +                    * (sazts + 0.5d0*delaz*s
     +                                  + 0.5d0*delza*(s*cone-mz2)) ) 
     +          *a(sig,tau)*conjg(a(sig,tau))*conjg(bornamp(sig,tau))

          endif
          enddo
      enddo
*
      end subroutine deltaself
*
**
** Subroutine for the calculation of the 
** contribution to deltavirt due to triangles on
** initial and final state state
**
** (pu)  u  \            / l+ 
**           \          /
**          Z|\________/    (etc)
**           |/        \
**           /          \
** (pud) u~ /            \ l-
**
**        s = (pu+pud)**2
*
*
* eq. (3.28) of Dittmaier-Huber 0911.2329
*
      subroutine deltavert(s,mqbig2,vert)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
      include'pwhg_math.h'
      include'PhysPars.h'
*
      real*8 s,mqbig2
      complex*16 vert
*
      complex*16 bornamp(0:1,0:1)
      common/bornamps/bornamp

      real*8 mlep2
      common/leptmass/mlep2
*
      complex*16 getdeltar
      external getdeltar
*
*
      complex*16 sazt0,saztmz2,saatp0,swtmw2
      complex*16 szztmz2,szztpmz2,swtpmw2
      complex*16 saatp0notop,saatmz2notop

      complex*16 delze,delswow,delza,delaz
      complex*16 delaa,delzz,deltamw2,deltamz2

      complex*16 delll_ct_g_alpha0(0:1),delll_ct_z_alpha0(0:1)
      complex*16 deluu_ct_g_alpha0(0:1),deluu_ct_z_alpha0(0:1)
      complex*16 deldd_ct_g_alpha0(0:1),deldd_ct_z_alpha0(0:1)
      complex*16 delbb_ct_g_alpha0(0:1),delbb_ct_z_alpha0(0:1)

      complex*16 deluu_ct_g(0:1),deluu_ct_z(0:1)
      complex*16 deldd_ct_g(0:1),deldd_ct_z(0:1)
      complex*16 delbb_ct_g(0:1),delbb_ct_z(0:1)

      complex*16 delll_ct_g(0:1),delll_ct_z(0:1)
      complex*16 delqq_ct_g(0:1),delqq_ct_z(0:1)
      complex*16 delta
      complex*16 deltatop
*
      complex*16 delguu(0:1),delgdd(0:1),delgbb(0:1)
      complex*16 delzuu(0:1),delzdd(0:1),delzbb(0:1)
      complex*16 delgll(0:1),delgqq(0:1)
      complex*16 delzll(0:1),delzqq(0:1)

      complex*16 fzqqg
      complex*16 fzllg

      complex*16 fzqqw(0:1),fzllw(0:1)
      complex*16 fgqqw(0:1),fgllw(0:1)

      real*8 bornm2
*
      integer ifirst
      data ifirst /0/
      save ifirst

      save delll_ct_g,delll_ct_z,
     +     deluu_ct_g,deluu_ct_z,deldd_ct_g,deldd_ct_z,
     +     delbb_ct_g,delbb_ct_z

      integer sig,tau

      logical calculatec0
      common/ifc0/calculatec0
*
      if (ifirst.eq.0) then
          ifirst = 1

          call sigmaazt(zero,sazt0)
          call sigmaazt(dble(mz2)*cone,saztmz2)
          call sigmaaatp(zero,saatp0,.true.)
          call sigmawt(dble(mw2)*cone,swtmw2)
          call sigmawtp(dble(mw2)*cone,swtpmw2)
          call sigmazzt(dble(mz2)*cone,szztmz2)
          call sigmazztp(dble(mz2)*cone,szztpmz2)

*
* eq. (3.28) of Dittmaier-Huber 0911.2329
*
          delza =   2d0/mz2*sazt0

          if (complexmasses) then
              delaz = - 2d0/mz2*saztmz2 + (mz2/dble(mz2)-cone)*delza
          else
              delaz = (- 2d0/mz2*saztmz2)
          endif

          delaa = - saatp0

          if (complexmasses) then
              delzz = - szztpmz2
          else
              delzz = - szztpmz2
          endif
*
** delze according to eq. (3.32) of ArXiv:0709.1075 (Denner Fortschritte)
*
          delze = 0.5d0*saatp0 - sw/cw * sazt0/mz2
*
* eq. (4.29) of ArXiv:0505042
*
          if (complexmasses) then
              deltamw2 = swtmw2 + ii*dimag(mw2)*swtpmw2
     +                          - ii*dimag(mw2)*alpha/pi

              deltamz2 = szztmz2 + ii*dimag(mz2)*szztpmz2
          else
              deltamw2 = swtmw2

              deltamz2 = szztmz2
          endif

*
* delsw according to eq. (3.35) of ArXiv:0709.1075 (Denner Fortschritte)
* and eq. (4.13) of ArXiv:0505042
*
          if (complexmasses) then
              delswow = -0.5d0*cw2/sw2*( deltamw2/mw2 - deltamz2/mz2 )
          else
              delswow = -0.5d0*cw2/sw2*dble(deltamw2/mw2 - deltamz2/mz2)
          endif
*
* eq. (3.32) of Dittmaier-Huber 
* 
          delgll(1) = -sw/cw*ql*(delze + 1d0/cw2*delswow)

          delgll(0) = i3l/(sw*cw)*(delze + (sw2-cw2)/cw2*delswow)
     +              + delgll(1)

          delguu(1) = -sw/cw*qu*(delze + 1d0/cw2*delswow)

          delguu(0) = 0.5d0/(sw*cw)*(delze + (sw2-cw2)/cw2*delswow)
     +              + delguu(1)

          delgdd(1) = -sw/cw*qd*(delze + 1d0/cw2*delswow)

          delgdd(0) = -0.5d0/(sw*cw)*(delze + (sw2-cw2)/cw2*delswow)
     +              + delgdd(1)
*
*
* delzfl: wave functions renormalization according to eq. (3.20) 
*         of ArXiv:0709.1075 (Denner Fortschritte). only 
*         diagonal terms are considered since ckm matrix = 1 
* here f= up and down quark (if we change final state this should 
*                            be changed accordingly)
*
          call deltazfl(qu*cone,gu(0),gu(1),zero,zero,delzuu(0))
          call deltazfr(qu*cone,gu(0),gu(1),zero,zero,delzuu(1))

          call deltazfl(qd*cone,gd(0),gd(1),zero,zero,delzdd(0))
          call deltazfr(qd*cone,gd(0),gd(1),zero,zero,delzdd(1))

          if (iftopinloop) then
              call deltazfl(qd*cone,gd(0),gd(1),zero,mqbig2*cone,
     +                      delzbb(0))
              call deltazfr(qd*cone,gd(0),gd(1),zero,mqbig2*cone,
     +                      delzbb(1))
          endif

          call deltazfl(ql*cone,gl(0),gl(1),mlep2*cone,zero,delzll(0))
          call deltazfr(ql*cone,gl(0),gl(1),mlep2*cone,zero,delzll(1))

          if (scheme.eq.0) then
              delta = zero
          elseif (scheme.eq.1) then
              call sigmaaat(dble(mz2)*cone,saatmz2notop,.false.)
              call sigmaaatp(zero,saatp0notop,.false.)
              delta = 0.5d0 * (saatp0notop - saatmz2notop/mz2)
          else
              delta = 0.5d0 * getdeltar()
          endif
*
** eq. (3.31) of Dittmaier-Huber 
* 
* l+ (gamma)
          delll_ct_g_alpha0(1) = delze + 0.5d0*delaa + delzll(1)
     +                         - 0.5d0*gl(1)/ql*delza
          delll_ct_g(1) = delll_ct_g_alpha0(1) - delta

* l- (gamma)
          delll_ct_g_alpha0(0) = delze + 0.5d0*delaa + delzll(0)
     +                         - 0.5d0*gl(0)/ql*delza
          delll_ct_g(0) = delll_ct_g_alpha0(0) - delta

* l+ (Z)
          delll_ct_z_alpha0(1) = delgll(1)/gl(1) + 0.5d0*delzz 
     +                         + delzll(1) - 0.5d0*ql/gl(1)*delaz
          delll_ct_z(1) = delll_ct_z_alpha0(1) - delta

* l- (Z)
          delll_ct_z_alpha0(0) = delgll(0)/gl(0) + 0.5d0*delzz 
     +                         + delzll(0) - 0.5d0*ql/gl(0)*delaz
          delll_ct_z(0) = delll_ct_z_alpha0(0) - delta



* u+ (gamma)
          deluu_ct_g_alpha0(1) = delze + 0.5d0*delaa + delzuu(1)
     +                         - 0.5d0*gu(1)/qu*delza
          deluu_ct_g(1) = deluu_ct_g_alpha0(1) - delta

* u- (gamma)
          deluu_ct_g_alpha0(0) = delze + 0.5d0*delaa + delzuu(0)
     +                         - 0.5d0*gu(0)/qu*delza
          deluu_ct_g(0) = deluu_ct_g_alpha0(0) - delta

* u+ (Z)
          deluu_ct_z_alpha0(1) = delguu(1)/gu(1) + 0.5d0*delzz 
     +                         + delzuu(1) - 0.5d0*qu/gu(1)*delaz
          deluu_ct_z(1) = deluu_ct_z_alpha0(1) - delta

* u- (Z)
          deluu_ct_z_alpha0(0) = delguu(0)/gu(0) + 0.5d0*delzz 
     +                         + delzuu(0) - 0.5d0*qu/gu(0)*delaz
          deluu_ct_z(0) = deluu_ct_z_alpha0(0) - delta



* d+ (gamma)
          deldd_ct_g_alpha0(1) = delze + 0.5d0*delaa + delzdd(1)
     +                         - 0.5d0*gd(1)/qd*delza
          deldd_ct_g(1) = deldd_ct_g_alpha0(1) - delta

* d- (gamma)
          deldd_ct_g_alpha0(0) = delze + 0.5d0*delaa + delzdd(0)
     +                         - 0.5d0*gd(0)/qd*delza
          deldd_ct_g(0) = deldd_ct_g_alpha0(0) - delta

* d+ (Z)
          deldd_ct_z_alpha0(1) = delgdd(1)/gd(1) + 0.5d0*delzz 
     +                         + delzdd(1) - 0.5d0*qd/gd(1)*delaz
          deldd_ct_z(1) = deldd_ct_z_alpha0(1) - delta

* d- (Z)
          deldd_ct_z_alpha0(0) = delgdd(0)/gd(0) + 0.5d0*delzz 
     +                         + delzdd(0) - 0.5d0*qd/gd(0)*delaz
          deldd_ct_z(0) = deldd_ct_z_alpha0(0) - delta

          if (iftopinloop) then
* b+ (gamma)
              delbb_ct_g_alpha0(1) = delze + 0.5d0*delaa + delzbb(1)
     +                             - 0.5d0*gd(1)/qd*delza
              delbb_ct_g(1) = delbb_ct_g_alpha0(1) - delta

* b- (gamma)
              delbb_ct_g_alpha0(0) = delze + 0.5d0*delaa + delzbb(0)
     +                             - 0.5d0*gd(0)/qd*delza
              delbb_ct_g(0) = delbb_ct_g_alpha0(0) - delta

* b+ (Z)
              delbb_ct_z_alpha0(1) = delgbb(1)/gd(1) + 0.5d0*delzz 
     +                             + delzbb(1) - 0.5d0*qd/gd(1)*delaz
              delbb_ct_z(1) = delbb_ct_z_alpha0(1) - delta

* b- (Z)
              delbb_ct_z_alpha0(0) = delgbb(0)/gd(0) + 0.5d0*delzz 
     +                             + delzbb(0) - 0.5d0*qd/gd(0)*delaz
              delbb_ct_z(0) = delbb_ct_z_alpha0(0) - delta
          endif


      endif
*
      if (qq.gt.0d0) then
          delqq_ct_g = deluu_ct_g
          delqq_ct_z = deluu_ct_z
      else
          if (mqbig2.gt.0d0) then
              delqq_ct_g = delbb_ct_g
              delqq_ct_z = delbb_ct_z
          else
              delqq_ct_g = deldd_ct_g
              delqq_ct_z = deldd_ct_z
          endif
      endif
*
      call zffg(s,qq,0d0,fzqqg)
      call zffg(s,ql,mlep2,fzllg)

      call ffwp(s,qq,fzqqw(1))
      call ffwp(s,ql,fzllw(1))
                                                                    
      fgqqw(1) = fzqqw(1)
      fgllw(1) = fzllw(1)

      calculatec0=.true.
      call zffwm(s,qq,i3q,fzqqw(0))
      calculatec0=.false.
      call zffwm(s,ql,i3l,fzllw(0))

      call gffwm(s,qq,i3q,fgqqw(0))
      call gffwm(s,ql,i3l,fgllw(0))

      if (mqbig2.gt.0d0) then

          call dgbbwm(s,qq,deltatop)
          fgqqw (0) = fgqqw(0) + deltatop

          call dzbbwm(s,qq,deltatop)
          fzqqw (0) = fzqqw(0) + deltatop

      endif

*
** Eq. (3.30) of Dittmaier-Huber - ArXiv:0911.2329
*
      fzqqw =  fzqqw + delqq_ct_z

      fgqqw =  fgqqw + delqq_ct_g

      fzllw =  fzllw + delll_ct_z

      fgllw =  fgllw + delll_ct_g

*
** Eq. (3.29) of Dittmaier-Huber - ArXiv:0911.2329
*
      vert = zero

      do sig=0,1
          do tau=0,1
          if (complexmasses) then
              vert =   vert 
     +             - el2_scheme * ( qq*ql/s*( fgqqw(sig) + fgllw(tau) )
     +                        + gq(sig)*gl(tau)/(s*cone-mz2) * 
     +                                  ( fzqqw(sig) + fzllw(tau) )  )
     +           *a(sig,tau)*conjg(a(sig,tau))*conjg(bornamp(sig,tau))
          else
              vert =   vert
     +             - el2_scheme * ( qq*ql/s*( fgqqw(sig) + fgllw(tau) )
     +                      + gq(sig)*gl(tau)/(s*cone-mz2+ii*ph_ZmZw) * 
     +                                  ( fzqqw(sig) + fzllw(tau) )  ) 
     +           *a(sig,tau)*conjg(a(sig,tau))*conjg(bornamp(sig,tau)) 
          endif

          enddo
      enddo
*
** photonic vertices
*
      bornm2 = 0d0
      do sig=0,1
          do tau=0,1
              bornm2 = bornm2 + 
     +                  a(sig,tau)*conjg(a(sig,tau))*
     +                    bornamp(sig,tau)*conjg(bornamp(sig,tau))
          enddo
      enddo
*
      vert = vert + ( qq**2*fzqqg + ql**2*fzllg ) * bornm2

      end subroutine deltavert
*
**
** Vertex photonic correction form factor
** Eq. (3.6) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine zffg(s,qf,mf2,fzffg)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
      include'pwhg_math.h'
*
      real*8 s,qf,mf2
      complex*16 fzffg
*
      complex*16 b0dmf20mf,b0dsmfmf,c0dmf2smf2mgmfmf
*
      call b0m0(mf2*cone,mf2*cone,b0dmf20mf)

      call b0reg(s*cone,mf2*cone,mf2*cone,b0dsmfmf)

      call c0ir_ffp(s*cone,mf2*cone,mf2*cone,c0dmf2smf2mgmfmf)

      fzffg =  alsu4pi*( 
     -                   - 2d0*cone
     -                   + 4d0*b0dmf20mf 
     -                   - 3d0*b0dsmfmf
     -                   - 2*s*c0dmf2smf2mgmfmf )

      end subroutine zffg
*
**
** Vertex weak correction form factor for Z/gamma and spin +
** Eq. (B.4) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine ffwp(s,qf,fffwp)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
*
      real*8 s,qf
      complex*16 fffwp
*
      complex*16 b0ds00,b0d00mz,c0d00s0mz0
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d00mz
*
      if (ifirst.eq.0) then

          ifirst = 1

          call b0p0m0(mz2,b0d00mz)

      endif

      call b0reg(s*cone,zero,zero,b0ds00)

      call c0reg2(s*cone,mz2,c0d00s0mz0)
      
      fffwp   = - alsu4pi * qf**2*sw2/cw2 * ( 
     +                      2d0*cone - 2d0*(mz2+2d0*s*cone)/s*b0d00mz
     +                    + ( 3d0*s*cone + 2d0*mz2 )/s*b0ds00
     +                    + 2d0*( mz2 + s*cone )**2/s*c0d00s0mz0 )
*
      end subroutine ffwp
*
**
** Vertex weak correction form factor for gamma and spin -
** Eq. (B.5) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine gffwm(s,qf,i3,fzffgwm)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
*
      real*8 s,qf
      complex*16 i3
      complex*16 fzffgwm
*
      complex*16 b0ds00,b0dsmwmw,b0d00mw,b0d00mz
      complex*16 c0d00s0mz0,c0d00s0mw0,c0d00smw0mw
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d00mw,b0d00mz

      common/forc0/b0ds00,b0dsmwmw,c0d00s0mz0,
     +             c0d00s0mw0,c0d00smw0mw

      logical calculatec0
      common/ifc0/calculatec0
*
      if (ifirst.eq.0) then

          ifirst = 1

          call b0reg(zero,zero,mw2,b0d00mw)
          call b0reg(zero,zero,mz2,b0d00mz)

      endif

      if (calculatec0) then
          call b0reg(s*cone,zero,zero,b0ds00)
          call b0reg(s*cone,mw2,mw2,b0dsmwmw)
    
          call c0reg2(s*cone,mz2,c0d00s0mz0)
          call c0reg2(s*cone,mw2,c0d00s0mw0)
          call c0reg3(s*cone,mw2,c0d00smw0mw)
      endif

      fzffgwm = alsu4pi/2d0*(
     +  1d0/(s*sw2*qf)*(  2d0*qf*(2d0*s*cone+mw2)*b0d00mw
     +       + (2.d0*i3 - qf*cone)*(   2d0*s*cone+(3d0*s*cone+2d0*mw2)
     +                                                           *b0ds00
     +                               + 2d0*(s*cone+mw2)**2*c0d00s0mw0 )
     +       - 2d0*i3*( (s*cone+2d0*mw2)*b0dsmwmw  
     +                 - 2d0*mw2*(2d0*s*cone+mw2)*c0d00smw0mw )
     +          )
     + + (i3-qf*sw2)**2/(s*cw2*sw2)*(
     +                    -4d0*s*cone + 4d0*( 2d0*s*cone + mz2 )*b0d00mz
     +                    -2d0*( 3d0*s*cone + 2d0*mz2 )*b0ds00
     +                    -4d0*( s*cone + mz2 )**2*c0d00s0mz0
     +                 )
     + )

      end subroutine gffwm
*
**
** Vertex weak correction form factor for Z and spin -
** Eq. (B.6) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine zffwm(s,qf,i3,fzffzwm)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
*
      real*8 s,qf
      complex*16 i3
      complex*16 fzffzwm
*
      complex*16 b0ds00,b0dsmwmw,b0d00mw,b0d00mz
      complex*16 c0d00s0mz0,c0d00s0mw0,c0d00smw0mw
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d00mw,b0d00mz

      common/forc0/b0ds00,b0dsmwmw,c0d00s0mz0,
     +             c0d00s0mw0,c0d00smw0mw

      logical calculatec0
      common/ifc0/calculatec0
*
      if (ifirst.eq.0) then

          ifirst = 1

          call b0p0m0(mw2,b0d00mw)
          call b0p0m0(mz2,b0d00mz)

      endif
*
      if (calculatec0) then
          call b0m0m0(s*cone,b0ds00)
          call b0reg(s*cone,mw2,mw2,b0dsmwmw)

          call c0reg2(s*cone,mz2,c0d00s0mz0)
          call c0reg2(s*cone,mw2,c0d00s0mw0)
          call c0reg3(s*cone,mw2,c0d00smw0mw)
      endif
*
      fzffzwm   = alsu4pi/2d0*(
     +     1d0/(s*sw2*(i3-qf*sw2))*(
     +           2d0*( i3*cone - qf*sw2 )*( 2d0*s*cone + mw2 )*b0d00mw
     +         +( i3*cw2 - i3*sw2 + qf*sw2 )*( 2d0*s*cone
     +                                + (3d0*s*cone+2d0*mw2)*b0ds00
     +                                + 2d0*(s*cone+mw2)**2*c0d00s0mw0 )
     +         - 2d0*cw2*i3*( ( s*cone + 2d0*mw2 )*b0dsmwmw
     +                      - 2d0*mw2*( 2d0*s*cone + mw2 )*c0d00smw0mw )
     +         )
     +  + (i3-qf*sw2)**2/(s*cw2*sw2)*(
     +        - 4d0*s*cone + 4d0*( 2d0*s*cone + mz2 )*b0d00mz
     +        - 2d0*( 3d0*s*cone + 2d0*mz2 )*b0ds00
     +        - 4d0*( s*cone + mz2 )**2*c0d00s0mz0
     +    )
     + )

      end subroutine zffwm
*
**
** Delta Vertex weak correction form factor for gamma and spin -
** for B quark
** Eq. (B.8) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine dgbbwm(s,qf,dfgbbwm)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
*
      real*8 s,qf
      complex*16 dfgbbwm
*
      complex*16 b0d00mw,b0d0mtmw,b0d0mwmw,b0dsmwmw
      complex*16 b0ds00,b0dsmtmt
      complex*16 c0d00s0mw0,c0d00smtmwmt
      complex*16 c0d00smw0mw,c0d00smwmtmw
*
      integer ifirst
      data ifirst/0/
      save ifirst
      save b0d00mw,b0d0mtmw,b0d0mwmw

      if (ifirst.eq.0) then

          ifirst = 1

          call b0p0m0(mw2,b0d00mw)
          call b0p0(mt2*cone,mw2,b0d0mtmw)
          call b0p0mm(mw2,b0d0mwmw)

      endif
*
      call b0m0m0(s*cone,b0ds00)
      call b0reg(s*cone,mt2*cone,mt2*cone,b0dsmtmt)
      call c0reg2(s*cone,mw2,c0d00s0mw0)
      call c0reg5(s*cone,mt2*cone,mw2,mt2*cone,c0d00smtmwmt)
      call c0reg3(s*cone,mw2,c0d00smw0mw)
      call c0reg5(s*cone,mw2,mt2*cone,mw2,c0d00smwmtmw)
      call b0reg(s*cone,mw2,mw2,b0dsmwmw)
*
      dfgbbwm   = alsu4pi/2d0*( 
     +            -2d0/(s*sw2)*(
     +                     (2d0*s*cone+mw2)*(b0d00mw-b0d0mtmw) 
     +                    +(3d0*s*cone+2d0*mw2)*(b0ds00-b0dsmtmt) 
     +                    +2d0*(s*cone+mw2)**2*(c0d00s0mw0-c0d00smtmwmt)
     +                    +3d0*mw2*(2d0*s*cone+mw2)*
     +                              (c0d00smw0mw-c0d00smwmtmw)
     +                   )
     +            +mt2/(s*sw2*mw2)*(
     +                     (mt2*cone+mw2)*
     +                              (b0d0mtmw+2d0*b0dsmtmt-3d0*b0dsmwmw)
     +                    +s*(b0dsmtmt-1.5d0*b0dsmwmw-2.5d0*cone)
     +                    +(mw2*(2d0*s*cone+3d0*mw2)-mt2*(mt2+s)*cone)*
     +                        (2d0*c0d00smtmwmt+3d0*c0d00smwmtmw)
     +                   )
     +               ) 

      end subroutine dgbbwm
*
**
** Delta Vertex weak correction form factor for Z and spin -
** for B quark
** Eq. (B.9) of Dittmaier-Huber - ArXiv:0911.2329
**
*
      subroutine dzbbwm(s,qf,dfzbbwm)
      implicit none
      include'pwhg_physpar.h'
      include'mathx.h'
*
      real*8 s,qf
      complex*16 dfzbbwm
*
      complex*16 b0d00mw,b0d0mtmw,b0ds00,b0dsmtmt
      complex*16 b0d0mwmw,b0dsmwmw
      complex*16 c0d00s0mw0,c0d00smtmwmt,c0d00smw0mw,c0d00smwmtmw
*
      integer ifirst
      data ifirst/0/
      save ifirst
      save b0d00mw,b0d0mtmw,b0d0mwmw

      if (ifirst.eq.0) then

          ifirst = 1

          call b0p0m0(mw2,b0d00mw)
          call b0p0(mt2*cone,mw2,b0d0mtmw)
          call b0p0mm(mw2,b0d0mwmw)

      endif

      call b0m0m0(s*cone,b0ds00)
      call b0reg(s*cone,mt2*cone,mt2*cone,b0dsmtmt)
      call c0reg2(s*cone,mw2,c0d00s0mw0)
      call c0reg5(s*cone,mt2*cone,mw2,mt2*cone,c0d00smtmwmt)
      call c0reg3(s*cone,mw2,c0d00smw0mw)
      call c0reg5(s*cone,mw2,mt2*cone,mw2,c0d00smwmtmw)
      call b0reg(s*cone,mw2,mw2,b0dsmwmw)


      dfzbbwm   = alsu4pi/2d0*(
     +          -1d0/(2d0*sw2-3d0*cone)/s/sw2*(
     +             2d0*(2d0*s+mw2)*(2d0*sw2-3d0*cone)*(b0d00mw-b0d0mtmw)
     +            +(3d0*s+2d0*mw2)*(4d0*sw2-3d0*cone)*(b0ds00-b0dsmtmt)
     +            -12d0*cw2*mw2*(2d0*s*cone+mw2)*
     +                          (c0d00smw0mw-c0d00smwmtmw)
     +            +2d0*(4d0*sw2-3d0*cone)*(s*cone+mw2)**2*
     +                          (c0d00s0mw0-c0d00smtmwmt)
     +           )
     +          -mt2/(2d0*sw2-3d0*cone)/s/sw2/mw2*(
     +            (mt2*cone-mw2)*(4d0*b0dsmtmt*sw2
     +              -3d0*b0dsmwmw*(sw2-cw2)+b0d0mtmw*(2d0*sw2-3d0*cone))
     +            +(2d0*s*sw2-6d0*mw2)*b0dsmtmt
     +            +s*(1.5d0-5d0*sw2)
     +            +(6d0*mw2-1.5d0*s*(sw2-cw2))*b0dsmwmw
     +            -12d0*s*mw2*c0d00smtmwmt
     +            -3d0*(4d0*mw2**2-2d0*mt2*mw2-mt2*s*cone)*
     +                      (c0d00smtmwmt+c0d00smwmtmw)
     +            -3d0*(2d0*mw2**2-mt2**2*cone)*c0d00smwmtmw
     +            +2d0*sw2*(mw2*(2d0*s+3d0*mw2)-mt2*(mt2+s)*cone)*
     +                        (2d0*c0d00smtmwmt+3d0*c0d00smwmtmw)
     +           )
     +         )

      end subroutine dzbbwm
*
*
*
**
** Subroutine for the calculation of the 
** contribution to deltavirt due to boxes
**
** (pu)  u  \          / l+  (pl)
**           \ ______ /
**            |      |   (etc)
**            |______| 
**           /        \
** (pud) d~ /          \ l- (pld)
**
**        s = (pu+pud)**2
**        t = (pu-pl)**2 = (pud-pld)**2
**        u = (pu-pld)**2 = (pud-pl)**2
*
* eq. (A.6) of Dittmaier-Kraemer PRD65 073007
*
      subroutine deltabox(s,t,u,mqbig2,box)
      implicit none
      include'pwhg_physpar.h'
      include'pwhg_math.h'
      include'mathx.h'
      include'PhysPars.h'
*
      real*8 s,t,u,mqbig2
      complex*16 box
*
      complex*16 bornamp(0:1,0:1)
      common/bornamps/bornamp

      real*8 mlep2
      common/leptmass/mlep2
*
      complex*16 c0d00smv0mvp,c0d00smvmqmvp,d0d0000stmvmqmvp0
      logical calculated0

      common/ford01/c0d00smv0mvp,c0d00smvmqmvp,d0d0000stmvmqmvp0
      common/ifd0/calculated0
*
      complex*16 tzz(0:1,0:1),tww(0:1,0:1),tgg(0:1,0:1),tzg(0:1,0:1)
      complex*16 uzz(0:1,0:1),uww(0:1,0:1),ugg(0:1,0:1),uzg(0:1,0:1)
      complex*16 fzz(0:1,0:1),fww(0:1,0:1),fgg(0:1,0:1),fzg(0:1,0:1)

      integer sig,tau
*
*
** ZZ
*
      calculated0=.true.
      call btpp(s*cone,t*cone,u*cone,mz2,mz2,zero,tzz(1,1))
      tzz(0,0) = tzz(1,1) 

      calculated0=.false.
      call btpm(s*cone,t*cone,u*cone,mz2,mz2,zero,tzz(1,0))
      tzz(0,1) = tzz(1,0) 

      calculated0=.true.
      call bupm(s*cone,t*cone,u*cone,mz2,mz2,zero,uzz(1,0))
      uzz(0,1) = uzz(1,0) 

      calculated0=.false.
      call bupp(s*cone,t*cone,u*cone,mz2,mz2,zero,uzz(1,1))
      uzz(0,0) = uzz(1,1) 
*
** WW
*

      tww(1,1) = zero
      tww(0,1) = zero
      tww(1,0) = zero

      calculated0=.true.
      if (qq.lt.0d0) then
          call btpp(s*cone,t*cone,u*cone,mw2,mw2,mqbig2*cone,tww(0,0))
      else
          call bupp(s*cone,t*cone,u*cone,mw2,mw2,mqbig2*cone,uww(0,0))
      endif

*
** gg
*

      calculated0=.true.
      call btpmg(s*cone,t*cone,u*cone,mlep2*cone,zero,tgg(1,0))
      calculated0=.false.
      call btppg(s*cone,t*cone,u*cone,mlep2*cone,zero,tgg(1,1))

      tgg(0,1) = tgg(1,0) 
      tgg(0,0) = tgg(1,1) 

      calculated0=.true.
      call buppg(s*cone,t*cone,u*cone,mlep2*cone,zero,ugg(1,1))
      calculated0=.false.
      call bupmg(s*cone,t*cone,u*cone,mlep2*cone,zero,ugg(1,0))

      ugg(0,1) = ugg(1,0) 
      ugg(0,0) = ugg(1,1) 

*
** zg
*

      calculated0=.true.
      call btpmgz(s*cone,t*cone,u*cone,mz2,mlep2*cone,zero,tzg(1,0))
      calculated0=.false.
      call btppgz(s*cone,t*cone,u*cone,mz2,mlep2*cone,zero,tzg(1,1))

      tzg(0,1) = tzg(1,0) 
      tzg(0,0) = tzg(1,1) 

      calculated0=.true.
      call buppgz(s*cone,t*cone,u*cone,mz2,mlep2*cone,zero,uzg(1,1))
      calculated0=.false.
      call bupmgz(s*cone,t*cone,u*cone,mz2,mlep2*cone,zero,uzg(1,0))

      uzg(0,1) = uzg(1,0) 
      uzg(0,0) = uzg(1,1) 
*
**
*
      fww(1,1) = zero
      fww(0,1) = zero
      fww(1,0) = zero

      if (qq.lt.0d0) then
          fww(0,0) = tww(0,0)/4d0/sw4
      else                    
          fww(0,0) = uww(0,0)/4d0/sw4
      endif
*
**
*
      do sig=0,1
          do tau=0,1
              fzz(sig,tau) = (gq(sig)*gl(tau))**2
     +                                 *( tzz(sig,tau) + uzz(sig,tau) )
              fgg(sig,tau) = (qq*ql)**2*( tgg(sig,tau) + ugg(sig,tau) )

              fzg(sig,tau) = 2d0*qq*ql*gq(sig)*gl(tau)
     +                                 *( tzg(sig,tau) + uzg(sig,tau) )
          enddo
      enddo

      box = zero

      do sig=0,1
          do tau=0,1
          if (complexmasses) then
              box = box + 
     +                    alpha*ph_alphaem*( 
     +                      + fww(sig,tau) + fzz(sig,tau) 
     +                      + fgg(sig,tau) + fzg(sig,tau) ) 
     +                    *a(sig,tau)
     +                    *conjg(a(sig,tau))*conjg(bornamp(sig,tau))
          else
              box = box + 
     +                    alpha*ph_alphaem*( 
     +                      + fww(sig,tau) + fzz(sig,tau) 
     +                      + fgg(sig,tau) + fzg(sig,tau) ) 
     +                    *a(sig,tau)
     +                    *conjg(a(sig,tau))*conjg(bornamp(sig,tau))
          endif
          enddo
      enddo

      end subroutine deltabox
*
** b functions for deltabox eq. (B.12-15) Dittmaier-Huber
** v, vp boson != gamma
*
      subroutine btpm(s,t,u,mv2,mvp2,mqbig2,bout)
* mqbig2 is only important for top quark and WW
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,mv2,mvp2,mqbig2
      complex*16 bout
*
      complex*16 c0d00smv0mvp,c0d00smvmqmvp,d0d0000stmvmqmvp0
      logical calculated0
*
      common/ford01/c0d00smv0mvp,c0d00smvmqmvp,d0d0000stmvmqmvp0
      common/ifd0/calculated0

      if (calculated0) then
          call c0reg1(s,mv2,mvp2,c0d00smv0mvp)

          call c0reg5(s,mv2,mqbig2,mvp2,c0d00smvmqmvp)

          call myd0reg(zero,zero,zero,zero,s,t,
     +                 mv2,mqbig2,mvp2,zero,
     +                 d0d0000stmvmqmvp0)
      endif
*
      bout = -2d0*(  c0d00smv0mvp + c0d00smvmqmvp 
     +             - (t-mqbig2)*d0d0000stmvmqmvp0 )
*
      end subroutine btpm

      subroutine btpp(s,t,u,mv2,mvp2,mqbig2,bout)
* mqbig2 is only important for top quark and WW
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,mv2,mvp2,mqbig2
      complex*16 bout
*
      complex*16 b0dsmvmvp,b0dtmq0
      complex*16 c0d00t0mvmq,c0d00t0mvpmq,c0d00smv0mvp,c0d00smvmqmvp
      complex*16 d0d0000stmvmqmvp0
      logical calculated0
* 
      common/ford01/c0d00smv0mvp,c0d00smvmqmvp,d0d0000stmvmqmvp0
      common/ifd0/calculated0

      call b0reg(s,mv2,mvp2,b0dsmvmvp)
      call b0m0(t,mqbig2,b0dtmq0)

      call c0reg4(t,mv2,mqbig2,c0d00t0mvmq)
      call c0reg4(t,mvp2,mqbig2,c0d00t0mvpmq)

      if (calculated0) then
          call c0reg1(s,mv2,mvp2,c0d00smv0mvp)
          call c0reg5(s,mv2,mqbig2,mvp2,c0d00smvmqmvp)

          call myd0reg(zero,zero,zero,zero,s,t,
     +                 mv2,mqbig2,mvp2,zero,
     +                 d0d0000stmvmqmvp0)
      endif

      bout = ( 2d0*u*(b0dsmvmvp-b0dtmq0) 
     +        - t*(mqbig2-mv2-mvp2-t+u)*(c0d00t0mvmq + c0d00t0mvpmq)
     +        - ( t**2 + u**2 + s*(mqbig2-mv2-mvp2) )*
     +                        (c0d00smv0mvp + c0d00smvmqmvp)
     +        + (   t*(mv2+mvp2-mqbig2-2d0*s)**2 
     +           + (u*(2d0*u-mqbig2)-2d0*s**2)*(mv2+mvp2-mqbig2-2d0*s)
     +           - 2d0*u*(u**2-mv2*mvp2) + u*s*(s-mqbig2) - s**3 )
     +          *d0d0000stmvmqmvp0
     +       ) / (u**2)

      end subroutine btpp

      subroutine bupm(s,t,u,mv2,mvp2,mqbig2,bout)
* mqbig2 is only important for top quark and WW
      implicit none
      complex*16 s,t,u,mv2,mvp2,mqbig2
      complex*16 bout
*
      call btpp(s,u,t,mv2,mvp2,mqbig2,bout)

      bout = - bout

      end subroutine bupm

      subroutine bupp(s,t,u,mv2,mvp2,mqbig2,bout)
* mqbig2 is only important for top quark and WW
      implicit none
      complex*16 s,t,u,mv2,mvp2,mqbig2
      complex*16 bout
*
      call btpm(s,u,t,mv2,mvp2,mqbig2,bout)

      bout = - bout

      end subroutine bupp
*
** b functions for deltabox eq. (B.12-15) Dittmaier-Huber
** v, vp boson = gamma
*
      subroutine btpmg(s,t,u,ml2,mq2,bout)
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,mq2,ml2
      complex*16 bout
*
      complex*16 c0dml2ml2s0ml0
      complex*16 c0dmq2mq2s0mq0
      complex*16 d0mq2mq2ml2ml2st0mq0ml
      logical calculated0

      common/ford02/c0dml2ml2s0ml0,c0dmq2mq2s0mq0,d0mq2mq2ml2ml2st0mq0ml
      common/ifd0/calculated0
*
      if (calculated0) then
          call c0reg6(s,ml2,c0dml2ml2s0ml0)
          call c0reg6(s,mq2,c0dmq2mq2s0mq0)
          call d02g(mq2,mq2,ml2,ml2,s,t,
     +              zero,mq2,zero,ml2,
     +              d0mq2mq2ml2ml2st0mq0ml)
      endif

      bout = -2d0*(  c0dml2ml2s0ml0 + c0dmq2mq2s0mq0
     +             - t*d0mq2mq2ml2ml2st0mq0ml ) 

      end subroutine btpmg
*
      subroutine btppg(s,t,u,ml2,mq2,bout)
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,ml2,mq2
      complex*16 bout
*
      complex*16 b0ds00,b0dt00
      complex*16 c0dml2mq2tmlmgmq,c0dml2ml2s0ml0,c0dmq2mq2s0mq0
      complex*16 d0mq2mq2ml2ml2st0mq0ml
      logical calculated0

      common/ford02/c0dml2ml2s0ml0,c0dmq2mq2s0mq0,d0mq2mq2ml2ml2st0mq0ml
      common/ifd0/calculated0
*      
      call b0m0m0(s,b0ds00)
      call b0m0m0(t,b0dt00)

      call c0lq(ml2,mq2,t,zero,c0dml2mq2tmlmgmq)

      if (calculated0) then
          call c0reg6(s,ml2,c0dml2ml2s0ml0)
          call c0reg6(s,mq2,c0dmq2mq2s0mq0)
          call d02g(mq2,mq2,ml2,ml2,s,t,
     +              zero,mq2,zero,ml2,
     +              d0mq2mq2ml2ml2st0mq0ml)
      endif
*
      bout = (   2d0*u*( b0ds00 - b0dt00 ) 
     +         - t*(-t+u)*( 2d0*c0dml2mq2tmlmgmq )
     +         - (t**2+u**2)*( c0dml2ml2s0ml0 + c0dmq2mq2s0mq0 )
     +         + (  t*(-2d0*s)**2 + (u*2d0*u-2d0*s**2)*(-2d0*s)
     +            - 2d0*u*u**2 + u*s*s - s**3 )
     +            *d0mq2mq2ml2ml2st0mq0ml
     +       ) / (u**2)

      end subroutine btppg
*
      subroutine bupmg(s,t,u,ml2,mq2,bout)
      implicit none
      complex*16 s,t,u,ml2,mq2
      complex*16 bout

      call btppg(s,u,t,ml2,mq2,bout)

      bout = - bout

      end subroutine bupmg
*
      subroutine buppg(s,t,u,ml2,mq2,bout)
      implicit none
      complex*16 s,t,u,ml2,mq2
      complex*16 bout

      call btpmg(s,u,t,ml2,mq2,bout)

      bout = - bout

      end subroutine buppg
*
** b functions for deltabox eq. (B.12-15) Dittmaier-Huber
** v = z, vp = gamma
*
      subroutine btpmgz(s,t,u,mv2,ml2,mq2,bout)
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,mv2,ml2,mq2
      complex*16 bout
*
      complex*16 c0dml2ml2smvml0,c0dmq2mq2smvmq0,d0d0000stmvmqmvp0
      logical calculated0

      common/ford03/c0dml2ml2smvml0,c0dmq2mq2smvmq0,d0d0000stmvmqmvp0
      common/ifd0/calculated0
*
      if (calculated0) then
          call c0fz(ml2,s,ml2,mv2,c0dml2ml2smvml0)
          call c0fz(mq2,s,mq2,mv2,c0dmq2mq2smvmq0)
          call d0sing(mq2,mq2,ml2,ml2,s,t,
     +                zero,mq2,mv2,ml2,
     +                d0d0000stmvmqmvp0)
      endif

      bout   = -2d0*(  c0dml2ml2smvml0 + c0dmq2mq2smvmq0
     +               - (t)*d0d0000stmvmqmvp0 )

      end subroutine btpmgz
*
      subroutine btppgz(s,t,u,mv2,ml2,mq2,bout)
      implicit none
      include'mathx.h'
*
      complex*16 s,t,u,mv2,ml2,mq2
      complex*16 bout
*
      complex*16 b0dsmv0,b0dt00
      complex*16 c0d00t0mv0,c0dml2mq2tmlmgmq,c0dml2ml2smvml0,
     +           c0dmq2mq2smvmq0
      complex*16 d0d0000stmvmqmvp0
      logical calculated0

      common/ford03/c0dml2ml2smvml0,c0dmq2mq2smvmq0,d0d0000stmvmqmvp0
      common/ifd0/calculated0
*      
      call b0m0(s,mv2,b0dsmv0)
      call b0m0m0(t,b0dt00)

      call c0reg2(t,mv2,c0d00t0mv0)
      call c0lq(ml2,mq2,t,zero,c0dml2mq2tmlmgmq)

      if (calculated0) then
          call c0fz(ml2,s,ml2,mv2,c0dml2ml2smvml0)
          call c0fz(mq2,s,mq2,mv2,c0dmq2mq2smvmq0)

          call d0sing(mq2,mq2,ml2,ml2,s,t,
     +                zero,mq2,mv2,ml2,
     +                d0d0000stmvmqmvp0)
      endif
*
      bout = (   2d0*u*(b0dsmv0-b0dt00) 
     +         - t*(-mv2-t+u)*(c0d00t0mv0 + c0dml2mq2tmlmgmq)
     +         - ( t**2 + u**2 - s*mv2 )*
     +                       (c0dml2ml2smvml0 + c0dmq2mq2smvmq0)
     +         + (   t*(mv2-2d0*s)**2 + 2d0*(u**2-s**2)*(mv2-2d0*s)
     +             - 2d0*u**3 + u*s**2 - s**3 )
     +           *d0d0000stmvmqmvp0
     +      ) / (u**2)
*
      end subroutine btppgz
*
      subroutine bupmgz(s,t,u,mv2,ml2,mq2,bout)
      implicit none
      complex*16 s,u,t,mv2,ml2,mq2
      complex*16 bout

      call btppgz(s,u,t,mv2,ml2,mq2,bout)

      bout = - bout

      end subroutine bupmgz
*
      subroutine buppgz(s,t,u,mv2,ml2,mq2,bout)
      implicit none
      complex*16 s,t,u,mv2,ml2,mq2
      complex*16 bout

      call btpmgz(s,u,t,mv2,ml2,mq2,bout)

      bout = - bout

      end subroutine buppgz

*
*************************************************************
*
* EXTERNAL WAVE FUNCTION COUNTERTERMS
*
* eq. (3.20) of ArXiv:0709.1075 (Denner Fortschritte) and 4.17 of
* ArXiv:0505042
*
      subroutine deltazfl(qf,gfm,gfp,mf2,mfp2,dzflout)
* only flavour diagonal terms are considered since ckm=1 always
*
      implicit none
      include'PhysPars.h'
*
      complex*16 qf,mf2,mfp2,gfm,gfp
      complex*16 dzflout
*
      complex*16 sfl,sflp,sfrp,sfsp
*
      call sigmafl (qf,gfm,gfp,mf2,mfp2,mf2,sfl )
      call sigmaflp(qf,gfm,gfp,mf2,mfp2,mf2,sflp)
      call sigmafrp(qf,gfm,gfp,mf2,mfp2,mf2,sfrp)
      call sigmafsp(qf,gfm,gfp,mf2,mfp2,mf2,sfsp)

      if (complexmasses) then
          dzflout = - sfl - mf2*( sflp + sfrp + 2.d0*sfsp )
      else
          dzflout = dble(- sfl - mf2*( sflp + sfrp + 2.d0*sfsp ))
      endif
 
      end subroutine deltazfl
*
* eq. (3.20) of ArXiv:0709.1075 (Denner Fortschritte) and 4.17 of
* ArXiv:0505042
*
      subroutine deltazfr(qf,gfm,gfp,mf2,mfp2,dzfrout)
* only flavour diagonal terms are considered since ckm=1 always
*
      implicit none
      include'PhysPars.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2
      complex*16 dzfrout
*
      complex*16 sfr,sflp,sfrp,sfsp
*
      call sigmafr (qf,gfm,gfp,mf2,mfp2,mf2,sfr )
      call sigmaflp(qf,gfm,gfp,mf2,mfp2,mf2,sflp)
      call sigmafrp(qf,gfm,gfp,mf2,mfp2,mf2,sfrp)
      call sigmafsp(qf,gfm,gfp,mf2,mfp2,mf2,sfsp)

      if (complexmasses) then
          dzfrout = - sfr - mf2*( sflp + sfrp + 2.d0*sfsp )
      else
          dzfrout = dble(- sfr - mf2*( sflp + sfrp + 2.d0*sfsp ))
      endif

      end subroutine deltazfr
*
* generic fermion self-energies
*
* eq. (B.6) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmafr(qf,gfm,gfp,mf2,mfp2,s,sfrout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2,s
      complex*16 sfrout
*
      complex*16 b1dsmfmg,b1dsmfmz,b1dsmfmh,b1dsmfpmw
*
      if(abs(qf).lt.epsilon) then

          b1dsmfmg = zero

      else

          call b1reg(s*cone,mf2*cone,zero,b1dsmfmg)

      endif

      call b1reg(s*cone,mf2*cone,mz2,b1dsmfmz)
      call b1reg(s*cone,mf2*cone,mh2,b1dsmfmh)
      call b1reg(s*cone,mfp2*cone,mw2,b1dsmfpmw)

      sfrout  = - alsu4pi * (  
     +               qf**2 * (2.d0*b1dsmfmg + cone) 
     +             + gfp**2 * (2.d0*b1dsmfmz + cone)
     +             + 0.5d0/sw2*mf2/2.d0/mw2*(b1dsmfmz + b1dsmfmh)
     +             + 0.5d0/sw2*sqrt(mf2)*sqrt(mfp2)/mw2*b1dsmfpmw  
     +          )

      end subroutine sigmafr
*
* eq. (B.6) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmafl(qf,gfm,gfp,mf2,mfp2,s,sflout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2,s
      complex*16 sflout
*
      complex*16 b1dsmfmg,b1dsmfmz,b1dsmfmh,b1dsmfpmw
*
      if(abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b1dsmfmg = zero

      else

          call b1reg(s*cone,mf2*cone,zero,b1dsmfmg)

      endif

      call b1reg(s*cone,mf2*cone,mz2,b1dsmfmz)
      call b1reg(s*cone,mf2*cone,mh2,b1dsmfmh)
      call b1reg(s*cone,mfp2*cone,mw2,b1dsmfpmw)

      sflout  = - alsu4pi * (   
     +               qf**2 * (2.d0*b1dsmfmg + cone) 
     +             + gfm**2 * (2.d0*b1dsmfmz + cone)
     +             + 0.5d0/sw2*mf2/2.d0/mw2*(b1dsmfmz + b1dsmfmh)
     +             + 0.5d0/sw2*((2.d0*cone+mfp2/mw2)*b1dsmfpmw+cone)  
     +            )

      end subroutine sigmafl
*
* generic fermion self-energies derivatives w.r.t. p^2
*
* derivative (w.r.t. p^2 (s)) of eq. (B.6) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmaflp(qf,gfm,gfp,mf2,mfp2,s,sflpout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2,s
      complex*16 sflpout
*

      complex*16 b1pdsmfmg,b1pdsmfmz,b1pdsmfmh,b1pdsmfpmw
*
      if (abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b1pdsmfmg = zero
      
      else

          call b1pir(s*cone,mf2*cone,b1pdsmfmg)

      endif

      call b1preg(s*cone,mf2*cone,mz2,b1pdsmfmz)
      call b1preg(s*cone,mf2*cone,mh2,b1pdsmfmh)
      call b1preg(s*cone,mfp2*cone,mw2,b1pdsmfpmw)

      sflpout  = - alsu4pi * ( 
     +                qf**2  * (2.d0*b1pdsmfmg) 
     +              + gfm**2 * (2.d0*b1pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b1pdsmfmz + b1pdsmfmh)
     +              + 0.5d0/sw2*(2.d0*cone+mfp2/mw2)*b1pdsmfpmw 
     +            )

      end subroutine sigmaflp
*
* derivative (w.r.t. p^2 (s)) of eq. (B.7) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmafrp(qf,gfm,gfp,mf2,mfp2,s,sfrpout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2,s
      complex*16 sfrpout
*      

      complex*16 b1pdsmfmg,b1pdsmfmz,b1pdsmfmh,b1pdsmfpmw
*
      if (abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b1pdsmfmg = zero
      
      else

          call b1pir(s*cone,mf2*cone,b1pdsmfmg)

      endif

      call b1preg(s*cone,mf2*cone,mz2,b1pdsmfmz)
      call b1preg(s*cone,mf2*cone,mh2,b1pdsmfmh)
      call b1preg(s*cone,mfp2*cone,mw2,b1pdsmfpmw)

      sfrpout  = - alsu4pi * (
     +                qf**2 * (2.d0*b1pdsmfmg) 
     +              + gfp**2 * (2.d0*b1pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b1pdsmfmz + b1pdsmfmh)
     +              + 0.5d0/sw2*mfp2/mw2*b1pdsmfpmw 
     +         )

      end subroutine sigmafrp
*
* derivative (w.r.t. p^2 (s)) of eq. (B.8) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmafsp(qf,gfm,gfp,mf2,mfp2,s,sfspout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 qf,gfm,gfp,mf2,mfp2,s
      complex*16 sfspout
*
      complex*16 b0pdsmfmg,b0pdsmfmz,b0pdsmfmh,b0pdsmfpmw
*
      if (abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b0pdsmfmg   = zero

      else

          call b0pir(s*cone,mf2*cone,b0pdsmfmg)

      endif

      call b0preg(s*cone,mf2*cone,mz2,b0pdsmfmz)
      call b0preg(s*cone,mf2*cone,mh2,b0pdsmfmh)
      call b0preg(s*cone,mfp2*cone,mw2,b0pdsmfpmw)

      sfspout  = - alsu4pi * ( 
     +                qf**2 * 4.d0*b0pdsmfmg 
     +              + gfp*gfm * (4.d0*b0pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b0pdsmfmz - b0pdsmfmh)
     +              + 0.5d0/sw2*mfp2/mw2*b0pdsmfpmw 
     +      )

      end subroutine sigmafsp
*
** end of external wavefunction conuterterms
*
************************************************************************
************************************************************************
*
** VECTOR BOSON SELF-ENERGIES AND THEIR DERIVATIVES W.R.T. P^2 WHEN USED
*
* eq. (B.2) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmaazt(s,saztout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 saztout
*
      complex*16 b0dsmeme,b0dsmmmm,b0dsmtlmtl,
     +           b0dsmumu,b0dsmdmd, 
     +           b0dsmcmc,b0dsmsms,
     +           b0dsmtmt,b0dsmbmb,
     +           b0dsmwmw,b0d0mwmw,
     +           b0d0meme,b0d0mmmm,b0d0mtlmtl,
     +           b0d0mumu,b0d0mdmd, 
     +           b0d0mcmc,b0d0msms,
     +           b0d0mtmt,b0d0mbmb
*
      integer ifirst
      data ifirst /0/
      save ifirst 
      save b0d0meme,b0d0mmmm,b0d0mtlmtl,
     +     b0d0mumu,b0d0mdmd,b0d0mcmc,
     +     b0d0msms,b0d0mtmt,b0d0mbmb,
     +     b0d0mwmw


      if (ifirst.eq.0) then
          ifirst = 1
* leptons
          call b0p0mm(me2*cone,b0d0meme)
          call b0p0mm(mm2*cone,b0d0mmmm)
          call b0p0mm(mtl2*cone,b0d0mtlmtl)

* quarks
          call b0p0mm(mu2*cone,b0d0mumu)
          call b0p0mm(md2*cone,b0d0mdmd)
          call b0p0mm(mc2*cone,b0d0mcmc)
          call b0p0mm(ms2*cone,b0d0msms)
          call b0p0mm(mt2*cone,b0d0mtmt)
          call b0p0mm(mb2*cone,b0d0mbmb)

* W
          call b0p0mm(mw2,b0d0mwmw)

      endif

* leptons
      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)

* quarks
      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)
      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)

* W
      call b0reg(s,mw2,mw2,b0dsmwmw)

***

      saztout  =  2.d0/3.d0*(
* sum over three charged leptons
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*me2*cone)*b0dsmeme
     +                     +2.d0*me2*b0d0meme + 1.d0/3.d0*s )
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*mm2*cone)*b0dsmmmm
     +                     +2.d0*mm2*b0d0mmmm + 1.d0/3.d0*s )
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*mtl2*cone)*b0dsmtlmtl
     +                     +2.d0*mtl2*b0d0mtlmtl + 1.d0/3.d0*s )
* sum over quarks
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mu2*cone)*b0dsmumu
     +                          +2.d0*mu2*b0d0mumu + 1.d0/3.d0*s )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*md2*cone)*b0dsmdmd
     +                          +2.d0*md2*b0d0mdmd + 1.d0/3.d0*s )
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mc2*cone)*b0dsmcmc
     +                          +2.d0*mc2*b0d0mcmc + 1.d0/3.d0*s )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*ms2*cone)*b0dsmsms
     +                          +2.d0*ms2*b0d0msms + 1.d0/3.d0*s )
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mt2*cone)*b0dsmtmt
     +                          +2.d0*mt2*b0d0mtmt + 1.d0/3.d0*s )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*mb2*cone)*b0dsmbmb
     +                          +2.d0*mb2*b0d0mbmb + 1.d0/3.d0*s )  )
* bosonic part
     +    -1.d0/3.d0/sw/cw*( ((9.d0*cw2+0.5d0*cone)*s
     +                        +(12.d0*cw2+4.d0*cone)*mw2)*b0dsmwmw
     +                      -(12.d0*cw2-2.d0*cone)*mw2*b0d0mwmw
     +                      +1.d0/3.d0*s )
      saztout  = - alsu4pi * saztout

      end subroutine sigmaazt

******************************************************************
*
* derivative w.r.t. p2 eq. (B.2) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmaaztp(s,saztpout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 saztpout
*
      complex*16 b0dsmeme,b0dsmmmm,b0dsmtlmtl,
     +           b0dsmumu,b0dsmdmd, 
     +           b0dsmcmc,b0dsmsms,
     +           b0dsmtmt,b0dsmbmb,
     +           b0dsmwmw,
     +           b0pdsmeme,b0pdsmmmm,b0pdsmtlmtl,
     +           b0pdsmumu,b0pdsmdmd, 
     +           b0pdsmcmc,b0pdsmsms,
     +           b0pdsmtmt,b0pdsmbmb,
     +           b0pdsmwmw

* leptons
      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)

      call b0preg(s,me2*cone,me2*cone,b0pdsmeme)
      call b0preg(s,mm2*cone,mm2*cone,b0pdsmmmm)
      call b0preg(s,mtl2*cone,mtl2*cone,b0pdsmtlmtl)

* quarks
      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)
      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)

      call b0preg(s,mu2*cone,mu2*cone,b0pdsmumu)
      call b0preg(s,md2*cone,md2*cone,b0pdsmdmd)
      call b0preg(s,mc2*cone,mc2*cone,b0pdsmcmc)
      call b0preg(s,ms2*cone,ms2*cone,b0pdsmsms)
      call b0preg(s,mt2*cone,mt2*cone,b0pdsmtmt)
      call b0preg(s,mb2*cone,mb2*cone,b0pdsmbmb)

* W
      call b0reg(s,mw2,mw2,b0dsmwmw)

      call b0preg(s,mw2,mw2,b0pdsmwmw)

***

      saztpout  =  2.d0/3.d0*(
* sum over three charged leptons
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*me2*cone)*b0pdsmeme
     +                     - b0dsmeme + cone/3.d0 )
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*mm2*cone)*b0pdsmmmm
     +                     - b0dsmmmm + cone/3.d0 )
     +     -ql*(gl(1)+gl(0))*( - (s+2.d0*mtl2*cone)*b0pdsmtlmtl
     +                     - b0dsmtlmtl + cone/3.d0 )
* sum over quarks
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mu2*cone)*b0pdsmumu
     +                          - b0dsmumu + cone/3.d0 )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*md2*cone)*b0pdsmdmd
     +                          - b0dsmdmd + cone/3.d0 )
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mc2*cone)*b0pdsmcmc
     +                          - b0dsmcmc + cone/3.d0 )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*ms2*cone)*b0pdsmsms
     +                          - b0dsmsms + cone/3.d0 )
     +     -3.d0*qu*(gu(1)+gu(0))*( - (s+2.d0*mt2*cone)*b0pdsmtmt
     +                          - b0dsmtmt + cone/3.d0 )
     +     -3.d0*qd*(gd(1)+gd(0))*( - (s+2.d0*mb2*cone)*b0pdsmbmb
     +                          - b0dsmbmb + cone/3.d0 )  )
* bosonic part
     +    -1.d0/3.d0/sw/cw*( ((9.d0*cw2+0.5d0*cone)*s
     +                        +(12.d0*cw2+4.d0*cone)*mw2)*b0pdsmwmw
     +                      +(9.d0*cw2+0.5d0*cone)*b0dsmwmw     
     +                      + cone/3.d0 )
      saztpout  = - alsu4pi * saztpout

      end subroutine sigmaaztp
*
*
******************************************************************
*
* eq. (B.3) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmazzt(s,szztout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 szztout
*
      complex*16 b0ds00,
     +           b0dsmeme,b0dsmmmm,b0dsmtlmtl,
     +           b0d0meme,b0d0mmmm,b0d0mtlmtl,
     +           b0dsmumu,b0dsmdmd, 
     +           b0dsmcmc,b0dsmsms,
     +           b0dsmtmt,b0dsmbmb,
     +           b0d0mumu,b0d0mdmd, 
     +           b0d0mcmc,b0d0msms,
     +           b0d0mtmt,b0d0mbmb,
     +           b0dsmwmw,b0d0mwmw,
     +           b0dsmzmh,b0d0mzmh,
     +           b0d0mzmz,b0d0mhmh

      integer ifirst
      data ifirst /0/
      save ifirst 
      save b0d0meme,b0d0mmmm,b0d0mtlmtl,
     +     b0d0mumu,b0d0mdmd,b0d0mcmc,
     +     b0d0msms,b0d0mtmt,b0d0mbmb,
     +     b0d0mzmz,b0d0mhmh,b0d0mwmw,b0d0mzmh

*
      if (ifirst.eq.0) then
          ifirst = 1

* leptons
          call b0p0mm(me2*cone,b0d0meme)
          call b0p0mm(mm2*cone,b0d0mmmm)
          call b0p0mm(mtl2*cone,b0d0mtlmtl)

* quarks
          call b0p0mm(mu2*cone,b0d0mumu)
          call b0p0mm(md2*cone,b0d0mdmd)
          call b0p0mm(mc2*cone,b0d0mcmc)
          call b0p0mm(ms2*cone,b0d0msms)
          call b0p0mm(mt2*cone,b0d0mtmt)
          call b0p0mm(mb2*cone,b0d0mbmb)

* bosons
          call b0p0mm(mz2,b0d0mzmz)
          call b0p0mm(mh2,b0d0mhmh)
          call b0p0mm(mw2,b0d0mwmw)
          call b0p0(mz2,mh2,b0d0mzmh)
      endif

* nu
      call b0m0m0(s,b0ds00)

* leptons
      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)

* quarks
      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)
      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)

* bosons
      call b0reg(s,mw2,mw2,b0dsmwmw)
      call b0reg(s,mz2,mh2,b0dsmzmh)

***

      szztout  = 2.d0/3.d0*(
* sum over three lepton families
     +      3.d0*(gn(1)**2+gn(0)**2)*( - (s)*b0ds00 + 1.d0/3.d0*s)
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*me2*cone)*b0dsmeme
     +                  +2.d0*me2*b0d0meme + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*me2*b0dsmeme
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*mm2*cone)*b0dsmmmm
     +                  +2.d0*mm2*b0d0mmmm + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mm2*b0dsmmmm
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*mtl2*cone)*b0dsmtlmtl
     +                  +2.d0*mtl2*b0d0mtlmtl + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mtl2*b0dsmtlmtl 
* sum over quarks
     +     +3.d0*( (gu(1)**2+gu(0)**2)*(- (s+2.d0*mu2*cone)*b0dsmumu
     +                  +2.d0*mu2*b0d0mumu + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*mu2*b0dsmumu 
     +           + (gd(1)**2+gd(0)**2)*(- (s+2.d0*md2*cone)*b0dsmdmd
     +                  +2.d0*md2*b0d0mdmd + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*md2*b0dsmdmd 
     +           + (gu(1)**2+gu(0)**2)*(- (s+2.d0*mc2*cone)*b0dsmcmc
     +                  +2.d0*mc2*b0d0mcmc + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*mc2*b0dsmcmc 
     +           + (gd(1)**2+gd(0)**2)*(- (s+2.d0*ms2*cone)*b0dsmsms
     +                  +2.d0*ms2*b0d0msms + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*ms2*b0dsmsms 
     +           + (gu(1)**2+gu(0)**2)*(- (s+2.d0*mt2*cone)*b0dsmtmt
     +                  +2.d0*mt2*b0d0mtmt + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mt2*b0dsmtmt 
     +           + (gd(1)**2+gd(0)**2)*(- (s+2.d0*mb2*cone)*b0dsmbmb
     +                  +2.d0*mb2*b0d0mbmb + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mb2*b0dsmbmb    )    )
* bosonic part
     +    + 1.d0/6.d0/sw2/cw2 * (
     +             ( (18.d0*cw4+2.d0*cw2-0.5d0*cone)*s
     +               +(24.d0*cw4+16.d0*cw2-10.d0*cone)*mw2 )*b0dsmwmw
     +             -(24.d0*cw4-8.d0*cw2+2.d0*cone)*mw2*b0d0mwmw
     +             +(4.d0*cw2-cone)/3.d0*s )
     +    + 1.d0/12.d0/sw2/cw2 * (
     +             (2.d0*mh2-10.d0*mz2-s)*b0dsmzmh
     +            -2.d0*mz2*b0d0mzmz - 2.d0*mh2*b0d0mhmh
     +            -(mz2-mh2)**2/s*(b0dsmzmh-b0d0mzmh)
     +            -2.d0/3.d0*s)
      szztout = - alsu4pi * szztout

      end subroutine sigmazzt
*
******************************************************************
*
* derivative of eq. (B.3) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmazztp(s,szztpout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 szztpout
*
      complex*16 b0ds00,
     +           b0dsmeme,b0dsmmmm,b0dsmtlmtl,
     +           b0dsmumu,b0dsmdmd, 
     +           b0dsmcmc,b0dsmsms,
     +           b0dsmtmt,b0dsmbmb,
     +           b0dsmwmw,
     +           b0dsmzmh

      complex*16 b0pds00,
     +           b0pdsmeme,b0pdsmmmm,b0pdsmtlmtl,
     +           b0pdsmumu,b0pdsmdmd,
     +           b0pdsmcmc,b0pdsmsms,
     +           b0pdsmtmt,b0pdsmbmb,
     +           b0pdsmwmw,
     +           b0pdsmzmh,
     +           b0d0mzmh

* nu
      call b0m0m0(s,b0ds00)
      call b0preg(s,zero,zero,b0pds00)

* leptons
      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0preg(s,me2*cone,me2*cone,b0pdsmeme)

      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0preg(s,mm2*cone,mm2*cone,b0pdsmmmm)

      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)
      call b0preg(s,mtl2*cone,mtl2*cone,b0pdsmtlmtl)

* quarks
      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0preg(s,mu2*cone,mu2*cone,b0pdsmumu)

      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0preg(s,md2*cone,md2*cone,b0pdsmdmd)

      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0preg(s,mc2*cone,mc2*cone,b0pdsmcmc)

      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0preg(s,ms2*cone,ms2*cone,b0pdsmsms)

      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)
      call b0preg(s,mt2*cone,mt2*cone,b0pdsmtmt)

      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)
      call b0preg(s,mb2*cone,mb2*cone,b0pdsmbmb)

* bosons
      call b0reg(s,mw2,mw2,b0dsmwmw)
      call b0preg(s,mw2,mw2,b0pdsmwmw)

      call b0reg(s,mz2,mh2,b0dsmzmh)
      call b0p0(mz2,mh2,b0d0mzmh)
      call b0preg(s,mz2,mh2,b0pdsmzmh)

***

      szztpout= 2.d0/3.d0*(
* sum over three neutrino families
     +      3.d0*(gn(1)**2+gn(0)**2)*(-s*b0pds00-b0ds00+cone/3.d0)
* sum over three lepton families
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*me2*cone)*b0pdsmeme
     +                  - (1.d0)*b0dsmeme + cone/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*me2*b0pdsmeme
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*mm2*cone)*b0pdsmmmm
     +                  - (1.d0)*b0dsmmmm + cone/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*mm2*b0pdsmmmm
     +     +(gl(1)**2+gl(0)**2)*(- (s+2.d0*mtl2*cone)*b0pdsmtlmtl
     +                  - (1.d0)*b0dsmtlmtl + cone/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*mtl2*b0pdsmtlmtl
* sum over quarks
     +     +3.d0*( (gu(1)**2+gu(0)**2)*(- (s+2.d0*mu2*cone)*b0pdsmumu
     +                  - b0dsmumu + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mu2*b0pdsmumu )
     +     +3.d0*( (gd(1)**2+gd(0)**2)*(- (s+2.d0*md2*cone)*b0pdsmdmd
     +                  - b0dsmdmd + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*md2*b0pdsmdmd )
     +     +3.d0*( (gu(1)**2+gu(0)**2)*(- (s+2.d0*mc2*cone)*b0pdsmcmc
     +                  - b0dsmcmc + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mc2*b0pdsmcmc )
     +     +3.d0*( (gd(1)**2+gd(0)**2)*(- (s+2.d0*ms2*cone)*b0pdsmsms
     +                  - b0dsmsms + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*ms2*b0pdsmsms )
     +     +3.d0*( (gu(1)**2+gu(0)**2)*(- (s+2.d0*mt2*cone)*b0pdsmtmt
     +                  - b0dsmtmt + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mt2*b0pdsmtmt )
     +     +3.d0*( (gd(1)**2+gd(0)**2)*(- (s+2.d0*mb2*cone)*b0pdsmbmb
     +                  - b0dsmbmb + cone/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mb2*b0pdsmbmb ) )
* bosonic part
     +    + 1.d0/6.d0/sw2/cw2 * (
     +             + ( (18.d0*cw4+2.d0*cw2-0.5d0*cone)*s
     +                +(24.d0*cw4+16.d0*cw2-10.d0*cone)*mw2 )*b0pdsmwmw
     +             + (18.d0*cw4+2.d0*cw2-0.5d0*cone)*b0dsmwmw
     +             + (4.d0*cw2-cone)/3.d0 )
     +    + 1.d0/12.d0/sw2/cw2 * (
     +             (2.d0*mh2-10.d0*mz2-s)*b0pdsmzmh
     +            +(-1.d0)*b0dsmzmh
     +            -(mz2-mh2)**2/s*(b0pdsmzmh)
     +            +(mz2-mh2)**2/s/s*(b0dsmzmh-b0d0mzmh)
     +            -2.d0/3.d0*cone )
      szztpout  = - alsu4pi * szztpout

*
      end subroutine sigmazztp
*
******************************************************************
*
* eq. (B.4) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmawt(s,swtout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 swtout
*
      complex*16 b0ds0me,b0d0meme,b0d00me,
     +           b0ds0mm,b0d0mmmm,b0d00mm,
     +           b0ds0mtl,b0d0mtlmtl,b0d00mtl,
     +           b0dsmumd,b0d0mumu,b0d0mdmd,b0d0mumd,
     +           b0dsmcms,b0d0mcmc,b0d0msms,b0d0mcms,
     +           b0dsmtmb,b0d0mtmt,b0d0mbmb,b0d0mtmb,
     +           b0dsmwmg,b0d0mwmw,b0d0mwmg,b0dsmwmz,
     +           b0d0mzmz,b0dsmwmh,b0d0mhmh,b0d0mwmz,b0d0mwmh
* derivative for s = 0
      complex*16 b0pds0me,b0pds0mm,b0pds0mtl,
     +           b0pdsmumd,b0pdsmcms,b0pdsmtmb,
     +           b0pdsmwmg,b0pdsmwmz,b0pdsmwmh
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d0meme,b0d00me,b0d0mmmm,b0d00mm,b0d0mtlmtl,b0d00mtl,
     +     b0d0mumu,b0d0mdmd,b0d0mumd,b0d0mcmc,b0d0msms,b0d0mcms,
     +     b0d0mtmt,b0d0mbmb,b0d0mtmb,
     +     b0d0mwmg,b0d0mwmz,b0d0mwmw,b0d0mzmz,b0d0mwmh,b0d0mhmh

      if (ifirst.eq.0) then
          ifirst = 1

* leptons
          call b0p0mm(me2*cone,b0d0meme)
          call b0p0m0(me2*cone,b0d00me)

          call b0p0mm(mm2*cone,b0d0mmmm)
          call b0p0m0(mm2*cone,b0d00mm)

          call b0p0mm(mtl2*cone,b0d0mtlmtl)
          call b0p0m0(mtl2*cone,b0d00mtl)

* quarks
          call b0p0mm(mu2*cone,b0d0mumu)
          call b0p0mm(md2*cone,b0d0mdmd)
          call b0p0(mu2*cone,md2*cone,b0d0mumd)

          call b0p0mm(mc2*cone,b0d0mcmc)
          call b0p0mm(ms2*cone,b0d0msms)
          call b0p0(mc2*cone,ms2*cone,b0d0mcms)

          call b0p0mm(mt2*cone,b0d0mtmt)
          call b0p0mm(mb2*cone,b0d0mbmb)
          call b0p0(mt2*cone,mb2*cone,b0d0mtmb)

* bosons
          call b0p0m0(mw2,b0d0mwmg)

          call b0p0(mw2,mz2,b0d0mwmz)

          call b0p0mm(mw2,b0d0mwmw)

          call b0p0mm(mz2,b0d0mzmz)

          call b0p0(mw2,mh2,b0d0mwmh)

          call b0p0mm(mh2,b0d0mhmh)

      endif

* leptons
      call b0m0(s,me2*cone,b0ds0me)
      call b0m0(s,mm2*cone,b0ds0mm)
      call b0m0(s,mtl2*cone,b0ds0mtl)

* quarks
      call b0reg(s,mu2*cone,md2*cone,b0dsmumd)
      call b0reg(s,mc2*cone,ms2*cone,b0dsmcms)
      call b0reg(s,mt2*cone,mb2*cone,b0dsmtmb)
  
* bosons
      call b0m0(s,mw2,b0dsmwmg)
      call b0reg(s,mw2,mz2,b0dsmwmz)
      call b0reg(s,mw2,mh2,b0dsmwmh)

      if (abs(s).gt.epsilon) then

          swtout= 2.d0/3.d0/2.d0/sw2*(
* sum over three charged leptons
     +         -(s-me2/2.d0*cone)*b0ds0me + 1.d0/3.d0*s
     +         +me2*b0d0meme + me2*me2/2.d0/s*(b0ds0me-b0d00me)
     +         -(s-mm2/2.d0*cone)*b0ds0mm + 1.d0/3.d0*s
     +         +mm2*b0d0mmmm + mm2*mm2/2.d0/s*(b0ds0mm-b0d00mm)
     +         -(s-mtl2/2.d0*cone)*b0ds0mtl + 1.d0/3.d0*s
     +         +mtl2*b0d0mtlmtl + mtl2*mtl2/2.d0/s*(b0ds0mtl-b0d00mtl) )
* sum over three quark families
     +          + 2.d0/3.d0/2.d0/sw2 * 3.d0 * (
     +         -(s-(mu2+md2)*cone/2.d0)*b0dsmumd + 1.d0/3.d0*s
     +         +mu2*b0d0mumu + md2*b0d0mdmd 
     +         + (mu2-md2)**2/2.d0/s * (b0dsmumd - b0d0mumd)
     +         -(s-(mc2+ms2)*cone/2.d0)*b0dsmcms + 1.d0/3.d0*s
     +         +mc2*b0d0mcmc + ms2*b0d0msms 
     +         + (mc2-ms2)**2/2.d0/s * (b0dsmcms - b0d0mcms)
     +         -(s-(mt2+mb2)/2.d0)*b0dsmtmb + 1.d0/3.d0*s
     +         +mt2*b0d0mtmt + mb2*b0d0mbmb 
     +         + (mt2-mb2)**2/2.d0/s * (b0dsmtmb - b0d0mtmb) )
* bosonic part
     +          + 2.d0/3.d0*(
     +            (2.d0*mw2 + 5.d0*s)*b0dsmwmg - 2.d0*mw2*b0d0mwmw
     +           -mw2**2/s*(b0dsmwmg-b0d0mwmg) + 1.d0/3.d0*s )
     +          + 1.d0/12.d0/sw2*( ((40.d0*cw2-cone)*s 
     +                   +(16.d0*cw2+54.d0*cone-10.d0/cw2)*mw2)*b0dsmwmz
     +                   -(16.d0*cw2+2.d0*cone)*(mw2*b0d0mwmw
     +                                          +mz2*b0d0mzmz)
     +                   +(4.d0*cw2-cone)*2.d0/3.d0*s
     +                   -(8.d0*cw2+cone)*(mw2-mz2)**2/s
     +                              *(b0dsmwmz - b0d0mwmz) )
     +          + 1.d0/12.d0/sw2*( (2.d0*mh2-10.d0*mw2-s)*b0dsmwmh
     +                        -2.d0*mw2*b0d0mwmw - 2.d0*mh2*b0d0mhmh
     +                        -(mw2-mh2)**2/s*(b0dsmwmh-b0d0mwmh)
     +                        -2.d0/3.d0*s )
          swtout  = - alsu4pi * swtout

***

      else ! ( s = 0 )

* leptons
          call b0preg(s,zero,me2*cone,b0pds0me)
          call b0preg(s,zero,mm2*cone,b0pds0mm)
          call b0preg(s,zero,mtl2*cone,b0pds0mtl)

* quarks
          call b0preg(s,mu2*cone,md2*cone,b0pdsmumd)
          call b0preg(s,mc2*cone,ms2*cone,b0pdsmcms)
          call b0preg(s,mt2*cone,mb2*cone,b0pdsmtmb)

* bosons
          call b0preg(s,mw2,mz2,b0pdsmwmz)
          call b0preg(s,mw2,mh2,b0pdsmwmh)
          call b0preg(s,mw2,zero,b0pdsmwmg)

***
* sum over three charged leptons
          swtout= 2.d0/3.d0/2.d0/sw2*(
     +             - (-me2/2.d0) * b0ds0me
     +             + me2*b0d0meme + me2*me2/2.d0*b0pds0me
     +             - (-mm2/2.d0)*b0ds0mm
     +             + mm2*b0d0mmmm + mm2*mm2/2.d0*b0pds0mm
     +             - (-mtl2/2.d0)*b0ds0mtl 
     +             + mtl2*b0d0mtlmtl + mtl2*mtl2/2.d0*b0pds0mtl )
* sum over three quark families
     +          + 2.d0/3.d0/2.d0/sw2 * 3.d0 * (
     +             - (-(mu2+md2)/2.d0)*b0dsmumd
     +             + mu2*b0d0mumu + md2*b0d0mdmd 
     +             + (mu2-md2)**2/2.d0*b0pdsmumd
     +             - (-(mc2+ms2)/2.d0)*b0dsmcms
     +             + mc2*b0d0mcmc + ms2*b0d0msms 
     +             + (mc2-ms2)**2/2.d0*b0pdsmcms
     +             - (-(mt2+mb2)/2.d0)*b0dsmtmb 
     +             + mt2*b0d0mtmt + mb2*b0d0mbmb 
     +             + (mt2-mb2)**2/2.d0*b0pdsmtmb )
* bosonic part
     +          + 2.d0/3.d0*(
     +                    (2.d0*mw2)*b0dsmwmg - 2.d0*mw2*b0d0mwmw
     +                   - mw2**2*b0pdsmwmg  
     +            )
     +          + 1.d0/12.d0/sw2*(  
     +                   + (16.d0*cw2+54.d0*cone-10.d0/cw2)*mw2*b0dsmwmz
     +                   - (16.d0*cw2+2.d0*cone)*(mw2*b0d0mwmw
     +                                           +mz2*b0d0mzmz)
     +                   - (8.d0*cw2+cone)*(mw2-mz2)**2*b0pdsmwmz 
     +            )
     +          + 1.d0/12.d0/sw2*( (2.d0*mh2-10.d0*mw2)*b0dsmwmh
     +                        -2.d0*mw2*b0d0mwmw - 2.d0*mh2*b0d0mhmh
     +                        -(mw2-mh2)**2*b0pdsmwmh )
          swtout  = - alsu4pi * swtout

      endif

      end subroutine sigmawt
*
******************************************************************
*
* eq. (B.1) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmaaat(s,saatout,ifheavy)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 saatout
      logical ifheavy
*
      complex*16 b0dsmeme,b0dsmmmm,b0dsmtlmtl
      complex*16 b0d0meme,b0d0mmmm,b0d0mtlmtl

      complex*16 b0dsmumu,b0dsmcmc,b0dsmtmt
      complex*16 b0d0mumu,b0d0mcmc,b0d0mtmt

      complex*16 b0dsmdmd,b0dsmsms,b0dsmbmb
      complex*16 b0d0mdmd,b0d0msms,b0d0mbmb

      complex*16 b0dsmwmw
      complex*16 b0d0mwmw
*
      integer ifirst
      data ifirst/0/
      save ifirst
      save b0d0meme,b0d0mmmm,b0d0mtlmtl
      save b0d0mumu,b0d0mcmc,b0d0mtmt
      save b0d0mdmd,b0d0msms,b0d0mbmb
      save b0d0mwmw
*
      if (ifirst.eq.0) then
          ifirst = 1

          call b0p0mm(me2*cone,b0d0meme)
          call b0p0mm(mm2*cone,b0d0mmmm)
          call b0p0mm(mtl2*cone,b0d0mtlmtl)

          call b0p0mm(mu2*cone,b0d0mumu)
          call b0p0mm(mc2*cone,b0d0mcmc)
          call b0p0mm(mt2*cone,b0d0mtmt)

          call b0p0mm(md2*cone,b0d0mdmd)
          call b0p0mm(ms2*cone,b0d0msms)
          call b0p0mm(mb2*cone,b0d0mbmb)

          call b0p0mm(mw2,b0d0mwmw)
      endif

      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)

      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)

      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)

      call b0reg(s,mw2*cone,mw2*cone,b0dsmwmw)
*
      saatout = 2.d0/3.d0*(
* sum over three charged leptons
     +             2.d0*ql**2*( - (s+2.d0*me2*cone)*b0dsmeme
     +                          + 2.d0*me2*b0d0meme + s/3.d0
     +                          - (s+2.d0*mm2*cone)*b0dsmmmm
     +                          + 2.d0*mm2*b0d0mmmm + s/3.d0
     +                          - (s+2.d0*mtl2*cone)*b0dsmtlmtl
     +                          + 2.d0*mtl2*b0d0mtlmtl + s/3.d0)
* sum over quarks
     +  + 3.d0*(   2.d0*qu**2*( - (s+2.d0*mu2*cone)*b0dsmumu
     +                          + 2.d0*mu2*b0d0mumu + s/3.d0
     +                          - (s+2.d0*mc2*cone)*b0dsmcmc
     +                          + 2.d0*mc2*b0d0mcmc + s/3.d0)
     +           + 2.d0*qd**2*( - (s+2.d0*md2*cone)*b0dsmdmd
     +                          + 2.d0*md2*b0d0mdmd + s/3.d0
     +                          - (s+2.d0*ms2*cone)*b0dsmsms
     +                          + 2.d0*ms2*b0d0msms + s/3.d0
     +                          - (s+2.d0*mb2*cone)*b0dsmbmb
     +                          + 2.d0*mb2*b0d0mbmb + s/3.d0)
     +         ) 
     +         )

      if (ifheavy) 
     +    saatout = saatout + 2d0/3d0 * 
     +              3d0*2d0*qu**2*(
     +                 - (s+2.d0*mt2*cone)*b0dsmtmt
     +                 + 2.d0*mt2*b0d0mtmt + s/3.d0
     +              )
* bosonic part
     +      + ( 3.d0*s + 4.d0*mw2 )*b0dsmwmw - 4.d0*mw2*b0d0mwmw

      saatout = - alsu4pi * saatout
 
      end subroutine sigmaaat
******************************************************************
*
* derivative (w.r.t. p^2 (s)) of eq. (B.1) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmaaatp(s,saatpout,ifheavy)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 saatpout
      logical ifheavy
*
      complex*16 b0dsmeme,b0pdsmeme,
     +           b0dsmmmm,b0pdsmmmm,
     +           b0dsmtlmtl,b0pdsmtlmtl,
     +           b0dsmumu,b0pdsmumu,
     +           b0dsmdmd,b0pdsmdmd, 
     +           b0dsmcmc,b0pdsmcmc,
     +           b0dsmsms,b0pdsmsms,
     +           b0dsmtmt,b0pdsmtmt,
     +           b0dsmbmb,b0pdsmbmb,
     +           b0dsmwmw,b0pdsmwmw

* leptons
      call b0reg(s,me2*cone,me2*cone,b0dsmeme)
      call b0preg(s,me2*cone,me2*cone,b0pdsmeme)

      call b0reg(s,mm2*cone,mm2*cone,b0dsmmmm)
      call b0preg(s,mm2*cone,mm2*cone,b0pdsmmmm)

      call b0reg(s,mtl2*cone,mtl2*cone,b0dsmtlmtl)
      call b0preg(s,mtl2*cone,mtl2*cone,b0pdsmtlmtl)

* quarks
      call b0reg(s,mu2*cone,mu2*cone,b0dsmumu)
      call b0preg(s,mu2*cone,mu2*cone,b0pdsmumu)

      call b0reg(s,md2*cone,md2*cone,b0dsmdmd)
      call b0preg(s,md2*cone,md2*cone,b0pdsmdmd)

      call b0reg(s,mc2*cone,mc2*cone,b0dsmcmc)
      call b0preg(s,mc2*cone,mc2*cone,b0pdsmcmc)

      call b0reg(s,ms2*cone,ms2*cone,b0dsmsms)
      call b0preg(s,ms2*cone,ms2*cone,b0pdsmsms)

      call b0reg(s,mt2*cone,mt2*cone,b0dsmtmt)
      call b0preg(s,mt2*cone,mt2*cone,b0pdsmtmt)

      call b0reg(s,mb2*cone,mb2*cone,b0dsmbmb)
      call b0preg(s,mb2*cone,mb2*cone,b0pdsmbmb)

* bosons
      call b0reg(s,mw2,mw2,b0dsmwmw)
      call b0preg(s,mw2,mw2,b0pdsmwmw)

***

      saatpout= 2.d0/3.d0*(
* sum over three charged leptons
     +      2.d0*ql*ql*( (-1.d0)*b0dsmeme - (s+2.d0*me2*cone)*b0pdsmeme
     +                  + cone/3.d0 )
     +     +2.d0*ql*ql*( (-1.d0)*b0dsmmmm - (s+2.d0*mm2*cone)*b0pdsmmmm
     +                  + cone/3.d0 )
     +     +2.d0*ql*ql*( (-1.d0)*b0dsmtlmtl
     +                    -(s+2.d0*mtl2*cone)*b0pdsmtlmtl
     +                  + cone/3.d0 )
* sum over quarks
     +     +3.d0* 2.d0*qu*qu*( (-1.d0)*b0dsmumu 
     +                          - (s+2.d0*mu2*cone)*b0pdsmumu
     +                        + cone/3.d0 )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmdmd 
     +                          - (s+2.d0*md2*cone)*b0pdsmdmd
     +                        + cone/3.d0 )
     +     +3.d0* 2.d0*qu*qu*( (-1.d0)*b0dsmcmc 
     +                          - (s+2.d0*mc2*cone)*b0pdsmcmc
     +                        + cone/3.d0 )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmsms 
     +                          - (s+2.d0*ms2*cone)*b0pdsmsms
     +                        + cone/3.d0 )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmbmb 
     +                          - (s+2.d0*mb2*cone)*b0pdsmbmb
     +                        + cone/3.d0 ) )
      if (ifheavy)
     +    saatpout = saatpout + 2d0/3d0*
     +    3.d0*2.d0*qu*qu*( (-1.d0)*b0dsmtmt 
     +                          - (s+2.d0*mt2*cone)*b0pdsmtmt
     +                        + 1.d0/3.d0*cone )
* bosonic part
     +    +3.d0 * b0dsmwmw + (3.d0*s + 4.d0*mw2) * b0pdsmwmw

      saatpout  = - alsu4pi * saatpout

      end subroutine sigmaaatp
*
******************************************************************
*
* derivative (w.r.t. p^2 (s)) of eq. (B.4) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmawtp(s,swtpout)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
      include'pwhg_physpar.h'
*
      complex*16 s
      complex*16 swtpout
*
      complex*16 b0ds0me,b0d0meme,b0d00me,b0ds0mm,b0d0mmmm,b0d00mm,
     +           b0ds0mtl,b0d0mtlmtl,b0d00mtl,
     +           b0dsmumd,b0d0mumu,b0d0mdmd,b0d0mumd,
     +           b0dsmcms,b0d0mcmc,b0d0msms,b0d0mcms,
     +           b0dsmtmb,b0d0mtmt,b0d0mbmb,b0d0mtmb,
     +           b0dsmwmg,b0d0mwmw,b0d0mwmg,b0dsmwmz,
     +           b0d0mzmz,b0dsmwmh,b0d0mhmh,b0d0mwmz,b0d0mwmh
      complex*16 b0pds0me,b0pds0mm,b0pds0mtl,b0pdsmumd,
     +           b0pdsmcms,b0pdsmtmb,b0pdsmwmg,b0pdsmwmh,b0pdsmwmz
*
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d0meme,b0d00me,b0d0mmmm,b0d00mm,b0d0mtlmtl,b0d00mtl,
     +     b0d0mumu,b0d0mdmd,b0d0mumd,b0d0mcmc,b0d0msms,b0d0mcms,
     +     b0d0mtmt,b0d0mbmb,b0d0mtmb,
     +     b0d0mwmg,b0d0mwmz,b0d0mwmh,b0d0mwmw,b0d0mzmz,b0d0mhmh


      if (ifirst.eq.0) then
        ifirst = 1

* leptons
        call b0p0mm(me2*cone,b0d0meme)
        call b0p0m0(me2*cone,b0d00me)

        call b0p0mm(mm2*cone,b0d0mmmm)
        call b0p0m0(mm2*cone,b0d00mm)

        call b0p0mm(mtl2*cone,b0d0mtlmtl)
        call b0p0m0(mtl2*cone,b0d00mtl)

* quarks
        call b0p0mm(mu2*cone,b0d0mumu)
        call b0p0mm(md2*cone,b0d0mdmd)
        call b0p0(mu2*cone,md*cone,b0d0mumd)
   
        call b0p0mm(mc2*cone,b0d0mcmc)
        call b0p0mm(ms2*cone,b0d0msms)
        call b0p0(mc2*cone,ms2*cone,b0d0mcms)

        call b0p0mm(mt2*cone,b0d0mtmt)
        call b0p0mm(mb2*cone,b0d0mbmb)
        call b0p0(mt2*cone,mb2*cone,b0d0mtmb)

* bosons
        call b0p0m0(mw2,b0d0mwmg)
        call b0p0(mw2,mz2,b0d0mwmz)
        call b0p0(mw2,mh2,b0d0mwmh)
        call b0p0mm(mw2,b0d0mwmw)
        call b0p0mm(mz2,b0d0mzmz)
        call b0p0mm(mh2,b0d0mhmh)

      endif

* leptons
      call b0m0(s,me2*cone,b0ds0me)
      call b0pm10(s,me2*cone,b0pds0me)

      call b0m0(s,mm2*cone,b0ds0mm)
      call b0pm10(s,mm2*cone,b0pds0mm)

      call b0m0(s,mtl2*cone,b0ds0mtl)
      call b0pm10(s,mtl2*cone,b0pds0mtl)

* quarks
      call b0reg(s,mu2*cone,md2*cone,b0dsmumd)
      call b0preg(s,mu2*cone,md2*cone,b0pdsmumd)

      call b0reg(s,mc2*cone,ms2*cone,b0dsmcms)
      call b0preg(s,mc2*cone,ms2*cone,b0pdsmcms)

      call b0reg(s,mt2*cone,mb2*cone,b0dsmtmb)
      call b0preg(s,mt2*cone,mb2*cone,b0pdsmtmb)

* bosons
      call b0m0(s,mw2,b0dsmwmg)
      call b0reg(s,mw2,mz2,b0dsmwmz)
      call b0reg(s,mw2,mh2,b0dsmwmh)

      if(abs(dble(s-mw2)).lt.epsilon
     +   .and.abs(dimag(s-mw2)).lt.epsilon)then
          call b0pir(s,mw2,b0pdsmwmg)
      else
          call b0preg(s,mw2,zero,b0pdsmwmg)
      endif

      
      call b0preg(s,mw2,mz2,b0pdsmwmz)
      call b0preg(s,mw2,mh2,b0pdsmwmh)

***

      swtpout= 2.d0/3.d0/2.d0/sw2*(
* sum over three charged leptons
     +        -b0ds0me -(s-me2/2.d0*cone)*b0pds0me 
     +                        + cone/3.d0
     +        + me2**2/2.d0/s*(b0pds0me)
     +        - me2**2/2.d0/s/s*(b0ds0me-b0d00me)
     +        -b0ds0mm -(s-mm2/2.d0*cone)*b0pds0mm 
     +                        + cone/3.d0
     +        + mm2**2/2.d0/s*(b0pds0mm)
     +        - mm2**2/2.d0/s/s*(b0ds0mm-b0d00mm)
     +        -b0ds0mtl-(s-mtl2/2.d0*cone)*b0pds0mtl
     +                        + cone/3.d0
     +        + mtl2**2/2.d0/s*(b0pds0mtl)
     +        - mtl2**2/2.d0/s/s*(b0ds0mtl-b0d00mtl) )
* sum over three quark families
     +      + 2.d0/3.d0/2.d0/sw2 * 3.d0 * (
     +        -b0dsmumd -(s-(mu2+md2)*cone/2.d0)*b0pdsmumd
     +        + cone/3.d0
     +        + (mu2-md2)**2/2.d0/s * (b0pdsmumd)
     +        - (mu2-md2)**2/2.d0/s/s * (b0dsmumd - b0d0mumd)
     +        -b0dsmcms -(s-(mc2+ms2)*cone/2.d0)*b0pdsmcms
     +        + cone/3.d0
     +        + (mc2-ms2)**2/2.d0/s * (b0pdsmcms)
     +        - (mc2-ms2)**2/2.d0/s/s * (b0dsmcms - b0d0mcms)
     +        -b0dsmtmb -(s-(mt2+mb2)*cone/2.d0)*b0pdsmtmb
     +        + cone/3.d0
     +        + (mt2-mb2)**2/2.d0/s * (b0pdsmtmb)
     +        - (mt2-mb2)**2/2.d0/s/s * (b0dsmtmb - b0d0mtmb))
* bosonic part
     +        + 2.d0/3.d0*(
     +             (5.d0)*b0dsmwmg
     +           + (2.d0*mw2 + 5.d0*s)*b0pdsmwmg
     +           -mw2**2/s*(b0pdsmwmg)
     +           +mw2**2/s/s*(b0dsmwmg-b0d0mwmg)
     +           + cone/3.d0)
     +      +1.d0/12.d0/sw2*( ((40.d0*cw2-cone)*s 
     +                  +(16.d0*cw2+54.d0*cone-10.d0/cw2)*mw2)*b0pdsmwmz
     +                  +(40.d0*cw2-cone)*b0dsmwmz
     +                  +(4.d0*cw2 - cone)*2.d0/3.d0*cone
     +                  -(8.d0*cw2 + cone)*(mw2-mz2)**2/s
     +                                  *(b0pdsmwmz)
     +                  +(8.d0*cw2 + cone)*(mw2-mz2)**2/s/s
     +                              *(b0dsmwmz - b0d0mwmz) )
     +       +1.d0/12.d0/sw2*( (2.d0*mh2-10.d0*mw2-s)*b0pdsmwmh
     +                        +(-1.d0)*b0dsmwmh
     +                        -(mw2-mh2)**2/s*(b0pdsmwmh)
     +                        +(mw2-mh2)**2/s/s*(b0dsmwmh-b0d0mwmh)
     +                        -2.d0/3.d0*cone )
      swtpout= - alsu4pi * swtpout
      
      end subroutine sigmawtp
*
** end of vector boson self energies
*
****************************************************************
****************************************************************
*
** SCALAR FUNCTIONS WITH SPECIAL CASES 
*
*
* subroutine for regular b0
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0reg(p2i,m02i,m12i,b0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02i,m12i
      complex*16 b0out
      complex*16 p2,m02,m12
*
      complex*16 m0m1
      complex*16 m1sum0
      complex*16 r1,r2
      complex*16 a,b,c
      complex*16 cquad1,cquad2
      complex*16 beta
      complex*16 tmp
*
      p2 = dble(p2i)*cone

      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif
      if (abs(m02i).gt.abs(m12i)) then
          tmp = m02
          m02 = m12
          m12 = tmp
      endif
*
** p2 = 0
*
      if (abs(p2).lt.2d0*epsilon) then

          call b0p0(m02,m12,b0out)
*
** m02 = 0
*
      elseif (dble(m02).le.2d0*epsilon) then

          if (dble(m12).le.2d0*epsilon) then

              call b0m0m0(p2,b0out)
     
          else

              call b0m0(p2,m12,b0out)

          endif
*
** m02 = p2, m12 >> m02
*
      elseif (abs(p2-m02).le.2d0*epsilon.and.abs(p2/m12).lt.1.d-4) then

          m0m1=sqrt(m12*m02)
          m1sum0=m12/m02
 
          a = cone
          b = sqrt(m1sum0)
          c = cone
     
          call qcroots(a,b,c,r2,r1,cquad1,cquad2)
          r1 = -r1
          r2 = -r2

          b0out = (2.d0+log(mudim2))*cone - log(m0m1)
     +            + cone*log(m1sum0)/2.d0
     +            - b*( b*log(m1sum0)/2.d0 
     +                + (r2-r1)*log(r1) )

*
** m12 = m02
** 5.18 of B.P.
*
      elseif (abs(m12-m02).le.2d0*epsilon) then

          beta = sqrt( 1d0 + 4d0*m02/(-p2-ii*epsilon) )

          b0out = 2d0*cone + log(mudim2/m02)
     +           - beta*log( (beta+1d0) / (beta-1d0) )

*
** general case
*
      else

          m0m1=sqrt(m02*m12)
     
          a = cone
          b = (m02+m12-p2)/m0m1
          c = cone
     
          call qcroots(a,b,c,r2,r1,cquad1,cquad2)
          r2 = -r2
          r1 = -r1
     
          b0out  =  ( 2.d0 + log(mudim2) )*cone - log(m0m1)
     +            + (m02-m12)/p2*log(m12/m02)/2.d0 
     +            - m0m1/p2*(r2-r1)*log(r1)

      endif

      end subroutine b0reg
*
* subroutine for regular b0(p2,m0,m1) when m0 = m1 = 0
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0m0m0(p2i,b0out)
      implicit none
      include 'mathx.h'
*
      complex*16 p2i
      complex*16 b0out
      complex*16 p2
*
      p2=dble(p2i)*cone
*
      if(abs(p2).lt.2d0*epsilon) then

          b0out = zero

      else

          b0out = 2.d0*cone - log(-p2/mudim2)

      endif

      end subroutine b0m0m0

*
* subroutine for regular b0(p2,m0,m1) when p2 = m0 = 0
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0p0m0(m12i,b0out)
      implicit none
      include 'mathx.h'
*
      complex*16 m12i
      complex*16 b0out
      complex*16 m12
*
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif

      if(dble(m12).le.2d0*epsilon)then

          b0out = zero

      else

          b0out = cone - log(m12/mudim2)

      endif

      end subroutine b0p0m0
*
* subroutine for regular b0(p2,m0,m1) when p2 = 0
* and m1 != m2      
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0p0(m02i,m12i,b0out)
      implicit none
      include 'mathx.h'
*
      complex*16 m02i,m12i
      complex*16 b0out
      complex*16 m02,m12
*
      complex*16 tmp
*
      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif
      if (abs(m02i).gt.abs(m12i)) then
          tmp = m02
          m02 = m12
          m12 = tmp
      endif
*
      if (dble(m02).le.2d0*epsilon) then

          if (dble(m12).le.2d0*epsilon) then

              b0out   = zero
              return

          endif

          call b0p0m0(m12,b0out)
      
      elseif (abs(m02-m12).le.2d0*epsilon) then

          call b0p0mm(m02,b0out)

      else

          b0out = cone
     +            + (m12*log(m12/mudim2) - m02*log(m02/mudim2))
     +              /(m02-m12)

      endif

      end subroutine b0p0
*
* subroutine for regular b0(p2,m0,m1) when p2 = 0
* and m1 = m2
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0p0mm(m02i,b0out)
      implicit none
      include 'mathx.h'
*
      complex*16 m02i
      complex*16 b0out
      complex*16 m02
*
      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif

      b0out = log(mudim2/m02)

      end subroutine b0p0mm
*
* subroutine for regular b0(p2,m0,m1) when m0 = 0
* according to Denner hep-ph/0709.1075 eq. (4.23)
*
      subroutine b0m0(p2i,m12i,b0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      complex*16 p2i,m12i
      complex*16 b0out
      complex*16 p2,m12
      complex*16 arglog
*
      p2  = dble(p2i)*cone

      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif
*
      if(abs(p2).lt.2d0*epsilon) then

          call b0p0m0(m12,b0out)

      elseif (abs(m12-p2).le.2d0*epsilon) then

          b0out = 2.d0*cone - log(p2/mudim2)

      elseif (dble(m12).lt.abs(p2)) then

          arglog = ( p2 - m12 )/p2

          b0out = 2.d0*cone - log(p2/mudim2)
     +          - m12/p2*log(m12/p2) 
     +          - (cone-m12/p2)*log(arglog)
     +          + ii*pi*(cone-m12/p2)
      
      else

          arglog = -( p2 - m12 )/p2

          b0out = 2.d0*cone - log(p2/mudim2)
     +          - m12/p2*log(m12/p2) 
     +          - (cone-m12/p2)*log(arglog)

      endif
            
      end subroutine b0m0
*
* subroutine for the derivative of the b0
* according to Denner hep-ph/0709.1075 eq. (4.25)
*
      subroutine b0preg(p2i,m02i,m12i,b0pout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02i,m12i
      complex*16 b0pout
      complex*16 p2,m02,m12
*
      complex*16 m0m1
      complex*16 r1,r2
      complex*16 a,b,c
      complex*16 la
      complex*16 cquad1,cquad2
      complex*16 tmp
*
      p2 = dble(p2i)*cone

      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif
      if (abs(m02i).gt.abs(m12i)) then
          tmp = m02
          m02 = m12
          m12 = tmp
      endif

*
* special cases: 
* p2 = 0
*
      if (abs(p2).lt.2d0*epsilon) then
*
** p2=m02=m12=0
*
          if(dble(m02).le.2d0*epsilon.and.dble(m12).le.2d0*epsilon)then

              b0pout = zero
*
** p2=m12=0
*
          elseif (dble(m02).le.2d0*epsilon) then

              b0pout = 1.d0/(2.d0*m12)
*
** p2=0, m02=m12
*
          elseif (abs(m12-m02).le.2d0*epsilon) then

              b0pout = 1.d0/(6.d0*m02)

          else
*
** p2=0, m02 << m12
*
              if (abs(m02/m12).gt.1.d-4) then

                  b0pout =-(  (2.d0*m12*log(-m02) - m02)*m02
     +                      -(2.d0*m02*log(-m12) - m12)*m12 )/
     +                    ( 2.d0*(m02-m12)**3)
*
** p2=0, m02 << m12
*
              else

                  b0pout = +1.d0/(2.d0*m12)
     +                - m02*(-3.d0*cone - 2.d0*ii*pi - 2.d0*log(m02)
     +                  + 2.d0*log(-m12))/(2.d0*m12**2)
     +                + m02**2*(-5.d0 - 6.d0*ii*pi - 6.d0*log(m02)
     +                  + 6.d0*log(-m12))/(2.d0*m12**3)
     +                + m02**3*(-7.d0 - 12.d0*ii*pi - 12.d0*log(m02)
     +                  + 12.d0*log(-m12))/(2.d0*m12**4)
     +                + m02**4*(-9.d0 - 20.d0*ii*pi - 20.d0*log(m02)
     +                  + 20.d0*log(-m12))/(2.d0*m12**5)
     +                + m02**5*(-11.d0 - 30.d0*ii*pi - 30.d0*log(m02)
     +                  + 30.d0*log(-m12))/(2.d0*m12**6)

              endif

          endif

*
* m02 = 0
*
      elseif (dble(m02).le.2d0*epsilon) then
          
          call b0pm10(p2,m12,b0pout)
*
* p2 = m02, p2 << m12
*
      elseif (abs(p2-m02).le.2d0*epsilon.and.abs(p2/m12).lt.1.d-4) then

          a  = p2/m12
          la = log(a)
          b0pout =
     +       (a*(1.d0/2.*cone
     +      - a*((-11.d0/6.d0    -       la)
     +      - a*((-89.d0/12.d0   - 5.d0* la)
     +      - a*((-589.d0/20.d0  - 21.d0*la)
     +      - a*((-1732.d0/15.d0 - 84.d0*la) 
     +      - a*(-18899.d0/42.d0 - 330.d0*la) )))))
     +      )/p2
*
* general cases
*
      else

          m0m1=sqrt(m02*m12)
         
          a = cone
          b = ( m02 + m12 - p2 )/m0m1
          c = cone
         
          call qcroots(a,b,c,r1,r2,cquad1,cquad2)
          r1 = -r1
          r2 = -r2
        
          b0pout = -(m02-m12)/p2/p2*log(m12/m02)/2.d0
     +          + m0m1/p2/p2*(r2-r1)*log(r1)
     +          -( cone + (r1**2+cone)/(r1**2-cone)*log(r1) )/p2

      endif

      end subroutine b0preg
*
* subroutine for the derivative of the b0(p2,m02,m12)
* when m12 = 0
* according to Denner hep-ph/0709.1075 eq. (4.25)
*
      subroutine b0pm10(p2i,m02i,b0pout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      complex*16 p2i,m02i
      complex*16 b0pout
      complex*16 p2,m02
      complex*16 arglog
*
      p2 = dble(p2i)*cone

      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
*
* m02 = m12 = 0
*
      if(dble(m02).le.2d0*epsilon) then

          b0pout = -1d0/p2
*
* m12 = 0, p2 << m02
*
      elseif (abs(p2/m02).lt.1.d-4) then

          b0pout = 
     +         (  60*p2**5
     +           + m02*(70*p2**4 
     +           + m02*(84*p2**3 
     +           + m02*(105.d0*p2**2 
     +           + m02*(140.d0*p2 
     +           + 210.d0*m02 ) ))))
     +           /(420.d0*m02**6)

      else
*
* m12 = 0
*
          arglog = ( m02 - p2 )/m02

          if (dble(m02).le.dble(p2)) then
   
              b0pout = conjg(- m02*log(arglog)/(p2**2)-cone/p2)

          else

              b0pout = - m02*log(arglog)/(p2**2)-cone/p2

          endif
      
      endif

      end subroutine b0pm10

*
* subroutine for the regular b1
* according to Denner hep-ph/0709.1075 eq. (b.9)
*
      subroutine b1reg(p2i,m02i,m12i,b1out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02i,m12i
      complex*16 b1out
      complex*16 p2,m02,m12
*
      complex*16 b0out,b00out,b0pout
      complex*16 tmp
*
      p2  = dble(p2i)*cone

      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif
*
* p2 = 0 || p2 << m0/1
*
      if(abs(p2).lt.2d0*epsilon.and.dble(m02).le.2d0*epsilon
     +   .and.dble(m12).le.2d0*epsilon) then

         b1out = zero

      elseif(abs(p2)/(abs(m12)+epsilon).lt.1.d-5
     +   .and.abs(m02)/(abs(m12)+epsilon).lt.1.d-5) then

          b1out = - 0.25d0*cone + 0.5d0*log(m12/mudim2)

      elseif (abs(p2).lt.2d0*epsilon.or.
     +    abs(p2/(m12+epsilon)).lt.1.d-4.or.
     +    abs(p2/(m02+epsilon)).lt.1.d-4)then

          call b0preg(p2,m02,m12,b0pout)
          call b0p0(m02,m12,b00out)

          b1out = (m12-m02)/2.d0*b0pout - b00out/2.d0

*
* general case
*
      else

         call b0reg(p2,m02,m12,b0out)
         call b0p0(m02,m12,b00out)

         b1out = (m12-m02)/(2.d0*p2)*(b0out - b00out) - b0out/2.d0

      endif

      end subroutine b1reg
*
* subroutine for the regular b1
* according to Denner hep-ph/0709.1075 eq. (b.10)
*
      subroutine b1preg(p2i,m02i,m12i,b1pout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02i,m12i
      complex*16 b1pout
      complex*16 p2,m02,m12
*
      complex*16 b0out,b00out,b0pout
*
      p2  = dble(p2i)*cone
      if (abs(dimag(m02i)).lt.epsilon) then
          m02 = m02i - ii*epsilon
      else
          m02 = m02i
      endif
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif

*
      if (abs(p2).lt.2d0*epsilon) then
*
** p2 = m02 = m12 = 0
*
          if(dble(m02).le.2d0*epsilon.and.dble(m12).le.2d0*epsilon)then

              b1pout = zero
*
** p2 = m02 = 0
*
          elseif (dble(m02).le.2d0*epsilon) then

              b1pout = -1.d0/(6.d0*m12)
*
** p2 = m12 = 0
*
          elseif (dble(m12).le.2d0*epsilon) then

              b1pout = -1.d0/(6.d0*m02)
*
** p2 = 0
*
          else

              b1pout =-(2*m02**3 + 3*m02**2*m12 - 6*m02*m12**2 + m12**3+
     +                  6*m02**2*m12*(-log(-m02) + log(-m12)))/
     +                 (6.*(m02 - m12)**4)

          endif

      elseif (dble(m12).le.2d0*epsilon) then
*
** m12 = 0, p2 << m02
*
          if (abs(p2/m02).lt.1.d-4) then

              b1pout = 1/(3.*m02) + p2/(4.*m02**2) + p2**2/(5.*m02**3) +
     +                 p2**3/(6.*m02**4) + p2**4/(7.*m02**5) + 
     +                 p2**5/(8.*m02**6) + p2**6/(9.*m02**7) + 
     +                 p2**7/(10.*m02**8) + p2**8/(11.*m02**9) + 
     +                 p2**9/(12.*m02**10) + p2**10/(13.*m02**11)
*
** m12 = 0
*
          else

              b1pout = - ( 2.d0*m02**2*(log(m02)-log(m02-p2))
     +                                  -p2*(2.d0*m02+p2)   )
     +                   /(2.d0*p2**3)

          endif

      elseif (dble(m02).le.2d0*epsilon) then
*
** m02 = 0, p2 << m12
*
          if (abs(p2/m12).lt.1.d-4) then

              b1pout = - ( 60060*m12**10 + 30030*m12**9*p2 
     +              + 18018*m12**8*p2**2 + 12012*m12**7*p2**3 
     +              + 8580*m12**6*p2**4 + 6435*m12**5*p2**5 
     +              + 5005*m12**4*p2**6 + 4004*m12**3*p2**7 
     +              + 3276*m12**2*p2**8 + 2730*m12*p2**9 + 2310*p2**10 )
     +              /(360360.*m12**11)
*
** m02 = 0
*
          else

              b1pout = - ( p2*(2.d0*m12-p2)
     +                    +2.d0*m12*(m12-p2)*(log(p2-m12)-log(-m12)))
     +                  /(2.d0*p2**3)

          endif

      elseif (abs(m02-p2).le.2d0*epsilon) then
*
** p2 = m02 << m12
*
          if (abs(p2/m12).lt.1.d-4) then

              b1pout = (-4620*m12**9 + 6930*m12**8*p2 + 
     +        35574*m12**7*p2**2 + 64218*m12**6*p2**3 + 
     +        198924*m12**5*p2**4 + 557073*m12**4*p2**5 + 
     +        1561439*m12**3*p2**6 + 4341711*m12**2*p2**7 + 
     +        12012782*m12*p2**8 + 33128851*p2**9 + 
     +        (0,27720)*p2**2*
     +         (m12**7 + 3*m12**6*p2 + 8*m12**5*p2**2 + 
     +           22*m12**4*p2**3 + 60*m12**3*p2**4 + 
     +           164*m12**2*p2**5 + 448*m12*p2**6 + 1224*p2**7)*pi
     +        + 27720*p2**2*
     +          (m12**7 + 3*m12**6*p2 + 8*m12**5*p2**2 + 
     +           22*m12**4*p2**3 + 60*m12**3*p2**4 + 
     +           164*m12**2*p2**5 + 448*m12*p2**6 + 1224*p2**7)*
     +         (-log(-m12) + log(p2)))/(27720.*m12**10)
*
** p2 = m02
*
          else

              b1pout = -((2*m12 - 3*p2)*p2 - 2*(m12 - 2*p2)*(m12 - p2)*
     +                 (log(-m12 + p2) - log(-m12 + 2*p2)))
     +                /(2.*p2**3)

          endif

      elseif (abs(m02)/(abs(m12)+epsilon).lt.1.d-5
     +        .and.abs(p2)/(abs(m12)+epsilon).lt.1.d-5) then

          b1pout = -1d0/(6d0*m12)
*
** general case
*
      else

          call b0reg(p2,m02,m12,b0out)
          call b0p0(m02,m12,b00out)
          call b0preg(p2,m02,m12,b0pout)

          b1pout  = - 0.5d0*( - (m12-m02)/(2.d0*p2**2)*(b0out-b00out)
     +                        + (m12-m02-p2)/(2.d0*p2)*b0pout )

      endif

      end subroutine b1preg
*
** subroutine for the singular b0p (p2,lambda,m12)
*
      subroutine b0pir(p2i,m12i,b0pout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m12i
      complex*16 b0pout
      complex*16 p2,m12
*
      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif

      if (dble(m12).le.2d0*epsilon) then
          b0pout = zero
          return
      endif

      p2  = dble(p2i)*cone
*
* eq. 5.51 bardin - passarino
*
      b0pout = - ( 2.d0*cone - log(m12/mudim2) )/(2d0*m12)

      end subroutine b0pir
*
** subroutine for the singular b1p (p2,m12,0d0)
*
      subroutine b1pir(p2i,m12i,b1pout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m12i
      complex*16 b1pout
      complex*16 p2,m12
*
      p2 = dble(p2i)*cone

      if (abs(dimag(m12i)).lt.epsilon) then
          m12 = m12i - ii*epsilon
      else
          m12 = m12i
      endif

*
      if(abs(p2).le.2d0*epsilon)then

          b1pout = zero

      else
*
* eq. 5.54 bardin - passarino
*
          b1pout = + ( ( 3.d0 )*cone - log(m12/mudim2) )/(2.d0*m12)

      endif

      end subroutine b1pir
*
* subroutine for the singular collinear C0 (and also threshold, 
* cured with complex W mass)
* according to Dittmaier-Kraemer PRD65 073007 
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (b.2)
*
* the sequence of external p2 and internal masses is understood as 
* m1^2,m1^2,s,0,m1^2,mw^2, n.b.: where m1 is the lepton mass
* this is not truely collinear singular because there is the lepton mass
*
* here mz used is complex
*
* from non abelian diagram
*
      subroutine c0cl(m12in,s_in,m22in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real*8 m12in,s_in
      complex*16 m22in
      complex*16 c0out
      real*8 m12,s
*
      complex*16 myli2
      external myli2
*
      complex*16 m22,s2bar
      complex*16 arglog1,arglog2
      complex*16 argli21,argli22
*
      m12= m12in
      s  = s_in
      m22= m22in

      s2bar = s*cone + ii*epsilon
*
* eq. (B.2) Dittmaier
*

      arglog1= ( m22 - s2bar )/m12
      arglog2= ( m22 - s2bar )/m22
      argli21= - s2bar/( m22 - s2bar )

      argli22 = s2bar/m22

      c0out = - log(arglog1)*log(arglog2)
     +         - 2.d0*myli2(argli21)
     +         - myli2(argli22)
      c0out = c0out/(-s)
*
      end subroutine c0cl

*
* subroutine for the singular collinear C0 (and also threshold, 
* cured with complex W mass)
* according to Dittmaier-Kraemer PRD65 073007 
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (B.2) and eq. (B.4)
*
* the sequence of external p2 and internal masses is understood as 
* m1^2,0,s,0,m1^2,mw^2, n.b.: where m1 is the quark mass
*
* this is truely collinear singular since the quark mass -> 0
* here mw used is complex
*
* from non-abelian diagram
*
*
      subroutine c0cq(m12in,s_in,m22in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real*8 m12in,s_in
      complex*16 m22in
      complex*16 c0out
      real*8 m12,s
*
      complex*16 myli2
      external myli2
*
      complex*16 m22,s2bar
      complex*16 arglog
*
      m12= m12in - ii*epsilon
      s  = s_in
      m22= m22in - ii*epsilon

      s2bar= s*cone + ii*epsilon
*
      arglog = ( m22 - s2bar )/m22

      c0out= - log(mudim2/m22) * log(arglog)
     +       + ( log(arglog) )**2
     +       + myli2(s2bar/m22)

      c0out= c0out/s
*
      end subroutine c0cq

* subroutine for the singular ir-coll C0
* according to Dittmaier-Kraemer PRD65 073007
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (B.5) and eq. (B.12)
*
* the sequence of external p2 and internal masses is understood as 
* p102in,p202in,p212in,p102in,mg2in,p202in, where p102in= mq^2, p202in= mlep^2
*
      subroutine c0ircoll_ql(p102in,p202in,p212in,mg2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      real*8 p102in,p202in,p212in
      real*8 mg2in
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      complex*16 m22,r,argli21
*
      m22= p202in - ii*epsilon
      r= p212in + ii*epsilon
*
* eq. (B.12) Dittmaier
*
* it is understood that m1= m_quark and m2= m_lepton
*
      argli21= r/(r-m22)

      c0out =  pi2/12d0
     +       - ( log(mudim2/m22) )**2 / 4d0
     +       + ( log(mudim2/(m22-r)) )**2 / 2d0
      c0out = c0out - myli2(argli21)
      c0out = c0out/(r - m22)
*
      end subroutine c0ircoll_ql

*
* subroutine for the singular ir-coll C0
* according to Dittmaier-Kraemer PRD65 073007 
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (B.5) and eq. (B.16)
*
* the sequence of external p2 and internal masses is understood as 
* p102in,p202in,p212in,p102in,mg2in,p202in, where p102in= mq^2, p202in= mqp^2
*
      subroutine c0ircoll_qqp(p102in,p202in,p212in,mg2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p102in,p202in,p212in
      complex*16 mg2in
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      real*8 r
      complex*16 rbar
*
      r   = p212in
      rbar= r + ii*epsilon
*
* eq. (B.16) Dittmaier
*
      c0out = ( log(-mudim2/rbar) )**2 / 2d0 / r
*
      end subroutine c0ircoll_qqp

*
* subroutine for the singular ir-coll C0
* according to Dittmaier NPB 675
* eq. (B.8/9) and eq. (B.12)
*
* the sequence of external p2 and internal masses is understood as 
* p102in,p202in,p212in,p102in,mg2in,p202in, where p102in= ml^2, p202in= mq^2
*
* mq -> 0
*
      subroutine c0lq(p102in,p202in,p212in,mg2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p102in,p202in,p212in,mg2in
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      complex*16 r,ml2
      complex*16 rbar

      r = p212in*cone 
      rbar = r + ii*epsilon

      ml2  = p102in*cone - ii*epsilon

      c0out = - myli2(rbar/(rbar-ml2))
     -        + pi2/12d0 - (log(mudim2/ml2))**2/4d0
     -        + (log(mudim2/(ml2-rbar)))**2/2d0

      c0out = c0out/(r-ml2)
*
      end subroutine c0lq
*
* subroutine for the singular ir-coll C0
* according to Dittmaier NPB 675
* eq. (B.2) and eq. (B.4)
*
* the sequence of external p2 and internal masses is understood as 
* s1,mf2,s2,mv2,mf2,0
* mf2 is lambda2 (soft scale)
*
      subroutine c0fz(s1,s2,mf2,mv2in,c0out)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
*
      complex*16 mf2,mv2in,s1,s2
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      complex*16 s1bar,s2bar,mv2

      s1bar = s1 + ii*epsilon
      s2bar = s2 + ii*epsilon

      if(abs(dimag(mv2in)).lt.epsilon) then
          mv2 = mv2in - ii*epsilon
      else
          mv2 = mv2in
      endif

      if (abs(mf2).gt.epsilon) then

          c0out = (   log( (mv2-s1bar)/mf2 )*log( (mv2-s1bar)/mv2 )
     +              - log( (mv2-s2bar)/mf2 )*log( (mv2-s2bar)/mv2 )
     +              - 2d0*myli2( (s1bar-s2bar)/(mv2-s2bar) )
     +              + myli2(s1bar/mv2) - myli2(s2bar/mv2)   )/(s1 - s2)

      else

          c0out =( (log((mv2-s1bar)/mv2))**2 - (log((mv2-s2bar)/mv2))**2
     +            + myli2(s1bar/mv2) - myli2(s2bar/mv2) 
     +            + log(mudim2/mv2)*log((mv2-s2bar)/(mv2-s1bar)) )
     +           /(s1-s2)

      endif
*
      end subroutine c0fz
*
* subroutine for the infrared c0
* according to Dittmaier NPB 675
* eq. (B.9) and eq. (B.6)
*
* the sequence of external p2 and internal masses is understood as 
* ml2,mlp2,s,ml2,lambda2,mlp2
*
      subroutine c0ir_ffp(s,ml2in,mlp2in,c0out)
      implicit none
      include'mathx.h'
      include'pwhg_math.h'
*
      complex*16 ml2in,mlp2in
      complex*16 s,ml2,mlp2
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      complex*16 sbar,m1,m2
      complex*16 xs

      if (abs(ml2in).lt.epsilon.and.abs(mlp2in).lt.epsilon) then
          call c0ircoll_qqp(ml2in,mlp2in,s,zero,c0out)
          return
      endif

      sbar = s + ii*epsilon

      ml2  = ml2in  - ii*epsilon
      mlp2 = mlp2in - ii*epsilon

      m1 = sqrt(ml2)
      m2 = sqrt(mlp2)

      if(abs(m1*m2/s).lt.1.d-6) then
         xs= -m1*m2/s*cone
      else
         xs= (sqrt(cone
     +        -4.d0*m1*m2/(s*cone+ii*epsilon
     +                           -cone*(m1-m2)**2))-cone)/
     +       (sqrt(cone
     +        -4.d0*m1*m2/(s*cone+ii*epsilon
     +                           -cone*(m1-m2)**2))+cone)
      endif
*
      c0out = - 0.5d0*(log(xs))**2 + 2d0*log(xs)*log(cone-xs*xs)
     +        + 0.5d0*(log(m1/m2))**2 - pi2/6d0 + myli2(xs*xs)
     +        + myli2(cone-xs*m1/m2) + myli2(cone - xs*m2/m1)
     +        - log(mudim2/m1/m2)*log(xs)

      c0out = c0out*xs/( m1*m2*(cone-xs*xs) )
*
      end subroutine c0ir_ffp
*
* subroutine for the regular C0 with generic internal masses (also complex)
* according to Denner arXiv:0709.1075 (Fortscrhitte), 
* Eqs. (4.26), (4.27), (4.28)
*
      subroutine c0reg(p102,p212,p202,
     +                 m02in,m12in,m22in,
     +                 c0out)
      implicit none
*
      real*8 p102,p202,p212
      complex*16 m02in,m12in,m22in
      complex*16 m02,m12,m22
      complex*16 c0out
*
      real*8 pi,pis
      common/pigreco/pi,pis

      real*8 epsilon
      common/epsi/epsilon

      complex*16 cone,zero,ii
      common/complexunit/cone,zero,ii
*
      complex*16 myli2
      external myli2

      complex*16 myeta
      external myeta

      complex*16 kallen
      external kallen
*
      complex*16 alpha
      complex*16 arg1,arg2,arg3
      complex*16 y00,y01,y02
      complex*16 y10,y20
      complex*16 alpha0,alpha1,alpha2
      complex*16 x0p,x0m,x1p,x1m,x2p,x2m
      complex*16 y0p,y0m,y1p,y1m,y2p,y2m
      complex*16 argli0p1,argli0p2,argli0m1,argli0m2
      complex*16 argli1p1,argli1p2,argli1m1,argli1m2
      complex*16 argli2p1,argli2p2,argli2m1,argli2m2
      complex*16 eta0p1,eta0p2,eta0m1,eta0m2
      complex*16 eta1p1,eta1p2,eta1m1,eta1m2
      complex*16 eta2p1,eta2p2,eta2m1,eta2m2
      complex*16 eta03,eta04
      real*8 theta21,theta0
      complex*16 eta13,eta14
      real*8 theta20,theta1
      complex*16 eta23,eta24
      real*8 theta10,theta2
*
      m02= m02in
      m12= m12in
      m22= m22in
*
      c0out= zero
*
* Eq. (4.27)
*
      alpha= kallen(p102*cone,p212*cone,p202*cone)
* i=0, j=1, k=2
      y00= 1.d0/2.d0/alpha/p212*(p212*( (p212-p202-p102)*cone 
     +                                 +2.d0*m02-m12-m22 )
     +                                 -(p202-p102)*(m12-m22)
     +                                 +alpha*(p212*cone-m12+m22))
      arg1= p212 * cone
      arg2= m12  * cone
      arg3= m22  * cone
      alpha0= kallen(arg1,arg2,arg3)*(cone + ii*epsilon*p212)
      x0p= 1.d0/2.d0/p212*(p212*cone - m12 + m22 + alpha0)
      x0m= 1.d0/2.d0/p212*(p212*cone - m12 + m22 - alpha0)
      y0p= y00 - x0p
      y0m= y00 - x0m
* i=1, j=2, k=0
      y01= 1.d0/2.d0/alpha/p202*(p202*( (p202-p102-p212)*cone 
     +                                 +2.d0*m12-m22-m02 )
     +                                 -(p102-p212)*(m22-m02)
     +                                 +alpha*(p202*cone-m22+m02))
      arg1= p202 * cone
      arg2= m22  * cone
      arg3= m02  * cone
      alpha1= kallen(arg1,arg2,arg3)*(cone + ii*epsilon*p202)
      x1p= 1.d0/2.d0/p202*(p202*cone - m22 + m02 + alpha1)
      x1m= 1.d0/2.d0/p202*(p202*cone - m22 + m02 - alpha1)
      y1p= y01 - x1p
      y1m= y01 - x1m
* i=2, j=0, k=1
      y02= 1.d0/2.d0/alpha/p102*(p102*( (p102-p212-p202)*cone 
     +                                 +2.d0*m22-m02-m12 )
     +                                 -(p212-p202)*(m02-m12)
     +                                 +alpha*(p102*cone-m02+m12))
      arg1= p102 * cone
      arg2= m02  * cone
      arg3= m12  * cone
      alpha2= kallen(arg1,arg2,arg3)*(cone + ii*epsilon*p102)
      x2p= 1.d0/2.d0/p102*(p102*cone - m02 + m12 + alpha2)
      x2m= 1.d0/2.d0/p102*(p102*cone - m02 + m12 - alpha2)
      y2p= y02 - x2p
      y2m= y02 - x2m
* i=1, j=2, k=0
      y10= y01
* i=2, j=0, k=1
      y20= y02
*
* Building blocks of Eq. (4.26)
*
*
* i= 0, sigma= p,m
*
      argli0p1= (y00-cone)/y0p
      argli0p2= y00/y0p

      argli0m1= (y00-cone)/y0m
      argli0m2= y00/y0m

      arg1= cone - x0p
      arg2= 1.d0/y0p
      eta0p1 = myeta(arg1,arg2)

      arg1= -x0p
      arg2= 1.d0/y0p
      eta0p2 = myeta(arg1,arg2)

      arg1= cone - x0m
      arg2= 1.d0/y0m
      eta0m1 = myeta(arg1,arg2)

      arg1= -x0m
      arg2= 1.d0/y0m
      eta0m2 = myeta(arg1,arg2)
*
* i= 1, sigma= p,m
*
      argli1p1= (y01-cone)/y1p
      argli1p2= y01/y1p

      argli1m1= (y01-cone)/y1m
      argli1m2= y01/y1m
      
      arg1= cone - x1p
      arg2= 1.d0/y1p
      eta1p1 = myeta(arg1,arg2)

      arg1= -x1p
      arg2= 1.d0/y1p
      eta1p2 = myeta(arg1,arg2)

      arg1= cone - x1m
      arg2= 1.d0/y1m
      eta1m1 = myeta(arg1,arg2)

      arg1= -x1m
      arg2= 1.d0/y1m
      eta1m2 = myeta(arg1,arg2)
*
* i= 2, sigma= p,m
*
      argli2p1= (y02-cone)/y2p
      argli2p2= y02/y2p

      argli2m1= (y02-cone)/y2m
      argli2m2= y02/y2m

      arg1= cone - x2p
      arg2= 1.d0/y2p
      eta2p1 = myeta(arg1,arg2)

      arg1= -x2p
      arg2= 1.d0/y2p
      eta2p2 = myeta(arg1,arg2)

      arg1= cone - x2m
      arg2= 1.d0/y2m
      eta2m1 = myeta(arg1,arg2)

      arg1= -x2m
      arg2= 1.d0/y2m
      eta2m2 = myeta(arg1,arg2)
*
* last two eta functions independent of sigma
*
*
* i= 0
*
      arg1= -x0p
      arg2= -x0m
      eta03 = myeta(arg1,arg2)

      arg1= y0p
      arg2= y0m
      eta04 = myeta(arg1,arg2)

      theta21= 0.d0
      if(-p212.gt.0.d0) theta21= 1.d0
      theta0= 0.d0
      if(-dimag(y0p*y0m).gt.0.d0) theta0= 1.d0
*
* i= 1
*
      arg1= -x1p
      arg2= -x1m
      eta13 = myeta(arg1,arg2)

      arg1= y1p
      arg2= y1m
      eta14 = myeta(arg1,arg2)

      theta20= 0.d0
      if(-p202.gt.0.d0) theta20= 1.d0
      theta1= 0.d0
      if(-dimag(y1p*y1m).gt.0.d0) theta1= 1.d0
*
* i= 2
*
      arg1= -x2p
      arg2= -x2m
      eta23 = myeta(arg1,arg2)

      arg1= y2p
      arg2= y2m
      eta24 = myeta(arg1,arg2)

      theta10= 0.d0
      if(-p102.gt.0.d0) theta10= 1.d0
      theta2= 0.d0
      if(-dimag(y2p*y2m).gt.0.d0) theta2= 1.d0
*
* Eq. (4.26)
*
      c0out= zero
* i= 0
     +     + myli2(argli0p1) - myli2(argli0p2)
     +     + eta0p1*log(argli0p1)
     +     - eta0p2*log(argli0p2)
     +     + myli2(argli0m1) - myli2(argli0m2)
     +     + eta0m1*log(argli0m1)
     +     - eta0m2*log(argli0m2)
     +     -(eta03-eta04-2.d0*pi*(0.d0,1.d0)*theta21*theta0)
     +                                       *log((cone-y00)/(-y00))
* i= 1
     +     + myli2(argli1p1) - myli2(argli1p2)
     +     + eta1p1*log(argli1p1)
     +     - eta1p2*log(argli1p2)
     +     + myli2(argli1m1) - myli2(argli1m2)
     +     + eta1m1*log(argli1m1)
     +     - eta1m2*log(argli1m2)
     +     -(eta13-eta14-2.d0*pi*(0.d0,1.d0)*theta20*theta1)
     +                                       *log((cone-y10)/(-y10))
* i= 2
     +     + myli2(argli2p1) - myli2(argli2p2)
     +     + eta2p1*log(argli2p1)
     +     - eta2p2*log(argli2p2)
     +     + myli2(argli2m1) - myli2(argli2m2)
     +     + eta2m1*log(argli2m1)
     +     - eta2m2*log(argli2m2)
     +     -(eta23-eta24-2.d0*pi*(0.d0,1.d0)*theta10*theta2)
     +                                       *log((cone-y20)/(-y20))
      c0out= 1.d0/alpha*c0out

      end subroutine c0reg

*
* eq. (5.66) of Bardin - Passarino
* C0(0,0,p2,m02,0,m12)
*
      subroutine c0reg1(p2i,m02in,m12in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02in,m12in
      complex*16 p2,m02,m12
      complex*16 c0out
*
      complex*16 kallen
      external kallen
*
      complex*16 x1,x2

* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone

      if (abs(dimag(m02in)).lt.epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif
      if (abs(dimag(m12in)).lt.epsilon) then
          m12 = m12in - ii*epsilon
      else
          m12 = m12in
      endif

      x1 = (p2+m02-m12-kallen(-p2,m02,m12))/(2.d0*p2)
      x2 = (p2+m02-m12+kallen(-p2,m02,m12))/(2.d0*p2)

      c0out  = -log(x2/(x2-cone))*log((x1-cone)/x1) / p2

      end subroutine c0reg1
*
* eq. (5.59) of Bardin - Passarino
* C0(0,0,p2,0,m2,0)
*
      subroutine c0reg2(p2i,m2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      complex*16 p2i,m2in
      complex*16 p2,m2
      complex*16 c0out
*
      complex*16 myli2
      external myli2

      if (abs(dimag(m2in)).lt.epsilon) then
          m2 = m2in - ii*epsilon
      else
          m2 = m2in
      endif
*
* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone
*
      c0out = - ( pi2/6d0 - myli2(cone-(p2/m2)) ) / p2

      end subroutine c0reg2

*
* eq. (5.67) of Bardin - Passarino
* C0(0,0,p2,m2,0,m2)
*
      subroutine c0reg3(p2i,m2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m2in
      complex*16 p2,m2
      complex*16 c0out
*
      complex*16 kallen
      external kallen
*
      complex*16 betaq

* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone

      if (abs(dimag(m2in)).lt.epsilon) then
          m2 = m2in - ii*epsilon
      else
          m2 = m2in
      endif

      betaq = sqrt(cone + 4d0*m2/p2)

      c0out  = - ( log((betaq+cone)/(betaq-cone)) )**2 / p2

      end subroutine c0reg3
*
* eq. (5.68) of Bardin - Passarino
* C0(0,0,p2,0,m02,m12)
*
      subroutine c0reg4(p2i,m02in,m12in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02in,m12in
      complex*16 p2,m02,m12
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      if (abs(m12in).lt.epsilon) then
          call c0reg2(p2i,m02in,c0out)
          return
      endif
* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone

      if (abs(dimag(m02in)).lt.epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif
      if (abs(dimag(m12in)).lt.epsilon) then
          m12 = m12in - ii*epsilon
      else
          m12 = m12in
      endif


      c0out  = - ( myli2(cone-m12/m02)-myli2(cone-(p2+m12)/m02) ) / p2

      end subroutine c0reg4
*
* eq. (5.64) of Bardin-Passarino
* C0(0,0,p2,m02,m12,m22)
*
      subroutine c0reg5(p2i,m02in,m12in,m22in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m02in,m12in,m22in
      complex*16 p2,m02,m12,m22
      complex*16 c0out
*
      complex*16 myli2
      external myli2

      complex*16 kallen
      external kallen
*
      complex*16 x0,x1,x2,x3

* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone

      if (abs(dimag(m02in)).lt.epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif
      if (abs(dimag(m12in)).lt.epsilon) then
          m12 = m12in - ii*epsilon
      else
          m12 = m12in
      endif
      if (abs(dimag(m22in)).lt.epsilon) then
          m22 = m22in - ii*epsilon
      else
          m22 = m22in
      endif

* eq. (5.65) of Bardin - Passarino

      x0 = cone + (m02-m12)/p2
      x1 = ( p2 + m02 - m22 - kallen(-p2,m02,m22) )/(2d0*p2)
      x2 = ( p2 + m02 - m22 + kallen(-p2,m02,m22) )/(2d0*p2)
      x3 = m22/(m22-m12)

      c0out  = - ( + myli2((x0-cone)/(x0-x1)) - myli2(x0/(x0-x1))
     +             + myli2((x0-cone)/(x0-x2)) - myli2(x0/(x0-x2))
     +             - myli2((x0-cone)/(x0-x3)) + myli2(x0/(x0-x3))
     +           ) /p2

      end subroutine c0reg5
*
* eq. (5.69) of Bardin-Passarino
* C0(m2,m2,p2,0,m2,0)
*
      subroutine c0reg6(p2i,m2in,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m2in
      complex*16 p2,m2
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      complex*16 y1,y2

      if (abs(m2in).lt.epsilon) then
          call c0ircoll_qqp(m2in,m2in,p2i,zero,c0out)
          return
      endif

* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone

      if (abs(dimag(m2in)).lt.epsilon) then
          m2 = m2in - ii*epsilon
      else
          m2 = m2in
      endif

      y1 = -p2/(2d0*m2)*( cone + sqrt(cone+4d0*m2/p2) )
      y2 = -p2/(2d0*m2)*( cone - sqrt(cone+4d0*m2/p2) )

      c0out  = -1d0/(m2*(y1-y2))*(
     +           + 2d0*myli2(1d0/y1) - 2d0*myli2(1d0/y2)
     +           + myli2(y1) - myli2(y2)
     +         )

      end subroutine c0reg6
*
* subroutine for the D0 with 2 photons
* according to Denner-Dittmaier ArXiv:1005.2076, 
* eqs.(4.11), (4.45)
*
* typical entries: p102= mq2, p212= mq2, p322= mlep2, p302= mlep2,
*                  p202= s, p312= t, 
*                  m02= lambda2, m12= mq2, m22= lambda2, m32= mlep2
*
      subroutine d02g(p102,p212,p322,p302,p202,p312,
     +                m02,m12in,m22,m32in,
     +                boxout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12in,m22,m32in
      complex*16 boxout
*
      complex*16 m12,m32
      complex*16 m1,m3
      complex*16 p202b,p312b,x31

      p202b = p202 + ii*epsilon
      p312b = p312 + ii*epsilon

      m12 = m12in - ii*epsilon
      m32 = m32in - ii*epsilon

      m1 = sqrt(m12)
      m3 = sqrt(m32)

      boxout = 2d0*log(mudim2/(-p202b))*log(sqrt(mudim2)*m3/(m32-p312b))
     +         - 5d0/6d0*pi2
     +         + pi2/3d0

      boxout = boxout/(p202*(p312-m32))

      end subroutine d02g
*
* subroutine for the infrared/collinear D0 
* according to Denner-Dittmaier ArXiv:1005.2076, 
* eqs.(4.4), (4.5), (4.19)
* only m2 is complex
*
* at the end we give also 
* the expression (A.8) of Dittmaier-Kraemer PRD 65 073007
*
* we use the relation 
* D0(mu2,ml2,0,0,u,s,mu,mg,ml,mw) = D0(mu2,0,0,ml2,s,u,mg,mu,mw,ml)
* the same is true also for u -> t
*
* typical entries: p102= mq2, p212= mqp2, p322= 0.d0, p302= mlep2,
*                  p202= s, p312= u, 
*                  m02= lambda2, m12= mq2, m22= mw2, m32= mlep2
*
      subroutine d0sing(p102,p212,p322,p302,p202,p312,
     +                  m02,m12,m22in,m32,
     +                  boxout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12,m32
      complex*16 m22in
      complex*16 boxout
*
      complex*16 myli2
      external myli2

      complex*16 li2cont
      external li2cont
*
      complex*16 m22
      real*8 m3
      complex*16 m2,x32
      complex*16 li21lp1,li21lm1,argli2
      complex*16 arg1,arg2
      complex*16 arglog
*
      if (dimag(m22in).lt.epsilon) then
          m22 = m22in - ii*epsilon
      else
          m22 = m22in 
      endif
      m2  = sqrt(m22)
      m3  = sqrt(m32)
      x32 = m3/m2
*
* eq. (4.19)
*
* l= +1
*
      arg1= (m32-p312)/(m22-p212*cone)
      arg2= m2/m3*x32
      li21lp1 = li2cont(arg1,arg2)
*
* l= -1
*
      arg1= (m32-p312)/(m22-p212*cone)
      arg2= m2/m3/x32
      li21lm1 = li2cont(arg1,arg2)
*
c      arglog = m22 - p202*cone 
c      argli2= (p212-p202)/(m22-p202) * cone
      if (.not.complexmasses.and.abs(m22-mz2).le.2d0*epsilon) then
          arglog = m22 - ii*ph_ZmZw - p202*cone
          argli2= (p212-p202)/(m22-ii*ph_ZmZw-p202) * cone
      else
          arglog = m22 - p202*cone 
          argli2= (p212-p202)/(m22-p202) * cone
      endif

      boxout =  2.d0 * log((m22-p212)/arglog)
     +               * log(sqrt(mudim2)*m3/(m32-p312))
     +        + (log(sqrt(mudim2)*m3/(m32-p312)))**2 
     +        - 2.d0*myli2(argli2)
     +        + li21lp1 + li21lm1 - pi2/6.d0*cone
      boxout = boxout + pi2/12d0
      if (complexmasses.and.abs(m22-mz2).le.2d0*epsilon) then
          boxout = boxout / (p202-m22)/(p312-m32)
      else
          boxout = boxout / (p202-m22+ii*ph_ZmZw)/(p312-m32)
      endif
*
      end subroutine d0sing

*
* subroutine for the regular box with two massless internal particles
* according to Denner-Dittmaier ArXiv:1005.2076, 
* eqs.(3.78), (3.49), (3.51), (3.52), (2.4),(2.5)
*
      subroutine d0reg(p102,p212,p322,p302,p202,p312,
     +                 m02in,m12,m22,m32in,
     +                 boxout)
      implicit none
      include 'mathx.h'
*
      real*8 p102,p212,p322,p302,p202,p312
      complex*16 m02in,m32in
      complex*16 m02,m12,m22,m32

      complex*16 boxout
*
      complex*16 li2cont
      external li2cont
*
      complex*16 y01,y02,y03,y12,y13,y23
      complex*16 a,b,c,d
      complex*16 disc,rsdisc
      complex*16 li21k1,li22k1,li23k1,li24k1
      complex*16 li21k2,li22k2,li23k2,li24k2
      complex*16 x1,x2,r031,r032
      complex*16 qdr1,qdr2,cquad1,cquad2
*
      if(abs(dimag(m02)).lt.epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif

      if(abs(dimag(m32)).lt.epsilon) then
          m32 = m32in - ii*epsilon
      else
          m32 = m32in
      endif
*
      y01 = m02 + m12 - cone*p102
      y02 = m02 + m22 - cone*p202
      y03 = m02 + m32 - cone*p302
      y12 = m12 + m22 - cone*p212
      y13 = m12 + m32 - cone*p312
      y23 = m22 + m32 - cone*p322
*
      a =  y13 * y23  -  m32 * y12
      b =  y02 * y13  +  y01 * y23  -  y03 * y12
      c =  y01 * y02  -  m02 * y12
      d =  y12

      disc   = b*b - (4.d0*a*c)
      rsdisc = sqrt(disc)

      x1 = (- b + rsdisc)/2.d0/a - ii*epsilon*d/rsdisc
      x2 = (- b - rsdisc)/2.d0/a + ii*epsilon*d/rsdisc

      call qcroots(m32,y03,m02,qdr1,qdr2,cquad1,cquad2)

      r031 = -1.d0/qdr1
      r032 = -1.d0/qdr2
*
* calculates the building blocks of eq. (3.78)
*
* k=1
*
      li21k1 = li2cont(-x1,y23/y02)

      li22k1 = li2cont(-x1,r031)

      li23k1 = li2cont(-x1,r032)

      li24k1 = li2cont(-x1,y13/y01)
*
* k=2
*
      li21k2 = li2cont(-x2,y23/y02)

      li22k2 = li2cont(-x2,r031)

      li23k2 = li2cont(-x2,r032)

      li24k2 = li2cont(-x2,y13/y01)
*
      boxout=  (li21k1 - li22k1 - li23k1 + li24k1)
     +        -log(-x1)*(log(y01/y12)+log(y02/m02))
     +        -((li21k2 - li22k2 - li23k2 + li24k2)
     +        -log(-x2)*(log(y01/y12)+log(y02/m02)))
*
      if(abs(x1-x2).lt.1.d-15) then
         print*,'x1-x2 in d0reg= ',x1-x2
         print*,'x1= ',x1
         print*,'x2= ',x2
         stop
      endif

      if(abs(a).lt.1.d-15) then
         print*,'a in d0reg= ',a
         stop
      endif

      boxout  = boxout/a/(x1-x2)

      end subroutine d0reg

*
* subroutine for the regular box with generic internal particles
* according to Denner-Dittmaier arXiv:1005.2076, 
* Eq.(3.38) ('t Hooft - Veltman method)
* Eq.(3.65) (Denner-Nierste-Scharf method)
*
      subroutine myd0reg(p102,p212,p322,p302,p202,p312,
     +                   m02,m12,m22,m32,
     +                   boxout)
      implicit none
      include 'mathx.h'
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12,m22,m32
      complex*16 boxout,boxeps_2,boxeps_1,boxeps_0
      complex*16 boxouttv,boxeps_2tv,boxeps_1tv,boxeps_0tv
      complex*16 boxoutdns,boxeps_2dns,boxeps_1dns,boxeps_0dns
*
* D0 according to Eq. (3.38)
*
c      call myd0regtv(p102,p212,p322,p302,p202,p312,
c     +               m02,m12,m22,m32,
c     +               boxouttv,boxeps_2tv,boxeps_1tv,boxeps_0tv)
c
c      boxout= boxouttv
*
* D0 according to Eq. (3.65)
*
c      if(dble(m12).gt.epsilon) then
c
c          call myd0regdns(p302,p322,p212,p102,p312,p202,
c     +                    m32,m02,m12,m22,
c     +                    boxoutdns)
c
c      else

c          call myd0regdnsI(p102,p312,p212,p202,p302,p322,
c     +                    m02,m32,m12,m22,
c     +                    boxoutdns)
          call d02zeri(p102,p312,p212,p202,p302,p322,
     +                 m02,m32,m12,m22,
     +                 boxoutdns)

c      endif

      boxout = boxoutdns
      
      end
*
* subroutine for the regular box with generic internal particles
* according to Denner-Dittmaier arXiv:1005.2076, 
* Eq.(3.38) ('t Hooft - Veltman method)
*
      subroutine myd0regtv(p102in,p212in,p322in,p302in,p202in,p312in,
     +                     m02in,m12in,m22in,m32in,
     +                     boxout,boxeps_2,boxeps_1,boxeps_0)
      implicit none
*
      complex*16 p102in,p212in,p322in,p302in,p202in,p312in
      complex*16 m02in,m12in,m22in,m32in
      complex*16 boxout,boxeps_2,boxeps_1,boxeps_0
*
      complex*16 cone,zero,ii
      common/complexunit/cone,zero,ii

      real*8 epsilon
      common/epsi/epsilon
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12,m22,m32
      complex*16 al,al1,al2,quad1,quad2
      complex*16 be11,be12,be21,be22,be31,be32
      complex*16 be(3)
      complex*16 L31x11p,L21x21p,L01x31p,L01x11pp 
      complex*16 L31x21pp,L21x31pp,L30x11,L23x21,L02x31
      complex*16 L31x12p,L21x22p,L01x32p,L01x12pp 
      complex*16 L31x22pp,L21x32pp,L30x12,L23x22,L02x32
      complex*16 m2big,m02big,m12big,m22big,m32big
*
      complex*16 Y(0:3,0:3)
      complex*16 y00,y01,y02,y03,
     +           y10,y11,y12,y13,
     +           y20,y21,y22,y23,
     +           y30,y31,y32,y33
      complex*16 a,b,c
      complex*16 x1,x2
      complex*16 arg1,arg2
      complex*16 qdr1,qdr2,cquad1,cquad2
      complex*16 alin,bein,x0in,x1in,x2in,x3in
      complex*16 pab1,pab2,pab3
      complex*16 a1,b1,c1
      complex*16 a2,b2,c2
      complex*16 a3,b3,c3
      complex*16 pa1,pa2,pa3,pa4,pa5
      complex*16 pb1,pb2,pb3,pb4,pb5
      complex*16 dety(3)
      complex*16 z(3,2),x(3,2),xp(3,2),xpp(3,2)
      complex*16 r011,r021,r031,
     +           r101,r121,r131,
     +           r201,r211,r231,
     +           r301,r311,r321
      complex*16 r012,r022,r032,
     +           r102,r122,r132,
     +           r202,r212,r232,
     +           r302,r312,r322
      integer i,j
      real*8 diff,diff1,diff2,diff3
      real*8 gamma
*
      p102= p102in
      p212= p212in
      p322= p322in
      p302= p302in
      p202= p202in
      p312= p312in
      if (abs(dimag(m02in)).lt.epsilon) then
          m02= m02in - ii*epsilon
      else
          m02= m02in 
      endif
      if (abs(dimag(m12in)).lt.epsilon) then
          m12= m12in - ii*epsilon
      else
          m12= m12in
      endif
      if (abs(dimag(m22in)).lt.epsilon) then
          m22= m22in - ii*epsilon
      else
          m22= m22in
      endif
      if (abs(dimag(m32in)).lt.epsilon) then
          m32= m32in - ii*epsilon
      else
          m32= m32in
      endif

*
      boxout= zero
      boxeps_2= zero
      boxeps_1= zero
      boxeps_0= zero
*
* Eq. (2.4) Denner-Dittmaier arXiv:1005.2076
*
      Y(0,0)= m02 + m02
      Y(0,1)= m02 + m12 - cone*p102
      Y(0,2)= m02 + m22 - cone*p202
      Y(0,3)= m02 + m32 - cone*p302
      Y(1,0)= Y(0,1)
      Y(1,1)= m12 + m12
      Y(1,2)= m12 + m22 - cone*p212
      Y(1,3)= m12 + m32 - cone*p312
      Y(2,0)= Y(0,2)
      Y(2,1)= Y(1,2)
      Y(2,2)= m22 + m22
      Y(2,3)= m22 + m32 - cone*p322
      Y(3,0)= Y(0,3)
      Y(3,1)= Y(1,3)
      Y(3,2)= Y(2,3)
      Y(3,3)= m32 + m32
*
* Eq. (2.5) Denner-Dittmaier arXiv:1005.2076
*
      call qcroots(m12,y(0,1),m02,qdr1,qdr2,cquad1,cquad2)
      r011= -1.d0/qdr1
      r012= -1.d0/qdr2

      call qcroots(m22,y(0,2),m02,qdr1,qdr2,cquad1,cquad2)
      r021= -1.d0/qdr1
      r022= -1.d0/qdr2

      call qcroots(m32,y(0,3),m02,qdr1,qdr2,cquad1,cquad2)
      r031= -1.d0/qdr1
      r032= -1.d0/qdr2

      call qcroots(m02,y(1,0),m12,qdr1,qdr2,cquad1,cquad2)
      r101= -1.d0/qdr1
      r102= -1.d0/qdr2

      call qcroots(m22,y(1,2),m12,qdr1,qdr2,cquad1,cquad2)
      r121= -1.d0/qdr1
      r122= -1.d0/qdr2

      call qcroots(m32,y(1,3),m12,qdr1,qdr2,cquad1,cquad2)
      r131= -1.d0/qdr1
      r132= -1.d0/qdr2

      call qcroots(m02,y(2,0),m22,qdr1,qdr2,cquad1,cquad2)
      r201= -1.d0/qdr1
      r202= -1.d0/qdr2

      call qcroots(m12,y(2,1),m22,qdr1,qdr2,cquad1,cquad2)
      r211= -1.d0/qdr1
      r212= -1.d0/qdr2

      call qcroots(m32,y(2,3),m22,qdr1,qdr2,cquad1,cquad2)
      r231= -1.d0/qdr1
      r232= -1.d0/qdr2

      call qcroots(m02,y(3,0),m32,qdr1,qdr2,cquad1,cquad2)
      r301= -1.d0/qdr1
      r302= -1.d0/qdr2

      call qcroots(m12,y(3,1),m32,qdr1,qdr2,cquad1,cquad2)
      r311= -1.d0/qdr1
      r312= -1.d0/qdr2

      call qcroots(m22,y(3,2),m32,qdr1,qdr2,cquad1,cquad2)
      r321= -1.d0/qdr1
      r322= -1.d0/qdr2
*
* Eq. (3.7) of Denner-Dittmaier arXiv:1005.2076
*
      a= p202
      b= p302 - p322 + p202
      c= p302
      call qcroots(a,b,c,al1,al2,quad1,quad2)
      al= al1
      if(abs(al2).lt.abs(al1)) then
         al= al2
      endif
*
* Eq. (3.16) of Denner-Dittmaier arXiv:1005.2076
*
      a= p302
      b= p312 - p102 + p302
      c= p312
      call qcroots(a,b,c,be11,be12,quad1,quad2)
      be(1)= be11
      if(abs(be12).lt.abs(be11)) then
         be(1)= be12
      endif
*
      a= p322
      b= p212 - p312 + p322
      c= p212
      call qcroots(a,b,c,be21,be22,quad1,quad2)
      be(2)= be21
      if(abs(be22).lt.abs(be21)) then
         be(2)= be22
      endif
*
      a= p202
      b= p102 - p212 + p202
      c= p102
      call qcroots(a,b,c,be31,be32,quad1,quad2)
      be(3)= be31
      if(abs(be32).lt.abs(be31)) then
         be(3)= be32
      endif
*
* Eq. (3.36) of Denner-Dittmaier arXiv:1005.2076
*
      m2big= -(1.d0+al)*m02 + al*m22 + m32
      m02big= m2big - p302*cone - al*p202*cone
      m12big= m2big + (1.d0+al)*p102*cone
     +              - al*p212*cone - p312*cone
      m22big= m2big + (1.d0+al)*p202*cone - p322*cone
      m32big= m2big + (1.d0+al)*p302*cone - al*p322*cone
*
*  Eq. (3.26) of Denner-Dittmaier arXiv:1005.2076
*
      if(abs(al).lt.1.d-14) then
         print*,'alpha = ',al
         stop
      endif
      if(abs(al+1.d0).lt.1.d-14) then
         print*,'alpha = ',al
         stop
      endif
* a1
      alin= al
      bein= al*(1.d0+be(1))
      call Pab(alin,bein,p102,p202,p302,p212,p312,p322,pab1)

      alin= al
      x0in= -1.d0
      x1in= 0.d0
      x2in= 1.d0
      x3in= 0.d0
      call Pa(alin,x0in,x1in,x2in,x3in,Y,pa1)

      bein= al*(1.d0+be(1))
      x0in= -1.d0
      x1in= 0.d0
      x2in= 1.d0
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb1)

      a1= p202*pab1*cone
     +   -pa1*pb1
     +   +(1.d0+be(1))*pa1*pa1
      a1= al*a1
* b1
      alin= al
      x0in= 1.d0+al
      x1in= 0.d0
      x2in= -al
      x3in= 0.d0
      call Pa(alin,x0in,x1in,x2in,x3in,Y,pa2)

      bein= al*(1.d0+be(1))
      x0in= 1.d0+al
      x1in= 0.d0
      x2in= -al
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb2)

      b1= (Y(2,0)-2.d0*m02-2.d0*al*p202*cone)*pab1
     +   -pa2*pb1
     +   -pa1*pb2
     +   +2.d0*(1.d0+be(1))*pa1*pa2
* c1
      c1= ((1.d0+al)*Y(3,0)-al*Y(2,3)-m32)*pab1
     +   -pa2*pb2
     +   +(1.d0+be(1))*pa2*pa2
      c1= 1.d0/al*c1
* a2
      alin= al
      bein= -1.d0-(1.d0+al)*be(2)
      call Pab(alin,bein,p102,p202,p302,p212,p312,p322,pab2)

      alin= al
      x0in= 1.d0
      x1in= 0.d0
      x2in= -1.d0
      x3in= 0.d0
      call Pa(alin,x0in,x1in,x2in,x3in,Y,pa3)

      bein= -1.d0-(1.d0+al)*be(2)
      x0in= 1.d0
      x1in= 0.d0
      x2in= -1.d0
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb3)

      alin= al
      x0in= 0.d0
      x1in= 0.d0
      x2in= 1.d0
      x3in= 0.d0
      call Pa(alin,x0in,x1in,x2in,x3in,Y,pa4)

      bein= -1.d0-(1.d0+al)*be(2)
      x0in= 0.d0
      x1in= 0.d0
      x2in= 1.d0
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb4)
*
      a2= -p202*pab2*cone
     +    +pa3*pb3
     +    +be(2)*pa3*pa3
      a2= (1.d0+al)*a2
* b2
      b2= -(Y(2,0)-2.d0*m22)*pab2
     +    +pa3*pb4
     +    +pa4*pb3
     +    +2.d0*be(2)*pa3*pa4
* c2
      c2= -m22*pab2
     +    +pa4*pb4
     +    +be(2)*pa4*pa4
      c2= c2/(1.d0+al)
* a3
      alin= al
      bein= be(3)
      call Pab(alin,bein,p102,p202,p302,p212,p312,p322,pab3)

      bein= be(3)
      x0in= -1.d0
      x1in= 0.d0
      x2in= 1.d0
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb4)

      alin= al
      x0in= 1.d0
      x1in= 0.d0
      x2in= 0.d0
      x3in= 0.d0
      call Pa(alin,x0in,x1in,x2in,x3in,Y,pa5)

      bein= be(3)
      x0in= 1.d0
      x1in= 0.d0
      x2in= 0.d0
      x3in= 0.d0
      call Pb(bein,x0in,x1in,x2in,x3in,Y,pb5)
*
      a3= p202*pab3*cone
     +   -pa1*pb4
* b3
      b3= (Y(2,0)-2.d0*m02)*pab3
     +   -pa5*pb4
     +   -pa1*pb5
* c3
      c3= m02*pab3
     +   -pa5*pb5
*
*  Eq. (3.27) of Denner-Dittmaier arXiv:1005.2076
*
      
      dety(1)= b1*b1-4.d0*a1*c1
      dety(2)= b2*b2-4.d0*a2*c2
      dety(3)= b3*b3-4.d0*a3*c3

* check
      diff1= abs(dety(1)-dety(2))
      diff2= abs(dety(2)-dety(3))
      diff3= abs(dety(1)-dety(3))
      diff= max(diff1,max(diff2,diff3))
*
*  Eq. (3.28) of Denner-Dittmaier arXiv:1005.2076
*
      z(1,1)= (-b1+sqrt(dety(1)))/2.d0/a1
      z(1,2)= (-b1-sqrt(dety(1)))/2.d0/a1

      z(2,1)= (-b2+sqrt(dety(1)))/2.d0/a2
      z(2,2)= (-b2-sqrt(dety(1)))/2.d0/a2

      z(3,1)= (-b3+sqrt(dety(1)))/2.d0/a3
      z(3,2)= (-b3-sqrt(dety(1)))/2.d0/a3
*
*  Eq. (3.33) of Denner-Dittmaier arXiv:1005.2076
*
      do i=1,3
         do j=1,2
            x(i,j)= z(i,j)/(cone-z(i,j))
            xp(i,j)= -z(i,j)/(z(i,j)+be(i))
            xpp(i,j)= (cone-z(i,j))/(z(i,j)+be(i))
         enddo
      enddo
*
* Eq. (3.38)
*
* k=1
*
      call Lij(m32big,m32,m12big,r311,r312,xp(1,1),L31x11p) 
      call Lij(m22big,m22,m12big,r211,r212,xp(2,1),L21x21p) 
      call Lij(m02big,m02,m12big,r011,r012,xp(3,1),L01x31p) 
      call Lij(m02big,m02,m12big,r011,r012,xpp(1,1),L01x11pp) 
      call Lij(m32big,m32,m12big,r311,r312,xpp(2,1),L31x21pp) 
      call Lij(m22big,m22,m12big,r211,r212,xpp(3,1),L21x31pp) 
      call Lij(m32big,m32,m02big,r301,r302,x(1,1),L30x11) 
      call Lij(m22big,m22,m32big,r231,r232,x(2,1),L23x21) 
      call Lij(m02big,m02,m22big,r021,r022,x(3,1),L02x31) 
*
* k=2
*
      call Lij(m32big,m32,m12big,r311,r312,xp(1,2),L31x12p) 
      call Lij(m22big,m22,m12big,r211,r212,xp(2,2),L21x22p) 
      call Lij(m02big,m02,m12big,r011,r012,xp(3,2),L01x32p)
      call Lij(m02big,m02,m12big,r011,r012,xpp(1,2),L01x12pp)
      call Lij(m32big,m32,m12big,r311,r312,xpp(2,2),L31x22pp) 
      call Lij(m22big,m22,m12big,r211,r212,xpp(3,2),L21x32pp)
      call Lij(m32big,m32,m02big,r301,r302,x(1,2),L30x12) 
      call Lij(m22big,m22,m32big,r231,r232,x(2,2),L23x22) 
      call Lij(m02big,m02,m22big,r021,r022,x(3,2),L02x32) 
*
      boxout= 1.d0/sqrt(dety(1)) 
     +            * ( + (  L31x11p + L21x21p + L01x31p 
     +                   - L01x11pp - L31x21pp - L21x31pp
     +                   - L30x11 - L23x21 - L02x31)
     +                - (  L31x12p + L21x22p + L01x32p
     +                   - L01x12pp - L31x22pp - L21x32pp
     +                   - L30x12 - L23x22 - L02x32) )
      boxeps_0= boxout
*
      end
*
* subroutine for the regular box with generic internal particles
* according to Denner-Dittmaier arXiv:1005.2076, 
* Eq.(3.65) (Denner-Nierste-Scharf method)
*
      subroutine myd0regdns(p102in,p212in,p322in,p302in,p202in,p312in,
     +                      m02in,m12in,m22in,m32in,
     +                      boxout)
      implicit none
      include 'mathx.h'
*
      complex*16 p102in,p212in,p322in,p302in,p202in,p312in
      complex*16 m02in,m12in,m22in,m32in
      complex*16 boxout
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12,m22,m32
*
      complex*16 Y(0:3,0:3)
      complex*16 y00,y01,y02,y03,
     +           y10,y11,y12,y13,
     +           y20,y21,y22,y23,
     +           y30,y31,y32,y33
      complex*16 a,b,c,d
      complex*16 arg1,arg2
      complex*16 qdr1,qdr2,cquad1,cquad2
      complex*16 dety,x1,x2
      complex*16 r011,r021,r031,
     +           r101,r121,r131,
     +           r201,r211,r231,
     +           r301,r311,r321
      complex*16 r012,r022,r032,
     +           r102,r122,r132,
     +           r202,r212,r232,
     +           r302,r312,r322
      complex*16 li2cont111,li2cont112,li2cont113,li2cont114
      complex*16 li2cont121,li2cont122,li2cont123,li2cont124
      complex*16 li2cont211,li2cont212,li2cont213,li2cont214
      complex*16 li2cont221,li2cont222,li2cont223,li2cont224
      complex*16 p1numout,p1denout,p2numout,p2denout
      complex*16 eta1r13,eta2r13
      complex*16 d01out
      complex*16 d02out
      integer i,j
      real*8 gamma
*
* choice of gamma (Eqs. (3.68) and (3.69))
*
* no t-channel-like invariant as input
*
c      if(p102.ge.0.d0.and.
c     +   p212.ge.0.d0.and.
c     +   p322.ge.0.d0.and.
c     +   p302.ge.0.d0.and.
c     +   p202.ge.0.d0.and.
c     +   p312.ge.0.d0) then
c      endif
*
* still to be implemented: tha case of two t-channel-like invariants
*
c         gamma= 1.d0/2.d0/p312*(p312+p302-p102
c     +            -sqrt((p312+p302-p102)**2-4.d0*p312*p302))
c      print*,'gamma= ',gamma
*
      gamma= (   p312in + p302in - p102in
     +         + sqrt((p312in+p302in-p102in)**2-4.d0*p312in*p302in) )
     +        /2.d0/p312in
      if(gamma.lt.0.d0.or.gamma.gt.1.d0) then
         print*,'gamma= ',gamma
         print*,'change permutation of input!!!'
         stop
      endif
*
* redefinition of inputs to compute the first D0 of Eq.(3.65)
*
c      p102= p302
c      p212= gamma*gamma*p312in
c      p322= gamma*(p212in - p102in) 
c     +    + (1.d0-gamma)*(p322in - p302in)
c      p302= p202in
c      p202= 0d0
c      p312= p322in
c
c      m02= m02in
c      m12= m32in
c      m22= gamma*(m12in-p102*cone)+(1.d0-gamma)*(m32in-p302in*cone)
c      m32= m22in

      p102= 0d0
      p212= gamma*(p212in - p102in) 
     +    + (1.d0-gamma)*(p322in - p302in)
      p322= p322in
      p302= p302in
      p202= p202in 
      p312= gamma**2*p312in
      m02= m02in
      m12= gamma*(m12in-p102*cone)+(1.d0-gamma)*(m32in-p302in*cone)
      m22= m22in
      m32= m32in
*
      call myd0regdnsI(p102,p212,p322,p302,p202,p312,
     +                 m02,m12,m22,m32,
     +                 d01out)
*
* redefinition of inputs to compute the second D0 of Eq.(3.65)
*
c      p102= p102in
c      p212= (1.d0-gamma)**2*p312in
c      p322= gamma*(p212in - p102in) 
c     +    + (1.d0-gamma)*(p322in - p302in)
c      p302= p202in
c      p202= 0.d0
c      p312= p212in
c      m02= m02in
c      m12= m12in
c      m22= gamma*(m12in-p102in*cone)+(1.d0-gamma)*(m32in-p302in*cone)
c      m32= m22in

      p102= p102in
      p212= p212in
      p322= gamma*(p212in - p102in) 
     +    + (1.d0-gamma)*(p322in - p302in)
      p302= 0d0
      p202= p202in
      p312= (1.d0-gamma)**2*p312in
      m02= m02in
      m12= m12in
      m22= m22in
      m32= gamma*(m12in-p102in*cone)+(1.d0-gamma)*(m32in-p302in*cone)

      call myd0regdnsI(p102,p212,p322,p302,p202,p312,
     +                 m02,m12,m22,m32,
     +                 d02out)
*
      boxout = gamma*d01out + (1.d0-gamma)*d02out
*
      end
*
**
*
      subroutine myd0regdnsI(p102in,p212in,p322in,p302in,p202in,p312in,
     +                       m02in,m12in,m22in,m32in,
     +                       d0iout)
      implicit none
      include 'mathx.h'
*
      complex*16 p102in,p212in,p322in,p302in,p202in,p312in
      complex*16 m02in,m12in,m22in,m32in
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02,m12,m22,m32
      complex*16 d0iout
*
      complex*16 li2cont
      external li2cont

      complex*16 myeta
      external myeta
*
      complex*16 Y(0:3,0:3)
      complex*16 y00,y01,y02,y03,
     +           y10,y11,y12,y13,
     +           y20,y21,y22,y23,
     +           y30,y31,y32,y33
      complex*16 a,b,c,d
      complex*16 arg1,arg2
      complex*16 qdr1,qdr2,cquad1,cquad2
      complex*16 dety,x1,x2
      complex*16 r011,r021,r031,
     +           r101,r121,r131,
     +           r201,r211,r231,
     +           r301,r311,r321
      complex*16 r012,r022,r032,
     +           r102,r122,r132,
     +           r202,r212,r232,
     +           r302,r312,r322
      complex*16 li2cont111,li2cont112,li2cont113,li2cont114
      complex*16 li2cont121,li2cont122,li2cont123,li2cont124
      complex*16 li2cont211,li2cont212,li2cont213,li2cont214
      complex*16 li2cont221,li2cont222,li2cont223,li2cont224
      complex*16 p1numout,p1denout,p2numout,p2denout
      complex*16 eta1r13,eta2r13
      integer i,j
      real*8 gamma
*
      p102= p102in + ii*epsilon
      p212= p212in + ii*epsilon
      p322= p322in + ii*epsilon
      p302= p302in + ii*epsilon
      p202= p202in + ii*epsilon
      p312= p312in + ii*epsilon

      if (abs(dimag(m02in)).lt.epsilon) then
          m02= m02in - ii*epsilon
      else
          m02= m02in
      endif
      if (abs(dimag(m12in)).lt.epsilon) then
          m12= m12in - ii*epsilon
      else
          m12= m12in
      endif
      if (abs(dimag(m22in)).lt.epsilon) then
          m22= m22in - ii*epsilon
      else
          m22= m22in
      endif
      if (abs(dimag(m32in)).lt.epsilon) then
          m32= m32in - ii*epsilon
      else
          m32= m32in 
      endif
c      print*,m12
*
* Eq. (2.4) Denner-Dittmaier arXiv:1005.2076
*
      Y(0,0)= m02 + m02
      Y(0,1)= m02 + m12 - cone*p102
      Y(0,2)= m02 + m22 - cone*p202
      Y(0,3)= m02 + m32 - cone*p302
      Y(1,0)= Y(0,1)
      Y(1,1)= m12 + m12
      Y(1,2)= m12 + m22 - cone*p212
      Y(1,3)= m12 + m32 - cone*p312
      Y(2,0)= Y(0,2)
      Y(2,1)= Y(1,2)
      Y(2,2)= m22 + m22
      Y(2,3)= m22 + m32 - cone*p322
      Y(3,0)= Y(0,3)
      Y(3,1)= Y(1,3)
      Y(3,2)= Y(2,3)
      Y(3,3)= m32 + m32
*
* Eq. (2.5) Denner-Dittmaier arXiv:1005.2076
*
      call qcroots(m12,Y(0,1),m02,qdr1,qdr2,cquad1,cquad2)
      r011= -1.d0/qdr1
      r012= -1.d0/qdr2

      call qcroots(m22,Y(0,2),m02,qdr1,qdr2,cquad1,cquad2)
      r021= -1.d0/qdr1
      r022= -1.d0/qdr2

      call qcroots(m32,Y(0,3),m02,qdr1,qdr2,cquad1,cquad2)
      r031= -1.d0/qdr1
      r032= -1.d0/qdr2

      call qcroots(m02,Y(1,0),m12,qdr1,qdr2,cquad1,cquad2)
      r101= -1.d0/qdr1
      r102= -1.d0/qdr2

      call qcroots(m22,Y(1,2),m12,qdr1,qdr2,cquad1,cquad2)
      r121= -1.d0/qdr1
      r122= -1.d0/qdr2

      call qcroots(m32,Y(1,3),m12,qdr1,qdr2,cquad1,cquad2)
      r131= -1.d0/qdr1
      r132= -1.d0/qdr2

      call qcroots(m02,Y(2,0),m22,qdr1,qdr2,cquad1,cquad2)
      r201= -1.d0/qdr1
      r202= -1.d0/qdr2

      call qcroots(m12,Y(2,1),m22,qdr1,qdr2,cquad1,cquad2)
      r211= -1.d0/qdr1
      r212= -1.d0/qdr2

      call qcroots(m32,Y(2,3),m22,qdr1,qdr2,cquad1,cquad2)
      r231= -1.d0/qdr1
      r232= -1.d0/qdr2

      call qcroots(m02,Y(3,0),m32,qdr1,qdr2,cquad1,cquad2)
      r301= -1.d0/qdr1
      r302= -1.d0/qdr2

      call qcroots(m12,Y(3,1),m32,qdr1,qdr2,cquad1,cquad2)
      r311= -1.d0/qdr1
      r312= -1.d0/qdr2

      call qcroots(m22,Y(3,2),m32,qdr1,qdr2,cquad1,cquad2)
      r321= -1.d0/qdr1
      r322= -1.d0/qdr2
*
* Eq. (3.50) of Denner-Dittmaier arXiv:1005.2076
*
      a = m12*r131*(Y(2,3)-Y(0,3)/r201) - m32*(Y(1,2) - Y(0,1)/r201)
      b =  (m22*r201 - m02/r201)*(m12*r131 - m32/r131)
     +  + Y(0,1)*Y(2,3) - Y(0,3)*Y(1,2)
      c = m22*r201*(Y(0,1)-Y(0,3)/r131) - m02*(Y(1,2)-Y(2,3)/r131)
      d = Y(1,2) - Y(0,1)/r201 - Y(2,3)/r131 + Y(0,3)/(r201*r131)

*
* Eq. (3.52)
*
      dety = b*b - 4.d0*a*c
      x1   = ( -b + sqrt(dety) )/2.d0/a - ii*epsilon*d/sqrt(dety)
      x2   = ( -b - sqrt(dety) )/2.d0/a + ii*epsilon*d/sqrt(dety)
*
* building blocks of Eq. (3.61)
*
c k=1
c l=1
      li2cont111= li2cont(-x1,r231/r201)
      li2cont112= li2cont(-x1,r031)
      li2cont113= li2cont(-x1*r131,r211/r201)
      li2cont114= li2cont(-x1*r131,r011)
*
c k=1
c l=2
      li2cont121= li2cont(-x1,r232/r201)
      li2cont122= li2cont(-x1,r032)
      li2cont123= li2cont(-x1*r131,r212/r201)
      li2cont124= li2cont(-x1*r131,r012)
*
c k=2
c l=1
      li2cont211= li2cont(-x2,r231/r201)
      li2cont212= li2cont(-x2,r031)
      li2cont213= li2cont(-x2*r131,r211/r201)
      li2cont214= li2cont(-x2*r131,r011)
*
c k=2
c l=2
      li2cont221= li2cont(-x2,r232/r201)
      li2cont222= li2cont(-x2,r032)
      li2cont223= li2cont(-x2*r131,r212/r201)
      li2cont224= li2cont(-x2*r131,r012)
*
      arg1= x1/r201
      call P(Y,m02,m12,m22,m32,zero,zero,cone,arg1,p1numout)
      arg1= x1
      call P(Y,m02,m12,m22,m32,cone,zero,zero,arg1,p1denout)

      arg1= x2/r201
      call P(Y,m02,m12,m22,m32,zero,zero,cone,arg1,p2numout)
      arg1= x2
      call P(Y,m02,m12,m22,m32,cone,zero,zero,arg1,p2denout)
*
      d0iout= 
     +         + li2cont111 - li2cont112 - li2cont113 + li2cont114
     +         + li2cont121 - li2cont122 - li2cont123 + li2cont124
     +      - (  li2cont211 - li2cont212 - li2cont213 + li2cont214
     +         + li2cont221 - li2cont222 - li2cont223 + li2cont224 )
     +      -myeta(-x1,r131)*(log( (p1numout-ii*epsilon/r201**2)
     +                    /(p1denout-ii*epsilon)) 
     +                -log(m22/m02))
     +      +myeta(-x2,r131)*(log( (p2numout-ii*epsilon/r201**2)
     +                    /(p2denout-ii*epsilon)) 
     +                -log(m22/m02))
      d0iout= d0iout/a/(x1-x2)
*
      end
*
** Eq. (3.78) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine d02zeri(p102,p212,p322,p302,p202,p312,
     +                   m02in,m12,m22,m32in,
     +                   boxout)
      implicit none
      include 'mathx.h'
*
      complex*16 p102,p212,p322,p302,p202,p312
      complex*16 m02in,m32in
      complex*16 m02,m12,m22,m32

      complex*16 boxout
*
      complex*16 li2cont
      external li2cont
*
      complex*16 y01,y02,y03,y12,y13,y23
      complex*16 a,b,c,d
      complex*16 disc,rsdisc
      complex*16 x1,x2,r031,r032
      complex*16 qdr1,qdr2,cquad1,cquad2
*
      if(abs(dimag(m02)).le.2d0*epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif

      if(abs(dimag(m32)).le.2d0*epsilon) then
          m32 = m32in - ii*epsilon
      else
          m32 = m32in
      endif
*
      y01 = + m02 - cone*p102
      y02 = + m02 - cone*p202
      y03 = + m02 + m32 - cone*p302
      y12 = - cone*p212
      y13 = + m32 - cone*p312
      y23 = + m32 - cone*p322
*
      a =  y13 * y23  -  m32 * y12
      b =  y02 * y13  +  y01 * y23  -  y03 * y12
      c =  y01 * y02  -  m02 * y12
      d =  y12
*
      disc   = b*b - (4d0*a*c)
      rsdisc = sqrt(disc)

      x1 = (- b + rsdisc)/2.d0/a - ii*epsilon*d/rsdisc
      x2 = (- b - rsdisc)/2.d0/a + ii*epsilon*d/rsdisc

      call qcroots(m32,y03,m02,qdr1,qdr2,cquad1,cquad2)

      r031 = -1.d0/qdr1
      r032 = -1.d0/qdr2
*
      boxout = ( + (  li2cont(-x1,y23/y02) - li2cont(-x1,r031)
     +              - li2cont(-x1,r032)    + li2cont(-x1,y13/y01)
     +              - log(-x1) * ( log(y01/y12)+log(y02/m02) )
     +             )
     +           - (  li2cont(-x2,y23/y02) - li2cont(-x2,r031)
     +              - li2cont(-x2,r032)    + li2cont(-x2,y13/y01)
     +              - log(-x2) * ( log(y01/y12)+log(y02/m02) )
     +             ) ) / (a*(x1-x2))

      end
*
* Eq. (3.39) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine Lij(mi2big,mi2,mj2big,rij1,rij2,x,Lijout)
      integer i,j
      complex*16 mi2big,mi2,mj2big,rij1,rij2,x
      complex*16 Lijout
*
      complex*16 li2cont
      external li2cont
*
      complex*16 arg1,arg2,cone
      complex*16 li2contmjmi,li2contrij1,li2contrij2
      complex*16 myli2,arg,li2opx
      external myli2
      real*8 zero
*
      zero= 1.d-15
*
      cone= (1.d0,0.d0)
      arg= cone+x
      li2opx= myli2(arg)
*
      if(abs(mi2big).lt.zero) then

         arg1= rij1
         arg2= -x
         li2contrij1 = li2cont(arg1,arg2)
*
         arg1= rij2
         arg2= -x
         li2contrij2 = li2cont(arg1,arg2)
*
         Lijout= -(log(mj2big)-log(mi2))*log(-x)
     +           -0.5d0*(log(-x))**2
     +           +li2opx
     +           -li2contrij1
     +           -li2contrij2
*
      else
*
         arg1= mj2big/mi2big
         arg2= -x
         li2contmjmi = li2cont(arg1,arg2)
*
         arg1= rij1
         arg2= -x
         li2contrij1 = li2cont(arg1,arg2)
*
         arg1= rij2
         arg2= -x
         li2contrij2 = li2cont(arg1,arg2)
*
         Lijout= -(log(mi2big)-log(mi2))*log(-x)
     +           +li2opx
     +           +li2contmjmi
     +           -li2contrij1
     +           -li2contrij2
*
      endif

      end
*
* Eq. (2.3) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine P(Y,m02,m12,m22,m32,x0,x1,x2,x3,pout)
      implicit none
*
      complex*16 Y(0:3,0:3)
      complex*16 m02,m12,m22,m32
      complex*16 x0,x1,x2,x3
      complex*16 pout
*
      pout= m02*x0**2+m12*x1**2+m22*x2**2+m32*x3**2
     +    + y(0,1)*x0*x1 + y(0,2)*x0*x2 + y(0,3)*x0*x3
     +    + y(1,2)*x1*x2 + y(1,3)*x1*x3
     +    + y(2,3)*x2*x3
*
      end
*
* Eq. (3.41) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine Q(Y,m02,m12,m22,m32,x0,x1,x2,x3,r20,qout)
      implicit none
*
      complex*16 Y(0:3,0:3)
      complex*16 m02,m12,m22,m32
      complex*16 x0,x1,x2,x3
      complex*16 r20
      complex*16 qout
*
      qout= -( 2d0*m02*x0 + y(0,1)*x1 + y(0,2)*x2 + y(0,3)*x3)/r20 +
     +         2d0*m22*x2 + y(0,2)*x0 + y(1,2)*x1 + y(2,3)*x2
*
      end

*
* Eq. (3.2) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine Pk(k,Y,x0,x1,x2,x3,pkout)
      implicit none
*
      integer k
      complex*16 Y(0:3,0:3)
      complex*16 x0,x1,x2,x3
      complex*16 pkout
*
      if(k.eq.0) then
         pkout= y(0,0)*x0 + y(1,0)*x1 + y(2,0)*x2 + y(3,0)*x3
      elseif(k.eq.1) then
         pkout= y(0,1)*x0 + y(1,1)*x1 + y(2,1)*x2 + y(3,1)*x3
      elseif(k.eq.2) then
         pkout= y(0,2)*x0 + y(1,2)*x1 + y(2,2)*x2 + y(3,2)*x3
      elseif(k.eq.3) then
         pkout= y(0,3)*x0 + y(1,3)*x1 + y(2,3)*x2 + y(3,3)*x3
      endif
*
      end
*
      subroutine Pa(al,x0,x1,x2,x3,Y,paout)
      implicit none
      complex*16 al,x0,x1,x2,x3
      complex*16 Y(0:3,0:3)
      complex*16 paout
*
      complex*16 p0,p2,p3
*
      call Pk(0,Y,x0,x1,x2,x3,p0)
      call Pk(2,Y,x0,x1,x2,x3,p2)
      call Pk(3,Y,x0,x1,x2,x3,p3)
      paout= -(1.d0+al)*p0 + al*p2 + p3
*
      end
*
      subroutine Pb(be,x0,x1,x2,x3,Y,pbout)
      implicit none
      complex*16 be,x0,x1,x2,x3
      complex*16 Y(0:3,0:3)
      complex*16 pbout
*
      complex*16 p0,p1,p2
*
      call Pk(0,Y,x0,x1,x2,x3,p0)
      call Pk(1,Y,x0,x1,x2,x3,p1)
      call Pk(2,Y,x0,x1,x2,x3,p2)
      pbout= -(1.d0+be)*p0 + p1 + be*p2
*
      end
*
* Eq. (3.3) of Denner-Dittmaier arXiv:1005.2076
*
      subroutine Pab(al,be,p102,p202,p302,p212,p312,p322,pabout)
      implicit none
*
      complex*16 al,be,p102,p202,p302,p212,p312,p322
      complex*16 pabout
*
      pabout= (1.d0+al)*p102
     +      +(al+be+2.d0*al*be)*p202
     +      +(1.d0+be)*p302
     +      -al*p212 - p312 - be*p322
*
      end

*
** end of scalar point functions
*
**********************************************************
**********************************************************
*
** other useful functions
*
** eq. (4.28) of ArXiv:0709.1075 (Denner Fortschritte)
*
      complex*16 function kallen(x,y,z)
      implicit none
      complex*16 x,y,z

      kallen= sqrt(x**2+y**2+z**2-2.d0*(x*y+x*z+y*z))

      end function kallen
*
***************************************************************
**  DILOGARITHM STUFF
*
** if |z| > 1 calculate li2(1/z) and convert
*
      complex*16  function myli2(zz)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 zz
*
      complex*16 myli22
      external myli22
*
      if(abs(zz-cone).lt.1.d-8) then
         myli2=cone*pi2/6.d0
         return
      endif

      if (abs(zz).le.1.d0) then
         myli2=myli22(zz)
      else  
         myli2=-log(-zz)**2/2d0-cone*pi2/6d0-myli22(cone/zz)
      endif

      end function myli2
*
** if |z|>1/2 calculate li2(1-z) and convert
*
      complex*16 function myli22(zz)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 zz
*
      complex*16 myli23
      external myli23
*
      if (abs(zz).lt.0.5d0) then
         myli22=myli23(zz)
      else
         myli22=cone*pi2/6d0-log(zz)*log(cone-zz)-myli23(cone-zz)
      endif
*
      end function myli22
*
**  main function for dilogarithms
*
      complex*16 function myli23(zz)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 zz
*
      integer dim
      parameter (dim = 25)
      real*8 b(dim)
      integer i,n
      complex*16 z,add
*
* even n Bernoulli numbers, already divided by (n + 1)!
*
      data b /
     $  2.777777777777777777777777777777777777777777777777777777778d-2,
     $ -2.777777777777777777777777777777777777777777777777777777778d-4,
     $  4.724111866969009826152683295540438397581254724111866969009d-6,
     $ -9.185773074661963550852439741328630217519106407995296884185d-8,
     $  1.897886998897099907200917301927402937503947604957705967806d-9,
     $ -4.064761645144225526805909386291966674547057127439707822288d-11,
     $  8.921691020456452555217987316752748851514283613049045147810d-13,
     $ -1.993929586072107568723644347793789705630694749653880147036d-14,
     $  4.518980029619918191650476552855593228396819014466618405199d-16,
     $ -1.035651761218124701448341154221865666596091238168650515964d-17,
     $  2.395218621026186745740283743000980381678949001942974256251d-19,
     $ -5.581785874325009336283074505625419905567054667644398095136d-21,
     $  1.309150755418321285812307399186592301749849838783303836854d-22,
     $ -3.087419802426740293242279764866462431595565256132745695326d-24,
     $  7.315975652702203420357905609252148591033401063690875035693d-26,
     $ -1.740845657234000740989055147759702545340841421754271264171d-27,
     $  4.157635644613899719617899620775226673488254159511563860825d-29,
     $ -9.962148488284622103194006702455838849854860017394488768062d-31,
     $  2.394034424896165300521167987893749562934279156932915750221d-32,
     $ -5.768347355367390084291793161877654244072332317926275110062d-34,
     $  1.393179479647007977827886603911548331732411625673399565806d-35,
     $ -3.372121965485089470468473635254930958979742891656539304386d-37,
     $  8.178208777562102621764777214872834267876189462495503276198d-39,
     $ -1.987010831152385925564820669234786567541858995824743201790d-40,
     $  4.835778518040550896287059373115378207694465369420827842732d-42
     $  /
*
      z=-log(cone-zz)

      myli23=z-z**2/4d0

      if (abs(myli23).gt.epsilon) then
          do i=1,dim
             n    = 2*i
             add  = z**(n+1)*b(i)
             if (abs(add/myli23).gt.epsilon) then
                 myli23 = myli23 + add
             else
                 exit
             endif
          enddo
      endif

      end function myli23
*
* analytically continued dilogarithm of two variables
* according to eq. (2.14) of Denner-Dittmaier ArXiv:1005.2076
*
      complex*16 function li2cont(arg1,arg2)
      implicit none
      include 'mathx.h'
*
      complex*16 arg1,arg2
*
      complex*16 myli2
      external myli2
      complex*16 myeta
      external myeta
*
      complex*16 arg
*
      arg     = cone - arg1*arg2
      li2cont = myli2(arg) + myeta(arg1,arg2)*log(arg)
*
      end function li2cont
*
* analytically continued dilogarithm of three variables
* according to eq. (2.13) of Denner-Dittmaier ArXiv:1005.2076
*
      complex*16 function li23arg(arg1,arg2,arg3)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 arg1,arg2,arg3
*
      complex*16 myli2
      external myli2
*
      complex*16 arg,arg123
      real*8 argif,theta
*
      arg123 = arg1*arg2*arg3
      arg    = cone - arg123

      argif= abs(arg123)-1.d0
      theta= 0.d0
      if(argif.gt.0.d0) theta= 1.d0

      li23arg = myli2(arg)
     +          + (log(arg)- theta*( log(-arg123)
     +                     -0.5d0*log(arg123)
     +                     -0.5d0*(log(arg1)+log(arg2)+log(arg3))))
     +           *(log(arg123)-log(arg1)-log(arg2)-log(arg3))
*
      end function li23arg
*
      complex*16 function myeta(arg1,arg2)
*
* computes eta(arg1,arg2)  (  log(a * b) = log(a) + log(b) + eta(a,b)  )
*
* remembering that eta(a,b)= 2 pi i (theta(-im a) theta(-im b) theta(im(a b))
*                                   -theta(im a) theta(im b) theta(-im(a b)))
*
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 arg1,arg2
*
      complex*16 arg12,eta1,eta2

      arg12= arg1*arg2
      eta1= zero
      if((dimag(-arg1)).ge.0.d0) then
         if((dimag(-arg2)).ge.0.d0) then
            if((dimag(arg12)).ge.0.d0) then
               eta1= cone
            endif
         endif
      endif
      eta2= zero
      if((dimag(arg1)).ge.0.d0) then
         if((dimag(arg2)).ge.0.d0) then
            if((dimag(-arg12)).ge.0.d0) then
               eta2= cone
            endif
         endif
      endif
      
      myeta= 2.d0*pi * ii * (eta1 - eta2)
*
      end function myeta
*
* subroutine for precise calculation of the roots of quadratic equations
*
      subroutine qcroots(qda,qdb,qdc,qdr1,qdr2,cquad1,cquad2)
      implicit none
      complex*16 qda,qdb,qdc,qdr1,qdr2,cquad1,cquad2
      complex*16 qa,qb,qc,qr1,qr2
      real*8 si
      complex*16 sqdisc,argsq
      complex*16 qtmp
*
      qa= qda
      qb= qdb
      qc= qdc
      if (abs(qda).gt.0.d0) then
         if (abs(qb).gt.0.d0) then
            argsq=1.d0 *
     +           (1.d0-2.d0*sqrt(qa*qc)/qb)*(1.d0+2.d0*sqrt(qa*qc)/qb)
            argsq = dcmplx(dble(argsq),dimag(- 4.d0*qa*qc/qb/qb))
            sqdisc = qb*sqrt(argsq)
            si = 1.d0
            if (dble(conjg(qb)*sqdisc).lt.0.d0) si = -1.d0
         else
            argsq=1.d0 *
     +           (qb-2.d0*sqrt(qa*qc))*(qb+2.d0*sqrt(qa*qc))
            argsq = dcmplx(dble(argsq),dimag(- 4.d0*qa*qc))
            sqdisc = sqrt(argsq)
            si = 1.d0
            if (dble(conjg(qb)*sqdisc).lt.0.d0) si = -1.d0
         endif

         qtmp = -0.5d0*(qb + si*sqdisc)
         qr1 = qtmp/qa
         qr2 = qc/qr1/qa
      else
         if (abs(qdb).gt.0.d0) then
            qr1 = -qdc/qdb
            qr2 = qr1
         else
            qr1 = (0.d0,0.d0)
            qr2 = (0.d0,0.d0)
         endif
      endif
      cquad1= qa*qr1*qr1 + qb*qr1 + qc
      cquad2= qa*qr2*qr2 + qb*qr2 + qc
      qdr1 = qr1
      qdr2 = qr2
      end subroutine qcroots


*
* subroutine for precise calculation of the roots of quadratic equations
* real arguments 
*
      subroutine qroots(qda,qdb,qdc,qdr1,qdr2,cquad1,cquad2)
      implicit none
      real*8 qda,qdb,qdc,qdr1,qdr2,cquad1,cquad2
      complex*16 qa,qb,qc,qr1,qr2,qr3,qcquad1,qcquad2
      real*8 aqb,si
      double precision ar,ai,br,bi,cr,ci,x1r,x1i,x2r,x2i,tmp
      complex*16 qdisc,sqdisc,argsq
      complex*16 qr1c,qr2c,qtmp,sol1,sol2,cq1,cq2
      integer i
*
      qa= qda*(1.d0,0.d0)
      qb= qdb*(1.d0,0.d0)
      qc= qdc*(1.d0,0.d0)
      if (abs(qda).gt.0.d0) then
         if (abs(qb).gt.0.d0) then
            argsq=1.d0 *
     $           (1.d0-2.d0*sqrt(qa*qc)/qb)*(1.d0+2.d0*sqrt(qa*qc)/qb)
!
            argsq = dcmplx(dble(argsq),dimag(- 4.d0*qa*qc/qb/qb))
            sqdisc = qb*sqrt(argsq)
            si = 1.d0
            if (dble(conjg(qb)*sqdisc).lt.0.d0) si = -1.d0
         else
            argsq=1.d0 *
     $           (qb-2.d0*sqrt(qa*qc))*(qb+2.d0*sqrt(qa*qc))
!
            argsq = dcmplx(dble(argsq),dimag(- 4.d0*qa*qc))
            sqdisc = sqrt(argsq)
            si = 1.d0
            if (dble(conjg(qb)*sqdisc).lt.0.d0) si = -1.d0
         endif

         qtmp = -0.5d0*(qb + si*sqdisc)
         qr1 = qtmp/qa
         qr2 = qc/qr1/qa
      else
         if (abs(qdb).gt.0.d0) then
            qr1 = -qdc/qdb
            qr2 = qr1
         else
            qr1 = (0.d0,0.d0)
            qr2 = (0.d0,0.d0)
         endif
      endif
      argsq= qa*qr1*qr1 + qb*qr1 + qc
      cquad1= dble(argsq)
      argsq= qa*qr2*qr2 + qb*qr2 + qc
      cquad2= dble(argsq)
      if(abs(dimag(qr1)).gt.1.d-15.or.abs(dimag(qr2)).gt.1.d-15) then
         print*,'alpha1 is not real: ',qr1
         print*,'alpha2 is not real: ',qr2
      endif
      qdr1 = dble(qr1)
      qdr2 = dble(qr2)
      end
*
** deltar according to eq. (8.14) of ArXiv:0709.1075 (Denner Fortschritte)
*
      complex*16 function getdeltar()
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
*
      complex*16 swtmw2
      complex*16 swtpmw2
      complex*16 swt0
      complex*16 szztmz2
      complex*16 szztpmz2
      complex*16 sazt0
      complex*16 saatp0
      complex*16 deltamw2
      complex*16 deltamz2
*
      call sigmawt(dble(mw2)*cone,swtmw2)
      call sigmawtp(dble(mw2)*cone,swtpmw2)
      call sigmawt(zero,swt0)
      call sigmazzt(dble(mz2)*cone,szztmz2)
      call sigmazztp(dble(mz2)*cone,szztpmz2)
      call sigmaazt(zero,sazt0)
      call sigmaaatp(zero,saatp0,.true.)

      deltamw2  = swtmw2 + ii*dimag(mw2)*swtpmw2
     +          - ii*dimag(mw2)*alpha/pi

      deltamz2  = szztmz2 + ii*dimag(mz2)*szztpmz2
*
      getdeltar = saatp0 
     +          - cw2/sw2*
     +            ( deltamz2/mz2 - deltamw2/mw2 )
     +          + ( swt0 - deltamw2 )/mw2 
     +          + 2.d0*cw/sw*sazt0/mz2
     +          + alsu4pi/sw2 
     +              * ( 6.d0*cone+(7.d0*cone-4.d0*sw2)
     +                              /2.d0/sw2*log(cw2) )
*
      end function getdeltar
