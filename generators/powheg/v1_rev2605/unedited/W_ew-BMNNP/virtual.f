c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual
      real * 8 virtual_st
      complex*16 virtual_ew
      real * 8 dummy(0:3,0:3)
      complex*16 born
      real *8 s,dotp
      external dotp

      s=2d0*dotp(p(0,1),p(0,2))
      call compborn(p,vflav,born,dummy)
      virtual_st = pi2 - 8 - 3*log(st_muren2/s) - log(st_muren2/s)**2

      call deltavirt(p,vflav,virtual_ew)

      virtual=( virtual_st*cf + 
     +          2d0*virtual_ew/(st_alpha/(2d0*pi)) )
     +         *born

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
** (pu) u \            / l+ (pl)
**         \----------/
**         |\________/    (etc)
**         |/   W+   \ 
**         /          \
** (pd) d~/            \ nu_l (pnu)
*

      subroutine deltavirt(p,flav,delvirt)
      implicit none
      include 'nlegborn.h'
      include 'mathx.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real* 8 p(0:3,nlegborn)
      integer flav(nlegborn)
      complex*16 delvirt
*
      real*8 dotp
      external dotp

      real * 8 chargeofparticle
      external chargeofparticle
*
      real*8 s,t,u
      complex*16 delww  
      complex*16 delwdu 
      complex*16 delwnul
      complex*16 delbox 

      real*8 pu(0:3),pd(0:3),pn(0:3),pl(0:3)
      real*8 ptmpu(0:3),ptmpt(0:3),ptmps(0:3)
      integer nu
*
      if (mod(abs(flav(1)),2).eq.0) then
          do nu=0,3
              pu(nu) = p(nu,1)
              pd(nu) = p(nu,2)
          enddo
      else
          do nu=0,3
              pd(nu) = p(nu,1)
              pu(nu) = p(nu,2)
          enddo
      endif
      do nu=0,3
          pl(nu) = p(nu,3)
          pn(nu) = p(nu,4)
      enddo

*
** mandelstam invariants
*
      do nu=0,3
          ptmps(nu) = pu(nu) + pd(nu)
          ptmpt(nu) = pu(nu) - pn(nu)
          ptmpu(nu) = pd(nu) - pn(nu)
      enddo
      s = dotp(ptmps,ptmps)
      t = dotp(ptmpt,ptmpt)
      u = dotp(ptmpu,ptmpu)
*
      complexmasses = .false.
*
      call deltaww  (s,delww)

      call deltawdu (s,delwdu)

      call deltawnul(s,delwnul)

      call deltabox (s,t,u,delbox)
*
* eq. (2.10) of Dittmaier-Kraemer PRD65 073007
*
      delvirt = delww + delwdu + delwnul + delbox

      return
      end subroutine deltavirt
*
** Subroutine for the calculation of the 
** contribution to delvirt due to W self-energy
**
** (pu) u  \            / l+ 
**          \    l     /
**           \___/\___/    (etc)
**           /   \/   \ 
**          /    l'~   \
** (pd) d~ /            \ nu_l
**
**        s = (pu+pd)**2
*
      subroutine deltaww(s,delww)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
*
      real*8 s
      complex*16 delww
*
      complex*16 swtmw2
      complex*16 swts
      complex*16 swtpmw2
      complex*16 deltamw2
      complex*16 deltazw
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save deltamw2
      save deltazw
*
      if (ifirst.eq.0) then
          ifirst = 1
*
          call sigmawt(dble(mw2)*cone,swtmw2)
          call sigmawtp(dble(mw2)*cone,swtpmw2)
*
* eq. (2.11) of Dittmaier-Kraemer + eq. (3.19) of ArXiv:0709.1075 
* (Denner Fortschritte) and eq. (4.29) of ArXiv:0505042
*
          deltamw2  =  swtmw2 + ii*dimag(mw2)*swtpmw2
     +              - ii*dimag(mw2)*alpha/pi
*
* eq. (3.19) of ArXiv:0709.1075 (Denner Fortschritte)
*
          deltazw  = - swtpmw2
*
      endif

      call sigmawt(s*cone,swts)

      delww   = - ( swts - deltamw2 ) / ( s*cone - mw2 ) - deltazw

      return
      end subroutine deltaww
*
**
** Subroutine for the calculation of the 
** contribution to deltavirt due to triangles on
** initial state
**
** (pu) u  \            / l+ 
**          \          /
**         W|\________/    (etc)
**          |/        \ 
**          /          \
** (pd) d~ /            \ nu_l
**
**        s = (pu+pd)**2
*
      subroutine deltawdu(s,delwdu)
*
* eq. (2.12) of Dittmaier-Kraemer PRD65 073007
*
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
*
      real*8 s
      complex*16 delwdu
*
      complex*16 getdeltar
      external getdeltar
*
      integer scheme
      common/sch/scheme
*
      complex*16 fwdu
      complex*16 delwdu_ct_alpha0
      complex*16 delwdu_ct
      complex*16 deltar
      complex*16 delze
      complex*16 delswow
      complex*16 delzw
      complex*16 delzdl
      complex*16 delzul
 
      complex*16 swtmw2
      complex*16 szztmz2
      complex*16 sazt0
      complex*16 saatp0
      complex*16 swtpmw2
      complex*16 szztpmz2

      complex*16 deltamw2
      complex*16 deltamz2
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save delwdu_ct

      if (ifirst.eq.0) then
          ifirst = 1
*
* q_f and q_f' must respect the following relation: q_f - q_f' = 2 i_3_f
*
          call sigmawt(dble(mw2)*cone,swtmw2)
          call sigmazzt(dble(mz2)*cone,szztmz2)
          call sigmazztp(dble(mz2)*cone,szztpmz2)
          call sigmaazt(zero,sazt0)
          call sigmaaatp(zero,saatp0)
          call sigmawtp(dble(mw2)*cone,swtpmw2)
*
* delze according to eq. (3.32) of ArXiv:0709.1075 (Denner Fortschritte)
* and eq. (4.9) of ArXiv:0505042
*
          delze  = 0.5d0*saatp0 - sw/cw * sazt0/mz2
*
* eq. (4.29) of ArXiv:0505042
*
          deltamw2  = swtmw2 + ii*dimag(mw2)*swtpmw2
     +              - ii*dimag(mw2)*alpha/pi

          deltamz2  = szztmz2 + ii*dimag(mz2)*szztpmz2
*
* delsw according to eq. (3.35) of ArXiv:0709.1075 (Denner Fortschritte)
* and eq. (4.8) of ArXiv:0505042
*
          delswow  = -0.5d0*cw2/sw2*( deltamw2/mw2 - deltamz2/mz2 )
*
* eq. (3.19) of ArXiv:0709.1075 (Denner Fortschritte) and 4.10 of ArXiv:0505042
*
          delzw = - swtpmw2
*
* delzfl: wave functions renormalization according to eq. (3.20) 
*         of ArXiv:0709.1075 (Denner Fortschritte). only 
*         diagonal terms are considered since ckm matrix = 1 
* here f= up and down quark (if we change final state this should 
*                            be changed accordingly)
*
          call deltazfl(qd,gdm,gdp,zero,zero,delzdl)
          call deltazfl(qu,gum,gup,zero,zero,delzul)
*
* eq. (2.13) of Dittmaier-Kraemer PRD65 073007
*
          delwdu_ct_alpha0 = delze-delswow+0.5d0*(delzw+delzdl+delzul)
*
* eq. (2.17) of Dittmaier-Kraemer PRD65 073007
*
          if(scheme.eq.0) then

              delwdu_ct = delwdu_ct_alpha0

          else
*
* deltar according to eq. (8.14) of ArXiv:0709.1075 (Denner Fortschritte)
*
              deltar = getdeltar()

              delwdu_ct = delwdu_ct_alpha0 - 0.5d0*deltar

          endif

      endif
*
      call wffp(s,qu,qd,0d0,0d0,fwdu)
*
* eq. (2.12) of Dittmaier-Kraemer
*
      delwdu = fwdu + delwdu_ct

      return
      end subroutine deltawdu
*
**
** Subroutine for the calculation of the 
** contribution to deltavirt due to triangles on
** final state
**
** (pu) u  \            / l+ 
**          \          /
**           \________/|   (etc)
**           /        \|W 
**          /          \
** (pd) d~ /            \ nu_l
**
**        s = (pu+pd)**2
*
      subroutine deltawnul(s,delwnul)
*
* eq. (2.12) of Dittmaier-Kraemer PRD65 073007
*
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
*
      real*8 s
      complex*16 delwnul
*
      complex*16 getdeltar
      external getdeltar
*
      real*8 mlep2
      common/leptmass/mlep2

      integer scheme
      common/sch/scheme
*
      complex*16 fwnul
      complex*16 delwnul_ct_alpha0
      complex*16 delwnul_ct
      complex*16 deltar
      complex*16 delze
      complex*16 delswow
      complex*16 delzw
      complex*16 delzll
      complex*16 delznul

      complex*16 swtmw2
      complex*16 szztmz2
      complex*16 swtpmw2
      complex*16 sazt0
      complex*16 saatp0
      complex*16 szztpmz2

      complex*16 deltamw2
      complex*16 deltamz2
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save delwnul_ct

      if (ifirst.eq.0) then
          ifirst = 1
*
* q_f and q_f' must respect the following relation: q_f - q_f' = 2 i_3_f
*
          call sigmawt(dble(mw2)*cone,swtmw2)
          call sigmazzt(dble(mz2)*cone,szztmz2)
          call sigmazztp(dble(mz2)*cone,szztpmz2)
          call sigmaazt(zero,sazt0)
          call sigmaaatp(zero,saatp0)
          call sigmawtp(dble(mw2)*cone,swtpmw2)
*
* delze according to eq. (3.32) of ArXiv:0709.1075 (Denner Fortschritte)
* and eq. (4.9) of ArXiv:0505042
*
          delze = 0.5d0*saatp0 - sw/cw * sazt0/mz2
*
* eq. (4.29) of ArXiv:0505042
*
          deltamw2 = swtmw2 + ii*dimag(mw2)*swtpmw2
     +              - ii*dimag(mw2)*alpha/pi

          deltamz2 = szztmz2 + ii*dimag(mz2)*szztpmz2
*
* delsw according to eq. (3.35) of ArXiv:0709.1075 (Denner Fortschritte)
* and eq. (4.8) of ArXiv:0505042
*
          delswow  = -0.5d0*cw2/sw2*( deltamw2/mw2 - deltamz2/mz2 )
*
*   eq. (3.19) of ArXiv:0709.1075 (Denner Fortschritte) 
* + eq. (4.10) of ArXiv:0505042
*
          delzw = - swtpmw2
*
* delzfl: wave functions renormalization according to eq. (3.20) 
*         of ArXiv:0709.1075 (Denner Fortschritte). only 
*         diagonal terms are considered since ckm matrix = 1 
* here f= electron and neutrino (if we change final state this should 
*                                be changed accordingly)
*
          call deltazfl(ql,glm,glp,mlep2*cone,zero,delzll)
          call deltazfl(qnu,gnm,gnp,zero,mlep2*cone,delznul)
*
* eq. (2.13) of Dittmaier-Kraemer PRD65 073007
*

          delwnul_ct_alpha0=delze-delswow+0.5d0*(delzw+delzll+delznul)
*
* eq. (2.17) of Dittmaier-Kraemer PRD65 073007
*
          if(scheme.eq.0) then

              delwnul_ct = delwnul_ct_alpha0

          else
*
* deltar according to eq. (8.14) of ArXiv:0709.1075 (Denner Fortschritte)
*
              deltar = getdeltar()

              delwnul_ct = delwnul_ct_alpha0 - 0.5d0*deltar

          endif

      endif
*
      call wffp(s,qnu,ql,0.d0,mlep2,fwnul)
*
* eq. (2.12) of Dittmaier-Kraemer
*
      delwnul = fwnul + delwnul_ct

      return
      end subroutine deltawnul
*
**
** Vertex correction form factor 
** Eq. (A.1) of Dittmaier-Kraemer PRD65 073007
**
*
      subroutine wffp(s,qf,qfp,mf2,mfp2,fwffp)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real*8 s,qf,qfp,mf2,mfp2
      complex*16 fwffp
*
      complex*16 b0d00mw,b0d00mz,b0ds0mw,b0dsmwmz,b0ds00,b0dmf20mf
      complex*16 b0dmfp20mfp

      complex*16 c0dmf2mfp2smfmgmfp,c0dmf20s0mfmw
      complex*16 c0dmfp20s0mfpmw,c0d00smw0mz
      complex*16 c0d00s0mz0 
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d00mw
      save b0d00mz
*
* q_f and q_f must respect the following relation: q_f - q_f' = 2 i_3_f
*
      if (ifirst.eq.0) then

          ifirst = 1
          call b0p0m0(mw2,b0d00mw)
          call b0p0m0(mz2,b0d00mz)

      endif

      call b0m0(s,mw2,b0ds0mw)
      call b0reg(s,mw2,mz2,b0dsmwmz)
      call b0m0m0(s,b0ds00)

      if (abs(qf).lt.epsilon) then
          b0dmf20mf = zero
      else
          call b0m0(mf2*cone,mf2*cone,b0dmf20mf)
      endif

      if (abs(qfp).lt.epsilon) then
          b0dmfp20mfp = zero
      else
          call b0m0(mfp2*cone,mfp2*cone,b0dmfp20mfp)
      endif
*
* C0 form factors
*
*
      if(abs(qf).lt.epsilon.or.abs(qfp).lt.epsilon) then
*        in this case the contribution is zero
*        (neutrino not coupled to gamma)
         c0dmf2mfp2smfmgmfp = zero
      elseif(abs(qf-qu).lt.epsilon) then
         call c0ircoll_qqp(mf2,mfp2,s,0d0,
     +                     c0dmf2mfp2smfmgmfp)
*
      else
         print*,'problems in quantum numbers in wffp'
      endif
*
      if(abs(qf).lt.epsilon) then
*        the contribution comes only from the charged lepton fp
*
* the sequence of external p2 and internal masses is understood as 
* mfp2,0,s,0,mfp2,mw^2
*
* here mw used is complex
*
         call c0cl(mfp2,s,mw2,c0dmfp20s0mfpmw)

         c0dmf20s0mfmw = zero

      elseif(abs(qfp).lt.epsilon) then
*        the contribution comes only from the charged lepton f
*
* the sequence of external p2 and internal masses is understood as 
* mf2,0,s,0,mf2,mW^2
*
* here mW used is complex
*
         call c0cl(mf2,s,mw2,c0dmf20s0mfmw)

         c0dmfp20s0mfpmw = zero

      else
* contributions from quark up and down
         call c0cq(mf2,s,mw2,c0dmf20s0mfmw)
*
         call c0cq(mfp2,s,mw2,c0dmfp20s0mfpmw)
*
      endif
*
      call c0reg1(s*cone,mw2,mz2,c0d00smw0mz)

      call c0reg2(s*cone,mz2,c0d00s0mz0)

***

      fwffp  = alsu4pi/4.d0*(
     +          2.d0*cone*(cone-2.d0*sw2-4.d0*qf*qfp*sw2)
     +                                          /cw2/sw2
     +       + 4.d0/sw2*(2.d0*cone+mw2/s)*b0d00mw
     +       + 2.d0/cw2/sw2*(cone-2.d0*sw2+2.d0*qf*qf*sw4
     +                       +2.d0*qfp*qfp*sw4)
     +                     *(2.d0*cone+mz2/s)*b0d00mz
     +       - 4.d0*(cone+mw2/s)*b0ds0mw
     +       - 4.d0*cw2/sw2*(cone+(mw2+mz2)/s)*b0dsmwmz
     +       + 1.d0/cw2/sw2*(2.d0*mz2/s*(qf-qfp-2.d0*qf*sw2)
     +                                 *(qf-qfp+2.d0*qfp*sw2)
     +               +3.d0*(cone-2.d0*sw2-4.d0*qf*qfp*sw2))*b0ds00
     +       + 8.d0*qf*qf*b0dmf20mf
     +       + 8.d0*qfp*qfp*b0dmfp20mfp
     +       - 8.d0*qf*qfp*s*c0dmf2mfp2smfmgmfp
     +       + 8.d0*qf*(qf-qfp)*mw2*c0dmf20s0mfmw
     +       - 8.d0*qfp*(qf-qfp)*mw2*c0dmfp20s0mfpmw
     +       + 8.d0*cw2/sw2*(mw2+mz2+mw2*mz2/s)*c0d00smw0mz
     +       + 2.d0/cw2/sw2*(qf-qfp-2.d0*qf*sw2)*(qf-qfp+2.d0*qfp*sw2)
     +                     *(cone+mz2/s)**2 *s*c0d00s0mz0 )

      return
      end subroutine wffp
*
**
** Subroutine for the calculation of the 
** contribution to deltavirt due to boxes
**
** (pu) u  \          / l+  (pl)
**          \ ______ /
**           |      |   (etc)
**           |______| 
**          /        \
** (pd) d~ /          \ nu_l (pn)
**
**        s = (pu+pd)**2
**        t = (pu-pl)**2 = (pd-pn)**2
**        u = (pu-pn)**2 = (pd-pl)**2
*
* eq. (A.6) of Dittmaier-Kraemer PRD65 073007
*
      subroutine deltabox(s,t,u,delbox)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real*8 s,t,u
      complex*16 delbox
*
      real*8 mlep2
      common/leptmass/mlep2
*
      complex*16 b0dt00,b0dsmwmz,b0ds0mw
      complex*16 c0d00smw0mz,c0d00t0mw0,c0d00t0mz0,c0dmd20s0mdmw,
     +           c0dml20s0mlmw,
     +           c0dmd2ml2tmdmgml,c0dmu20s0mumw
      complex*16 d0d0000ts0mw0mz,d0d0000us0mw0mz,
     +           d0dmd2ml200tsmdmgmlmw,d0dmu2ml200usmumgmlmw

      complex*16 mw2in,mz2in
*
*  remember: s+t+u=0! (if m=0)
*
      call b0m0m0(t*cone,b0dt00)
      call b0reg(s*cone,mw2,mz2,b0dsmwmz)
      call b0m0(s*cone,mw2,b0ds0mw)
     
*
* regular C0
*
      call c0reg1(s*cone,mw2,mz2,c0d00smw0mz)
*
      call c0reg2(t*cone,mw2,c0d00t0mw0)
*
      call c0reg2(t*cone,mz2,c0d00t0mz0)
*
* singular C0
*
*
* the sequence of external p2 and internal masses is understood as 
* mlep2,0,s,0,mlep2,mw^2
*
* here mw used is complex
*
      call c0cl(mlep2,s,mw2,c0dml20s0mlmw)
*
* the sequence of external p2 and internal masses is understood as 
* md2,0,s,0,md2,mw^2
*
* here mw used is complex
*
      call c0cq(0d0,s,mw2,c0dmd20s0mdmw)
*
* the sequence of external p2 and internal masses is understood as 
* mu2,0,s,0,mu2,mw^2
*
* here mw used is complex
*
      call c0cq(0d0,s,mw2,c0dmu20s0mumw)
*
* the sequence of external p2 and internal masses is understood as 
*     md2,mlep2,t,md2,lambda2,mlep2
*
      call c0ircoll_ql(0d0,mlep2,t,0d0,c0dmd2ml2tmdmgml)
*
* regular D0
*
* D0(p102,p212,p322,p302,p202,p312,m02,m12,m22,m32) = 
* D0(p102,p202,p322,p312,p212,p302,m12,m02,m22,m32)
*
* d0reg works only with complex masses
*
      call d0reg(0.d0,t,0.d0,s,0.d0,0.d0,
     +           mw2,zero,zero,mz2,
     +           d0d0000ts0mw0mz)

*
* regular D0
*
      call d0reg(0.d0,u,0.d0,s,0.d0,0.d0,
     +           mw2,zero,zero,mz2,
     +           d0d0000us0mw0mz)
*
* singular D0
* D0(md2,ml2,0,0,t,s,md,mg,ml,mw)
*    we use the relation 
* D0(md2,ml2,0,0,t,s,md,mg,ml,mw) = D0(md2,0,0,ml2,s,t,mg,md,mw,ml)
*
      call d0sing(0d0,0.d0,0.d0,mlep2,s,t,
     +            0d0,0d0,mw2,mlep2,
     +            d0dmd2ml200tsmdmgmlmw)
*
* singular D0
* D0(mu2,ml2,0,0,u,s,mu,mg,ml,mw)
*    we use the relation 
* D0(mu2,ml2,0,0,u,s,mu,mg,ml,mw) = D0(mu2,0,0,ml2,s,u,mg,mu,mw,ml)
*
      call d0sing(0d0,0.d0,0.d0,mlep2,s,u,
     +            0d0,0d0,mw2,mlep2,
     +            d0dmu2ml200usmumgmlmw)
*
      delbox= alsu4pi * (s*cone - mw2)
     +       * ( (gdm*glm + gnm*gum) * (2.d0/u*(b0dt00-b0dsmwmz)
     +                                   +4.d0*c0d00smw0mz
     +              -(mw2+mz2+s*cone+2.d0*t*cone)/u/u
     +               *(t*c0d00t0mw0+t*c0d00t0mz0
     +                                      +2.d0*s*c0d00smw0mz)
     +             -(t*(mw2+mz2+s*cone+2.d0*t*cone)**2
     +               +2.d0*(mw2+t*cone)*(mz2+t*cone)*u)
     +                                  /u/u*d0d0000ts0mw0mz)
     +          -2.d0*(gdm*gnm+glm*gum) * (2.d0*c0d00smw0mz
     +                                      -u*d0d0000us0mw0mz)
     +          +qd*ql*(2.d0/u*(b0dt00-b0ds0mw)+2.d0*c0dmd20s0mdmw
     +                                         +2.d0*c0dml20s0mlmw
     +             -(mw2+s*cone+2.d0*t*cone)/u/u*(s*c0dmd20s0mdmw
     +                                 +s*c0dml20s0mlmw
     +                                 +t*c0d00t0mw0
     +                                 +t*c0dmd2ml2tmdmgml)
     +             -t*(cone+(mw2+t*cone)**2/u/u)*d0dmd2ml200tsmdmgmlmw)
     +          -2.d0*qu*ql*(c0dml20s0mlmw + c0dmu20s0mumw 
     +                       -u*d0dmu2ml200usmumgmlmw) )

      return
      end subroutine deltabox

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
      include 'mathx.h'
*
      real*8 qf,mf2,mfp2
      complex*16 gfm,gfp
      complex*16 dzflout
*
      complex*16 sfl,sflp,sfrp,sfsp
*
      call sigmafl (qf,gfm,gfp,mf2,mfp2,mf2,sfl )
      call sigmaflp(qf,gfm,gfp,mf2,mfp2,mf2,sflp)
      call sigmafrp(qf,gfm,gfp,mf2,mfp2,mf2,sfrp)
      call sigmafsp(qf,gfm,gfp,mf2,mfp2,mf2,sfsp)

      dzflout = - sfl - mf2*( sflp + sfrp + 2.d0*sfsp )

      return
      end subroutine deltazfl
*
* generic fermion self-energies
*
* eq. (B.6) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmafl(qf,gfm,gfp,mf2,mfp2,s,sflout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real*8 qf,mf2,mfp2,s
      complex*16 gfm,gfp
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

      sflout  = - alsu4pi * (   qf**2 * (2.d0*b1dsmfmg + cone)
     +             + gfm**2 * (2.d0*b1dsmfmz + cone)
     +             + 0.5d0/sw2*mf2/2.d0/mw2*(b1dsmfmz + b1dsmfmh)
     +             + 0.5d0/sw2*((2.d0*cone+mfp2/mw2)*b1dsmfpmw+cone)   )

***
      return
      end subroutine sigmafl
*
* generic fermion self-energies derivatives w.r.t. p^2
*
* derivative (w.r.t. p^2 (s)) of eq. (B.6) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmaflp(qf,gfm,gfp,mf2,mfp2,s,sflpout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real*8 qf,mf2,mfp2,s
      complex*16 gfm,gfp
      complex*16 sflpout
*
      complex*16 b1pdsmfmg,b1pdsmfmz,b1pdsmfmh,b1pdsmfpmw

      if (abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b1pdsmfmg = zero
      
      else

          call b1pir(s*cone,mf2*cone,b1pdsmfmg)

      endif

      call b1preg(s*cone,mf2*cone,mz2,b1pdsmfmz)
      call b1preg(s*cone,mf2*cone,mh2,b1pdsmfmh)
      call b1preg(s*cone,mfp2*cone,mw2,b1pdsmfpmw)

      sflpout  = - alsu4pi * ( qf**2  * (2.d0*b1pdsmfmg)
     +              + gfm**2 * (2.d0*b1pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b1pdsmfmz + b1pdsmfmh)
     +              + 0.5d0/sw2*(2.d0*cone+mfp2/mw2)*b1pdsmfpmw )

      return
      end subroutine sigmaflp
*
* derivative (w.r.t. p^2 (s)) of eq. (B.7) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmafrp(qf,gfm,gfp,mf2,mfp2,s,sfrpout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real*8 qf,mf2,mfp2,s
      complex*16 gfm,gfp
      complex*16 sfrpout
*      
      complex*16 b1pdsmfmg,b1pdsmfmz,b1pdsmfmh,b1pdsmfpmw

      if (abs(qf).lt.epsilon) then
      ! terms only multiplied by qf
          b1pdsmfmg = zero
      
      else

          call b1pir(s*cone,mf2*cone,b1pdsmfmg)

      endif

      call b1preg(s*cone,mf2*cone,mz2,b1pdsmfmz)
      call b1preg(s*cone,mf2*cone,mh2,b1pdsmfmh)
      call b1preg(s*cone,mfp2*cone,mw2,b1pdsmfpmw)

      sfrpout  = - alsu4pi * ( qf**2 * (2.d0*b1pdsmfmg)
     +              + gfp**2 * (2.d0*b1pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b1pdsmfmz + b1pdsmfmh)
     +              + 0.5d0/sw2*mfp2/mw2*b1pdsmfpmw )

      return
      end subroutine sigmafrp
*
* derivative (w.r.t. p^2 (s)) of eq. (B.8) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmafsp(qf,gfm,gfp,mf2,mfp2,s,sfspout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real*8 qf,mf2,mfp2,s
      complex*16 gfm,gfp
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

      sfspout  = - alsu4pi * ( qf**2 * 4.d0*b0pdsmfmg
     +              + gfp*gfm * (4.d0*b0pdsmfmz)
     +              + 0.5d0/sw2*mf2/2.d0/mw2*(b0pdsmfmz - b0pdsmfmh)
     +              + 0.5d0/sw2*mfp2/mw2*b0pdsmfpmw )

      return
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
      include 'mathx.h'
      include 'pwhg_physpar.h'
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
      save b0d0meme,
     +     b0d0mmmm,
     +     b0d0mtlmtl,
     +     b0d0mumu,
     +     b0d0mdmd,
     +     b0d0mcmc,
     +     b0d0msms,
     +     b0d0mtmt,
     +     b0d0mbmb,
     +     b0d0mwmw


      if (ifirst.eq.0) then
          ifirst=1
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
     +     -ql*(glp+glm)*( - (s+2.d0*me2*cone)*b0dsmeme
     +                     +2.d0*me2*b0d0meme + 1.d0/3.d0*s 
     +                     - (s+2.d0*mm2*cone)*b0dsmmmm
     +                     +2.d0*mm2*b0d0mmmm + 1.d0/3.d0*s 
     +                     - (s+2.d0*mtl2*cone)*b0dsmtlmtl
     +                     +2.d0*mtl2*b0d0mtlmtl + 1.d0/3.d0*s )
* sum over quarks
     +     -3.d0*qu*(gup+gum)*( - (s+2.d0*mu2*cone)*b0dsmumu
     +                          +2.d0*mu2*b0d0mumu + 1.d0/3.d0*s 
     +                          - (s+2.d0*mc2*cone)*b0dsmcmc
     +                          +2.d0*mc2*b0d0mcmc + 1.d0/3.d0*s 
     +                          - (s+2.d0*mt2*cone)*b0dsmtmt
     +                          +2.d0*mt2*b0d0mtmt + 1.d0/3.d0*s )
     +     -3.d0*qd*(gdp+gdm)*( - (s+2.d0*md2*cone)*b0dsmdmd
     +                          +2.d0*md2*b0d0mdmd + 1.d0/3.d0*s 
     +                          - (s+2.d0*ms2*cone)*b0dsmsms
     +                          +2.d0*ms2*b0d0msms + 1.d0/3.d0*s 
     +                          - (s+2.d0*mb2*cone)*b0dsmbmb
     +                          +2.d0*mb2*b0d0mbmb + 1.d0/3.d0*s )  )
* bosonic part
     +    -1.d0/3.d0/sw/cw*( ((9.d0*cw2+0.5d0*cone)*s
     +                        +(12.d0*cw2+4.d0*cone)*mw2)*b0dsmwmw
     +                      -(12.d0*cw2-2.d0*cone)*mw2*b0d0mwmw
     +                      +1.d0/3.d0*s )
      saztout  = - alsu4pi * saztout

***

      return
      end subroutine sigmaazt
*
******************************************************************
*
* eq. (B.3) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmazzt(s,szztout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
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
      save b0d0meme,b0d0mmmm,b0d0mtlmtl,b0d0mumu,b0d0mdmd,
     +     b0d0mcmc,b0d0msms,b0d0mtmt,b0d0mbmb,b0d0mzmz,
     +     b0d0mhmh,b0d0mwmw,b0d0mzmh

*
      if (ifirst.eq.0) then
          ifirst=1

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
     +      3.d0*(gnp**2+gnm**2)*( - (s)*b0ds00 + 1.d0/3.d0*s)
     +     +(glp**2+glm**2)*(- (s+2.d0*me2*cone)*b0dsmeme
     +                  +2.d0*me2*b0d0meme + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*me2*b0dsmeme
     +     +(glp**2+glm**2)*(- (s+2.d0*mm2*cone)*b0dsmmmm
     +                  +2.d0*mm2*b0d0mmmm + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mm2*b0dsmmmm
     +     +(glp**2+glm**2)*(- (s+2.d0*mtl2*cone)*b0dsmtlmtl
     +                  +2.d0*mtl2*b0d0mtlmtl + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mtl2*b0dsmtlmtl 
* sum over quarks
     +     +3.d0*( (gup**2+gum**2)*(- (s+2.d0*mu2*cone)*b0dsmumu
     +                  +2.d0*mu2*b0d0mumu + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*mu2*b0dsmumu 
     +           + (gdp**2+gdm**2)*(- (s+2.d0*md2*cone)*b0dsmdmd
     +                  +2.d0*md2*b0d0mdmd + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*md2*b0dsmdmd 
     +           + (gup**2+gum**2)*(- (s+2.d0*mc2*cone)*b0dsmcmc
     +                  +2.d0*mc2*b0d0mcmc + 1.d0/3.d0*s) 
     +                  +3.d0/4.d0/sw2/cw2*mc2*b0dsmcmc 
     +           + (gdp**2+gdm**2)*(- (s+2.d0*ms2*cone)*b0dsmsms
     +                  +2.d0*ms2*b0d0msms + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*ms2*b0dsmsms 
     +           + (gup**2+gum**2)*(- (s+2.d0*mt2*cone)*b0dsmtmt
     +                  +2.d0*mt2*b0d0mtmt + 1.d0/3.d0*s)
     +                  +3.d0/4.d0/sw2/cw2*mt2*b0dsmtmt 
     +           + (gdp**2+gdm**2)*(- (s+2.d0*mb2*cone)*b0dsmbmb
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

***

      return
      end subroutine sigmazzt
*
******************************************************************
*
* derivative of eq. (B.3) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmazztp(s,szztpout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
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
      complex*16 b0pds00,b0pdsmeme,b0pdsmmmm,b0pdsmtlmtl,b0pdsmumu,
     +           b0pdsmdmd,b0pdsmcmc,b0pdsmsms,b0pdsmtmt,b0pdsmbmb,
     +           b0pdsmwmw,b0pdsmzmh,b0d0mzmh

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
     +      3.d0*(gnp**2+gnm**2)*(-s*b0pds00-b0ds00+1.d0/3.d0)
* sum over three lepton families
     +     +(glp**2+glm**2)*(- (s+2.d0*me2*cone)*b0pdsmeme
     +                  - (1.d0)*b0dsmeme + 1.d0/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*me2*b0pdsmeme
     +     +(glp**2+glm**2)*(- (s+2.d0*mm2*cone)*b0pdsmmmm
     +                  - (1.d0)*b0dsmmmm + 1.d0/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*mm2*b0pdsmmmm
     +     +(glp**2+glm**2)*(- (s+2.d0*mtl2*cone)*b0pdsmtlmtl
     +                  - (1.d0)*b0dsmtlmtl + 1.d0/3.d0)
     +                  +3.d0/4.d0/sw2/cw2*mtl2*b0pdsmtlmtl
* sum over quarks
     +     +3.d0*( (gup**2+gum**2)*(- (s+2.d0*mu2*cone)*b0pdsmumu
     +                  - b0dsmumu + 1.d0/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mu2*b0pdsmumu )
     +     +3.d0*( (gdp**2+gdm**2)*(- (s+2.d0*md2*cone)*b0pdsmdmd
     +                  - b0dsmdmd + 1.d0/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*md2*b0pdsmdmd )
     +     +3.d0*( (gup**2+gum**2)*(- (s+2.d0*mc2*cone)*b0pdsmcmc
     +                  - b0dsmcmc + 1.d0/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mc2*b0pdsmcmc )
     +     +3.d0*( (gdp**2+gdm**2)*(- (s+2.d0*ms2*cone)*b0pdsmsms
     +                  - b0dsmsms + 1.d0/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*ms2*b0pdsmsms )
     +     +3.d0*( (gup**2+gum**2)*(- (s+2.d0*mt2*cone)*b0pdsmtmt
     +                  - b0dsmtmt + 1.d0/3.d0) 
     +                  +3.d0/4.d0/sw2/cw2*mt2*b0pdsmtmt )
     +     +3.d0*( (gdp**2+gdm**2)*(- (s+2.d0*mb2*cone)*b0pdsmbmb
     +                  - b0dsmbmb + 1.d0/3.d0) 
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
     +            -2.d0/3.d0 )
      szztpout  = - alsu4pi * szztpout

***

      return
      end subroutine sigmazztp
*
******************************************************************
*
* eq. (B.4) of ArXiv:0709.1075 (Denner Fortschritte)
*
      subroutine sigmawt(s,swtout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
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
      complex*16 b0pds0me,b0pds0mm,b0pds0mtl,b0pdsmumd,
     +           b0pdsmcms,b0pdsmtmb,b0pdsmwmg,b0pdsmwmz,
     +           b0pdsmwmh
*
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d0meme,b0d00me,b0d0mmmm,b0d00mm,b0d0mtlmtl,b0d00mtl,
     +     b0d0mumu,b0d0mdmd,b0d0mumd,b0d0mcmc,b0d0msms,b0d0mcms,
     +     b0d0mtmt,b0d0mbmb,b0d0mtmb,b0d0mwmg,b0d0mwmz,b0d0mwmw,
     +     b0d0mzmz,b0d0mwmh,b0d0mhmh

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

***

      endif

      return
      end subroutine sigmawt
*
******************************************************************
*
* derivative (w.r.t. p^2 (s)) of eq. (B.1) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmaaatp(s,saatpout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      complex*16 s
      complex*16 saatpout
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
     +                  + 1.d0/3.d0*cone )
     +     +2.d0*ql*ql*( (-1.d0)*b0dsmmmm - (s+2.d0*mm2*cone)*b0pdsmmmm
     +                  + 1.d0/3.d0*cone )
     +     +2.d0*ql*ql*( (-1.d0)*b0dsmtlmtl
     +                    -(s+2.d0*mtl2*cone)*b0pdsmtlmtl
     +                  + 1.d0/3.d0*cone )
* sum over quarks
     +     +3.d0* 2.d0*qu*qu*( (-1.d0)*b0dsmumu 
     +                          - (s+2.d0*mu2*cone)*b0pdsmumu
     +                        + 1.d0/3.d0*cone )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmdmd 
     +                          - (s+2.d0*md2*cone)*b0pdsmdmd
     +                        + 1.d0/3.d0*cone )
     +     +3.d0* 2.d0*qu*qu*( (-1.d0)*b0dsmcmc 
     +                          - (s+2.d0*mc2*cone)*b0pdsmcmc
     +                        + 1.d0/3.d0*cone )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmsms 
     +                          - (s+2.d0*ms2*cone)*b0pdsmsms
     +                        + 1.d0/3.d0*cone )
     +     +3.d0* 2.d0*qu*qu*( (-1.d0)*b0dsmtmt 
     +                          - (s+2.d0*mt2*cone)*b0pdsmtmt
     +                        + 1.d0/3.d0*cone )
     +     +3.d0* 2.d0*qd*qd*( (-1.d0)*b0dsmbmb 
     +                          - (s+2.d0*mb2*cone)*b0pdsmbmb
     +                        + 1.d0/3.d0*cone ) )
* bosonic part
     +    +3.d0 * b0dsmwmw + (3.d0*s + 4.d0*mw2) * b0pdsmwmw
      saatpout  = - alsu4pi * saatpout

***

      return

      end subroutine sigmaaatp
*
******************************************************************
*
* derivative (w.r.t. p^2 (s)) of eq. (B.4) of ArXiv:0709.1075 
* (Denner Fortschritte)
*
      subroutine sigmawtp(s,swtpout)
      implicit none
      include 'mathx.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
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
      integer ifirst
      data ifirst /0/
      save ifirst
      save b0d0meme,b0d00me,b0d0mmmm,b0d00mm,b0d0mtlmtl,b0d00mtl,
     +     b0d0mumu,b0d0mdmd,b0d0mumd,b0d0mcmc,b0d0msms,b0d0mcms,
     +     b0d0mtmt,b0d0mbmb,b0d0mtmb,b0d0mwmg,b0d0mwmz,b0d0mwmh,
     +     b0d0mwmw,b0d0mzmz,b0d0mhmh


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

      if (.not.complexmasses) then

          call b0preg(s,mw2,zero,b0pdsmwmg)
          
      else

          if (abs(dimag(mw2)).lt.epsilon.and.abs(mw2-s).lt.epsilon) then

              call b0pir(s,mw2,b0pdsmwmg)
              print*,"W on-shell, divergence not subtracted"

          else

              call b0preg(s,mw2,zero,b0pdsmwmg)

          endif

      endif
      
      call b0preg(s,mw2,mz2,b0pdsmwmz)
      call b0preg(s,mw2,mh2,b0pdsmwmh)

***

      swtpout= 2.d0/3.d0/2.d0/sw2*(
* sum over three charged leptons
     +        -b0ds0me -(s-me2/2.d0*cone)*b0pds0me 
     +                        + 1.d0/3.d0*cone
     +        + me2**2/2.d0/s*(b0pds0me)
     +        - me2**2/2.d0/s/s*(b0ds0me-b0d00me)
     +        -b0ds0mm -(s-mm2/2.d0*cone)*b0pds0mm 
     +                        + 1.d0/3.d0*cone
     +        + mm2**2/2.d0/s*(b0pds0mm)
     +        - mm2**2/2.d0/s/s*(b0ds0mm-b0d00mm)
     +        -b0ds0mtl-(s-mtl2/2.d0*cone)*b0pds0mtl
     +                        + 1.d0/3.d0*cone
     +        + mtl2**2/2.d0/s*(b0pds0mtl)
     +        - mtl2**2/2.d0/s/s*(b0ds0mtl-b0d00mtl) )
* sum over three quark families
     +      + 2.d0/3.d0/2.d0/sw2 * 3.d0 * (
     +        -b0dsmumd -(s-(mu2+md2)*cone/2.d0)*b0pdsmumd
     +        + 1.d0/3.d0*cone
     +        + (mu2-md2)**2/2.d0/s * (b0pdsmumd)
     +        - (mu2-md2)**2/2.d0/s/s * (b0dsmumd - b0d0mumd)
     +        -b0dsmcms -(s-(mc2+ms2)*cone/2.d0)*b0pdsmcms
     +        + 1.d0/3.d0*cone
     +        + (mc2-ms2)**2/2.d0/s * (b0pdsmcms)
     +        - (mc2-ms2)**2/2.d0/s/s * (b0dsmcms - b0d0mcms)
     +        -b0dsmtmb -(s-(mt2+mb2)*cone/2.d0)*b0pdsmtmb
     +        + 1.d0/3.d0*cone
     +        + (mt2-mb2)**2/2.d0/s * (b0pdsmtmb)
     +        - (mt2-mb2)**2/2.d0/s/s * (b0dsmtmb - b0d0mtmb))
* bosonic part
     +        + 2.d0/3.d0*(
     +             (5.d0)*b0dsmwmg
     +           + (2.d0*mw2 + 5.d0*s)*b0pdsmwmg
     +           -mw2**2/s*(b0pdsmwmg)
     +           +mw2**2/s/s*(b0dsmwmg-b0d0mwmg)
     +           + 1.d0/3.d0*cone)
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
      
      return
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
      if (abs(p2).le.2d0*epsilon) then

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

      return
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
      if(abs(p2).le.2d0*epsilon) then

          b0out = zero

      else

          b0out = 2.d0*cone - log(-p2/mudim2)

      endif

      return
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

      return
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

      return
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

      return
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
*
      real*8 mlep2
      common/leptmass/mlep2
*
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
      if(abs(p2).le.2d0*epsilon) then

          call b0p0m0(m12,b0out)

      elseif (abs(m12-p2).le.2d0*epsilon) then

          b0out = 2.d0*cone - log(p2/mudim2)

      elseif (dble(m12).lt.abs(p2)) then

          if (.not.complexmasses.and.abs(m12-mw2).le.epsilon.and.
     +        abs(p2-mlep2).gt.epsilon) then
              arglog = ( p2 - m12 + ii*ph_WmWw )/p2
          else
              arglog = ( p2 - m12 )/p2
          endif

          b0out = 2.d0*cone - log(p2/mudim2)
     +          - m12/p2*log(m12/p2) 
     +          - (cone-m12/p2)*log(arglog)
     +          + ii*pi*(cone-m12/p2)
      
      else

          if (.not.complexmasses.and.abs(m12-mw2).le.epsilon.and.
     +        abs(p2-mlep2).gt.epsilon) then
              arglog = -( p2 - m12 + ii*ph_WmWw )/p2
          else
              arglog = -( p2 - m12 )/p2
          endif

          b0out = 2.d0*cone - log(p2/mudim2)
     +          - m12/p2*log(m12/p2) 
     +          - (cone-m12/p2)*log(arglog)

      endif
          
      return
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
      if (abs(p2).le.2d0*epsilon) then
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

      return
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
          if (.not.complexmasses.and.abs(m02-mw2).lt.2d0*epsilon) then
              arglog = ( m02 - ii*ph_WmWw - p2 )/m02
          else
              arglog = ( m02 - p2 )/m02
          endif

          if (dble(m02).le.dble(p2)) then
   
              b0pout = conjg(- m02*log(arglog)/(p2**2)-cone/p2)

          else

              b0pout = - m02*log(arglog)/(p2**2)-cone/p2

          endif
      
      endif

      return
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

      return
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

      return
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

      return
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
          return

      else
*
* eq. 5.54 bardin - passarino
*
          b1pout = + ( ( 3.d0 )*cone - log(m12/mudim2) )/(2.d0*m12)

      endif

      return
      end subroutine b1pir

* subroutine for the singular collinear C0 (and also threshold, 
* cured with complex w mass)
* according to Dittmaier-Kraemer PRD65 073007 
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (b.2)
*
* the sequence of external p2 and internal masses is understood as 
* m1^2,0,s,0,m1^2,m2^2, m2->mw n.b.: where m1 is the lepton mass
* this is not truely collinear singular because there is the lepton mass
*
* here mW used is complex
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

      if (.not.complexmasses.and.abs(m22-mw2).lt.2d0*epsilon) then
          arglog1= ( m22 - ii*ph_WmWw - s2bar )/m12
          arglog2= ( m22 - ii*ph_WmWw - s2bar )/m22
          argli21= - s2bar/( m22 - ii*ph_WmWw - s2bar )
      else
          arglog1= ( m22 - s2bar )/m12
          arglog2= ( m22 - s2bar )/m22
          argli21= - s2bar/( m22 - s2bar )
      endif

      argli22 = s2bar/m22

      c0out = - log(arglog1)*log(arglog2)
     +         - 2.d0*myli2(argli21)
     +         - myli2(argli22)
      c0out = c0out/(-s)
*
      return
      end subroutine c0cl
*
* subroutine for the singular collinear C0 (and also threshold, 
* cured with complex W mass)
* according to Dittmaier-Kraemer PRD65 073007 
* eq. (A.3) and Dittmaier, NPB 675 (2003) 447, eq. (B.2) and eq. (B.4)
*
* the sequence of external p2 and internal masses is understood as 
* m1^2,0,s,0,m1^2,m2^2, m2->mw n.b.: where m1 is the quark mass
*
* this is truely collinear singular since the quark mass -> 0
* here mw used is complex
*
* from non-abelian diagram
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
      if (.not.complexmasses.and.abs(m22-mw2).le.2d0*epsilon) then
          arglog = ( m22 - ii*ph_WmWw - s2bar )/m22
      else
          arglog = ( m22 - s2bar )/m22
      endif

      c0out= - log(mudim2/m22) * log(arglog)
     +       + ( log(arglog) )**2
     +       + myli2(s2bar/m22)

      c0out= c0out/s
*
      return
      end subroutine c0cq
*
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

      return
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
      real*8 p102in,p202in,p212in
      real*8 mg2in
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
      real*8 r
      complex*16 rbar
*
      r    = p212in
      rbar = r + ii*epsilon

*
* eq. (B.16) Dittmaier
*
      c0out = ( log(-mudim2/rbar) )**2 / 2d0 / r

      return
      end subroutine c0ircoll_qqp
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

      return
      end subroutine c0reg1
*
* eq. (5.59) of Bardin - Passarino
* C0(0,0,p2,0,m2,0)
*
      subroutine c0reg2(p2i,m2,c0out)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
*
      complex*16 p2i,m2
      complex*16 p2
      complex*16 c0out
*
      complex*16 myli2
      external myli2
*
* "-" is because B-P use a different metric 
      p2 = -dble(p2i)*cone
*
      c0out = - ( pi2/6d0 - myli2(cone-(p2/m2)) ) / p2

      return
      end subroutine c0reg2
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
*                  m02= 0d0, m12= mq2, m22= mw2, m32= mlep2
*
      subroutine d0sing(p102,p212,p322,p302,p202,p312,
     +                  m02,m12,m22,m32,
     +                  boxout)
      implicit none
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
*
      real*8 p102,p212,p322,p302,p202,p312
      real*8 m02,m12,m32
      complex*16 m22
      complex*16 boxout
*
      complex*16 myli2
      external myli2

      complex*16 li2cont
      external li2cont
*
      real*8 m3
      complex*16 m2,x32
      complex*16 li21lp1,li21lm1,argli2
      complex*16 arg1,arg2
      complex*16 arglog
*
      boxout= zero
*
      m2= sqrt(m22)
      m3= sqrt(m32)
      x32= m3/m2
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
      if (.not.complexmasses.and.abs(m22-mw2).le.2d0*epsilon) then
          arglog = m22 - ii*ph_WmWw - p202*cone
          argli2= (p212-p202)/(m22-ii*ph_WmWw-p202) * cone
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
      boxout = boxout / (p202-m22)/(p312-m32)
*
      return
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
      if(abs(dimag(m02in)).lt.epsilon) then
          m02 = m02in - ii*epsilon
      else
          m02 = m02in
      endif

      if(abs(dimag(m32in)).lt.epsilon) then
          m32 = m32in - ii*epsilon
      else
          m32 = m32in
      endif
*
      y01 = m02 - cone*p102
      y02 = m02 - cone*p202
      y03 = m02 + m32 - cone*p302
      y12 = - cone*p212
      y13 = m32 - cone*p312
      y23 = m32 - cone*p322
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

      return
      end subroutine d0reg
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

      return
      end function kallen
*
      real*8 function rkallen(x,y,z)
      implicit none
      real*8 x,y,z

      rkallen= sqrt(x**2+y**2+z**2-2.d0*(x*y+x*z+y*z))

      return
      end function rkallen
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

      return
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
      return
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

      do i=1,dim
         n    = 2*i
         add  = z**(n+1)*b(i)
         if (abs(add/myli23).gt.epsilon) then
             myli23 = myli23 + add
         else
             exit
         endif
      enddo

      return  
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
      arg= cone - arg1*arg2
      li2cont= myli2(arg) + myeta(arg1,arg2)*log(arg)
*
      return
      end function li2cont
*
* analytically continued dilogarithm of three variables
* according to eq. (2.13) of Denner-Dittmaier ArXiv:1005.2076
*
      complex*16 function li23arg(arg1,arg2,arg3)
      implicit none
      include 'mathx.h'
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
      return
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
      return
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
            if (dble(dconjg(qb)*sqdisc).lt.0.d0) si = -1.d0
         else
            argsq=1.d0 *
     +           (qb-2.d0*sqrt(qa*qc))*(qb+2.d0*sqrt(qa*qc))
            argsq = dcmplx(dble(argsq),dimag(- 4.d0*qa*qc))
            sqdisc = sqrt(argsq)
            si = 1.d0
            if (dble(dconjg(qb)*sqdisc).lt.0.d0) si = -1.d0
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
      return
      end subroutine qcroots
*
**
*
      real*8 function factorial(nl)
      integer n,nl
*
      n=nl
      factorial = 1.d0
      do while(n.gt.0)
         factorial = 1.d0*n * factorial
         n = n - 1
      enddo

      return
      end function factorial
*
**
*
      real*8 function binomial(n,k)
      integer n,k
*
      real*8 factorial
      external factorial

      binomial = factorial(n)/factorial(k)/factorial(n-k)

      return
      end function binomial
*
**
*
      real*8 function bernoulli(n)
      integer n
*
      real*8 binomial
      external binomial

      integer k,r
      
      if (n.eq.0) then
          bernoulli = 1d0
          return
      else
          bernoulli = 0d0
          k=0
          do while (k.le.n)
              add = 0d0
              r=0
              do while (r.le.k)
                  add = add + (-1)**r*binomial(k,r)*r**n
                  r=r+1
              enddo
              bernoulli = bernoulli + add/(k+1)
              k=k+1
          enddo
      endif

      return
      end function bernoulli

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
      call sigmaaatp(zero,saatp0)

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
