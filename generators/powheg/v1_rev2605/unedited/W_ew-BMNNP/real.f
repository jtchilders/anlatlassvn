      subroutine setreal(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
*
      real * 8 p(0:3,nlegreal)
      integer fermion_flav(nlegreal)
      real * 8 amp2
*
      if(fermion_flav(nlegreal).eq.22) then

          call setreal_ew(p,fermion_flav,amp2)

      else

          call setreal_st(p,fermion_flav,amp2)

      endif
*
      end


      subroutine setreal_st(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
c -*- Fortran -*-
c      character *2 flav(-5:5)
      real * 8 charge(-5:5)
c      data (charge(ijkh),ijkh=-5,5) 
c      data (flav(ijkh),ijkh=-5,5) 
c      data flav
c     #     /'b~','c~','s~','u~','d~','g','d','u','s','c','b'/
      data charge
     #     / 0.33333333333333333333d0, !   1d0/3
     #      -0.66666666666666666667d0, !  -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #      -0.66666666666666666667d0, !   -2d0/3
     #       0.33333333333333333333d0, !   1d0/3 
     #       0d0,                      !   0d0   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3   
     #      -0.33333333333333333333d0, !   -1d0/3
     #       0.66666666666666666667d0, !   2d0/3 
     #      -0.33333333333333333333d0/ !   -1d0/3
c      include 'QuarkFlavs.h'
      include 'PhysPars.h'
      integer nleg
      parameter (nleg=nlegreal)
      real * 8 p(0:3,nleg)
      integer fermion_flav(nleg)
      real * 8 amp2
      integer ferm_type(nleg)
      real * 8 ferm_charge(nleg)
      integer i,j,k,l,count,tmp_type
      real *8 tmp_charge
c     vector boson id and decay
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode  

c     i is the flavour index of first incoming parton
c     j is the flavour index of second incoming parton
c     k is the flavour of outgoing parton in the order particle,antiparticle,gluon
c     with the convention:
c     
c      -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6                    
c      t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t                    
      
      i = fermion_flav(1)
      j = fermion_flav(2)
      k = fermion_flav(5)
      ferm_charge(1) = charge(i)
      ferm_charge(2) = charge(j)
      ferm_charge(5) = charge(k)
      

      if (i.eq.0) then
         ferm_type(1) = 0
      else 
         ferm_type(1) = i/abs(i)
      endif 
      if (j.eq.0) then
         ferm_type(2) = 0
      else 
         ferm_type(2) = j/abs(j)
      endif   
      if (k.eq.0) then
         ferm_type(5) = 0
      else 
         ferm_type(5) = k/abs(k)
      endif   

c     antilepton-neutrino from W decay
      ferm_type(3) = fermion_flav(3)/abs(fermion_flav(3))
      ferm_charge(3) = ferm_type(3)*(-1d0)
      ferm_type(4) = -ferm_type(3)
      ferm_charge(4) = 0d0

      if(idvecbos.eq.24) then
         if (i.eq.0) then
c     g q -> W+ qp
            call g_aqp_to_al_vl_aq(p,ferm_type,ferm_charge,amp2) 
         elseif ((i.ne.0).and.(j.ne.0)) then
c     q aqp -> W+ g
            call q_aqp_to_al_vl_g(p,ferm_type,ferm_charge,amp2)
         elseif (j.eq.0) then
c     q g -> W+ qp
            call q_g_to_al_vl_qp(p,ferm_type,ferm_charge,amp2)
         else
            amp2 = 0d0
         endif
      elseif(idvecbos.eq.-24) then
         if (i.eq.0) then
c     g q -> W- qp
            call g_aqp_to_l_avl_aq(p,ferm_type,ferm_charge,amp2) 
         elseif ((i.ne.0).and.(j.ne.0)) then
c     q aqp -> W- g
            call q_aqp_to_l_avl_g(p,ferm_type,ferm_charge,amp2)
         elseif (j.eq.0) then
c     q g -> W- qp
            call q_g_to_l_avl_qp(p,ferm_type,ferm_charge,amp2)
         else
            amp2 = 0d0
         endif

      else
         write(*,*) 'ERROR: this subroutine deals only with W+ or W- '
         call exit(1)
      endif

      if (i.eq.0) i=abs(k)
      if (j.eq.0) j=abs(k)
      if(mod(abs(i),2).eq.0) then
         amp2=amp2*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      elseif(mod(abs(i),2).eq.1) then   
         amp2=amp2*ph_CKM(abs(j)/2,(abs(i)+1)/2)**2
      endif
c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2*pi))
      
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C AMPLITUDES RELATES BY CROSSING:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine g_aqp_to_al_vl_aq(pphy,fermion_type,fermion_charge,
     #     amp2)
      implicit none
      include 'nlegborn.h'
      integer nleg
      parameter (nleg=nlegreal)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 pp(0:3,nleg), ferm_charge(nleg)
      integer ferm_type(nleg)
      real * 8 amp2
      integer mu,i

      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

      do mu=0,3
c     exchange initial gluon <-> final quark
         pp(mu,5) = -pphy(mu,1)
         pp(mu,1) = -pphy(mu,5)
      enddo

c no useful information is in ferm_type(1) or ferm_charge(1), 
c since it's the gluon, BEFORE the following exchange
      ferm_type(1) = -ferm_type(5)
c NOTE the MINUS sign     !!!
      ferm_charge(1) = -ferm_charge(5)

c     if the following two lines are missing 
      ferm_type(5)=0
      ferm_charge(5)=0d0 
c     ferm_type(5) and ferm_charge(5) don't contain
c     their correct values, but this does not affect 
c     the correct call of

       call q_aqp_to_al_vl_g(pp,ferm_type,ferm_charge,
     #     amp2)

c     correct for color average
      amp2 = amp2 * 3d0/8d0
      
      end

ccccccccccccccccccccccccccccccccccccccccccccc

      subroutine q_g_to_al_vl_qp(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
      integer nleg
      parameter (nleg=5)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 pp(0:3,nleg), ferm_charge(nleg)
      integer ferm_type(nleg)
      real * 8 amp2
      integer mu,i

      do i = 1,nleg
         do mu=0,3
            pp(mu,i) = pphy(mu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

      do mu=0,3
c     exchange initial gluon <-> final quark
         pp(mu,5) = -pphy(mu,2)
         pp(mu,2) = -pphy(mu,5)
      enddo

c no useful information is in ferm_type(2) or ferm_charge(2), 
c since it's the gluon, BEFORE the following exchange
      ferm_type(2) = -ferm_type(5)
c NOTE the MINUS sign     !!!
      ferm_charge(2) = -ferm_charge(5)

c     if the following two lines are missing 
      ferm_type(5)=0
      ferm_charge(5)=0d0 
c     ferm_type(5) and ferm_charge(5) don't contain
c     their correct values, but this does not affect 
c     the correct call of

       call q_aqp_to_al_vl_g(pp,ferm_type,ferm_charge,
     #     amp2)

c     correct for color average
      amp2 = amp2 * 3d0/8d0
      
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine q_aqp_to_l_avl_g(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
      include 'nlegborn.h'
c the 5 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer nleg
      parameter (nleg=nlegreal)
      integer fermion_type(nleg),i
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 amp2
      real * 8 ferm_charge(nleg)
      integer ferm_type(nleg)

       if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with W- decay'
         stop
      endif

      do i=1,nleg
         ferm_charge(i) = -fermion_charge(i)
         ferm_type(i) = -fermion_type(i)
      enddo
            
      
      call q_aqp_to_al_vl_g(pphy,ferm_type,ferm_charge,
     #     amp2)

      end

ccccccccccccccccccccccccccccccccccccccccccccc

       subroutine g_aqp_to_l_avl_aq(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
      include 'nlegborn.h'
c the 5 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer nleg
      parameter (nleg=nlegreal)
      integer fermion_type(nleg),i
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 amp2
      real * 8 ferm_charge(nleg)
      integer ferm_type(nleg)

       if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with W- decay'
         stop
      endif

      do i=1,nleg
         ferm_charge(i) = -fermion_charge(i)
         ferm_type(i) = -fermion_type(i)
      enddo
      
      call g_aqp_to_al_vl_aq(pphy,ferm_type,ferm_charge,
     #     amp2)

      end

ccccccccccccccccccccccccccccccccccccccccccccc

       subroutine q_g_to_l_avl_qp(pphy,fermion_type,fermion_charge,
     #     amp2)
   
      implicit none
c the 5 4-momentum vectors
c p(i,1) is the i-th component of vector p1...   
      integer nleg
      parameter (nleg=5)
      integer fermion_type(nleg),i
      real * 8 fermion_charge(nleg)
      real * 8 pphy(0:3,nleg)
      real * 8 amp2
      real * 8 ferm_charge(nleg)
      integer ferm_type(nleg)

      if ((fermion_type(3).ne.1).and.(fermion_type(4).ne.-1)) then
         write(*,*) 'ERROR: this subroutine deals only with W- decay'
         stop
      endif

      do i=1,nleg
      
         ferm_charge(i) = -fermion_charge(i)
         ferm_type(i) = -fermion_type(i)
      enddo
      
      call q_g_to_al_vl_qp(pphy,ferm_type,ferm_charge,
     #     amp2)

      end
*
** Subroutine for strong real part
** mu = md = 0
**
**(pu) u \            / l+ (pl)
**        \---gl     /
**         \________/
**         /   W+   \       (gl can be emitted by both initial leg)
**        /          \
**(pd) d~/            \ nu_l (pn)
*
      subroutine q_aqp_to_al_vl_g(pphy,fermion_type,fermion_charge,amp2)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'pwhg_physpar.h'
*
      integer nleg
      parameter (nleg=nlegreal)
      real * 8 pphy(0:3,nleg)
      integer fermion_type(nleg)
      real * 8 fermion_charge(nleg)
      real * 8 amp2
*
      real*8 dotp
      external dotp
*
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
*
      real*8 p(0:3,nleg)
      integer ferm_type(nleg)
      real*8 ferm_charge(nleg)
      real*8 pu(0:3),pd(0:3),pn(0:3),pl(0:3),k(0:3)
      real*8 sp

      integer i,j,nu

      real*8 pupd,pupe,pupn,puk
      real*8 pdpe,pdpn,pdk
      real*8 pnpe,pnk
      real*8 pek
      complex*16 mw2m1,mw2m2,mw2dm1
      real*8 tmp
      real*8 ml2,ml4
      complex*16 gs2
      real*8 densp2

      integer ifirst
      data ifirst/0/
      save ifirst,ml2,ml4,mw2m1,mw2dm1,mw2m2

      if (ifirst.eq.0) then
          ifirst = 0

          ml2 = kn_masses(3)*kn_masses(3)
          ml4 = ml2*ml2
    
          mw2m1 = 1d0/mw2
          mw2dm1= conjg(mw2m1)
          mw2m2 = mw2m1*mw2dm1
      endif

      gs2 = 4d0*pi*st_alpha
*
c  copy of local variables
      do i=1,nlegreal
         do nu=0,3
            p(nu,i) = pphy(nu,i)
         enddo
         ferm_charge(i) = fermion_charge(i)
         ferm_type(i) = fermion_type(i)
      enddo

c     exchance particle 1 and 2
      if (ferm_type(1).eq.-1) then
         if (ferm_type(2).eq.1) then
            call exchange_momenta(p(0,1),p(0,2))
            tmp = ferm_charge(1)
            ferm_charge(1)=ferm_charge(2)
            ferm_charge(2)=tmp
            tmp = ferm_type(1)
            ferm_type(1)=ferm_type(2)
            ferm_type(2)=tmp
         else
            write(*,*) 'Error in the type of the quark 1-2'
            stop
         endif
      endif

      if (ferm_charge(1)*ferm_type(1).lt.0d0) then
          do nu=0,3
              pu(nu) = p(nu,2)
              pd(nu) = p(nu,1)
          enddo
      else
          do nu=0,3
              pu(nu) = p(nu,1)
              pd(nu) = p(nu,2)
          enddo
      endif

      do nu=0,3
          pl(nu) = p(nu,3)
          pn(nu) = p(nu,4)
          k (nu) = p(nu,5)
      enddo

*
      pupd = dotp(pu,pd)
      pupe = dotp(pu,pl)
      pupn = dotp(pu,pn)
      puk  = dotp(pu,k)
*
      pdpe = dotp(pd,pl)
      pdpn = dotp(pd,pn)
      pdk  = dotp(pd,k)
*
      pnpe = dotp(pl,pn)
      pnk  = dotp(pn,k)
*
      pek  = dotp(pl,k)
*
      sp = 2.d0*pnpe + ml2
*
      densp2 = 1d0/(sp*cone-ph_Wmass2+ii*ph_WmWw)/
     -             (sp*cone-ph_Wmass2-ii*ph_WmWw)

      amp2 =
     -   (densp2*(pdpe*pdpn*puk - pdpn*pek*pupd + 
     -      (-(pdpn*puk) + pnk*puk + 2*pdpn*pupd - pnk*pupd)*pupe + 
     -      pdk*(pdpn*(pek - pupe) + pupe*pupn)))/(pdk*puk)

      amp2 = amp2 * 16 
     +            * g2*conjg(g2) * gs2/4d0
     +            * CF/4d0/nc
* for initial state gluons
      if(dble(amp2).lt.0d0) amp2 = -amp2

      end

*
** Subroutine for real part
** mu = md = 0
**
**(pu) u \            / l+ (pl)
**        \          /--- g (k)
**         \________/
**         /   W+   \       (g can be emitted by every charged leg)
**        /          \
**(pd) d~/            \ nu_l (pn)
*
*
**
*

      subroutine setreal_ew(p,fermion_flav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'mathx.h'
      include 'pwhg_physpar.h'
*
      real * 8 p(0:3,nlegreal)
      integer fermion_flav(nlegreal)
      real * 8 amp2
*
      complex*16 e_
      external e_

      real*8 dotp
      external dotp
*
      real*8 mlep2
      common/leptmass/mlep2

      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode

*
      real*8 qu2,qd2

      real*8 pu(0:3),pd(0:3),pn(0:3),pl(0:3),k(0:3)
      real*8 s,sp

      integer i,j,nu

      real*8 pupd,pupe,pupn,puk
      real*8 pdpe,pdpn,pdk
      real*8 pnpe,pnk
      real*8 pek
      complex*16 mw2m1,mw2m2,mw2dm1
      complex*16 epupdpek,epupdpnk,epupnpek,epdpnpek,epupdpnpe
      real*8 ml2

      complex*16 dens,densc,densp,denspc,dens2,densp2
*
      integer ifirst
      data ifirst/0/
      save ifirst,ml2,mw2m1,mw2dm1,mw2m2

      if (ifirst.eq.0) then
          ifirst = 1

          ml2 = mlep2
    
          mw2m1 = 1d0/mw2
          mw2dm1= conjg(mw2m1)
          mw2m2 = mw2m1*mw2dm1
      endif

      i = fermion_flav(1)
      j = fermion_flav(2)

      if (mod(abs(i),2).eq.0) then
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
          k (nu) = p(nu,5)
      enddo
*
      ! CP symmetry for W-
      if (idvecbos.lt.0) then
          call invertspace(pu,pu)
          call invertspace(pd,pd)
          call invertspace(pl,pl)
          call invertspace(pn,pn)
          call invertspace(k,k)
      endif

*
      pupd = dotp(pu,pd)
      pupe = dotp(pu,pl)
      pupn = dotp(pu,pn)
      puk  = dotp(pu,k)
*
      pdpe = dotp(pd,pl)
      pdpn = dotp(pd,pn)
      pdk  = dotp(pd,k)
*
      pnpe = dotp(pl,pn)
      pnk  = dotp(pn,k)
*
      pek  = dotp(pl,k)
*
      epupdpek  = e_(pu,pd,pl,k)
      epupdpnk  = e_(pu,pd,pn,k)
      epupnpek  = e_(pu,pn,pl,k)
      epdpnpek  = e_(pd,pn,pl,k)
      epupdpnpe = e_(pu,pd,pn,pl)
*
      s  = 2.d0*pupd
      sp = 2.d0*pnpe + ml2

      dens   = 1d0/(s*cone - ph_Wmass2 + ii*ph_WmWw)
      densc  = dconjg(dens)
      densp  = 1d0/(sp*cone - ph_Wmass2 + ii*ph_WmWw)
      denspc = dconjg(densp)
      dens2  = dens*densc
      densp2 = densp*denspc
*
      amp2 =        (-4*epupnpek*(densp*puk*qd*
     -        (densc*(2*pdpe - 
     -             ml2*(denspc*mw2dm1*pek*(pdk + pupd) + 
     -                mw2m1*(pdpn - denspc*pek*(pdk + pupd))))
     -            + denspc*ml2*(-mw2dm1 + mw2m1)*pek*qd) + 
     -       dens*puk*(2*densc*(densp - denspc)*pdk*pdpe + 
     -          denspc*ml2*
     -           (densp*mw2m1*pek*(pdk + pupd) + 
     -             mw2dm1*(pdpn - densp*pek*(pdk + pupd)))*qd)
     -         + dens*denspc*pdk*
     -        (ml2*mw2dm1*(pdpe + pdpn) + 
     -          densp*ml2*
     -           (mw2m1*pdk + mw2dm1*(-pdk + pdpe + pdpn))*pek
     -            + 2*densp*pek*(pdpn + pupd))*qu - 
     -       densp*(densc*pdk*
     -           (ml2*mw2m1*(pdpe + pdpn) + 
     -             denspc*ml2*
     -              (mw2dm1*pdk + mw2m1*(-pdk + pdpe + pdpn))*
     -              pek + 2*denspc*pek*(pdpn + pupd)) + 
     -          2*denspc*ml2*(-mw2dm1 + mw2m1)*pek*pupd*qd)*qu
     -       ))/(pdk*pek*puk) + 
     -  (4*epdpnpek*(dens*densc*pdk*puk*
     -        (-(denspc*ml2*mw2dm1*puk) + 
     -          densp*ml2*mw2m1*puk - 2*densp*pupe + 
     -          2*denspc*pupe) + 
     -       denspc*puk*(densc*densp*pek*
     -           ((-2 + ml2*mw2m1)*pupe + ml2*mw2m1*pupn) - 
     -          dens*((-2 + densp*(-2 + ml2*mw2dm1)*pek)*
     -              pupe + densp*ml2*mw2dm1*pek*pupn))*qd + 
     -       ml2*(dens*denspc*mw2dm1 - densc*densp*mw2m1)*pdk*
     -        (puk + pupe)*qu + 
     -       densp*denspc*ml2*(-mw2dm1 + mw2m1)*pdk*pek*qu**2)
     -     )/(pdk*pek*puk) + 
     -  (4*epupdpek*(densc*densp*
     -        (-((-2 + ml2*mw2m1)*pdpn*(1 + 2*denspc*pek)) + 
     -          2*(1 + denspc*ml2*mw2m2*(pdpe + pdpn)*pek)*
     -           pnpe)*puk*qd - 
     -       densp*(densc*pdk*
     -           (ml2*mw2m1*pnpe + 
     -             denspc*pek*
     -              (2*pupn + 
     -                ml2*
     -                 (mw2dm1*(pdpn + pnk) + 
     -                   2*mw2m2*pnpe*(pupe + pupn) - 
     -                   mw2m1*(pdpn + pnk + 2*pupn)))) + 
     -          denspc*ml2*(mw2dm1 - mw2m1)*pek*(pdpn - pupn)*
     -           qd)*qu + 
     -       dens*(densc*pdk*puk*
     -           (denspc*((-2 + ml2*mw2dm1)*pdpn - 2*pnpe) + 
     -             densp*
     -              (pdpn*
     -                 (2 - ml2*mw2m1 + 
     -                   2*denspc*ml2*(mw2dm1 - mw2m1)*pek) + 
     -                2*pnpe + 
     -                denspc*ml2*(mw2dm1 - mw2m1)*pek*
     -                 (pnk + pnpe + 2*pupn))) + 
     -          denspc*((-2 + ml2*mw2dm1)*pdpn*
     -              (1 + 2*densp*pek) - 
     -             2*densp*ml2*mw2m2*(pdpe + pdpn)*pek*pnpe)*
     -           puk*qd + 
     -          denspc*pdk*
     -           (ml2*mw2dm1*pnpe + 
     -             densp*pek*
     -              (2*pupn + 
     -                ml2*
     -                 (mw2m1*(pdpn + pnk) + 
     -                   2*mw2m2*pnpe*(pupe + pupn) - 
     -                   mw2dm1*(pdpn + pnk + 2*pupn))))*qu)))
     -    /(pdk*pek*puk) + 
     -  (4*epupdpnpe*(dens*
     -        (densc*pdk*puk*
     -           (denspc*
     -              (-2*(puk + pupe) + 
     -                ml2*
     -                 (-1 + 
     -                   mw2dm1*
     -                    (pdpn + pek + pnk - puk + pupn))) + 
     -             densp*
     -              (2*(puk + pupe) - 
     -                ml2*
     -                 (-1 - 
     -                   denspc*mw2dm1*pek*
     -                    (pdpe + pdpn + 2*pek + 2*pnk - 
     -                     2*pupd + pupe + pupn) + 
     -                   mw2m1*
     -                    (pdpn + denspc*pdpn*pek + pnk - 
     -                     puk + pupn + 
     -                     pek*
     -                     (1 + 
     -                     denspc*
     -                     (pdpe + 2*pek + 2*pnk - 2*pupd + 
     -                     pupe + pupn)))))) + 
     -          denspc*puk*
     -           (-2*densp*pek*(pdk + pdpn - pnk) + 
     -             ml2*mw2dm1*
     -              (pek - pdpe*(1 + densp*pek) + pnk + 
     -                densp*pek*(-pdpn + pek + pnk)) + 
     -             densp*ml2*(-mw2dm1 + mw2m1)*pek*pupd)*qd + 
     -          denspc*pdk*
     -           (ml2*mw2dm1*(puk + pupe) + 
     -             densp*pek*
     -              (2*(pek - pupe) + 
     -                ml2*
     -                 (mw2m1*(pdk - pupd) + 
     -                   mw2dm1*
     -                    (-pdk - pek - pnk + pupd + pupe + 
     -                     pupn))))*qu) + 
     -       densp*(densc*puk*
     -           (2*pek*(1 + denspc*(pdk + pdpn - pnk)) + 
     -             ml2*mw2m1*
     -              (pdpe + denspc*pdpe*pek - pnk - 
     -                pek*(1 + denspc*(-pdpn + pek + pnk))) + 
     -             denspc*ml2*(-mw2dm1 + mw2m1)*pek*pupd)*qd
     -           - (densc*pdk*
     -              (ml2*mw2m1*(puk + pupe) + 
     -                denspc*pek*
     -                 (2*(pek - pupe) + 
     -                   ml2*
     -                    (mw2dm1*(pdk - pupd) + 
     -                     mw2m1*
     -                     (-pdk - pek - pnk + pupd + pupe + 
     -                     pupn)))) + 
     -             2*denspc*ml2*(mw2dm1 - mw2m1)*pek*
     -              (pdk - pupd)*qd)*qu)))/(pdk*pek*puk) + 
     -  (8*densp*(densc*(-(denspc*pdpn**2*pek*puk*pupe*qd) + 
     -          pdpn*puk*
     -           (-(pdk*pupe) + 
     -             pek*((1 + 2*denspc*pek - denspc*pnpe)*
     -                 pupd + pupe + 
     -                denspc*(-pdk + pek + pnk - 2*pupd)*pupe)
     -             )*qd - 
     -          puk*(pdk*pnpe*pupe + 
     -             denspc*pek*
     -              (-(pnk*pupd*(pnpe + 2*pupe)) + 
     -                pdk*(-2*pnk*pupe + pnpe*(pupd + pupe))))
     -            *qd - pdpe*puk*
     -           (pdpn*puk + 2*denspc*pdpn*pek*puk + 
     -             2*pdpn*pupe - pnk*pupe + 
     -             denspc*pek*pnk*pupe - 
     -             denspc*pek*(pdk + pdpn - pnk)*pupn)*qd + 
     -          pdk*pdpn*
     -           (2*pupe*(puk + pupe) + 
     -             pek*(puk*(-2 + denspc*(pnpe + pupe)) + 
     -                pupe*
     -                 (-2 + denspc*(-pnk + 2*pupd + pupe)))
     -              - denspc*pek**2*
     -              (2*puk + 2*pupd + pupe - pupn))*qu - 
     -          denspc*pdk*pdpe*pek*(-pek + puk + pupe)*pupn*
     -           qu + denspc*pdk*pek*
     -           (-(pek*pnpe*pupd) - 2*pnk*pupd*pupe + 
     -             pnpe*pupd*(puk + pupe) + 2*pdk*pupe*pupn)*
     -           qu) + 2*denspc*pek*
     -        (qu*((pdpe*pdpn*puk - 
     -                pdpn*
     -                 (pek*pupd + (pdk + puk - 2*pupd)*pupe)
     -                 + pdk*pupe*pupn)*qd + pdk*pdpn*pek*qu)
     -           + pnk*pupe*qd*(puk*qd - pupd*qu))))/
     -   (pdk*pek*puk) + 
     -  dens*((8*densc*(pdpn*
     -           (puk*(2 + 
     -                denspc*(pdpe + 2*pek - 2*pnpe - pupe))
     -              + denspc*
     -              (-(pupe*(pdk + 2*pnpe + pupe)) + 
     -                pek*(pupd + pupe))) + 
     -          denspc*(2*pdpe*pnk*pupe - 
     -             pnpe*(2*pdk*pupe + 
     -                pupd*(-pek + puk + pupe)) + 
     -             pdpe*(-pek + puk + pupe)*pupn)))/pek + 
     -     (-8*denspc*puk*
     -         (pdpn*(pdpe*puk - pek*pupd) + 
     -           (2*pdpe*pdpn - pdpn*pek - pdpe*pnk + 
     -              pdk*(pdpn + pnpe))*pupe)*qd - 
     -        16*denspc*pdk*pdpn*(pek - pupe)*(puk + pupe)*qu)
     -       /(pdk*pek*puk) + 
     -     densp*(8*densc*
     -         (pnpe*pupd - (pdk*(pdpn + 2*pnpe)*pupe)/pek + 
     -           pdpn*(pupd + pupe) - pdpe*pupn + 
     -           2*denspc*
     -            (-(pdk*pnpe*pupd) + pdpe*pnpe*pupd + 
     -              pek*pnpe*pupd + pnk*pnpe*pupd - 
     -              pnpe**2*pupd - pnpe*puk*pupd - 
     -              pnpe*pupd**2 + pdk*pnk*pupe - 
     -              2*pdk*pnpe*pupe + pdpe*pnpe*pupe + 
     -              2*pnk*pupd*pupe - 2*pnpe*pupd*pupe + 
     -              pdpn*
     -               (pdpe*pupd + (-pdk + pdpe + pnk)*pupe + 
     -                 pek*(puk + 2*pupd + pupe) - 
     -                 pupe*(puk + pupd - pupn) - 
     -                 pnpe*(2*puk + 2*pupd + pupe - pupn)) + 
     -              (pdk*pdpe - pdpe**2 + 
     -                 pdpe*
     -                  (-pek - pnk + pnpe + puk + pupd) + 
     -                 pupd*(pnpe + pupe))*pupn - pdpe*pupn**2
     -              ) + (puk*
     -              (-(pnpe*pupd) - 
     -                pdpn*(-2*pek + 2*pnpe + pupe) + 
     -                pdpe*(pdpn + pupn)))/pek - 
     -           (pupe*(pnpe*pupd + pdpn*(2*pnpe + pupe) - 
     -                pdpe*(2*pnk + pupn)))/pek) + 
     -        (8*denspc*(-(puk*
     -                ((-2*pdpn*pek + 
     -                     (pdk + pdpn - pnk)*pnpe)*pupd + 
     -                  (pdpn**2 + 
     -                     pdk*(pdpn - 2*pnk + pnpe) - 
     -                     pdpn*(pek + pnk - 2*pupd) - 
     -                     2*pnk*pupd)*pupe + 
     -                  pdpe*
     -                   (2*pdpn*puk + pnk*pupe - 
     -                     (pdk + pdpn - pnk)*pupn))*qd) + 
     -             pdk*(pupd*
     -                 (-2*pnk*pupe + pnpe*(puk + pupe)) + 
     -                pdpn*
     -                 (pnpe*puk + 
     -                   pupe*(-pnk + puk + 2*pupd + pupe) - 
     -                   pek*(2*puk + 2*pupd + pupe - pupn))
     -                 - (-2*pdk*pupe + pdpe*(puk + pupe))*
     -                 pupn + pek*(-(pnpe*pupd) + pdpe*pupn))*
     -              qu))/(pdk*puk))) - 
     -  (2*epupdpnk*(dens*
     -        (2*densc*ml2*pdk*puk*
     -           (denspc*(-1 + mw2dm1*(ml2 - pdpn + pnpe)) + 
     -             densp*
     -              (1 - 
     -                ml2*
     -                 (mw2m1 + denspc*(-mw2dm1 + mw2m1)*pek)
     -                 + denspc*mw2dm1*pek*
     -                 (2*pdpe + pek + pnpe + 2*pupe) + 
     -                mw2m1*
     -                 (pdpn - pnpe - 
     -                   denspc*pek*
     -                    (2*pdpe + pek + pnpe + 2*pupe)))) + 
     -          denspc*ml2*
     -           (ml2*mw2dm1 - 
     -             2*mw2dm1*
     -              (pdpn + 2*densp*pdpn*pek - pnpe) - 
     -             4*densp*(pdpe + pdpn)*pek*
     -              (mw2m1 - mw2m2*pnpe))*puk*qd + 
     -          denspc*pdk*
     -           (ml2**2*mw2dm1 - 4*densp*pek*pupe + 
     -             2*densp*ml2*pek*
     -              (-(mw2dm1*(pdpe + pek - 2*pupn)) - 
     -                2*mw2m2*pnpe*(pupe + pupn) + 
     -                mw2m1*(pdpe + pek + 2*(pupe + pupn))))*
     -           qu) - densp*
     -        (2*denspc*ml2*(mw2dm1 - mw2m1)*pek*
     -           (pdpe - pupe)*qd*qu + 
     -          densc*(ml2*
     -              (-4 + ml2*mw2m1 - 
     -                2*mw2m1*
     -                 (pdpn + 2*denspc*pdpn*pek - pnpe) - 
     -                4*denspc*(pdpe + pdpn)*pek*
     -                 (mw2dm1 - mw2m2*pnpe))*puk*qd + 
     -             pdk*(ml2**2*mw2m1 - 4*denspc*pek*pupe + 
     -                2*denspc*ml2*pek*
     -                 (-(mw2m1*(pdpe + pek - 2*pupn)) - 
     -                   2*mw2m2*pnpe*(pupe + pupn) + 
     -                   mw2dm1*(pdpe + pek + 2*(pupe + pupn))
     -                   ))*qu))))/(pdk*pek*puk) + 
     -  (2*ml2**2*(dens*(2*densc*pdk*puk*
     -           (densp*mw2m1*
     -              (pdpn*puk - pnk*pupd + pdk*pupn) + 
     -             denspc*
     -              (mw2dm1*
     -                 (-((1 + densp*pek)*pnk*pupd) + 
     -                   pdpn*
     -                    (puk + densp*pek*puk + 
     -                     densp*pek*pupd) + pdk*pupn + 
     -                   densp*pek*(pdk + pupd)*pupn) + 
     -                densp*pek*
     -                 (2*mw2m2*pnpe*(-(pdk*puk) + pupd**2) + 
     -                   mw2m1*
     -                    (-(pnk*pupd) + pdpn*(puk + pupd) + 
     -                     (pdk + pupd)*pupn)))) + 
     -          denspc*puk*
     -           (-2*densp*mw2m2*pek*pnpe*
     -              (pdk*puk + (puk - pupd)*pupd) + 
     -             mw2dm1*
     -              (pdpn*(puk + 2*densp*pek*pupd) - 
     -                (1 + 2*densp*pek)*(pnk*pupd - pdk*pupn))
     -             )*qd - 
     -          denspc*pdk*
     -           (mw2dm1*(1 + 2*densp*pek)*
     -              (pdpn*puk - pnk*pupd) - 
     -             2*densp*mw2m2*pek*pnpe*
     -              (-pupd**2 + pdk*(puk + pupd)) + 
     -             mw2dm1*(pdk + 2*densp*pek*pupd)*pupn)*qu)
     -        + densp*(densc*puk*
     -           (-2*denspc*mw2m2*pek*pnpe*
     -              (pdk*puk + (puk - pupd)*pupd) + 
     -             mw2m1*
     -              (pdpn*(puk + 2*denspc*pek*pupd) - 
     -                (1 + 2*denspc*pek)*(pnk*pupd - pdk*pupn)
     -                ))*qd - 
     -          densc*pdk*
     -           (mw2m1*(1 + 2*denspc*pek)*
     -              (pdpn*puk - pnk*pupd) - 
     -             2*denspc*mw2m2*pek*pnpe*
     -              (-pupd**2 + pdk*(puk + pupd)) + 
     -             mw2m1*(pdk + 2*denspc*pek*pupd)*pupn)*qu - 
     -          2*denspc*mw2m2*pek*pnpe*
     -           (puk**2*qd**2 - 
     -             2*(pdk + puk - pupd)*pupd*qd*qu + 
     -             pdk**2*qu**2))))/(pdk*pek*puk) + 
     -  ml2*(densp*((4*denspc*
     -           (-(puk*(mw2dm1*
     -                   (-(pnpe*puk) + pnk*pupe + pek*pupn + 
     -                     2*pnk*pupn) + 
     -                  mw2m1*
     -                   (-(pnpe*puk) + pnk*pupe + pek*pupn + 
     -                     2*pnk*pupn) + 
     -                  2*mw2m2*pnpe*
     -                   (pnpe*puk - 
     -                     (pek + pnk)*(pupe + pupn)))*qd**2)
     -              + (mw2dm1*
     -                 (-2*pdpn**2*puk - 2*pdk*pnpe*pupd - 
     -                   2*pnpe*puk*pupd + 2*pnpe*pupd**2 + 
     -                   pnk*pupd*pupe + pek*pupd*pupn + 
     -                   2*pnk*pupd*pupn - 2*pdk*pupe*pupn - 
     -                   2*pdk*pupn**2 + 
     -                   pdpe*
     -                    (-2*pdpn*puk + pnk*pupd + 
     -                     (pdk + puk - 2*pupd)*pupn) + 
     -                   pdpn*
     -                    (pek*pupd + 2*pnk*pupd + 
     -                     (pdk + puk - 2*pupd)*
     -                     (pupe + 2*pupn))) + 
     -                mw2m1*
     -                 (-2*pdpn**2*puk - 2*pdk*pnpe*pupd - 
     -                   2*pnpe*puk*pupd + 2*pnpe*pupd**2 + 
     -                   pnk*pupd*pupe + pek*pupd*pupn + 
     -                   2*pnk*pupd*pupn - 2*pdk*pupe*pupn - 
     -                   2*pdk*pupn**2 + 
     -                   pdpe*
     -                    (-2*pdpn*puk + pnk*pupd + 
     -                     (pdk + puk - 2*pupd)*pupn) + 
     -                   pdpn*
     -                    (pek*pupd + 2*pnk*pupd + 
     -                     (pdk + puk - 2*pupd)*
     -                     (pupe + 2*pupn))) - 
     -                2*mw2m2*pnpe*
     -                 (-(pdpe**2*puk) - pdpn**2*puk - 
     -                   2*pdk*pnpe*pupd - 2*pnpe*puk*pupd + 
     -                   2*pnpe*pupd**2 + pek*pupd*pupe + 
     -                   pnk*pupd*pupe - pdk*pupe**2 + 
     -                   pek*pupd*pupn + pnk*pupd*pupn - 
     -                   2*pdk*pupe*pupn - pdk*pupn**2 + 
     -                   pdpn*
     -                    (pek*pupd + pnk*pupd + 
     -                     (pdk + puk - 2*pupd)*(pupe + pupn))
     -                     + pdpe*
     -                    (-2*pdpn*puk + 
     -                     (pdk + puk)*(pupe + pupn) + 
     -                     pupd*(pek + pnk - 2*(pupe + pupn)))
     -                   ))*qd*qu - 
     -             pdk*((mw2dm1 + mw2m1)*
     -                 (pdpe*pnk + pdpn*(pek + 2*pnk)) - 
     -                ((mw2dm1 + mw2m1)*pdk + 
     -                   2*mw2m2*(pdpe + pdpn)*(pek + pnk))*
     -                 pnpe + 2*mw2m2*pdk*pnpe**2)*qu**2))/
     -         (pdk*puk) + 
     -        densc*((8*pdpn*qu)/pek + 
     -           (4*mw2m1*
     -              (puk*
     -                 (pdpn**2*puk + pdpe**2*pupn + 
     -                   pnpe*(pek*pupd + pdk*pupn) - 
     -                   pdpn*
     -                    (pnk*pupd - pdk*(pupe + pupn) + 
     -                     pek*(pupd + pupe + pupn)) + 
     -                   pdpe*
     -                    (-(pnpe*pupd) - (pek + pnk)*pupn + 
     -                     pdpn*(puk + pupe + 2*pupn)))*qd + 
     -                pdk*
     -                 (pdpe*pnk*puk + 
     -                   pnpe*
     -                    (-(pdk*puk) + 
     -                     pupd*(-pek + puk + pupe)) - 
     -                   pdpe*(-pek + puk + pupe)*pupn + 
     -                   pdpn*
     -                    (2*pnk*puk - pnpe*puk + pnk*pupe - 
     -                     puk*pupe - pupe**2 - 
     -                     2*(puk + pupe)*pupn + 
     -                     pek*(puk + pupe + pupn)))*qu))/
     -            (pdk*pek*puk) + 
     -           (4*denspc*
     -              (mw2m1*
     -                 (pdpn**2*puk*(2*puk - pupe)*qd + 
     -                   pdpe**2*puk*pupn*qd + 
     -                   pdpe*puk*
     -                    (2*pdpn*puk - pnpe*pupd - 
     -                     pdpn*pupe + pnk*pupe + 
     -                     (pdpn - pek + pupd)*pupn)*qd + 
     -                   pdpn*puk*
     -                    (pnk*(-2*pupd + pupe) - 
     -                     pek*(2*pupd + pupn) + 
     -                     pupd*(pnpe + pupe + 2*pupn))*qd + 
     -                   puk*
     -                    (pdk*
     -                     (-(pnk*pupe) - pek*pupn - 
     -                     2*pnk*pupn + 
     -                     pnpe*(puk - pupe + pupn)) - 
     -                     pupd*
     -                     (pnpe*(-puk + pupd) + 
     -                     pek*(-pnpe + pupn) + 
     -                     pnk*(pnpe + pupe + 2*pupn)))*qd + 
     -                   pdk*pdpe*
     -                    (pnpe*puk + 
     -                     pnk*(puk + pupd - pupe - pupn) + 
     -                     pupn*(-pupd + pupe + pupn))*qu + 
     -                   pdk*pdpn*
     -                    (-(pnpe*puk) + 2*pnk*(puk + pupd) - 
     -                     pupe*(pupd + pupe) - 
     -                     (2*pupd + pupe)*pupn + 
     -                     pek*(puk + pupd + pupe + pupn))*qu
     -                    - pdk*
     -                    (-(pnpe*pupd*
     -                     (-pek + pnk + pupd + pupe)) + 
     -                     (-2*(pek + pnk) + pnpe)*pupd*
     -                    pupn + 
     -                     pdk*
     -                     (pnpe*(puk + pupd) + 
     -                     2*pupn*(pupe + pupn)))*qu) + 
     -                mw2dm1*
     -                 (puk*
     -                    (2*pdpn**2*puk + pdk*pnpe*puk + 
     -                     pnpe*puk*pupd - pnpe*pupd**2 - 
     -                     pdk*pnk*pupe - pnk*pupd*pupe - 
     -                     (pek + 2*pnk)*(pdk + pupd)*pupn + 
     -                     pdpn*pupd*
     -                     (-2*pnk + pupe + 2*pupn) + 
     -                     pdpe*
     -                     (2*pdpn*puk + pupd*(-2*pnk + pupn))
     -                     )*qd + 
     -                   pdk*
     -                    (-(pdk*pnpe*puk) - pdk*pnpe*pupd + 
     -                     pnpe*pupd**2 + 
     -                     pdpe*pnk*(puk + pupd) + 
     -                     2*pnk*pupd*pupe - pdpe*pupd*pupn + 
     -                     2*pnk*pupd*pupn - 
     -                     2*pdk*pupe*pupn - 2*pdk*pupn**2 + 
     -                     pdpn*
     -                     (pek*(puk + pupd) + 
     -                     2*pnk*(puk + pupd) - 
     -                     pupd*(pupe + 2*pupn)))*qu) + 
     -                2*(puk*
     -                    (-(pdpn*pupd*qd) + pnk*pupd*qd + 
     -                     pdk*pdpn*qu) + 
     -                   mw2m2*pnpe*
     -                    (-(pdpe**2*puk**2*qd) + 
     -                     pdpe*puk*
     -                     (-2*pdpn*puk + 
     -                     pupd*(pek + pnk - pupe - pupn))*qd
     -                     + puk*
     -                     (-(pdpn**2*puk) + 
     -                     pdpn*pupd*
     -                     (pek + pnk - pupe - pupn) + 
     -                     pdk*
     -                     (-(pnpe*puk) + 
     -                     (pek + pnk)*(pupe + pupn)) + 
     -                     pupd*
     -                     (pnpe*(-puk + pupd) + 
     -                     (pek + pnk)*(pupe + pupn)))*qd - 
     -                     pdk*pdpe*
     -                     (pek*(puk + pupd) + 
     -                     pnk*(puk + pupd) - 
     -                     pupd*(pupe + pupn))*qu + 
     -                     pdk*
     -                     (-(pupd*
     -                     (pnpe*pupd + 
     -                     (pek + pnk)*(pupe + pupn))) - 
     -                     pdpn*
     -                     (pek*(puk + pupd) + 
     -                     pnk*(puk + pupd) - 
     -                     pupd*(pupe + pupn)) + 
     -                     pdk*
     -                     (pnpe*(puk + pupd) + 
     -                     (pupe + pupn)**2))*qu))))/(pdk*puk)
     -           )) + dens*
     -      ((-4*densc*(denspc*mw2dm1*pdpn**2*pek*
     -              (puk + pupe) + 
     -             pdpn*(4*pupe + 
     -                puk*
     -                 (4 + 
     -                   denspc*pek*
     -                    (5 + 
     -                     mw2dm1*
     -                     (pdpe + pek + 2*pnk - pnpe + pupe))
     -                   ) + 
     -                denspc*pek*
     -                 ((2 + mw2dm1*(pek + pnk - pnpe))*
     -                    pupd + 
     -                   (3 + mw2dm1*(-pdk + pek + pnk))*
     -                    pupe - 
     -                   (2 + mw2dm1*(pdk + pdpe + pupe))*pupn
     -                   )) + 
     -             denspc*pek*
     -              ((pnk + pnpe)*pupd - (pdk + pdpe)*pupn + 
     -                mw2dm1*
     -                 (-(pnpe*
     -                     (pupd*(pek - puk + pupn) + 
     -                     pdk*(puk + pupn))) + 
     -                   pdpe*
     -                    (pupn*(pek - puk + pupn) + 
     -                     pnk*(puk + pupn))))))/pek**2 + 
     -        (4*denspc*(2*pdk*pdpn*puk*qu + 
     -             mw2dm1*
     -              (puk*
     -                 (pdpn**2*puk + pdpe**2*pupn + 
     -                   pnpe*(pek*pupd + pdk*pupn) - 
     -                   pdpn*
     -                    (pnk*pupd - pdk*(pupe + pupn) + 
     -                     pek*(pupd + pupe + pupn)) + 
     -                   pdpe*
     -                    (-(pnpe*pupd) - (pek + pnk)*pupn + 
     -                     pdpn*(puk + pupe + 2*pupn)))*qd + 
     -                pdk*
     -                 (pdpe*pnk*puk + 
     -                   pnpe*
     -                    (-(pdk*puk) + 
     -                     pupd*(-pek + puk + pupe)) - 
     -                   pdpe*(-pek + puk + pupe)*pupn + 
     -                   pdpn*
     -                    (2*pnk*puk - pnpe*puk + pnk*pupe - 
     -                     puk*pupe - pupe**2 - 
     -                     2*(puk + pupe)*pupn + 
     -                     pek*(puk + pupe + pupn)))*qu)))/
     -         (pdk*pek*puk) + 
     -        densp*(densc*
     -            ((-20*pdpn*puk)/pek + (4*pdk*pupn)/pek + 
     -              (-4*(2*pdpn + pnk + pnpe)*pupd - 
     -                 12*pdpn*pupe + 4*(pdpe + 2*pdpn)*pupn)/
     -               pek + 
     -              (4*mw2m1*
     -                 (-(pdpn**2*(puk + pupe)) - 
     -                   pdpn*
     -                    (-(pnpe*(puk + pupd)) + 
     -                     (-pdk + puk)*pupe + 
     -                     pek*(puk + pupd + pupe) + 
     -                     pnk*(2*puk + pupd + pupe)) + 
     -                   pdpn*(pdk + pupe)*pupn - 
     -                   pdpe*
     -                    ((pdpn + pnk)*puk + 
     -                     (-pdpn + pek + pnk - puk)*pupn + 
     -                     pupn**2) + 
     -                   pnpe*
     -                    (pupd*(pek - puk + pupn) + 
     -                     pdk*(puk + pupn))))/pek + 
     -              4*denspc*
     -               (-6*pdpn*puk + 2*pdk*pupn + 
     -                 2*
     -                  (-(pdpn*(4*pupd + pupe - 2*pupn)) + 
     -                    pdpe*pupn + 
     -                    pupd*(pnk - 2*pnpe + 2*pupn)) + 
     -                 2*mw2m2*pnpe*
     -                  (((pdpe + pdpn)*(pek + pnk) - 
     -                     2*pdk*pnpe)*puk + 
     -                    pdk*(pek + pnk)*(pupe + pupn) + 
     -                    2*pek*pupd*
     -                     (pdpe + pdpn + pupe + pupn) + 
     -                    pupd*
     -                     (pdpe**2 + pdpn**2 + 2*pdpn*pnk + 
     -                     2*pdpe*(pdpn + pnk) + 
     -                     2*pnpe*pupd + 
     -                     (pupe + pupn)*(2*pnk + pupe + pupn)
     -                     )) + 
     -                 mw2dm1*
     -                  (-((pdpn*(pek + 2*pnk - pnpe) - 
     -                     2*pdk*pnpe + pdpe*(pnk + pnpe))*puk
     -                     ) - pdpn**2*(2*pupd + pupe) + 
     -                    pek*pupd*
     -                    (-2*pdpn + pnpe - 2*pupn) + 
     -                    pdpe**2*pupn - 
     -                    pdk*
     -                     ((pnk + pnpe)*pupe + 
     -                     (pek + 2*pnk - pnpe)*pupn) + 
     -                    pdpn*
     -                     (-4*pnk*pupd + pnpe*pupd + 
     -                     pupe*(pupe + pupn)) - 
     -                    pdpe*
     -                     (2*pnk*pupd + pnpe*pupd + 
     -                     pdpn*(2*pupd + pupe - pupn) + 
     -                     pupn*(pupe + pupn)) - 
     -                    pupd*
     -                     (pnpe*(2*pupd + pupe - pupn) + 
     -                     2*pupn*(pupe + pupn) + 
     -                     pnk*(pnpe + 2*pupe + 4*pupn))) + 
     -                 mw2m1*
     -                  (-((pdpn*(pek + 2*pnk - pnpe) - 
     -                     2*pdk*pnpe + pdpe*(pnk + pnpe))*puk
     -                     ) - pdpn**2*(2*pupd + pupe) + 
     -                    pek*pupd*
     -                    (-2*pdpn + pnpe - 2*pupn) + 
     -                    pdpe**2*pupn - 
     -                    pdk*
     -                     ((pnk + pnpe)*pupe + 
     -                     (pek + 2*pnk - pnpe)*pupn) + 
     -                    pdpn*
     -                     (-4*pnk*pupd + pnpe*pupd + 
     -                     pupe*(pupe + pupn)) - 
     -                    pdpe*
     -                     (2*pnk*pupd + pnpe*pupd + 
     -                     pdpn*(2*pupd + pupe - pupn) + 
     -                     pupn*(pupe + pupn)) - 
     -                    pupd*
     -                     (pnpe*(2*pupd + pupe - pupn) + 
     -                     2*pupn*(pupe + pupn) + 
     -                     pnk*(pnpe + 2*pupe + 4*pupn))))) + 
     -           (4*denspc*
     -              (mw2dm1*
     -                 (pdpn**2*puk*(2*puk - pupe)*qd + 
     -                   pdpe**2*puk*pupn*qd + 
     -                   pdpe*puk*
     -                    (2*pdpn*puk - pnpe*pupd - 
     -                     pdpn*pupe + pnk*pupe + 
     -                     (pdpn - pek + pupd)*pupn)*qd + 
     -                   pdpn*puk*
     -                    (pnk*(-2*pupd + pupe) - 
     -                     pek*(2*pupd + pupn) + 
     -                     pupd*(pnpe + pupe + 2*pupn))*qd + 
     -                   puk*
     -                    (pdk*
     -                     (-(pnk*pupe) - pek*pupn - 
     -                     2*pnk*pupn + 
     -                     pnpe*(puk - pupe + pupn)) - 
     -                     pupd*
     -                     (pnpe*(-puk + pupd) + 
     -                     pek*(-pnpe + pupn) + 
     -                     pnk*(pnpe + pupe + 2*pupn)))*qd + 
     -                   pdk*pdpe*
     -                    (pnpe*puk + 
     -                     pnk*(puk + pupd - pupe - pupn) + 
     -                     pupn*(-pupd + pupe + pupn))*qu + 
     -                   pdk*pdpn*
     -                    (-(pnpe*puk) + 2*pnk*(puk + pupd) - 
     -                     pupe*(pupd + pupe) - 
     -                     (2*pupd + pupe)*pupn + 
     -                     pek*(puk + pupd + pupe + pupn))*qu
     -                    - pdk*
     -                    (-(pnpe*pupd*
     -                     (-pek + pnk + pupd + pupe)) + 
     -                     (-2*(pek + pnk) + pnpe)*pupd*
     -                    pupn + 
     -                     pdk*
     -                     (pnpe*(puk + pupd) + 
     -                     2*pupn*(pupe + pupn)))*qu) + 
     -                mw2m1*
     -                 (puk*
     -                    (2*pdpn**2*puk + pdk*pnpe*puk + 
     -                     pnpe*puk*pupd - pnpe*pupd**2 - 
     -                     pdk*pnk*pupe - pnk*pupd*pupe - 
     -                     (pek + 2*pnk)*(pdk + pupd)*pupn + 
     -                     pdpn*pupd*
     -                     (-2*pnk + pupe + 2*pupn) + 
     -                     pdpe*
     -                     (2*pdpn*puk + pupd*(-2*pnk + pupn))
     -                     )*qd + 
     -                   pdk*
     -                    (-(pdk*pnpe*puk) - pdk*pnpe*pupd + 
     -                     pnpe*pupd**2 + 
     -                     pdpe*pnk*(puk + pupd) + 
     -                     2*pnk*pupd*pupe - pdpe*pupd*pupn + 
     -                     2*pnk*pupd*pupn - 
     -                     2*pdk*pupe*pupn - 2*pdk*pupn**2 + 
     -                     pdpn*
     -                     (pek*(puk + pupd) + 
     -                     2*pnk*(puk + pupd) - 
     -                     pupd*(pupe + 2*pupn)))*qu) + 
     -                2*(puk*
     -                    (-(pdpn*pupd*qd) + pnk*pupd*qd + 
     -                     pdk*pdpn*qu) + 
     -                   mw2m2*pnpe*
     -                    (-(pdpe**2*puk**2*qd) + 
     -                     pdpe*puk*
     -                     (-2*pdpn*puk + 
     -                     pupd*(pek + pnk - pupe - pupn))*qd
     -                     + puk*
     -                     (-(pdpn**2*puk) + 
     -                     pdpn*pupd*
     -                     (pek + pnk - pupe - pupn) + 
     -                     pdk*
     -                     (-(pnpe*puk) + 
     -                     (pek + pnk)*(pupe + pupn)) + 
     -                     pupd*
     -                     (pnpe*(-puk + pupd) + 
     -                     (pek + pnk)*(pupe + pupn)))*qd - 
     -                     pdk*pdpe*
     -                     (pek*(puk + pupd) + 
     -                     pnk*(puk + pupd) - 
     -                     pupd*(pupe + pupn))*qu + 
     -                     pdk*
     -                     (-(pupd*
     -                     (pnpe*pupd + 
     -                     (pek + pnk)*(pupe + pupn))) - 
     -                     pdpn*
     -                     (pek*(puk + pupd) + 
     -                     pnk*(puk + pupd) - 
     -                     pupd*(pupe + pupn)) + 
     -                     pdk*
     -                     (pnpe*(puk + pupd) + 
     -                     (pupe + pupn)**2))*qu))))/(pdk*puk)
     -           )))



      amp2 = amp2*
     +       g2*conjg(g2)*el2/4d0
     +       /4d0/3
     +       /(st_alpha/(2*pi))

      if(mod(abs(i),2).eq.0) then
         amp2=amp2*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      elseif(mod(abs(i),2).eq.1) then   
         amp2=amp2*ph_CKM(abs(j)/2,(abs(i)+1)/2)**2
      endif

      return
      end 
*
**
*
      complex*16 function e_(q1,q2,q3,q4)
      implicit none
      real*8 q1(0:3),q2(0:3),q3(0:3),q4(0:3)

      complex*16 esign

      esign = -(0.d0,1.d0)

      e_ = q1(1)*q2(2)*q3(3)*q4(0) - q1(1)*q2(2)*q3(0)*q4(3) +
     .     q1(1)*q2(3)*q3(0)*q4(2) - q1(1)*q2(3)*q3(2)*q4(0) +
     .     q1(1)*q2(0)*q3(2)*q4(3) - q1(1)*q2(0)*q3(3)*q4(2) +
     .     q1(2)*q2(1)*q3(0)*q4(3) - q1(2)*q2(1)*q3(3)*q4(0) +
     .     q1(2)*q2(3)*q3(1)*q4(0) - q1(2)*q2(3)*q3(0)*q4(1) +
     .     q1(2)*q2(0)*q3(3)*q4(1) - q1(2)*q2(0)*q3(1)*q4(3) +
     .     q1(3)*q2(1)*q3(2)*q4(0) - q1(3)*q2(1)*q3(0)*q4(2) +
     .     q1(3)*q2(2)*q3(0)*q4(1) - q1(3)*q2(2)*q3(1)*q4(0) +
     .     q1(3)*q2(0)*q3(1)*q4(2) - q1(3)*q2(0)*q3(2)*q4(1) +
     .     q1(0)*q2(1)*q3(3)*q4(2) - q1(0)*q2(1)*q3(2)*q4(3) +
     .     q1(0)*q2(2)*q3(1)*q4(3) - q1(0)*q2(2)*q3(3)*q4(1) +
     .     q1(0)*q2(3)*q3(2)*q4(1) - q1(0)*q2(3)*q3(1)*q4(2)

      e_ = esign * e_ 

      return
      end function e_
*
**
*
      subroutine invertspace(p,pp)
      implicit none
      real*8 p(0:3),pp(0:3)
      integer nu

      pp(0)=p(0)
      do nu=1,3
          pp(nu)=-p(nu)
      enddo

      end
