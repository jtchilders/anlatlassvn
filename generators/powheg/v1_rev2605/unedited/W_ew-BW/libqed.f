c*******************************************************************
      subroutine mathard_cc(i,matri)
c matri: MES and averaged (summed) over IS (FS) spin and colour DOF 
c of hard photon emission in the W+ production process at the Tevatron and LHC:
c i (=quark with charge qi, mass mi and four momentum pi,1) +
c i'(=antiquark,qis,mis,pis,2) --> W+  -->
c f (=neutrino,qf,mf,pf,4) + f'(=lepton+,qfs,mfs,pfs,3) + gamma(momentum k)
c i=1,2,3: initial state radiation, final state radiation and interference.
c j=1,2: U(x1)D~(x2), U(x2)D~(x1)
      implicit real*8(a-z)
      integer i,ii,j,jj
      real*8 matri(2),dotp
      real*8 sinv(5,5)
      real*8 pi,pi2
      external dotp
      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
c     choice of representation:
      integer i1(2),i2(2)
      data i1,i2/1,2,2,1/
      pi=4d0*datan(1d0)
      pi2 = pi**2
c     alpha representation:
      if (rep.eq.1) then
      sig0s = 64d0*(pi*alpha/sw2)**2*pi2/3d0
c     G_mu representation:
      else
      sig0s = 64d0*(w2*gfermi*xmw)**2*pi2/3d0
      endif

c itwasme - BEGIN we need rl_sinv here
      do ii = 1,5
        do jj = 1,5
        sinv(ii,jj) = rl_sinv(ii,jj)
        enddo
      enddo
c itwasme END
      do j=1,2
         shat = sinv(1,2)
         kq = (sinv(1,5)+sinv(2,5))/2d0
         kpi = sinv(i1(j),5)/2d0
         kpis = sinv(i2(j),5)/2d0
         kpf = sinv(4,5)/2d0
         kpfs = sinv(3,5)/2d0
         pfpis = sinv(i2(j),4)/2d0
         pfspi = sinv(i1(j),3)/2d0
         pfpi = sinv(i1(j),4)/2d0
         pfspis = sinv(i2(j),3)/2d0
      
         dw = shat-xmw
         mat = 0d0

c        initial state:
         if (i.eq.1) then 
c         mat = alpha0/pi/((dw-2d0*kq)**2+xmw*gamw**2)*
c     $        (-(qi*pfpis*mmi/kpi)**2-(qis*pfspi*mmis/kpis)**2+
c     $        (pfspi**2+pfpis**2)*(-qi**2/kpi-qis**2/kpis+
c     $        (qi-qis)**2/kq+shat/2d0*(qi*qis/kpi/kpis+
c     $        qi*(qi-qis)/kpi/kq-
c     $        qis*(qi-qis)/kpis/kq-(qi-qis)**2/kq**2)))
            if (wopt.eq.1) then
            mat = alpha0/pi/((dw-2d0*kq)**2+xmw*gamw**2)*(
c     $            (-(qi*pfpis*mmi/kpi)**2-(qis*pfspi*mmis/kpis)**2+
     $            (pfspi**2+pfpis**2)*(qi**2/kpi+qis**2/kpis-
     $            (qi-qis)**2/kq)*(shat/2d0/kq-1d0))
            else if (wopt.eq.2) then
            mat = alpha0/pi/((dw-2d0*kq)**2+
     $            ((shat-2d0*kq)/mw*gamw)**2)*(
c     $            (-(qi*pfpis*mmi/kpi)**2-(qis*pfspi*mmis/kpis)**2+
     $            (pfspi**2+pfpis**2)*(qi**2/kpi+qis**2/kpis-
     $            (qi-qis)**2/kq)*(shat/2d0/kq-1d0))
            end if
         end if
c final state:
         if (i.eq.2) then 
            if (wopt.eq.1) then
            mat = alpha0/pi/(dw**2+xmw*gamw**2)*(
     $           -(qf*pfspi*mmf/kpf)**2-(qfs*pfpis*mmfs/kpfs)**2+
     $           (pfspi**2+pfpis**2)*shat/2d0/kq*(qf**2/kpf+
     $           qfs**2/kpfs-(qf-qfs)**2/kq))
            else if (wopt.eq.2) then
            mat = alpha0/pi/(dw**2+(shat/mw*gamw)**2)*(
     $           -(qf*pfspi*mmf/kpf)**2-(qfs*pfpis*mmfs/kpfs)**2+
     $           (pfspi**2+pfpis**2)*shat/2d0/kq*(qf**2/kpf+
     $           qfs**2/kpfs-(qf-qfs)**2/kq))
            end if
         end if

c interference:
         if (i.eq.3) then
            if (wopt.eq.1) then
            mat = alpha0/pi/2d0/(dw**2+xmw*gamw**2)*
     $              (dw*(dw-2d0*kq)+xmw*gamw**2)/
     $              ((dw-2d0*kq)**2+xmw*gamw**2)*(pfspi**2+pfpis**2)*
     $              ((pfpi+pfspis)*(qi*qf/kpi/kpf+qis*qfs/kpis/kpfs)-
     $              (pfspi+pfpis)*(qi*qfs/kpfs/kpi+qis*qf/kpf/kpis)-
     $              (shat/kq-1d0)*(qi*(qf-qfs)/kpi-qis*(qf-qfs)/kpis+
     $              qf*(qi-qis)/kpf-qfs*(qi-qis)/kpfs-
     $              2d0*(qi-qis)*(qf-qfs)/kq))
            else if (wopt.eq.2) then
            mat = alpha0/pi/2d0/(dw**2+(shat/mw*gamw)**2)*
     $              (dw*(dw-2d0*kq)+shat*(shat-2d0*kq)/xmw*gamw**2)/
     $              ((dw-2d0*kq)**2+((shat-2d0*kq)/mw*gamw)**2)*
     $              (pfspi**2+pfpis**2)*
     $              ((pfpi+pfspis)*(qi*qf/kpi/kpf+qis*qfs/kpis/kpfs)-
     $              (pfspi+pfpis)*(qi*qfs/kpfs/kpi+qis*qf/kpf/kpis)-
     $              (shat/kq-1d0)*(qi*(qf-qfs)/kpi-qis*(qf-qfs)/kpis+
     $              qf*(qi-qis)/kpf-qfs*(qi-qis)/kpfs-
     $              2d0*(qi-qis)*(qf-qfs)/kq))
            end if
        end if
      matri(j) = sig0s*mat
      end do
      return
      end
c********************************************************************
      subroutine mathard_qcd(matri)

c matri: matrix element squared and averaged (summed) over 
c initial state (final state) spin and colour degrees
c of freedom to hard gluon emission in the
c W+ production process at the Tevatron and LHC:
c i(=quark with charge qi, mass mi and four momentum pi,1) +
c i'(=antiquark,qis,mis,pis,2) -> W+ 
c -> f(=neutrino,qf,mf,pf,4) + f'(=lepton+,qfs,mfs,pfs,3) + gluon(momentum k)
c j =1,2: U(x1)D~(x2), U(x2)D~(x1)

      implicit real*8(a-z)
      integer i,ii,j,jj
      real*8 matri(2)
      include 'pwhg_wzgrad.h'
      real*8 sinv(5,5)
      real*8 pi,pi2

      integer i1(2),i2(2)
      data i1,i2/1,2,2,1/

      pi = 4d0*datan(1d0)
      pi = pi**2
c alpha representation:
      if (rep.eq.1) then
         sig0s = 64d0*(pi*alpha/sw2)**2*pi2/3d0
c G_mu representation:
      else
         sig0s = 64d0*(w2*gfermi*xmw)**2*pi2/3d0
      endif
c itwasme - BEGIN we need rl_sinv here
c      do ii = 1,5
c        do jj = 1,5
c        sinv(ii,jj) = rl_sinv(ii,jj)
c        enddo
c      enddo
c itwasme END
      do j=1,2
         shat = sinv(1,2)
         pipis=shat/2d0
         kq = (sinv(1,5)+sinv(2,5))/2d0
         kpi = sinv(i1(j),5)/2d0
         kpis = sinv(i2(j),5)/2d0
         kpf = sinv(4,5)/2d0
         kpfs = sinv(3,5)/2d0
         pfpis = sinv(i2(j),4)/2d0
         pfspi = sinv(i1(j),3)/2d0
         pfpi = sinv(i1(j),4)/2d0
         pfspis = sinv(i2(j),3)/2d0
      
         dw = shat-xmw

         mat = 0d0

         cf=4d0/3d0
         ci=kpfs*pfpis+pfpi*pfspi-pfpis*pfspi
         cis=kpf*pfspi-pfpis*pfspi+pfpis*pfspis
         ciis=-((kpfs*pfpis+(kpf-2d0*pfpis)*pfspi)*pipis)

         if (wopt.eq.1) then
            mat = alphas/pi/((dw-2d0*kq)**2+xmw*gamw**2)*cf*
     $           (cis/kpis+ci/kpi+ciis/kpi/kpis)
         else if (wopt.eq.2) then
            mat = alphas/pi/((dw-2d0*kq)**2+
     $           ((shat-2d0*kq)/mw*gamw)**2)*cf*
     $           (cis/kpis+ci/kpi+ciis/kpi/kpis)
         end if

         matri(j) = sig0s*mat
         
      end do
      return
      end
c**********************************************************************
      subroutine softcollqed_cc(i,fqed)
c virtual+soft photon contribution, eps = upper limit of photon energy
c sig(0+1) = sig0(g0->g(0+1))*(1+fqed)
c final-state collinear contribution and PDF CT added
      implicit real*8(a-z)
      integer i,k,j,l
      real*8 fqed(2)
      real*8 sinv(4,4)
      complex*16 spence

      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
      real*8 pi
      integer i1(2)
      data i1/1,2/
      pi = 4d0*datan(1d0)
c     itwasme BEGIN
      do k = 1,4
      do l = 1,4
      sinv(k,l) = br_sinv(k,l)
      enddo
      enddo
c     itwasme END  
      do j=1,2!main do-loop
         shat = sinv(1,2)
         that = -sinv(i1(j),4)
         eps = deltas*dsqrt(shat)/2d0
         le = dlog(deltas)
         dw = shat-xmw
         if (dw.eq.0d0) then
            phase = 0d0
            lw = dlog(gamw**2/(gamw**2+4d0*eps**2))/2d0
            lws = dlog(xmw/(gamw**2+4d0*eps**2))/2d0
         else
            phase = dw/mw/gamw*(datan(dw/mw/gamw)+
     &           datan((2d0*dsqrt(shat)*eps-dw)/mw/gamw))
            lw = dlog((dw**2+xmw*gamw**2)/(xmw*gamw**2+
     &           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
            lws = dlog(xmw**2/(xmw*gamw**2+
     &           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
         endif
c-----i = 1------------------------------------------------------------
c        initial state:
         if (i.eq.1) then
c unsubtracted (with shift in w-propagator):
c            call beta(i,shat,that,beti)
c            call delta(i,shat,that,delti)
c            fqed(j) = beti*(le+lw+phase)+2d0*delti
c minimal factorisation: mi replaced by mu_f (and
c shift in propagator neglected):

c            fqed(j) = alpha0/pi*((qi**2+qis**2)*((le+3d0/4d0)*
c     $           dlog(shat/mu_f**2)+pi**2/6d0-1d0)-
c     $           le+3d0/2d0+pi**2/24d0)
c pdf counterterm added (and shift in propagator neglected):
c       what are all these variables in fqed(1) and where
c       are they defined in pwhg? where are they defined in 
c       wzgrad? are they defined the same way?
            fvps = 9d0+2d0/3d0*pi**2+3d0*le-2d0*le**2
        
            fqed(j) = alpha0/pi*((qi**2+qis**2)*((le+3d0/4d0)*
     $           dlog(shat/mu_f**2)+pi**2/6d0-2d0+le**2+lfc/4d0*fvps)-
     $           le+3d0/2d0+pi**2/24d0)
         end if
c-----i = 2------------------------------------------------------------
c        final state:
         if (i.eq.2) then
c the following expressions are only valid for a leptonic final state
            if (test(4).eq.1) then
               if (collcut.eq.0) then
               call beta(i,shat,that,betf)
               call delta(i,shat,that,deltf)
c               fqed(j) = betf*le+2d0*deltf
c with terms prop. to deltas:
               fqed(j) = betf*le
     &                 + 2d0*deltf
     &                 + betf*deltas*(deltas-4d0)/4d0
     &                 + alpha0/pi*(-dreal(spence(deltas+ieps))
     &                 + (1d0-deltas)*(3d0-deltas)/4d0*dlog(1d0-deltas)
     &                 + deltas*(18d0-deltas)/24d0)

               else
c v+s+collinear part added (plus terms prop. deltas from soft part):
                  fqed(j) = alpha0/pi*(-(1d0+dlog(deltac/2d0))*
     $                 (le+(1d0-deltas)*(3d0-deltas)/4d0)-
     $                 pi**2/8d0+dreal(spence(deltas+ieps))-
     $                 (1d0-deltas)*(3d0-deltas)/4d0*
     $                 dlog(1d0-deltas)+65d0/24d0+
     $                 (1d0-deltas)*(19d0+deltas)/24d0)
               end if
            elseif (test(4).eq.2) then
               dqedf = alpha0/pi*(77d0/24d0-pi**2/3d0)
               if (collcut.eq.0) then
                  call beta(i,mw**2,that,betf)
                  dqedh = -betf*dlog(deltas)+alpha0/pi*
     $                 (qf**2*(-3d0/2d0*dlog(mw/mmf)-pi**2/6d0+
     $                 11d0/8d0)+qfs**2*(-3d0/2d0*dlog(mw/mmfs)-
     $                 pi**2/6d0+11d0/8d0)+5d0/6d0)
c terms prop. deltas added:
c               dqedh = alpha0/pi*((2d0-dlog(xmw/mmfs**2))*
c     $              (dlog(deltas)+(1d0-deltas)*(3d0-deltas)/4d0)-
c     $              pi**2/6d0+dreal(spence(deltas+ieps))-
c     $              (1d0-deltas)*(3d0-deltas)/4d0*
c     $              dlog(1d0-deltas)+
c     $              (1d0-deltas)*(17d0-deltas)/24d0)
               fqed(j) = dqedf-dqedh
               else
c hard part-collinear part(= colli. part added to v+s):
               dqedhc = alpha0/pi*((1d0+dlog(deltac/2d0*shat/xmw))*
     $                 (dlog(deltas)+(1d0-deltas)*(3d0-deltas)/4d0)+
     $                 pi**2/6d0-dreal(spence(deltas+ieps))+
     $                 (1d0-deltas)*(3d0-deltas)/4d0*
     $                 dlog(1d0-deltas)-
     $                 (1d0-deltas)*(19d0+deltas)/24d0)
               fqed(j) = dqedf-dqedhc
               end if
            else
            print*,'wrong choice for QED final state calc.'
            fqed(j) = 0d0
            end if
         end if
c-----i = 3------------------------------------------------------------
c interference:
         if (i.eq.3) then
         call beta(i,shat,that,betif)
         fqed(j) = betif*(le+lws)
         end if
      end do
      return
      end
c*********************************************************************
      subroutine virtsoftqcd_cc(fqcd)
c virtual+soft gluon contribution, eps = upper limit of photon energy
c sig(0+1) = sig0(g0->g(0+1))*(1+fqcd)
c soft part of PDF CT added
      implicit real*8(a-z)
      integer i,ii,j,jj
      real*8 fqcd(2)
      real*8 sinv(4,4)
      include 'pwhg_wzgrad.h'
c      include 'config.inc'
      real*8 pi,pi2
      integer i1(2)
      data i1/1,2/
      pi = 4d0*datan(1d0)
      pi = pi**2
      fqcd(1)=0d0
      fqcd(2)=0d0
      cf=3d0/4d0
c itwasme BEGIN
      do ii = 1,4
      do jj = 1,4
      sinv(ii,jj) = br_sinv(ii,jj)
      enddo
      enddo
c itwasme END 
      do j=1,2
         shat = sinv(1,2)
         that = -sinv(i1(j),4)
         eps = deltas*dsqrt(shat)/2d0
         le = dlog(deltas)
         dw = shat-xmw
         if (dw.eq.0d0) then
            phase = 0d0
            lw = dlog(gamw**2/(gamw**2+4d0*eps**2))/2d0
            lws = dlog(xmw/(gamw**2+4d0*eps**2))/2d0
         else
            phase = dw/mw/gamw*(datan(dw/mw/gamw)+
     $           datan((2d0*dsqrt(shat)*eps-dw)/mw/gamw))
            lw = dlog((dw**2+xmw*gamw**2)/(xmw*gamw**2+
     $           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
            lws = dlog(xmw**2/(xmw*gamw**2+
     $           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
         endif
c pdf counterterm added (and shift in propagator neglected):
c MSbar scheme
         fqcd(j) = alphas/pi*cf*(2d0*(le+3d0/4d0)*
     $        dlog(shat/mu_f**2)+pi2/3d0-7d0/2d0+2d0*le**2)
      end do
      return
      end
c******************************************************************
      subroutine beta(i,s,t,bet) 

      implicit real*8(a-z)
      integer i

      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs

      pi = 4d0*datan(1d0)
      u = -s-t
      s2 = s**2
      t2 = t**2
      u2 = u**2
      mi2 = mmi**2
      mis2 = mmis**2
      mf2 = mmf**2
      mfs2 = mmfs**2
      bet = 0d0

      if (i.eq.1) then 
         bet = alpha0/pi*(qi**2*(dlog(s/mi2)-1d0)+
     $        qis**2*(dlog(s/mis2)-1d0)-1d0)
      end if

      if (i.eq.2) then
         bet = alpha0/pi*(qf**2*(dlog(s/mf2)-1d0)+
     $        qfs**2*(dlog(s/mfs2)-1d0)-1d0)
      end if

      if (i.eq.3) then
         bet = alpha0/pi*((qi*qf+qis*qfs)*dlog(t2/s2)-
     $        (qis*qf+qi*qfs)*dlog(u2/s2)+2d0)
      end if

      return
      end
c********************************************************************
      subroutine delta(i,s,t,delt)
      implicit real*8(a-z)
      integer i
      complex*16 spence
      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
      real*8 pi,pi2
      u = -s-t
      spu = dreal(spence(1d0+s/u+ieps))
      spt = dreal(spence(1d0+s/t+ieps))
      mi2 = mmi**2
      mis2 = mmis**2
      mf2 = mmf**2
      mfs2 = mmfs**2
      s2 = s**2
      t2 = t**2
      u2 = u**2
      delt = 0d0
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      if (i.eq.1) then
         delt = alpha0/4d0/pi*(qi**2*(3d0/2d0*dlog(s/mi2)+
     &        pi2/3d0-2d0)+qis**2*(3d0/2d0*dlog(s/mis2)+
     &        pi2/3d0-2d0)+3d0+pi2/12d0)
      end if
      if (i.eq.2) then
         delt = alpha0/4d0/pi*(qf**2*(3d0/2d0*dlog(s/mf2) 
     &                                + pi2/3d0-2d0) 
     &                         + qfs**2*(3d0/2d0*dlog(s/mfs2)
     &                                + pi2/3d0-2d0)+3d0+pi2/12d0)

      end if
        
      if (i.eq.3) then
         delt = alpha0/4d0/pi*((qf*qi+qfs*qis)*(-2d0*spt+
     &        dlog(t2/s2)/2d0-dlog(t2/s2)**2/4d0)-
     &        (qfs*qi+qf*qis)*(-2d0*spu+
     &        dlog(u2/s2)/2d0-dlog(u2/s2)**2/4d0)-6d0-7d0/6d0*pi2)
      end if
      return
      end
*-------------------------------------------------------------------
      subroutine fvpureg(s,xq,xqs,xm,xms,form)
*
c pure photonic contribution:
c formfactors fvin(i)=fain(i) and fvfin(i) describing the 
c ith contribution to the initial or final state vertex correction.
c vertex: i lamda_mu(s) = i*e/2/w2/sw*g_mu*(1-g5)*fvertex(s)
c on-shell sing. terms have been subtracted!

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      real*8 pi,pi2
      complex*16 b0,c0,b0m,b0ms,b0s,spence
      complex*16 fv(5),factor,form,lws,cdlog
      common/renorm/mue,mue2,mue4
      q=xq
      qs=xqs
      m=xm
      ms=xms
      pi = 4d0*datan(1d0)
      pi2 = pi**2

      if(dabs(s-xmw).lt.1d-5)then
         lws=0d0
      else
         lws=(s-xmw)*cdlog((xmw-s)/xmw-ieps)
      endif
c g-exchange:      
      c0 = dcmplx(-dlog(s/m/ms)*dlog(lambda**2/s)-dlog(s/m**2)**2/4d0-
     $     dlog(s/ms**2)**2/4d0-2d0/3d0*pi2,
     $     pi*dlog(lambda**2/s))/s
      b0m=dcmplx(-dlog(m**2/mue2)+2d0,0d0)
      b0ms=dcmplx(-dlog(ms**2/mue2)+2d0,0d0)
      b0s=dcmplx(-dlog(s/mue2)+2d0,-pi)
      factor = -2d0*s*c0+2d0*b0m+2d0*b0ms-3d0*b0s-2d0
      fv(1) = alpha0/4d0/pi*q*qs*factor

c wgw-vertex:
      c0 = (-spence(1d0-xmw/s+ieps)-dlog(s/xmw)**2/2d0-pi**2/6d0)/s
      b0=-dlog(xmw/mue2)+2d0-lws/s
      factor = q*(2d0*xmw*c0+2d0*b0m+(2d0+xmw/s)*(-dlog(xmw/mue2)+1d0)-
     $     (1d0+xmw/s)*b0-2d0*dlog(s/m**2)*lws/s)-
     $     qs*(2d0*xmw*c0+2d0*b0ms+(2d0+xmw/s)*(-dlog(xmw/mue2)+1d0)-
     $     (1d0+xmw/s)*b0-2d0*dlog(s/ms**2)*lws/s)
      fv(2) = alpha0/4d0/pi*factor

c add on-shell singular terms:
c      fv(2)=fv(2)-alpha0/4d0/pi*2d0*(q*dlog(s/m**2)-qs*dlog(s/ms**2))*
c     $     cdlog(mw**2/dcmplx(mw**2-s,-mw*2.1d0))
*
c self energy of the external fermions:
c g:
      fv(3) = -1d0/2d0*alpha0/4d0/pi*q**2*
     $     dcmplx(-dlog(s/mue2)+3d0*dlog(s/m**2)+4d0+
     $     2d0*dlog(lambda**2/s),0d0)
      fv(4) = -1d0/2d0*alpha0/4d0/pi*qs**2*
     $     dcmplx(-dlog(s/mue2)+3d0*dlog(s/ms**2)+4d0+
     $     2d0*dlog(lambda**2/s),0d0)
*
c add contribution from W selfenergy:
      fv(5)=alpha0/4d0/pi/2d0*(-10d0/3d0*dlog(xmw/mue2)+68d0/9d0+
     $     2d0/3d0*(1d0-xmw/s)*(lws/s-1d0))
c add on-shell sing. term:
c      fv(5)=fv(5)+alpha0/4d0/pi/2d0*
c     $     4d0*cdlog(mw**2/dcmplx(mw**2-s,-mw*2.1d0))

c pure photonic contribution to the formfactors:
      form = dcmplx(0d0,0d0)
      do i=1,5
         form = form+fv(i)
      end do
      return
      end
c******************************************************************
      subroutine fvboxwg(s,t,form)

c pure photonic contribution to t/u-channel box diagram:
c on-shell sing. terms have been subtracted!

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
      real*8 pi,pi2
      complex*16 b0,b0t,b0u,d0,fvt,fvu,spence
      complex*16 factor,form,lws,cdlog

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2

      if(dabs(s-xmw).lt.1d-5)then
         lws=0d0
      else
         lws=(s-xmw)*cdlog((xmw-s)/xmw-ieps)
      endif
c t channel:
      d0=dlog(t**2/mmi**2/mmf**2)*dlog(lambda/mw)+
     $     dlog(mw/mmi)**2+dlog(mw/mmf)**2+spence(1d0+xmw/t+ieps)+
     $     pi2/3d0
      b0t=dcmplx(-dlog(dabs(t)/mue2)+2d0,0d0)
      b0=-dlog(xmw/mue2)+2d0-lws/s
      fvt = 2d0*(s+t)*(b0-b0t)+(2d0*t+s+xmw)*
     $     (spence(1d0+t/xmw+ieps)-spence(1d0+xmw/t+ieps)-
     $     2d0*spence(1d0-xmw/s+ieps)-dlog(t**2/s**2)**2/8d0-
     $     dlog(s/xmw)**2/2d0-
     $     dlog(t**2/xmw**2)*dlog(s/xmw)/2d0-2d0/3d0*pi2)+
     $     2d0*(t**2-s*xmw)/s*(-2d0*spence(1d0-xmw/s+ieps)-
     $     dlog(s/xmw)**2-pi2/3d0)
      factor = qi*qf*(2d0*d0+1d0/(s+t)**2*((s-xmw)*fvt+
     $     ((2d0*t+s+xmw)*dlog(s**2*t**2/mmf**4/mmi**4)+
     $     2d0*(t**2-s*xmw)/s*dlog(s**2/mmf**2/mmi**2))*lws))

      d0 = dlog(t**2/mmis**2/mmfs**2)*dlog(lambda/mw)+
     $     dlog(mw/mmis)**2+dlog(mw/mmfs)**2+spence(1d0+xmw/t+ieps)+
     $     pi2/3d0
      factor = factor+qis*qfs*(2d0*d0+1d0/(s+t)**2*((s-xmw)*fvt+
     $     ((2d0*t+s+xmw)*dlog(s**2*t**2/mmfs**4/mmis**4)+
     $     2d0*(t**2-s*xmw)/s*dlog(s**2/mmfs**2/mmis**2))*lws))

c u channel:
      u=-s-t
      d0=dlog(u**2/mmi**2/mmfs**2)*dlog(lambda/mw)+
     $     dlog(mw/mmi)**2+dlog(mw/mmfs)**2+spence(1d0+xmw/u+ieps)+
     $     pi2/3d0
      b0u=dcmplx(-dlog(dabs(u)/mue2)+2d0,0d0)
      fvu = 2d0*(s+u)*(b0-b0u)+(2d0*u+s+xmw)*
     $     (spence(1d0+u/xmw+ieps)-spence(1d0+xmw/u+ieps)-
     $     2d0*spence(1d0-xmw/s+ieps)-
     $     dlog(u**2/s**2)**2/8d0-dlog(s/xmw)**2/2d0-
     $     dlog(u**2/xmw**2)*dlog(s/xmw)/2d0-2d0/3d0*pi2)+
     $     2d0*(u**2-s*xmw)/s*(-2d0*spence(1d0-xmw/s+ieps)-
     $     dlog(s/xmw)**2-pi2/3d0)
      factor = factor-qi*qfs*(2d0*d0+1d0/(s+u)**2*((s-xmw)*fvu+
     $     ((2d0*u+s+xmw)*dlog(s**2*u**2/mmfs**4/mmi**4)+
     $     2d0*(u**2-s*xmw)/s*dlog(s**2/mmfs**2/mmi**2))*lws))

      d0 = dlog(u**2/mmis**2/mmf**2)*dlog(lambda/mw)+
     $     dlog(mw/mmis)**2+dlog(mw/mmf)**2+spence(1d0+xmw/u+ieps)+
     $     pi2/3d0
      factor = factor-qis*qf*(2d0*d0+1d0/(s+u)**2*((s-xmw)*fvu+
     $     ((2d0*u+s+xmw)*dlog(s**2*u**2/mmf**4/mmis**4)+
     $     2d0*(u**2-s*xmw)/s*dlog(s**2/mmf**2/mmis**2))*lws))
c add on-shell sing. terms:
c      factor=factor+2d0*(qi*qf*dlog(t**2/mmf**2/mmi**2)+
c     $     qis*qfs*dlog(t**2/mmfs**2/mmis**2)-
c     $     qi*qfs*dlog(u**2/mmfs**2/mmi**2)-
c     $     qis*qf*dlog(u**2/mmf**2/mmis**2))*
c     $     cdlog(mw**2/dcmplx(mw**2-s,-mw*2.1d0))

      form = alpha0/4d0/pi*factor
      return
      end
c******************************************************************
      subroutine fvyfs(s,t,form)
c photonic contribution, only modified YFS part
c sum of initial, final and interf.
c when adding soft BR contribution and the on-shell
c and mass sing. parts of the IR finite contribution, 
c it should coincide with the sofqed formfactors.

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      complex*16 spence

c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
      common/renorm/mue,mue2,mue4

      do i=1,1
c         mue=1d-1**(i**2)
c         mue2=mue**2
c         mue4=mue2**2

      u=-s-t
      form=0d0
c IR singular parts from diagrams I,II,V:
      factor = -dlog(s/mue2)+2d0*dlog(s/mmi/mmis)*(1d0+
     $     dlog(lambda**2/s))+2d0+dlog(s/mmi**2)**2/2d0+
     $     dlog(s/mmis**2)**2/2d0+4d0/3d0*pi2
      form = alpha0/4d0/pi*qi*qis*factor
      factor = -dlog(s/mue2)+2d0*dlog(s/mmf/mmfs)*(1d0+
     $     dlog(lambda**2/s))+2d0+dlog(s/mmf**2)**2/2d0+
     $     dlog(s/mmfs**2)**2/2d0+4d0/3d0*pi2
      form = form+alpha0/4d0/pi*qf*qfs*factor
      factor=-dlog(s/mue2)+3d0*dlog(s/mmi**2)+4d0+
     $     2d0*dlog(lambda**2/s)
      form = form-1d0/2d0*alpha0/4d0/pi*qi**2*factor
      factor=-dlog(s/mue2)+3d0*dlog(s/mmis**2)+4d0+
     $     2d0*dlog(lambda**2/s)
      form = form-1d0/2d0*alpha0/4d0/pi*qis**2*factor
      factor=-dlog(s/mue2)+3d0*dlog(s/mmf**2)+4d0+
     $     2d0*dlog(lambda**2/s)
      form = form-1d0/2d0*alpha0/4d0/pi*qf**2*factor
      factor=-dlog(s/mue2)+3d0*dlog(s/mmfs**2)+4d0+
     $     2d0*dlog(lambda**2/s)
      form = form-1d0/2d0*alpha0/4d0/pi*qfs**2*factor
      factor=-dlog(s/mue2)+dlog(t**2/mmf**2/mmi**2)*dlog(lambda**2/s)+
     $     dlog(s/mmf**2)**2/2d0+dlog(s/mmi**2)**2/2d0-
     $     dlog(t**2/s**2)**2/4d0+dlog(t**2/s**2)/2d0+
     $     2d0*dlog(s/mmf/mmi)+2d0+pi2/3d0
      form = form+alpha0/4d0/pi*qf*qi*factor
      factor=-dlog(s/mue2)+
     $     dlog(t**2/mmfs**2/mmis**2)*dlog(lambda**2/s)+
     $     dlog(s/mmfs**2)**2/2d0+dlog(s/mmis**2)**2/2d0-
     $     dlog(t**2/s**2)**2/4d0+dlog(t**2/s**2)/2d0+
     $     2d0*dlog(s/mmfs/mmis)+2d0+pi2/3d0
      form = form+alpha0/4d0/pi*qfs*qis*factor
      factor=-dlog(s/mue2)+dlog(u**2/mmfs**2/mmi**2)*dlog(lambda**2/s)+
     $     dlog(s/mmfs**2)**2/2d0+dlog(s/mmi**2)**2/2d0-
     $     dlog(u**2/s**2)**2/4d0+dlog(u**2/s**2)/2d0+
     $     2d0*dlog(s/mmfs/mmi)+2d0+pi2/3d0
      form = form-alpha0/4d0/pi*qfs*qi*factor
      factor=-dlog(s/mue2)+dlog(u**2/mmf**2/mmis**2)*dlog(lambda**2/s)+
     $     dlog(s/mmf**2)**2/2d0+dlog(s/mmis**2)**2/2d0-
     $     dlog(u**2/s**2)**2/4d0+dlog(u**2/s**2)/2d0+
     $     2d0*dlog(s/mmf/mmis)+2d0+pi2/3d0
      form = form-alpha0/4d0/pi*qf*qis*factor

c add terms subtracted from finite parts:
c mass singularities and term prop. q_j**2:
      form=form+alpha0/4d0/pi*(qi*2d0*dlog(s/mmi**2)-
     $     qis*2d0*dlog(s/mmis**2)+
     $     qf*2d0*dlog(s/mmf**2)-qfs*2d0*dlog(s/mmfs**2))
      form=form-2d0*(qi**2+qis**2+qf**2+qfs**2)*alpha0/4d0/pi
      form=form+alpha0/4d0/pi*(qi*qis*2d0*dlog(s/mmi/mmis)+
     $     qf*qfs*2d0*dlog(s/mmf/mmfs))
      form=form+alpha0/4d0/pi*(-2d0*qf*qi*dlog(s/mmf/mmi)-
     $     2d0*qfs*qis*dlog(s/mmfs/mmis)+2d0*qfs*qi*dlog(s/mmfs/mmi)+
     $     2d0*qf*qis*dlog(s/mmf/mmis))
c      write(6,*)'sum1',form
c (modified)YFS factors from PhD thesis:
c initial state
      form=alpha0/4d0/pi*(qi**2*((dlog(s/mmi**2)-1d0)*dlog(lambda**2/s)+
     $     dlog(s/mmi**2)**2/2d0+dlog(s/mmi**2)/2d0+2d0/3d0*pi2-2d0)+
     $     qis**2*((dlog(s/mmis**2)-1d0)*dlog(lambda**2/s)+
     $     dlog(s/mmis**2)**2/2d0+dlog(s/mmis**2)/2d0+2d0/3d0*pi2-2d0)+
     $     1d0+pi2/12d0-dlog(lambda**2/s))
c final state
      form=form+
     $     alpha0/4d0/pi*(qf**2*((dlog(s/mmf**2)-1d0)*dlog(lambda**2/s)+
     $     dlog(s/mmf**2)**2/2d0+dlog(s/mmf**2)/2d0+2d0/3d0*pi2-2d0)+
     $     qfs**2*((dlog(s/mmfs**2)-1d0)*dlog(lambda**2/s)+
     $     dlog(s/mmfs**2)**2/2d0+dlog(s/mmfs**2)/2d0+2d0/3d0*pi2-2d0)+
     $     1d0+pi2/12d0-dlog(lambda**2/s))
c interference:
      form=form+
     $     alpha0/4d0/pi*((qi*qf+qis*qfs)*
     $     (-dlog(t**2/s**2)**2/4d0+dlog(t**2/s**2)*
     $     (dlog(lambda**2/s)+1d0/2d0))-
     $     (qi*qfs+qis*qf)*
     $     (-dlog(u**2/s**2)**2/4d0+dlog(u**2/s**2)*
     $     (dlog(lambda**2/s)+1d0/2d0))-
     $     2d0-7d0/6d0*pi2+2d0*dlog(lambda**2/s))
c      write(6,*)'sum2',form
c subtract delta_interf._v+s:
      form=form-alpha0/4d0/pi*((qi*qf+qis*qfs)*
     $     (-dlog(t**2/s**2)**2/4d0+dlog(t**2/s**2)/2d0-
     $     2d0*dreal(spence(1d0+s/t+ieps)))-
     $     (qi*qfs+qis*qf)*
     $     (-dlog(u**2/s**2)**2/4d0+dlog(u**2/s**2)/2d0-
     $     2d0*dreal(spence(1d0+s/u+ieps)))-6d0-7d0/6d0*pi2)
c      write(6,*)mue,form
      enddo
      return
      end
c*******************************************************************
      subroutine fvsoftcoll(i,fqed)
c soft photon contribution, eps = upper limit of photon energy
c sig(0+1) = sig0(g0->g(0+1))*(1+fqed)
c collinear contribution added
      implicit real*8(a-z)
      integer i,ii,j,jj
      real*8 fqed(2)
      real*8 sinv(4,4)      
      complex*16 spence,cdlog
      include 'pwhg_wzgrad.h'
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
      real*8 pi,pi2
      integer i1(2)
      data i1/1,2/
      pi = 4d0*datan(1d0)
      pi2 = pi**2
c itwasme BEGIN
      do ii = 1,4
      do jj = 1,4
      sinv(ii,jj) = br_sinv(ii,jj)
      enddo
      enddo
c itwasme END  
      do j=1,2
         shat = sinv(1,2)
         that = -sinv(i1(j),4)
         uhat = -shat-that
         eps = deltas*dsqrt(shat)/2d0
         le = dlog(deltas)
         ll = dlog(dsqrt(shat)/lambda)
         dw = shat-xmw
         if (dw.eq.0d0) then
            phase = 0d0
            lw = dlog(gamw**2/(gamw**2+4d0*eps**2))/2d0
            lws = dlog(xmw/(gamw**2+4d0*eps**2))/2d0
         else
            phase = dw/mw/gamw*(datan(dw/mw/gamw)+
     &           datan((2d0*dsqrt(shat)*eps-dw)/mw/gamw))
            lw = dlog((dw**2+xmw*gamw**2)/(xmw*gamw**2+
     &           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
            lws = dlog(xmw**2/(xmw*gamw**2+
     &           (dw-2d0*dsqrt(shat)*eps)**2))/2d0
         endif
     
c initial state:

         if (i.eq.1) then

c pdf counterterm added (and shift in propagator neglected):
            
            fvps = 9d0+2d0/3d0*pi**2+3d0*le-2d0*le**2
            
            fqed(j) = alpha0/pi*(qi**2*((le+3d0/4d0)*
     $           dlog(mmi**2/mu_f**2)-1d0+le+le**2+lfc/4d0*fvps)+
     $           qis**2*((le+3d0/4d0)*
     $           dlog(mmis**2/mu_f**2)-1d0+le+le**2+lfc/4d0*fvps)+
     $           qi**2*((dlog(shat/mmi**2)-1d0)*(le+ll+
     $           0d0*(lw+phase))+
     $           1d0/4d0*dlog(shat/mmi**2)*(2d0-dlog(shat/mmi**2))-
     $           pi2/6d0)+
     $           qis**2*((dlog(shat/mmis**2)-1d0)*(le+ll+
     $           0d0*(lw+phase))+
     $           1d0/4d0*dlog(shat/mmis**2)*(2d0-dlog(shat/mmis**2))-
     $           pi2/6d0)-(le+ll+0d0*(lw+phase))+1d0)

c            fqed(j) = alpha0/pi*(qi**2*((dlog(shat/mmi**2)-1d0)*(le+ll+
c     $           (lw+phase))+
c     $           1d0/4d0*dlog(shat/mmi**2)*(2d0-dlog(shat/mmi**2))-
c     $           pi2/6d0)+
c     $           qis**2*((dlog(shat/mmis**2)-1d0)*(le+ll+
c     $           (lw+phase))+
c     $           1d0/4d0*dlog(shat/mmis**2)*(2d0-dlog(shat/mmis**2))-
c     $           pi2/6d0)-(le+ll+(lw+phase))+1d0)
         end if
*
c final state:
*
         if (i.eq.2) then

c the following expressions are only valid for a leptonic final state
            
            if (collcut.eq.0) then

               fqed(j) = alpha0/pi*(qf**2*((dlog(shat/mmf**2)-1d0)*
     $              (le+ll)+1d0/4d0*dlog(shat/mmf**2)*
     $              (2d0-dlog(shat/mmf**2))-pi2/6d0)+
     $              qfs**2*((dlog(shat/mmfs**2)-1d0)*
     $              (le+ll)+1d0/4d0*dlog(shat/mmfs**2)*
     $              (2d0-dlog(shat/mmfs**2))-pi2/6d0)-(le+ll)+1d0)

            else
c collinear part added (plus terms prop. deltas from soft part):

               fqed(j) = alpha0/pi*(qf**2*((dlog(shat/mmf**2)-1d0)*
     $              (le+ll)+1d0/4d0*dlog(shat/mmf**2)*
     $              (2d0-dlog(shat/mmf**2))-pi2/6d0)+
     $              qfs**2*((dlog(shat/mmfs**2)-1d0)*
     $              (le+ll)+1d0/4d0*dlog(shat/mmfs**2)*
     $              (2d0-dlog(shat/mmfs**2))-pi2/6d0)-(le+ll)+1d0)

               fqed(j) = fqed(j)+alpha0/pi*(9d0/4d0-pi2/3d0+le-
     $              deltas*(deltas+2d0)/4d0+deltas*(18d0-deltas)/24d0-
     $              (3d0-4d0*deltas+deltas**2)*dlog(1d0-deltas)/4d0+
     $              dreal(spence(deltas+ieps))+
     $              (-3d0/4d0+0d0*(deltas-deltas**2/4d0)-le)*
     $              (qf**2*dlog(shat/mmf**2)+
     $              qfs**2*dlog(shat/mmfs**2))+
     $              (-3d0/4d0+deltas-deltas**2/4d0-le)*
     $              (qf**2+qfs**2)*dlog(deltac/2d0))

            end if
         end if
*
c interference:
c on-shell sing. terms have been canceled by 
c on-shell sing. parts of F_V^gamma(s) !
*
         if (i.eq.3) then

            fqed(j) = alpha0/pi*((qi*qf+qis*qfs)*
     $           (dlog(that**2/shat**2)*(lws+le+ll)-
     $           dreal(spence(1d0+shat/that+ieps)))-
     $           (qis*qf+qi*qfs)*
     $           (dlog(uhat**2/shat**2)*(lws+le+ll)-
     $           dreal(spence(1d0+shat/uhat+ieps)))+
     $           2d0*(lws+le+ll-1d0))
c add on-shell sing. terms:
c            fqed(j)=fqed(j)+alpha0/pi*((qi*qf+qis*qfs)*
c     $           (dlog(that**2/shat**2))-
c     $           (qis*qf+qi*qfs)*
c     $           (dlog(uhat**2/shat**2))+
c     $           2d0)*dreal(cdlog(dcmplx(mw**2-shat,-mw*2.1d0)/mw**2))
         end if
*
      end do
      return
      end






