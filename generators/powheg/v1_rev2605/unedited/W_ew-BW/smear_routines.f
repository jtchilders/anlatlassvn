c***********************************************************************
      subroutine smear(p1,p2,p3,sp1,sp2,sp3)
c.....
c     this subroutine smears the FS momenta 
c     input: p1,p2,p3 
c     output: sp1,sp2,sp3 
c.....
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      character*1 lep
      real*8 p1(4),p2(4),p3(4)
      real*8 sp1(4),sp2(4),sp3(4),scalet
      common/scalet/scalet
      real*8 powheginput
      external powheginput
      logical flg_saverand
      common/pwhg_flg_saverand/flg_saverand
!     =========================================
!     electron or muon ?
!     =========================================
      if(powheginput('vdecaymode').eq.1)lep='e'
      if(powheginput('vdecaymode').eq.2)lep='m'
!     ==========================================
!     call smearing routines 
!     ==========================================
      scalet=dsqrt(p1(1)**2+p1(2)**2)
      if(powheginput('ih2').eq.-1)then !TEVATRON
      call d0upgrsmear(1,p1,sp1,lep) !lepton
      call metsmear(2,p2,sp2)        !neutrino
      call d0upgrsmear(3,p3,sp3,'e') !photon
      elseif(powheginput('ih2').eq.1)then !LHC
      call atlsmear(1,p1,sp1,lep)    !lepton
      call atlsmear(2,p2,sp2,'n')    !neutrino
      call atlsmear(3,p3,sp3,'e')    !photon
      else
      print*, 'unknown type for hadron 2, error in real_EW.f'
      stop
      endif
      flg_saverand=.false.
      end
c***********************************************************************
      subroutine d0upgrsmear(ip,p,sp,lep)
c.....
c     this subroutine smear the 3-momentum component of the 4-vector p
c     and is suitable for Tevatron 
c     input: p, 4-vector to be smeared
c     lep  : particle type, 'e' for electron/photon, 'h' for hadron,
c            'n' for missing pt (neutrino)
c     output: sp, smeared 4-vector
c.....
      implicit none
      integer i,j,ip
      character*1 lep
      real*8 scalet
      common/scalet/scalet
      real*8 p(4),sp(4)
      real*8 ratio,eta,c,s,en,et,a,b,ret,del,dpx,dpy,dp

c     =======================================
c     calculate eta and Et - used to calc del
c     =======================================
      if(p(4).ne.p(3))then
         ratio = (p(4)+p(3))/(p(4)-p(3))
         if(ratio.gt.0d0)then
         eta=dabs(0.5d0*dlog((p(4)+p(3))/(p(4)-p(3))))
         else
         eta=10d0
         endif
      else
      eta=10d0
      endif
      et=dsqrt(p(1)**2+p(2)**2)
c     =====================================
c     calculate del for electron
c     =====================================
      if(lep.eq.'e')then !setup
         c=0.003d0
         if(eta.lt.1.3d0)then
         s=0.14d0
         en=0.14d0
         elseif(eta.ge.1.3d0.and.eta.lt.4d0)then
         s=0.157d0
         en=0.290d0
         endif
         if(eta.lt.4d0)then
         del=dsqrt((en/p(4))**2+s**2/p(4)+c**2)*p(4)
         elseif(eta.ge.4d0)then
         del=0d0
         endif
      goto 9 !now smear
      endif
c     =====================================
c     calculate del for muon
c     =====================================
      if(lep.eq.'m')then !setup
         a=0.015d0
         b=0.0016d0
         if(eta.lt.1.6d0)then
         del=dsqrt((b*et)**2+a**2)*et
         else
         del=0d0
         endif
      goto 9 !now smear
      endif
c     =====================================
c     calculate del for jets
c     =====================================
      if(lep.eq.'h')then !setup
      del=dsqrt((2.16d0/et)**2+0.74d0**2/et+0.01d0**2)*et
      goto 9 !now smear
      endif
c     =====================================
c     perform smearing for neutrinos
c     =====================================
      if(lep.eq.'n')then
      ret=scalet
      del=1.08d0+0.019d0*ret
      call rgauss(ip,1,del,dpx)
      call rgauss(ip,2,del,dpy)
      sp(1)=p(1)+dpx
      sp(2)=p(2)+dpy
      sp(3)=p(3)
      sp(4)=dsqrt(sp(1)**2+sp(2)**2+sp(3)**2)
      goto 15
      endif
c     =====================================
c     perform smearing for everthing else
c     =====================================
    9 call rgauss(ip,0,del,dp) 
      if(p(4).gt.0d0.and.dp.gt.-p(4))then
      sp(4)=p(4)+dp
         do i=1,3
         sp(i)=p(i)*(1d0+dp/p(4))
         enddo
      else
      sp(4)=1d-6*p(4)
         do j=1,3
         sp(j)=1d-6*p(j)
         enddo
      endif
   15 return
      end
c***********************************************************************
      subroutine metsmear(ip,p,sp)

c     this subroutine smear the pt-momentum component of the 4-vector p
c     of a neutrino

c     input:
c     p    : 4-vector to be smeared

c     output: 
c     sp    : smeared 4-vector
      implicit none
      integer i,ip
      real*8 p(4),sp(4),pi,x,et,del,dp,random

      pi=4d0*datan(1d0)
      et=dsqrt(p(1)**2+p(2)**2)
      del=5d0

      call rgauss(ip,0,del,dp)

c     if et=0, then smear in a random direction
c if calo=1, also bare=1 so that ET=0 
c is always cut, so it is ok not to save the random number here
      if(et.eq.0d0)then
      x=random()
      sp(1)=dp*dcos(x*pi)
      sp(2)=dp*dsin(x*pi)
      else
         do i=1,2
         sp(i)=p(i)*(1d0+dp/et)
         enddo
      endif
    
      sp(3)=p(3)
      sp(4)=dsqrt(sp(3)**2+sp(2)**2+sp(1)**2)

      end
c***********************************************************************
      subroutine atlsmear(ip,p,sp,lep)
c.....
c     this subroutine smear the 3-momentum component of the 4-vector p
c     and is suitable for Tevatron 
c     input: p, 4-vector to be smeared
c     lep  : particle type, 'e' for electron/photon, 'h' for hadron,
c            'n' for missing pt (neutrino)
c     output: sp, smeared 4-vector
c.....

      implicit none
      character*1 lep
      integer i,j,ip
      real*8 scalet
      common/scalet/scalet
      real*8 p(4),sp(4)
      real*8 diff,sum,eta,arg,et,del,ret,dpx,dpy,dp
      
      diff=p(4)-p(3)
      sum=p(4)+p(3)

      if(dabs(diff).lt.1d-5.or.dabs(sum).lt.1d-5)then
      eta=20d0
      else
      arg=sum/diff
         if(arg.le.0d0)then
         eta=20d0
         elseif(arg.gt.0d0)then
         eta=dabs(0.5d0*dlog((p(4)+p(3))/(p(4)-p(3))))
         endif
      endif

      et=dsqrt(p(1)**2+p(2)**2)

      if(lep.eq.'e')then  !electron
         if(eta.lt.2.5d0)then
         del=dsqrt(0.095d0**2/p(4)+0.005d0**2)*p(4)
         elseif(eta.ge.2.5d0)then
         del=0d0
         endif 
      endif

      if(lep.eq.'m')then  !muon
         if(eta.lt.2d0)then
         del=dsqrt((5d-4*et)**2+0.012d0**2)*et*
     &       (1d0+eta**10/7000d0)
         elseif(eta.ge.2d0)then
         del=0d0
         endif
      endif

      if(lep.eq.'n')then  !neutrino
      ret=scalet
      del=0.45d0*dsqrt(ret)
      call rgauss(ip,1,del,dpx)
      call rgauss(ip,2,del,dpy)
      sp(1)=p(1)+dpx
      sp(2)=p(2)+dpy
      sp(3)=p(3)
      sp(4)=dsqrt(sp(1)**2+sp(2)**2+sp(3)**2)
      goto 15
      endif

    9 call rgauss(ip,0,del,dp)
      
       if(p(4).gt.0d0.and.dp.gt.-p(4))then
       sp(4)=p(4)+dp
          do i=1,3
          sp(i)=p(i)*(1d0+dp/p(4))
          enddo
       else
       sp(4)=1d-15*p(4)
          do j=1,3
          sp(j)=1d-15*p(j)
          enddo
       endif
   15  return  
       end
c***********************************************************************
      subroutine rgauss(ip,icomp,del,dp)

c     yields gaussian distribution of random numbers
      implicit none
      integer ip,icomp
      real*8 x(2),random,xsave(3,0:2,2)
      real*8 del,dp,gauss,test 
      common/saverandom/xsave
      logical flg_saverand
      common/pwhg_flg_saverand/flg_saverand

      if(del.le.0d0)then !del must be positive 
      dp=0d0
      return
      endif

    5 x(1) = random()
      x(2) = random()
c the same random numbers are used for smearing of momenta
c belonging to one event
      if(flg_saverand)then
         xsave(ip,icomp,1)=x(1)
         xsave(ip,icomp,2)=x(2)
      else
         x(1)=xsave(ip,icomp,1)
         x(2)=xsave(ip,icomp,2)
      endif

      dp=3d0*del*(-1d0+2d0*x(1))
      gauss=dexp(-dp**2/2d0/del**2)
      test=1d0*x(2)-gauss

      if(test.le.0d0)then
      return
      else
         if(.not.flg_saverand)then
            write(6,*)'test failed in smearing routine'
            stop
         endif
      goto 5  
      endif
      end
c***********************************************************************
      subroutine delR(sp1,sp2,sp3,dR)
      implicit none
      integer i
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      
      real*8 phil,phip,del_phi,pi,ptl,ptp
      real*8 thl,thp,yl,yp,dR,sp1(4),sp2(4),sp3(4)

      integer flag1,flag2
      common/leptonidf/flag1,flag2

      logical flg_inbtilde,flg_inequiv
      common/pwhg_flg_EW/flg_inbtilde,flg_inequiv

      real*8 powheginput
      external powheginput

      pi = 4d0*datan(1d0)
      !define delta_R
      phil = datan2(sp1(2),sp1(1))
      phip = datan2(sp3(2),sp3(1))
      if(phil.lt.0.d0)phil=phil+2d0*pi
      if(phip.lt.0.d0)phip=phip+2d0*pi

      del_phi = dabs(phil-phip)
      if(del_phi.gt.pi)del_phi=2d0*pi-del_phi

      ptl = dsqrt(sp1(1)**2+sp1(2)**2)
      ptp = dsqrt(sp3(1)**2+sp3(2)**2)

      thl = datan2(ptl,sp1(3))      
      thp = datan2(ptp,sp3(3))      

      yl = -dlog(dtan(thl/2d0))
      yp = -dlog(dtan(thp/2d0))
 
      dR = dsqrt(del_phi**2+(yl-yp)**2)

cimpose delta_R cuts only if calculating real EW MES in Btilde
      flag2=0
      if(flg_inbtilde)then
         if(powheginput('vdecaymode').eq.1)then
c only recombine when a hard photon is present
            if(dR.lt.0.1d0.and.flag1.eq.0)then
               do i=1,4
                  sp1(i) = sp1(i)+sp3(i)
                  sp3(i) = 0d0
               end do
            elseif((dR.ge.0.1d0).and.(dR.le.0.4d0))then
               if(sp3(4).gt.0.1d0*sp1(4))then
                  flag2=1
               endif
            endif
         elseif(powheginput('vdecaymode').eq.2)then
            if((dR.lt.0.1d0).and.(sp3(4).gt.2d0))then
            flag2=1
            endif
            if((dR.ge.0.1d0).and.(dR.le.0.4d0))then
               if(sp3(4).gt.0.1d0*sp1(4))then
               flag2=1
               endif
            endif
         endif
      endif
      end  
c***********************************************************************
      function dotpr(p1,p2)
      implicit none
      real * 8 dotpr,p1(0:3),p2(0:3)
      dotpr = (p1(0)*p2(0) - p1(3)*p2(3)) - p1(1)*p2(1) - p1(2)*p2(2)
      end

