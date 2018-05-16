c*******************************************************************
      subroutine formweak(i,fvweak)
c calculation of the (modified) weak 1-loop corrections 
c to the w-exchange process: 
c    i(p4)+i'(p3)->w+(s=q^2)->f(p2)+f'(p1).
c i=1: initial state contribution, i=2: final state
      implicit none
      include 'pwhg_wzgrad.h'
      complex*16 fvweak,dweak
      real*8 g0,dqed,dqcd,dr,gwi,gwf,g0f,g0i,gw1f,gw1i,ncf
      integer i,ff,ffp,ii,iip
*
      common/split/g0,dqed,dweak,dqcd
      common/drit/dr 
*
c QED-substracted partial W-width 
*
      if (i.eq.0) fvweak = (0d0,0d0)
      if (i.eq.1) then
         ii = 7
         iip = 10
         call gwffp(ii,iip,g0i,gw1i)
         gwi = g0i*(1d0+dreal(dweak)-(rep-1)*dr)
         fvweak = dweak-(rep-1)*dr
      end if
      if (i.eq.2) then
         ff = 1
         ffp = 4
         ncf = 1d0
         call gwffp(ff,ffp,g0f,gw1f)
         gwf = g0f*(1d0+dreal(dweak)-(rep-1)*dr+(ncf-1d0)/2d0*dqcd)
         fvweak = dweak-(rep-1)*dr
      end if
      return
      end
c******************************************************************
      real*8 function a2qqw(xs,xt,xu,sigborn)
c     full (resonant+non-resonant) EW 1-loop corrections to the 
c     W-exchange process:  i(p1)+i'(p2)->W+(s=q^2)->f=vu(p4)+f'(p3).
      implicit none
      include 'pwhg_wzgrad.h'
      integer i
      real*8 boxq,boxt_wz,boxt_zw,boxu_wz,boxu_zw
      real*8 boxt_v1v2,boxu_v1v2,sigborn,sigq
      real*8 t1,t2,sp,mmff,mni,xs,xt,xu
      real*8 xvf,xaf,xvi,xai,vap,vam
      real*8 vf12,vi34,af12,ai34
      real*8 xvfs,xafs,xvis,xais,xm,xnc,xqm
      real*8 fvweak,fvgam,denw,dz2w,ftfinite,fufinite
      real*8 fvu,fvt,w,b0,b0t,b0u,lws
      complex*16 swwos,sww,dswwos,dsww,fvwff,fvwii
      complex*16 spence,cdlog
      real*8 pi,pi2
      real*8 mue,mue2,mue4
      common/renorm/mue,mue2,mue4
      real*8 dr
      common/drit/dr 
c      real*8 mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      a2qqw=0d0
      sigq=0d0
      boxq=0d0
      sp=xs
      t1=xt
      t2=xu
      denw=sp-xmw
c v,a couplings of fermions to the Z boson
c u db -> W^+ ->v_l l:
      call prop(7,xm,xqm,xvi,xai,xnc)
      call prop(10,xm,xqm,xvis,xais,xnc)
      call prop(1,xm,xqm,xvf,xaf,xnc)
      call prop(4,xm,xqm,xvfs,xafs,xnc)
c      if (mz.gt.0d0) goto 99
c-------------------------------------------------------------
c matrix element squared with weak vertex and self energy corrections
c to W+ production: full s-dependence !
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
c pure weak contribution evaluated:
c check of mu-dependence
c      do i=1,3
c         mue=1d1**(i**2)
c         mue2=mue*mue
c         mue4=mue2*mue2
c      if(mz.gt.0d0) goto 99
c weak part to W vertex and selfenergy corrections
      call wselfw(xmw,swwos,dswwos)
      call wselfw(sp,sww,dsww)
      call fvpurew(sp,xvf,xaf,xvfs,xafs,fvwff)
      call fvpurew(sp,xvi,xai,xvis,xais,fvwii)
      if(dabs(denw).lt.1d-5) then
         fvweak = dreal(fvwff+fvwii-dswwos)
      else
         fvweak = dreal(fvwff+fvwii-(sww-swwos)/denw)
      end if
 99   continue
c IR finite (virtual) photon contribution: 
c full virtual photon contribution where the (modified) YFS factors
c have been subtracted, plus the
c delta^interf_vs term from the YFS prescription, 
c which in the resonant part 
c had been absorbed into the modifed weak contribution.
c 2*delta gamma_rem=fvgam(s=Mw**2):
      fvgam=alpha0/4d0/pi*(-3d0*pi2/2d0+68d0/9d0-
     $     25d0/3d0*dlog(xmw/mue2))
      if(dabs(denw).lt.1d-5)then
         lws=0d0
      else
         lws=(sp-xmw)*dreal(cdlog((xmw-sp)/xmw-ieps))
      endif
      w=xmw/sp
c term proportional to (sp-xmw)/(s+t)^2 of the finite
c box contribution f_V,t:
      b0t=-dlog(dabs(t1)/mue2)+2d0
      b0=-dlog(xmw/mue2)+2d0-lws/sp
      fvt = 2d0*(sp+t1)*(b0-b0t)+(2d0*t1+sp+xmw)*
     $     (dreal(spence(1d0+t1/xmw+ieps))-
     $     dreal(spence(1d0+xmw/t1+ieps))-
     $     2d0*dreal(spence(1d0-xmw/sp+ieps))-dlog(t1**2/sp**2)**2/8d0-
     $     dlog(sp/xmw)**2/2d0-
     $     dlog(t1**2/xmw**2)*dlog(sp/xmw)/2d0-2d0/3d0*pi2)+
     $     2d0*(t1**2-sp*xmw)/sp*(-2d0*dreal(spence(1d0-xmw/sp+ieps))-
     $     dlog(sp/xmw)**2-pi2/3d0)
      ftfinite = 1d0/(sp+t1)**2*((sp-xmw)*(qi*qf+qis*qfs)*fvt+
     $     qi*qf*((2d0*t1+sp+xmw)*dlog(sp**2*t1**2/mmf**4/mmi**4)+
     $     2d0*(t1**2-sp*xmw)/sp*dlog(sp**2/mmf**2/mmi**2))*lws+
     $     qis*qfs*((2d0*t1+sp+xmw)*dlog(sp**2*t1**2/mmfs**4/mmis**4)+
     $     2d0*(t1**2-sp*xmw)/sp*dlog(sp**2/mmfs**2/mmis**2))*lws)
c
c term proportional to (sp-xmw)/(s+u)^2 of the finite
c box contribution f_V,u:
      b0u=-dlog(dabs(t2)/mue2)+2d0
      fvu = 2d0*(sp+t2)*(b0-b0u)+(2d0*t2+sp+xmw)*
     $     (dreal(spence(1d0+t2/xmw+ieps))-
     $     dreal(spence(1d0+xmw/t2+ieps))-
     $     2d0*dreal(spence(1d0-xmw/sp+ieps))-
     $     dlog(t2**2/sp**2)**2/8d0-dlog(sp/xmw)**2/2d0-
     $     dlog(t2**2/xmw**2)*dlog(sp/xmw)/2d0-2d0/3d0*pi2)+
     $     2d0*(t2**2-sp*xmw)/sp*(-2d0*dreal(spence(1d0-xmw/sp+ieps))-
     $     dlog(sp/xmw)**2-pi2/3d0)
      fufinite = -1d0/(sp+t2)**2*((sp-xmw)*(qi*qfs+qis*qf)*fvu+
     $     qi*qfs*((2d0*t2+sp+xmw)*dlog(sp**2*t2**2/mmfs**4/mmi**4)+
     $     2d0*(t2**2-sp*xmw)/sp*dlog(sp**2/mmfs**2/mmi**2))*lws+
     $     qis*qf*((2d0*t2+sp+xmw)*dlog(sp**2*t2**2/mmf**4/mmis**4)+
     $     2d0*(t2**2-sp*xmw)/sp*dlog(sp**2/mmf**2/mmis**2))*lws)

      fvgam=fvgam+alpha0/4d0/pi*
     $     ((qi*qf+qis*qfs)*(2d0*dreal(spence(1d0+xmw/t1+ieps)-
     $     spence(1d0+sp/t1+ieps))+
     $     dlog(w)*(dlog(w)-dlog(t1**2/sp**2)))+ftfinite-
     $     (qis*qf+qi*qfs)*(2d0*dreal(spence(1d0+xmw/t2+ieps)-
     $     spence(1d0+sp/t2+ieps))+
     $     dlog(w)*(dlog(w)-dlog(t2**2/sp**2)))+fufinite+
     $     2d0*((1d0-w)/3d0*(2d0+pi2)+2d0/3d0*(2d0+w)*lws/sp-
     $     2d0*w*dreal(spence(1d0-w+ieps))-w*dlog(w)**2-dlog(w)-
     $     (qf*dlog(sp/mmf**2)-qfs*dlog(sp/mmfs**2)+
     $     qi*dlog(sp/mmi**2)-qis*dlog(sp/mmis**2))*lws/sp)+
     $     5d0*dlog(w))
c
c add the W self energy CT (photon+weak part):
 999  continue
      call wrenorm(dz2w)
      sigq=2d0*(fvgam+fvweak-dz2w-(rep-1)*dr)*sigborn
c      write(6,*)mue,fvgam+fvweak-dz2w-(rep-1)*dr
c      enddo
c      if (mz.gt.0d0) goto 9999
c-------------------------------------------------------------
c matrix element squared for W,Z boxes to W production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
      mni=1d-2
      mmff=1d-2
c WZ box:
      boxt_wz=0d0
      vf12=xvfs+xafs
      vi34=xvis+xais
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxt_wz=boxt_v1v2(t1,sp,mmff,mni,mmff,mw,mni,mz,vap,vam)
c ZW box:
      vf12=xvf+xaf
      vi34=xvi+xai
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxt_zw=boxt_v1v2(t1,sp,mmff,mni,mmff,mz,mni,mw,vap,vam)
c WZ box (crossed):
      vf12=xvf+xaf
      vi34=xvis+xais
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxu_wz=boxu_v1v2(t2,sp,mmff,mni,mmff,mw,mni,mz,vap,vam)
c ZW box (crossed):
      vf12=xvfs+xafs
      vi34=xvi+xai
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxu_zw=boxu_v1v2(t2,sp,mmff,mni,mmff,mz,mni,mw,vap,vam)
c sum of box diagrams
      boxq=alpha0/4d0/pi*8d0*(sp-xmw)*
     $     (boxt_wz+boxt_zw+boxu_wz+boxu_zw)
c sum:
 9999 continue
      a2qqw=sigq+boxq

      return
      end
c******************************************************************
      real*8 function a2qqw_check(xs,xt,xu,sigborn)
c     full (resonant+non-resonant)
c     electroweak 1-loop corrections to the 
c     W-exchange process:  i(p4)+i'(p3)->W+(s=q^2)->f=vu(p2)+f'(p1).

      implicit none
      include 'pwhg_wzgrad.h'
      integer i
      real*8 boxq,boxt_wz,boxt_zw,boxu_wz,boxu_zw
      real*8 boxt_v1v2,boxu_v1v2,sigborn,sigq
      real*8 t1,t2,sp,mmff,mmii,xs,xt,xu
      real*8 xvf,xaf,xvi,xai,vap,vam
      real*8 vf12,vi34,af12,ai34
      real*8 xvfs,xafs,xvis,xais,xm,xnc,xqm
      real*8 fvweak,fvg,denw,dz2w
      complex*16 swwos,sww,dswwos,dsww,fvwff,fvwii
      complex*16 fvgff,fvgii,fvbox
      complex*16 spence
      real*8 pi,pi2
      real*8 mue,mue2,mue4
      common/renorm/mue,mue2,mue4
      real*8 dr
      common/drit/dr 
c      real*8 mmi,mmis,mmf,mmfs,qi,qis,qf,qfs
c      common/ferm/mmi,mmis,mmf,mmfs,qi,qis,qf,qfs

      pi=4d0*datan(1d0)
      pi2 = pi**2
      a2qqw_check=0d0
      sigq=0d0
      boxq=0d0
      sp=xs
      t1=xt
      t2=xu
      denw=sp-xmw
c v,a couplings of fermions to the Z boson
c u db -> W^+ ->v_l l:
      call prop(7,xm,xqm,xvi,xai,xnc)
      call prop(10,xm,xqm,xvis,xais,xnc)
      call prop(1,xm,xqm,xvf,xaf,xnc)
      call prop(4,xm,xqm,xvfs,xafs,xnc)
c-------------------------------------------------------------
c matrix element squared with weak vertex and self energy corrections
c to W+ production: full s-dependence !
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
c pure weak contribution evaluated at shat:
c check of mue-dependence
c      do i=1,3
c         mue=1d1**(i**2)
c         mue2=mue*mue
c         mue4=mue2*mue2
c weak part to W vertex and selfenergy corrections
      call wselfw(xmw,swwos,dswwos)
      call wselfw(sp,sww,dsww)
      call fvpurew(sp,xvf,xaf,xvfs,xafs,fvwff)
      call fvpurew(sp,xvi,xai,xvis,xais,fvwii)
      if(dabs(denw).lt.1d-5) then
         fvweak = dreal(fvwff+fvwii-dswwos)
      else
         fvweak = dreal(fvwff+fvwii-(sww-swwos)/denw)
      end if
c pure photonic contribution:
      call fvpureg(sp,qf,qfs,mmf,mmfs,fvgff)      
      call fvpureg(sp,qi,qis,mmi,mmis,fvgii)
c photonic box contribution:
      call fvboxwg(sp,t1,fvbox)
      fvg = dreal(fvgff+fvgii+fvbox)
c add the W self energy CT (photon+weak part):
 999  continue
      call wrenorm(dz2w)
      sigq=2d0*(fvg+fvweak-dz2w-(rep-1)*dr)*sigborn
c      write(6,*)mue,sigq
c      enddo
c      if (mz.gt.0d0) goto 9999
c-------------------------------------------------------------
c matrix element squared for W,Z boxes to W production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
      mmii=1d-2
      mmff=1d-2
c WZ box:
      boxt_wz=0d0
      vf12=xvfs+xafs
      vi34=xvis+xais
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxt_wz=boxt_v1v2(t1,sp,mmff,mmii,mmff,mw,mmii,mz,vap,vam)
c ZW box:
      vf12=xvf+xaf
      vi34=xvi+xai
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxt_zw=boxt_v1v2(t1,sp,mmff,mmii,mmff,mz,mmii,mw,vap,vam)
c WZ box (crossed):
      vf12=xvf+xaf
      vi34=xvis+xais
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxu_wz=boxu_v1v2(t2,sp,mmff,mmii,mmff,mw,mmii,mz,vap,vam)
c ZW box (crossed):
      vf12=xvfs+xafs
      vi34=xvi+xai
      af12=vf12
      ai34=vi34
      vap=vf12*vi34+af12*ai34
      vam=0d0
      boxu_zw=boxu_v1v2(t2,sp,mmff,mmii,mmff,mz,mmii,mw,vap,vam)
c sum of box diagrams
      boxq=alpha0/4d0/pi*16d0*(sp-xmw)*
     $     (boxt_wz+boxt_zw+boxu_wz+boxu_zw)
c sum:
 9999 continue
      a2qqw_check=sigq+boxq
      return
      end
c******************************************************************
      real*8 function a2qqz_weak(xs,xt,xu,xvf,xaf,xqf,xvi,xai,xqi)
c     calculation of the weak 1-loop corrections to the z,gamma-exchange
c     process:  i(p4)+i'(p3)->gamma,z(s=q^2)->f(p2)+f'(p1).

      implicit none
      include 'pwhg_wzgrad.h'
      integer i,mix
      real*8 boxq,box_ww,boxt_zz,boxu_zz
      real*8 boxt_v1v2,boxu_v1v2
      real*8 t1,t2,sp,den,denz,mmff,mmii,xs,xt,xu
      real*8 xvf,xaf,xqf,xvi,xai,xqi,vap,vam
      real*8 vf12,vi34,af12,ai34,v12,a12,v34,a34
      real*8 sigz_weak,sigg_weak,siggz_weak,sigq_weak,sigq_weak_all
      real*8 xvfs,xafs,xqfs,xvis,xais,xqis,xm,xnc
      real*8 alaq1,alaq2,alaq3,blbq1,blbq2,blbq3
      real*8 mue,mue2,mue4
      real*8 dz2z,dz2g,dz2gz,dsigzz
      real*8 pi,pi2
      real*8 pig0,sgz0
      real*8 born,a2qqz
      complex*16 zero
      complex*16 rsiggg,rsiggz
      complex*16 dsiggz,dsiggg,sigzz,siggg,siggz,sigzos
      complex*16 piggs,pigzs,pizs,pizzs
      complex*16 fvzf,fvzi,fvgf,fvgi,gazf,gazi,gagf,gagi
      complex*16 dgzmix
      complex*16 vfg,afg,vig,aig,vfz,afz,viz,aiz
      complex*16 vfgc,afgc,vigc,aigc,vfzc,afzc,vizc,aizc
      complex*16 denzc

      common/renorm/mue,mue2,mue4
      common/zgrenorm/dz2z,dz2g,dz2gz
      common/sigzero/pig0,sgz0
      common/par_prop/den
      common/sigzon/sigzos
      common/derivmix/dgzmix
      common/par_cprop/denzc
      common/gzmix/mix
                
      pi=4d0*datan(1d0)
      pi2 = pi**2
      sigq_weak_all=0d0
      boxq=0d0

      sp=xs
      t1=xt
      t2=xu
c-------------------------------------------------------------
c matrix element squared with weak vertex and self energy corrections
c to Z,gamma production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
      ieps=dcmplx(0d0,1d-8)
      zero=dcmplx(0d0,0d0)
*
c renormalised self energies:
      rsiggg=siggg(sp)+sp*dz2g
      rsiggz=siggz(sp)-sgz0-sp*dz2gz
      if (dabs(sp).lt.1d-10) then
         piggs=0d0
         pigzs=dsiggz(0d0)-dz2gz
      else
         piggs=rsiggg/sp
         pigzs=rsiggz/sp
      end if
      denz=sp-xmz
c if mix=0, then pizzs is ok, but if
c mix=1, then pizs must be used.
      if(dabs(denz).lt.0.001d0) then
         pizzs=dsigzz(xmz)+dz2z
         if(mix.ne.0)then
            pizs=pizzs+dgzmix
         else
            pizs=pizzs
         endif
      else
         pizzs=(sigzz(sp)-dreal(sigzos)+denz*dz2z)/denz
         if(mix.ne.0)then
            pizs=pizzs-rsiggz**2/(sp+rsiggg)/denz
         else
            pizs=pizzs
         endif
      end if
c check
c      do i=1,100
c         sp=xmz+i*0.001d0
c         rsiggg=siggg(sp)+sp*dz2g
c         rsiggz=siggz(sp)-sgz0-sp*dz2gz
c         denz=sp-xmz
c         pizzs=dsigzz(xmz)+dz2z
c         pizs=pizzs+dgzmix
c         write(6,*)'1',denz,pizzs,pizs
c         pizzs=(sigzz(sp)-dreal(sigzos)+denz*dz2z)/denz
c         pizs=pizzs-rsiggz**2/(sp+rsiggg)/denz
c         write(6,*)'2',denz,pizzs,pizs
c      enddo
c      stop
*
c renormalised vertex functions:
      call prop(1,xm,xqfs,xvfs,xafs,xnc)
      if (xqi.gt.0d0)
     $     call prop(10,xm,xqis,xvis,xais,xnc)
      if (xqi.lt.0d0)
     $     call prop(7,xm,xqis,xvis,xais,xnc)
      call formgz(sp,xvf,xaf,xqf,xvfs,xafs,xqfs,fvzf,fvgf,gazf,gagf)
      call formgz(sp,xvi,xai,xqi,xvis,xais,xqis,fvzi,fvgi,gazi,gagi)

      if(mz.gt.0d0) goto 90
c
c strict 1-loop, i.e. no quadratic terms:
c      piggs=0d0
c      pigzs=0d0
c      pizs=0d0
c      gazf=0d0
c      gazi=0d0
c      fvzf=0d0
c      fvzi=0d0
c      gagf=0d0
c      gagi=0d0
c      fvgi=0d0
c      fvgf=0d0

      alaq1=(xai**2+xvi**2)*(xaf*dreal(gazf)+xvf*dreal(fvzf))+
     $     (xaf**2+xvf**2)*(xai*dreal(gazi)+xvi*dreal(fvzi))-
     $     ((xaf**2+xvf**2)*(xai**2+xvi**2)*dreal(pizzs)-
     $     ((xaf**2+xvf**2)*xvi*xqi+
     $     (xai**2+xvi**2)*xvf*xqf)*dreal(pigzs))
      blbq1=-2d0*xvi*xai*(xvf*dreal(gazf)+xaf*dreal(fvzf))-
     $     2d0*xaf*xvf*(xvi*dreal(gazi)+xai*dreal(fvzi))-
     $     (-4d0*xvf*xaf*xvi*xai*dreal(pizzs)+
     $     (2d0*xvf*xaf*xai*xqi+2d0*xai*xvi*xaf*xqf)*dreal(pigzs))
      sigz_weak=16.d0*(alaq1*(t1**2+t2**2)+blbq1*(t1**2-t2**2))
      alaq2=-xqi*xqf*(xqi*dreal(fvgf)+xqf*dreal(fvgi))-
     $     (xqi*xqf)**2*dreal(piggs)
      blbq2=0d0
      sigg_weak=16.d0*(alaq2*(t1**2+t2**2)+blbq2*(t1**2-t2**2))
      alaq3=-xvi*xqi*(xaf*dreal(gagf)+xvf*dreal(fvgf))-
     $     xvf*xqf*(xai*dreal(gagi)+xvi*dreal(fvgi))+
     $     xqi*xqf*(xvi*dreal(fvzf)+xvf*dreal(fvzi))-
     $     (xqi*xqf*xvi*xvf*dreal(piggs)+xvf*xvi*xqi*xqf*dreal(pizzs)-
     $     xqi*xqf*(xvf*xqi+xvi*xqf)*dreal(pigzs))
      blbq3=xai*xqi*(xvf*dreal(gagf)+xaf*dreal(fvgf))+
     $     xaf*xqf*(xvi*dreal(gagi)+xai*dreal(fvgi))-
     $     xqi*xqf*(xai*dreal(gazf)+xaf*dreal(gazi))+
     $     xqi*xqf*xai*xaf*dreal(piggs)+xaf*xai*xqi*xqf*dreal(pizzs)
      siggz_weak=16.d0*(alaq3*(t1**2+t2**2)+blbq3*(t1**2-t2**2))
*
      sigq_weak_all=16d0*alpha0**2*(sigz_weak/den+sigg_weak/sp**2+
     $     (sp-mz**2)/den/sp*siggz_weak)
 90   continue
c      if(mz.gt.0d0) goto 99
*
c factorised weak corrections:
c only vertex corrections:
c      pizs=zero
c      pigzs=zero
c      gazf=zero
c      gazi=zero
c      fvzf=zero
c      fvzi=zero
c      gagf=zero
c      gagi=zero
c      fvgi=zero
c      fvgf=zero

      afz=xaf+gazf
      aiz=xai+gazi
      vfz=xvf+fvzf+xqf*pigzs/(1d0+piggs)
      viz=xvi+fvzi+xqi*pigzs/(1d0+piggs)
      afg=-gagf
      aig=-gagi
      vfg=xqf-fvgf
      vig=xqi-fvgi
c conjugate complex:
      afzc=dcmplx(dreal(afz),-dimag(afz))
      aizc=dcmplx(dreal(aiz),-dimag(aiz))
      vfzc=dcmplx(dreal(vfz),-dimag(vfz))
      vizc=dcmplx(dreal(viz),-dimag(viz))
      afgc=dcmplx(dreal(afg),-dimag(afg))
      aigc=dcmplx(dreal(aig),-dimag(aig))
      vfgc=dcmplx(dreal(vfg),-dimag(vfg))
      vigc=dcmplx(dreal(vig),-dimag(vig))
*
      alaq1= (cdabs(afz)**2+cdabs(vfz)**2)*
     $     (cdabs(aiz)**2+cdabs(viz)**2)
      blbq1=-4d0*dreal(vfz*afzc)*dreal(viz*aizc)
      sigz_weak=8.d0*(alaq1*(t1**2+t2**2)+blbq1*(t1**2-t2**2))
      alaq2=(cdabs(afg)**2+cdabs(vfg)**2)*
     $     (cdabs(aig)**2+cdabs(vig)**2)
      blbq2=-4d0*dreal(vfg*afgc)*dreal(vig*aigc)
      sigg_weak=8.d0*(alaq2*(t1**2+t2**2)+blbq2*(t1**2-t2**2))
      alaq3=2d0*dreal(denzc*(afg*afzc+vfg*vfzc)*(aig*aizc+vig*vizc))
      blbq3=-2d0*dreal(denzc*(afg*vfzc+vfg*afzc)*(aig*vizc+vig*aizc))
      siggz_weak=8.d0*(alaq3*(t1**2+t2**2)+blbq3*(t1**2-t2**2))
*
      sigq_weak_all=16d0*alpha0**2*
     $     (sigz_weak/den/(1d0+dreal(pizs))**2+
     $     sigg_weak/sp**2/(1d0+dreal(piggs))**2+
     $     1d0/den/sp/(1d0+dreal(pizs))/(1d0+dreal(piggs))*siggz_weak)
 99   continue
c      if (mz.gt.0d0) goto 9999
c-------------------------------------------------------------
c matrix element squared for weak boxes to Z,gamma production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
      mmii=1d-2
      mmff=1d-2
c WW box:
      box_ww=0d0
c      if (mz.gt.0d0) goto 999
      vf12=2d0*(1d0/2d0/w2/sw)**2*(xvf+xaf)
      vi34=2d0*(1d0/2d0/w2/sw)**2*(xvi+xai)
      af12=vf12
      ai34=vi34
      v12=2d0*(1d0/2d0/w2/sw)**2
      v34=v12
      a12=v12
      a34=v12
      vap=(vf12*vi34+af12*ai34)*(sp-mz**2)/den+xqf*xqi/sp*
     $     (v12*v34+a12*a34)
      vam=(vf12*vi34-af12*ai34)*(sp-mz**2)/den+xqf*xqi/sp*
     $     (v12*v34-a12*a34)
      if (xqi.lt.0d0)
     $     box_ww=boxt_v1v2(t1,sp,mmff,mmii,1d-5,mw,mmii,mw,vap,vam)
      if (xqi.gt.0d0)
     $     box_ww=boxu_v1v2(t2,sp,mmff,mmii,1d-5,mw,mmii,mw,vap,vam)
c 999  continue
c ZZ box:
      boxt_zz=0d0
      boxu_zz=0d0
c      if (mz.gt.0d0) goto 999
      vf12=(xvf**2+xaf**2)*xvf+2d0*xvf*xaf**2
      vi34=(xvi**2+xai**2)*xvi+2d0*xvi*xai**2
      af12=(xvf**2+xaf**2)*xaf+2d0*xvf**2*xaf
      ai34=(xvi**2+xai**2)*xai+2d0*xvi**2*xai
      v12=xvf**2+xaf**2
      v34=xvi**2+xai**2
      a12=2d0*xvf*xaf
      a34=2d0*xvi*xai
      vap=(vf12*vi34+af12*ai34)*(sp-mz**2)/den+xqf*xqi/sp*
     $     (v12*v34+a12*a34)
      vam=(vf12*vi34-af12*ai34)*(sp-mz**2)/den+xqf*xqi/sp*
     $     (v12*v34-a12*a34)
      boxt_zz=boxt_v1v2(t1,sp,mmff,mmii,mmff,mz,mmii,mz,vap,vam)
c ZZ box (crossed):
      boxu_zz=boxu_v1v2(t2,sp,mmff,mmii,mmff,mz,mmii,mz,vap,vam)
c 999  continue
c sum of box diagrams
      boxq=alpha0**3/pi*128d0*(box_ww+boxt_zz+boxu_zz)
c sum:
 9999 continue
      a2qqz_weak=sigq_weak_all+boxq
      return
      end
*-------------------------------------------------------------------
      subroutine gztot(g0tot,gamz1)
*
c determination of the total Z-width
c (in alpha (rep=1) or G_mu representation (rep=2),
c including qcd corrections (qcd=1) in the limit of massless 
c decay products,
c mzero = 0 means the fermion masses have been taken into account
c in the tree level expression (kap is not 0 )):
*
      implicit real*8(a-z)
      integer f,fp
*
      gz0l = 0d0
      gz1l = 0d0
      do f = 1,3
         call gzffp(f,f+3,gz0,gz1)
         gz0l = gz0l+gz0
         gz1l = gz1l+gz1
c         write(6,*)'nu',gz0,gz1
      end do
      do f = 4,5
         call gzffp(f,f-3,gz0,gz1)
         gz0l = gz0l+gz0
         gz1l = gz1l+gz1
         if (f.eq.4) gz1e=gz1
c      write(6,*)'l',gz0,gz1
      end do
c modification for Z->tau tau:
      gz1t=gz1e-0.18d-3
      call gzffp(6,3,gz0,gz1)
c      write(6,*)gz0,gz1t
      gz0l = gz0l+gz0
      gz1l = gz1l+gz1t
      gz0q = 0d0
      gz1q = 0d0
      do f = 7,8
         call gzffp(f,f+3,gz0,gz1)
         gz0q = gz0q+gz0
         gz1q = gz1q+gz1
c         write(6,*)gz0,gz1
      end do
      do f = 10,11
         call gzffp(f,f-3,gz0,gz1)
         gz0q = gz0q+gz0
         gz1q = gz1q+gz1
         if (f.eq.10) gz1d=gz1
c         write(6,*)gz0,gz1
      end do
c modification for Z->bb:
      gz1b=gz1d-8.8d-3
      call gzffp(12,9,gz0,gz1)
c      write(6,*)gz0,gz1b
      gz0q = gz0q+gz0
      gz1q = gz1q+gz1b
c total Z-width:
      g0tot = gz0l+gz0q
      gamz1 = gz1l+gz1q
*
      return
      end
c******************************************************************
      subroutine formgz(s,v,a,q,vs,as,qs,fvz,fvg,gaz,gag)
c     formfactors describing the weak contribution to the (Z,gamma)
c     ffbar vertex correction.
c     vertex: i lamda_mu(s) = i*e*g_mu*(fv(s)-ga(s) g5)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      real*8 pi,pi2
      complex*16 b0,b1,c0,spence,cscal
      complex*16 factor,fvz,fvg,gaz,gag
      complex*16 vfz,vfw,vfws,vfww
      common/renorm/mue,mue2,mue4
      common/sigzero/pig0,sgz0

      pi=4d0*datan(1d0)
      pi2 = pi**2
c     z-exchange:
      w=xmz/s
      c0=dcmplx(dreal(spence(-1d0/w+ieps))-dlog(w)*dlog(1d0+1d0/w),
     $     -pi*dlog(1d0+1d0/w))
      factor=dcmplx(-dlog(xmz/mue2)-(2d0*w+3d0)*dlog(w)-2d0*w-
     $     4d0,-pi*(2d0*w+3d0))-2d0*(1d0+w)**2*c0
      vfz=alpha0/4d0/pi*factor
c self energy of the external fermions:
      vfz=vfz+alpha0/4d0/pi*dcmplx(dlog(xmz/mue2)+1d0/2d0,0d0)
c w-exchange:
      w=xmw/s
      c0=dcmplx(dreal(spence(-1d0/w+ieps))-dlog(w)*dlog(1d0+1d0/w),
     $     -pi*dlog(1d0+1d0/w))
      factor=dcmplx(-dlog(xmw/mue2)-(2d0*w+3d0)*dlog(w)-2d0*w-
     $     4d0,-pi*(2d0*w+3d0))-2d0*(1d0+w)**2*c0
      vfw=alpha0/4d0/pi*factor
c self energy of the external fermions:
      vfws=alpha0/4d0/pi*dcmplx((dlog(xmw/mue2)+1d0/2d0),0d0)
c z(g)ww-vertex:
      c0=cscal(s,1d-1,mw,mw,0d0)
      call bquer(s,mw,mw,b0,b1)
      factor=-3d0*dlog(xmw/mue2)+2d0*w+4d0-(2d0*w+1d0)*b0+
     $     2d0*xmw*(2d0+w)*c0
      vfww=alpha0/4d0/pi*factor
*
      t3f=2d0*cw*sw*a
      fvz=(v*(v**2+a**2)+a*2d0*v*a)*vfz+
     $     (vs+as)/4d0/sw**2*vfw+(v+a)/4d0/sw**2*vfws-
     $     (-2d0*t3f)*cw/sw/4d0/sw**2*vfww-cw/sw*a*sgz0/xmz
      fvg=-q*(v**2+a**2)*vfz-qs/4d0/sw**2*vfw+(-2d0*t3f)*
     $     1d0/4d0/sw**2*vfww-q/4d0/sw**2*vfws+a*sgz0/xmz
      gaz=(a*(v**2+a**2)+v*2d0*v*a)*vfz+
     $     (vs+as)/4d0/sw**2*vfw+(v+a)/4d0/sw**2*vfws-
     $     (-2d0*t3f)*cw/sw/4d0/sw**2*vfww-cw/sw*a*sgz0/xmz
      gag=-q*2d0*v*a*vfz-qs/4d0/sw**2*vfw-q/4d0/sw**2*vfws+
     $     (-2d0*t3f)/4d0/sw**2*vfww+a*sgz0/xmz
*
      return
      end
c******************************************************************
      subroutine gzffp(f,fp,gammatree,gammaloop)
c     calculation of the fermionic partial Z-width Z->f fbar: tree level,
c     1-loop EW and higher order QCD (+resummation); 
c     gammaloop = gammatree*(1+corr.)
c     as taken from Wolf's PhD thesis !

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer f,fp
      complex*16 fvzf,gazf,fvgf,gagf,zm,gvf,gaf
      
      common/renorm/mue,mue2,mue4
      common/drit/deltar
      common/zgrenorm/dz2z,dz2g,dz2gz
      common/waverenorm/zz,zm
      common/s2eff/sl2eff,su2eff,sd2eff

      pi=4d0*datan(1d0)
      pi2 = pi**2

      if (f.eq.9) then
         gammatree = 0d0
         gammaloop = 0d0
         return
      end if
*
c choice of the representation: rep=0,1: alpha-repr. and
c rep=2: gmu-repr.:
*
      call prop(f,m,q,v,a,nc)
      if (rep.eq.2) then
         g0z = mz**3/3d0/pi*w2*gfermi*sw2*cw2
      else
         g0z = alpha0*mz/3d0
      end if
      m2 = m**2
      kap = m2/xmz
*
      if (mzero.eq.1) kap = 0d0
      gammatree = g0z*nc*dsqrt(1d0-4d0*kap)*((1d0+2d0*kap)*v**2+
     $     (1d0-4d0*kap)*a**2)
*
c QED contribution (virt.+real(soft+hard)):
      dqedz = alpha0/pi*3d0/4d0*q**2
c weak contribution:
c renormalised vertex form factors:
      call prop(fp,ms,qs,vs,as,nc)
      call formgz(xmz,v,a,q,vs,as,qs,fvzf,fvgf,gazf,gagf)
*
c      if (f.eq.12) goto 99

      gvf=v+fvzf-q*zm
      gaf=a+gazf
      if (f.eq.4)
     $     sl2eff=1d0/(4d0*dabs(q))*(1d0-dreal(gvf)/dreal(gaf))
      if (f.eq.7)
     $     su2eff=1d0/(4d0*dabs(q))*(1d0-dreal(gvf)/dreal(gaf))
      if (f.eq.10)
     $     sd2eff=1d0/(4d0*dabs(q))*(1d0-dreal(gvf)/dreal(gaf))

      dweakz = zz*((1d0+2d0*kap)*cdabs(gvf)**2+
     $     (1d0-4d0*kap)*cdabs(gaf)**2)
*
      gammaloop = nc*g0z*dsqrt(1d0-4d0*kap)*dweakz*(1d0+dqedz)
c      if (rep.eq.2) gammaloop = gammaloop*(1d0-deltar)
*
c QCD correction:
      dqcdz = 0d0
      if ((qcd.eq.1).and.(f.ge.7)) then 
         dqcdz = alphas/pi*(1d0+1.405d0*alphas/pi-
     $        12.8d0*(alphas/pi)**2-q**2/4d0/pi*alpha0)
         gammaloop = gammaloop*(1d0+dqcdz)
      end if
*
 99   continue
c      if (f.eq.12) then
c         ll = dlog(xmz/m**2)
c         zt = xmz/4d0/mtop**2
c         aspi=alphas/pi
c         ii = -9.25d0+1.037d0*zt+0.0632d0*zt**2
c     $        +6d0*log(2d0*sqrt(zt))
c         rv = 12d0*m**2/xmz*(aspi+(6.07d0-2d0*ll)*aspi**2
c     $        +(2.38d0-24.29*ll+0.083d0*ll**2)*aspi**3)
c         ra = 6d0*m**2/xmz*((2d0*ll-1d0)*aspi
c     $        +(17.96d0+log(mtop**2/xmz)+14.14d0*ll
c     $        -0.083d0*ll**2)*aspi**2)
c     $        +ii/3d0*aspi**2
c         wqcd  = nc*g0z*zz*(rv*cdabs(v+fvzf-q*zm)**2+
c     $        ra*cdabs(a+gazf)**2)
c         gammaloop=nc*g0z*zz*(cdabs(v+fvzf-q*zm)**2+
c     $        cdabs(a+gazf)**2)*(1d0+dqcdz)*(1d0+dqedz)+wqcd
c         write(6,*)wqcd,gammaloop
c      end if
      return
      end

*-----------------------------------------------------------------
      double precision function boxt_v1v2(xt,xs,xmf,xmi,xmfs,
     $     xmv1,xmis,xmv2,xvap,xvam)
*-----------------------------------------------------------------
c t-channel box contribution with two vector bosons v1,v2 and
c two fermions i,f. f and i are considered to be massless !
c mi=mf=0 !
*-----------------------------------------------------------------
      implicit none
      real*8 xt,xs,xmv1,xmv2,xvap,xvam,xmf,xmi,xmfs,xmis
      real*8 t,u,s,mv1,mv2,mf,mi,mfs,mis
      real*8 d1(0:3),d2(0:6),d3(0:3,0:3),ch1(0:2)
      real*8 d0,d11,d12,d13,d20,d21,d22,d23,d212,d213,d223
      s=xs
      t=xt
      u=-s-t
      mf=xmf
      mi=xmi
      mv1=xmv1
      mv2=xmv2
      mfs=xmfs
      mis=xmis
      call dmn(mf**2,0d0,0d0,t,u,mf**2,mfs,mv1,mis,mv2,d1,d2,d3,ch1)
      d0 = d1(0)
      d11 = d1(1)
      d12 = d1(2)
      d13 = d1(3)
      d20 = d2(0)
      d21 = d2(1)
      d22 = d2(2)
      d23 = d2(3)
      d212 = d2(4)
      d213 = d2(5)
      d223 = d2(6)
      boxt_v1v2=xvap*u**2*(2d0*d20+t*(d11+d12+d13+d22+d223+d212))+
     $     xvam*t**2*(8d0*d20+t*(d11+2d0*d12+d13+
     $     2d0*(d22+d223+d212))-2d0*s*d213)
      return
      end
*-----------------------------------------------------------------
      double precision function boxu_v1v2(xu,xs,xmf,xmi,xmfs,
     $     xmv1,xmis,xmv2,xvap,xvam)
*-----------------------------------------------------------------
c u-channel box contribution with two vector bosons v1,v2 and
c two fermions i,f. f and i are considered to be massless !
c mi=mf=0 !
*-----------------------------------------------------------------
      implicit none
      real*8 xu,xs,xmv1,xmv2,xvap,xvam,xmf,xmi,xmis,xmfs
      real*8 t,u,s,mv1,mv2,mf,mi,mis,mfs
      real*8 d1(0:3),d2(0:6),d3(0:3,0:3),ch1(0:2)
      real*8 d0,d11,d12,d13,d20,d21,d22,d23,d212,d213,d223
      s=xs
      u=xu
      t=-s-u
      mf=xmf
      mi=xmi
      mv1=xmv1
      mv2=xmv2
      mfs=xmfs
      mis=xmis
      call dmn(mf**2,0d0,0d0,u,t,mf**2,mfs,mv1,mis,mv2,d1,d2,d3,ch1)
      d0 = d1(0)
      d11 = d1(1)
      d12 = d1(2)
      d13 = d1(3)
      d20 = d2(0)
      d21 = d2(1)
      d22 = d2(2)
      d23 = d2(3)
      d212 = d2(4)
      d213 = d2(5)
      d223 = d2(6)
      boxu_v1v2=-(xvam*t**2*(2d0*d20+u*(d11+d12+d13+d22+d223+d212))+
     $     xvap*u**2*(8d0*d20+u*(d11+2d0*d12+d13+
     $     2d0*(d22+d223+d212))-2d0*s*d213))
      return
      end
*-------------------------------------------------------------------
      subroutine gwtot(g0tot,gamw1)
*
c determination of the total W-width
c (in alpha (rep =1) or G_mu representation (rep = 2),
c including qcd corrections (qcd=1) in the limit of massless 
c decay products,
c mzero = 0 means the fermion masses have been taken into account
c in the tree level expression (kap is not 0 )):
*
      implicit real*8(a-z)
      integer f,fp
*
      gw0l = 0d0
      gw1l = 0d0
      do f = 1,3
         call gwffp(f,f+3,gw0,gw1)
         gw0l = gw0l+gw0
         gw1l = gw1l+gw1
      end do
      gw0q = 0d0
      gw1q = 0d0
      do f = 7,8
         do fp = 10,12
            call gwffp(f,fp,gw0,gw1)
            gw0q = gw0q+gw0
            gw1q = gw1q+gw1
         end do
      end do
c total W-width:
      g0tot = gw0l+gw0q
      gamw1 = gw1l+gw1q
*
      return
      end
c******************************************************************
      subroutine fvpurew(s,v,a,vs,as,form)
c     formfactors fvin(i)=fain(i) and fvfin(i) describing the 
c     ith contribution to the initial or final state vertex correction.
c     vertex: i lamda_mu(s) = i*e/2/w2/sw*g_mu*(1-g5)*fvertex(s)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      complex*16 b0,b1,c0,spence,cscal
      complex*16 fv(5),factor,form

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
c z-exchange:      
      w = xmz/s
      c0 = dcmplx(dreal(spence(-1d0/w+ieps))-dlog(w)*dlog(1d0+1d0/w),
     $     pi*dlog(1d0+1d0/w))/s
      factor = dcmplx(-dlog(xmz/mue2)-(2d0*w+3d0)*dlog(w)-2d0*w-
     $     4d0,-pi*(2d0*w+3d0))-2d0*(1d0+w)**2*s*c0
      fv(1) = alpha/4d0/pi*(v+a)*(vs+as)*factor
*
c wzw-vertex:
      c0 = cscal(s,1d-1,mz,mw,0d0)
      call bquer(s,mz,mw,b0,b1)
      factor = -3d0/2d0*dlog(xmz*xmw/mue4)-1d0/2d0*(xmw-xmz)/s*
     $     dlog(xmw/xmz)+4d0+(xmw+xmz)/s-((xmw+xmz)/s+1d0)*b0+
     $     2d0*(xmz+xmw+xmz*xmw/s)*c0
      fv(2) = alpha/4d0/pi*cw/sw*(v+a-vs-as)*factor
*
c self energy of the external fermions:
c Z:
      fv(3) = alpha/4d0/pi*((v+a)**2+(vs+as)**2)*
     $     dcmplx((dlog(xmz/mue2)+1d0/2d0)/2d0,0d0)
c w:
      fv(4) = alpha/4d0/pi*dcmplx((dlog(xmw/mue2)+
     $     1d0/2d0)/2d0/sw2,0d0)
*
c counterterm: dz1_w-dz2_w:
      fv(5) = alpha/4d0/pi*dcmplx(2d0/sw2*dlog(xmw/mue2),0d0)
*
c pure weak contribution to the formfactors:
      form = dcmplx(0d0,0d0)
      do i=1,5
         form = form+fv(i)
      end do
*
      return
      end
c******************************************************************
      complex*16 function sigzz(s)
c     Z-self energy (has been checked by comparison with Wolf)
        
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 azferm,azneut,azboson,cf,b0h,b1h,b20h,b0w,b1w,b20w

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      azferm = dcmplx(0d0,0d0)
      azneut = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cf = dcmplx(f(s,m,m),g(s,m,m))
         azferm = azferm+
     $        4d0/3d0*nc*((v**2+a**2)*(-s*dlog(m2/mue2)+
     $        (2d0*m2+s)*cf-s/3d0)-a**2*6d0*m2*(-dlog(m2/mue2)+cf))
      end do
      do ff=1,3
         call prop(ff,m,q,v,a,nc)
         azneut = azneut+
     $        8d0/3d0*nc*a**2*s*dcmplx(-dlog(s/mue2)+5d0/3d0,pi)
      end do
      call bquer(s,mw,mw,b0w,b1w)
      b20w = (xmw*b0w+s/2d0*b1w-7d0/24d0*s+3d0/2d0*xmw)/3d0
      b20w = b20w-(6d0*xmw-s)/12d0*dlog(xmw/mue2)
      b0w = b0w-dlog(xmw/mue2)
      a0w = xmw*(-dlog(xmw/mue2)+1d0)
      azboson = -cw2/sw2*(10d0*b20w+2d0*a0w+
     $     2d0*(xmw+2d0*s)*b0w-4d0*xmw+2d0/3d0*s)
      azboson = azboson+2d0*cw2/sw2*(3d0*a0w-2d0*xmw)
      azboson = azboson+xmw*sw2/cw2*2d0*b0w
      azboson = azboson+(sw2-cw2)**2/2d0/sw2/cw2*(a0w-2d0*b20w)
      azboson = azboson+cw2/sw2*b20w*2d0
      call bquer(s,mz,mh,b0h,b1h)
      b20h = (xmz*b0h+(s+xmz-xmh)/2d0*b1h-7d0/24d0*s+3d0/8d0*
     $     (xmz+3d0*xmh))/3d0
      b20h = b20h-(xmz+3d0*xmh-s)/12d0*dlog(xmh/mue2)-
     $     xmz/6d0*dlog(xmz/mue2)
      b0h = b0h-dlog(xmh*xmz/mue4)/2d0
      a0h = xmh*(-dlog(xmh/mue2)+1d0)
      a0z = xmz*(-dlog(xmz/mue2)+1d0)
      azboson = azboson-(b20h-xmz*b0h-a0h/4d0-
     $     a0z/4d0)/cw2/sw2
*
      sigzz = alpha/4d0/pi*(azferm+azneut+azboson)
      return
      end
c******************************************************************
      subroutine wselfw(s,sigww,dsigww)
c     W-selfenergy : (has been checked by comparison with Wolf)
c     only weak part!
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 b0z,b1z,b20z,db20z,b0h,b1h,b20h,db20h,b0f,b1f,
     &     b20f,db20f
      complex*16 awlep,dawlep,awquark,dawquark,awboson,dawboson,
     &     sigww,dsigww

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
c weak contribution to A_ww:    
      awlep = dcmplx(0d0,0d0)
      dawlep = dcmplx(0d0,0d0)
      do ff=4,6
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         m4 = m2**2
         if (s.eq.0d0) then
            awlep = awlep+
     $           nc*m2*dcmplx(1d0/2d0*dlog(m2/mue2)-1d0/4d0,0d0)/sw2
            dawlep = 0d0
         else
            awlep = awlep+
     $           nc*dcmplx(-(s-3d0/2d0*m2-m4/s)*dlog(m2/mue2)+
     $           (s-m2/2d0-m4/s/2d0)*(m2/s-1d0)*dlog(dabs(1d0-s/m2))+
     $           5d0/3d0*s-m2-3d0/2d0*m4/s,-pi*theta(s-m2)*(s-m2/2d0-
     $           m4/s/2d0)*(m2/s-1d0))/3d0/sw2
            dawlep = dawlep+
     $           nc*dcmplx(-(1d0+m4/s**2)*dlog(m2/mue2)+
     $           ((m2/s)**2/2d0+(m2/s)**3-1d0)*dlog(dabs(1d0-s/m2))+
     $           2d0/3d0+2d0*m4/s**2+m2/2d0/s,-pi*theta(s-m2)*(1d0-
     $           (m2/s)**3))/3d0/sw2
         endif
      end do
      awquark = dcmplx(0d0,0d0)
      dawquark = dcmplx(0d0,0d0)
      do  ff=7,9
         call prop(ff,m,q,v,a,nc)
         call prop(ff+3,ms,qs,vs,as,ncs)
         m2 = m**2
         ms2 = ms**2
         call bquer(s,m,ms,b0f,b1f)
         b20f = (m2*b0f+(s+m2-ms2)/2d0*b1f-7d0/24d0*s+3d0/8d0*
     $        (m2+3d0*ms2))/3d0
         b20f = b20f-(m2+3d0*ms2-s)/12d0*dlog(ms2/mue2)-
     $        m2/6d0*dlog(m2/mue2)
         b0f = b0f-dlog(m2*ms2/mue4)/2d0
         a0f = m2*(-dlog(m2/mue2)+1d0)
         a0fs = ms2*(-dlog(ms2/mue2)+1d0)
         awquark = awquark+nc*(2d0*b20f-a0f/2d0
     $        -a0fs/2d0-(m2+ms2-s)*b0f/2d0)/sw2
*
         if (s.eq.0d0) then
            dawquark = 0d0
            goto 100
         else
            call bquer1(s,m,ms,p0,p1,db0f,db1f)
            call bquer(s,m,ms,b0f,b1f)
            db20f = (m2*db0f+(s+m2-ms2)/2d0*db1f+b1f/2d0-
     $           7d0/24d0)/3d0
            db20f = db20f+dlog(ms2/mue2)/12d0
            b0f = b0f-dlog(m2*ms2/mue4)/2d0
            dawquark = dawquark+nc*(2d0*db20f-(m2+ms2-s)*db0f/2d0+
     $           b0f/2d0)/sw2
         end if
 100     continue
      end do
      call bquer(s,mz,mw,b0z,b1z)
      call bquer(s,mw,mh,b0h,b1h)
      b20z = (xmz*b0z+(s+xmz-xmw)/2d0*b1z-7d0/24d0*s+3d0/8d0*
     $     (xmz+3d0*xmw))/3d0
      b20h = (xmw*b0h+(s+xmw-xmh)/2d0*b1h-7d0/24d0*s+3d0/8d0*
     $     (xmw+3d0*xmh))/3d0
      b20z = b20z-(xmz+3d0*xmw-s)/12d0*dlog(xmw/mue2)-
     $     xmz/6d0*dlog(xmz/mue2)
      b20h = b20h-(xmw+3d0*xmh-s)/12d0*dlog(xmh/mue2)-
     $     xmw/6d0*dlog(xmw/mue2)
      a0z = xmz*(-dlog(xmz/mue2)+1d0)
      a0w = xmw*(-dlog(xmw/mue2)+1d0)
      a0h = xmh*(-dlog(xmh/mue2)+1d0)
      b0z = b0z-dlog(xmw*xmz/mue4)/2d0
      b0h = b0h-dlog(xmw*xmh/mue4)/2d0
      awboson = -cw2/sw2*(10d0*b20z-2d0*xmz-2d0*xmw+2d0/3d0*s+
     $     a0z+a0w+(xmw+xmz+4d0*s)*b0z)
      awboson = awboson+xmw*sw2/cw2*b0z
      awboson = awboson-b20z/sw2
      awboson = awboson+cw2/sw2*b20z*2d0
      awboson = awboson+(xmw*b0h-b20h)/sw2
      awboson = awboson+cw2/sw2*(3d0*a0z-2d0*xmz)
      awboson = awboson+(3d0*a0w-2d0*xmw)/sw2
      awboson = awboson+(a0h+a0z+2d0*a0w)/4d0/sw2
*
      if (s.eq.0d0) then
         dawboson = 0d0
         goto 200
      else
         call bquer1(s,mz,mw,p0,p1,db0z,db1z)
         call bquer1(s,mw,mh,p0,p1,db0h,db1h)
         call bquer(s,mz,mw,b0z,b1z)
         call bquer(s,mw,mh,b0h,b1h)
         db20z = (xmz*db0z+(s+xmz-xmw)/2d0*db1z+b1z/2d0-
     $        7d0/24d0)/3d0
         db20h = (xmw*db0h+(s+xmw-xmh)/2d0*db1h+b1h/2d0-
     $        7d0/24d0)/3d0
         db20z = db20z+dlog(xmw/mue2)/12d0
         db20h = db20h+dlog(xmh/mue2)/12d0
         b0z = b0z-dlog(xmw*xmz/mue4)/2d0
         dawboson = -cw2/sw2*(10d0*db20z+2d0/3d0+
     $        (4d0*s+xmw+xmz)*db0z+4d0*b0z)
         dawboson = dawboson+xmw*sw2/cw2*db0z
         dawboson = dawboson-db20z/sw2
         dawboson = dawboson+cw2/sw2*db20z*2d0
         dawboson = dawboson+(xmw*db0h-db20h)/sw2
      endif
 200  continue
*
      sigww = alpha/4d0/pi*(awlep+awquark+awboson)
      dsigww = alpha/4d0/pi*(dawlep+dawquark+dawboson)
*
      return
      end
c******************************************************************
      complex*16 function sigwwp(s)
c     W-selfenergy : (has been checked by comparison with Wolf)
c     only photonic part!

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      complex*16 awp
      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
c photon contribution to A_ww:
      dw = s-xmw
      if (dw.eq.0d0) then
         awp = dcmplx(19d0*dlog(xmw/mue2)-89d0/3d0,0d0)*xmw/3d0
      else
         if (s.eq.0d0) then
            awp = xmw*dcmplx(3d0*dlog(xmw/mue2)-2d0,0d0)
         else
            w = xmw/s
            w2 = w**2
            w3 = w2*w
            lw = dlog(dabs(dw)/xmw)
            awp = dcmplx((10d0/w+9d0)*dlog(xmw/mue2)-
     $           11d0-62d0/3d0/w+2d0*w+(10d0/w+2d0*w2-6d0*w-6d0)*lw,
     $           -pi*theta(dw)*(10d0/w+2d0*w2-6d0*w-6d0))*xmw/3d0
         endif
      endif
      sigwwp = alpha/4d0/pi*awp
*
      return
      end
c******************************************************************
      complex*16 function siggg(s)
c     photon-self energy (has been checked by comparison with Wolf)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 agferm,cff,cfw,agboson
      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      agferm = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cff = dcmplx(f(s,m,m),g(s,m,m))
         agferm = agferm+
     $        4d0/3d0*nc*q**2*(-s*dlog(m2/mue2)-s/3d0+
     $        (2d0*m2+s)*cff)
      end do
      cfw = dcmplx(f(s,mw,mw),g(s,mw,mw))
      agboson = 3d0*s*dlog(xmw/mue2)-(4d0*xmw+3d0*s)*cfw
*
      siggg = alpha/4d0/pi*(agferm+agboson)
*
      return
      end
c******************************************************************
      complex*16 function dsiggg(s)

c     photon-self energy (has been checked by comparison with Wolf)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 dagferm,cff,cfw,dagboson
      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      dagferm = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cff = dcmplx(f(s,m,m),g(s,m,m))
         dff = df(s,m,m)
         dagferm = dagferm+
     $        4d0/3d0*nc*q**2*(-dlog(m2/mue2)-1d0/3d0+
     $        (2d0*m2+s)*dff+cff)
      end do
      cfw = dcmplx(f(s,mw,mw),g(s,mw,mw))
      dfw = df(s,mw,mw)
      dagboson = 3d0*dlog(xmw/mue2)-(4d0*xmw+3d0*s)*dfw-3d0*cfw
*
      dsiggg = alpha/4d0/pi*(dagferm+dagboson)
*
      return
      end
c******************************************************************
      complex*16 function siggz(s)
c     Photon-Z-self energy (has been checked by comparison with Wolf)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 agzferm,agzboson,cff,cfw

      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      agzferm = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cff = dcmplx(f(s,m,m),g(s,m,m))
         agzferm = agzferm+4d0/3d0*nc*q*v*
     $        (s*dlog(m2/mue2)+s/3d0-(2d0*m2+s)*cff)
      end do
      cfw = dcmplx(f(s,mw,mw),g(s,mw,mw))
      agzboson = (-(2d0*xmw+(1d0/6d0+3d0*cw2)*s)*dlog(xmw/mue2)+
     $     ((1d0/6d0+3d0*cw2)*s+4d0*xmw*(1d0/3d0+cw2))*cfw+
     $     s/9d0)/cw/sw
*
      siggz = alpha/4d0/pi*(agzferm+agzboson)
*
      return
      end
c******************************************************************
c     derivative of the Photon-Z-self energy
      complex*16 function dsiggz(s)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 dagzferm,dagzboson,cff,cfw

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      dagzferm = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cff = dcmplx(f(s,m,m),g(s,m,m))
         dff = df(s,m,m)
         dagzferm = dagzferm+4d0/3d0*nc*q*v*
     $        (dlog(m2/mue2)+1d0/3d0-(2d0*m2+s)*dff-cff)
      end do
      cfw = dcmplx(f(s,mw,mw),g(s,mw,mw))
      dfw = df(s,mw,mw)
      dagzboson = (-(1d0/6d0+3d0*cw2)*dlog(xmw/mue2)+
     $     ((1d0/6d0+3d0*cw2)*s+4d0*xmw*(1d0/3d0+cw2))*dfw+
     $     (1d0/6d0+3d0*cw2)*cfw+1d0/9d0)/cw/sw
*
      dsiggz = alpha/4d0/pi*(dagzferm+dagzboson)
*
      return
      end
*------------------------------------------------------------------
      double precision function theta(xx)
c definition of the theta-function:
      implicit real*8(a-z)
      if (xx.ge.0d0) then
         theta = 1d0
      else
         theta = 0d0
      endif
      return
      end
c******************************************************************
      subroutine deltar(dr)
c     determination of delta_r:

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      complex*16 sww0,dsww0,swwos,dswwos,swwp0,swwpos,sigwwp

      pi = 4d0*datan(1d0)
      pi2 = pi**2
c     calculation of the renormalisation constant dz2w:
      call wrenorm(dz2w)

      call wselfw(xmw,swwos,dswwos)
      call wselfw(0d0,sww0,dsww0)
      swwp0 = sigwwp(0d0)
      swwpos = sigwwp(xmw)
      rsigw0 = dreal(sww0+swwp0)-dreal(swwos+swwpos)-xmw*dz2w
*
      dr = rsigw0/xmw+alpha/4d0/pi/sw2*(6d0+
     $     (7d0-4d0*sw2)/2d0/sw2*dlog(cw2))
*
      return
      end
*------------------------------------------------------------------
      real*8 function rho2(x)
      implicit real*8(a-z)
      data p1/-0.74141d0/,p2/ -11.483d0  /,p3/  9.6577d0/,
     $     p4/ -6.7270d0/,p5/  3.0659d0  /,p6/-0.82053d0/,
     $     p7/ 0.11659d0/,p8/-0.67712d-02/
      pi=4d0*datan(1d0)
      pi2=pi*pi
      if(x.le.4d0) then
         rho2=p1+p2*x+p3*x**2+p4*x**3+p5*x**4+p6*x**5+p7*x**6+p8*x**7
      else
         rbth=1/x**2
         alrb=log(rbth)
         rho2=49d0/4d0+pi2+27d0/2d0*alrb+3d0/2d0*alrb**2
     $        +rbth/3d0*(2d0-12d0*pi2+12d0*alrb-27d0*alrb**2)
     $        +rbth**2/48d0*(1613d0-240d0*pi2-1500d0*alrb
     $        -720d0*alrb**2)
      endif
      return
      end
c******************************************************************
      subroutine mwiter(mwerror)
c     determination of mw from deltar:

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ii,mwerror

      common/drit/dr 

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      mwold = mw
c     calculating mw(dr(mw)):
      ashort = pi*alpha/w2/gfermi
      do ii = 0,12
         call deltar(dr)
         mw = mz*dsqrt(1d0+
     $        dsqrt(1d0-4d0*ashort/(1d0-dr)/xmz))/w2
         xmw = mw**2
         cw2 = xmw/xmz
         cw = dsqrt(cw2)
         sw2 = 1d0-cw2
         sw = dsqrt(sw2)
         if (dabs(mw-mwold).le.1d-5) goto 21
         mwold = mw
      end do
      mwerror = 1
 21   continue
      return
      end
c******************************************************************
c calculation of the renormalisation constant dz2w
      subroutine wrenorm(dz2w)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff,i,mix,hqcd,twol,had
      complex*16 swwos,dswwos,swwpos,sigzos,sigzz
     $     ,sgz,sgzz,rsgz,rsgzz,itszz,hold,sigwwp,siggg,siggz
      complex*16 dgzmix,dsiggg,dsiggz,dsgz,dsgzz

      common/renorm/mue,mue2,mue4
      common/gzmix/mix
      common/higher/hqcd,twol
      common/renself/rsgzz,rsgz
      common/derivmix/dgzmix
      common/dalpha/dhad,dferm,disint
      common/jeger/had
      common/inv/alpzinv

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      pig0f = 0d0
      do ff=4,12
          call prop(ff,m,q,v,a,nc)
          pig0f = pig0f-alpha/4d0/pi*
     $        4d0/3d0*nc*q**2*dlog(m**2/mue2)
      end do
*
      disint = 1d0-alpha*alpzinv-dferm-dhad
c calculation of the hadronic contribution to pig0
c by using the dispersion integral result:
      if (had.eq.1) pig0f = pig0f+dhad+disint
*
      pig0 = pig0f+alpha/4d0/pi*(3d0*dlog(xmw/mue2)-2d0/3d0)
      sgz0 = -alpha/4d0/pi*2d0*xmw/sw/cw*dlog(xmw/mue2)
      sigzos = sigzz(xmz)
*
c weak contribution:
      call wselfw(xmw,swwos,dswwos)
c photonic contribution:
      swwpos = sigwwp(xmw)
*
      dzww = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
*
c inclusion of higher order terms (QCD and two-loop):
c qcd corrections :
      kf = 12d0*pi/23d0
      alfst = kf/(2d0*dlog(mtop/mz)+kf/alphas)
      qcd = alfst/pi*(pi2/3d0+1d0)*gfermi*mtop**2/(4d0*pi2*w2)
c add higher order qcd corrections (taken from Wolf) :
      qcd = alfst/pi*2d0*(pi**2/3d0+1d0)
     $     *gfermi*mtop**2/(8d0*pi**2*w2)
     $     +(alfst/pi)**2*3d0*pi**2*(2.155165d0-5d0*0.180981d0)
     $     *gfermi*mtop**2/(8d0*pi**2*w2)
      if (hqcd.eq.0) qcd = 0d0
c two-loop-corrections :
      xt = gfermi*mtop**2/8d0/pi2/w2
      tl = -3d0*xt**2*rho2(mh/mtop)
      if (twol.eq.0) tl = 0d0
*
      dzww = dzww-qcd-tl
      dz2w = -pig0-2d0*cw/sw/xmz*sgz0+cw2/sw2*dzww
*
c inclusion of the shift in the renormalisation condition for
c delta M_Z due to photon-Z-mixing:
*
      dzwit = dzww
      sgz = siggg(xmz)
      sgzz = siggz(xmz)
      rsgz = sgz-xmz*pig0
      rsgzz = sgzz+sgz0-cw/sw*dzwit*xmz
*
      dsgzz = dsiggz(xmz)
      dsgz = dsiggg(xmz)
      dgzmix = -2d0*rsgzz*(dsgzz+2d0*sgz0/xmz-cw/sw*(dzwit))
     $     /(xmz+rsgz)
     $     +(1d0+dsgz-pig0)*rsgzz**2/(xmz+rsgz)**2
*
c      write(6,*)1d0/alpha,rsgzz,rsgz
      if (mix.eq.0) goto 99
*
      itszz = sigzos
      hold = rsgzz
      do i=1,12
         sigzos = itszz-rsgzz**2/(xmz+rsgz)
         dzwit = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
         dzwit = dzwit-qcd-tl
         rsgzz = sgzz+sgz0-cw/sw*dzwit*xmz
         if (cdabs(rsgzz-hold).lt.1d-5) goto 999
         hold = rsgzz
      end do
 999  continue
      sigzos = itszz-rsgzz**2/(xmz+rsgz)
      dzwit = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
      dzwit = dzwit-qcd-tl
*
      dz2w = -pig0-2d0*cw/sw/xmz*sgz0+cw2/sw2*dzwit
*
      dgzmix = -2d0*rsgzz*(dsgzz+2d0*sgz0/xmz-cw/sw*(dzwit))
     $     /(xmz+rsgz)
     $     +(1d0+dsgz-pig0)*rsgzz**2/(xmz+rsgz)**2
*
 99   continue
*
      return
      end
c******************************************************************
      subroutine gwffp(f,fp,gammatree,gammaloop)
c     calculation of the fermionic partial W-width W->f fp: tree level,
c     1-loop EW and higher order QCD; 
c     gammaloop = gammatree*(1+corr.)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer f,fp,hqcd,twol,mix
      real*8 ckm2(7:8,10:12)
      complex*16 fvwffp,swwos,dswwos,dawweak,dweak
      common/renorm/mue,mue2,mue4
      common/kobaya/ckm2
      common/drit/dr
      common/split/g0,dqed,dweak,dqcd
      common/gzmix/mix
      common/higher/hqcd,twol

      pi = 4d0*datan(1d0)
      pi2 = pi**2
c choice of the representation: rep=1: alpha-repr. and
c rep=2: gmu-repr.:
      if (rep.eq.1) then
         g0 = alpha*mw/sw2/12d0
      else 
         g0 = w2*gfermi/pi*xmw*mw/12d0
      end if
      call prop(f,m,q,v,a,nc)
      call prop(fp,mp,qp,vp,ap,ncp)
      m2 = m**2
      mp2 = mp**2
      kap = dsqrt((xmw-m2-mp2)**2-4d0*m2*mp2)/xmw*
     &     (1d0-(m2+mp2)/2d0/xmw-(m2-mp2)**2/2d0/xmw**2)
      if (mzero.eq.1) kap = 1d0
      gammatree = g0*kap*nc
      if (f.ge.7) gammatree = gammatree*ckm2(f,fp)
c gauge invariant separation of the complete 1-loop corrections 
c into a QED-like and a modified weak part:
c QED-like contribution (virt.+real(soft+hard)):
      dqed = alpha/pi*(3d0/8d0*(q**2+qp**2)+7d0/3d0+pi2/24d0)
c pure weak contribution:
      call wselfw(xmw,swwos,dswwos)
      call fvpurew(xmw,v,a,vp,ap,fvwffp)
c remainder of the photon contribution:
      dgamrem = alpha/4d0/pi*(-25d0/3d0*dlog(xmw/mue2)+68d0/9d0-
     $     3d0/2d0*pi2)-0d0*dz2wp
c wavefunction renorm.:
      call wrenorm(dz2w)
*
c      dweak = 2d0*dreal(fvwffp-dawweak/2d0)+dgamrem
      dweak = 2d0*(fvwffp-dswwos/2d0)+dgamrem-dz2w
*
      if (rep.eq.1) then
         gammaloop = gammatree*(1d0+dqed+dreal(dweak))
      else
         gammaloop = gammatree*(1d0+dqed+dreal(dweak)-dr)
      end if
c QCD correction:
      dqcd = 0d0
      if ((qcd.eq.1).and.(f.ge.7)) then 
c         dqcd = alphas/pi*(1d0+1.40932d0*alphas/pi-
c     $        12.76706d0*(alphas/pi)**2)
         dqcd = alphas/pi*(1d0+1.405d0*alphas/pi-
     $        12.8d0*(alphas/pi)**2)
c         dqcd=alphas/pi
         gammaloop = gammaloop+gammatree*dqcd
      end if
*
      return
      end
c******************************************************************
      double precision function hadvac(s)
c     calcof the hadronic vacuum polarisation for 5 active flavours

      implicit real*8(a-z)
      integer ff

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      pig0f = 0d0
      agferm = 0d0
      do ff=7,12
         if (ff.eq.9) goto 99
          call prop(ff,m,q,v,a,nc)
          pig0f = pig0f-alpha/4d0/pi*
     $        4d0/3d0*nc*q**2*dlog(m**2/mue2)
          agferm = agferm+
     $         4d0/3d0*nc*q**2*(-s*dlog(m**2/mue2)-s/3d0+
     $         (2d0*m**2+s)*f(s,m,m))
 99       continue
      end do
      siggg = alpha/4d0/pi*agferm
*
      hadvac = (siggg-s*pig0f)/s
*
      return
      end
c******************************************************************
c     calculation of the renormalisation constant dz2z:
      subroutine zrenorm(dz2z)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff,i,mix,hqcd,twol,had
      complex*16 sigzos,swwos,swwpos,dswwos,sigwwp,sigzz
     &     ,sgz,sgzz,rsgz,rsgzz,itszz,hold,siggg,siggz

      common/renorm/mue,mue2,mue4
      common/higher/hqcd,twol
      common/gzmix/mix
      common/dalpha/dhad,dferm,disint
      common/jeger/had
      common/inv/alpzinv

      pi = 4d0*datan(1d0)
      pi2 = pi**2

      pig0f = 0d0
      do ff=4,12
          call prop(ff,m,q,v,a,nc)
          pig0f = pig0f-alpha/4d0/pi*
     &        4d0/3d0*nc*q**2*dlog(m**2/mue2)
      end do

      disint = 1d0-alpha*alpzinv-dferm-dhad
c calculation of the hadronic contribution to pig0
c by using the dispersion integral result:
      if (had.eq.1) pig0f = pig0f+dhad+disint

      pig0 = pig0f+alpha/4d0/pi*(3d0*dlog(xmw/mue2)-2d0/3d0)
      sgz0 = -alpha/4d0/pi*2d0*xmw/sw/cw*dlog(xmw/mue2)
      sigzos = sigzz(xmz)
c weak contribution:
      call wselfw(xmw,swwos,dswwos)
c photonic contribution:
      swwpos = sigwwp(xmw)
*
      dzww = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
*
c inclusion of higher order terms (QCD and two-loop):
c qcd corrections :
      kf = 12d0*pi/23d0
      alfst = kf/(2d0*dlog(mtop/mz)+kf/alphas)
      qcd = alfst/pi*(pi2/3d0+1d0)*gfermi*mtop**2/(4d0*pi2*w2)
c add higher order qcd corrections (taken from Wolf) :
      qcd = alfst/pi*2d0*(pi**2/3d0+1d0)
     $     *gfermi*mtop**2/(8d0*pi**2*w2)
     $     +(alfst/pi)**2*3d0*pi**2*(2.155165d0-5d0*0.180981d0)
     $     *gfermi*mtop**2/(8d0*pi**2*w2)
      if (hqcd.eq.0) qcd = 0d0
c two-loop-corrections :
      xt = gfermi*mtop**2/8d0/pi2/w2
      tl = -3d0*xt**2*rho2(mh/mtop)
      if (twol.eq.0) tl = 0d0
*
      dzww = dzww-qcd-tl
*
      dz2z = -pig0-2d0*(cw2-sw2)/cw/sw/xmz*sgz0+(cw2-sw2)/sw2*dzww
*
c inclusion of the shift in the renormalisation condition for
c delta M_Z due to photon-Z-mixing:
*
      if (mix.eq.0) goto 99
*
      dzwit = dzww
      sgz = siggg(xmz)
      sgzz = siggz(xmz)
      rsgz = sgz-xmz*pig0
      rsgzz = sgzz+sgz0-cw/sw*dzwit*xmz
*
      itszz = sigzos
      hold = rsgzz
      do i=1,12
         sigzos = itszz-rsgzz**2/(xmz+rsgz)
         dzwit = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
         dzwit = dzwit-qcd-tl
         rsgzz = sgzz+sgz0-cw/sw*dzwit*xmz
         if (cdabs(rsgzz-hold).lt.1d-5) goto 999
         hold = rsgzz
      end do
 999  continue
      sigzos = itszz-rsgzz**2/(xmz+rsgz)
      dzwit = dreal(sigzos)/xmz-dreal(swwos+swwpos)/xmw
      dzwit = dzwit-qcd-tl
*
      dz2z = -pig0-2d0*(cw2-sw2)/cw/sw/xmz*sgz0+(cw2-sw2)/sw2*dzwit
 99   continue
*
      return
      end
c******************************************************************
      double precision function dsigzz(s)
c     derivative of the Z-self energy (only real part):

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 b0h,b1h,b0w,b1w

      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      dazferm = 0d0
      dazneut = 0d0
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         dcf = df(s,m,m)
         cf = f(s,m,m)
         dazferm = dazferm+
     $        4d0/3d0*nc*((v**2+a**2)*(-dlog(m2/mue2)+
     $        (2d0*m2+s)*dcf-1d0/3d0+cf)-a**2*6d0*m2*dcf)
      end do
      do ff=1,3
         call prop(ff,m,q,v,a,nc)
         dazneut = dazneut+
     $        8d0/3d0*nc*a**2*(-dlog(s/mue2)+5d0/3d0-1d0)
      end do
      call bquer(s,mw,mw,b0w,b1w)
      call bquer1(s,mw,mw,p0,p1,db0w,db1w)
      db20w = (xmw*db0w+s/2d0*db1w+dreal(b1w)/2d0-7d0/24d0)/3d0
      db20w = db20w+1d0/12d0*dlog(xmw/mue2)
      b0w = b0w-dlog(xmw/mue2)
      dazboson = -cw2/sw2*(10d0*db20w+
     $     2d0*(xmw+2d0*s)*db0w+4d0*dreal(b0w)+2d0/3d0)
      dazboson = dazboson+xmw*sw2/cw2*2d0*db0w
      dazboson = dazboson-(sw2-cw2)**2/2d0/sw2/cw2*2d0*db20w
      dazboson = dazboson+cw2/sw2*db20w*2d0
      call bquer(s,mz,mh,b0h,b1h)
      call bquer1(s,mz,mh,p0,p1,db0h,db1h)
      db20h = (xmz*db0h+(s+xmz-xmh)/2d0*db1h+
     &     dreal(b1h)/2d0-7d0/24d0)/3d0
      db20h = db20h+1d0/12d0*dlog(xmh/mue2)
      dazboson = dazboson-(db20h-xmz*db0h)/cw2/sw2

      dsigzz = alpha/4d0/pi*(dazferm+dazneut+dazboson)
      return
      end
c******************************************************************
      real*8 function a2qqw2_check(xs,xt,xu,xvf,xaf,xvfs,xafs,
     &     xvi,xai,xvis,xais,xqi,xqis,xqfs,xmi,xmis,xmfs,sigborn)

c     full EW 1-loop corrections to the 
c     W-exchange process:  f=nu(p4)+i(p3)->f'(p2)+i'(p1).
c     from nutev code

      implicit none
      include 'pwhg_wzgrad.h'
      integer i,iqqbar
      real*8 boxq,boxt_wz,boxt_zw,boxu_wz,boxu_zw
      real*8 boxt2_v1v2,boxu2_v1v2,sigborn,sigq
      real*8 mi,mis,mfs
      real*8 tp,up,sp,xs,xt,xu
      real*8 xvf,xaf,xvi,xai,kap1,kap2
      real*8 xvfs,xafs,xvis,xais,qqis,qqi,qqfs
      real*8 xqi,xqis,xqfs,xmi,xmis,xmfs
      real*8 fvweak,fvg,denw,dz2w
      complex*16 swwos,sww,dswwos,dsww,fvwff,fvwii
      complex*16 sigwwp,swwp,swwpos,sigzz,szzos
      complex*16 fvgff,fvgii,fvbox
      complex*16 spence,zero
      real*8 den
      common/par_prop/den
      real*8 mue,mue2,mue4
      common/renorm/mue,mue2,mue4
      real*8 deltar
      common/drit/deltar 
      real*8 pig0,sgz0
      common/sigzero/pig0,sgz0
      real*8 pi,pi2

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      a2qqw2_check=0d0
      sigq=0d0
      boxq=0d0
      sp=xs
      tp=xt
      up=xu
      qqis=xqis
      qqi=xqi
      qqfs=xqfs
      mis=xmis
      mi=xmi
      mfs=xmfs
      zero=dcmplx(0d0,0d0)
c-------------------------------------------------------------
c matrix element squared with weak vertex and self energy corrections
c to W production: 
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
c pure weak contribution evaluated at that:
c check of mue-dependence
c      do i=1,3
c         mue=1d1**(i**2)
c         mue2=mue*mue
c         mue4=mue2*mue2
c weak part to W vertex corrections
c weak part to W self energy
c W propagator
      denw=tp-xmw
      call wselfw(xmw,swwos,dswwos)
      call wselfw(tp,sww,dsww)
      call fvpurew(tp,xvf,xaf,xvfs,xafs,fvwff)
      call fvpurew(tp,xvis,xais,xvi,xai,fvwii)
      fvweak = dreal(fvwff+fvwii-(sww-swwos)/denw)
c pure photonic vertex contribution:
c photonic part to W self energy is included in form factors
      call fvpureg(tp,0d0,qqfs,1d-5,mfs,fvgff)      
      call fvpureg(tp,qqis,qqi,mis,mi,fvgii)
      fvg = dreal(fvgff+fvgii)
c photonic box contribution:
      call fvboxwg_new(iqqbar,tp,sp,qqi,qqis,qqfs,mi,mis,mfs,fvbox)
c add the W self energy CT (photonic+weak part):
      call wrenorm(dz2w)
      sigq=2d0*(fvg+fvweak-dz2w-(rep-1)*deltar)*sigborn+
     $     4d0*denw*dreal(fvbox)/8d0/sw2/sw2
*
c      write(6,*)mue,fvg+fvweak-dz2w
c      enddo
c      stop
c      if (mw.gt.0d0) goto 9999
c-------------------------------------------------------------
c matrix element squared for W,Z boxes to W production:
c (all fermions are considered to be massless !)
c-------------------------------------------------------------
c WZ box:
      kap1=(xvfs+xafs)*(xvis+xais)/8d0/sw2/sw2
      kap2=0d0
      if(iqqbar.eq.1)
     $     boxt_wz=boxt2_v1v2(tp,sp,mw,mz,kap1,kap2)
      if(iqqbar.eq.2)
     $     boxt_wz=boxt2_v1v2(tp,sp,mw,mz,kap2,kap1)
c ZW box:
      kap1=(xvf+xaf)*(xvi+xai)/8d0/sw2/sw2
      kap2=0d0
      if(iqqbar.eq.1)
     $     boxt_zw=boxt2_v1v2(tp,sp,mz,mw,kap1,kap2)
      if(iqqbar.eq.2)
     $     boxt_zw=boxt2_v1v2(tp,sp,mz,mw,kap2,kap1)
c WZ box (crossed):
      kap1=0d0
      kap2=(xvfs+xafs)*(xvi+xai)/8d0/sw2/sw2      
      if(iqqbar.eq.1)
     $     boxu_wz=boxu2_v1v2(up,sp,mw,mz,kap1,kap2)
      if(iqqbar.eq.2)
     $     boxu_wz=boxu2_v1v2(up,sp,mw,mz,kap2,kap1)
c ZW box (crossed):
      kap1=0d0
      kap2=(xvf+xaf)*(xvis+xais)/8d0/sw2/sw2
      if(iqqbar.eq.1)
     $     boxu_zw=boxu2_v1v2(up,sp,mz,mw,kap1,kap2)
      if(iqqbar.eq.2)
     $     boxu_zw=boxu2_v1v2(up,sp,mz,mw,kap2,kap1)
c sum of box diagrams
      boxq=alpha0/4d0/pi*32d0*denw*
     $     (boxt_wz+boxt_zw+boxu_wz+boxu_zw)
c sum:
 9999 continue
      a2qqw2_check=sigq+boxq
      return
      end
*-----------------------------------------------------------------
      double precision function boxt2_v1v2(xt,xs,xmv1,xmv2,xkap1,xkap2)
*-----------------------------------------------------------------
c t-channel box contribution with two vector bosons v1,v2 and
c two fermions. All fermions are considered to be massless !
*-----------------------------------------------------------------
      implicit none
      real*8 xt,xs,xu,xmv1,xmv2,xkap1,xkap2
      real*8 t,u,s,mv1,mv2
      complex*16 d0,d11,d12,d13,d20,d21,d22,d23,d212,d213,d223
      complex*16 c0,C0_
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      s=xs
      t=xt
      u=-s-t
      mv1=xmv1
      mv2=xmv2
      call ddfunc(s,0d0,0d0,s/2d0,s/2d0,-t/2d0,
     $     0d0,0d0,mv2,mv1,d1,d2,d3,c1,c2)
c      c0=C0_(0d0,t,0d0,0d0,mv2,mv1,0)                            
      c0 = c1(0,1)
      d0 = d1(0)
      d11 = d1(1)
      d12 = d1(2)
      d13 = d1(3)
      d20 = d2(0)
      d21 = d2(1)
      d22 = d2(2)
      d23 = d2(3)
      d212 = d2(4)
      d213 = d2(5)
      d223 = d2(6)

      boxt2_v1v2=16d0*dreal(xkap1*(-s**3*(2d0*d11+d12+d13)-2d0*s**2*c0)+
     $     xkap2*(-u**3*d223-u*s**2*(d21+d212+d213+d223)+
     $     u**2*s*(d11+d12+d13-2d0*d223)+u**2*(c0-2d0*d20)+
     $     s*u*(c0-4d0*d20)))
      return
      end
*-----------------------------------------------------------------
      double precision function boxu2_v1v2(xu,xs,xmv1,xmv2,xkap1,xkap2)
*-----------------------------------------------------------------
c u-channel box contribution with two vector bosons v1,v2 and
c two fermions. All fermions are considered to be massless !
*-----------------------------------------------------------------
      implicit none
      real*8 xu,xs,xmv1,xmv2,xkap1,xkap2
      real*8 t,u,s,mv1,mv2
      complex*16 d0,d11,d12,d13,d20,d21,d22,d23,d212,d213,d223
      complex*16 c0,C0_
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      s=xs
      u=xu
      t=-s-u
      mv1=xmv1
      mv2=xmv2
      call ddfunc(u,0d0,0d0,u/2d0,u/2d0,-t/2d0,
     $     0d0,0d0,mv2,mv1,d1,d2,d3,c1,c2)
c      c0=C0_(0d0,u,0d0,0d0,mv2,mv1,0)                            
      c0 = c1(0,1)
      d0 = d1(0)
      d11 = d1(1)
      d12 = d1(2)
      d13 = d1(3)
      d20 = d2(0)
      d21 = d2(1)
      d22 = d2(2)
      d23 = d2(3)
      d212 = d2(4)
      d213 = d2(5)
      d223 = d2(6)
      boxu2_v1v2=16d0*dreal(xkap1*(-u**3*(2d0*d11+d12+d13)-2d0*u**2*c0)+
     $     xkap2*(-s**3*d223-s*u**2*(d21+d212+d213+d223)+
     $     s**2*u*(d11+d12+d13-2d0*d223)+s**2*(c0-2d0*d20)+
     $     s*u*(c0-4d0*d20)))
      return
      end
c******************************************************************
      subroutine formgz2(xs,xv,xa,xq,xvs,xas,xqs,fvz,fvg,gaz,gag)

c     formfactors describing the weak contribution to the (Z,gamma)
c     ffbar vertex correction.
c     vertex: i lamda_mu(s) = i*e*g_mu*(fv(s)-ga(s) g5)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      complex*16 b0,b1,c0,C0_,spence,cscal
      complex*16 factor,fvz,fvg,gaz,gag
      complex*16 vfz,vfw,vfws,vfww,sc,wc

      common/renorm/mue,mue2,mue4
      common/sigzero/pig0,sgz0

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      s=xs
      sc=s+ieps
      v=xv
      a=xa
      q=xq
      vs=xvs
      as=xas
      qs=xqs
c
      fvz=dcmplx(0d0,0d0)
      fvg=dcmplx(0d0,0d0)
      gaz=dcmplx(0d0,0d0)
      gag=dcmplx(0d0,0d0)
c z-exchange:
      w=xmz/s
      wc=xmz/sc
c c0 is already multiplied with s:
c      c0=dcmplx(dreal(spence(-1d0/w+ieps))-dlog(dabs(w))*
c     $     dlog(1d0+1d0/w),-pi*dlog(1d0+1d0/w))
      c0=-(spence(1d0+s/(xmz+ieps))-pi2/6d0)
      factor=-dlog(xmz/mue2)-(2d0*w+3d0)*cdlog(-wc)-2d0*w-
     $     4d0-2d0*(1d0+w)**2*c0
      vfz=alpha0/4d0/pi*factor
c self energy of the external fermions:
      vfz=vfz+alpha0/4d0/pi*dcmplx(dlog(xmz/mue2)+1d0/2d0,0d0)
c w-exchange:
      w=xmw/s
      wc=xmw/sc
c c0 is already multiplied with s:
c      c0=dcmplx(dreal(spence(-1d0/w+ieps))-dlog(dabs(w))*
c     $     dlog(1d0+1d0/w),-pi*dlog(1d0+1d0/w))
      c0=-(spence(1d0+s/(xmw+ieps))-pi2/6d0)
      factor=-dlog(xmw/mue2)-(2d0*w+3d0)*cdlog(-wc)-2d0*w-
     $     4d0-2d0*(1d0+w)**2*c0
      vfw=alpha0/4d0/pi*factor      
c self energy of the external fermions:
      vfws=alpha0/4d0/pi*dcmplx((dlog(xmw/mue2)+1d0/2d0),0d0)
c
c z(g)ww-vertex:
c      c0=cscal(s,1d-1,mw,mw,0d0)
      c0=C0_(0d0,s,0d0,0d0,mw,mw,0)
      call bquer(s,mw,mw,b0,b1)
      factor=-3d0*dlog(xmw/mue2)+2d0*w+4d0-(2d0*w+1d0)*b0+
     $     2d0*xmw*(2d0+w)*c0
      vfww=alpha0/4d0/pi*factor
      t3f=2d0*cw*sw*a
      fvz=(v*(v**2+a**2)+a*2d0*v*a)*vfz+
     $     (vs+as)/4d0/sw**2*vfw+(v+a)/4d0/sw**2*vfws-
     $     (-2d0*t3f)*cw/sw/4d0/sw**2*vfww-cw/sw*a*sgz0/xmz
      fvg=-q*(v**2+a**2)*vfz-qs/4d0/sw**2*vfw+(-2d0*t3f)*
     $     1d0/4d0/sw**2*vfww-q/4d0/sw**2*vfws+a*sgz0/xmz
      gaz=(a*(v**2+a**2)+v*2d0*v*a)*vfz+
     $     (vs+as)/4d0/sw**2*vfw+(v+a)/4d0/sw**2*vfws-
     $     (-2d0*t3f)*cw/sw/4d0/sw**2*vfww-cw/sw*a*sgz0/xmz
      gag=-q*2d0*v*a*vfz-qs/4d0/sw**2*vfw-q/4d0/sw**2*vfws+
     $     (-2d0*t3f)/4d0/sw**2*vfww+a*sgz0/xmz
*
      return
      end
c******************************************************************
      subroutine formzqed(xt,xq,xm,form)

c     pure photonic contribution:
c     formfactors fvin(i)=fain(i) and fvfin(i) describing the 
c     ith contribution to the initial state vertex correction.
c     vertex: i e lamda_mu(s) = i*e*g_mu*(vi-ai g5)*fvertex(s)
c     no contribution to the final state vertex for nuN scattering.

      implicit none
      include 'pwhg_wzgrad.h'
      integer i
      real*8 xt,t,xq,q,xm,m
      real*8 pi,pi2
      real*8 mue,mue2,mue4
      complex*16 b0,c0,C0_,b0m,b0ms,b0s,spence
      complex*16 fv(2),factor,form,cdlog,tc

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)  
      pi2 = pi**2
      q=xq
      m=xm
      t=xt
      tc=t+ieps
*
c g-exchange:
c c0 is already multiplied with t:
c      c0 = dcmplx(-cdlog(tc/m/ms)*cdlog(lambda**2/tc)-
c     $     cdlog(tc/m**2)**2/4d0-
c     $     cdlog(tc/ms**2)**2/4d0-2d0/3d0*pi2,
c     $     pi*dlog(lambda**2/tc))
      c0=-(dlog(-m**2/t)*dlog(m**2/lambda**2)-dlog(-m**2/t)**2/2d0+
     $     pi2/6d0)
      b0m=dcmplx(-dlog(m**2/mue2)+2d0,0d0)
      b0s=-dlog(-t/mue2)+2d0
      factor = -2d0*c0+4d0*b0m-3d0*b0s-2d0
      fv(1) = alpha0/4d0/pi*q**2*factor
c      write(6,*)'photon exch.:',s,c0/s,fv(1)
c      stop
*
c self energy of the external fermions:
c g:
      factor=-1d0/2d0*(-dlog(m**2/mue2)+4d0+
     $     2d0*dlog(lambda**2/m**2))
      fv(2)=2d0*alpha0/4d0/pi*q**2*factor
*
c pure photonic contribution to the formfactors:
      form = dcmplx(0d0,0d0)
      do i=1,2
         form = form+fv(i)
      end do
      return
      end
c******************************************************************
      subroutine fvpurew2(xt,xv,xa,xvs,xas,form)

c     weak contribution:
c     formfactors fvin(i)=fain(i) and fvfin(i) describing the 
c     ith contribution to the initial or final state vertex correction.
c     vertex: i lamda_mu(s) = i*e/2/w2/sw*g_mu*(1-g5)*fvertex(s)

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      complex*16 b0,b1,c0,C0_,spence,cscal
      complex*16 fv(5),factor,form,tc,wc
      complex*16 c11,c12,c21,c22,c20,c23,
     &     b01,b11,b20

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      t=xt
      tc=t+ieps
      v=xv
      a=xa
      vs=xvs
      as=xas
c z-exchange:      
      w=xmz/t
      wc=xmz/tc
c c0 is already multiplied with t:
c      c0 = dcmplx(dreal(spence(-1d0/w+ieps))-
c     $     dlog(w)*dlog(1d0+1d0/w),
c     $     pi*dlog(1d0+1d0/w))
      c0=-(spence(1d0+t/(xmz+ieps))-pi2/6d0)
      factor = -dlog(xmz/mue2)-(2d0*w+3d0)*dlog(-w)-2d0*w-
     $     4d0-2d0*(1d0+w)**2*c0
c      write(6,*)c0/t,factor
c alternative calculation:
c      b01=-cdlog(tc/mue2)+2d0+dcmplx(0d0,1d0)*pi
c      call ccfunc2(0d0,0d0,-t/2d0,mz,0d0,0d0,
c     $     c0,c11,c12,c21,c22,c20,c23)
c      factor=b01-2d0-xmz*(c11+c12)-2d0*s*(c0+c11+c12)
c
      fv(1) = alpha0/4d0/pi*(v+a)*(vs+as)*factor
*
c wzw-vertex:
c      c0 = cscal(t,1d-1,mz,mw,0d0)
      c0=C0_(0d0,t,0d0,0d0,mz,mw,0)
      call bquer(t,mz,mw,b0,b1)
      factor = -3d0/2d0*dlog(xmz*xmw/mue4)-1d0/2d0*(xmw-xmz)/t*
     $     dlog(xmw/xmz)+4d0+(xmw+xmz)/t-((xmw+xmz)/t+1d0)*b0+
     $     2d0*(xmz+xmw+xmz*xmw/t)*c0
      fv(2) = alpha0/4d0/pi*cw/sw*(v+a-vs-as)*factor
c      write(6,*)'WZW',c0,factor
c check - ok
c      call ccfunc2(0d0,0d0,-t/2d0,0d0,mz,mw,
c     $     c0,c11,c12,c21,c22,c20,c23)
c      write(6,*)'WZW',c0,2d0*(b0+2d0*c20-0.5d0-t*(c11+c12)-
c     $     1d0/2d0*dlog(xmz*xmw/mue4))
*
c self energy of the external fermions:
c Z:
      fv(3) = alpha0/4d0/pi*((v+a)**2+(vs+as)**2)*
     $     dcmplx((dlog(xmz/mue2)+1d0/2d0)/2d0,0d0)
c W:
      fv(4) = alpha0/4d0/pi*dcmplx((dlog(xmw/mue2)+
     $     1d0/2d0)/2d0/sw2,0d0)
*
c counterterm: dz1_w-dz2_w:
      fv(5) = alpha0/4d0/pi*dcmplx(2d0/sw2*dlog(xmw/mue2),0d0)
*
c pure weak contribution to the formfactors:
      form = dcmplx(0d0,0d0)
      do i=1,5
         form = form+fv(i)
      end do
*
      return
      end
c******************************************************************
      subroutine fvpureg2(xt,xq,xqs,xm,xms,form)

c     pure photonic contribution:
c     formfactors fvin(i)=fain(i) and fvfin(i) describing the 
c     ith contribution to the initial or final state vertex correction.
c     vertex: i e lamda_mu(s) = i*e/2/w2/sw*g_mu*(1-g5)*fvertex(s)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer i
      complex*16 b0,c0,C0_,b0m,b0ms,b0s,spence
      complex*16 fv(5),factor,form,cdlog,tc

      common/renorm/mue,mue2,mue4

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      q=xq
      qs=xqs
      m=xm
      ms=xms
      t=xt
      tc=t+ieps
      lwt=(t-xmw)*dlog(1d0-t/xmw)
c g-exchange:
c c0 is already multiplied with t:
      c0=-(dlog(-m*ms/t)*dlog(m*ms/lambda**2)-dlog(-m*ms/t)**2/2d0+
     $     dlog(m/ms)**2/2d0+pi2/6d0)
      b0m=dcmplx(-dlog(m**2/mue2)+2d0,0d0)
      b0ms=dcmplx(-dlog(ms**2/mue2)+2d0,0d0)
      b0s=-dlog(-t/mue2)+2d0
      factor=-2d0*c0+2d0*b0m+2d0*b0ms-3d0*b0s-2d0
      fv(1)=alpha0/4d0/pi*q*qs*factor
c      write(6,*)'photon exch.:',t,c0/t,fv(1)
c      stop
*
c wgw-vertex:
c collinear singularity has been subtracted and added to factor:
      c0 = -(spence(1d0-xmw/t+ieps)+
     $     cdlog(tc/xmw)**2/2d0+pi2/6d0)/t
c check (add back the coll. singularity)
c      write(6,*)c0-cdlog(m**2/tc)*dlog(1d0-t/xmw)/t
c      c0=C0_(m**2,0d0,t,0d0,m,mw,0)
c      write(6,*)'c0',c0
c      stop
      b0=-dlog(xmw/mue2)+2d0-lwt/t
      factor=q*(2d0*xmw*c0+2d0*b0m+(2d0+xmw/t)*(-dlog(xmw/mue2)+1d0)-
     $     (1d0+xmw/t)*b0+2d0*xmw/t*cdlog(tc/m**2)*dlog(1d0-t/xmw))-
     $     qs*(2d0*xmw*c0+2d0*b0ms+(2d0+xmw/t)*(-dlog(xmw/mue2)+1d0)-
     $     (1d0+xmw/t)*b0+2d0*xmw/t*cdlog(tc/ms**2)*dlog(1d0-t/xmw))
      fv(2)=alpha0/4d0/pi*factor
c      write(6,*)'wgw vertex:',q,alpha0/4d0/pi*
c     $     dsqrt(4d0*pi*alpha0)/2d0/w2/sw*
c     $     q*(2d0*xmw*c0+2d0*b0m+(2d0+xmw/t)*(-dlog(xmw/mue2)+1d0)-
c     $     (1d0+xmw/t)*b0+2d0*xmw/t*cdlog(tc/m**2)*dlog(1d0-t/xmw))
c      write(6,*)'wgw vertex:',qs,alpha0/4d0/pi*
c     $     dsqrt(4d0*pi*alpha0)/2d0/w2/sw*
c     $     (-qs)*(2d0*xmw*c0+2d0*b0ms+(2d0+xmw/t)*(-dlog(xmw/mue2)+1d0)-
c     $     (1d0+xmw/t)*b0+2d0*xmw/t*cdlog(tc/ms**2)*dlog(1d0-t/xmw))
*
c self energy of the external fermions:
c g:
      fv(3) = -1d0/2d0*alpha0/4d0/pi*q**2*
     $     (-dlog(m**2/mue2)+4d0+2d0*dlog(lambda**2/m**2))
      fv(4) = -1d0/2d0*alpha0/4d0/pi*qs**2*
     $     (-dlog(ms**2/mue2)+4d0+2d0*dlog(lambda**2/ms**2))
*
c add contribution from W selfenergy:
      fv(5)=alpha0/4d0/pi/2d0*(-10d0/3d0*dlog(xmw/mue2)+62d0/9d0+
     $    2d0/3d0*xmw/t-4d0*dlog(1d0-t/xmw)+2d0/3d0*(1d0-xmw/t)*lwt/t)
*
c pure photonic contribution to the formfactors:
      form = dcmplx(0d0,0d0)
      do i=1,5
         form = form+fv(i)
      end do
      return
      end
c******************************************************************
      subroutine fvboxwg_new(iqqbar,xt,xs,xqi,xqis,xqfs,xmi,xmis,
     &                                                  xmfs,form)
c     pure photonic contribution to W exchange box diagrams:
      implicit none
      include 'pwhg_wzgrad.h'
      integer iqqbar
      real*8 xt,xs,t,s,u,lwt,qqis,qqi,qqfs
      real*8 xqi,xqis,xqfs,xmi,xmis,xmfs,mi,mis,mfs
      real*8 pi,pi2
      complex*16 d0,D0_,C0_,spence,c0234,c0134,c0123,c0124
      complex*16 factor(2),form,cdlog,tc,sps,spt,spu,spu2,sps2
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)

      pi = 4d0*datan(1d0)
      pi2 = pi**2

      t=xt
      s=xs
      u=-s-t

      tc=t+ieps

      qqis=xqis
      qqi=xqi
      qqfs=xqfs

      mis=xmis
      mi=xmi
      mfs=xmfs
      
      lwt=(t-xmw)*dlog(1d0-t/xmw)
      sps=spence(1d0+xmw/s+ieps)
      spt=spence(1d0-xmw/t+ieps)
      spu=spence(1d0+xmw/u+ieps)
      spu2=spence(1d0+u/xmw+ieps)
      sps2=spence(1d0+s/xmw+ieps)
c W gamma box:
      d0=-(dlog(s**2/mis**2/mfs**2)*dlog(mw*lambda/(xmw-t))+
     $     dlog(mw/mis)**2+dlog(mw/mfs)**2+sps+
     $     pi2/3d0)/s/(t-xmw)
c      write(6,*)d0
c      d0=d0_(0d0,mfs**2,mis**2,0d0,t,s,mw,mfs,0d0,mis,0)
c      write(6,*)'d0',d0
      c0234=-1d0/t*(spt+pi2/6d0+cdlog(xmw/tc)**2/2d0+
     $     cdlog(mis**2/tc)*dlog(1d0-t/xmw))
      c0134=-1d0/t*(spt+pi2/6d0+cdlog(xmw/tc)**2/2d0+
     $     cdlog(mfs**2/tc)*dlog(1d0-t/xmw))
      c0123=-1d0/s*(dlog(s/mis/mfs)*dlog(lambda**2/s)+2d0*pi2/3d0+
     $     dlog(s/mfs**2)**2/4d0+dlog(s/mis**2)**2/4d0)
      c0124=-1d0/s*(sps2-pi2/6d0)
c      write(6,*)'c0',c0124
c      c0124=C0_(s,0d0,0d0,0d0,0d0,mw,0)
c      write(6,*)c0124
c      stop
      if(iqqbar.eq.1)
     $     factor(1) = 128d0*qqis*qqfs*s**2*(s*d0-c0234-c0134)
      if(iqqbar.eq.2)
     $     factor(1) = 64d0*qqis*qqfs*
     $     (-s*(s**2+u**2+2d0*xmw*s+xmw**2)*d0+
     $     (s**2+u**2+xmw*s+xmw*u)*(c0134+c0234)+
     $     (s*u-s**2-xmw*s)*(c0123+c0124)+2d0*u*(dlog(xmw/s)+lwt/t))
c crossed W gamma box:
      d0=-(dlog(u**2/mi**2/mfs**2)*dlog(mw*lambda/(xmw-t))+
     $     dlog(mw/mi)**2+dlog(mw/mfs)**2+spu+
     $     pi2/3d0)/u/(t-xmw)
c      d0=d0_(0d0,mfs**2,mi**2,0d0,t,u,mw,mfs,0d0,mi,0)
c      write(6,*)d0
      c0234=-1d0/t*(spt+pi2/6d0+cdlog(xmw/tc)**2/2d0+
     $     cdlog(mi**2/tc)*dlog(1d0-t/xmw))
      c0134=-1d0/t*(spt+pi2/6d0+cdlog(xmw/tc)**2/2d0+
     $     cdlog(mfs**2/tc)*dlog(1d0-t/xmw))
      c0123=-1d0/u*(dlog(-mi*mfs/u)*dlog(mi*mfs/lambda**2)+pi2/6d0-
     $     dlog(-mi*mfs/u)**2/2d0+dlog(mi/mfs)**2/2d0)
c      write(6,*)c0123
c      c0123=C0_(mi**2,u,mfs**2,0d0,mi,mfs,0)
c      write(6,*)c0123
c      stop
      c0124=-1d0/u*(spu2-pi2/6d0)
c      write(6,*)'c0',c0124
c      c0124=C0_(u,0d0,0d0,0d0,0d0,mw,0)
c      write(6,*)c0124
c      stop
      if(iqqbar.eq.1)
     $     factor(2)=64d0*qqi*qqfs*(-u*(s**2+u**2+2d0*xmw*u+xmw**2)*d0+
     $     (s**2+u**2+xmw*s+xmw*u)*(c0134+c0234)+
     $     (s*u-u**2-xmw*u)*(c0123+c0124)+2d0*s*(dlog(-xmw/u)+lwt/t))
      if(iqqbar.eq.2) factor(2)=128d0*qqi*qqfs*u**2*(u*d0-c0234-c0134)
*
      form = alpha0/4d0/pi*(factor(1)+factor(2))
      return
      end
c******************************************************************
      complex*16 function sigzz2(xs)
c     Z-self energy (has been checked by comparison with Wolf)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 azferm,azneut,azboson,cf,b0h,b1h,b20h,b0w,b1w,b20w
      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      s=xs
      sc=s+ieps
      azferm = dcmplx(0d0,0d0)
      azneut = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cf = dcmplx(f(s,m,m),g(s,m,m))
         azferm = azferm+
     $        4d0/3d0*nc*((v**2+a**2)*(-s*dlog(m2/mue2)+
     $        (2d0*m2+s)*cf-s/3d0)-a**2*6d0*m2*(-dlog(m2/mue2)+cf))
      end do
      do ff=1,3
         call prop(ff,m,q,v,a,nc)
c change on 12/14/99: -pi-> pi !
         azneut = azneut+
c     $        8d0/3d0*nc*a**2*s*dcmplx(-cdlog(sc/mue2)+5d0/3d0,pi)
     $        8d0/3d0*nc*a**2*s*dcmplx(-dlog(sc/mue2)+5d0/3d0,pi)
      end do
      call bfunc(s,mw,mw,b0w,b1w,b20w)
c      call bquer(s,mw,mw,b0w,b1w)
c      b20w = (xmw*b0w+s/2d0*b1w-7d0/24d0*s+3d0/2d0*xmw)/3d0
c      b20w = b20w-(6d0*xmw-s)/12d0*dlog(xmw/mue2)
c      b0w = b0w-dlog(xmw/mue2)
      a0w = xmw*(-dlog(xmw/mue2)+1d0)
      azboson = -cw2/sw2*(10d0*b20w+2d0*a0w+
     $     2d0*(xmw+2d0*s)*b0w-4d0*xmw+2d0/3d0*s)
      azboson = azboson+2d0*cw2/sw2*(3d0*a0w-2d0*xmw)
      azboson = azboson+xmw*sw2/cw2*2d0*b0w
      azboson = azboson+(sw2-cw2)**2/2d0/sw2/cw2*(a0w-2d0*b20w)
      azboson = azboson+cw2/sw2*b20w*2d0
      call bfunc(s,mz,mh,b0h,b1h,b20h)
c      call bquer(s,mz,mh,b0h,b1h)
c      b20h = (xmz*b0h+(s+xmz-xmh)/2d0*b1h-7d0/24d0*s+3d0/8d0*
c     $     (xmz+3d0*xmh))/3d0
c      b20h = b20h-(xmz+3d0*xmh-s)/12d0*dlog(xmh/mue2)-
c     $     xmz/6d0*dlog(xmz/mue2)
c      b0h = b0h-dlog(xmh*xmz/mue4)/2d0
      a0h = xmh*(-dlog(xmh/mue2)+1d0)
      a0z = xmz*(-dlog(xmz/mue2)+1d0)
      azboson = azboson-(b20h-xmz*b0h-a0h/4d0-
     $     a0z/4d0)/cw2/sw2
*
      sigzz2 = alpha0/4d0/pi*(azferm+azneut+azboson)
*
c      write(6,*)sigzz,alpha0/4d0/pi*azboson,
c     $     alpha0/4d0/pi*(azferm+azneut)
c      stop
      return
      end
c******************************************************************
      subroutine wselfw2(xs,sigww,dsigww)
c     W-selfenergy : (has been checked by comparison with Wolf)
c     only weak part!

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'     
      integer ff
      complex*16 b0z,b1z,b20z,db20z,b0h,b1h,b20h,db20h,b0f,b1f,
     &     b20f,db20f
      complex*16 awlep,dawlep,awquark,dawquark,awboson,dawboson,
     &     sigww,dsigww
      complex*16 sc

      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2

      s=xs
      sc=s+ieps
c     weak contribution to A_ww:    
      awlep = dcmplx(0d0,0d0)
      dawlep = dcmplx(0d0,0d0)
      do ff=4,6
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         m4 = m2**2
         if (dabs(s).eq.0d0) then
            awlep = awlep+
     $           nc*m2*dcmplx(1d0/2d0*dlog(m2/mue2)-1d0/4d0,0d0)/sw2
            dawlep = 0d0
         else
            awlep = awlep+
     $           nc*dcmplx(-(s-3d0/2d0*m2-m4/s)*dlog(m2/mue2)+
     $           (s-m2/2d0-m4/s/2d0)*(m2/s-1d0)*dlog(dabs(1d0-s/m2))+
     $           5d0/3d0*s-m2-3d0/2d0*m4/s,-pi*theta(s-m2)*(s-m2/2d0-
     $           m4/s/2d0)*(m2/s-1d0))/3d0/sw2
            dawlep = dawlep+
     $           nc*dcmplx(-(1d0+m4/s**2)*dlog(m2/mue2)+
     $           ((m2/s)**2/2d0+(m2/s)**3-1d0)*dlog(dabs(1d0-s/m2))+
     $           2d0/3d0+2d0*m4/s**2+m2/2d0/s,-pi*theta(s-m2)*(1d0-
     $           (m2/s)**3))/3d0/sw2
         endif
      end do
      awquark = dcmplx(0d0,0d0)
      dawquark = dcmplx(0d0,0d0)
      do  ff=7,9
         call prop(ff,m,q,v,a,nc)
         call prop(ff+3,ms,qs,vs,as,ncs)
         m2 = m**2
         ms2 = ms**2
         call bquer(s,m,ms,b0f,b1f)
         b20f = (m2*b0f+(s+m2-ms2)/2d0*b1f-7d0/24d0*s+3d0/8d0*
     $        (m2+3d0*ms2))/3d0
         b20f = b20f-(m2+3d0*ms2-s)/12d0*dlog(ms2/mue2)-
     $        m2/6d0*dlog(m2/mue2)
         b0f = b0f-dlog(m2*ms2/mue4)/2d0
         a0f = m2*(-dlog(m2/mue2)+1d0)
         a0fs = ms2*(-dlog(ms2/mue2)+1d0)
         awquark = awquark+nc*(2d0*b20f-a0f/2d0
     $        -a0fs/2d0-(m2+ms2-s)*b0f/2d0)/sw2
*
         if (s.eq.0d0) then
            dawquark = 0d0
            goto 100
         else
            call bquer1(s,m,ms,p0,p1,db0f,db1f)
            call bquer(s,m,ms,b0f,b1f)
            db20f = (m2*db0f+(s+m2-ms2)/2d0*db1f+b1f/2d0-
     $           7d0/24d0)/3d0
            db20f = db20f+dlog(ms2/mue2)/12d0
            b0f = b0f-dlog(m2*ms2/mue4)/2d0
            dawquark = dawquark+nc*(2d0*db20f-(m2+ms2-s)*db0f/2d0+
     $           b0f/2d0)/sw2
         end if
 100     continue
      end do
      call bquer(s,mz,mw,b0z,b1z)
      call bquer(s,mw,mh,b0h,b1h)
      b20z = (xmz*b0z+(s+xmz-xmw)/2d0*b1z-7d0/24d0*s+3d0/8d0*
     $     (xmz+3d0*xmw))/3d0
      b20h = (xmw*b0h+(s+xmw-xmh)/2d0*b1h-7d0/24d0*s+3d0/8d0*
     $     (xmw+3d0*xmh))/3d0
      b20z = b20z-(xmz+3d0*xmw-s)/12d0*dlog(xmw/mue2)-
     $     xmz/6d0*dlog(xmz/mue2)
      b20h = b20h-(xmw+3d0*xmh-s)/12d0*dlog(xmh/mue2)-
     $     xmw/6d0*dlog(xmw/mue2)
      a0z = xmz*(-dlog(xmz/mue2)+1d0)
      a0w = xmw*(-dlog(xmw/mue2)+1d0)
      a0h = xmh*(-dlog(xmh/mue2)+1d0)
      b0z = b0z-dlog(xmw*xmz/mue4)/2d0
      b0h = b0h-dlog(xmw*xmh/mue4)/2d0
      awboson = -cw2/sw2*(10d0*b20z-2d0*xmz-2d0*xmw+2d0/3d0*s+
     $     a0z+a0w+(xmw+xmz+4d0*s)*b0z)
      awboson = awboson+xmw*sw2/cw2*b0z
      awboson = awboson-b20z/sw2
      awboson = awboson+cw2/sw2*b20z*2d0
      awboson = awboson+(xmw*b0h-b20h)/sw2
      awboson = awboson+cw2/sw2*(3d0*a0z-2d0*xmz)
      awboson = awboson+(3d0*a0w-2d0*xmw)/sw2
      awboson = awboson+(a0h+a0z+2d0*a0w)/4d0/sw2
*
      if (dabs(s).eq.0d0) then
         dawboson = 0d0
         goto 200
      else
         call bquer1(s,mz,mw,p0,p1,db0z,db1z)
         call bquer1(s,mw,mh,p0,p1,db0h,db1h)
         call bquer(s,mz,mw,b0z,b1z)
         call bquer(s,mw,mh,b0h,b1h)
         db20z = (xmz*db0z+(s+xmz-xmw)/2d0*db1z+b1z/2d0-
     $        7d0/24d0)/3d0
         db20h = (xmw*db0h+(s+xmw-xmh)/2d0*db1h+b1h/2d0-
     $        7d0/24d0)/3d0
         db20z = db20z+dlog(xmw/mue2)/12d0
         db20h = db20h+dlog(xmh/mue2)/12d0
         b0z = b0z-dlog(xmw*xmz/mue4)/2d0
         dawboson = -cw2/sw2*(10d0*db20z+2d0/3d0+
     $        (4d0*s+xmw+xmz)*db0z+4d0*b0z)
         dawboson = dawboson+xmw*sw2/cw2*db0z
         dawboson = dawboson-db20z/sw2
         dawboson = dawboson+cw2/sw2*db20z*2d0
         dawboson = dawboson+(xmw*db0h-db20h)/sw2
      endif
 200  continue
*
      sigww = alpha0/4d0/pi*(awlep+awquark+awboson)
      dsigww = alpha0/4d0/pi*(dawlep+dawquark+dawboson)
*
c      write(6,*)sigww,alpha0/4d0/pi*(awlep+awquark),
c     $     alpha0/4d0/pi*awboson
c      stop
*
      end
c******************************************************************
      complex*16 function siggg2(xs)
c photon-self energy (has been checked by comparison with Wolf)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff
      complex*16 agferm,cff,cfw,agboson

      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      s=xs
      agferm = dcmplx(0d0,0d0)
      do ff=4,12
         call prop(ff,m,q,v,a,nc)
         m2 = m**2
         cff = dcmplx(f(s,m,m),g(s,m,m))
         agferm = agferm+
     $        4d0/3d0*nc*q**2*(-s*dlog(m2/mue2)-s/3d0+
     $        (2d0*m2+s)*cff)
      end do
      cfw = dcmplx(f(s,mw,mw),g(s,mw,mw))
      agboson = 3d0*s*dlog(xmw/mue2)-(4d0*xmw+3d0*s)*cfw
*
      siggg2 = alpha0/4d0/pi*(agferm+agboson)
*
      end
c******************************************************************
      subroutine fdeltar(dr)
c     determination of delta_r:

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      complex*16 sww0,dsww0,swwos,dswwos,swwp0,swwpos,sigwwp

      pi = 4d0*datan(1d0)
      pi2 = pi**2
c     calculation of the renormalisation constant dz2w:
      call wrenorm(dz2w)

      call wselfw(xmw,swwos,dswwos)
      call wselfw(0d0,sww0,dsww0)
      swwp0 = sigwwp(0d0)
      swwpos = sigwwp(xmw)
      rsigw0 = dreal(sww0+swwp0)-dreal(swwos+swwpos)-xmw*dz2w
*
      dr = rsigw0/xmw+alpha0/4d0/pi/sw2*(6d0+
     $     (7d0-4d0*sw2)/2d0/sw2*dlog(cw2))
*
      end
c******************************************************************
      double precision function fermvac(xs)
c     calculation of the vacuum polarisation:

      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
      integer ff

      common/renorm/mue,mue2,mue4
      pi = 4d0*datan(1d0)
      pi2 = pi**2
      s=xs
      pig0f = 0d0
      agferm = 0d0
      do ff=7,12
          call prop(ff,m,q,v,a,nc)
          pig0f = pig0f-alpha0/4d0/pi*
     $        4d0/3d0*nc*q**2*dlog(m**2/mue2)
          agferm = agferm+
     $         4d0/3d0*nc*q**2*(-s*dlog(m**2/mue2)-s/3d0+
     $         (2d0*m**2+s)*f(s,m,m))
      end do
c add leptons:
      do ff=4,6
          call prop(ff,m,q,v,a,nc)
          pig0f = pig0f-alpha0/4d0/pi*
     $        4d0/3d0*nc*q**2*dlog(m**2/mue2)
          agferm = agferm+
     $         4d0/3d0*nc*q**2*(-s*dlog(m**2/mue2)-s/3d0+
     $         (2d0*m**2+s)*f(s,m,m))
      end do
      siggg = alpha0/4d0/pi*agferm
*
      fermvac = (siggg-s*pig0f)/s
*
      return
      end
c******************************************************************
      subroutine higgsdr(mwerror)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
c     calculation of M_higgs from delta_r with M_w as input
c     (using regula falsi)
      integer mwerror
      parameter (wtol = 1d-5)

      common/cdfvalue/mwexp
      common/drit/deltar 
      common/hbounds/hlow,hhigh

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      a = hlow
      b = hhigh
      mh = a
      xmh = mh**2
      call mwiter(mwerror)
      fa = mw-mwexp
      mh = b
      xmh = mh**2
      call mwiter(mwerror)
      fb = mw-mwexp
      if ((fa*fb).gt.0d0) then
         if (fa.gt.0d0) then
            if (fa.lt.fb) then
               mh = a
               xmh = mh**2
               call mwiter(mwerror)
               goto 21
            else
               mh = b
               xmh = mh**2
               call mwiter(mwerror)
               goto 21
            end if
         else
            if (fa.lt.fb) then
               mh = b
               xmh = mh**2
               call mwiter(mwerror)
               goto 21
            else
               mh = a
               xmh = mh**2
               call mwiter(mwerror)
               goto 21
            end if
         end if
      end if
 20   continue
      x = (a*fb-b*fa)/(fb-fa)
      mh = x
      xmh = mh**2
      call mwiter(mwerror)
      fx = mw-mwexp
      if (dabs(fx).lt.wtol) goto 21
      if (fa*fx.lt.0d0) then
         b = x
         fb = fx
         fa = fa/2d0
         goto 20
      else 
         a = x
         fa = fx
         fb = fb/2d0
         goto 20
      end if
 21   continue
      end
c******************************************************************
      subroutine higgs_sirlin(mwerror)
      implicit real*8(a-z)
      include 'pwhg_wzgrad.h'
c     calculation of M_higgs from aprox. delta r with M_w as input
c     (using regula falsi)
      integer mwerror
      parameter (wtol = 1d-5)

      common/cdfvalue/mwexp
      common/drit/deltar 
      common/hbounds/hlow,hhigh

      pi = 4d0*datan(1d0)
      pi2 = pi**2
      a = hlow
      b = hhigh
      mh = a
      xmh = mh**2
      mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $     0.0078d0*dlog(mh/100d0)**2
      fa = mw-mwexp
      mh = b
      xmh = mh**2
      mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $     0.0078d0*dlog(mh/100d0)**2
      fb = mw-mwexp
      if ((fa*fb).gt.0d0) then
         if (fa.gt.0d0) then
            if (fa.lt.fb) then
               mh = a
               xmh = mh**2
               mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $              0.0078d0*dlog(mh/100d0)**2
               goto 21
            else
               mh = b
               xmh = mh**2
               mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $              0.0078d0*dlog(mh/100d0)**2
               goto 21
            end if
         else
            if (fa.lt.fb) then
               mh = b
               xmh = mh**2
               mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $              0.0078d0*dlog(mh/100d0)**2
               goto 21
            else
               mh = a
               xmh = mh**2
               mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $              0.0078d0*dlog(mh/100d0)**2
               goto 21
            end if
         end if
      end if
 20   continue
      x = (a*fb-b*fa)/(fb-fa)
      mh = x
      xmh = mh**2
      mw=80.3805d0-0.0581d0*dlog(mh/100d0)-
     $     0.0078d0*dlog(mh/100d0)**2
      fx = mw-mwexp
      if (dabs(fx).lt.wtol) goto 21
      if (fa*fx.lt.0d0) then
         b = x
         fb = fx
         fa = fa/2d0
         goto 20
      else 
         a = x
         fa = fx
         fb = fb/2d0
         goto 20
      end if
 21   continue
      end
