c******************************************************************        
        subroutine prop(f,mf,qqf,vf,af,ncf)

        implicit real*8(a-z)
        include 'pwhg_wzgrad.h'
	real*8 m(12),q(12),t3(12),v(12),a(12),nc(12)
        integer i,f

        do i=1,3
           q(i) = 0d0
           t3(i) = 1d0/2d0
           m(i) = 1d-5
           nc(i) = 1d0
        end do
        do i=4,6
           q(i) = -1d0
           t3(i) = -1d0/2d0
           nc(i) = 1d0
        end do
        do i=7,9
           q(i) = 2d0/3d0
           t3(i) = 1d0/2d0
           nc(i) = 3d0
        end do
        do i=10,12
           q(i) = -1d0/3d0
           t3(i) = -1d0/2d0
           nc(i) = 3d0
        end do
        qqf = q(f)
        ncf = nc(f)

        do i=1,12
           v(i) = (t3(i)-2d0*q(i)*sw2)/2d0/sw/cw
           a(i) = t3(i)/2d0/sw/cw
        end do
        vf = v(f)
        af = a(f)

        m(4) = me
        m(5) = mmu
        m(6) = 1.77684d0
        m(7) = meff
        m(8) = 1.2d0
        m(9) = mtop
        m(10) = meff+1d-5
        m(11) = 150d-3 
        m(12) = 4.6d0

        mf = m(f)

        return
        end
c***********************************************************************
        subroutine winput
c       weak input:
        implicit none
        include 'pwhg_wzgrad.h'
        integer i,ff,mwerror,mix,hqcd,twol,had
	real*8 mue,mue2,mue4,deltar,mw0,s2eff0,s2eff
	real*8 mwexp,pi
	real*8 f,hadvac,mf,qqf,vf,af,ncf
	real*8 pig0,sgz0,dz2z,dz2g,dz2gz
	real*8 dhad,dferm,dtop,dlep,disint,alpzinv,alpmz,hlow,hhigh
	real*8 ckm2(7:8,10:12)
        common/renorm/mue,mue2,mue4
        common/kobaya/ckm2
        common/drit/deltar 
        common/gzmix/mix
        common/higher/hqcd,twol
        common/dalpha/dhad,dferm,disint
        common/jeger/had
        common/inv/alpzinv
        common/cdfvalue/mwexp
        common/hbounds/hlow,hhigh
        common/sigzero/pig0,sgz0
        common/zgrenorm/dz2z,dz2g,dz2gz

        mue=10d0
        mue2=mue**2
        mue4=mue2**2
        pi=4d0*datan(1d0)
        mwexp=mw

c       mw=mwexp
c       Jegerlehners dispersion integral and the resulting effective 
c       light quark mass:
        disint=0.027572d0
c       meff=0.0465d0

c use value for Milano W mass workshop studies:
        meff=69.83d-3
        dhad=disint
        alpzinv=128.93d0

c approximation of rad. corr. to MW by Sirlin et al.
c       mw0=80.3805d0
c       mw=mw0-0.0581d0*dlog(mh/100d0)-0.0078d0*dlog(mh/100d0)**2-
c       1    0.518d0*(dhad/0.028d0-1d0)+
c       2    0.537d0*((mtop/175d0)**2-1d0)-
c       2    0.085d0*(alphas/0.118d0-1d0)

c approximation of rad. corr. to sl2eff by Sirlin et al.
        s2eff0=0.23154d0
        s2eff=s2eff0+0.000526d0*dlog(mh/100d0)+
     1    0.00986d0*(dhad/0.028d0-1d0)-
     2    0.00268d0*((mtop/175d0)**2-1d0)+
     2    0.00044d0*(alphas/0.118d0-1d0)
c       write(6,*)mh,mw,s2eff

c Kobayashi-Maskawa matrix elements squared:
c only used in the calculation of the W width
        ckm2(7,10) = 0.975d0**2
        ckm2(7,11) = 0.222d0**2
        ckm2(7,12) = 0d0
        ckm2(8,10) = 0.222d0**2
        ckm2(8,11) = 0.975d0**2
        ckm2(8,12) = 0d0

c calculation of pigamma(0) using dispersion integral (had = 1):
        had=0

c calculation of the renormalisation constant dz2w
c (with inclusion of higher order terms: 
c corr. to rho param.:QCD and two-loop)
c and gZ-mixing in Z renormalis. constant; (0: off, 1: on):
        hqcd=0
        twol=0
        mix=1

c calculation of meff from delta alpha_hadronic:
        dhad=hadvac(xmz)

c calculation of alpha(xmz):
        dlep=0d0
        do ff=4,6
           call prop(ff,mf,qqf,vf,af,ncf)
           dlep=dlep+alpha/pi/3d0*((2d0*mf**2/xmz+1d0)*
     1       f(xmz,mf,mf)-1d0/3d0)
        end do
        dtop=alpha/pi*4d0/9d0*((2d0*mtop**2/xmz+1d0)*
     1    f(xmz,mtop,mtop)-1d0/3d0)
        dferm=-dlep-dhad-dtop
        alpmz=alpha/(1d0-dferm)
c       write(6,*)'hadvac',dhad,dferm,1d0/alpmz

c determination of MH with MH as input:
        mwerror=0

c       call mwiter(mwerror)
c       if (mwerror.ne.0) then
c          write(6,*)'mw not found for mhiggs=',
c       1       mh,' and mtop=',mtop 
c       endif
c        xmw = mw**2
c        cw2 = xmw/xmz
c        cw = dsqrt(cw2)
c        sw2 = 1d0-cw2
c        sw = dsqrt(sw2)
c        write(6,*)deltar,mtop,mw,mh
c calculation of delta r:

        call fdeltar(deltar)
c       write(6,*)'delta_r=',deltar

c calculation of the Z renormalisation constants:
        call zrenorm(dz2z)
        dz2g=-pig0
        dz2gz=cw*sw/(cw2-sw2)*(dz2z-dz2g)
        return
        end
c**********************************************************************
