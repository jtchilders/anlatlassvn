      subroutine setvirtual_EW(vflav,born,virt_soft)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'pwhg_wzgrad.h'
      integer vflav(nlegborn),j,i
      real*8 virt_soft,born,CKM,dotp,shat,o_fac,sig012(2),a2qqw
      real*8 t1,t2,sig_virt12,sig_virt21,sig_virt,k_fac(2),k_faci(2)
      real*8 sig_soft12,sig_soft21,sig_soft,soft,virt
      complex*16 fvw,fvwp(2),prop,prop2
      common/fvwpfactors/fvwp
      logical fvwpini
      common/flg_fvwp/fvwpini
      real*8 powheginput
      external powheginput,dotp
      logical flg_inbtilde,flg_inequiv
      common/pwhg_flg_EW/flg_inbtilde,flg_inequiv
      if(flg_inbtilde)then
      call init_phys_EW
         if(flg_wgrad2.eq.1)then
!         call CKMs(vflav,CKM)
         shat=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2))
!        choice of constant of s-dependent width
         if(wopt.eq.1)then
         prop=1d0/dcmplx(shat-xmw,ph_WmWw)
         elseif(wopt.eq.2)then
         prop=1d0/dcmplx(shat-xmw,ph_Wwidth*shat/mw)
         endif
         prop2=prop*dconjg(prop)
!	 choice of representation
         if(rep.eq.1)then
         o_fac=(pi*ph_alphaem/ph_sthw2)**2/3d0
         elseif(rep.eq.2)then
         o_fac=(w2*xmw*gfermi)**2/3d0
         endif
         if(powheginput('idvecbos').eq.24)then
         sig012(1)=4d0*(2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,3)))**2 
         sig012(2)=4d0*(2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,3)))**2 
         elseif(powheginput('idvecbos').eq.-24)then
         sig012(1)=4d0*(2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,4)))**2 
         sig012(2)=4d0*(2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,4)))**2 
         endif
!        ----------------------virtual finite piece--------------------
c weak fromfactors only need to be calculated once
         if(fvwpini)then
            do j=1,2
               call formweak(j,fvw)
               fvwp(j)=fvw
            enddo
            fvwpini=.false.
         endif
         if(qnonr.eq.0)then     !resonant weak virtual contrib
            virt=born*dreal(fvwp(1)+fvwp(2))*2d0*pi/st_alpha
         elseif(qnonr.eq.1)then !nonresonant weak virtual contrib
               if(dabs(shat-ph_Wmass2).le.1d0)then
               virt=born*dreal(fvwp(1)+fvwp(2))*2d0*pi/st_alpha 
               else
                  if(powheginput('idvecbos').eq.24)then
                  t1=-2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,4))
                  t2=-2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,4))
                  CKM=ph_CKM(1,3)
                  elseif(powheginput('idvecbos').eq.-24)then
                  t1=-2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,3))
                  t2=-2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,3))
                  CKM=ph_CKM(2,1)
                  endif
               sig_virt12=o_fac*prop2*a2qqw(shat,t1,t2,sig012(1))
               sig_virt21=o_fac*prop2*a2qqw(shat,t2,t1,sig012(2))
                  if(vflav(1).lt.0)then
                  sig_virt=sig_virt21
                  else
                  sig_virt=sig_virt12
                  endif
               virt=CKM**2*sig_virt*2d0*pi/st_alpha
               endif
            endif
!        ----------------------soft finite piece--------------------
            if(test(3).le.3)then
            call softcollqed_cc(test(3),k_fac)
            elseif(test(3).eq.4)then
            k_fac(1)=0d0
            k_fac(2)=0d0
               do i=1,3
               call softcollqed_cc(i,k_faci)
               k_fac(1)=k_fac(1)+k_faci(1)
               k_fac(2)=k_fac(2)+k_faci(2)
               enddo
            endif

            if(powheginput('idvecbos').eq.24)then
            sig_soft12=o_fac*k_fac(1)*sig012(1)*prop2
            sig_soft21=o_fac*k_fac(2)*sig012(2)*prop2
            CKM=ph_CKM(1,3)
            elseif(powheginput('idvecbos').eq.-24)then
            sig_soft12=o_fac*k_fac(2)*sig012(1)*prop2
            sig_soft21=o_fac*k_fac(1)*sig012(2)*prop2
            CKM=ph_CKM(2,1)
            endif

            if(vflav(1).lt.0)then
            sig_soft=sig_soft21
            else
            sig_soft=sig_soft12
            endif
         soft=CKM**2*sig_soft*2d0*pi/st_alpha
         else
         virt_soft=0d0
         endif
!        ----------------------final result--------------------
      virt_soft = virt + soft 
c      virt_soft = soft !WZGRAD TEST EDIT
      else 
      virt_soft=0d0
      endif

      end
c***********************************************************************
      subroutine CKMs(vflav,CKM)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      integer i,j,vflav(nlegborn)
      real*8 CKM
         do i=1,2
            do j=2,1,-1
            if(vflav(i).eq.-5)then

               if(vflav(j).eq.2)then
               CKM=ph_CKM(1,3)
               elseif(vflav(j).eq.4)then
               CKM=ph_CKM(2,3)
               endif

            elseif(vflav(i).eq.-3)then

               if(vflav(j).eq.2)then
               CKM=ph_CKM(1,2)
               elseif(vflav(j).eq.4)then
               CKM=ph_CKM(2,2)
               endif

            elseif(vflav(i).eq.-1)then

               if(vflav(j).eq.2)then
               CKM=ph_CKM(1,1)
               elseif(vflav(j).eq.4)then
               CKM=ph_CKM(2,1)
               endif

            endif
            enddo
         enddo
      end
c***********************************************************************
      subroutine init_phys_EW
      implicit none

      include 'nlegborn.h'
      include 'pwhg_wzgrad.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'

      real * 8 powheginput,dotp
      external powheginput,dotp
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass

      !choice to include EW corrections
      flg_wgrad2 = powheginput('wgrad2')
      !coupling constants
      alpha0=ph_alphaem !used in libqed.f, libweak.f
      alphas= 0.12!used in libqed.f, libweak.f
      !masses
      me=0.51099891d-3
      mmu=0.1056583668d0

      mff=1d-5 !neutrino mass used only for QED part
      if(powheginput('vdecaymode').eq.1)then
         mffs=lepmass(1) 
      elseif(powheginput('vdecaymode').eq.2)then
         mffs=lepmass(2)
      else
         write(6,*)'only vdecaymode=1,2 is allowed'
         stop
      endif

      mii=1d-5 !common mass for IS q's-used to regulate QED coll sing
      mmi = mii
      mmis = mmi
      mmf = 1d-5
      mmfs = mffs
      qi = 2d0/3d0
      qis = -1d0/3d0
      qf = 0d0
      qfs = -1d0

      mtop=171.2d0
      mh=115d0

      mw=ph_Wmass
      mz=ph_Zmass
      mv=mw !for now only W production
      xmh=mh**2
      xmw=mw**2
      xmz=mz**2

      !widths
      wopt=1 !switch for width (0: s depdt, 1: constant width) 
      gamw=ph_Wwidth
      gamz=ph_Zwidth
      gamv=gamw

      !mixing angles
      cw=ph_cthw
      sw=ph_sthw
      sw2=sw**2
      cw2=cw**2

      !factorization and renormalization scales set in 
      !subroutine set_fac_ren_scales in W/Born_phsp.f
      !because momenta are needed to define it ito 
      !M_invmass

      !PSS parameters
      deltas=0.01d0 
      deltac=0.005d0 

      !miscellaneous
      gfermi=1.16637d-5
      w2=dsqrt(2d0)

      !options
      rep=1
         if(rep.eq.2)then
         alpha=xmw*sw2/pi*w2*gfermi
         elseif(rep.eq.1)then
         alpha=alpha0
         endif
      qcd=0 !inclusion of QCD in W width(0:n, 1:y)
      mzero=1 !all fermions to be massless(0:n, 1:y)
      !choice of factorizations scheme - 0: MSbar    1:DIS
      lfc=powheginput('lfc') 
      !factorization and renormalization scales
      call set_fac_ren_scales(mu_f,mu_r)
      !collinear FSR cut - 0:no cut   1:with cut
      !if massless fermions, collcut must be 1 - masses 
      !regulate collinear singularities
      collcut=0 
      !fictitious photon mass
      lambda=1d0
      !small imaginary part 
      ieps=dcmplx(0d0,1d-18) 
      !QED corrxns - 1: ISR   2:FSR   3:interference   4:all
      test(3)=powheginput('QED') 
      !W-resonance contributions - 0:res only   1:nonres contr
      qnonr=powheginput('qnonr') 
      if(qnonr.eq.1)test(3)=4
      !FS radiation - 1:full FSR   2:FSR a la Berends etal.
      test(4)=1
      call winput !further prep for input.f, libweak.f and libqed.f 

      !calculate invariants for use in libqed.f
      br_sinv(1,2)=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,2))
      br_sinv(1,3)=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,3))
      br_sinv(2,3)=2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,3))
      br_sinv(1,4)=2d0*dotp(kn_cmpborn(0,1),kn_cmpborn(0,4))
      br_sinv(2,4)=2d0*dotp(kn_cmpborn(0,2),kn_cmpborn(0,4))
      br_sinv(3,3)=2d0*dotp(kn_cmpborn(0,3),kn_cmpborn(0,3))
      br_sinv(3,4)=2d0*dotp(kn_cmpborn(0,3),kn_cmpborn(0,4))
      br_sinv(4,4)=2d0*dotp(kn_cmpborn(0,4),kn_cmpborn(0,4))
      br_sinv(2,1)=br_sinv(1,2)
      br_sinv(3,1)=br_sinv(1,3)
      br_sinv(4,1)=br_sinv(1,4)
      br_sinv(3,2)=br_sinv(2,3)
      br_sinv(4,2)=br_sinv(2,4)
      br_sinv(4,3)=br_sinv(3,4)
      end
!***********************************************************************
      subroutine CKMconvert(CKM) !WZGRAD EDIT
      implicit none
      include 'PhysPars.h'
      real*8 CKM(12)
      character*20 pwgprefix
      real*8 powheginput
      external powheginput

      if(powheginput('idvecbos').eq.24)then
      CKM(1) = ph_CKM(1,3)
      CKM(2) = ph_CKM(2,3)
      CKM(3) = ph_CKM(1,2)
      CKM(4) = ph_CKM(2,2)
      CKM(5) = ph_CKM(1,1)
      CKM(6) = ph_CKM(2,1)
      CKM(7) = ph_CKM(1,3)
      CKM(8) = ph_CKM(1,2)
      CKM(9) = ph_CKM(1,1)
      CKM(10) = ph_CKM(2,3)
      CKM(11) = ph_CKM(2,2)
      CKM(12) = ph_CKM(2,1)
      elseif(powheginput('idvecbos').eq.-24)then
      CKM(1) = ph_CKM(2,1)
      CKM(2) = ph_CKM(2,2)
      CKM(3) = ph_CKM(2,3)
      CKM(4) = ph_CKM(1,1)
      CKM(5) = ph_CKM(1,2)
      CKM(6) = ph_CKM(1,3)
      CKM(7) = ph_CKM(2,1)
      CKM(8) = ph_CKM(1,1)
      CKM(9) = ph_CKM(2,2)
      CKM(10) = ph_CKM(1,2)
      CKM(11) = ph_CKM(2,3)
      CKM(12) = ph_CKM(1,3)
      endif

      end
!******************************************************************






















