      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'mathx.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_physpar.h'
      include 'LesHouches.h'
      real * 8 masswindow_low,masswindow_high
      real * 8 mass_low,mass_high

      real * 8 powheginput
      external powheginput
      logical verbose
      parameter(verbose=.true.)
      real *8 decmass
      common/clepmass/decmass
      real*8 mlep2
      common/leptmass/mlep2
c     renormalization scheme
c     0 -> alpha(0)
c     1 -> alpha(mz)
c     2 -> gmu
      integer cmpmass
      integer iftoptmp
      real * 8 cmass, bmass
      integer j
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/137.03599911d0

      physpar_aem = ph_alphaem

c Boson masses and witdhs
      ph_Zmass = powheginput("#Zmass")
      if (ph_Zmass.le.0d0) ph_Zmass  = 91.1876d0     
      ph_Zwidth = powheginput("#Zwidth")
      if (ph_Zwidth.le.0d0) ph_Zwidth =  2.4924d0

      ph_Wmass = powheginput("#Wmass")
      if (ph_Wmass.le.0d0) ph_Wmass  = 80.37399d0     
      ph_Wwidth = powheginput("#Wwidth")
      if (ph_Wwidth.le.0d0) ph_Wwidth =  2.0836d0
      ph_Hmass = powheginput("#Hmass")

      if (ph_Hmass.le.0d0) ph_Hmass  = 125.d0
      ph_Hwidth = powheginput("#Hwidth")
      if (ph_Hwidth.le.0d0) ph_Hwidth =  0.01d0

c Quark masses
      ph_mTop = powheginput("#Tmass")
      if (ph_mTop.le.0d0) ph_mTop = 174d0

      ph_mBot = powheginput("#Bmass")
      if (ph_mBot.le.0d0) ph_mBot = 4.6d0

      ph_mCha = powheginput("#Cmass")
      if (ph_mCha.le.0d0) ph_mCha = 1.2d0

      ph_mStr = powheginput("#Smass")
      if (ph_mStr.le.0d0) ph_mStr = 0.15d0

      ph_mUp = powheginput("#Umass")
      if (ph_mUp.le.0d0) ph_mUp = 0.06983d0

      ph_mDown = powheginput("#Dmass")
      if (ph_mDown.le.0d0) ph_mDown = 0.06984d0


c Lepton masses
      ph_mEl = powheginput("#Elmass")
      if (ph_mEl.le.0d0) ph_mEl = 0.51099892d-3

      ph_mMuon = powheginput("#Mumass")
      if (ph_mMuon.le.0d0) ph_mMuon = 0.105658369d0

      ph_mTau = powheginput("#Taumass")
      if (ph_mTau.le.0d0) ph_mTau = 1.77699d0


      physpar_ml(1)=ph_mEl
      physpar_ml(2)=ph_mMuon
      physpar_ml(3)=ph_mTau

c     Set here lepton and quark masses for momentum reshuffle in the LHE event file

c     number of light flavors
      st_nlight = 5

      do j=1,st_nlight         
         physpar_mq(j)=0d0
      enddo
c     read eventual c and b masses from the input file
      cmass=powheginput("#cmass_lhe")
      if (cmass.gt.0d0) physpar_mq(4)=cmass
      bmass=powheginput("#bmass_lhe")
      if (bmass.gt.0d0) physpar_mq(5)=bmass


      scheme = powheginput("#scheme")
      if (scheme.lt.0d0) scheme = 0

      iftoptmp = powheginput("#iftopinloop")
      if (iftoptmp.le.0d0) iftoptmp = 0
      if (iftoptmp.eq.0) then
          iftopinloop = .false.
      else
          iftopinloop = .true.
      endif

      ph_gmu = powheginput("#gmu")
      if (ph_gmu.le.0) ph_gmu = 1.16637d-5

      ph_alpha_mz = powheginput("#alphaem_z")
      if (ph_alpha_mz.le.0d0) ph_alpha_mz = 1d0/128.93d0

c     number of light flavors
      st_nlight = 5

c     mass window
      mass_low = powheginput("#mass_low")
      if (mass_low.lt.0d0) mass_low=-1d0
      mass_high = powheginput("#mass_high")
      if (mass_high.lt.0d0) mass_high=-1d0    

      masswindow_low = powheginput("#masswindow_low")
      if (masswindow_low.le.0d0) masswindow_low=30d0
      masswindow_high = powheginput("#masswindow_high")
      if (masswindow_high.le.0d0) masswindow_high=30d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


*
      epsilon = 1d-40

*
** fermion masses
*
      me  = ph_mEl
      mm  = ph_mMuon
      mtl = ph_mTau
      me2 = me*me
      mm2 = mm*mm
      mtl2= mtl*mtl

      mu= ph_mUp
      md= ph_mDown
      mu2= mu*mu
      md2= md*md

      mc= ph_mCha
      ms= ph_mStr
      mc2= mc*mc
      ms2= ms*ms

      mt= ph_mTop
      mb= ph_mBot
      mt2= mt*mt
      mb2= mb*mb

*  mass of final state lepton
      if(lprup(1).eq.10011) then
         decmass=ph_mEl
      elseif(lprup(1).eq.10013) then
         decmass=ph_mMuon
      elseif(lprup(1).eq.10015) then
         decmass=ph_mTau
      endif   

      mlep2=decmass*decmass

*
** couplings / boson masses (cw = mw/mz)
*
      cmpmass = powheginput("#complexmasses")
      if (cmpmass.lt.0) cmpmass = 0

      if (cmpmass.eq.0) then
          complexmasses = .false.
          mw2 = ph_Wmass**2
          mz2 = ph_Zmass**2
      else
          complexmasses = .true.
          mw2 = ph_Wmass**2 - ii*ph_Wmass*ph_Wwidth
          mz2 = ph_Zmass**2 - ii*ph_Zmass*ph_Zwidth
      endif

      mh2= ph_Hmass**2

      mw = sqrt(mw2)
      mz = sqrt(mz2)
      mh = sqrt(mh2)

      cw = mw/mz
      cw2= cw*cw
      sw2= 1.d0-cw2
      sw = sqrt(sw2)

      sw4= sw2*sw2
      cw4= cw2*cw2

      ph_sthw2= 1d0 - mw2/mz2
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass*ph_Zmass
      ph_Wmass2 = ph_Wmass*ph_Wmass
      ph_WmWw = ph_Wmass * ph_Wwidth
      ph_ZmZw = ph_Zmass * ph_Zwidth


c     set mass windows around Z-mass peak in unit of ph_Zwidth
c     It is used in the generation of the Born phase space

      if(mass_low.ge.0d0) then
         ph_Zmass2low=mass_low
      else
         ph_Zmass2low=max(0d0,ph_Zmass-masswindow_low*ph_Zwidth)
      endif
      if(mass_high.gt.0d0) then
         ph_Zmass2high=mass_high
      else
         ph_Zmass2high=ph_Zmass+masswindow_high*ph_Zwidth
      endif

      if (ph_Zmass2low.lt.1d0) then
         write(*,*) '*************************************'
         write(*,*) 'WARNING: Z virtuality cutoff at 1 GeV'
         write(*,*) '         to avoid the photon pole    '
         write(*,*) '*************************************'
         ph_Zmass2low=1d0
      endif
      ph_Zmass2low=ph_Zmass2low**2
      ph_Zmass2high=min(kn_sbeams,ph_Zmass2high**2)

      if( ph_Zmass2low.ge.ph_Zmass2high ) then
         write(*,*) "Error in init_couplings: mass_low >= mass_high"
         call exit(1)
      endif

*
      el2 = ph_alphaem*4.d0*pi
      alsu4pi = ph_alphaem/4d0/pi
      ph_unit_e = sqrt(4d0*pi*ph_alphaem)

      if (scheme.eq.0) then
          alpha   = dcmplx(ph_alphaem)
          el2_scheme = alpha*4d0*pi
      elseif (scheme.eq.1) then
          alpha   = dcmplx(ph_alpha_mz)
          el2_scheme = alpha*4d0*pi
      else
          el2_scheme = ph_gmu * 8d0/sqrt(2d0) * mw2 * sw2
          alpha   = el2_scheme/4d0/pi
      endif
*
* couplings of vectors to fermions
*
* given the expression (gv - ga gam_5) it is decomposed as 
*
* (gv - ga gam_5)= gm omega- + gp omega+ (see eqs. (a.14) and (a.15) 
* and comment after eq. (6.13) of arxiv:0709.1075 (denner fortschritte), 
* where omega- = (1 - gam_5)/2 and omega+ = (1 + gam_5)/2 (see eq. (2.9) 
* of arxiv:0709.1075). as a consequence
*
* gv= (gm+gp)/2
* ga= (gm-gp)/2
*
* gm= gv + ga
* gp= gv - ga
*
* using eq. (a.14) of arxiv:0709.1075 (Denner Fortschritte) 
* (in agreement (for the gfm couplings) with Dittmaier-Kramer prd65 073007) 
*
      qu =  2d0/3d0
      qd = -1d0/3d0
      ql = -1d0

      gl(0)= (-0.5d0*cone - sw2 * ql)/sw/cw  !charged leptons
      gl(1) = -sw/cw*ql
      gn(0)= +0.5d0/sw/cw                    !neutrinos
      gn(1) = 0.d0
      gu(0)= (+0.5d0*cone - sw2 * qu)/sw/cw  !quarks
      gu(1) = -sw/cw*qu
      gd(0)= (-0.5d0*cone - sw2 * qd)/sw/cw
      gd(1) = -sw/cw*qd

      if(verbose) then
          write(*,*) '*************************************'
          write(*,*) 'Z mass = ',ph_Zmass
          write(*,*) 'Z width = ',ph_Zwidth
          write(*,*) 'W mass = ',ph_Wmass
          write(*,*) 'W width = ',ph_Wwidth
          write(*,*) '1/alphaem = ',1d0/ph_alphaem
          write(*,*) 'sthw2 = ',ph_sthw2
          write(*,*) '(unit_e)^2 = ',ph_unit_e**2   
          write(*,*) '(g_w)^2 = ',ph_unit_e*ph_unit_e/ph_sthw2   
    
          write(*,*) '*************************************'
          write(*,*)
          write(*,*) '*************************************'
          write(*,*) sqrt(ph_Zmass2low),'< M_Z <',sqrt(ph_Zmass2high)
          write(*,*) '*************************************'
      endif
      end

