      subroutine init_couplings
      implicit none
      include 'LesHouches.h'
      include 'PhysPars.h'
      include 'Flags.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'    
      include 'nlegborn.h'      
      include 'pwhg_kn.h'      
      real * 8 powheginput
      external powheginput
      real * 8 ewamplitude,amplew
      common /ew/ewamplitude
      external amplew

c     number of light flavors
      st_nlight = 5

      massren = int(powheginput('scheme'))
      if (massren.eq.0) then
         massrenb = 0
         write(*,*) 'Unknow renormalization scheme selected. Abor
     $ing.'
c     Original OS scheme disabled
         stop
         write(*,*) 'On-shell renormalization scheme selected.'
      else if (massren.eq.1) then
         massrenb = 1
         write(*,*) 'DRBar renormalization scheme selected'
      else if (massren.eq.2) then
         massrenb = 2
         write(*,*) 'Mixed renormalization scheme selected'
      else
         write(*,*) 'Unknown renormalization scheme selected.
     $Aborting.'
         stop
      endif
      
      write(*,*) 'MSSM selected'
      afer = 2
      flg_noscalars = int(powheginput('#nosquarks'))
      if (flg_noscalars.le.0) then
         asca = 4
         flg_lhscalars = int(powheginput('#lhsquarks'))
         if (flg_lhscalars.gt.0) then
            write(*,*) 'Light Higgs limit for scalars enabled'
            flg_lhscalars = 1
         else
            flg_lhscalars = 0
         endif
      else
         asca = 0
      endif
c     This is a common factor which actually cancel out in the MSSM amplitudes and it
c     was in the BDV formalism to make the coupling adimensional
      mmaa = 256d0
      ih = int(powheginput('higgstype'))
      if ((ih.ne.1).and.(ih.ne.2).and.(ih.ne.3)) then
        write(*,*) 'Unrecognized value for higgstype'
        write(*,*) '1 -> Light CP-even neutral Higgs'
        write(*,*) '2 -> Heavy CP-even neutral Higgs'
        write(*,*) '3 -> CP-odd neutral Higgs'
        stop
      endif

!      flg_hHweightevents = int(powheginput('#flg_hHweightevents'))
!      if (flg_hHweightevents.gt.0) then
!         write(*,*) 'Weighted events generation enabled'
!         idwtup = -4
!         flg_hHweightevents = 1
!      endif


      call genmssmvar()
      call init_higgs()

c     EW corrections
      flg_ew = int(powheginput('ew'))

      if (flg_ew.le.0) then
         flg_ew = 0
         write(*,*) '2-loops EW corrections disabled'
      else
         if(ih.eq.3) then
            write (*,*) 'Warning: there no EW corrections for the pseudo
     $scalar. Disabling. Please check your input file'
            flg_ew = 0
         else
            flg_ew = 1
            write(*,*) '2-loops EW corrections enabled'

            flg_fast_ew = int(powheginput('#fastew'))
            if (flg_fast_ew.gt.0) then
             write (*,*) 'Enable fast ew corrections evaluations
     $ by mass sampling'
              flg_fast_ew = 1
             call init_cached_ew_corr
            endif
         endif
      endif

      flg_nnlo = int(powheginput('#nnlo'))
      if (flg_nnlo.le.0) then
         flg_nnlo = 0
         write(*,*) 'NNLO rescaling disabled'
         nnlorescfactor = 1d0
      else
         open(24,file='nnlo.dat')
         read (24,*) nnlorescfactor
         write(*,*) 'NNLO rescaling factor set to', nnlorescfactor
      endif

      write(*,*) 'Active fermions'
      write(*,*) afer
      write(*,*) 'fermion masses'
      write(*,*) mfer
      write(*,*) 'fermion couplings'
      write(*,*) lambdafer
      write(*,*) 'Active scalars'
      write(*,*) asca
      write(*,*) 'Scalar masses'
      write(*,*) msca
!     These are only valid for h/H
      if(ih.ne.3) then
        write(*,*) 'Scalar couplings'
         write(*,*) lambdasca
      endif
      write(*,*) 'mumass2'
      write(*,*) mumass2
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'GF = ',ph_GF

      end

c     If fast_ew is enabled, create the array with the campionazied deltaew value.
      subroutine init_cached_ew_corr
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      integer i,imh
      real * 8 m12,m0,y12,y0,value
      complex * 16 ampl, reduced,x12,x0,aux,bornqcd
      common /bornampl/ampl,bornqcd
      real * 8 deltaew
      external reduced,deltaew

      write(*,*) 'Starting cached EW corrections array initialization'
      open(unit=33,file="ew.dat")
c     Cannot loop on real, see fortran standards
      do imh=0,2500
c first bin caveat
        if(imh.eq.0) then
            mh = 0.5d0
        else
            mh=REAL(imh)
        endif
c        mh=REAL(imh)
         value = 0d0
c      write(*,*) ''
c      write(*,*) 'imh',imh
c      write(*,*) 'mh', mh
c      write(*,*) 'nint(mh)', nint(mh)
         ampl = dcmplx(0d0)
         do i=1,afer
            m12=mfer(i)
            y12=m12**2/mh**2
            x12 = reduced(1d0/y12)
            aux = lambdafer(i)*trfer(i)*
     &           (-4d0)*y12*(2d0-(1d0-4d0*y12)*0.5d0*log(x12)**2)
            ampl = ampl+aux
         end do
         if (flg_lhscalars.eq.0) then
            do i=1,asca
               m0 = msca(i)
               y0 = m0**2/mh**2
               x0 = reduced(1d0/y0)
               ampl = ampl + lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &              4d0*y0*(1d0+2d0*y0*0.5d0*log(x0)**2)
            end do
         else
            do i=1,asca
               m0 = msca(i)
               y0 = m0**2/mh**2
               ampl = ampl + lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &              (-1d0/3d0-8d0/(45d0*4d0*y0))
            end do
         endif
         cached_ew_corr(imh) = deltaew(ampl)
         write(33,*) 'deltaew',mh, imh, cached_ew_corr(imh),
     &REAL(ampl),AIMAG(ampl)
      end do
      close(33)
      write(*,*) 'Finished cached EW corrections array initialization'
      end subroutine init_cached_ew_corr

c     Initialize Higgs mass
      subroutine init_higgs()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      real * 8 masswindow,masswindow_high,masswindow_low
      real * 8 powheginput
      external powheginput
      logical verbose
      parameter(verbose=.true.)

      flg_onshell = int(powheginput('zerowidth'))
      call gethiggswidth()
      if (flg_onshell.eq.1) then
c     Higgs is produced on-shell
         ph_Hmass2low = ph_Hmass2
         ph_Hmass2high = ph_Hmass2low
         ph_HmHw = 0d0
      else


c     set mass windows around H-mass peak in unit of ph_Hwidth
c     It is used in the generation of the Born phase space
c     masswindow is an optonal  parameter passed by the user
c     the default vale is 10
c     Keep compatibility with old input file
      masswindow = powheginput("#masswindow")
      if(masswindow.gt.0d0) then
        masswindow_low = masswindow
        masswindow_high = masswindow
      else
        masswindow_low = powheginput("#masswindow_low")
        if(masswindow_low.lt.0d0) masswindow_low=10d0
        masswindow_high = powheginput("#masswindow_high")
        if(masswindow_high.lt.0d0) masswindow_high=10d0
      endif

      ph_Hmass2low=max(0.5d0,ph_Hmass-masswindow_low*ph_Hwidth)
      ph_Hmass2low= ph_Hmass2low**2
      ph_Hmass2high=ph_Hmass+masswindow_high*ph_Hwidth
      ph_Hmass2high= min(kn_sbeams,ph_Hmass2high**2)
      ph_HmHw = ph_Hmass * ph_Hwidth
      endif

      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(verbose) then
         write(*,*) '*************************************'
         write(*,*) 'H mass = ',ph_Hmass
         if (flg_onshell.ne.1) then
            write(*,*) 'H width = ',ph_Hwidth
         else
            write(*,*) 'H is produced on-shell'
         endif
         if (flg_onshell.ne.1) then
            write(*,*) '*************************************'
            write(*,*) sqrt(ph_Hmass2low),' < M_H <',sqrt(ph_Hmass2high)
         endif
         write(*,*) '*************************************'
      endif

      end subroutine


      subroutine gethiggswidth()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

      flg_hdecay = int(powheginput('#hdecaywidth'))
      flg_fhdecay = int(powheginput('#fhdecaywidth'))
      if ((flg_hdecay.gt.0).and.(flg_fhdecay.le.0)) then
         write(*,*) 'Higgs width from hdecay'
         call hdecayparser(ph_Hmass,ph_Hwidth)
      else if ((flg_hdecay.le.0).and.(flg_fhdecay.gt.0)) then
                write(*,*) 'Higgs width from FeynHiggs'
                ph_Hwidth = fh_Hwidth
      else if ((flg_hdecay.gt.0).and.(flg_fhdecay.gt.0)) then
           write(*,*) 'Error: conflicting options for the decay width.'
           write(*,*) '       Only one of hdecaywidth and hdecaywidth'
           write(*,*) '       can be enabled at the same time. Aborting'
           stop
      else
         write(*,*) 'Higgs width from powheg.input'
         ph_Hwidth = powheginput('hwidth')
      endif

      end subroutine

      subroutine getweakbosonswidth()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

      ph_Zwidth = powheginput("#Zwidth")
      if (ph_Zwidth.le.0d0) ph_Zwidth =  2.4952d0
      ph_Wwidth = powheginput("#Wwidth")
      if (ph_Wwidth.le.0d0) ph_Wwidth =  2.085d0

      end subroutine


c     MSSM DRBAR
      subroutine genmssmvar()
      implicit none
      include 'Flags.h'
#include "PhysPars.h"
#include "pwhg_st.h"
      real * 8 qtop,i3top,lambdat,lambdab
      real * 8 qbot,i3bot,hmixfact1,hmixfact2
      real * 8 c2b,c2t,cosa,s2b,s2t,sina,sin2b
      real * 8 Xb,Xt,fact1,fact2
      real * 8 powheginput
      external powheginput

c     We fill the common block with the physical parameters according
c     to the chosen renormalization scheme.
      if ((massren.eq.0).or.(massren.eq.2)) then
         call OSsusyparameter()
      else if (massren.eq.1) then
         call DRBARsusyparameter()
      else
         write(*,*) 'Error: unknown mass renormalization scheme seelcted
     $'
      end if

      s2t = 2d0*ph_ctt*ph_stt
      c2t = ph_ctt**2 - ph_stt**2
      s2b = 2d0*ph_ctb*ph_stb
      c2b = ph_ctb**2 - ph_stb**2
      if (massren.eq.1) then
         Xt = (ph_t1_2-ph_t2_2)*s2t/2d0/ph_topmass
         At = Xt - ph_mumssm/ph_tanb           ! note our sign convention for mu
         Xb = (ph_b1_2-ph_b2_2)*s2b/2d0/ph_mbsb
         Ab = Xb - ph_mumssm*ph_tanb           ! note our sign convention for m
      endif
c     Couplings (lambdas) calculation
c     Physcal parameters needed to calculate lambdafer e lambdasca
      qtop = 2d0/3d0
      i3top = 1d0/2d0
      qbot = -1d0/3d0
      i3bot = -1d0/2d0

      sina = sin(ph_alpha)
      cosa = cos(ph_alpha)

c     sin(2beta)
      sin2b = 2d0 * ph_sb * ph_cb
      if ((ih.eq.1).or.(ih.eq.2)) then
        lambdat = 1d0/ph_sb
        lambdab = 1d0/ph_cb
c     pseudoscalar
      else
        lambdat = 1/ph_tanb
        lambdab = ph_tanb
      endif

c     Light Higgs
      if (ih.eq.1) then
         hmixfact1 = -sina
         hmixfact2 = cosa
c     Heavy Higgs
      else if (ih.eq.2) then
         hmixfact1 = cosa
         hmixfact2 = sina
c     Pseudoscalar
      else
         hmixfact1 = 1d0
         hmixfact2 = 1d0
      endif

c     fermion 1  (top)
      trfer(1) = 1d0/2d0
      mfer(1) = ph_topmass
      lambdafer(1) = lambdat * hmixfact2
      if (massren.eq.1) then
         ferlogmratio(1) = log(mfer(1)**2/q**2)-1d0/3d0
      end if

c     fermion 2  (bottom)
      trfer(2) = 1d0/2d0
      if (flg_mssm_q_mh.eq.1) then
         mfer(2) = ph_bottommass_mh
      else
         mfer(2) = ph_bottommass
      end if
      lambdafer(2) = lambdab * hmixfact1 * botcouplfac
      if (massren.eq.1) then
         if (flg_mssm_q_mh.eq.1) then
            ferlogmratio(2) = log(mfer(2)**2/ph_Hmass**2)-1d0/3d0
         else
            ferlogmratio(2) = log(mfer(2)**2/q**2)-1d0/3d0
         end if
      end if

c     scalar 1  (stop 1)
      trsca(1) = 1d0/2d0
      msca(1) = ph_t1

      fact1 =   1d0/2d0 * s2t * ph_mumssm * mfer(1)
     $        + 1d0/4d0 * i3top * ph_Zmass**2 * sin2b * (1+c2t)
     $        - 1d0/2d0 * c2t * ph_Zmass**2 * qtop * sin2b * ph_sthw2

      fact2 =   mfer(1)**2
     $        + 1d0/2d0 * At * mfer(1) * s2t
     $        - 1d0/2d0 * i3top * ph_Zmass**2 * ph_sb**2 * (1+c2t)
     $        + c2t * ph_Zmass**2 * qtop * ph_sb**2 * ph_sthw2

      lambdasca(1) = (lambdat/mmaa**2) * (
     $hmixfact1 * fact1 +
     $hmixfact2 * fact2
     $)

c     scalar 2  (stop 2)
      trsca(2) = 1d0/2d0
      msca(2) = ph_t2

      fact1 = - 1d0/2d0 * s2t * ph_mumssm * mfer(1)
     $        + 1d0/4d0 * i3top * ph_Zmass**2 * sin2b * (1-c2t)
     $        + 1d0/2d0 * c2t * ph_Zmass**2 * qtop * sin2b * ph_sthw2

      fact2 =  mfer(1)**2
     $       - 1d0/2d0 * At * mfer(1) * s2t
     $       - 1d0/2d0 * i3top * ph_Zmass**2 * ph_sb**2 * (1-c2t)
     $       - c2t * ph_Zmass**2 * qtop * ph_sb**2 * ph_sthw2

      lambdasca(2) = (lambdat/mmaa**2) * (
     $hmixfact1 * fact1 +
     $hmixfact2 * fact2
     $)

c     scalar 1  (sbottom 1)
      trsca(3) = 1d0/2d0
      msca(3) = ph_b1

      fact1 =   ph_mbsb**2
     $        + 1d0/2d0 * i3bot * ph_Zmass**2 * ph_cb**2
     $        + 1d0/2d0 * i3bot * ph_Zmass**2 * ph_cb**2 * c2b
     $        + 1d0/2d0 * Ab * ph_mbsb * s2b
     $        - c2b * ph_cb**2 * ph_Zmass**2 * qbot * ph_sthw2

      fact2 =   1d0/2d0 * s2b * ph_mumssm * ph_mbsb
     $        - 1d0/4d0 * i3bot * ph_Zmass**2 * sin2b * (1+c2b)
     $        + 1d0/2d0 * c2b * ph_Zmass**2 * qbot * sin2b * ph_sthw2

      lambdasca(3) = (lambdab/mmaa**2) * (
     $hmixfact1 * fact1 +
     $hmixfact2 * fact2
     $)

c     scalar 2  (sbottom 2)
      trsca(4) = 1d0/2d0
      msca(4) = ph_b2

      fact1 =   ph_mbsb**2 
     $        - 1d0/2d0 * Ab * ph_mbsb * s2b
     $        + 1d0/2d0 * ph_cb**2 * i3bot * ph_Zmass**2 * (1-c2b)
     $        + c2b * ph_cb**2 * ph_Zmass**2 * qbot * ph_sthw2

      fact2 = - 1d0/2d0 * s2b * ph_mumssm * ph_mbsb
     $        - 1d0/4d0 * i3bot * ph_Zmass**2 * sin2b * (1-c2b)
     $        - 1d0/2d0 * c2b * ph_Zmass**2 * qbot * sin2b * ph_sthw2

      lambdasca(4) = (lambdab/mmaa**2) * (
     $hmixfact1 * fact1 +
     $hmixfact2 * fact2
     $)

      end

      subroutine DRBARsusyparameter()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 v,mh_1,mh_2,alpha
      real * 8 g,gp,msbr,msql,mstr
      real * 8 q_g,q_yb,q_yt,qsoft,q_ad,q_au
      real * 8 yb,yt
      real * 8 checkvalue
      real * 8 mstop2(2),msbot2(2)
      real * 8 powheginput
      real * 8 q_yb_mh,v_mh,yb_mh
      real * 8 ma_slha
      external powheginput
      character*50 valname
      integer error
c     SLHALib definition
#include "SLHA.h"
      complex * 16 slhadata(nslhadata)
c     First implementation of quarks mass renormalised at mh scale.
c     These are the only two model parameter we actually retrive from the powheg.input file
      ph_GF= powheginput('gfermi')

c     The others are read from the SLHA file
      call SLHAclear(slhadata)
      call SLHARead(error, slhadata, "mssm-param.slha", 0)
      if(error.ne.0) stop "Read error from mssm-param.slha"
      write(*,*) "** SLHA file input parameters **"
c     CP-Even neutral Higgs double masses
      valname = "Light CP-even Higgs mass"
      mh_1 = checkvalue(Mass_Mh0,valname)
      valname = "Heavy CP-even Higgs mass"
      mh_2 = checkvalue(Mass_MHH,valname)
c     CP-ODD neutral Higgs mass
      valname = 'CP-odd Higgs mass'
      ph_ma_pole = checkvalue(Mass_MA0,valname)
c     Higgsino mass parameters
      valname = "Higgsino mass parameter"
      ph_mumssm = checkvalue(HMix_MUE,valname)
c     Tan beta
      valname = "Tan(beta)"
      ph_tanb = checkvalue(HMix_TB,valname)
c     Vev
      valname = "Higgs vev"
      v = checkvalue(HMix_VEV,valname)
c     MSOFT Block
      valname = "Soft Block renormalization scale"
      qsoft = checkvalue(MSoft_Q,valname)
      valname = "Gluino mass"
      ph_mg  = checkvalue(MSoft_M3,valname)
      valname = "mqL3"
      msql = checkvalue(MSoft_MSQ(3),valname)
      valname = "mtR"
      mstr = checkvalue(MSoft_MSU(3),valname)
      valname = "mbR"
      msbr = checkvalue(MSoft_MSD(3),valname)
c     AU Block
      valname = "Au block renormalization scale"
      q_au = checkvalue(Au_Q,valname)
      valname = "At"
      At   = checkvalue(Au_At,valname)
c     AD Block
      valname = "Ad block renormalization scale"
      q_ad = checkvalue(Ad_Q,valname)
      valname = "Ab"
      Ab   = checkvalue(Ad_Ab,valname)
c     YU Block
      valname = "Yu block renormalization scale"
      q_yt = checkvalue(Yu_Q,valname)
      valname = "Yt"
      yt = checkvalue(Yu_Yt,valname)
c     YD Block
      valname = "Yd block renormazation scale"
      q_yb = checkvalue(Yd_Q,valname)
      valname = "Yb"
      yb = checkvalue(Yd_Yb,valname)
c     Mixing angle in the CP-Even Higgs sector
      valname = "CP-even Higgs sector mixing angle"
      alpha = checkvalue(Alpha_Alpha,valname)
c     GAUGE Block
      valname = "Gauge block renormalization scale"
      q_g = checkvalue(Gauge_Q,valname)
      valname = "g'"
      gp  = checkvalue(Gauge_g1,valname)
      valname = "g "
      g  = checkvalue(Gauge_g2,valname)

      if ((qsoft.eq.q_au).and.(q_au.eq.q_ad).and.(q_ad.eq.q_yt)
     $     .and.(q_yt.eq.q_yb).and.(q_yb.eq.q_g)) then
         q = qsoft
      else
         write(*,*) 'Error: inconsistency in the
     $renormalization scales. Aborting.'
         stop
      endif
      ph_mumssm = -ph_mumssm                  ! we use the opposite convention

c     DR-BAR
c     compute the running masses for quarks, squarks and gauge bosons

      write(*,*) 'Computing the running masses for quarks, squarks
     $ and gauge bosons'

      ph_cb = 1d0/sqrt(1+ph_tanb**2)
      ph_sb = ph_tanb*ph_cb

      ph_topmass = yt*v/sqrt(2d0)*ph_sb
      ph_bottommass = yb*v/sqrt(2d0)*ph_cb

      ph_Wmass = g*v/2d0
      ph_Zmass = sqrt(g**2+gp**2)*v/2d0
      call getweakbosonswidth()


      ph_sthw2 = 1-(ph_Wmass/ph_Zmass)**2
      ph_alphaem = g**2*gp**2/(g**2+gp**2)/16d0/atan(1d0) ! note it's at Q=MS!!!

      write(*,*) 'Running top mass ', ph_topmass
      write(*,*) 'Running bottom mass ', ph_bottommass
      write(*,*) 'Running Z mass ', ph_Zmass
      write(*,*) 'Running W mass ', ph_Wmass
      write(*,*) 'Running sin(thetaw)', ph_sthw2

c     compute the squark masses

      call diagonalize(ph_topmass,ph_Wmass,ph_Zmass,msql,mstr,At,
     $     ph_mumssm,ph_tanb,1,0,mstop2,ph_stt,ph_ctt)

      ph_t1 = sqrt(mstop2(1))
      ph_t1_2 = mstop2(1)
      ph_t2 = sqrt(mstop2(2))
      ph_t2_2 = mstop2(2)

      write(*,*) 'stop masses and sin(tht)',ph_t1,ph_t2,ph_stt

      call diagonalize(ph_bottommass,ph_Wmass,ph_Zmass,msql,msbr,Ab,
     $     ph_mumssm,ph_tanb,2,0,msbot2,ph_stb,ph_ctb)

      ph_b1 = sqrt(msbot2(1))
      ph_b1_2 = msbot2(1)
      ph_b2 = sqrt(msbot2(2))
      ph_b2_2 = msbot2(2)

      write(*,*) 'sbottom masses and sin(thb)',ph_b1,ph_b2,ph_stb

      mumass2 = q**2
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw2 = 1-ph_sthw2
      ph_cthw = sqrt(ph_cthw2)

      ph_alpha = alpha
c     In the DRBAR scheme we fix these to the DRBAR values
      ph_mbsb = ph_bottommass
      ph_mbpole = ph_bottommass
      ph_mbrun = ph_bottommass
      botcouplfac = 1d0

      if (ih.eq.1) then
         write(*,*) 'Light neutral CP-even Higgs selected'
         ph_Hmass = mh_1
      elseif (ih.eq.2) then
         write(*,*) 'Heavy neutral CP-even Higgs selected'
         ph_Hmass = mh_2
      else
         write(*,*) 'Neutral CP-odd Higgs selected'
         ph_Hmass = ph_ma_pole
      endif
      ph_Hmass2 = ph_Hmass**2

      flg_mssm_q_mh = int(powheginput('#mssmmbatmh'))
      if (flg_mssm_q_mh.eq.1) then
         call SLHAclear(slhadata)
         call SLHARead(error, slhadata, "mssm-param-mh.slha", 0)
         if(error.ne.0) stop "Read error from mssm-param-mh.slha"
         write(*,*) "** SLHA file input parameters at scale MH **"
c     Tan beta@mh
         valname = "Tan(beta)@mh"
         ph_tanb_mh = checkvalue(HMix_TB,valname)
c     Vev
         valname = "Higgs vev@mh"
         v_mh = checkvalue(HMix_VEV,valname)
c     YD Block
         valname = "Yd block renormazation scale == mh?"
         q_yb_mh = checkvalue(Yd_Q,valname)
         valname = "Yb@mh"
         yb_mh = checkvalue(Yd_Yb,valname)

         ph_cb_mh = 1d0/sqrt(1+ph_tanb_mh**2)

         ph_bottommass_mh = yb_mh*v_mh/sqrt(2d0)*ph_cb_mh
         write(*,*) 'mb(mh) is', ph_bottommass_mh
      else
         ph_bottommass_mh = 0d0
      end if

      end

      real * 8 function checkvalue(val,valname)
      complex * 16 val
      character*50 valname
c     SLHALib definition
#include "SLHA.h"
      if (val.eq.invalid) then
         write(*,*) 'Invalid/Absent value in SLHA file for ', valname
         stop
      endif
      checkvalue = REAL(val)
      write(*,*) valname,checkvalue
      end

c     MSSM OS
c     This subroutine read the SUSY input parameters required for OS calculation
      subroutine OSsusyparameter()
      implicit none
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'Flags.h'
      real * 8 M2,M3,mstl,msbr,mstr,msf,Af,asq
      real * 8 mstop2(2),msbot2(2)
      real * 8 mh_int(2),v,asb
c     TODO: REMOVE
      real * 8 mbpole,mtrun,mbrun,msbl
      real * 8 powheginput,asf
      real * 8 pwhg_alphas
      external pwhg_alphas
      external powheginput
      external asf

c     From Feynhiggs we get the everything
      call runFH(ph_ma_mz,M2,ph_mg,mstl,mstr,msbr,msf,Af,
     $ ph_mumssm,ph_tanb_mz,mh_int,ph_alpha)
      call getweakbosonswidth()
c     Compute various physical quantities
      v = 1d0/sqrt(sqrt(2d0)*ph_GF)
      ph_sthw2 = 1-(ph_Wmass/ph_Zmass)**2
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw2 = 1-ph_sthw2
      ph_cthw = sqrt(ph_cthw2)

      ph_alphaem = ph_Wmass**2*ph_sthw2/v**2/pi
      ph_alphaemmz = ph_alphaem


      asb = pwhg_alphas(ph_mbmb,st_lambda5MSB,st_nlight)
      ph_bottommass = powheginput('#bottommass')
      if(ph_bottommass.le.0) then
        write(*,*) 'Calculating mb(OS)'
        write(*,*) 'alphas(mb) = ',asb
        ph_bottommass = ph_mbmb*(1 + 4d0/3d0*asb/pi)
      endif
c      ph_bottommass = ph_mbmb*(1 + 4d0/3d0*asb/pi + 9.2778d0*(asb/pi)**2
c     $     + 94.4182d0*(asb/pi)**3)

c      ph_bottommass = 4.75d0

      write(*,*) 'mb(mb) = ',ph_mbmb
      write(*,*) 'mb_pole = ', ph_bottommass

      ph_tanb = ph_tanb_mz
      ph_cb = 1d0/sqrt(1+ph_tanb**2)
      ph_sb = ph_tanb*ph_cb

      ph_mbsb = sqrt(ph_ctb**2*ph_b1_2 + ph_stb**2*ph_b2_2 - msbl2
     $     - ph_Zmass**2*cos(2.d0*atan(ph_tanb))*
     $     (-1/2d0+1/3d0*ph_sthw2))

      write(*,*) 'mbsb',ph_mbsb
      write(*,*) 'MH from FeynHiggs', mh_int
      write(*,*) 'alpha from FeynHiggs', ph_alpha

      if (ih.eq.1) then
         write(*,*) 'Light neutral CP-even Higgs selected'
         ph_Hmass = mh_int(1)
      elseif (ih.eq.2) then
         write(*,*) 'Heavy neutral CP-even Higgs selected'
         ph_Hmass = mh_int(2)
      else
        write(*,*) 'Neutral CP-odd Higgs selected'
         ph_Hmass = ph_ma_mz
         ph_ma_pole = ph_ma_mz
      end if
      ph_Hmass2 = ph_Hmass**2

      q = ph_topmass
      mtrun = ph_topmass

      botcouplfac = 1d0

      if(massren.eq.2) then
         ph_mbrun = ph_bottommass/(1+ph_deltab) ! NOTE: use "effective" mass
         massrenb = 4                ! force "gluglu_bot" routine
c         OS = 1                 ! force "twoloop" routine

c$$$         if(i_h.eq.1) then
c$$$            mbrun = mbpole/(1+deltab)*(1-deltab/tan(alpha)/tanb)
c$$$         else
c$$$            mbrun = mbpole/(1+deltab)*(1+deltab*tan(alpha)/tanb)
c$$$         endif
c$$$         OSb = 5                ! force "gluglu_bot" routine
c$$$         OS = 1                 ! force "twoloop" routine
         botcouplfac = ph_mbrun/ph_bottommass
         write(*,*) 'mbrun(mh)=',ph_mbrun,' Delta_=',ph_deltab,'fac=',
     $ botcouplfac
      ph_bottommass_mh=ph_mbrun
      endif

      end subroutine


      subroutine diagonalize(mq,mw,mz,msql,msqr,Aq,mu,tb,iq,lim,
     $     msq2,sth,cth)

      implicit none

      integer iq,lim
      double precision mq,mw,mz,msql,msqr,Aq,mu,tb,msq2(2),sth,cth
      double precision mq2,mz2,mw2,c2b,xq,yq,zq,tth,dx,dy

      mq2 = mq**2
      mz2 = mz**2
      mw2 = mw**2
      c2b = (1-tb**2)/(1+tb**2)

      if(lim.eq.0) then
         if(iq.eq.1) then
            zq = mq*(Aq+mu/tb)
            dx = 1d0/4d0*mz2*c2b
            dy = 1d0/12d0*(8*mw2-5*mz2)*c2b
         elseif(iq.eq.2) then
            zq = mq*(Aq+mu*tb)
            dx = -1d0/4d0*mz2*c2b
            dy = -1d0/12d0*(4*mw2-mz2)*c2b 
         else
            write(*,*) 'ERROR: iq out of range'
         endif
      else
c$$$         if(iq.eq.1) then
c$$$            zq = mq*Aq
c$$$            dx = -1d0/4d0*mz2
c$$$            dy = -1d0/12d0*(8*mw2-5*mz2)
c$$$         elseif(iq.eq.2) then
c$$$            zq = mq*mu*sqrt(1+tb**2)
c$$$            dx = 1d0/4d0*mz2
c$$$            dy = 1d0/12d0*(4*mw2-mz2)
c$$$         else
c$$$            write(*,*) 'ERROR: iq out of range'
c$$$         endif
         if(iq.eq.1) then
            zq = mq*Aq
            dx = 0d0!-1d0/4d0*mz2
            dy = 0d0!-1d0/12d0*(8*mw2-5*mz2)
         elseif(iq.eq.2) then
            mq2 = 0d0
            zq = mq*mu*sqrt(1+tb**2)
            dx = 0d0!1d0/4d0*mz2
            dy = 0d0!1d0/12d0*(4*mw2-mz2)
         else
            write(*,*) 'ERROR: iq out of range'
         endif
      endif

      xq = mq2 + 1d0/2d0*(msql**2+msqr**2) + dx
      yq = 1d0/2d0*(msql**2-msqr**2) + dy

      msq2(1) = xq + sqrt(yq**2+zq**2)
      msq2(2) = xq - sqrt(yq**2+zq**2)

      if(zq.eq.0.and.yq.ge.0)then
         cth=1d0
         sth=0d0
      elseif(zq.eq.0.and.yq.lt.0)then
         cth=0d0
         sth=1d0
      else
         tth = 1/zq*(sqrt(yq**2+zq**2)-yq)
         cth = 1/sqrt(1+tth**2)
         sth = cth*tth
      endif

      return
      end
