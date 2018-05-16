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
         write(*,*) 'On-shell renormalization scheme selected.'
      else if (massren.eq.1) then
         write(*,*) 'MSBAR renormalization scheme selected'
      else if (massren.eq.2) then
         write(*,*) 'DRBAR renormalization scheme selected'
      else
         write(*,*) 'Unknown renormalization scheme selected.
     $Aborting.'
         stop
      endif
      
      flg_2HDMtype=int(powheginput('2HDMtype'))
      if (flg_2HDMtype.eq.1) then
        write(*,*) 'Type I 2HDM selected'
      else if (flg_2HDMtype.eq.2) then
        write(*,*) 'Type II 2HDM selected'
      else if (flg_2HDMtype.eq.3) then
        write(*,*) 'Lepton specific 2HDM selected'
      else if (flg_2HDMtype.eq.4) then
        write(*,*) 'Flipped 2HDM selected'
      else
        write(*,*) 'Error: unrecognized 2HDM type,', flg_2HDMtype
        write(*,*) 'Please check your input file.'
        write(*,*) 'Aborting.'
        stop
      endif


c     This is a common factor which actually cancel out in the MSSM amplitudes and it
c     was in the BDV formalism to make the coupling adimensional
      mmaa = 256d0


      flg_hHweightevents = int(powheginput('#flg_hHweightevents'))
      if (flg_hHweightevents.gt.0) then
         write(*,*) 'Weighted events generation enabled'
         idwtup = -4
         flg_hHweightevents = 1
      endif

      call init_ind()
      call init_fermions()
      call gen2HDMvar()
      call init_higgs()
      call init_ew()


      write(*,*) 'Active fermions'
      write(*,*) afer
      write(*,*) 'fermion masses'
      write(*,*) mfer
      write(*,*) 'fermion couplings'
      write(*,*) lambdafer
      write(*,*) 'mumass2'
      write(*,*) mumass2
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'GF = ',ph_GF
      end


c     Indipendent variables for SM
      subroutine init_ind()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

      ph_Hmass = powheginput('hmass')
      ph_Hmass2 = ph_Hmass**2

      flg_hdecay = int(powheginput('#hdecaywidth'))
      call gethiggswidth()
      ih = int(powheginput('higgstype'))
      if ((ih.ne.1).and.(ih.ne.2).and.(ih.ne.3)) then
        write(*,*) 'Unrecognized value for higgstype'
        write(*,*) '1 -> Light CP-even neutral Higgs'
        write(*,*) '2 -> Heavy CP-even neutral Higgs'
        write(*,*) '3 -> CP-odd neutral Higgs'
        stop
      endif

      call getweakbosonspars()

      ph_sthw2 = abs(1d0-(ph_Wmass/ph_Zmass)**2)
      ph_GF= powheginput('gfermi')

      ph_sthw = sqrt(ph_sthw2)
      ph_cthw2 = 1-ph_sthw2
      ph_cthw = sqrt(ph_cthw2)

c     NNLO normalization rescaling factor
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

      end subroutine

      subroutine init_ew()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

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
      end subroutine

c     If fast_ew is enabled, create the array with the campionazied deltaew value.
      subroutine init_cached_ew_corr
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      integer i,imh
      real * 8 m12,y12,value
      complex * 16 ampl, reduced,x12,aux,bornqcd
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
         ampl = dcmplx(0d0)
         do i=1,afer
            m12=mfer(i)
            y12=m12**2/mh**2
            x12 = reduced(1d0/y12)
            aux = lambdafer(i)*trfer(i)*
     &           (-4d0)*y12*(2d0-(1d0-4d0*y12)*0.5d0*log(x12)**2)
            ampl = ampl+aux
         end do
         cached_ew_corr(imh) = deltaew(ampl)
         write(33,*) 'deltaew',mh, cached_ew_corr(imh),
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
      real * 8 masswindow
      real * 8 powheginput
      external powheginput
      logical verbose
      parameter(verbose=.true.)

      flg_onshell = int(powheginput('zerowidth'))

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
         masswindow = powheginput("#masswindow")
         if(masswindow.lt.0d0) masswindow=10d0
c     ph_Hmass2low=(ph_Hmass-masswindow*ph_Hwidth)^2
         ph_Hmass2low=max(0.5d0,ph_Hmass-masswindow*ph_Hwidth)
         ph_Hmass2low= ph_Hmass2low**2
c     ph_Hmass2high=(ph_Hmass+masswindow*ph_Hwidth)^2
         ph_Hmass2high=min(kn_sbeams,(ph_Hmass+masswindow*ph_Hwidth)**2)
         ph_HmHw = ph_Hmass * ph_Hwidth
         ph_unit_e = sqrt(4*pi*ph_alphaem)
      endif
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
      if (flg_hdecay.gt.0) then
         write(*,*) 'Higgs width from hdecay'
         call hdecayparser(ph_Hmass,ph_Hwidth)
      else
         write(*,*) 'Higgs width from powheg.input'
         ph_Hwidth = powheginput('hwidth')
      endif

      end subroutine

c     2HDM model couplings
      subroutine gen2HDMvar()
      implicit none
      include 'Flags.h'
#include "PhysPars.h"
      real * 8 sina,cosa
      real * 8 powheginput
      external powheginput

      ph_alpha = powheginput('alpha')
      sina = sin(ph_alpha)
      cosa = cos(ph_alpha)

      ph_tanb = powheginput('tanb')
      ph_cb = cos(atan(ph_tanb))
      ph_sb = sin(atan(ph_tanb))

c     h higgs
      if(ih.eq.1) then
          if (flg_2HDMtype.eq.1) then
c     Type I 2HDM selected
              lambdafer(1) = cosa/ph_sb
              lambdafer(2) = cosa/ph_sb
c     Type II 2HDM selected
          else if (flg_2HDMtype.eq.2) then
              lambdafer(1) = cosa/ph_sb
              lambdafer(2) = -sina/ph_cb
c     Lepton specific 2HDM selected
          else if (flg_2HDMtype.eq.3) then
              lambdafer(1) = cosa/ph_sb
              lambdafer(2) = cosa/ph_sb
c     Flipped 2HDM selected
          else if (flg_2HDMtype.eq.4) then
              lambdafer(1) = cosa/ph_sb
              lambdafer(2) = -sina/ph_cb
          endif
c     H higgs
      else if(ih.eq.2) then
          if (flg_2HDMtype.eq.1) then
c     Type I 2HDM selected
              lambdafer(1) = sina/ph_sb
              lambdafer(2) = sina/ph_sb
c     Type II 2HDM selected
          else if (flg_2HDMtype.eq.2) then
              lambdafer(1) = sina/ph_sb
              lambdafer(2) = cosa/ph_cb
c     Lepton specific 2HDM selected
          else if (flg_2HDMtype.eq.3) then
              lambdafer(1) = sina/ph_sb
              lambdafer(2) = sina/ph_sb
c     Flipped 2HDM selected
          else if (flg_2HDMtype.eq.4) then
              lambdafer(1) = sina/ph_sb
              lambdafer(2) = -cosa/ph_cb
          endif
c     A higgs
      else
           if (flg_2HDMtype.eq.1) then
c     Type I 2HDM selected
              lambdafer(1) = 1d0/ph_tanb
              lambdafer(2) = -1d0/ph_tanb
c     Type II 2HDM selected
          else if (flg_2HDMtype.eq.2) then
              lambdafer(1) = 1d0/ph_tanb
              lambdafer(2) = ph_tanb
c     Lepton specific 2HDM selected
          else if (flg_2HDMtype.eq.3) then
              lambdafer(1) = 1d0/ph_tanb
              lambdafer(2) = -1d0/ph_tanb
c     Flipped 2HDM selected
          else if (flg_2HDMtype.eq.4) then
              lambdafer(1) = 1d0/ph_tanb
              lambdafer(2) = ph_tanb
          endif
      endif
      end subroutine


      subroutine getweakbosonspars()
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

      ph_Zmass = powheginput("#Zmass")
      if (ph_Zmass.le.0d0) ph_Zmass = 91.1876d0
      ph_Wmass = powheginput("#Wmass")
      if (ph_Wmass.le.0d0) ph_Wmass = 80.385d0
      ph_Zwidth = powheginput("#Zwidth")
      if (ph_Zwidth.le.0d0) ph_Zwidth =  2.4952d0
      ph_Wwidth = powheginput("#Wwidth")
      if (ph_Wwidth.le.0d0) ph_Wwidth =  2.085d0
      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/137.035999679d0
      ph_alphaemmz = powheginput("#alphaemmz")
      if (ph_alphaemmz.le.0d0) ph_alphaemmz = 1d0/128d0
      end subroutine


      subroutine init_fermions()
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'Flags.h'
      real * 8 powheginput,getRSmass,muren2
      real * 8 pwhg_alphas,alpha_s
      external pwhg_alphas
      external powheginput

c     TODO: for the final release this should be changed
      muren2 = ph_Hmass2
      alpha_s =  pwhg_alphas(ph_Hmass2,st_lambda5MSB,st_nlight)
      write(*,*) 'lambda', st_lambda5MSB, st_nlight

c     For the moment we force both top and bottom present
      afer = 0
      ph_topmass = powheginput('topmass')
      if (ph_topmass.ne.-1000000d0) then
         write(*,*) "Top Quark enabled"
         afer = afer + 1
         trfer(afer) = 1d0/2d0
         mfer(afer) = getRSmass(ph_topmass,massren,ph_Hmass,alpha_s)
         lambdafer(afer) = 1d0
         ferlogmratio(afer) = log(mfer(afer)**2/muren2)
         if (massren.eq.2) then
            ferlogmratio(afer) = ferlogmratio(afer) - 1d0/3d0
         endif
      endif
      ph_bottommass = powheginput('#bottommass')
      if (ph_bottommass.ne.-1000000d0) then
         write(*,*) "Bottom Quark enabled"
         afer = afer+1
         trfer(afer) = 1d0/2d0
         mfer(afer) = getRSmass(ph_bottommass,massren,ph_Hmass,alpha_s)
         lambdafer(afer) = 1d0
         ferlogmratio(afer) = log(mfer(afer)**2/muren2)
         if (massren.eq.2) then
            ferlogmratio(afer) = ferlogmratio(afer) - 1d0/3d0
         endif
      endif

      end subroutine


      real * 8 function getRSmass(osmass,RS,q,alpha_s)
      implicit none
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 osmass,q,alpha_s
      integer RS

c     OS
      if (RS.eq.0) then
         getRSmass = osmass
c     MSBAR
      else if (RS.eq.1) then
         getRSmass = osmass*(1-alpha_s/pi*(log(q**2/osmass**2)
     $                  +4d0/3d0))
c     DRBAR
      else if (RS.eq.2) then
         getRSmass = osmass*(1-alpha_s/pi*(log(q**2/osmass**2)
     $                  +5d0/3d0))
      else
         write(*,*) 'Error: unknown renormalization scheme! Exiting'
         stop
      endif

c      write(*,*) 'st_alpha', alpha_s

      end function
