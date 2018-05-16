c     input of physical parameters
c     2012-07 Andreas von Manteuffel


c     running alpha_em at a scale Q2
c     this is an on-shell quantity, defined to include 1-loop running
c     for 5 quarks and 3 charged leptons in small mass limit
      real *8 function alemeff(Q2)
      implicit none
#include "PhysPars.h"
      include 'pwhg_math.h'
      real *8 powheginput
      external powheginput
      real *8 dalq2,Q2
      if (ph_runningemscale) then
        dalq2 = (ph_alphaem0/(3*pi))
     &            *(3 * 1 + 3 * 1/3d0 + 2 * 4/3d0)*log(q2/ph_Zmass2)
     &        + ph_dalphaemmz
        alemeff = ph_alphaem0 / (1 - dalq2)
      else
        alemeff = ph_alphaem0
      endif
      end



c     extracts mass and mixing matrix row
c     from slhadata according to particle id,
c     assumes slhadata to be filled already
      subroutine set_slepton_params(id,mass,Uxl,Uxr)
      implicit none
#include "PhysPars.h"
      integer id
      real *8 mass, norm
      complex *16 Uxl, Uxr

      ! invalid default values
      mass = invalid
      Uxl = DCMPLX(invalid, invalid)
      Uxr = DCMPLX(invalid, invalid)
      ! extract parameters from slhadata
      select case(abs(id))
      case (1000011) ! selectron left
        mass = dreal(Mass_MSf(1,2,1))
        Uxl = 1d0
        Uxr = 0d0
      case (1000013) ! smuon left
        mass = dreal(Mass_MSf(1,2,2))
        Uxl = 1d0
        Uxr = 0d0
      case (1000015) ! stau 1
        mass = dreal(Mass_MSf(1,2,3))
        Uxl = StauMix_USf(1,1)
        Uxr = StauMix_USf(1,2)
      case (2000011) ! selectron right
        mass = dreal(Mass_MSf(2,2,1))
        Uxl = 0d0
        Uxr = 1d0
      case (2000013) ! smuon right
        mass = dreal(Mass_MSf(2,2,2))
        Uxl = 0d0
        Uxr = 1d0
      case (2000015) ! stau 2
        mass = dreal(Mass_MSf(2,2,3))
        Uxl = StauMix_USf(2,1)
        Uxr = StauMix_USf(2,2)
      case default
        write(*,*) 'encountered unknown slepton ID ', id
        stop
      end select

      ! check that parameters are set correctly
      if ((mass.eq.invalid).or.
     &    (dreal(Uxl).eq.invalid).or.
     &    (aimag(Uxl).eq.invalid).or.
     &    (dreal(Uxr).eq.invalid).or.
     &    (aimag(Uxr).eq.invalid)) then
        stop "mass or mixing parameters are not set"
      endif
      ! check that mixing matrix row has norm 1
      norm = sqrt( abs(Uxl)**2 + abs(Uxr)**2 )
      if ((norm.lt.0.999d0).or.(norm.gt.1.001d0)) then
         write(*,*) 'slepton mixing matrix not unitary'
         write(*,*) 'check input values'
         stop
      endif

      end



c     setup SM and MSSM parameters in PhysPars.h
c     requires encoded process ID in lprup(1)
c     please check also setup_PYTHIA_parameters if you modify this
      subroutine init_couplings
      implicit none
#include "PhysPars.h"
      include 'LesHouches.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      real *8 powheginput
      external powheginput
      real *8 almzms,dalmzms,s2thw,mw2rough,mz2
      real *8 almz,alemeff
      logical verbose
      parameter (verbose=.true.)
      ! for reading parameters from SHLA file
      integer error
      character *11 slhafilename
      data slhafilename/"params.slha"/
      ! IDs
      integer id_sla, id_slb

      ph_alphaem0 = invalid
      ph_dalphaemmz = invalid
      ph_Zmass = invalid
      ph_Zwidth = invalid
      ph_Wmass = invalid
      ph_Wwidth = invalid
      ph_sthw = invalid
      ph_cthw = invalid
      ph_gf = invalid
      ph_Zmass2 = invalid
      ph_Wmass2 = invalid
      ph_sthw2 = invalid
      ph_runningemscale = .false.

      ! hard coded input parameters
      ph_alphaem0 = 1d0/137.035999074 ! alpha_em^ONS(0) PDG 2012
      ph_Zwidth = 0d0
      ph_Wwidth = 0d0

      ! read in parameters from SLHA file
      call SLHARead(error, slhadata, slhafilename, 1)
      if (error.ne.0)  stop "error reading parameters from SLHA file"
      ! SM parameters from SLHA
      write(*,*) "note: ignoring alpha_s value given in SLHA file"
      if ((SMInputs_MZ.ne.invalid).and.
     &    (SMInputs_invAlfaMZ.ne.invalid).and.
     &    (SMInputs_GF.ne.invalid).and.
     &    (ph_alphaem0.ne.invalid)) then
        write(*,*) "calculating SM parameters from SLHA"
        mz2 = dreal(SMInputs_MZ)**2
        mw2rough = 80.385d0**2  ! used for higher order radcor only
        if (Mass_MW.ne.invalid) then
          mw2rough = (dreal(Mass_MW))**2 ! prefer SLHA
        endif
        ! QED running coupling: almzms = al0 / (1 - dalmzms)
        almzms = 1d0 / dreal(SMInputs_invAlfaMZ)  ! alpha_em^MSbar(mz)
        dalmzms = 1d0 - ph_alphaem0 / almzms     ! Delta alpha_em^MSbar(mz)
        ! convert MSbar to ONS (see PDG2012 review, Phys.Rev.D86,010001)
        ph_dalphaemmz = dalmzms - (ph_alphaem0/pi)*( ! to ONS
     &                   100/27d0 ! fermions
     &                  -1/6d0-7/4d0*log(mz2/mw2rough) ) ! W^+-
        ! weak parameters
        ph_gf = dreal(SMInputs_GF)
        almz = ph_alphaem0 / (1d0 - ph_dalphaemmz)
        s2thw = sqrt(2*sqrt(2d0)*pi*almz/(ph_gf*mz2))
        ph_sthw2 = (sin(asin(s2thw)/2d0))**2
        ! Z, W mass
        ph_Zmass = dreal(SMInputs_MZ)
        if (Mass_MZ.ne.invalid)  ph_Zmass = dreal(Mass_MZ) ! prefer SLHA
        ph_Wmass  = ph_Zmass * sqrt(1d0 - ph_sthw2)
        if (Mass_MW.ne.invalid)  ph_Wmass = dreal(Mass_MW) ! prefer SLHA
      endif

      ! extract slepton types from process id
      call decode_sfermion_pair(lprup(1),id_sla,id_slb)
      write(*,*) 'in init_couplings'
      write(*,*) 'extracted from process ID: ', lprup(1)
      write(*,*) 'the sfermion pair: ', id_sla, ', ', id_slb
      ! set slepton masses and mixing parameters in own common block
      call set_slepton_params(id_sla,ph_slepton3mass,ph_c_U3l,ph_c_U3r)
      call set_slepton_params(id_slb,ph_slepton4mass,ph_c_U4l,ph_c_U4r)

      ! running scale for alpha_em
      if (powheginput('runningemscale').ne.0) then
        ph_runningemscale = .true. ! default
      else
        ph_runningemscale = .false.
      endif
      ! override SLHA parameters by explicit parameters from powheg.input
      write(*,*) "checking for parameters in powheg.input"
      write(*,*) "(if given, they override SLHA SM parameters)"
      if (powheginput('#invalphaem0').gt.0)
     &                  ph_alphaem0 = 1d0/powheginput('invalphaem0')
      if (powheginput('#dalphaemmz').ge.0)
     &               ph_dalphaemmz = powheginput('dalphaemmz')
      if (powheginput('#Zmass').ge.0)  ph_Zmass=powheginput('Zmass')
      if (powheginput('#Wmass').ge.0)  ph_Wmass=powheginput('Wmass')
      if (powheginput('#Zwidth').ge.0) ph_Zwidth=powheginput('Zwidth')
      if (powheginput('#Wwidth').ge.0) ph_Wwidth=powheginput('Wwidth')
      if (powheginput('#sthw2').gt.0)  ph_sthw2=powheginput('sthw2')
      if (powheginput('#gf').gt.0)     ph_gf=powheginput('gf')

      ! derived quantities
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1d0 - ph_sthw2)
      ph_Zmass2 = ph_Zmass**2
      ph_Wmass2 = ph_Wmass**2

      ! set number of light flavors
      st_nlight = 5

      almz = alemeff(ph_Zmass2)

      if(verbose) then
      write(*,*) '****************************************************'
      write(*,*) 'parameters set in init_couplings:'
      write(*,*)
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'W mass = ',ph_Wmass
      write(*,*) 'Z width = ',ph_Zwidth
      write(*,*) '1/alphaem0 = ',1d0/ph_alphaem0
      write(*,*) '1/alphaem(mz^2) = ', 1d0 / almz
      write(*,*) 'sin(thetha_W)^2 = ',ph_sthw2
      write(*,*)
      write(*,*) 'U3l = ', ph_c_U3l
      write(*,*) 'U3r = ', ph_c_U3r
      write(*,*) 'U4l = ', ph_c_U4l
      write(*,*) 'U4r = ', ph_c_U4r
      write(*,*) 'slepton 3 mass = ', ph_slepton3mass
      write(*,*) 'slepton 4 mass = ', ph_slepton4mass
      write(*,*)
      if (ph_runningemscale) then
        write(*,*) 'with alpha_em(Q2)'
      else
        write(*,*) 'with fixed alpha_em'
      endif
      write(*,*) '****************************************************'
      endif

      if ((ph_alphaem0.eq.invalid).or.
     &   ((ph_dalphaemmz.eq.invalid).and.ph_runningemscale).or.
     &    (ph_Zmass.eq.invalid).or.
     &    (ph_Zwidth.eq.invalid).or.
     &    (ph_Wmass.eq.invalid).or.
     &    (ph_Wwidth.eq.invalid).or.
     &    (ph_sthw.eq.invalid).or.
     &    (ph_cthw.eq.invalid).or.
     &    (ph_gf.eq.invalid).or.
     &    (ph_Zmass2.eq.invalid).or.
     &    (ph_Wmass2.eq.invalid).or.
     &    (ph_sthw2.eq.invalid)) then
        stop "one or more parameters are not set"
      endif

      end
