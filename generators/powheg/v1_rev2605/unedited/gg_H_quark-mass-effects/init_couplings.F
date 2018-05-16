      subroutine init_couplings
      implicit none
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

      model = int(powheginput('model'))
      if (model.ne.0) then
         write(*,*) 'Other model than the SM not supported in 
     $this version'
         stop
      endif

      massren = int(powheginput('massren'))
      if (massren.eq.0) then
         write(*,*) 'On-shell mass renormalization scheme selected.'
      else if (massren.eq.1) then
         write(*,*) 'MSBar mass renormalization scheme selected'
      else if (massren.eq.2) then
         write(*,*) 'DRBar mass renormalization scheme selected'
      else
         write(*,*) 'Unknow mass renormalization scheme selected.
     $Aborting.'
         stop
      endif
      
      write(*,*) 'SM selected'
      call init_ind()
      call init_fermions()
      call init_higgs()
      call init_ew()

      write(*,*) 'Active fermions'
      write(*,*) afer
      write(*,*) 'fermion masses'
      write(*,*) mfer
      write(*,*) 'fermion couplings'
      write(*,*) lambdafer
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) 'GF = ',ph_GF

      end subroutine

c     EW corrections initialization
      subroutine init_ew()
      implicit none
      include 'Flags.h'
      real * 8 powheginput
      external powheginput

      flg_ew = int(powheginput('ew'))
      if (flg_ew.le.0) then
         flg_ew = 0
         write(*,*) '2-loops EW corrections disabled'
         flg_fast_ew = 0
         flg_passarino = 0
      else
         write(*,*) '2-loops EW corrections enabled'
c     For release only passarino is supported, with mass sampling
         flg_ew = 1
         flg_fast_ew = 1
         flg_passarino = 1
         call init_cached_ew_corr
      endif

      end subroutine

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

      if ((flg_hdecay.gt.0).and.(model.eq.0)) then
         call hdecayparser(ph_Hmass,ph_Hwidth)
      else
         ph_Hwidth = powheginput('hwidth')
      endif

      ph_Zmass  = 91.1876d0
      ph_Zwidth =  2.4952d0
      ph_Wmass  = 80.398d0
      ph_Wwidth =  2.141d0

      ph_alphaem = 1d0/137.035999679d0
      ph_alphaemmz = 1d0/128d0
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


c     If fast_ew is enabled, create the array with the campionazied deltaew value.
      subroutine init_cached_ew_corr
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      integer i
      real * 8 m12,m0,y12,y0,value,ph_Hmasshigh
      complex * 16 ampl, reduced,x12,x0,aux,bornqcd
      common /bornampl/ampl,bornqcd
      real * 8 deltaew
      external reduced,deltaew

      ph_Hmasshigh = sqrt(ph_Hmass2high)

      write(*,*) 'Starting cached EW corrections array initialization'
      if ((ph_Hmass.le.700d0).and.(ph_Hmasshigh.gt.1800d0)) then
              write(*,*) 'The POWHEG-BOX aborts: 
     $For m_H less or equal than 700 GeV, EW corrections available
     $for virtualities up to 1800GeV. Please fix the maxwindow parameter 
     $accordingly.'
              call exit(1)
      else if (((ph_Hmass.gt.700d0).and.(ph_Hmass.le.1000d0)).and.
     $          (ph_Hmasshigh.gt.2000d0)) then
              write(*,*) 'The POWHEG-BOX aborts: 
     $For m_H greater than 700 and less than 1000, EW corrections
     $available for virtualities up to 2000GeV.
     $Please fix the maxwindow parameter accordingly.'
         call exit(1)
      else if (ph_Hmass.gt.1000d0) then
              write(*,*) 'The POWHEG-BOX aborts: 
     $For m_H greater than 1000, EW corrections are NOT available. 
     $Please disable them.' 
         call exit(1)
      end if

c     2500 is the end of Passarino's grid
      do mh=0,2500
         value = 0d0
         if (flg_passarino.eq.1) then
            call exa_ewgrid(mh,value)
c            write(*,*) 'mh',mh, value
            cached_ew_corr(nint(mh)) = value/100d0
         else
c     write(*,*) 'mh', mh
            ampl = dcmplx(0d0)
            do i=1,afer
               m12=mfer(i)
               y12=m12**2/mh**2
               x12 = reduced(1d0/y12)
               aux = lambdafer(i)*trfer(i)*
     &              (-4d0)*y12*(2d0-(1d0-4d0*y12)*0.5d0*log(x12)**2)
               ampl = ampl+aux
            end do
            cached_ew_corr(int(mh)) = deltaew(ampl)
         end if
         write(33,*) 'deltaew',mh, cached_ew_corr(nint(mh))
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

c     Fermions initialization subroutines used by SM and MW
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

      afer = 0
      ph_topmass = powheginput('#topmass')
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
      ph_charmmass = powheginput('#charmmass')
      if (ph_charmmass.ne.-1000000d0) then
         write(*,*) "Charm Quark enabled"
         afer = afer+1
         trfer(afer) = 1d0/2d0
         mfer(afer) = getRSmass(ph_charmmass,massren,ph_Hmass,alpha_s)
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
