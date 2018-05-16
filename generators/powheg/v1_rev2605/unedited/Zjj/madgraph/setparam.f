      subroutine set_madgraph_parameters
c***********************************************************************
c     !: This subroutine performs the assignation of POWHEG parameters
c     to variables in the common blocks included by MADGRAPH
c***********************************************************************
      implicit none
      include '../nlegborn.h'
      include '../../include/pwhg_flst.h'
      include '../../include/pwhg_math.h'
      include '../../include/pwhg_rad.h'
      include '../PhysPars.h'
      include 'coupl.inc'
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
ccccccccccccccccccccccccccccccccc
c     common bl. originally present in setrun.f, needed
c     by HELAS subroutines. 
C
C     BEAM POLARIZATION 
C
      double precision POL(2)
      common/to_polarization/ POL
      data POL/1d0,1d0/
C
c     common bl. originally present in driver.f, needed
c     by HELAS subroutines. 
C     isum_hel       ===> Exact helicity sum (0 yes, n = number/event)
C     multi_channel  ===> Suppress amplitude (0 no, 1 yes)
C     
      integer          isum_hel
      logical                    multi_channel
      common/to_matrix/isum_hel, multi_channel
      DATA isum_hel/0/
      DATA multi_channel/.false./
ccccccccccccccccccccccccccccccc

      double precision www

      write(*,*) 'POWHEG: set_madgraph_parameters called'

      alpha=ph_alphaem
      zmass=ph_Zmass
      zwidth=ph_Zwidth
      wmass=ph_Wmass
      wwidth=ph_Wwidth
      
      tmass= 171.3
      twidth= 1.745
      mtMS=tmass
      bmass=0d0
      mbMS=bmass
      cmass=0d0
      mcMS=cmass
      mtaMS=1.777d0

c     !ER: this is not needed for the Z case
c$$$      Vud=ph_CKM(1,1)
c$$$      Vus=ph_CKM(1,2)
c$$$      Vub=ph_CKM(1,3)
c$$$      Vcd=ph_CKM(2,1)
c$$$      Vcs=ph_CKM(2,2)
c$$$      Vcb=ph_CKM(2,3)
c$$$      Vtd=ph_CKM(3,1)
c$$$      Vts=ph_CKM(3,2)
c$$$      Vtb=ph_CKM(3,3)

c     Setting of wm MadGraph parameter. This is used only to
c     calculate the g_w (weak coupling) used in HELAS subroutines.
c     To have the same coupling of POWHEG, the following ad-hoc definition
c     of Madgraph gfermi is mandatory.
c     The following is an inversion of the assignment formula
c     for wm (see mad_setparam subroutine).
      www=zmass*sqrt(1-ph_sthw2)
      gfermi=pi*zmass**2*alpha/sqrt(2.)
      gfermi=gfermi/(zmass**2*www**2 - www**4)

c     setting of other remaining couplings is done by mad_setparam 
c     event by event 
      call mad_setparam

      write(*,*) 'Relevant couplings used by Madgraph'
      write(*,*) 'e/(sinw*sqrt(2))     = ',gwf(1)
      write(*,*) 'MCFM gwsq (from Mad) = ',(gwf(1)*sqrt(2.))**2


      end


      subroutine mad_setparam
c***********************************************************************
c     !: This subroutine is similar to setpara (Mad v4) and sets up
c     the HELAS couplings of the STANDARD MODEL without reading a card.
c     Original subroutine left in file couplings.f.orig
c***********************************************************************
      implicit none
c
c local
c
!:      character*(*) param_name
      logical readlha
      integer i
      double precision dum
c
c     common file with the couplings
c
      include 'coupl.inc'

c
c     local
c
      double precision  v
      double precision  ee, ee2, ez, ey, sw, cw, sc2, sin2w, wm
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  alphas, alfa, alfaw, mfrun
      external          alphas, alfa, alfaw, mfrun
c
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Assign to MADGRAPH routines the 
C  POWHEG a_s value event by event
C
Cccccccccccccccccccccccccccccccccccccccccccc
      include '../../include/pwhg_st.h'
      alfas=st_alpha
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c------------------------------------------
c Start calculating the couplings for HELAS
c------------------------------------------
c
      G = SQRT(4d0*PI*ALFAS)   ! use setting of the param_card.dat @ NLO

c     
c     Strong coupling
c
c     As a rule we first check if a pdf has been chosen in the    
c     run_card.dat (which has been already read at this stage).
c     If there pdfs in the initial state, then the alpha_s(MZ) used
c     is set to the corresponding value.  

      GG(1) = -G
      GG(2) = -G     
c
c auxiliary local values
c
      wm = sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      sin2w  = One-(wm/zmass)**2
      cw  = sqrt( One - sin2w )
      ee2 = alpha * Fourpi
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee   ! the wmass is used to calculate v
      lambda = hmass**2 / (Two * v**2)
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c gauge & higgs boson coupling constants
c
      gwwh  = cmplx( ee2/sin2w*Half*v, Zero )
      gzzh  = cmplx( ee2/sc2*Half*v, Zero )
      ghhh  = cmplx( -hmass**2/v*Three, Zero )
      gwwhh = cmplx( ee2/sin2w*Half, Zero )
      gzzhh = cmplx( ee2/sc2*Half, Zero)
      ghhhh = ghhh/v
c
c fermion-fermion-vector couplings
c
      gal(1) = cmplx(  ee          , Zero )
      gal(2) = cmplx(  ee          , Zero )
      gau(1) = cmplx( -ee*Two/Three, Zero )
      gau(2) = cmplx( -ee*Two/Three, Zero )
      gad(1) = cmplx(  ee/Three    , Zero )
      gad(2) = cmplx(  ee/Three    , Zero )

      gwf(1) = cmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = cmplx(  Zero              , Zero )

      

      gzn(1) = cmplx( -ez*Half                     , Zero )
      gzn(2) = cmplx(  Zero                        , Zero )
      gzl(1) = cmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = cmplx( -ey                          , Zero )
      gzu(1) = cmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = cmplx(  ey*Two/Three                , Zero )
      gzd(1) = cmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = cmplx( -ey/Three                    , Zero )
c
c fermion-fermion-Higgs couplings (complex) hff(2)
c
c NOTE: the running mass is evaluated @ the same order 
c nloop of alpha_s set by the PDF choice
c 

      if(mtMS.gt.1d0) then
         ghtop(1) = cmplx( -mtMS/v, Zero )
      else
         ghtop(1) = cmplx( Zero,Zero)
      endif
      ghtop(2) = ghtop(1)

      if(mbMS.gt.1d0) then
         ghbot(1) = cmplx( -mbMS/v, Zero )
      else
         ghbot(1) = cmplx( Zero, Zero )
      endif
      ghbot(2) = ghbot(1)
      
      if(mcMS.gt.1d0) then
         ghcha(1) = cmplx( -mcMS/v, Zero )
      else
         ghcha(1) = cmplx( Zero, Zero )
      endif
      ghcha(2) = ghcha(1)

      ghtau(1) = cmplx( -mtaMS/v, Zero )
      ghtau(2) = ghtau(1)
c
c     CKM matrix: 
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd Vcb=Vub=0
c
c     >>>>>>>>>>>>>>>***** NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter only when interaction_CKM.dat
c     is used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal
c     CKM matrix.


c     CKM matrix: 
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd Vcb=Vub=0
c
c     >>>>>>>>>>>>>>>***** NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter only when interaction_CKM.dat
c     is used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal
c     CKM matrix.

c$$$      do i=1,2  !ER: cancellare....
c$$$         gwfud(i) = gwf(i)*Vud
c$$$         gwfus(i) = gwf(i)*Vus
c$$$         gwfub(i) = gwf(i)*Vub
c$$$         gwfcd(i) = gwf(i)*Vcd
c$$$         gwfcs(i) = gwf(i)*Vcs
c$$$         gwfcb(i) = gwf(i)*Vcb
c$$$         gwftd(i) = gwf(i)*Vtd
c$$$         gwfts(i) = gwf(i)*Vts
c$$$         gwftb(i) = gwf(i)*Vtb
c$$$      enddo

c---------------------------------------------------------
c Set Photon Width to Zero, used by symmetry optimization
c---------------------------------------------------------

      awidth = 0d0
            
c----------------------------
c end subroutine coupsm
c----------------------------

      return
      end
