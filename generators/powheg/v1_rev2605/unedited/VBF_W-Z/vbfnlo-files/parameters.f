********************************************************************************
********************************************************************************
**                                                                           ***
******************************************************************           ***
*** sophy@particle.uni-karlsruhe.de                            ***           ***
***                                                            ***           ***
*** This set of subroutines sets all of the parameters used    ***           ***
*** in the calculation of electroweak corrections in the code: ***           ***
*** i.e. it calls FeynHiggs, and calculates Higgs masses and   ***           ***
*** renormalisation constants.                                 ***           ***
***                                                            ***           ***
*** VERY IMPORTANT NOTE!!!                                     ***           ***
*** The implementation of the user input EWSCHEME was changed  ***           ***
*** for version 2.5.0 of VBFNLO.  In order to switch back to   ***           ***
*** the old implementation, comment out the NEW subroutine     ***           ***
*** 'setEWpara' and uncomment the OLD subroutine 'setEWpara',  ***           ***
*** which are the first and second subroutines in this file    ***           ***
*** respectively.                                              ***           ***
*** For more details, see the wiki.                            ***           ***
******************************************************************           ***
**                                                                           ***
********************************************************************************
********************************************************************************

*** NEW v 2.5.0 subroutine that sets electroweak parameters according to 
*** EWSCHEME, and e and g2, which control the couplings.

      subroutine setEWpara(e,g2,s,c,z,w,q,g)

      implicit none

      double precision e, g2, s, c, z, w, q, g, xmwsq


      include "koppln_ew.inc"
! #include "global.inc"
      include "mssm.inc"

      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)
** ewscheme = 1: MW2 and SW2 are set according to the inputs alfa, gfermi and
** MZ.  NOTE ALSO:  using this scheme can lead to odd values of the W mass.  
** A warning is output if this is the case.
      if (ewscheme .eq. 1) then
         XMWSQ = XMZ**2/2.d0 + SQRT(XMZ**4/4.d0 - PI*ALFA*XMZ**2/
     -        SQRT(2.d0)/GF) 
         SIN2W = 1.d0 - XMWSQ/XMZ**2
         XMW = XMZ*SQRT(1.d0 - SIN2W)

** ewscheme = 2: alfa_qed and W mass are set according to the inputs gfermi,
** MZ and SW2
      else if (ewscheme .eq. 2) then
         XMW = XMZ*SQRT(1.d0 - SIN2W)

** ewschemes 3, 5 and 6: SW2 is set according to the input MW and MZ
      else if ((ewscheme .eq. 3) .or. (ewscheme .ge. 5)) then
        SIN2W = 1.d0 - (XMW/XMZ)**2
      end if


      G = SQRT(4*PI*ALFAS)
      S = SQRT(SIN2W)
      C = SQRT(1.d0 - SIN2W)
   

** G2 = E/S IS CALCULATED FROM GFERMI AND THE Z MASS
      G2 = SQRT(8.d0*GF/SQRT(2.d0))*XMZ*C
 
** For ewscheme = 1,5,6, all couplings are set according to the input alfa.
      if ((ewscheme .eq. 1) .or. (ewscheme .ge. 5)) then
         E = SQRT(4.d0*PI*ALFA)
         G2 = E/S

** For ewscheme 4, photon couplings are set according to the input alfa, and 
** all other couplings are set according to the input gfermi.  NOTE: This can
** lead to problems with gauge invariance if alfa and gfermi are not consistent!
      else if (ewscheme .eq. 4) then
         E = SQRT(4.d0*PI*ALFA)

** For ewscheme 2 and 3, all couplings are set according to the input gfermi
      else if ((ewscheme .eq. 2) .or. (ewscheme .eq. 3)) then
         E = G2*S
         ALFA = E**2/(4.d0*PI)
      end if


      Z = G2/4.d0/C
      W = G2/SQRT(8.d0)
      Q = G2*C


* call to subroutine setting parameters 
*      if ((HiggsType .ne. 0) .or. ewcor_switch) then
      call setparams(G2)
*      end if


** WARNING MESSAGES:
      if (ewscheme .eq. 4) then
         write(*,*)'  '
         write(*,*)'WARNING! Note that the inputs:'
         write(*,*)'  ALFA, FERMI_CONST, WMASS, ZMASS and SIN2W '
         write(*,*)'are not independent quantities.  If they are not'
         write(*,*)'consistent, problems with gauge invariance may'
         write(*,*)'arise.'
         write(*,*)'In this scheme, the photon couplings are set'
         write(*,*)'according to the input ALFA, and all other'
         write(*,*)'couplings are set according to FERMI_CONST.'
         write(*,*)'  '
      end if

      if (XMW .lt. 80D0) then
         write(*,*)'  '
         write(*,*)'WARNING!! MW =', XMW
         write(*,*)'This is very low!'
         write(*,*)'  '
      end if
    

      end



********************************************************************************
********************************************************************************

c$$$*** OLD (pre v 2.5.0) subroutine that sets electroweak parameters according 
c$$$*** to EWSCHEME, and e and g2, which control the couplings.
c$$$
c$$$      subroutine setEWpara(e,g2,s,c,z,w,q,g)
c$$$
c$$$      implicit none
c$$$
c$$$      double precision e, g2, s, c, z, w, q, g, xmwsq
c$$$
c$$$** Flag that will warn if there are problems with consistency - i.e. fermion
c$$$** and higgs couplings set according to different inputs
c$$$      integer iflag
c$$$
c$$$#include "koppln.inc"
c$$$#include "global.inc"
c$$$#include "mssm.inc"
c$$$
c$$$
c$$$      iflag = 0
c$$$
c$$$
c$$$* for 1/alpha_QED < 136 the Weak angle is calculated from the input
c$$$* values of alpha(assumed at scale m_W**2), m_Z and G_Fermi
c$$$
c$$$** setting MW and sin_theta_w according to the input variable ewscheme
c$$$      if (ewscheme .eq. 1) then
c$$$         if (alfa .gt. 1d0/136d0) then
c$$$            XMWSQ = XMZ**2/2.d0 + SQRT(XMZ**4/4.d0 - 
c$$$     -           PI*ALFA*XMZ**2/SQRT(2.d0)/GF) 
c$$$            SIN2W = 1.d0 - XMWSQ/XMZ**2
c$$$            XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$         else
c$$$** NOTE: if alfa_qed is input as less than or equal to alfa(0),sin^2 theta_w
c$$$**       is fixed at 0.2312.
c$$$            SIN2W = 0.2312d0
c$$$            XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$         end if
c$$$      else if (ewscheme .eq. 2) then
c$$$         XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$      else if ((ewscheme .eq. 3) .or. (ewscheme .ge. 5)) then
c$$$        SIN2W = 1.d0 - (XMW/XMZ)**2
c$$$      end if
c$$$
c$$$
c$$$      G = SQRT(4*PI*ALFAS)
c$$$      S = SQRT(SIN2W)
c$$$      C = SQRT(1.d0 - SIN2W)
c$$$   
c$$$c
c$$$c G2 = E/S IS CALCULATED FROM GFERMI AND THE Z MASS
c$$$c
c$$$      G2 = SQRT(8.d0*GF/SQRT(2.d0))*XMZ*C
c$$$      IF ( ALFA.GT.1.d0/136.d0 .or. ewscheme .ge. 5) THEN
c$$$** ewscheme 1 or 4 if alfa > alfa(0), and ewscheme 5 and 6
c$$$        E = SQRT(4.d0*PI*ALFA)
c$$$        G2 = E/S
c$$$      ELSEIF ( ALFA.GT.1.d0/138.d0 ) THEN
c$$$** ewscheme 1 or 4 if alfa = alfa(0): e calculated from (input) alfa, G2 from 
c$$$** g_fermi. WARNING! gauge cancellations may no longer work!
c$$$        E = SQRT(4.d0*PI*ALFA)
c$$$        iflag = 1
c$$$      ELSE
c$$$** ewscheme 2 or 3 or if alfa < alfa(0): e, G2 and alfa all calculated from 
c$$$** g_fermi. 
c$$$        E = G2*S
c$$$        ALFA = E**2/(4.d0*PI)
c$$$        if (ewscheme .eq. 1 .or. ewscheme .eq. 4) then
c$$$           iflag = 2
c$$$        end if
c$$$      ENDIF
c$$$
c$$$
c$$$      Z = G2/4.d0/C
c$$$      W = G2/SQRT(8.d0)
c$$$      Q = G2*C
c$$$
c$$$
c$$$* call to subroutine setting parameters 
c$$$*      if ((HiggsType .ne. 0) .or. ewcor_switch) then
c$$$      call setparams(G2)
c$$$*      end if
c$$$
c$$$
c$$$      if (iflag .eq. 1) then
c$$$         write(*,*)"WARNING! The Higgs coupling is set from G_FERMI,"
c$$$         write(*,*)"and the electron charge is set from the input ALFA."
c$$$         write(*,*)"Gauge cancellations may no longer work"
c$$$      else if (iflag .eq. 2) then
c$$$         write(*,*)"NOTE that all couplings are set from G_FERMI, not"
c$$$         write(*,*)"from the input ALFA"
c$$$      end if
c$$$
c$$$
c$$$      end



********************************************************************************
********************************************************************************

***  Subroutine that sets all (non-ewscheme) parameters 

      subroutine setparams(G2)

      implicit none


*** The following three parameters are used by FeynHiggs.  They are set in the
*** subroutine 'bench_set':

** This is a FeynHiggs flag: t1_CplxApprox determines how the two-loop 
** corrections are treated in the presence of complex parameters
      integer t1_CplxApprox

** The soft-SUSY breaking parameters in the sfermion sector
      double precision MSusy

** Scale information for FeynHiggs
      double precision scalefactor


** These do-loop counters are used to set squark and other SUSY parameters 
      integer sfe1, sfe2, gen, type
      integer neu1, neu2, cha1, cha2

** Flag identifying which of the neutral/charged currents is being considered
      integer cur

** Variables used for choosing the charge renormalisation scheme
      integer charge
      double precision Alfa0, AlfaMZ, AlfaGF
      double precision delA0, delAmz, delAgf
      double precision G2

** Flags for the PDF set which is used
      integer pdf

** Flag for whether to include bottom Yukawa coupling, and value of Delta MB
** from FeynHiggs
      logical delMB_switch
      double complex DeltaMB


** By (annoying) convention, the sign of the su(2) covariant derivative is 
** different in the SM and MSSM. This factor takes care of this fact, and 
** ensures everything's consistent.
      double precision su2sign


** Dummy parameters used in calculation of the tree level higgs sector
      double precision avgMH2, deltaMH2
      double precision Mh02tree, MHH2tree, MA02tree, MHp2tree


** Dummy parameters: these are used to calculate the SUSY parameters (sfermion,
** chargino, neutralino sectors) if they haven't been set by either a SLHA file
** or FeynHiggs
      logical first
      double precision MSfcalc(2,4,3), MChaCalc(2), MNeuCalc(4)
      double complex USfcalc(2,2,4,3), UChaCalc(2,2), VChaCalc(2,2)
      double complex ZNeuCalc(4,4)


** Function which determines whether a parameter has already been set with a 
** reasonable value (for real and complex parameters)
      logical replace
      external replace

** Debugging flag
      integer debug
      
      
** Parameter declarations
      include "koppln_ew.inc"
! #include "global.inc"
      include "mssm.inc"
      double precision Pi
      parameter (Pi = 3.1415926535897932384626433832795029D0)
      ! set debug
      debug = 0


** Setting gauge boson parameters to match input (read either from SLHA file - 
** in slha_read.F, or vbfnlo.dat - in koppln.F)
      MZ = DBLE(xmz)
      MW = DBLE(xmw)

      MZ2 = MZ**2D0
      MW2 = MW**2D0

      SW2 = Sin2W
      SW = sqrt(SW2)
      CW2 = 1 - SW2
      CW = sqrt(CW2)

*      CW2 = MW2/MZ2
*      SW2 = 1 - CW2
*      CW = sqrt(CW2)
*      SW = sqrt(SW2)



   
** Light fermion masses - set according to the PDG if they haven't been set in 
** the slha file.  c,, t and b masses will have been set by this point, either 
** by a slha file or via the .dat input files
      if (replace(Mf(2,1),0)) then
         ME = 0.00051099891D0
      else
         ME = Mf(2,1)
      end if
      if (replace(Mf(3,1),0)) then
         MsU = 0.0024D0
      else
         MsU = Mf(3,1)
      end if
      if (replace(Mf(4,1),0)) then
         MD = 0.00475D0
      else
         MD = Mf(4,1)
      end if
      if (replace(Mf(2,2),0)) then
         MM = 0.105658367D0
      else
         MM = Mf(2,2)
      end if
      if (replace(Mf(4,2),0)) then
         MS = 0.104D0
      else
         MS = Mf(4,2)
      end if

** b,t,c and tau masses - input (either slha or vbfnlo.dat)
      MC = DBLE(xmc)
      MB = DBLE(xmb)
      MT = DBLE(xmt)
      ML = DBLE(xmtau)

** Quark masses - PDG
*      MC = 1.27D0 
*      MT = 172.6D0   
*      MB = 4.2D0  
** Matching 0804.2676
*      MT = 170.9D0
*      MB = 4.7D0
** Matching 0710.4749
*      MsU = 0.0066D0
*      MC = 1.2D0
*      MT = 174.3D0      
*      MD = 0.0066D0
*      MS = 0.150D0
*      MB = 4.3D0
*      ME = 0.000510998927D0
*      MM = 0.105658369D0
*      ML = 1.77699D0

      MU2 = MsU**2D0
      MC2 = MC**2D0
      MT2 = MT**2D0
      MD2 = MD**2D0
      MS2 = MS**2D0
      MB2 = MB**2D0


** Choosing charge renormalisation scheme based on value of AlfaQED and ewscheme
*      AlfaQED = (G2**2D0)*SW2/(4D0*Pi)
      AlfaQED = Alfa
      AlfaQED2 = AlfaQED**2D0
      gfermi = GF


      EL = sqrt(4D0*Pi*AlfaQED)
*      write(*,*)'Alfa =', AlfaQED


** These are the values that the charges and L/R axial vector couplings can
** take for up-type and down-type fermions.  The variable 'su2sign' takes care
** of the always-fun sin-theta-w conventions
!       if (model .eq. 2) then
!          su2sign = -1D0
!       else
         su2sign = 1D0
!       end if

      Qu = 2/3D0
      Qd = -1/3D0
      
      gL3 = su2sign*((1/2D0)-(2*SW2/3D0))/(SW*CW)
      gR3 = -su2sign*(2*SW/(3D0*CW))
      gL4 = su2sign*(-(1/2D0)+(SW2/3D0))/(SW*CW)
      gR4 = su2sign*SW/(3D0*CW)


** Setting renormalisation scales
      scalefactor = 1D0
      massQ = 1D-5





         MH = XMH
         MH2 = MH**2D0


    
** Deubgging messages
      if (debug .ne. 0) then
         write(*,*)'    '
         write(*,*)'SM parameters ...'
         write(*,"(A,F7.3)")'MS = ', MS
         write(*,"(A,F7.3)")'MC = ', MC
         write(*,"(A,F7.3)")'MB = ', MB
         write(*,"(A,F7.3)")'MT = ', MT
         write(*,"(A,F7.3)")'MW = ', MW
         write(*,"(A,F7.3)")'MZ = ', MZ
      
!          if (model .eq. 2) then
!             write(*,*)'      '
!             write(*,*)'SUSY parameters ...'
!             write(*,"(A,F6.3)") "Tan[beta] = ", TB
!             write(*,"(A,F8.3)") "MA0 = ", MA0
!             write(*,"(A,F8.3)") "MHp =", MHp
!             write(*,"(A,'('F8.3,',',F8.3')')") "MUE = ", MUE
!             write(*,"(A,'('F8.3,',',F8.3')')") "At = ", Af(3,3,3)
!             write(*,"(A,F8.3)") "MSusy = ", MSusy
!             write(*,"(A,'('F8.3,',',F8.3')')") "M_1 = ", M_1
!             write(*,"(A,'('F8.3,',',F8.3')')") "M_2 = ", M_2
!             write(*,"(A,'('F8.3,',',F8.3')')") "M_3 = ", M_3
!          end if
      end if


** Quark masses:  Here I reset all first and second generation quark masses to
** a single parameter massQ, which controls the collinear divergences.  massQ
** is set in the file renormalisation.F.
      MsU = massQ
      MC = massQ
      MD = massQ
      MS = massQ
 
      MU2 = MsU**2D0
      MC2 = MC**2D0

      MD2 = MD**2D0
      MS2 = MS**2D0

      ME2 = ME**2D0
      MM2 = MM**2D0
      ML2 = ML**2D0


** Setting the fermion masses
      Mf(1,1) = 0D0
      Mf(1,2) = 0D0
      Mf(1,3) = 0D0

      Mf(2,1) = ME
      Mf(2,2) = MM
      Mf(2,3) = ML

      Mf(3,1) = MsU
      Mf(3,2) = MC
      Mf(3,3) = MT

      Mf(4,1) = MD
      Mf(4,2) = MS
      Mf(4,3) = MB

      do type = 1, 4
         do gen = 1, 3
            Mf2(type, gen) = Mf(type, gen)**2D0
         end do
      end do



      end 



********************************************************************************
********************************************************************************

*** A function which determines whether a parameter needs to be set from the 
*** .dat input files, or via an internal calculation.  

      logical function replace(para,calc)

      implicit none

** The parameter which (perhaps) needs to be set
      double precision para

** What type of parameter it is:
**   0: an input parameter (4-3-11 for full list)
**   1: sfermion masses and mixings
      integer calc 

      
      include "koppln_ew.inc"
! #include "global.inc"
      include "mssm.inc"

      
      replace = .false.


      if (calc .eq. 0) then

            replace = .true.

      elseif (calc .eq. 3) then ! Higgs widths


         if (para .eq. -999d0) replace = .true.

      else

         write(*,*)'SOPHY! You have a problem with REPLACE.'
         stop

      end if         


      end



********************************************************************************
********************************************************************************

*** Subroutine to initialise branching ratios and total widths to -999

      subroutine clearwidths

      implicit none


      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT


** widths:
      xgw = -999d0
      xgz = -999d0
      xgh = -999d0
      gamt = -999d0


** w branching ratios
      bwne = -999d0
      bwud = -999d0
      bwtb = -999d0


** z branching ratios
      bznn = -999d0
      bzee = -999d0
      bzuu = -999d0
      bzdd = -999d0
      bztt = -999d0


** higgs branching ratios
      bhww = -999d0
      bhzz = -999d0
      bhgg = -999d0
      bhtt = -999d0
      bhbb = -999d0
      bhcc = -999d0
      bhtau = -999d0
      bhmu = -999d0
      bhgam = -999d0
      bhgamz = -999d0


      end


********************************************************************************
********************************************************************************

*** This subroutine checks that the relevant branching ratio for a particular
*** process is not greater than one.  Only resonant Higgs bosons are considered.

      subroutine checkBR

      implicit none
! 


      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT



         if (BHGAM .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> 2 photons'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHGAM
            stop
         end if

         if (BHMU .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> 2 muons'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHMU
            stop
         end if

         if (BHTAU .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> tau+ tau-'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHTAU
            stop
         end if

         if (BHBB .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> b bbar'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHBB
            stop
         end if

         if (BHWW .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> W+ W-'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHWW
            stop
         end if

         if (BHZZ .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> Z Z'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHZZ
            stop
         end if



      end 
