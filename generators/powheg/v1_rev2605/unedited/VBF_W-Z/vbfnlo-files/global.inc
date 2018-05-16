****************************************************************************
*  This is an include-file of the VBFNLO main program.
*  It contains declarations of global variables.
c***************************************************************************
c  Author(s) : G. Klaemke
c***************************************************************************
*  Last modified: 15.09.2006
****************************************************************************

c the process identification number
      integer procID

c the number of iterations of the integration
      integer N_iterations1, N_iterations2
c the number of points/shots in the last iteration = 2**N_points
c an array is used because there are several loops (LO, virtual contrib., real emissions)
c -1       : real emissions loop
c  0       : LO loop
c  1,2,... : the various virtual contributions
      integer N_points(-1:5)	

c the center of mass energy of the collider
      double precision ecm	

c phasespace dimensionality - this is process dependend!!

      integer*4 PS_dimension

c length of the monaco random number array. (Must always be larger than PS_dimension.)

      integer max_PS_dim
      parameter (max_PS_dim = 30)

c the maximum number of 4-momenta stored in 'partons' 
c here, the partons are the incoming quarks/gluons + the 2 tagging jets
c + the real emission (for NLO only)

      integer max_p
      parameter (max_p=6)

c the maximum number of additional particles like a higgs boson, or any decay products.
c additional/decay particles means: everything which is not a parton involved in 
c the basic QCD 2->2(+1)+X process.

      integer max_v
      parameter (max_v=10)

c number of internal x parameters
      integer nx
      parameter( nx=2 )	

c the maximum number of jets 
      integer max_jets
      parameter (max_jets=max_p+max_v-nx)

c max number of momentum configurations for real emission processes
      integer max_kin
      parameter (max_kin=25)

c max number of different phase spaces
      integer max_ps
      parameter (max_ps = 5)

c actual number of entries in partons
      integer n_p

c actual number of entries in decay_particles
      integer n_v

c actual number of kinematics
      integer n_kin, n_qed

c     The particle identification code according to the Monte Carlo particle numbering scheme.
c     reminder :  e- =  11, nue    =  12, mu- =  13, numu    =  14, tau- =  15, nutau  =  16,
c                 e+ = -11, nuebar = -12, mu+ = -13, numubar = -14, tau+ = -15, nutaub = -16,
c             photon = 22,      Z0 = 23,   W+ = 24,       h0 = 25,    A0 = 36,      H+ = 37 
c     in addition, we identify a general parton (quark/gluon) with id=93
      integer particle_IDs(max_v)
c mass of the corresponding particle
      double precision particle_mass(max_v)
		
c file names for input and output of monaco grid information
      character*50 gridname1(1:max_ps), gridname2(1:max_ps), gridname3(1:max_ps)

c     proton/antiproton beams
      integer sign1, sign2

c     is NLO implemented for the process?
      logical hasNLO	

c     switches for born/virtual contribution/real emission
      logical doBorn, doVirtuals, doEmissions	


** Flags, switches, parameters for running in SUSY and with ew
** corrections
* switches: 
*  loopSqr_switch = whether |M_{loop}|^{2} is calculated
*  ewcor_switch = whether ew corns are included
*  FeynH_switch =  whether FeynHiggs is used
*  slha_switch = whether a SLHA file is used
c      logical loopSqr_switch, ewcor_switch, FeynH_switch, slha_switch
* flags: 
*  model = SM (1) and MSSM (2)
*  Higgs type  = SM (0), h0 (1), HH (2), A0 (3)
*  ewapprox = which approximations the user has selected for the ew corns
*             (0 = no approximations)
*  sector = which particles are included in the ew corns  
*  ewscheme = method by which LO em coupling is calculated       
c      integer model, HiggsType, ewapprox, sector, ewscheme
     

c     number of loops for different virtual contributions
c      integer NLO_Loops

c     number of loops for different phase space generators
c      integer PS_Loops

c put everything in global common blocks: split into integer, logical etc.
c      common /cglobalD/ ecm
c      common /cglobalI/ ewapprox, model, HiggsType, ewscheme, sector
c      common /cglobalI/ N_iterations1, N_iterations2, N_points, n_p, n_v
c      common /cglobalI/ sign1, sign2, NLO_Loops, PS_dimension, particle_IDs
c      common /cglobalI/ PS_Loops
c      common /cglobalL/ FeynH_switch, slha_switch, ewcor_switch, loopSqr_switch
c      common /cglobalL/ hasNLO, doBorn, doVirtuals, doEmissions
c      common /cglobalC/ gridname1, gridname2, gridname3 
	
      common /cglobal2/ n_kin, n_qed, procID

c distribute finite terms for Born*virtual+real emission = c*|Born|**2 + ....
c with c = -9/2 - pi^2/3 = creal + cvirt
      double precision creal, cvirt,  pi2o3p9h
      parameter ( pi2o3p9h = pi**2/3d0+9d0/2d0)
      parameter (cvirt = pi**2/3d0-7d0, creal = -pi2o3p9h-cvirt) ! for carlo

      integer  np, nv
c

      parameter( np=5 )         !number of external partons
      parameter( nv=2 )         !number of decay leptons (via decay bosons)

c fix number of external quark flavors to sum over; for CC processes both 
c up- and down-type quarks of the highest generation need to be allowed
c (i.e. the largest even number smaller or equal to nfl is used)
c nfl=5 gives the option to count initial and final b-quarks 
      integer nfl
      parameter (nfl=5)


c qsqAmin is the minimum virtuality for t-channel photon exchange. 
c Subamplitudes will be set to zero if photon virtuality falls below
c qsqAmin
      double precision qsqAmin
      parameter (qsqAmin=4d0)

c flag for warnings in matrixelement routines
      logical lwarn
      parameter (lwarn = .false.)

      integer igauge
      parameter (igauge=2)
c igauge allows to select different ways to insure EW gauge invariance
c in the face of finite W and Z width effects for EW W and Z production:
c  
c   igauge = 1       overall factor scheme
c            2       complex mass scheme, m**2-->m**2-i*m*Gamma
c                    in all propagators
c            3       for ibos=2 only: naive Z Breit Wigner
c           -1       width = 0, for testing purposes only

c      logical with_anom 	
c      common /anom_switch/ with_anom
c      logical with_anomHiggs
c      common /anomHiggs/ with_anomHiggs
c      logical with_kk 	
c      common /kk_switch/ with_kk
c      logical with_anomwgaj
c      common /anomwgaj/ with_anomwgaj

c instabilities for pentagons
      integer pentcountm, pentcount2m,pentcountps(3),pentcount2ps(3)
      common /instpent/ pentcountm, pentcount2m, pentcountps, pentcount2ps
      save /instpent/

c default values from GramDets (GGFLO)
      real*8 GRAMCTO0,GRAMDTO0
      common/GRAMTO/GRAMCTO0,GRAMDTO0

