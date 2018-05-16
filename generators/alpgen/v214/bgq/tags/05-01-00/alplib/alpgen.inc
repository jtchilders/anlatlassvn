c
c hard process specification
      integer ihrd
      common/prcnum/ihrd
      integer itopprc
      common/singletop/itopprc
c numerical parameters
      real *8 pbfac,pi
      parameter (pbfac=3.8937d8,pi=3.14159265358979324d0)
c particle parameters
      real *8 amass(24),mlep(3),mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
      common/pparam/amass,mlep,mc,mb,mt,mw,wwid,mz,zwid,mh,hwid
c couplings-related EW parameters
      integer nw,nz,nh,nph,iewopt
      real *8 aem,gfermi,gw,g1,gbar,stw,ctw,ccab2,scab2,ggh
      common/alpcpl/aem,gfermi,gw,g1,gbar,stw,ctw,ccab2,scab2,ggh,iewopt
      common/ngabos/nw,nz,nh,nph
      double precision cosvma
      integer itdec
      common/topcoup/cosvma,itdec
c alpha internal parameters
      double precision apar
      common/alphapar/apar(100)
c parameters for imode=2
      integer iWdecmode,iZdecmode,itdecmode
      common/decpar/iWdecmode,iZdecmode,itdecmode
c weight information
      real *8 maxwgt,avgwgt,totwgt,maxev
      integer unwev
      common/weights/maxwgt,avgwgt,totwgt,maxev,unwev
      logical nowgt,dump,evgen,colgen,evopt
      common/logcls/nowgt,dump,evgen,colgen,evopt
      integer niter
      double precision nopt
      common/optiter/nopt,niter
c colour/spin averaging common
      integer navg
      common/cspavg/navg
      integer nprtns
      common/npartons/nprtns
c
c     generation mode:
c     0: generate weighted events
c     1: generate weighted events, write events to file for later unweighting
c     2: read events from file for unweighting or processing
      integer imode
c     fname.stat: run statistics 
c     fname.top : topdrawer plots
c     fname.par : file with run parameters
c     fname.wgt : file with weighted events
c     fname.unw : file with unweighted events
c     fname_mpirank has the MPI rank string added to it
      character*140 fname
      character*512 fname_mpirank
      common/evtout/imode,fname,fname_mpirank
c topdrawer output file
      character *50 topfile
      common/topout/topfile
c grids I/O
      integer igrid
      common/gridio/igrid
c
c event kinematics: p(i,j)
c i=1-3: space components of momentum, i=4: energy, i=5: mass
c j=1,2: incoming partons, with x1,x2 beam momentum fractions
c j=3,4,...,npart: final state particles
      real *8 p
      integer ifl,icu
      integer maxpar, npart
      parameter (maxpar=10)
      common/event/p(5,maxpar),npart,ifl(maxpar),icu(2,maxpar)
c                                     
c beams, parton densities, alphas
      real *8 x1,x2,qsq,xlam,as,qfac,ebeam,roots,s,resc
      integer nflpdf
      real f1(-5:5),f2(-5:5)
      integer iqopt,ndns,nmnr,nloop,ih1,ih2       
      character*2 sche
      character*25 pdftyp
      common/beams/ebeam,roots,s,ih1,ih2
c iqopt: specifies the functional form of Q2; see code for mapping
c E.g. WBBGEN: iqopt=0: Q0^2=MW^2;   iqopt=1: Q0^2=MW^2+sum(Etjets)^2
c qfac:  Q=qfac*Q . Typically qfac=0.5, 1, 2
      common/pdfvar/x1,x2,qsq,qfac,xlam,as,resc,iqopt,ndns,nmnr,
     +   nloop,sche,pdftyp
      common/pdflum/f1,f2
c resonance and gauge invariance contorl parameters
      character*2 wmode
      character resonance
      real*8 winsize
      common/gauinv/winsize,resonance,wmode
c clustering parameters
      integer mxcpar
      parameter (mxcpar=20)
c     kres= resolution parameter (gev^2) -- not used
c     ickkw= if=1 => ckkw alphas rescaling
c     ifcol= appears in flavour-driven clustering, not used
c     avgrew= average alphas reweighting factor
c     xlclu= lambdaQCD used for reweighting
c     lpclu= loop order used in alphas for sudakovs and reweighting
c     clp=cluster momentum
c     icst= status: 0=unclustered, 1=clustered
c     icfl= cluster flavour
c     icfs= ini/final state (-1/+1)
c     icdau= cluster daughters
c     icmot= mother cluster
c     clkt= 1: kt at vtx where cluster/particle originates 
c           2: kt at vtx where cluster branches (kres for external states)
c     nclext= number of external particles to include in clustering
c     ncltot= total number of external particles plus clusters
c     
      double precision clp,clkt,kres,xlclu,avgrew
      integer ifcol,mrgopt,cluopt,icst,icfl,icfs,icdau,icmot,iccol,
     $  nclext,ncltot,lpclu,ickkw
      common/alclus/kres,avgrew,xlclu,clp(4,mxcpar),clkt(2,mxcpar),
     $     icst(mxcpar),icfl(mxcpar),icfs(mxcpar),icdau(2,mxcpar),
     $     icmot(mxcpar),iccol(2,mxcpar),nclext,ncltot,ifcol,ickkw,
     $     mrgopt,cluopt,lpclu
      double precision asmax,ktfac
      integer naspow
      common/asckkw/asmax,ktfac,naspow
c
c     I/O unit numbers
      integer niosta,niopar,niowgt,niounw
      common/iounits/niosta,niopar,niowgt,niounw
      integer iseed,iseed2
      common/seeds/iseed(2),iseed2(2)
c
c overall parameters
      integer npar,nprocs,nactprc
      parameter (npar=200,nprocs=20)
      integer parlen,partyp,paruse
      character chpar*8,chpdes*70,chprc*10
      double precision parval
      common/alpars/parval(npar),chpar(npar),chpdes(npar),parlen(npar),
     +  paruse(npar,nprocs),partyp(npar)
      common/alpprc/chprc(nprocs),nactprc
c hidden parameters
      integer impunw
      common/alhidp/impunw
c global event cuts
      real *8       ptjmin,ptjmax,etajmax,drjmin,
     +        ptj1min,ptj1max,
     +        ptQmin,ptQmax,etaQmax,drQmin,
     +        ptQ2min,ptQ2max,etaQ2max,drQ2min,
     +        ptbmin,ptbmax,etabmax,drbmin,
     +        ptcmin,ptcmax,etacmax,drcmin,
     +        ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +        ptlmin,etalmax,drlmin,metmin,mllmin,mllmax,
     +        ptcen,etagap,pthrmin,pthrmax,ptphmax

      integer irapgap
      common/glocut/ptjmin,ptjmax,etajmax,drjmin,
     +        ptj1min,ptj1max,
     +        ptbmin,ptbmax,etabmax,drbmin,
     +        ptcmin,ptcmax,etacmax,drcmin,
     +        ptQmin,ptQmax,etaQmax,drQmin,
     +        ptQ2min,ptQ2max,etaQ2max,drQ2min,
     +        ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +        ptlmin,etalmax,drlmin,metmin,mllmin,mllmax,
     +        ptcen,etagap,pthrmin,pthrmax,ptphmax,irapgap


c
c
c hard process specification
      double precision mq,mq2
      integer njets,ihvy,ihvy2,ilep,zin,zfstate
      common/hrdprc/mq,mq2,njets,ihvy,ihvy2,ilep,zin,zfstate
