c     ttjdecay knows if the current event is a born or a real by means
c     of nup value in the LesHouches common block. To perform the decay,
c     only massless momenta must be used (no reshuffling) because an
c     evaluation of Madgraph decayed matrix element is required and that
c     has to be done with the same momenta used to calculate POWHEG
c     matrix elements, where all particles (apart the top) are massless.
c     The reshuffling of momenta, in order to give masses to top (and W)
c     decay products, is done after the decay, by the routine
c     reshufflemom
      
c     This is different wrt the POWHEG-hvq code, and its POWHEGBOX
c     implementation under hvq/, where a nonvanishing bmass value was
c     used.  Here I prefer to avoid this, because bmass enters also in
c     the 'undecayed part' of the Madgraph full matrix element, being b
c     a light parton that can be produced together with the ttbar
c     pair. Moreover, at variance with single-top processes under ST_*/,
c     here it is not possible to spread the top and antitop virtuality
c     along a Breit-Wigner before calculating the matrix elements,
c     because they assume the two masses of the top and anti-top to be
c     equal. So the (anti-)top mass is kept fixed until the call to the
c     twvirts routine below. Inside there, new values for tops and W's
c     virtuality are generated, taking into account phase space and
c     luminosity correction factors too. Alternatively, it is always
c     possible to enforce the zero width limit by adding zerowidth 1 to
c     the POWHEG input card.  In any case, it is mandatory to always set
c     the top and W widths to a non zero value in MADGRAPH routines for
c     decayed matrix elements, otherwise the top and W propagators
c     present there may cause divergences (MADGRAPH propagators do
c     include the width).  There is still an ambiguity on which width
c     enters the undecayed matrix elements, two competing possibilities
c     are available:

c     1) twidth,wwidth=0 in undecayed matrix elements. This is advisable
c     since it corresponds to what is present in NLO and POWHEG calculation.
c
c     2) twidth,wwidth!=0 in undecayed matrix elements. This corresponds
c     instead to what is used in decayed matrix element.

c     In the following I adopt approach 1) , without removing from
c     MADGRAPH routine for decayed processes all the widths in the
c     propagators that are not resonant (e.g. in every top that is an
c     internal line). I prefer to not include this last trick, since the
c     improvement in the upper bounding evaluation does not justify the
c     work needed to do it.  The problem related to having a QCD
c     emission from the b coming from top decay is also avoided by
c     requiring explicitly MADGRAPH to generate routines which do
c     not include these contributions.



      subroutine ttjdecay
c From a : parton parton -> t tbar parton or: parton parton -> t tbar
c parton parton event on the Les Houches interface, adds the t-tbar
c decay products on the Les Houches interface.  The products of W decays
c are decided by the subroutine
c
c pickwdecays(iwp1,mdecwp1,iwp2,mdecwp2,iwm1,mdecwm1,iwm2,mdecwm2) 
c
c that returns the pdg codes and masses of the decay products of the W+
c and W-.  pickwdecays is initialized separately, at the beginning of
c the run (see inside init_couplings)
      implicit none
      include 'LesHouches.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'madgraph/coupl.inc'
      double precision pp(4,10),pps(4,10),ppss(4,10)
      double precision ppstd(0:3,10)
      double precision beta(3),betacm(3)
      integer itop,itbr,j,k,iwp,iwp1,iwp2,iwm,iwm1,iwm2,mu
      parameter (itop=3,itbr=4)
      double precision mdecwp1,mdecwp2,mdecwm1,mdecwm2
      double precision tpvirt2,tmvirt2,wpvirt2,ebp,wmvirt2,ebm,res_dec
     $     ,res_undec,rat,tten0,cmen0,krtop2,kpair,tten,cmen
      double precision xboundb
      double precision boundnorm,ubound
      integer ini
      data ini/0/
      data boundnorm,ubound/5d0,0d0/
      save ini,boundnorm,ubound
      double precision random,powheginput
      external random,powheginput
      double precision topdecamp,wdecamp,tvirt2,wvirt2,totbr
      logical debug,zerowidth,pwhg_isfinite
      external pwhg_isfinite
      parameter (debug=.true.)
      save zerowidth
c    
      double precision alpha,sin2w
c spin summmed and averaged amplitude squared for top decay
      topdecamp(tvirt2,wvirt2)=-(2*wvirt2**2-tvirt2*wvirt2 -bmass**2
     $     *wvirt2-tvirt2**2+2*bmass**2*tvirt2-bmass**4)/wvirt2 *4*pi
     $     *alpha/(2*sin2w) /2
c spin summmed and averaged amplitude squared for W decay
      wdecamp(wvirt2)=(2*wvirt2) / 3 * 4*pi*alpha/(2*sin2w)
c initialization
      if(ini.eq.2) then
c performs the decay only if a particular channel has been selected
         return
      elseif(ini.eq.0) then
         if(powheginput('#topdecaymode').le.0d0) then
            ini=2
            return
         endif
         zerowidth=powheginput('#zerowidth').eq.1
         if(zerowidth) then
            write(*,*) 
     $    "POWHEG: top and W on mass shell during decay"
         endif
         ini=1
         call set_madgraph_parameters   
      endif
c
      pp(1:4,1:nup)=pup(1:4,1:nup)
      pp(1:4,nup+1:10)=0d0
c pps shall not change from now on, in case one need to recompute
c things because of hit and miss
      pps(1:4,1:nup)=pup(1:4,1:nup)
      pps(1:4,nup+1:10)=0d0
      if(debug) call momcons(2,pp,nup-2,pp(1,3))
c exchange momenta in pp to put them in standard form (E,px,py,pz)
      call momstdform(pp,nup,ppstd)
c local variables 
      alpha=ph_alphaem
      sin2w=ph_sthw2
c set top and w width to zero to evaluate undecayed matrix elements
      twidth=0d0
      wwidth=0d0
c compute undecayed matrix elements for normalization
      if (nup.eq.6) then
         call real_undecayed(ppstd,res_undec)   
      elseif(nup.eq.5) then   
         call born_undecayed(ppstd,res_undec)
      else
         write(*,*)' error: ttdecay, invalid process nup'
         stop
      endif
c check that they are sensible
      if(.not.pwhg_isfinite(res_undec)) then
         write(*,*) "Error in ttjdecay : undecayed matrix elements NaN" 
         call exit(1) 
      endif
c find decay mode of W+ and W-; iwp1, iwp2 refer to the pdg id' of
c W decay products 1 and 2; by convention 1 is down type (e,mu,tau,d,s
c for W-, their antiparticles for W+) and 2 is up type.
c mdecw%% are the corresponding masses of the decay products.
      call pickwdecays(iwp1,mdecwp1,iwp2,mdecwp2, iwm1,mdecwm1,iwm2
     $     ,mdecwm2,totbr)

      rat=0
c Main hit-and-miss return point
 1    continue
c restore saved pp momenta in case of a miss
      pp(1:4,:)=pps(1:4,:)
c boost to partonic CM
      betacm(1:3)=-(pp(1:3,1)+pp(1:3,2))/(pp(4,1)+pp(4,2))
      call boostm(betacm,nup,pp(1,1),pp(1,1))
c ttbar energy in partonic CM (tops on shell!) 
      tten0=pp(4,3)+pp(4,4) 
c ttbar 3-momentum in partonic CM
      kpair=0
      do mu=1,3
         kpair=kpair+(pp(mu,3)+pp(mu,4))**2
      enddo
      kpair=sqrt(kpair)
c total final state energy in partonic CM for on-shell tops
      cmen0=tten0+pp(4,5)
      if (nup.eq.6) cmen0=cmen0+pp(4,6)
c boost heavy quarks in their rest frame 
      beta(1:3)=-(pp(1:3,3)+pp(1:3,4))/tten0
      call boostm(beta,2,pp(1,3),pp(1,3))
c krtop: relative top momentum (i.e. in the top-antitop CM)
      krtop2=pp(4,3)**2-ph_topmass**2
c generate tops and W's virtualities 
      if(.not.zerowidth) then
         call twvirts(krtop2,kpair,cmen0,
     $        tpvirt2,wpvirt2,tmvirt2,wmvirt2)
      else
         tpvirt2=ph_topmass**2
         tmvirt2=ph_topmass**2
         wpvirt2=wmass**2
         wmvirt2=wmass**2
      endif
c rescale top energies according to new invariant masses
      pp(4,3)=sqrt(tpvirt2+krtop2)
      pp(4,4)=sqrt(tmvirt2+krtop2)
c boost heavy quarks back to partonic CM
      if(nup.eq.5) then
         beta(1:3)=-pp(1:3,5)/sqrt(kpair**2+(pp(4,3)
     $        +pp(4,4))**2)
      else
         beta(1:3)=-(pp(1:3,5)+pp(1:3,6))/sqrt(kpair**2+(pp(4,3)
     $        +pp(4,4))**2)
      endif
      call boostm(beta,2,pp(1,3),pp(1,3))
c
c t tbar energy in partonic CM (tops are now off-shell if zerowidth=false)
      tten=pp(4,3)+pp(4,4)
c total final state energy in partonic CM 
      cmen=tten+pp(4,5)
      if (nup.eq.6) cmen=cmen+pp(4,6)
c rescale initial momenta to conserve energy
      pp(1:4,1:2)=pp(1:4,1:2)*cmen/cmen0
c boost back to lab frame
      betacm(:)=-betacm(:)
      call boostm(betacm,nup,pp(1,1),pp(1,1))
c to avoid bugs in HELAS, restore exact masslessness of incoming partons
      pp(4,1:2)=abs(pp(3,1:2))
c reset the values of the Les Houches common block, to account for
c off-shell t and tbar.  We need to do this now since pp is used later
c to store decay products momenta
      pup(1:4,1:nup)=pp(1:4,1:nup)
      pup(5,3)=sqrt(tpvirt2)
      pup(5,4)=sqrt(tmvirt2)
      if(debug) call momcons(2,pp,nup-2,pp(1,3))
c fill the first two entries of the full phase space array
c ppss has the order (parton ,parton -> b,l+,vl,bbar,l-,avl,j,j)
      ppss(1:4,1:2)=pp(1:4,1:2)
c extra partons in the final state fill the last entries
      ppss(1:4,9)=pp(1:4,5)
      if(nup.eq.6) then
         ppss(1:4,10)=pp(1:4,6)
      else
         ppss(1:4,10)=0
      endif
c add t and tbar decay products to ppss
      ebm=(tpvirt2+bmass**2-wpvirt2)/(2*sqrt(tpvirt2))
      ebp=(tmvirt2+bmass**2-wmvirt2)/(2*sqrt(tmvirt2))
      call decaytop(pup(1,3),wpvirt2,ebm,ppss(1,3))
      call decaytop(pup(1,4),wmvirt2,ebp,ppss(1,6))
c check momentum conservation
      if(debug) call momcons(1,pup(1,3),3,ppss(1,3))
      if(debug) call momcons(1,pup(1,4),3,ppss(1,6))
      if(debug) call momcons(2,ppss,nup-4+6,ppss(1,3))
c exchange momenta to get standard form
      call momstdform(ppss,nup-2+6,ppstd) 
c set top and w width different from zero for evaluating decayed matrix
c elements
      twidth=ph_topwidth
      wwidth=ph_wwidth
c compute decayed matrix elements 
      if(nup.eq.5) then
         call born_decayed(ppstd,res_dec)
      else
         call real_decayed(ppstd,res_dec)         
      endif
c check that they are sensible
      if(.not.pwhg_isfinite(res_dec)) then
         res_dec=0
      endif
c evaluate upper bound
      xboundb=1
      xboundb=xboundb/(((tpvirt2-ph_topmass**2)**2+ph_topmass**2
     $     *ph_topwidth**2)* ((wpvirt2-ph_wmass**2)**2+ph_wmass**2
     $     *ph_wwidth**2)* ((tmvirt2-ph_topmass**2)**2+ph_topmass**2
     $     *ph_topwidth**2)* ((wmvirt2-ph_wmass**2)**2+ph_wmass**2
     $     *ph_wwidth**2))
c amplitude squared for top decay
      xboundb=xboundb*topdecamp(tpvirt2,wpvirt2)*topdecamp(tmvirt2
     $     ,wmvirt2)*wdecamp(wpvirt2)*wdecamp(wmvirt2)
c
      if ((xboundb.le.0d0).or.(res_undec).le.0d0) then
         write(*,*) 'subprocess',(idup(j),j=1,nup)   
         write(*,*) 'Error zeron in denominator'
         write(*,*) 'xbound,undec,dec'
         write(*,*) xboundb,res_undec,res_dec
         write(*,*) 'top mass and width'
         write(*,*) ph_topmass,ph_topwidth
         write(*,*) 'W mass and width'
         write(*,*) ph_Wmass,ph_Wwidth
         write(*,*) 'current value of xi and y:',
     $           kn_csitilde,kn_y
         write(*,*) 'current value of tpvirt and wpvirt',
     $           sqrt(tpvirt2),sqrt(wpvirt2)
         write(*,*) 'current value of tmvirt and wmvirt',
     $           sqrt(tmvirt2),sqrt(wmvirt2)
         stop
         endif

      rat=res_dec/(xboundb*res_undec)
      if(rat.gt.boundnorm) then
         call  increasecnt ('ttjdecay upper bound violations')
         if(debug) then
         write(*,*) 'subprocess',(idup(j),j=1,nup)   
         write(*,*) 'topdecay upper bound violation: rat/bound= ',
     $           rat/boundnorm
         write(*,*) 'new decay upper bound is = ',rat
         write(*,*) 'current value of xi and y:',
     $           kn_csitilde,kn_y
         write(*,*) 'current value of tpvirt and wpvirt',
     $           sqrt(tpvirt2),sqrt(wpvirt2)
         write(*,*) 'current value of tmvirt and wmvirt',
     $           sqrt(tmvirt2),sqrt(wmvirt2)
         endif
         boundnorm=rat
      endif
cccc      ubound=max(rat,ubound)
      if(boundnorm*random().gt.rat) then
         call increasecnt('ttjdecay vetoed configurations')
         goto 1
      endif
c End of Hit and Miss loop
c
c check mom. conservation
      if(debug) call momcons(2,ppss,nup-4+6,ppss(1,3))
c write to pp the final momenta 
      pp(:,:)=ppss(:,:)
c Add the decay products on the Les Houches interface.
c Set istup=2 for t,tbar, W+ and W- (intermediate resonance, preserve mass)
      istup(itop)=2
      istup(itbr)=2
c b
      nup=nup+1
      idup(nup)=5
      pup(1:4,nup)=pp(1:4,3)
      pup(5,nup) = bmass
      istup(nup) = 1
      mothup(1,nup)=itop
      mothup(2,nup)=itop
      icolup(1,nup)=icolup(1,itop)
      icolup(2,nup)=icolup(2,itop)
      spinup(nup)=-1

c W+
      nup=nup+1
      iwp=nup
      idup(nup)=24
      pup(1:4,nup)=pp(1:4,4)+pp(1:4,5)
      pup(5,nup)=sqrt(wpvirt2)
      istup(nup)=2
      mothup(1,nup)=itop
      mothup(2,nup)=itop
      icolup(1,nup)=0
      icolup(2,nup)=0
      spinup(nup)=9

c bbar
      nup=nup+1
      idup(nup)=-5
      pup(1:4,nup)=pp(1:4,6)
      pup(5,nup) = bmass
      istup(nup) = 1
      mothup(1,nup)=itbr
      mothup(2,nup)=itbr
      icolup(1,nup)=icolup(1,itbr)
      icolup(2,nup)=icolup(2,itbr)
      spinup(nup)=1


c W-
      nup=nup+1      
      iwm=nup
      idup(nup)=-24
      pup(1:4,nup)=pp(1:4,7)+pp(1:4,8)
      pup(5,nup)=sqrt(wmvirt2)
      istup(nup)=2
      mothup(1,nup)=itbr
      mothup(2,nup)=itbr
      icolup(1,nup)=0
      icolup(2,nup)=0
      spinup(nup)=9

c 1st W+ decay product
      nup=nup+1
      idup(nup)=iwp1
      pup(1:4,nup)=pp(1:4,4)
      pup(5,nup) = mdecwp1
      istup(nup) = 1
      mothup(1,nup)=iwp
      mothup(2,nup)=iwp
      spinup(nup)=1
c setup colours
      if(abs(iwp1).lt.11) then
c find free colour tag
         k=1
         do j=1,nup
            k=max(icolup(1,j),icolup(2,j),k)
         enddo
         k=k+1
c tag found
         if(iwp1.gt.0) then
            icolup(1,nup)=k
            icolup(2,nup)=0
         else
            icolup(1,nup)=0
            icolup(2,nup)=k
         endif
      else
         icolup(1,nup)=0
         icolup(2,nup)=0
      endif

c 2nd W+ decay product
      nup=nup+1
      idup(nup)=iwp2
      pup(1:4,nup)=pp(1:4,5)
      pup(5,nup) = mdecwp2
      istup(nup) = 1
      mothup(1,nup)=iwp
      mothup(2,nup)=iwp
      icolup(1,nup)=icolup(2,nup-1)
      icolup(2,nup)=icolup(1,nup-1)
      spinup(nup)=-1

c 1st W- decay product
      nup=nup+1
      idup(nup)=iwm1
      pup(1:4,nup)=pp(1:4,7)
      pup(5,nup) = mdecwm1
      istup(nup) = 1
      mothup(1,nup)=iwm
      mothup(2,nup)=iwm
      spinup(nup)=-1
c setup colours
      if(abs(iwm1).lt.11) then
c find free colour tag
         k=1
         do j=1,nup
            k=max(icolup(1,j),icolup(2,j),k)
         enddo
         k=k+1
c tag found
         if(iwm1.gt.0) then
            icolup(1,nup)=k
            icolup(2,nup)=0
         else
            icolup(1,nup)=0
            icolup(2,nup)=k
         endif
      else
         icolup(1,nup)=0
         icolup(2,nup)=0
      endif

c 2nd W- decay product
      nup=nup+1
      idup(nup)=iwm2
      pup(1:4,nup)=pp(1:4,8)
      pup(5,nup) = mdecwp2
      istup(nup) = 1
      mothup(1,nup)=iwm
      mothup(2,nup)=iwm
      icolup(1,nup)=icolup(2,nup-1)
      icolup(2,nup)=icolup(1,nup-1)
      spinup(nup)=1
c
      do k=1,nup
         vtimup(k)=0
      enddo
      end


      subroutine born_undecayed(p,res)
      implicit none
      include 'LesHouches.h'
      double precision p(0:3,*),res
      double precision dummy1(nup,nup),dummy2(0:3,0:3,nup)
      integer flav(nup)
      call flavstdform(idup,nup,flav)
      call compborn(p,flav,res,dummy1,dummy2)
      end

      subroutine real_undecayed(p,res)
      implicit none
      include 'LesHouches.h'
      double precision p(0:3,*),res
      integer flav(nup)
      call flavstdform(idup,nup,flav)
      call compreal(p,flav,res) 
      end

      subroutine born_decayed(p,res)
      implicit none
      include 'LesHouches.h'
      double precision p(0:3,*),res
      integer flav(9)
      call flavstdform(idup,nup,flav(1))
      flav(9)=flav(5)
      flav(3)=5
      flav(4)=-11
      flav(5)=12
      flav(6)=-5
      flav(7)=11
      flav(8)=-12
      call compborn_decay(p,flav,res) 
      end

      subroutine real_decayed(p,res)
      implicit none
      include 'LesHouches.h'
      double precision p(0:3,*),res
      integer flav(10)
      call flavstdform(idup,nup,flav(1))
      flav(9)=flav(5)
      flav(10)=flav(6)
      flav(3)=5
      flav(4)=-11
      flav(5)=12
      flav(6)=-5
      flav(7)=11
      flav(8)=-12
      call compreal_decay(p,flav,res) 
      end

      function virt2(mass,gam,xx)
      implicit none
      double precision virt2,mass,gam,xx,xxM
      integer nwidth
      parameter (nwidth=10)
      include "pwhg_math.h"
c generate squared virtuality according to Breit-Wigner and random number xx inside nwidth widths
      xxM=datan(1d0/2d0/nwidth)/pi
      xxM=xxM+xx*(1d0-2d0*xxM)
      virt2=mass-gam/2d0/tan(pi*xxM)
      virt2=virt2*virt2
      return
      end

      subroutine 
     $     twvirts(krtop2,kpair,fsen,tpvirt2,wpvirt2,tmvirt2,wmvirt2)
c Computes the random virtualities for t+, t-, W+ and W-, according to
c their Breit-Wigner shape, to phase space, and to (small) luminosity
c changes due to t and tbar going off the pole mass. The scheme adopted
c here is that the relative 3-momentum of the top in the t-tbar rest
c system, and the 3 momentum of the radiated light partons in the
c partonic CM system are kept fixed. The values of the top and W mass
c and widths, and the value of the b mass, are taken from common blocks.
c
c Input: 
c krtop2: relative 3-momentum squared of t (i.e. t momentum in t+ t-
c rest frame) 
c kpair: t+t- total 3-momentum in partonic CM system
c fsen: total final state energy in partonic CM system    
c
c Output:
c tpvirt2,wpvirt2,tmvirt2,wmvirt2: t+,W+,t-,W- virtualities squared 
c
      implicit none
      include 'PhysPars.h'
      include 'madgraph/coupl.inc'
      double precision krtop2,kpair,fsen,tpvirt2,wpvirt2,tmvirt2,wmvirt2
     $     ,tpvirt,tmvirt,wpvirt,wmvirt
      double precision bmass2,tmass2,wmass2,eb0,ttmass0,tten0,cfac
     $     ,ttmass,kb0,wpkl0,wmkl0,tten,cmen,cmen0,ran(5)
      double precision  pdfcorr,decmom
      external pdfcorr,decmom
      integer j
      double precision pi
      parameter (pi=3.141592653589793d0)
      double precision virt2
      external virt2
      double precision cfacbound
      data cfacbound/2d0/
      save cfacbound
c check consistency
      if(tmass.ne.ph_topmass) then
         write(*,*)
     $ "Error in twvirts: different values in top mass assignment"
         call exit(-1)
      endif
c      
      bmass2=bmass**2
      tmass2=tmass**2
      wmass2=wmass**2
c b energy in t rest frame in zero width limit
      eb0=(tmass2+bmass2-wmass2)/(2*tmass)
c momentum of decay products in t rest frame
      kb0=decmom(tmass,wmass,bmass)
c momentum of decay products in W 
      wpkl0=decmom(wmass,0d0,0d0)
      wmkl0=decmom(wmass,0d0,0d0)
c invariant mass of t+t- system (tops on shell)
      ttmass0=2*sqrt(tmass2+krtop2)
      tten0=sqrt(ttmass0**2+kpair**2)
c total final state energy (tops on shell)
      cmen0=fsen
c Begin loop for hit and miss
 1    continue
c get 5 random numbers
      call rm48(ran,5)
c First get virtualities according to BW;
c avoid too extreme values 
c(default 10 widths around the pole mass defined in virt2)
      tpvirt2=virt2(ph_topmass,ph_topwidth,ran(1))
      tpvirt=sqrt(tpvirt2)
      if((abs(tpvirt-ph_topmass)/ph_topwidth).gt.10) stop "error tpvirt"
      tmvirt2=virt2(ph_topmass,ph_topwidth,ran(2))
      tmvirt=sqrt(tmvirt2)
      if((abs(tmvirt-ph_topmass)/ph_topwidth).gt.10) stop "error tmvirt"
      tmvirt2=virt2(ph_topmass,ph_topwidth,ran(2))
      wpvirt2=virt2(ph_Wmass,ph_Wwidth,ran(3))
      wpvirt=sqrt(wpvirt2)
      if((abs(wpvirt-ph_Wmass)/ph_Wwidth).gt.10) stop "error wpvirt"
      wmvirt2=virt2(ph_Wmass,ph_Wwidth,ran(4))
      wmvirt=sqrt(wmvirt2)
      if((abs(wmvirt-ph_Wmass)/ph_Wwidth).gt.10) stop "error wmvirt"
c estimate phase space factor and do hit and miss with it
c this is from missing 1/(k0(top) * k0(antitop)) in top-antitop rest frame
      cfac=(tmass2+krtop2)/sqrt((tpvirt2+krtop2)*(tmvirt2+krtop2))
c this is from phase space in t->W+b
      cfac=cfac*(decmom(tpvirt,wpvirt,bmass)/tpvirt)/(kb0/tmass)
c this is from phase space in tbar->W- + b
      cfac=cfac*(decmom(tmvirt,wmvirt,bmass)/tmvirt)/(kb0/tmass)
c this is from phase space in W+ decay
      cfac=cfac*(decmom(wpvirt,0d0,0d0)/wpvirt)/(wpkl0/wmass)
c this is from phase space in W- decay
      cfac=cfac*(decmom(wmvirt,0d0,0d0)/wmvirt)/(wmkl0/wmass)
c this is from ttmass/tten*sqrt(tau) factor
      ttmass=sqrt(tpvirt2+krtop2)+sqrt(tmvirt2+krtop2)
      tten=sqrt(ttmass**2+kpair**2)
c total final state energy with off shell tops
      cmen=tten+(fsen-tten0)            
      cfac=cfac * (cmen*ttmass/tten) / (cmen0*ttmass0/tten0)
      cfac=cfac*pdfcorr(cmen/cmen0)
c we optimize for cfac at most = 2 (no prize for larger values)
      if(cfac.gt.cfacbound) then
         call  increasecnt ('twvirts upper bound violations')
         cfacbound=cfac
      endif
      if(cfacbound*ran(5).gt.cfac) then
         call increasecnt('twvirts vetoed configurations')
         goto 1
      endif
      end


      function pdfcorr(corfac)
c compute the correction factor to the luminosity,
c due to the rescaling of the initial momenta when
c the top and anti-top invariant masses are off the pole mass
      implicit none
      include 'LesHouches.h'
      include 'pwhg_pdf.h'
      include 'PhysPars.h'
      double precision pdfcorr,corfac
      double precision mufc2,x(2),xc(2),r
      double precision fx(-6:6)
      double precision powheginput
      external powheginput
      integer ndns(2),ipart(2),ih(2),k
      ndns(1)=pdf_ndns1
      ndns(2)=pdf_ndns2
      ih(1)=pdf_ih1
      ih(2)=pdf_ih2
c Factorization scale fixed to top mass.  Choice is not unambiguos here,
c depending also if a ISR or a FSR occurred. I preferred to use a safe value.
      mufc2=ph_topmass**2

      do k=1,2
         x(k)=pup(4,k)/ebmup(k)
         xc(k)=x(k)*corfac
         if(xc(k).gt.1) then
            pdfcorr=0
            return
         endif
      enddo
      do k=1,2
         if(idbmup(k).eq.2212) then
            ih(k)=1
         elseif(idbmup(k).eq.-2212) then
            ih(k)=-1
         else
            write(*,*) ' pdfcorr: cannot handle incoming hadrons',
     $           idbmup(1),idbmup(2)
            stop
         endif
      enddo
c parton types, from pdg to generic pdf
      do k=1,2
         if(idup(k).eq.21) then
            ipart(k)=0
         else
            ipart(k)=idup(k)
         endif
      enddo
      r=1
      do k=1,2
         call genericpdf(ndns(k),ih(k),mufc2,x(k),fx)
         if(fx(ipart(k)).ne.0) then
            r=r/fx(ipart(k))
            call genericpdf(ndns(k),ih(k),mufc2,xc(k),fx)
            r=r*fx(ipart(k))
         endif
      enddo
      pdfcorr=r
      end

      function decmom(m,m1,m2)
c compute the momentum of the two body decay products, from mass m-> m1,m2
      implicit none
      double precision decmom,m,m1,m2
      if(m.lt.m1+m2) then
         decmom=0
      else
         decmom=sqrt((m-m1-m2)*(m-m1+m2)*(m+m1-m2)*(m+m1+m2))/(2*m)
      endif
      end


      subroutine decaytop(ptop,wvirt2,eb,decprod)
c Put momenta of top decay in decprod, distributed according phase space,
c for a given W virtuality
c Input (momenta as p(1:4)=(p1,p2,p3,E)
c ptop: top 4 momentum
c wvirt2:    W virtuality squared
c eb:   b energy in top rest frame
c m1:   mass of first W decay product
c m2:   mass of second W decay product
c output:
c decprod(1:4,1): bottom momentum
c decprod(1:4,2): W first decay product
c decprod(1:4,3): W second decay product

      implicit none
      include 'madgraph/coupl.inc'
      double precision ptop(4),wvirt2,eb,decprod(4,3)
      double precision wvirt,pb,ew,pl,beta(3)
      wvirt=sqrt(wvirt2)
c modulus of b momentum
      pb=sqrt(eb**2-bmass**2)
c W energy in top rest frame
      ew=sqrt(pb**2+wvirt2)
c generate random directed b momentum in t rest frame
      call rn3vec(decprod(1,1),pb)
      decprod(4,1)=eb
c 1st W decay product in W rest frame
      pl=wvirt/2
      call rn3vec(decprod(1,2),pl)
      decprod(4,2)=pl
c 2nd W decay product in W rest frame
      decprod(1:3,3)=-decprod(1:3,2)
      decprod(4,3)=pl
c boost W decay products with W velocity in t rest frame
c (its momentum is opposite to the b momentum)
      beta(1:3)=-decprod(1:3,1)/ew
      call boostm(beta,2,decprod(1,2),decprod(1,2))
c boost all decay products along top velocity
      beta(1:3)=ptop(1:3)/ptop(4)
      call boostm(beta,3,decprod(1,1),decprod(1,1))
      end


      subroutine pickwdecays(iwp1,mdecwp1,iwp2,mdecwp2,iwm1,mdecwm1
     $     ,iwm2,mdecwm2,totbr)
c Finds which decays to choose with correct probability, according
c to an integer of 5 digits that are either 0, 1 or 2, representing in
c the order the maximum number of the following particles(antiparticles)
c in the final state:
c          e  mu tau up charm
c For example
c 22222    All decays (up to 2 units of everything)
c 20000    both top go into b l nu (with the appropriate signs)
c 10011    one top goes into electron (or positron), the other into (any) hadrons
c 00022    Fully hadronic
c 00002    Fully hadronic with two charms
c 00011    Fully hadronic with a single charm
c 00012    Fully hadronic with at least one charm
      implicit none
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
      integer iwp1,iwp2, iwm1,iwm2, iwp(5,2),iwa(4)
      double precision mdecwp1,mdecwp2,mdecwm1,mdecwm2
      double precision probs(5,5),prbs(25),pr(5),mass(16),sin2cabibbo
     $     ,ebr,hbr,r,totbr
      equivalence (probs(1,1),prbs(1))
      integer ini,ii(5),j,k,imode
      logical semileptonic
      double precision random,powheginput
      external random,powheginput
      data ini/0/
c pdg id's of 1st and 2nd W+ decay products for e,mu,tau,up and charm decays (ignoring CKM)
      data ((iwp(j,k),k=1,2),j=1,5) /-11,12, -13,14, -15,16, -1,2, -3,4/
      save ini,probs,iwp,mass,sin2cabibbo,semileptonic
      if(ini.eq.2) return
      if(ini.eq.0) then
         ini=1
c on first run look for decay mode in powheginput
         imode=powheginput('#topdecaymode')
         semileptonic=powheginput('#semileptonic').eq.1
         if(imode.le.0) then
            ini=2
            return
         endif
         ii(1)=imode/10000
         imode=imode-ii(1)*10000
         ii(2)=imode/1000
         imode=imode-ii(2)*1000
         ii(3)=imode/100
         imode=imode-ii(3)*100
         ii(4)=imode/10
         imode=imode-ii(4)*10
         ii(5)=imode
c     load from input card the branching t->(b l vl) (only one lepton flavour)
         ebr=powheginput('#elbranching')
         if(ebr.lt.0) ebr=0.108
c     from ebr calculates the hadronic branching t->(b u d)
         hbr=(1-3*ebr)/2
c     Probabilities for top decay
         do j=1,5
            if(ii(j).eq.0) then
               pr(j)=0
            else
               if(j.le.3) then
                  pr(j)=ebr
               else
                  pr(j)=hbr
               endif
            endif
         enddo
         do j=1,5
            do k=1,5
               if(j.eq.k.and.ii(k).lt.2) then
                  probs(j,k)=0
               else
                  if(semileptonic.and.( (j.gt.3.and.k.gt.3)
     1                 .or.(j.le.3.and.k.le.3) ) ) then
                     probs(j,k)=0
                  else
                     probs(j,k)=pr(j)*pr(k)
                  endif
               endif
            enddo
         enddo
         do j=2,25
            prbs(j)=prbs(j)+prbs(j-1)
         enddo
         totbr=prbs(25)
         if(totbr.eq.0) then
            write(*,*) 'pickwdecays: the input parameters are such'
            write(*,*) 'that no decays are possible:'
            write(*,*) 'topdecaymode=',nint(powheginput('topdecaymode'))
            if(semileptonic) then
               write(*,*) 'semileptonic=1'
            endif
            write(*,*) ' Halting execution'
            stop
         endif
c     mass of decay products. For internal consistency, here one should
c     use the masses assumed by the shower. Leptonic W decay products
c     masses have to be assigned here. The 3 light quarks are assumed
c     massless. 
         mass(11)=physpar_ml(1)
         mass(13)=physpar_ml(2)
         mass(15)=physpar_ml(3)
         mass(12)=0
         mass(14)=0
         mass(16)=0
         mass(1)=physpar_mq(1)
         mass(2)=physpar_mq(2)
         mass(3)=physpar_mq(3)
         mass(4)=physpar_mq(4)  
         mass(5)=physpar_mq(5)  
CAVEAT: only 2 generation mixing in W decay implemented
         sin2cabibbo=ph_CKM(1,2)**2 
         return
      endif
c end initialization
      r=random()*probs(5,5)
      do k=1,5
         do j=1,5
            if(r.lt.probs(j,k)) goto 1
         enddo
      enddo
 1    continue
c now we have j,k decay mode
      iwa(1)=iwp(j,1)
      iwa(2)=iwp(j,2)
      iwa(3)=-iwp(k,1)
      iwa(4)=-iwp(k,2)
c if any is down or strange, it may turn to
c strange/down with a probability sin^2 theta
      do j=1,4,1
         if(abs(iwa(j)).eq.1) then
            if(random().lt.sin2cabibbo) then
               iwa(j)=sign(3,iwa(j))
            endif
         elseif(abs(iwa(j)).eq.3) then
            if(random().lt.sin2cabibbo) then
               iwa(j)=sign(1,iwa(j))
            endif
         endif
      enddo
      iwp1=iwa(1)
      iwp2=iwa(2)
      iwm1=iwa(3)
      iwm2=iwa(4)
      mdecwp1=mass(abs(iwp1))
      mdecwp2=mass(abs(iwp2))
      mdecwm1=mass(abs(iwm1))
      mdecwm2=mass(abs(iwm2))
      end

      subroutine rn3vec(vec,r)
c Generates a 3d vector in unit sphere
      implicit none
      double precision vec(3),r,r02,norm
      integer j
c generate in unit cube -1<x,y,z<1
 1    call rm48(vec,3)
      do j=1,3
         vec(j)=1-2*vec(j)
      enddo
c generate in unit sphere by hit and miss
      r02=vec(1)**2+vec(2)**2+vec(3)**2
      if(r02.gt.1) goto 1
c normalize
      norm=r/sqrt(r02)
      do j=1,3
         vec(j)=vec(j)*norm
      enddo
      end

c boosts the m 4 vectors vin, output in vout, boost velocity beta
c vin and vout may be the same
      subroutine boostm(beta,m,vin,vout)
      implicit none
      integer m
      double precision beta(3),vin(4,m),vout(4,m)
      double precision betav,gamma
      double precision vdotb
      integer ipart,idim
      betav=sqrt(beta(1)**2 + beta(2)**2 + beta(3)**2)
      gamma=1/sqrt(1-betav**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*beta(1)+vin(2,ipart)*beta(2)+vin(3,ipart)
     $        *beta(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)+beta(idim)/betav*((gamma
     $           -1)*vdotb/betav +gamma*betav*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb)
      enddo
      end

      subroutine flavstdform(fl,n,flav)
c     turns id 21 for gluons into id 0      
      implicit none
      integer n,fl(*),flav(*),k
      do k=1,n
         if (fl(k).eq.21) then
            flav(k)=0
         else
            flav(k)=fl(k)
         endif
      enddo
      end

      subroutine momstdform(pin,n,pout)
c     puts momenta given as (px,py,pz,E) in the standard order
c     (E,px,py,pz), for evaluating matrix elements
      implicit none
      integer n,mu,j
      double precision pin(4,n),pout(0:3,n) 
      pout(0,:)=pin(4,:)
      pout(1:3,:)=pin(1:3,:)
      end

      subroutine momcons(n1,p1,n2,p2)
c check that sum of n1 momenta in p1
c equals the sum of n2 momenta in p2
      implicit none
      integer n1,n2
      double precision p1(4,*),p2(4,*)
      double precision res(4),ini(4),fin(4)
      integer j,k,i
      res=0
      ini=0
      fin=0
      do k=1,n1
      ini=ini+p1(1:4,k)
      enddo
      do k=1,n2
      fin=fin+p2(1:4,k)
      enddo
      res=ini-fin
      do j=1,4
         if(abs(res(j)).gt.1d-7) then
            write(*,*) 'mom. check failed'
            write(*,'(4f15.8)') (res(k),k=1,4)
            write(*,*) 'initial momenta:'
            write(*,'(4f15.8,a)') ((p1(k,i),k=1,4),'\n',i=1,n1)
            write(*,'(a,4f15.8)') 'sum:',ini
            write(*,*) 'final momenta:'
            write(*,'(4f15.8,a)') ((p2(k,i),k=1,4),'\n',i=1,n2)
            write(*,'(a,4f15.8)') 'sum:',fin
            stop " Halting program in ttj decay"
            return
         endif
      enddo
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routines useful for debugging
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine printmom(n,p)
      implicit none
      integer n
      double precision p(4,*)
      double precision sum(4),ini(4),fin(4)
      integer j,k,i
      sum=0
      ini=0
      fin=0
      do k=1,2
      ini=ini+p(1:4,k)
      enddo
      do k=3,n
      fin=fin+p(1:4,k)
      enddo
      sum=ini-fin
      write(*,'(a)') 'momenta:'
      write(*,'(4f15.8,a)') ((p(k,i),k=1,4),'\n',i=1,n)
      write(*,'(a,4f15.8)') 'ini sum:',ini
      write(*,'(a,4f15.8)') 'fin sum:',fin
      write(*,'(a,4f15.8)') 'tot sum:',sum
      end




      subroutine invmasschain(pp)
c print the invarianr masses of the W's and tops in the event
c (for debugging)
      implicit none
      double precision pp(4,10),w(2),t(2),tt
      integer k,i,j
c w virtualities
      do k=1,2
         w(k)=0
c i=0,3
         i=3*(k-1)
         do j=1,3
            w(k)=w(k)-(pp(j,3+i)+pp(j,4+i))**2
         enddo
         w(k)=w(k)+(pp(4,3+i)+pp(4,4+i))**2
      enddo
c t virtualities
      do k=1,2
         t(k)=0
         i=3*(k-1)
         do j=1,3
            t(k)=t(k)-(pp(j,3+i)+pp(j,4+i)+pp(j,5+i))**2
         enddo
         t(k)=t(k)+(pp(4,3+i)+pp(4,4+i)+pp(4,5+i))**2
      enddo
      tt=0
      do j=1,3
         tt=tt-(pp(j,3)+pp(j,4)+pp(j,5)+pp(j,6)+pp(j,7)+pp(j,8))**2
      enddo
      tt=tt+(pp(4,3)+pp(4,4)+pp(4,5)+pp(4,6)+pp(4,7)+pp(4,8))**2
      write(*,*) sqrt(w(1)),sqrt(w(2)),sqrt(t(1)),sqrt(t(2)),sqrt(tt)
      end
