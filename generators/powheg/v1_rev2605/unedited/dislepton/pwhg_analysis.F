c     open some histograms and fill them with data
c     + several auxiliary subroutines

c     book all histograms, will be filled later
      subroutine init_hist
      implicit none
      include 'LesHouches.h'
      include '../pwhg_book.h'
      include 'pwhg_math.h'
      integer diag
      integer max_diag
      parameter (max_diag=200)
      real *8 binsize(max_diag)
      common/pwhghistcommon/binsize
      ! tell this analysis file which program is running it
      character *6 WHCPRG
      common/cWHCPRG/WHCPRG

      call pwhginihist

      ! total sigma
      diag=1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'total cross section','LOG',
     &     binsize(diag),0d0,1d0)

      ! invariant mass of the slepton system
      diag=2
      binsize(diag) = 10d0
      call pwhgbookup(diag,'invariant mass slepton pair','LIN',
     &     binsize(diag),350d0,800d0)

      ! pt(slepton A)
      diag=3
      binsize(diag) = 10d0
      call pwhgbookup(diag,'pt slepton A','LIN',binsize(diag),0d0
     &     ,700d0)

      ! eta(slepton A)
      diag=4
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'eta slepton A','LOG',
     &     binsize(diag),-4d0,4d0)

      ! number of jets
      diag=5
      binsize(diag) = 1d0
      call pwhgbookup(diag,'number of jets','LIN',
     &     binsize(diag),-0.5d0,34.5d0)

      ! number of jets with pt>20GeV
      diag=6
      binsize(diag) = 1d0
      call pwhgbookup(diag,'number of jets cut pt>20GeV','LIN',
     &     binsize(diag),-0.5d0,34.5d0)

      ! pt(jet)
      diag=7
      binsize(diag) = 5d0
      call pwhgbookup(diag,'pt jet 1','LOG',
     &     binsize(diag),0d0,200d0)

      ! y(jet) with pt>20GeV
      diag=8
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y jet 1 cut pt>20GeV','LOG',
     &     binsize(diag),-4d0,4d0)

      if(WHCPRG.eq.'PYTHIA') then
         ! pt(lepton A)
         diag=9
         binsize(diag) = 5d0
         call pwhgbookup(diag,'pt lepton A','LIN',
     &        binsize(diag),0d0,500d0)

         ! pt(lepton A) cuts 1
         diag=10
         binsize(diag) = 5d0
         call pwhgbookup(diag,'pt lepton A cuts 1','LIN',
     &        binsize(diag),0d0,500d0)

         ! eta(lepton A)
         diag=11
         binsize(diag) = 0.2d0
         call pwhgbookup(diag,'eta lepton A','LOG',
     &        binsize(diag),-4d0,4d0)

         ! eta(lepton A) cuts 1
         diag=12
         binsize(diag) = 0.2d0
         call pwhgbookup(diag,'eta lepton A cuts 1','LOG',
     &        binsize(diag),-4d0,4d0)

         ! azimuthal difference of the leptons
         ! in [-pi,pi] measured from lepton A
         diag=13
         binsize(diag) = pi/50d0
         call pwhgbookup(diag,'azimuth difference lepton pair','LIN',
     &        binsize(diag),-pi,pi)

#ifdef HAVE_MT2CODE
         ! MT2
         diag=14
         binsize(diag) = 2.5d0
         call pwhgbookup(diag,'MT2','LIN',
     &        binsize(diag),50d0,200d0)

         ! MT2 cuts 1
         diag=15
         binsize(diag) = 2.5d0
         call pwhgbookup(diag,'MT2 cuts 1','LIN',
     &        binsize(diag),50d0,200d0)

         ! MT2 observed cuts 1
         diag=16
         binsize(diag) = 2.5d0
         call pwhgbookup(diag,'MT2 cuts 1, pt(miss,obs)','LIN',
     &        binsize(diag),50d0,200d0)

         ! MT2 difference observed - ideal
         diag=17
         binsize(diag) = 1d0
         call pwhgbookup(diag,'MT2obs - MT2','LIN',
     &        binsize(diag),-20d0,20d0)
#endif
      endif

      end



c     extract all data required for the histograms
c     calculate quantities that shall be plotted
c     fill histograms
      subroutine analysis(dsig)
      implicit none
      include 'pwhg_math.h'
      include 'LesHouches.h'
      real *8 dsig ! total cross section
      ! tells if subroutine is called for the first time
      logical ini
      data ini/.true./
      save ini
      ! tell this analysis file which program is running it
      character *6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO'/ ! default, should be set when analysis is called
      ! for histograms
      integer diag, max_diag
      parameter (max_diag=200)
      real *8 binsize(max_diag)
      common/pwhghistcommon/binsize
      ! kinematical quantities
      real *8 p_sla(0:3), p_slb(0:3), p_slsl(0:3),
     &        p_la(0:3), p_lb(0:3)
      real *8 pt_sla, pt_la, pt_lb,
     &        eta_sla, eta_la, eta_lb,
     &        y_la, phi_la, phi_lb, phi_ll,
     &        m_slsl
      integer number_n
c     &        pt_slb, pt_slsl, y_slsl, delphi_slsl, mt_slsl,
c     &        p_ll(0:3), m_ll, y_lb, y_ll, pt_nn
      ! jet parameters
      integer j, mjets, num_hard_jets, maxjet
      parameter (maxjet=2048)
      real *8 kt(maxjet), eta(maxjet), rap(maxjet),
     &        phi(maxjet), pj(4,maxjet), ptrel(maxjet)
      real *8 ptmin, yjmax, R
      common /jetcuts/ptmin, yjmax, R
#ifdef HAVE_MT2CODE
      ! for MT2
      real *8 p_nn(0:3)
      real *8 pmissobs(0:2), mt2pa(0:2), mt2pb(0:2)
      real *8 calc_mt2, MT2, MT2obs, m_l, m_n
#endif
      ! IDs
      integer id_sla, id_slb, id_la, id_lb
      ! indices
      integer mu, ij
      ! if allowed decay took place
      logical decayed, isvaliddecayspectrum
      save decayed
      ! for cuts
      logical cuts1passed

      ! at first call
      if (ini) then
         ! check program that runs analysis
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) 'fixed-order NLO analysis'
         elseif(WHCPRG.eq.'LHE') then
            write(*,*) 'LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
            write(*,*) 'not implemented analysis'
            write(*,*) 'no plots will be present at the end of the run'
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write (*,*) '********************************************'

         ! extract slepton types from process id
         call decode_sfermion_pair(lprup(1),id_sla,id_slb)
         write(*,*) 'extracted from process ID: ', lprup(1)
         write(*,*) 'the sfermion pair: ', id_sla, ', ', id_slb

         ! determine lepton IDs in decay slepton -> lepton + neutralino
         id_la=mod(id_sla,100)
         id_lb=mod(id_slb,100)

         ! check sleptons in final state
         if (.not.( mod(id_sla,100).eq.mod(-id_slb,100) .and.
     &        (id_sla.eq.1000011 .or. id_sla.eq.2000011 .or.
     &         id_sla.eq.1000013 .or. id_sla.eq.2000013 .or.
     &         id_sla.eq.1000015 .or. id_sla.eq.2000015))) then
            write(*,*) 'analysis: final state not supported'
            call exit(1)
         endif

         ! write out slepton and lepton IDs
         write (*,*) '********************************************'
         write (*,*) 'in analysis routines'
         write (*,*) 'slepton A:', id_sla
         write (*,*) 'slepton B:', id_slb
         write (*,*) 'lepton A: ', id_la,
     &               ' (from decay of slepton A)'
         write (*,*) 'lepton B: ', id_lb,
     &               ' (from decay of slepton B)'
         write (*,*) '********************************************'

         ! find out if allowed decay took place
         decayed=.true.
         ! determine number of lightest neutralinos
         call getnumber(1000022,number_n)
         ! allowed decay requires lightest neutralinos in final state
         if (number_n.eq.0) then
            decayed=.false.
            write(*,*) 'no lightest neutralino found'
            write(*,*) '-> no allowed slepton decay took place'
         endif
         ! only if we already know that a decay took place
         if (decayed) then
            ! check spectrum if only allowed decays are possible
            if (isvaliddecayspectrum(id_sla) .eqv. .false.) then
               decayed=.false.
            endif
         endif
      endif

      ! from pico to femto
      dsig=dsig*1000

      ! find slepton 4-momenta
      call get4momentum_firstfound(id_sla,p_sla)
      call get4momentum_firstfound(id_slb,p_slb)

      if(decayed) then
         ! decay products
         ! find lepton 4-momenta
         ! undecayed lepton of the same kind as lepton A with largest pt
         call get4momentum_largestpt(id_la,p_la)
         ! undecayed lepton of the same kind as lepton B with largest pt
         call get4momentum_largestpt(id_lb,p_lb)

#ifdef HAVE_MT2CODE
         ! find lepton mass
         call getmass(id_la,m_l)
         ! determine sum of 4-momenta of outgoing neutralinos
         call getsum4momenta(1000022,p_nn)
         ! find neutralino mass
         call getmass(1000022,m_n)
#endif
      endif

      ! calculate slepton kinematical quantities:
      ! slepton A, B transverse momentum
      pt_sla=sqrt(p_sla(1)**2 + p_sla(2)**2)
c      pt_slb=sqrt(p_slb(1)**2 + p_slb(2)**2)
      ! slepton A pseudorapidity
      call get_pseudorap(p_sla,eta_sla)
      ! invariant mass of the slepton system
      do mu=0,3
         p_slsl(mu)=p_slb(mu)+p_sla(mu)
      enddo
      call getinvmass(p_slsl,m_slsl)
      ! transverse momentum of the slepton system
c      pt_slsl=sqrt((p_slb(1)+p_sla(1))**2 +
c     &             (p_slb(2)+p_sla(2))**2)
      ! rapidity of the slepton system
c      call getrapidity(p_slsl,y_slsl)
      ! azimuthal separation between sleptons
c      delphi_slsl = dabs(atan2(p_slb(2),p_slb(1)) -
c     &                   atan2(p_sla(2),p_sla(1)))
c      delphi_slsl=min(delphi_slsl,2d0*pi-delphi_slsl)
      ! transverse mass of the slepton system
c      mt_slsl=sqrt(2*pt_slb*pt_sla*(1d0-dcos(delphi_slsl)))

      ! jets
      ! jet cuts
      ptmin=1d0
      yjmax=4.5d0
      R=0.4d0
      ! initialize
      do j = 1,maxjet
        kt(j) = 0d0
      enddo
      mjets=30
      call buildjets(mjets,kt,eta,rap,phi,pj,ptrel)

      if(decayed) then
         ! calculate lepton kinematical quantities:
         ! lepton A, B transverse momentum
         pt_la=sqrt(p_la(1)**2 + p_la(2)**2)
         pt_lb=sqrt(p_lb(1)**2 + p_lb(2)**2)
         ! rapidity of lepton A, B
         call getrapidity(p_la,y_la)
c         call getrapidity(p_lb,y_lb)
         ! rapidity difference of the leptons
c         y_ll=y_la-y_lb
         ! pseudorapidity of lepton A, B
         call get_pseudorap(p_la,eta_la)
         call get_pseudorap(p_lb,eta_lb)
         ! azimuth of lepton A, B
         call getazimuth(p_la,phi_la)
         call getazimuth(p_lb,phi_lb)
         ! azimuthal difference of the leptons
         phi_ll=phi_lb-phi_la
         if(phi_ll.LT.(-pi)) phi_ll=phi_ll+2d0*pi
         if(phi_ll.GT.pi) phi_ll=phi_ll-2d0*pi
         ! invariant mass of the lepton system
c         do mu=0,3
c            p_ll(mu)=p_lb(mu)+p_la(mu)
c         enddo
c         call getinvmass(p_ll,m_ll)

         ! calculate neutralino kinematical quantities:
         ! transverse momentum of the neutralino system
c         pt_nn=sqrt(p_nn(1)**2 + p_nn(2)**2)

      endif

#ifdef HAVE_MT2CODE
      ! calculate MT2, defined in hep-ph/9906349
      if(decayed) then
         ! construct observed missing pt (ignores soft jets)
         ! assumption: only leptons and jets are observed !
         pmissobs(0) = 0d0 ! ignored
         do mu=1,2
            pmissobs(mu) = - p_la(mu) - p_lb(mu)
         enddo
         do ij = 1,mjets
            if ((kt(ij).gt.10d0).and.(abs(eta(ij)).lt.3)) then
               do mu=1,2
                  pmissobs(mu) = pmissobs(mu) - pj(mu,ij)
               enddo
            endif
         enddo

         mt2pa(0)=m_l
         mt2pb(0)=m_l
         do mu=1,2
            mt2pa(mu)=p_la(mu)
            mt2pb(mu)=p_lb(mu)
         enddo
         ! mt2pa(0:2),mt2pb(0:2),mt2pmiss(0:2)
         ! mt2p...(3) contains (mass,px,py)
         ! for the visible particles and the missing momentum. mt2pmiss(1) is not used.
         ! mn mass of the invisible particle
         ! here we use sum of neutralino momenta for ptmiss:
         MT2    = calc_mt2(mt2pa(0:2), mt2pb(0:2), p_nn(0:2), m_n)
         ! this tries to be more realistic for an observable ptmiss:
         MT2obs = calc_mt2(mt2pa(0:2), mt2pb(0:2), pmissobs(0:2), m_n)
c         write(*,*) "***************** MT2 check:"
c         write(*,*) "la: ", mt2pa(0), "  ", mt2pa(1), "  ", mt2pa(2)
c         write(*,*) "lb: ", mt2pb(0), "  ", mt2pb(1), "  ", mt2pb(2)
c         write(*,*) "n: ", m_n, "  ", p_nn(1), "  ", p_nn(2)
c         write(*,*) "MT2: ", MT2
      endif
#endif

      ! cuts 1
      cuts1passed=.true.
      ! demand decay
      cuts1passed=decayed
      ! lepton pseudorapidity
      if(abs(eta_la).GT.2.5d0) cuts1passed=.false.
      if(abs(eta_lb).GT.2.5d0) cuts1passed=.false.
      ! lepton pt
      if(pt_la.LT.20d0) cuts1passed=.false.
      if(pt_lb.LT.20d0) cuts1passed=.false.

      ! fill histograms
      ! total sigma
      diag=1
      call pwhgfill(diag,0.5d0,dsig/binsize(diag))

      ! invariant mass of the slepton system
      diag=2
      call pwhgfill(diag,m_slsl,dsig/binsize(diag))

      ! pt(slepton A)
      diag=3
      call pwhgfill(diag,pt_sla,dsig/binsize(diag))

      ! eta(slepton A)
      diag=4
      call pwhgfill(diag,eta_sla,dsig/binsize(diag))

      ! number of jets
      diag=5
      call pwhgfill(diag,dble(mjets),1d0)

      ! number of jets with pt>20GeV
      diag=6
      num_hard_jets=0
      do ij = 1,mjets
        if (kt(ij).gt.20d0) num_hard_jets=num_hard_jets+1
      enddo
      call pwhgfill(diag,dble(num_hard_jets),1d0)

      ! pt(jet)
      diag=7
      if(mjets.ge.1) call pwhgfill(diag,kt(1),dsig/binsize(diag))

      ! y(jet) with pt>20GeV
      diag=8
      if((mjets.ge.1).and.(kt(1).gt.20)) then
        call pwhgfill(diag,rap(1),dsig/binsize(diag))
      endif

      if(decayed) then
         ! pt(lepton A)
         diag=9
         call pwhgfill(diag,pt_la,dsig/binsize(diag))

         ! pt(lepton A) cuts 1
         diag=10
         if(cuts1passed) call pwhgfill(diag,pt_la,dsig/binsize(diag))

         ! eta(lepton A)
         diag=11
         call pwhgfill(diag,eta_la,dsig/binsize(diag))

         ! eta(lepton A) cuts 1
         diag=12
         if(cuts1passed) call pwhgfill(diag,eta_la,dsig/binsize(diag))

         ! azimuthal difference of the leptons
         ! in [-pi,pi] measured from lepton A
         diag=13
         call pwhgfill(diag,phi_ll,dsig/binsize(diag))

#ifdef HAVE_MT2CODE
         ! MT2
         diag=14
         if(m_n.ne.0d0) call pwhgfill(diag,MT2,dsig/binsize(diag))

         ! MT2 cuts 1
         diag=15
         if((m_n.ne.0d0).and.cuts1passed) then
            call pwhgfill(diag,MT2,dsig/binsize(diag))
         endif

         ! MT2 observed cuts 1
         diag=16
         if((m_n.ne.0d0).and.cuts1passed) then
            call pwhgfill(diag,MT2obs,dsig/binsize(diag))
         endif

         ! MT2 difference observed - ideal
         diag=17
         if(m_n.ne.0d0) call pwhgfill(diag,MT2obs-MT2,1d0)
#endif
      else
         if (ini) then
            write(*,*) 'no histograms for decay products are generated'
         endif
      endif

      ini=.false.

      end



c     takes PDG ID
c     returns number of undecayed particles with this ID in n
      subroutine getnumber(PDG_ID,n)
      implicit none
      include 'hepevt.h'
      integer PDG_ID, n, ihep

      n=0
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.PDG_ID) then
               n=n+1
            endif
         endif
      enddo

      end



c     takes ihep number of a particle in the HEPEVT common block
c     puts its 4-momentum into p
      subroutine get4momentum_ihep(ihep,p)
      implicit none
      include 'hepevt.h'
      integer ihep, mu
      real *8 p(0:3)

      do mu=0,3
         p(mu)=0d0
      enddo

      if(ihep.lt.0 .or. ihep.gt.nhep) then
         write(*,*) 'no particle with ihep number', ihep,
     &              'in current event'
         write(*,*) 'could not determine four-momentum'
      else
         p(0)=phep(4,ihep)
         do mu=1,3
            p(mu)=phep(mu,ihep)
         enddo
      endif

      end



c     takes particle data group identity code
c     puts 4-momentum of first particle in HEPEVT common block with this ID into p
      subroutine get4momentum_firstfound(PDG_ID,p)
      implicit none
      include 'hepevt.h'
      integer PDG_ID, ihep, mu
      real *8 p(0:3)
      logical found

      do mu=0,3
         p(mu)=0d0
      enddo
      found=.false.

      do ihep=1,nhep
         if(found .eqv. .false.) then
            if(idhep(ihep).eq.PDG_ID) then
               call get4momentum_ihep(ihep,p)
               found=.true.
            endif
         endif
      enddo

      if(found .eqv. .false.) then
         write(*,*) 'particle with PDG ID', PDG_ID, 'not found'
         write(*,*) 'could not determine four-momentum'
      endif

      end



c     takes PDG ID
c     returns four-momentum of undecayed particle with this ID and largest pt
      subroutine get4momentum_largestpt(PDG_ID,p)
      implicit none
      integer PDG_ID, n, nmax, mu
      parameter (nmax=200)
      integer list(nmax)
      real *8 p(0:3)

      do mu=0,3
         p(mu)=0d0
      enddo

      call getptsortedlist(PDG_ID,n,list)
      if (list(1).ne.0) then  ! only if particle with PDG_ID is found
         ! take four-momentum of first particle in list (largest pt)
         call get4momentum_ihep(list(1),p)
      endif

      end



c     takes PDG ID
c     returns number of undecayed particles with this ID in n
c     returns list with ihep numbers of these particles sorted by pt in list
      subroutine getptsortedlist(PDG_ID,n,list)
      implicit none
      include 'hepevt.h'
      integer PDG_ID, n, nmax, ihep
      parameter (nmax=200)
      integer list(nmax)

      do n=1,nmax
         list(n)=0
      enddo

      ! determine number of undecayed particles with PDG_ID in event
      n=0
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.PDG_ID) then
               n=n+1
            endif
         endif
      enddo
      ! check that maximal number of allowed particles is not exceeded
      if(n.gt.nmax) then
         write(*,*) 'more than', nmax, 'particles with ID', PDG_ID
         write(*,*) 'could not fill pt sorted list'
      else
         ! fill list with found particles
         n=0
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.PDG_ID) then
                  n=n+1
                  list(n)=ihep
               endif
            endif
         enddo
         ! sort by pt
         call sortbypt(n,list)
      endif

      end



c     takes particle data group identity code
c     puts sum of all undecayed particles with this ID
c     within the HEPEVT common block into p
      subroutine getsum4momenta(PDG_ID,p)
      implicit none
      include 'hepevt.h'
      integer PDG_ID, ihep, mu
      real *8 p(0:3)
      logical found

      do mu=0,3
         p(mu)=0d0
      enddo
      found=.false.

      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.PDG_ID) then
               p(0)=p(0)+phep(4,ihep)
               do mu=1,3
                  p(mu)=p(mu)+phep(mu,ihep)
               enddo
               found=.true.
            endif
         endif
      enddo

      if(found .eqv. .false.) then
         write(*,*) 'no undecayed particle with PDG ID', PDG_ID
         write(*,*) 'sum of four-momenta 0'
      endif

      end



c     takes particle data group identity code
c     puts mass of first particle in HEPEVT common block with this ID into m
      subroutine getmass(PDG_ID,m)
      implicit none
      include 'hepevt.h'
      integer PDG_ID, ihep
      real *8 m
      logical found

      m=0d0
      found=.false.

      do ihep=1,nhep
         if(found .eqv. .false.) then
            if(idhep(ihep).eq.PDG_ID) then
               m=phep(5,ihep)
               found=.true.
            endif
         endif
      enddo

      if(found .eqv. .false.) then
         write(*,*) 'particle with PDG ID', PDG_ID, 'not found'
         write(*,*) 'could not determine mass'
      endif

      end



c     checks spectrum to ensure that the only possible decay for sleptons is
c     slepton -> lepton + lightest neutralino
c     assumes slhadata to be filled already
c     electron and muon masses neglected, tau mass explicitly set here
      logical function isvaliddecayspectrum(id_sl)
      implicit none
#include "PhysPars.h"
      ! IDs
      integer id_sl
      ! masses
      real *8 m_sl, m_n_light, m_n_next, m_charg, m_l
      ! mixing matrix row (dummy, not used)
      complex *16 Uxl, Uxr

      isvaliddecayspectrum = .true.
      ! invalid default values
      m_sl = invalid
      m_n_light = invalid
      m_n_next = invalid
      m_charg = invalid

      ! extract masses from slhadata
      ! extract slepton mass
      call set_slepton_params(id_sl,m_sl,Uxl,Uxr)
      ! extract mass of lightest neutralino
      m_n_light = real(Mass_MNeu(1))
      ! extract minimum of other neutralino masses
      m_n_next = min(real(Mass_MNeu(2)),
     &               real(Mass_MNeu(3)),
     &               real(Mass_MNeu(4)))
      ! extract minimum of chargino masses
      m_charg = min(real(Mass_MCha(1)),real(Mass_MCha(2)))

      ! check that masses are set correctly
      if ((m_sl.eq.invalid).or.
     &    (m_n_light.eq.invalid).or.
     &    (m_n_next.eq.invalid).or.
     &    (m_charg.eq.invalid)) then
         write(*,*) 'mass parameters for decay are not set'
         isvaliddecayspectrum = .false.
      endif

      ! set charged lepton mass
      select case(abs(id_sl))
      case (1000011) ! selectron left
        m_l = 0d0        ! electron mass neglected
      case (1000013) ! smuon left
        m_l = 0d0        ! muon mass neglected
      case (1000015) ! stau 1
        m_l = 1.77682d0  ! tau mass (PDG 2012)
      case (2000011) ! selectron right
        m_l = 0d0        ! electron mass neglected
      case (2000013) ! smuon right
        m_l = 0d0        ! muon mass neglected
      case (2000015) ! stau 2
        m_l = 1.77682d0  ! tau mass (PDG 2012)
      case default
        write(*,*) 'encountered unknown slepton ID ', id_sl
        stop
      end select

      ! check that decay involving lightest neutralino is possible
      if (m_n_light.gt.(m_sl-m_l)) then
         write(*,*) 'lightest neutralino too heavy for slepton decay'
         isvaliddecayspectrum = .false.
      endif
      ! check that decay involving other neutralinos is not possible
      if (m_n_next.lt.(m_sl-m_l)) then
         write(*,*) 'slepton decay to other neutralinos possible'
         write(*,*) 'unsupported decay channel'
         isvaliddecayspectrum = .false.
      endif
      ! check that decay involving charginos is not possible
      if (m_charg.lt.m_sl) then
         write(*,*) 'slepton decay to charginos possible'
         write(*,*) 'unsupported decay channel'
         isvaliddecayspectrum = .false.
      endif

      ! if spectrum failed check write message
      if (isvaliddecayspectrum .eqv. .false.) then
         write(*,*) 'spectrum not suitable for allowed decay'
         write(*,*) 'actual decays in events generated by PYTHIA'
      endif

      end



c     calculates the separation between p1 and p2
      real *8 function rsep(p1,p2)
      implicit none
      real *8 p1(4),p2(4),r2,dely,delphi,e1,e2

      e1=dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      e2=dsqrt(p2(1)**2+p2(2)**2+p2(3)**2)

      dely = (e1+p1(3))*(e2-p2(3))/
     &      ((e2+p2(3))*(e1-p1(3)))
      dely = 0.5d0*dlog(dely)

      r2= (p1(1)*p2(1)+p1(2)*p2(2))
     &     /dsqrt((p1(1)**2+p1(2)**2)*(p2(1)**2+p2(2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1d0
      if (r2 .lt. -0.9999999D0) r2=-1d0
      delphi=dacos(r2)

      rsep=dsqrt(dely**2+delphi**2)

      end



c     calculates azimuth phi in [0,2pi] from x-axis
      subroutine getazimuth(p,phi)
      implicit none
      include 'pwhg_math.h'
      real *8 p(0:3),phi

      phi=0d0
      if(p(1).GT.0d0) then
         if(p(2).GE.0d0) phi=atan(p(2)/p(1))
         if(p(2).LT.0d0) phi=atan(p(2)/p(1))+2d0*pi
      elseif(p(1).LT.0d0) then
         phi=atan(p(2)/p(1))+pi
      elseif(p(1).EQ.0d0) then
         if(p(2).GT.0d0) phi=pi/2d0
         if(p(2).LT.0d0) phi=-pi/2d0
      endif

      end



      subroutine getrapidity(p,y)
      implicit none
      real *8 p(0:3),y
      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end



      subroutine getinvmass(p,m)
      implicit none
      real *8 p(0:3),m
      m=sqrt(abs((p(0)**2-p(3)**2)-p(1)**2-p(2)**2))
      end



      subroutine get_pseudorap(p,eta)
      implicit none
      real *8 p(0:3),eta,pt,th
      real *8 tiny
      parameter (tiny=1.d-5)

      pt=sqrt(p(1)**2+p(2)**2)
      if(pt.lt.tiny.and.abs(p(3)).lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      elseif(pt.lt.tiny) then   !: added this elseif
         eta=sign(1.d0,p(3))*1.d8
      else
         th=atan2(pt,p(3))
         eta=-log(tan(th/2.d0))
      endif
      end



c     sorts array iarr of ihep numbers by pt (decreasing)
      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include '../include/hepevt.h'
      integer j,k
      real *8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
      ! bubble sort
      touched=.true.
      do while(touched)
         touched=.false.
         do j=1,n-1
            if(pt(j).lt.pt(j+1)) then
               k=iarr(j)
               iarr(j)=iarr(j+1)
               iarr(j+1)=k
               tmp=pt(j)
               pt(j)=pt(j+1)
               pt(j+1)=tmp
               touched=.true.
            endif
         enddo
      enddo
      end



c     jets are orderd by decreasing pt
      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ptrel)
      implicit none
      integer mjets
      real *8  kt(mjets),eta(mjets),rap(mjets),phi(mjets),
     &          pjet(4,mjets),ptrel(mjets)
      include   '../include/hepevt.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu
      real *8 palg,pp,tmp
      real *8 ptmin,yjmax,R
      common /jetcuts/ptmin,yjmax,R
      logical isjetparticle

      ! Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo
      ntracks=0
      do j=1,mjets
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
      ! Extract final state particles to feed to jet finder
      do j=1,nhep
         if (isthep(j).eq.1 .and. isjetparticle(idhep(j))) then
            if(ntracks.eq.maxtrack) then
               write(*,*) 'analyze: need to increase maxtrack!'
               write(*,*) 'ntracks: ',ntracks
               stop
            endif
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,j)
            enddo
            itrackhep(ntracks)=j
         endif
      enddo
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
      ! --------------------------------------------------------------------
      ! -  Inclusive jet pT and Y spectra are to be compared to CDF data:  -
      ! --------------------------------------------------------------------

      ! note: ptmin and R are in common block "jetcuts"

      ! this is for kt-algorithm:
c      call fastjetktwhich(ptrack,ntracks,ptmin,R,
c     &     pjet,njets,jetvec)

      ! anti-kT (c.f. ZZ code):
      palg=-1
c      r=rr
c      ptmin=1d0
      ! jets ordered by decreasing pt (wrapper calls FJ's sorted_by_pt())
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     &                        jetvec)

      mjets=min(mjets,njets)
      if(njets.eq.0) return
      ! check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
      ! end check consistency

      ! ---------------------------------------------------------------------
      ! - Computing arrays of useful kinematics quantities for hardest jets -
      ! ---------------------------------------------------------------------
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

      call computeptrel(ptrack,ntracks,rap,kt,phi,mjets,jetvec,ptrel)

      end



c     determines if particle contributes to a jet
      logical function isjetparticle(PDG_ID)
      implicit none
      integer PDG_ID
      logical injetwhitelist, injetblacklist

      isjetparticle=.false.

      ! find out if particle is in
      ! jetwhitelist: particles that contribute to a jet
      ! jetblacklist: particles that do not contribute to a jet
      if (injetwhitelist(PDG_ID)) then
         isjetparticle=.true.
      else if (injetblacklist(PDG_ID)) then
         isjetparticle=.false.
      else
         write(*,*) 'could not determine if particle with ID', PDG_ID
         write(*,*) 'contributes to a jet'
         stop
      endif

      end



c     determines if particle is in whitelist
c     of particles that contribute to a jet
      logical function injetwhitelist(PDG_ID)
      implicit none
      integer PDG_ID, id

      injetwhitelist=.false.

      id=iabs(PDG_ID)

      ! the following particles contribute to a jet
      ! quarks
      if (1.le.id .and. id.le.6) injetwhitelist=.true.
      ! gluon
      if (id.eq.21) injetwhitelist=.true.
      if (id.eq.0) injetwhitelist=.true.  ! POWHEG numbering
      ! diquarks
      if ((1000.le.id .and. id.le.6000)
     &    .and. (mod(id,100).eq.1 .or. mod(id,100).eq.3)) then
         injetwhitelist=.true.
      endif
      ! mesons
      if ((111.le.id .and. id.le.600)
     &    .or. (10000.le.id .and. id.le.999999)
     &    .or. (9000000.le.id .and. id.le.9100000)) then
         injetwhitelist=.true.
      endif
      ! baryons
      if ((1000.le.id .and. id.le.6000)
     &    .and. (mod(id,10).eq.2 .or. mod(id,10).eq.4)) then
         injetwhitelist=.true.
      endif

      end



c     determines if particle is in blacklist
c     of particles that do not contribute to a jet
      logical function injetblacklist(PDG_ID)
      implicit none
      integer PDG_ID, id

      injetblacklist=.false.

      id=iabs(PDG_ID)

      ! the following particles may occur
      ! but do not contribute to a jet
      ! leptons
      if (id.eq.11) injetblacklist=.true.  ! e
      if (id.eq.12) injetblacklist=.true.  ! nue
      if (id.eq.13) injetblacklist=.true.  ! mu
      if (id.eq.14) injetblacklist=.true.  ! numu
      if (id.eq.15) injetblacklist=.true.  ! tau
      if (id.eq.16) injetblacklist=.true.  ! nu_tau
      ! sleptons
      if (id.eq.1000011) injetblacklist=.true.  ! ~e_L
      if (id.eq.1000012) injetblacklist=.true.  ! ~snue_L
      if (id.eq.1000013) injetblacklist=.true.  ! ~mu_L
      if (id.eq.1000014) injetblacklist=.true.  ! ~snumu_L
      if (id.eq.1000015) injetblacklist=.true.  ! ~stau_1
      if (id.eq.1000016) injetblacklist=.true.  ! ~snu_tau_L
      if (id.eq.2000011) injetblacklist=.true.  ! ~e_R
      if (id.eq.2000013) injetblacklist=.true.  ! ~mu_R
      if (id.eq.2000015) injetblacklist=.true.  ! ~stau_2
      ! lightest neutralino
      if (id.eq.1000022) injetblacklist=.true.  ! ~neutralino(1)

      end



c     pt rel is sum of the ptrack momentum projection ortogonal to the jet
c     momentum in the frame where the jet has zero rapidity
      subroutine computeptrel(ptracks,ntracks,rapjets,ktjets,phijets,
     &     njets,jetvec,ptrel)
      implicit none
      integer ntracks,njets,jetvec(ntracks)
      real *8 ptracks(4,ntracks),rapjets(njets),
     &     ktjets(njets),phijets(njets),ptrel(njets)
      integer j,i
      real *8 yj,kj1,kj2,y,pt(3)

      do j=1,njets
         ptrel(j)=0
      enddo
      do i=1,ntracks
         j=jetvec(i)
         if(j.gt.0.and.j.le.njets) then
            ! Track i belongs to jet j
            yj=rapjets(j)
            kj1=ktjets(j)*cos(phijets(j))
            kj2=ktjets(j)*sin(phijets(j))
            ! rapidity of track i
            y=0.5d0*log((ptracks(4,i)+ptracks(3,i))
     &                 /(ptracks(4,i)-ptracks(3,i)))
            ! rapidity of track i in frame where the jet has zero rapidity
            y=y-yj
            ! find momentum of track i in frame where the jet has zero rapidity
            pt(1)=ptracks(1,i)
            pt(2)=ptracks(2,i)
            pt(3)=sqrt(pt(1)**2+pt(2)**2)*sinh(y)
            ! pt rel is sum of the ptrack momentum projection ortogonal to the jet
            ! momentum in the frame where the jet has zero rapidity
            ptrel(j)=sqrt(((pt(1)*kj2-pt(2)*kj1)**2+
     &                     (         -pt(3)*kj2)**2+
     &                     (pt(3)*kj1          )**2)/
     &                     (kj1**2+kj2**2)) + ptrel(j)
         endif
      enddo

      end
