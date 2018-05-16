c
c returns Born (nlo=0) or 
c finite part of virtual contributions (NOT sum of Born + virtuals!) 
c      
      subroutine qqwwqq_vonly_channel(pbar,sign,nlo,L,k,ans)
      implicit none
c
c	Initial version by Dieter Zeppenfeld (vbfnlo)
c	Last modified for POWHEG: 2012 Dec. 
C
C  QQWWQQ calculates the matrix elements**2 for electroweak
c  weak boson pair production by quark quark scattering
C
C        q1 q3    ---->   q2 q4 W+W-,   W+ ---> f5-bar f6, W- ---> f7-bar f8
C
c  where the final state leptons are 2 charged leptons + 2 neutrinos. 
C  Crossing related processes can be computed as well. Pauli interference 
c  terms for identical fermions are neglected. In particular, only the
c  t-channel exchange of electroweak bosons is considered. s-channel
c  production of 3 weak bosons is NOT implemented.
C
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi**2)
      include '../boxfiles-pre2-1/pwhg_st.h'

      include 'global.inc'
      include 'tensor.inc'
      include '../higgs_graphs.h'
c
c electroweak couplings:
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
c
c alfas, scales etc
      include 'scales.inc'
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+nv), fac, musq
      double precision res(6),resv(6)
      double precision ans
      integer kl
      double precision  p(0:3,4+nv), p21(0:4), p43(0:4), pwp(0:4),
     1                  pwm(0:4), pww(0:4)
      integer  sign(4+nv), fs(4+nv), nlo, mu, i, j, jj, k, kk, id,
     1         isig, isig1, isig3
      integer  ifl(4,6), js1, js3, L, is1, is3
      double complex prop21(4), prop43(4), propwp(4)
      double complex mat(6,-1:1,-1:1,9), matv(6,-1:1,-1:1,9)
      double complex mm(6,-1:1,-1:1), 
     1               mv12(6,-1:1,-1:1), mv34(6,-1:1,-1:1)
      double complex maa, maz, mza, mzz, mww5, mww6
      double complex  m1, m2, ma(2), mz(2), mpm(4), mmp(4), m5(3,3:4)
      double complex mv1,mv2,mva(2),mvz(2), mv5(3,3:4)
      double complex psi(2,-1:1,4), jqq(0:5,-1:1,2), eps(0:5)
      double complex psiwp(2,4), psiwm(2,4), bkjqq(2,-1:1,-1:1,4),
     1               jwp(0:5,-1:1,4), jwm(0:5,-1:1,4),
     1               jvwp(0:5,-1:1,-1:1,4), jvwm(0:5,-1:1,-1:1,4),
     2               qepswm(6), qepswp(6)
      double complex ja(0:5,-1:1,-1:1,4), jz(0:5,-1:1,-1:1,4), 
     1               jva(0:5,-1:1,-1:1,4), jvz(0:5,-1:1,-1:1,4),
     2               psia(2,-1:1,-1:1,4), psiz(2,-1:1,-1:1,4),
     3               mvpm(4,2), mvmp(4,2)
      double complex j5pm(0:3,3,-1:1,2), j5mp(0:3,3,-1:1,2)
      double complex j5pma(0:3,3,-1:1,2), j5mpa(0:3,3,-1:1,2)
      double complex j5pmb(0:3,3,-1:1,2), j5mpb(0:3,3,-1:1,2)
      double complex jbornp(0:5), jbornm(0:5)
      double complex epsa(0:3), epsz(0:3), 
     1               epsNCwp(0:5,-1:1,3:4,2), epsNCwm(0:5,-1:1,3:4,2),
     2               epsCCwp(0:5,-1:1,3:4,2), epsCCwm(0:5,-1:1,3:4,2)
      double precision fqp(0:4,4), fqm(0:4,4), fq(0:4,4),
     1                 qpm(0:4), qmp(0:4), dummy(0:4)
      double complex jj21m,jj21p,jj43m,jj43p, zpm, zmp, zm2i(2:4)
      double complex z1,zp(4:5),zm(4:5), qepsCCwp(2), qepsCCwm(2)
      double complex propt(-1:1,-1:1,6,2), prop_pm(3), prop_mp(3)
      double complex contract_Tjj, dotcc, dotrc, dotqj, s1c
      external contract_Tjj, dotcc, dotrc, dotqj, s1c

      logical nc_type,cc_type

      logical ldebug
      parameter (ldebug = .false.)

      logical linit,lgc(4)
      data linit /.true./
      data lgc /4*.false./
      data lfs /4*.true./

      integer*8 icb1, icount1, icb2, icount2
      data icb1/0/,icount1/0/,icb2/0/,icount2/0/
      double precision xv1, xv2, xgc1, xgc2
      save ifl, zm2i, linit, lgc, icb1, icb2, icount1, icount2
      save ja,jz,jva,jvz,psia,psiz,fq,j5pm,j5mp
      double complex wpx(0:5), wmx(0:5), z2, zero
      parameter (zero = (0d0,0d0) )

      integer icountmax
      common/gauge_checks/icountmax 

      real *8 powheginput
      external powheginput 
      integer ncall2
      logical, save :: firsttime = .true.   
c
c variables for powheg:
      double precision q2_up,q2_lo,rup,rlo,lrup,lrlo
      double precision cvirtl   

c identify "bad" points (low qsq):
      logical qdamp      
c
c variables for virtual corrections:
	logical bad_gauge,bad_gauge_sin
	common /gauge / bad_gauge,bad_gauge_sin
c
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpent, lpt, lpq
      lnlo = NLO.ne.0    ! include some virtual stuff if T
      lbox = NLO.eq.1 .or. NLO.eq.-4 .or. NLO.eq.5    ! call boxline if T
      lpt  = NLO.eq.1 .or. NLO.le.-5                  ! true pentagon contributions on
      lpq  = NLO.eq.5 .or. NLO.eq.-4                  ! q^mu/m_V terms of Pentagon contributions only
      lpent= lpt .or. lpq                             ! include pentagons if T
c
      bad_gauge = .false. ! set F at beginning of each run
c
c----------------------------------------------------------
c----------------------------------------------------------
c	
      if (firsttime) then
        ncall2=powheginput('ncall2')
        firsttime = .false.
      endif

c define flavors of external quarks for the 4 NC and 2 CC subprocesses
c
      if (linit) then
         linit = .false.
c
         print*,'Born/virtual amplitudes:'
         if (higgs_only) print*,'H->WW graphs only'
         print*, 'minimum virtuality for t-channel photon exchange'
         print*, 'qsqA_min = ',qsqAmin,'GeV**2'
         print*, 'damping factor of 1d-20 below qsqAmin '
c
         kl = 1                  ! uucc
         ifl(1,kl) = 3
         ifl(2,kl) = 3
         ifl(3,kl) = 3
         ifl(4,kl) = 3
         kl = 2                  ! uuss
         ifl(1,kl) = 3
         ifl(2,kl) = 3
         ifl(3,kl) = 4
         ifl(4,kl) = 4
         kl = 3                  ! ddcc
         ifl(1,kl) = 4
         ifl(2,kl) = 4
         ifl(3,kl) = 3
         ifl(4,kl) = 3
         kl = 4                  ! ddss
         ifl(1,kl) = 4
         ifl(2,kl) = 4
         ifl(3,kl) = 4
         ifl(4,kl) = 4
         kl = 5                  ! udsc
         ifl(1,kl) = 3
         ifl(2,kl) = 4
         ifl(3,kl) = 4
         ifl(4,kl) = 3
         kl =6                   ! ducs
         ifl(1,kl) = 4
         ifl(2,kl) = 3
         ifl(3,kl) = 3
         ifl(4,kl) = 4
         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))	 

         if (ncall2.gt.10000) then 
            icountmax = ncall2
         else
            icountmax = 100000
	 endif

      endif

      if (k.le.4) then
         nc_type = .true.
         cc_type = .false.
      else
         cc_type = .true.
         nc_type = .false.
      endif

      do kl = 1,6
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               do i = 1,9
                  mat(kl,isig1,isig3,i) = 0
                  matv(kl,isig1,isig3,i) = 0
               enddo
            enddo
         enddo
      enddo
c
c identify fermion line sign factors
c
      is1 = sign(1)
      is3 = sign(3)
      js1 = (3+sign(1))/2       ! 1 for sign1=-1,2 for sign1=+1
      js3 = (7+sign(3))/2       ! 3 for sign3=-1,4 for sign3=+1
c
c define the internal momenta
c
      do mu = 0,3
         do i = 1,4+nv
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         p21(mu) = p(mu,2) - p(mu,1)
         p43(mu) = p(mu,4) - p(mu,3)
         pwp(mu) = p(mu,6) - p(mu,5)
         pwm(mu) = p(mu,8) - p(mu,7)
         pww(mu) = pwp(mu) + pwm(mu)
      enddo
      p21(4) = p21(0)**2 - p21(1)**2 - p21(2)**2 - p21(3)**2
      p43(4) = p43(0)**2 - p43(1)**2 - p43(2)**2 - p43(3)**2
      pwp(4) = pwp(0)**2 - pwp(1)**2 - pwp(2)**2 - pwp(3)**2
      pwm(4) = pwm(0)**2 - pwm(1)**2 - pwm(2)**2 - pwm(3)**2
      pww(4) = pww(0)**2 - pww(1)**2 - pww(2)**2 - pww(3)**2   

      qdamp = .false.
c eliminate contributions from too low qsq in t-channel:
      if ( (abs(p21(4)).lt.qsqAmin).or.(abs(p43(4)).lt.qsqAmin)) then
      	qdamp = .true.
      endif
      
c
c get the vector boson propagator factors
c
      prop21(1) = 1/p21(4)
      prop21(2) = 1/dcmplx(p21(4)-xm2(2),xmg(2))
      prop21(3) = 1/dcmplx(p21(4)-xm2(3),xmg(3))
      prop21(4) = prop21(3)

      prop43(1) = 1/p43(4)
      prop43(2) = 1/dcmplx(p43(4)-xm2(2),xmg(2))
      prop43(3) = 1/dcmplx(p43(4)-xm2(3),xmg(3))
      prop43(4) = prop43(3)
c
c for box-box and BV graphs we need the propagators for t-channel bosons between quark lines
c as seen from upper line these W momenta are INCOMING. They are OUTGOINg as seen from lower 
c line
      do mu = 0,3
         qpm(mu) = pwp(mu)+p(mu,2) - p(mu,1)    !W+ emitted on upper line
         qmp(mu) = pwm(mu)+p(mu,2) - p(mu,1)    !W- emitted on upper line
      enddo
      qpm(4) = qpm(0)**2-qpm(1)**2-qpm(2)**2-qpm(3)**2
      qmp(4) = qmp(0)**2-qmp(1)**2-qmp(2)**2-qmp(3)**2

      prop_pm(1) = clr(3,3,-1)**2/qpm(4)
      prop_pm(2) = clr(3,3,-1)**2/dcmplx(qpm(4)-xm2(2),xmg(2))
      prop_pm(3) = clr(3,3,-1)**4/dcmplx(qpm(4)-xm2(3),xmg(3))
      prop_mp(1) = clr(3,3,-1)**2/qmp(4)
      prop_mp(2) = clr(3,3,-1)**2/dcmplx(qmp(4)-xm2(2),xmg(2))
      prop_mp(3) = clr(3,3,-1)**4/dcmplx(qmp(4)-xm2(3),xmg(3))
c
c get the external quark spinors (including factor sqrt(2E) )
c
      call psi0m(4,pbar(0,1),sign(1),psi)
c
c get the f-fbar currents J21^mu=jqq(mu,*,1), J43^mu=jqq(mu,*,2) 
c
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
      call curr6(1,psi(1,-1,4),p(0,4),psi(1,-1,3),p(0,3),jqq(0,-1,2))
c
c --------------------------------------------
c
c contract with vvtoww tensors to get Vertex-Vertex scattering diagrams
c
      if (nc_type) then
      do isig1 = -1,1,2
         do isig3 = -1,1,2
            maa = contract_Tjj(aaww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            maz = contract_Tjj(azww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
            mza = contract_Tjj(zaww(0,0,L),jqq(0,isig3,2),jqq(0,isig1,1))
            mzz = contract_Tjj(zzww(0,0,L),jqq(0,isig1,1),jqq(0,isig3,2))
c            do k = 1,4
               mat(k,isig1,isig3,1) = 
     1              maa*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),1,isig3)
     2            + maz*clr(ifl(1,k),1,isig1)*clr(ifl(3,k),2,isig3)
     3            + mza*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),1,isig3)
     4            + mzz*clr(ifl(1,k),2,isig1)*clr(ifl(3,k),2,isig3)
               matv(k,isig1,isig3,1) = (0d0,0d0)
c            enddo
            
         enddo
      enddo
      else ! cc_type
         if (k.eq.5) then
            mww5 = contract_Tjj(wwww5(0,0,L),jqq(0,-1,2),jqq(0,-1,1))
            mat(5,-1,-1,1)  = mww5*clr(3,3,-1)**2
            matv(5,-1,-1,1) = (0d0,0d0)     
         elseif (k.eq.6) then   
            mww6 = contract_Tjj(wwww6(0,0,L),jqq(0,-1,1),jqq(0,-1,2))
            mat(6,-1,-1,1)  = mww6*clr(3,3,-1)**2
            matv(6,-1,-1,1) = (0d0,0d0)
         endif   
      endif !nc/cc

      if (higgs_only) goto 999
c
c --------------------------------------------      
c
c  prepare box diagrams: attach W+ and W- to external spinors
c 
c      isig = -1   : lefthanded spinors only coupling to W's
c
c  Notation for virtual 2-component spinors and momenta
c
c  W+ attached to quark number i: psiwp(*,i) with momentum fqp(mu,i)
c  W- attached to quark number i: psiwm(*,i) with momentum fqm(mu,i)
c  
c  the fermion current corresponding to a quark line with the real emitted W+
c  attached next to quark number i is stored in jwp(mu,isig,i). Similarly
c  jwm(mu,isig,i) is the corresponding current for a W- attached next to quark i
c 
c  For the virtual amplitudes the notation, e.g.     jvwp(mu,isig,is,i)
c  is used for the boxline correction to a quark line with one W+ attached next to
c  quark #i and a free Lorentz index mu for the second attached EW boson. is=+-1
c  refers to the sign factor of this quark (vs.antiquark line). They are recalculated 
c  only if this quark line sign has not been calculated yet for this phase space point
c  (i.e. lfs = .true.)  Otherwise they are taken from saved previous calculation
 
      do i = 1,3,2
         call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qp,wp,
     1              psiwp(1,i),fqp(0,i))
         call ket2c(psi(1,-1,i),.true.,p(0,i),-1,qm,wm,
     1              psiwm(1,i),fqm(0,i))
         call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qp,wp,
     1              psiwp(1,i+1),fqp(0,i+1))
         call bra2c(psi(1,-1,i+1),.true.,p(0,i+1),-1,qm,wm,
     1              psiwm(1,i+1),fqm(0,i+1))
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1                 psiwp(1,i),fqp(0,i), jwp(0,-1,i)   )
         call curr6(-1,psiwp(1,i+1),fqp(0,i+1),
     1                 psi(1,-1,i),p(0,i), jwp(0,-1,i+1) )
         call curr6(-1,psi(1,-1,i+1),p(0,i+1),
     1                 psiwm(1,i),fqm(0,i), jwm(0,-1,i)   )
         call curr6(-1,psiwm(1,i+1),fqm(0,i+1),
     1                 psi(1,-1,i),p(0,i), jwm(0,-1,i+1) )
      enddo
      if (lbox .and. lfs(js1).and.(.not.qdamp)) then
         do i = 1,2
            call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., wp,3-i,   jwp(0,-1,i), 
     2                    jvwp(0,-1,is1,i) )
            call boxlinec(-1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                    .true., wm,3-i,   jwm(0,-1,i), 
     2                    jvwm(0,-1,is1,i) )
         enddo
      endif
      if (lbox .and. lfs(js3).and.(.not.qdamp)) then
         do i = 3,4
            call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., wp,5-i,   jwp(0,-1,i), 
     2                    jvwp(0,-1,is3,i) )
            call boxlinec(-1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                    .true., wm,5-i,   jwm(0,-1,i), 
     2                    jvwm(0,-1,is3,i) )
         enddo
      endif
c
c now calculate the Vertex-box diagrams; get t-channel W currents first
      do kl = 1,2
         kk = 3-kl
c kl=1 and kk=2 is for "box correction" to upper line
c kl=2 and kk=1 is for "box correction" to lower line
         if (kl.eq.1) then
            zp(4) = -dcmplx(qmp(0),qmp(3))
            zp(5) = -dcmplx(qmp(1),qmp(2))
            zm(4) = -dcmplx(qpm(0),qpm(3))
            zm(5) = -dcmplx(qpm(1),qpm(2))
         else
            zm(4) = dcmplx(qmp(0),qmp(3))
            zm(5) = dcmplx(qmp(1),qmp(2))
            zp(4) = dcmplx(qpm(0),qpm(3))
            zp(5) = dcmplx(qpm(1),qpm(2))
         endif
         if (nc_type) then 
         do isig = -1,1,2
            call contract_T2j(NCwpa(0,0,kl,L),jqq(0,isig,kk),   epsa) !for W+ current
            call contract_T2j(NCwpz(0,0,kl,L),jqq(0,isig,kk),   epsz)
            do mu = 0,3
               epsNCwp(mu,isig,3,kl) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwp(mu,isig,4,kl) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo
            call contract_T2j(NCwma(0,0,kl,L),jqq(0,isig,kk),   epsa) !for W- current
            call contract_T2j(NCwmz(0,0,kl,L),jqq(0,isig,kk),   epsz)
            do mu = 0,3 
               epsNCwm(mu,isig,3,kl) =
     1              epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
               epsNCwm(mu,isig,4,kl) =
     1              epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
            enddo
            do mu = 4,5           ! add momentum info to the currents
               epsNCwp(mu,isig,3,kl) = zp(mu)
               epsNCwp(mu,isig,4,kl) = zp(mu)
               epsNCwm(mu,isig,3,kl) = zm(mu)
               epsNCwm(mu,isig,4,kl) = zm(mu)
            enddo
         enddo

         endif !nc

c and same for the CC processes (W attached to j43 or j21 current)
         if (cc_type) then        
         isig = -1
         call contract_T1j(CCwpa(0,0,kl,L),jqq(0,isig,kk),   epsa) !for W+ current
         call contract_T1j(CCwpz(0,0,kl,L),jqq(0,isig,kk),   epsz)
         do mu = 0,3 
            epsCCwp(mu,isig,3,kl) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwp(mu,isig,4,kl) =
     1           epsa(mu)*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
         enddo
         if (kl.eq.1) then
            qepsCCwp(kl) = -dotrc(qmp,epsz)*zm2i(2)
         else
            qepsCCwp(kl) = dotrc(qpm,epsz)*zm2i(2)
         endif
         call contract_T1j(CCwma(0,0,kl,L),jqq(0,isig,kk),   epsa) !for W- current
         call contract_T1j(CCwmz(0,0,kl,L),jqq(0,isig,kk),   epsz)
         do mu = 0,3 
            epsCCwm(mu,isig,3,kl) =
     1           epsa(mu)*clr(3,1,isig)+epsz(mu)*clr(3,2,isig)
            epsCCwm(mu,isig,4,kl) =
     1           epsa(mu )*clr(4,1,isig)+epsz(mu)*clr(4,2,isig)
         enddo
         do mu = 4,5            ! add momentum info to the currents
            epsCCwp(mu,isig,3,kl) = zp(mu)
            epsCCwp(mu,isig,4,kl) = zp(mu)
            epsCCwm(mu,isig,3,kl) = zm(mu)
            epsCCwm(mu,isig,4,kl) = zm(mu)
         enddo 
         if (kl.eq.1) then
            qepsCCwm(kl) = -dotrc(qpm,epsz)*zm2i(2)
         else
            qepsCCwm(kl) = dotrc(qmp,epsz)*zm2i(2)
         endif
      endif !cc   
      enddo !kl
      
      jj21m = dotcc(wm,jqq(0,-1,1))
      jj21p = dotcc(wp,jqq(0,-1,1))
      jj43m = dotcc(wm,jqq(0,-1,2))
      jj43p = dotcc(wp,jqq(0,-1,2))

c now construct the contribution to the amplitude by current contraction
c virtual contributions are assembled in subroutine boxline
      if (nc_type) then   
      do isig = -1,1,2            ! 2 bosons attached to 12 line
         do kl=3,4
            qepswm(kl) = -dotrc(qpm,epsNCwm(0,isig,kl,1))*zm2i(3)
            qepswp(kl) = -dotrc(qmp,epsNCwp(0,isig,kl,1))*zm2i(3)
         enddo
c         do k = 1,2                !uucc, uuss
         if (k.eq.1.or.k.eq.2) then
            m1 = dotcc(jwp(0,-1,1),epsNCwm(0,isig,ifl(3,k),1))      !uW+(tW-)u   
     1         - qepswm(ifl(3,k))*jj21p
            m2 = dotcc(jwm(0,-1,2),epsNCwp(0,isig,ifl(3,k),1))      !u(tW+)W-u
     1         + qepswp(ifl(3,k))*jj21m
            mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2)
            if (lbox.and.(.not.qdamp)) then
               mv1 = dotcc(jvwp(0,-1,is1,1),epsNCwm(0,isig,ifl(3,k),1)) !uW+(tW-)u   
               mv2 = dotcc(jvwm(0,-1,is1,2),epsNCwp(0,isig,ifl(3,k),1)) !u(tW+)W-u
               matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
c         enddo
         elseif (k.eq.3.or.k.eq.4) then   
c         do k = 3,4             !ddcc, ddss
            m1 = dotcc(jwp(0,-1,2),epsNCwm(0,isig,ifl(3,k),1))      !d(tW-)W+d
     1         + qepswm(ifl(3,k))*jj21p
            m2 = dotcc(jwm(0,-1,1),epsNCwp(0,isig,ifl(3,k),1))      !dW-(tW+)d
     1         - qepswp(ifl(3,k))*jj21m
            mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2)
            if (lbox.and.(.not.qdamp)) then 
               mv1 = dotcc(jvwp(0,-1,is1,2),epsNCwm(0,isig,ifl(3,k),1)) !d(tW-)W+d
               mv2 = dotcc(jvwm(0,-1,is1,1),epsNCwp(0,isig,ifl(3,k),1)) !dW-(tW+)d   
               matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
c         enddo
          endif !k  
      enddo !isig

      else !CC
      isig = -1
c      k = 5                     !udsc
      if (k.eq.5) then
      m1 = dotcc(jwp(0,-1,2),epsCCwm(0,isig,ifl(1,k),1))            !u(tAZ)W+d
      m2 = dotcc(jwp(0,-1,1),epsCCwm(0,isig,ifl(2,k),1))            !uW+(tAZ)d
      z1 =qepsCCwm(1)*jj21p*(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
      mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2+z1)
      if (lbox.and.(.not.qdamp)) then
         mv1 = dotcc(jvwp(0,-1,is1,2),epsCCwm(0,isig,ifl(1,k),1)) !u(tAZ)W+d
         mv2 = dotcc(jvwp(0,-1,is1,1),epsCCwm(0,isig,ifl(2,k),1)) !uW+(tAZ)d
         matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
c      k = 6                     !ducs
      elseif (k.eq.6) then 
      m1 = dotcc(jwm(0,-1,2),epsCCwp(0,isig,ifl(1,k),1))            !d(tAZ)W-u   
      m2 = dotcc(jwm(0,-1,1),epsCCwp(0,isig,ifl(2,k),1))            !dW-(tAZ)u
      z1 =qepsCCwp(1)*jj21m*(clr(ifl(1,k),2,isig)-clr(ifl(2,k),2,isig))
      mat(k,-1,isig,2) = clr(3,3,-1)**2 * (m1+m2+z1)
      if (lbox.and.(.not.qdamp)) then
         mv1 = dotcc(jvwm(0,-1,is1,2),epsCCwp(0,isig,ifl(1,k),1)) !d(tAZ)W-u   
         mv2 = dotcc(jvwm(0,-1,is1,1),epsCCwp(0,isig,ifl(2,k),1)) !dW-(tAZ)u
         matv(k,-1,isig,2) = clr(3,3,-1)**2 * (mv1+mv2)
      endif
      endif !k
      endif !nc/cc
c
c repeat the same for 2 bosons attached to 34 line

      if (nc_type) then  
c
      do isig = -1,1,2
         do kl=3,4
            qepswm(kl) = dotrc(qmp,epsNCwm(0,isig,kl,2))*zm2i(3)
            qepswp(kl) = dotrc(qpm,epsNCwp(0,isig,kl,2))*zm2i(3)
         enddo
         if (k.eq.1.or.k.eq.3) then 
c         do k = 1,3,2                  !uucc, ddcc
            m1 = dotcc(jwp(0,-1,3),epsNCwm(0,isig,ifl(1,k),2)) !cW+(tW-)c   
     1         - qepswm(ifl(1,k))*jj43p
            m2 = dotcc(jwm(0,-1,4),epsNCwp(0,isig,ifl(1,k),2)) !c(tW+)W-c
     1         + qepswp(ifl(1,k))*jj43m
            mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2)
            if (lbox.and.(.not.qdamp)) then
               mv1 = dotcc(jvwp(0,-1,is3,3),epsNCwm(0,isig,ifl(1,k),2)) !cW+(tW-)c   
               mv2 = dotcc(jvwm(0,-1,is3,4),epsNCwp(0,isig,ifl(1,k),2)) !c(tW+)W-c
               matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
c         enddo
         elseif (k.eq.2.or.k.eq.4) then   
c         do k = 2,4,2                  !uuss, ddss
            m1 = dotcc(jwp(0,-1,4),epsNCwm(0,isig,ifl(1,k),2))      !s(tW-)W+s
     1         + qepswm(ifl(1,k))*jj43p
            m2 = dotcc(jwm(0,-1,3),epsNCwp(0,isig,ifl(1,k),2))      !sW-(tW+)s
     1         - qepswp(ifl(1,k))*jj43m
            mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2)
            if (lbox.and.(.not.qdamp)) then
               mv1 = dotcc(jvwp(0,-1,is3,4),epsNCwm(0,isig,ifl(1,k),2)) !s(tW-)W+s
               mv2 = dotcc(jvwm(0,-1,is3,3),epsNCwp(0,isig,ifl(1,k),2)) !sW-(tW+)s
               matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
            endif
c         enddo
         endif !k   
      enddo !isig

      else !cc
      isig = -1
c      k = 5                     !udsc
      if (k.eq.5) then
         m1 = dotcc(jwm(0,-1,4),epsCCwp(0,isig,ifl(3,k),2)) !s(tAZ)W-c
         m2 = dotcc(jwm(0,-1,3),epsCCwp(0,isig,ifl(4,k),2)) !sW-(tAZ)c
         z1 =qepsCCwp(2)*jj43m*(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
         mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2+z1)
         if (lbox.and.(.not.qdamp)) then
            mv1 = dotcc(jvwm(0,-1,is3,4),epsCCwp(0,isig,ifl(3,k),2)) !s(tAZ)W-c
            mv2 = dotcc(jvwm(0,-1,is3,3),epsCCwp(0,isig,ifl(4,k),2)) !sW-(tAZ)c
            matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
         endif
c      k = 6                     !ducs
      elseif (k.eq.6) then    
         m1 = dotcc(jwp(0,-1,4),epsCCwm(0,isig,ifl(3,k),2)) !c(tAZ)W+s   
         m2 = dotcc(jwp(0,-1,3),epsCCwm(0,isig,ifl(4,k),2)) !cW+(tAZ)s
         z1 =qepsCCwm(2)*jj43p*(clr(ifl(3,k),2,isig)-clr(ifl(4,k),2,isig))
         mat(k,isig,-1,3) = clr(3,3,-1)**2 * (m1+m2+z1)
         if (lbox.and.(.not.qdamp)) then
            mv1 = dotcc(jvwp(0,-1,is3,4),epsCCwm(0,isig,ifl(3,k),2)) !c(tAZ)W+s   
            mv2 = dotcc(jvwp(0,-1,is3,3),epsCCwm(0,isig,ifl(4,k),2)) !cW+(tAZ)s
            matv(k,isig,-1,3) = clr(3,3,-1)**2 * (mv1+mv2)
         endif
      endif !k
      endif !nc/cc 
c
c ============================================================ 
c
c next come the A/Z-->WW currents attached to the quark lines. 
c The most effective structure is the contraction of two polarization 
c vectors with one fermion line. 
c First build these effective polarization vectors
c from the currents aww(mu) and azz(mu)
c
c NOTE: the aww and azz currents are conserved. Hence there is no 
c       need to consider q^mu * q^nu/m_Z^2  terms in the Z boson propagator
      if (lfs(js1)) then
         i = 1
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is1,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is1,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is1,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is1,i+1),fq(0,i+1))
         enddo
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is1,i),fq(0,i), ja(0,-1,is1,i)   )
         call curr6(1,psia(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is1,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is1,i),fq(0,i), jz(0,-1,is1,i)   )
         call curr6(1,psiz(1,-1,is1,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is1,i+1) )
         if (lbox.and.(.not.qdamp)) then
            do i = 1,2
               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                       .true., aww,3-i,   ja(0,-1,is1,i), 
     2                       jva(0,-1,is1,i) )
               call boxlinec(1,psi(1,-1,1),psi(1,-1,2),p(0,1),p(0,2),
     1                       .false., zww,3-i,  jz(0,-1,is1,i), 
     2                       jvz(0,-1,is1,i) )
            enddo
         endif
      endif
      if (lfs(js3)) then
         i = 3
         do isig = -1,1,2
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,aww,
     1                 psia(1,isig,is3,i),fq(0,i))
            call ket2c(psi(1,isig,i),.true.,p(0,i),isig,pww,zww,
     1                 psiz(1,isig,is3,i),fq(0,i))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,aww,
     1                 psia(1,isig,is3,i+1),fq(0,i+1))
            call bra2c(psi(1,isig,i+1),.true.,p(0,i+1),isig,pww,zww,
     1                 psiz(1,isig,is3,i+1),fq(0,i+1))
         enddo
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psia(1,-1,is3,i),fq(0,i), ja(0,-1,is3,i)   )
         call curr6(1,psia(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), ja(0,-1,is3,i+1) )
         call curr6(1,psi(1,-1,i+1),p(0,i+1),
     1              psiz(1,-1,is3,i),fq(0,i), jz(0,-1,is3,i)   )
         call curr6(1,psiz(1,-1,is3,i+1),fq(0,i+1),
     1              psi(1,-1,i),p(0,i), jz(0,-1,is3,i+1) )
         if (lbox.and.(.not.qdamp)) then
            do i = 3,4
               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                       .true., aww,5-i,   ja(0,-1,is3,i), 
     2                       jva(0,-1,is3,i) )
               call boxlinec(1,psi(1,-1,3),psi(1,-1,4),p(0,3),p(0,4),
     1                       .false., zww,5-i,   
     2                       jz(0,-1,is3,i), jvz(0,-1,is3,i) )
            enddo
         endif
      endif

      if (nc_type) then 
      do isig1 = -1,1,2
         do isig3 = -1,1,2
c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2)=j43 with momentum    p43        and
c    aww/zww(mu)         with momentum    pww

            ma(1) =      dotcc(ja(0,isig1,is1,2),jqq(0,isig3,2))
            ma(2) =      dotcc(ja(0,isig1,is1,1),jqq(0,isig3,2))
            mz(1) =      dotcc(jz(0,isig1,is1,2),jqq(0,isig3,2))
            mz(2) =      dotcc(jz(0,isig1,is1,1),jqq(0,isig3,2))
            if (lbox.and.(.not.qdamp)) then
               mva(1) = dotcc(jva(0,isig1,is1,2),jqq(0,isig3,2))
               mva(2) = dotcc(jva(0,isig1,is1,1),jqq(0,isig3,2))
               mvz(1) = dotcc(jvz(0,isig1,is1,2),jqq(0,isig3,2))
               mvz(2) = dotcc(jvz(0,isig1,is1,1),jqq(0,isig3,2))
            endif
            do kl = 1,4
               propt(isig1,isig3,kl,2) = 
     1            clr(ifl(1,kl),1,isig1)*clr(ifl(3,kl),1,isig3)*prop43(1) 
     2          + clr(ifl(1,kl),2,isig1)*clr(ifl(3,kl),2,isig3)*prop43(2)
            enddo !kl
   
               mat(k,isig1,isig3,4) = propt(isig1,isig3,k,2) *
     1            ( (ma(1)+ma(2))*clr(ifl(1,k),1,isig1) +
     2              (mz(1)+mz(2))*clr(ifl(1,k),2,isig1)   )
               if (lbox.and.(.not.qdamp)) matv(k,isig1,isig3,4) = propt(isig1,isig3,k,2)*
     1            ( (mva(1)+mva(2))*clr(ifl(1,k),1,isig1) +
     2              (mvz(1)+mvz(2))*clr(ifl(1,k),2,isig1)   )
           
c box correction to lower line: polarization vectors are 
c    jqq(mu,isig1,1)=j21 with momentum    p21        and
c    aww/zww(mu)         with momentum    pww

            ma(1) = dotcc(ja(0,isig3,is3,4),jqq(0,isig1,1))
            ma(2) = dotcc(ja(0,isig3,is3,3),jqq(0,isig1,1))
            mz(1) = dotcc(jz(0,isig3,is3,4),jqq(0,isig1,1))
            mz(2) = dotcc(jz(0,isig3,is3,3),jqq(0,isig1,1))
            if (lbox.and.(.not.qdamp)) then
               mva(1) = dotcc(jva(0,isig3,is3,4),jqq(0,isig1,1))
               mva(2) = dotcc(jva(0,isig3,is3,3),jqq(0,isig1,1))
               mvz(1) = dotcc(jvz(0,isig3,is3,4),jqq(0,isig1,1))
               mvz(2) = dotcc(jvz(0,isig3,is3,3),jqq(0,isig1,1))
            endif

            do kl = 1,4
               propt(isig1,isig3,kl,1) = 
     1            clr(ifl(1,kl),1,isig1)*clr(ifl(3,kl),1,isig3)*prop21(1)
     2          + clr(ifl(1,kl),2,isig1)*clr(ifl(3,kl),2,isig3)*prop21(2)
            enddo !kl   
               mat(k,isig1,isig3,5) = propt(isig1,isig3,k,1) *
     1            ( (ma(1)+ma(2))*clr(ifl(3,k),1,isig3) +
     2              (mz(1)+mz(2))*clr(ifl(3,k),2,isig3)   )
               if (lbox) matv(k,isig1,isig3,5) = propt(isig1,isig3,k,1)*
     1            ( (mva(1)+mva(2))*clr(ifl(3,k),1,isig3) +
     2              (mvz(1)+mvz(2))*clr(ifl(3,k),2,isig3)   )
           
         enddo !isig3
      enddo !isig1

      else !CC
      isig1 = -1
      isig3 = -1
c box correction to upper line: polarization vectors are 
c    jqq(mu,isig3,2)=j43 with momentum    p43        and
c    aww/zww(mu)         with momentum    pww
            ma(1) =      dotcc(ja(0,isig1,is1,2),jqq(0,isig3,2))
            ma(2) =      dotcc(ja(0,isig1,is1,1),jqq(0,isig3,2))
            mz(1) =      dotcc(jz(0,isig1,is1,2),jqq(0,isig3,2))
            mz(2) =      dotcc(jz(0,isig1,is1,1),jqq(0,isig3,2))
            if (lbox.and.(.not.qdamp)) then
               mva(1) = dotcc(jva(0,isig1,is1,2),jqq(0,isig3,2))
               mva(2) = dotcc(jva(0,isig1,is1,1),jqq(0,isig3,2))
               mvz(1) = dotcc(jvz(0,isig1,is1,2),jqq(0,isig3,2))
               mvz(2) = dotcc(jvz(0,isig1,is1,1),jqq(0,isig3,2))
            endif
           
c               do k = 5,6
                  propt(-1,-1,k,2) = clr(ifl(1,k),3,-1)**2*prop43(3)
                  mat(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( ma(1)*clr(ifl(2,k),1,-1) + ma(2)*clr(ifl(1,k),1,-1) 
     2       + mz(1)*clr(ifl(2,k),2,-1) + mz(2)*clr(ifl(1,k),2,-1) )
                  if (lbox.and.(.not.qdamp)) matv(k,-1,-1,4) = propt(-1,-1,k,2) *
     1       ( mva(1)*clr(ifl(2,k),1,-1) + mva(2)*clr(ifl(1,k),1,-1) 
     2       + mvz(1)*clr(ifl(2,k),2,-1) + mvz(2)*clr(ifl(1,k),2,-1) )
c               enddo

c box correction to lower line: polarization vectors are 
c    jqq(mu,isig1,1)=j21 with momentum    p21        and
c    aww/zww(mu)         with momentum    pww

            ma(1) = dotcc(ja(0,isig3,is3,4),jqq(0,isig1,1))
            ma(2) = dotcc(ja(0,isig3,is3,3),jqq(0,isig1,1))
            mz(1) = dotcc(jz(0,isig3,is3,4),jqq(0,isig1,1))
            mz(2) = dotcc(jz(0,isig3,is3,3),jqq(0,isig1,1))
            if (lbox.and.(.not.qdamp)) then
               mva(1) = dotcc(jva(0,isig3,is3,4),jqq(0,isig1,1))
               mva(2) = dotcc(jva(0,isig3,is3,3),jqq(0,isig1,1))
               mvz(1) = dotcc(jvz(0,isig3,is3,4),jqq(0,isig1,1))
               mvz(2) = dotcc(jvz(0,isig3,is3,3),jqq(0,isig1,1))
            endif

c               do k = 5,6
                  propt(-1,-1,k,1) = clr(ifl(1,k),3,-1)**2*prop21(3)
                  mat(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( ma(1)*clr(ifl(4,k),1,-1) + ma(2)*clr(ifl(3,k),1,-1) 
     2       + mz(1)*clr(ifl(4,k),2,-1) + mz(2)*clr(ifl(3,k),2,-1) )
                  if (lbox.and.(.not.qdamp)) matv(k,-1,-1,5) = propt(-1,-1,k,1) *
     1       ( mva(1)*clr(ifl(4,k),1,-1) + mva(2)*clr(ifl(3,k),1,-1) 
     2       + mvz(1)*clr(ifl(4,k),2,-1) + mvz(2)*clr(ifl(3,k),2,-1) )
c               enddo

      endif !cc            
c
c ==============================================
c
c  next do the box-box graphs with one W emitted from the upper and the 
c  other from the lower line. These are only possible for lefthanded quarks
c  on both lines, i.e. isig1 = -1 and isig3 = -1

c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,3)
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,1)
      mpm(1) = -dotcc(jwp(0,-1,1),jwm(0,-1,3))  !- zpm
      mmp(1) = -dotcc(jwm(0,-1,1),jwp(0,-1,3))  ! -zmp
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
      mpm(2) = -dotcc(jwp(0,-1,1),jwm(0,-1,4)) ! + zpm
      mmp(2) = -dotcc(jwm(0,-1,1),jwp(0,-1,4)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
      mpm(3) = -dotcc(jwp(0,-1,2),jwm(0,-1,3)) ! + zpm 
      mmp(3) = -dotcc(jwm(0,-1,2),jwp(0,-1,3)) ! + zmp
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
      mpm(4) = -dotcc(jwp(0,-1,2),jwm(0,-1,4)) ! - zpm
      mmp(4) = -dotcc(jwm(0,-1,2),jwp(0,-1,4)) ! - zmp

      if (lbox.and.(.not.qdamp)) then
c 1 = upper line box: eps1 = wp, eps2 = jwm(mu,-1,3)
c 2 = lower line box: eps1 = wm, eps2 = jwp(mu,-1,1)
         mvpm(1,1) = -dotcc(jvwp(0,-1,is1,1),jwm(0,-1,3))
         mvmp(1,1) = -dotcc(jvwm(0,-1,is1,1),jwp(0,-1,3))
         mvpm(1,2) = -dotcc(jwp(0,-1,1),jvwm(0,-1,is3,3))
         mvmp(1,2) = -dotcc(jwm(0,-1,1),jvwp(0,-1,is3,3))
         
c  upper line box: eps1 = wp, eps2 = jwm(mu,-1,4)
c  lower line box: eps1 = jwp(mu,-1,1), eps2 = wm
         mvpm(2,1) = -dotcc(jvwp(0,-1,is1,1),jwm(0,-1,4))
         mvmp(2,1) = -dotcc(jvwm(0,-1,is1,1),jwp(0,-1,4))
         mvpm(2,2) = -dotcc(jwp(0,-1,1),jvwm(0,-1,is3,4))
         mvmp(2,2) = -dotcc(jwm(0,-1,1),jvwp(0,-1,is3,4))
         
c  upper line box: eps1 = jwm(mu,-1,3), eps2 = wp
c  lower line box: eps1 = wm, eps2 = jwp(mu,-1,2)
         mvpm(3,1) = -dotcc(jvwp(0,-1,is1,2),jwm(0,-1,3)) 
         mvmp(3,1) = -dotcc(jvwm(0,-1,is1,2),jwp(0,-1,3))
         mvpm(3,2) = -dotcc(jwp(0,-1,2),jvwm(0,-1,is3,3)) 
         mvmp(3,2) = -dotcc(jwm(0,-1,2),jvwp(0,-1,is3,3))
         
c  upper line box: eps1 = jwm(mu,-1,4), eps2 = wp
c  lower line box: eps1 = jwp(mu,-1,2), eps2 = wm
         mvpm(4,1) = -dotcc(jvwp(0,-1,is1,2),jwm(0,-1,4))
         mvmp(4,1) = -dotcc(jvwm(0,-1,is1,2),jwp(0,-1,4))
         mvpm(4,2) = -dotcc(jwp(0,-1,2),jvwm(0,-1,is3,4))
         mvmp(4,2) = -dotcc(jwm(0,-1,2),jvwp(0,-1,is3,4))
         
      endif
c  for the q^mu*q^nu/M_V^2 terms in the gauge boson propagators we need
c      zpm = jj21p*jj43m*zm2i(3)
c      zmp = jj21m*jj43p*zm2i(3)
      if (nc_type) then
        zpm = jj21p*jj43m*zm2i(3)
        zmp = jj21m*jj43p*zm2i(3)
      else
        zpm = jj21p*jj43m*zm2i(2)
        zmp = jj21m*jj43p*zm2i(2)
      endif
c
      if (k.eq.1) then
      mat(k,-1,-1,6) = (mpm(2)+zpm)*prop_pm(3) + (mmp(3)+zmp)*prop_mp(3)
      elseif (k.eq.2) then
      mat(k,-1,-1,6) = (mpm(1)-zpm)*prop_pm(3) + (mmp(4)-zmp)*prop_mp(3) 
      elseif (k.eq.3) then
      mat(k,-1,-1,6) = (mpm(4)-zpm)*prop_pm(3) + (mmp(1)-zmp)*prop_mp(3) 
      elseif(k.eq.4) then
      mat(k,-1,-1,6) = (mpm(3)+zpm)*prop_pm(3) + (mmp(2)+zmp)*prop_mp(3)
      elseif (k.eq.5) then
      mat(k,-1,-1,6) = 
     1    mpm(1)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mpm(2)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mpm(3)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mpm(4)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3    + zpm*prop_pm(2)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4                     (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))
      
      elseif(k.eq.6) then
      mat(k,-1,-1,6) = 
     1    mmp(1)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mmp(2)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mmp(3)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mmp(4)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
     3    + zmp*prop_mp(2)*(clr(ifl(2,k),2,-1)-clr(ifl(1,k),2,-1))*
     4                     (clr(ifl(3,k),2,-1)-clr(ifl(4,k),2,-1))
      
      endif !k

      if (lbox.and.(.not.qdamp)) then
         do i = 1,2      ! 1 is for upper line, 2 for lower line QCD correction
            jj = 3 + 3*i ! stored in matv(...,6) and matv(...,9) respectively 
            if (k.eq.1) then
            matv(k,-1,-1,jj) = mvpm(2,i)*prop_pm(3)+mvmp(3,i)*prop_mp(3)
            elseif (k.eq.2) then
            matv(k,-1,-1,jj) = mvpm(1,i)*prop_pm(3)+mvmp(4,i)*prop_mp(3) 
            elseif (k.eq.3) then
            matv(k,-1,-1,jj) = mvpm(4,i)*prop_pm(3)+mvmp(1,i)*prop_mp(3) 
            elseif (k.eq.4) then
            matv(k,-1,-1,jj) = mvpm(3,i)*prop_pm(3)+mvmp(2,i)*prop_mp(3)
            elseif (k.eq.5) then
            matv(k,-1,-1,jj) = 
     1    mvpm(1,i)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvpm(2,i)*(prop_pm(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mvpm(3,i)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvpm(4,i)*(prop_pm(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_pm(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
            elseif (k.eq.6) then
            matv(k,-1,-1,jj) = 
     1    mvmp(1,i)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvmp(2,i)*(prop_mp(1)*clr(ifl(2,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(2,k),2,-1)*clr(ifl(3,k),2,-1) ) +
     1    mvmp(3,i)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(4,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(4,k),2,-1) ) +
     1    mvmp(4,i)*(prop_mp(1)*clr(ifl(1,k),1,-1)*clr(ifl(3,k),1,-1) + 
     2            prop_mp(2)*clr(ifl(1,k),2,-1)*clr(ifl(3,k),2,-1) ) 
            endif !k
         enddo
      endif
c
c -------------------------------------------------------------------
c
c and now, finally, the pentagon contributions, i.e. two W's emitted from the
c  same quark line
c
c prerequisites for virtual corrections:	
      if (lpent.and.icount1.eq.-1) then
c	 print*, 'pent. gauge check criterion chosen: mv5/mvc.gt.0.1d0'
         icount1 = icount1+1
      endif !icount1
      if (lpt) then
	 if (icountmax.eq.1) icountmax = 10	 
      endif !lpt   
	
      if (lpent .and. lfs(js1).and.(.not.qdamp)) then     ! need new pentagon graphs for 12 line
         musq = -p21(4)
         if (NLO.eq.1) then
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wp, wm, zero, zero, kl, -5, j5pm(0,kl,is1,1))
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wm, wp, zero, zero, kl, -5, j5mp(0,kl,is1,1))
            enddo
            icount1= icount1+1
         elseif (nlo.eq.-5) then
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wpp, wmp, xp, xm, kl, -5, j5pm(0,kl,is1,1) )
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wmp, wpp, xm, xp, kl, -5, j5mp(0,kl,is1,1) )
            enddo
            icount1= icount1+1
         else
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wpp, wmp, xp, xm, kl, 5, j5pm(0,kl,is1,1) )
               call J_virtual_tri_box_pent45(psi,p,1,2,sign, musq,
     1              wmp, wpp, xm, xp, kl, 5, j5mp(0,kl,is1,1) )
            enddo
         endif

         lgc(js1) = .false.
c         lfs(js1) = .false.
         if (lpt) then
	 
c  check numerical stability by gauge invariance
            if (.true.) then
               call gauge_check_5(j5pm(0,1,is1,1),p43,1,psi,
     %                            p,1,2,wpp,wmp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
               call gauge_check_5(j5mp(0,1,is1,1),p43,1,psi,
     %                            p,1,2,wmp,wpp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
            endif
	    if(bad_gauge) then
               icb1 = icb1+1
               lgc(js1) = .true.
            endif
c            if (mod(icount1,icountmax).eq.0 ) then
c               print*,' bad gauge check fraction ',real(icb1)/real(icount1), js1
c               print*,' bad gauge check params ',real(icb1),real(icount1)
c            endif
         endif   ! lpt
      endif

      if (lpent .and. lfs(js3).and.(.not.qdamp)) then     ! need new pentagon graphs for 34 line
         musq = -p43(4)
         if (NLO.eq.1) then
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wp, wm, zero, zero, kl, -5, j5pm(0,kl,is3,2))
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wm, wp, zero, zero, kl, -5, j5mp(0,kl,is3,2))
            enddo
            icount2= icount2+1
         elseif (nlo.eq.-5) then
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wpp, wmp, xp, xm, kl, -5, j5pm(0,kl,is3,2) )
                call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wmp, wpp, xm, xp, kl, -5, j5mp(0,kl,is3,2) )
            enddo
            icount2= icount2+1
         else
            do kl = 1,3
               call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wpp, wmp, xp, xm, kl, 5, j5pm(0,kl,is3,2) )
              call J_virtual_tri_box_pent45(psi,p,3,4,sign, musq,
     1              wmp, wpp, xm, xp, kl, 5, j5mp(0,kl,is3,2) )
            enddo
         endif

         lgc(js3) = .false.
c         lfs(js3) = .false.
         if (lpt) then
            if (.true.) then
               call gauge_check_5(j5pm(0,1,is3,2),p21,1,psi,
     %              p,3,4,wpp,wmp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
               call gauge_check_5(j5mp(0,1,is3,2),p21,1,psi,
     %                p,3,4,wmp,wpp,0d0)
	         if(bad_gauge_sin) bad_gauge = .true.
            endif
	    if(bad_gauge) then
               icb2 = icb2+1
               lgc(js3) = .true.
            endif
c            if (mod(icount2,icountmax).eq.0 ) then
c               print*,' bad gauge check fraction ',real(icb2)/real(icount2), js3
c            endif
         endif  !lpt
      endif

      fac = clr(3,3,-1)**2
c upper line:
      if (nc_type) then
      do isig3 = 1,-1,-2
         call ket2c(psi(1,-1,1),.true.,p(0,1),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,1),dummy)
         call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,2),dummy)
c  eps1=j43,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,2),wp,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wp,eps2=j43,eps3=wm
         m5(2,3) = -s1c(psiwm(1,2),jqq(0,isig3,2),.true.,-1,psiwp(1,1))
c  eps1=wp,eps2=wm,eps3=j43
         m5(3,3) = -s1c(bkjqq(1,-1,isig3,2),wm,.true.,-1,psiwp(1,1))
c  eps1=j43,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,2),wm,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wm,eps2=j43,eps3=wp
         m5(2,4) = -s1c(psiwp(1,2),jqq(0,isig3,2),.true.,-1,psiwm(1,1))
c  eps1=wm,eps2=wp,eps3=j43
         m5(3,4) = -s1c(bkjqq(1,-1,isig3,2),wp,.true.,-1,psiwm(1,1))
         if (lpent.and.(.not.qdamp)) then
            do j = 1,3
               mv5(j,3) = -dotcc(j5pm(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif

c         do k = 1,4
            kk = mod(k+2,4)
            if (k.eq.2) kk = 4 
            mat(k,-1,isig3,7) = fac*m5(2,ifl(1,k))*propt(-1,isig3,kk,2)+
     1         fac*( m5(1,ifl(2,k))+m5(3,ifl(1,k)) )*propt(-1,isig3,k,2)
            if (lpent.and.(.not.qdamp)) matv(k,-1,isig3,7) = 
     1         fac*mv5(2,ifl(1,k))*propt(-1,isig3,kk,2)+
     2         fac*( mv5(1,ifl(2,k))+mv5(3,ifl(1,k)) )*propt(-1,isig3,k,2)
c         enddo

      enddo !isig3
      
      else !cc
      isig3 = -1      
         call ket2c(psi(1,-1,1),.true.,p(0,1),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,1),dummy)
         call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     1        p43,jqq(0,isig3,2),bkjqq(1,-1,isig3,2),dummy)
c  eps1=j43,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,2),wp,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wp,eps2=wm,eps3=j43
         m5(3,3) = -s1c(bkjqq(1,-1,isig3,2),wm,.true.,-1,psiwp(1,1))
c  eps1=j43,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,2),wm,.true.,-1,bkjqq(1,-1,isig3,1))
c  eps1=wm,eps2=wp,eps3=j43
         m5(3,4) = -s1c(bkjqq(1,-1,isig3,2),wp,.true.,-1,psiwm(1,1))
         if (lpent.and.(.not.qdamp)) then
            do j = 1,3,2
               mv5(j,3) = -dotcc(j5pm(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is1,1),jqq(0,isig3,2))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif

c      do k = 5,6
         mat(k,-1,-1,7) = propt(-1,-1,k,2) * fac *
     1                    ( m5(1,ifl(2,k))+m5(3,ifl(1,k)) )
         if (lpent.and.(.not.qdamp)) matv(k,-1,-1,7) = propt(-1,-1,k,2) * fac *
     1                    ( mv5(1,ifl(2,k))+mv5(3,ifl(1,k)) )
         
c      enddo

      endif !nc/cc

c ======================================

c lower line:
      if (nc_type) then
      do isig1 = 1,-1,-2
         call ket2c(psi(1,-1,3),.true.,p(0,3),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,3),dummy)
         call bra2c(psi(1,-1,4),.true.,p(0,4),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,4),dummy)
c  eps1=j21,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,4),wp,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wp,eps2=j21,eps3=wm
         m5(2,3) = -s1c(psiwm(1,4),jqq(0,isig1,1),.true.,-1,psiwp(1,3))
c  eps1=wp,eps2=wm,eps3=j21
         m5(3,3) = -s1c(bkjqq(1,isig1,-1,4),wm,.true.,-1,psiwp(1,3))
c  eps1=j21,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,4),wm,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wm,eps2=j21,eps3=wp
         m5(2,4) = -s1c(psiwp(1,4),jqq(0,isig1,1),.true.,-1,psiwm(1,3))
c  eps1=wm,eps2=wp,eps3=j21
         m5(3,4) = -s1c(bkjqq(1,isig1,-1,4),wp,.true.,-1,psiwm(1,3))
         if (lpent.and.(.not.qdamp)) then
            do j = 1,3
               mv5(j,3) = -dotcc(j5pm(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif
c         do k = 1,4
            kk = k+1  - 2*mod(k+1,2)
            mat(k,isig1,-1,8) = fac*m5(2,ifl(3,k))*propt(isig1,-1,kk,1)+
     1         fac*( m5(1,ifl(4,k))+m5(3,ifl(3,k)) )*propt(isig1,-1,k,1)
            if (lpent.and.(.not.qdamp)) matv(k,isig1,-1,8) = 
     1         fac*mv5(2,ifl(3,k))*propt(isig1,-1,kk,1)+
     2         fac*( mv5(1,ifl(4,k))+mv5(3,ifl(3,k)) )*propt(isig1,-1,k,1)
c         enddo
      enddo !isig1

      else !cc
      
      isig1 = -1
         call ket2c(psi(1,-1,3),.true.,p(0,3),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,3),dummy)
         call bra2c(psi(1,-1,4),.true.,p(0,4),-1,
     1        p21,jqq(0,isig1,1),bkjqq(1,isig1,-1,4),dummy)
c  eps1=j21,eps2=wp,eps3=wm
         m5(1,3) = -s1c(psiwm(1,4),wp,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wp,eps2=wm,eps3=j21
         m5(3,3) = -s1c(bkjqq(1,isig1,-1,4),wm,.true.,-1,psiwp(1,3))
c  eps1=j21,eps2=wm,eps3=wp
         m5(1,4) = -s1c(psiwp(1,4),wm,.true.,-1,bkjqq(1,isig1,-1,3))
c  eps1=wm,eps2=wp,eps3=j21
         m5(3,4) = -s1c(bkjqq(1,isig1,-1,4),wp,.true.,-1,psiwm(1,3))
         if (lpent.and.(.not.qdamp)) then
            do j = 1,3,2
               mv5(j,3) = -dotcc(j5pm(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,3)
               mv5(j,4) = -dotcc(j5mp(0,j,is3,2),jqq(0,isig1,1))
c     1                    -cvirtc*m5(j,4)
            enddo
         endif
c      do k = 5,6
         mat(k,-1,-1,8) = propt(-1,-1,k,1) * fac *
     1                    ( m5(1,ifl(4,k))+m5(3,ifl(3,k)) )
         if (lpent.and.(.not.qdamp)) matv(k,-1,-1,8) = propt(-1,-1,k,1) * fac *
     1                    ( mv5(1,ifl(4,k))+mv5(3,ifl(3,k)) )
        
c      enddo
       
      endif !nc/cc   

c
c ------------------------------------------------------------
c
 999  continue
c
c sum the graphs, square them and map them onto uucc, uuss etc.

      if (lpt) then
         xgc1 = real(icount1+1)/real(icount1-icb1+1)
         xgc2 = real(icount2+1)/real(icount2-icb2+1)
      endif

	if (higgs_only) then
	  do kl = 1,6
	  do isig1 = -1,1,2
	  do isig3 = -1,1,2
	      mat(kl,isig1,isig3,2:9) = 0d0
	      matv(kl,isig1,isig3,2:9) = 0d0
	  enddo !isig3
          enddo !isig1
	  enddo !kl   
	endif

c      do k = 1,6
         res(k) = 0
         resv(k) = 0
         do isig1 = -1,1,2
            do isig3 = -1,1,2
               mm(k,isig1,isig3) = 0
               do i = 1,8
                  mm(k,isig1,isig3) = mm(k,isig1,isig3) + 
     1                                mat(k,isig1,isig3,i)
               enddo !i
               res(k) = res(k) + dreal(mm(k,isig1,isig3))**2
     &                         + dimag(mm(k,isig1,isig3))**2
               if (lnlo.and.(.not.qdamp)) then
                  mv12(k,isig1,isig3) = 
     1         +  matv(k,isig1,isig3,2) + matv(k,isig1,isig3,4) + 
     2            matv(k,isig1,isig3,6) + matv(k,isig1,isig3,7) 
                  mv34(k,isig1,isig3) =
     1          + matv(k,isig1,isig3,3) + matv(k,isig1,isig3,5) +
     2            matv(k,isig1,isig3,9) + matv(k,isig1,isig3,8)

                  if (lpt .and. isig1.eq.-1) then   
                  !Pentagon contributes for q_L only
                     if (lgc(js1)) then
                        mv12(k,isig1,isig3)=0
                     else
                        mv12(k,isig1,isig3)=mv12(k,isig1,isig3)*xgc1
                     endif !lgc
                  endif !lpt
                  if (lpt .and. isig3.eq.-1) then   
                  ! Pentagon contributes for q_L only
                     if (lgc(js3)) then
                        mv34(k,isig1,isig3)=0
                     else
                        mv34(k,isig1,isig3)=mv34(k,isig1,isig3)*xgc2
                     endif !lgc
                  endif !lpt

c  add Born type term and multiply by F_q = alphas*C_2/4pi;
c  the factor pi^2/3+9/2 for the born term is 
c  after adding the subtraction term
c  and the counter term for the renormalization of the pdfs
c
c  comply with POWHEG normalization:
c  vertex corrections in dreg:
c
c   V^nu = (4*Pi)^ep * Gamma(1+ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7)*B^nu
c
c         = (4*Pi)^ep / Gamma(1-ep) * CF * as/(4*Pi) * 
c           (-2/ep^2+(-2*ln(r)-3)/ep-ln(r)^2-3*ln(r)+Pi^2/3-7-Pi^2/3)*B^nu
c
c     The factor  (4*Pi)^ep/Gamma(1-ep) IS NOT RETURNED by this subroutine
c     and it's thought as factorized in front of the real counterterms too.
c
c     if (4*Pi)^ep / Gamma(1-ep) is collected in front then cvirt:
c     set parameter in dred (cvirtl = -8d0)
c
c     squared momentum of the weak boson connected with the upper line
          q2_up = p21(4)
c     squared momentum of the weak boson connected with the lower line
          q2_lo = p43(4)

          rup = st_muren2/(-q2_up)
          if (rup.lt.0d0.and..not.qdamp) then
              write(*,*) 'Error in setvirtual: q2_up should be < 0!!'
              write(*,*) 'q2_up = ', q2_up
              write(*,*) 'st_muren2 = ', st_muren2
              ! set to dummy value:
              	rup = 1d0
c              stop
          endif      
          rlo = st_muren2/(-q2_lo)
          if (rlo.lt.0d0.and..not.qdamp) then
              write(*,*) 'Error in setvirtual: q2_lo should be < 0!!'
              write(*,*) 'q2_lo = ', q2_lo
              write(*,*) 'st_muren2 = ', st_muren2
              ! set to dummy value:
              	rlo = 1d0
c              stop
          endif      
          lrup = log(rup)
          lrlo = log(rlo) 

          cvirtl = -8d0

                  if (nlo.gt.0) then
                     mv12(k,isig1,isig3) = als(1,1)*c2o4pi*
     1                ( mv12(k,isig1,isig3) + 
     2                  mm(k,isig1,isig3)*(-lrup**2-3*lrup+cvirtl) )
                     mv34(k,isig1,isig3) = als(2,1)*c2o4pi*
     1                ( mv34(k,isig1,isig3) + 
     2                  mm(k,isig1,isig3)*(-lrlo**2-3*lrlo+cvirtl) )
                  else
                     mv12(k,isig1,isig3) = 
     1                    als(1,1)*c2o4pi*mv12(k,isig1,isig3)
                     mv34(k,isig1,isig3) = 
     1                    als(2,1)*c2o4pi*mv34(k,isig1,isig3)
                  endif
                  resv(k) = resv(k) + 2*dreal(
     1                 mm(k,isig1,isig3)   *
     1              conjg( mv12(k,isig1,isig3)+mv34(k,isig1,isig3) )  )
               endif !lnlo
            enddo !isig3
         enddo  !isig1

         if (nlo.eq.0) then
            res(k) = res(k)*9d0
         elseif (nlo.gt.0) then
            if (qdamp) then
               !fakevirt:
                res(k) = res(k)*9d0   ! Born-type
                res(k) = res(k)*0.2d0 ! 'fakevirt' factor
            else !no damping  
c virt only (without Born):
               res(k) = (resv(k))*9d0 ! 9 is the color sum factor
            endif !qdamp
         else !nlo.lt.0
            if (qdamp) then
               !fakevirt:
                res(k) = res(k)*9d0   ! Born-type
                res(k) = res(k)*0.02d0 ! 'fakevirt' for boxes, pentagons
            else !no damping  
c virt only (without Born):
               res(k) = resv(k)*9d0 ! 9 is the color sum factor
            endif !qdamp
         endif !lnlo
c      enddo !k
           
cc eliminate processes with photon virtuality below cutoff
      if (qdamp) then 
          res(k) = res(k)*1d-20
      endif

      ans = res(k)

      return
      end


