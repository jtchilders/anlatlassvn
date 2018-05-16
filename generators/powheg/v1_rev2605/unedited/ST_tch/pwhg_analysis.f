c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include '../pwhg_book.h'
      include 'pwhg_math.h' 
      integer numplots
      real * 8 binsize(100)
      common/pwhghistcommon/binsize,numplots
      real *8 cut_etabj,cut_ptbj,cut_etal,cut_ptl,
     $     cut_etasj,cut_ptsj,cut_ptv
      common/ccuts/cut_etabj,cut_ptbj,cut_etal,cut_ptl,
     $     cut_etasj,cut_ptsj,cut_ptv
      real *8 jet_ktRadius,jet_ktptmin
      common/cjetsparam/jet_ktRadius,jet_ktptmin
    
      integer diag
      character * 10 cut
      integer icut

c     jet parameters
      jet_ktRadius = 0.7d0          
      jet_ktptmin  = 0d0

      jet_ktRadius = 1d0          
      jet_ktptmin  = 10d0


      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '                JET PARAMETERS                    '
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '   inclusive kt (fastkt implementation): '
      write(*,*) '      jet radius ',  jet_ktRadius
      write(*,*) '      jet ptmin ',   jet_ktptmin
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'


c     cut_1 set (as in Motylinski, apart from the cut on
c     the spectator jet, that is still as in 0702198)
      cut_etabj = 2d0     !2.5d0
      cut_ptbj  = 20d0    !50d0
      cut_etal = 2.5d0    !2.5d0
      cut_ptl  = 10d0     !20d0
      cut_ptv  = 20d0     !20d0

      cut_etasj= 2.5d0
      cut_ptsj = 20d0
c$$$      cut_etasj= 2.5d0
c$$$      cut_ptsj = 50d0

      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '                ANALYSIS CUTS                     '
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'
      write(*,*) '   CUTS: '
      write(*,*) '      eta b ',          cut_etabj 
      write(*,*) '      pt b ',           cut_ptbj  
      write(*,*) '      eta l ',          cut_etal 
      write(*,*) '      pt l ',           cut_ptl  
      write(*,*) '      pt vl ',          cut_ptv  
      write(*,*) '      eta spect-jet ',  cut_etasj         
      write(*,*) '      pt spect-jet ',   cut_ptsj 
      write(*,*) '**************************************************'
      write(*,*) '**************************************************'

      call pwhginihist
      
      icut=-1
      diag=0

 111  continue
      icut=icut+1

      if (icut.eq.0) then
         cut = ' nocut'
      elseif (icut.eq.1) then
         numplots=diag
         cut = ' cut_1'
      elseif (icut.eq.2) then
         goto 222
      else
         write(*,*) 'Error in init_hist, icut ',icut
         call exit(1)
      endif


c     top
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'t pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1.2d0
      call pwhgbookup(diag,'t pt zoom'//cut,'LOG',binsize(diag),
     $     0d0,60d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'t pt tail'//cut,'LOG',binsize(diag),
     $     0d0,1000d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'t y'//cut,'LOG',binsize(diag),-6d0,6d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'t eta'//cut,'LOG',binsize(diag),-6d0,6d0)

      diag=diag+1
      binsize(diag) = 0.8d0
      call pwhgbookup(diag,'t minv'//cut,'LOG',binsize(diag),
     $     175d0-20d0,175d0+20d0)


c     W
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'W pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1.2d0
      call pwhgbookup(diag,'W pt zoom'//cut,'LOG',binsize(diag),
     $     0d0,60d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'W pt tail'//cut,'LOG',binsize(diag),
     $     0d0,1000d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'W y'//cut,'LOG',binsize(diag),-6d0,6d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'W eta'//cut,'LOG',binsize(diag),-6d0,6d0)

      diag=diag+1
      binsize(diag) = 0.8d0
      call pwhgbookup(diag,'W minv'//cut,'LOG',binsize(diag),
     $     80d0-20d0,80d0+20d0)

c     hardest b
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'b pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1.2d0
      call pwhgbookup(diag,'b pt zoom'//cut,'LOG',binsize(diag),
     $     0d0,60d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'b pt tail'//cut,'LOG',binsize(diag),
     $     0d0,1000d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'b y'//cut,'LOG',binsize(diag),-6d0,6d0)

c     hardest bbar
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'bbar pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 1.2d0
      call pwhgbookup(diag,'bbar pt zoom'//cut,'LOG',binsize(diag),
     $     0d0,60d0)

      diag=diag+1
      binsize(diag) = 20d0
      call pwhgbookup(diag,'bbar pt tail'//cut,'LOG',binsize(diag),
     $     0d0,1000d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'bbar y'//cut,'LOG',binsize(diag),-6d0,6d0)

c     J1
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'j1 pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'j1 eta'//cut,'LOG',binsize(diag),-6d0,6d0)

c     J2
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'j2 pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'j2 eta'//cut,'LOG',binsize(diag),-6d0,6d0)

c     pt(t,j1)
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'tj1 pt'//cut,'LOG',binsize(diag),0d0,200d0)

c     dphi t-j1
      diag=diag+1
      binsize(diag) = pi/50d0
      call pwhgbookup(diag,'phitj1'//cut,'LIN',binsize(diag),0d0,pi)

c     hardest charged lepton
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'lep pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'lep eta'//cut,'LOG',binsize(diag),-6d0,6d0)

c     hardest neutrino
      diag=diag+1
      binsize(diag) = 4d0
      call pwhgbookup(diag,'v pt'//cut,'LOG',binsize(diag),0d0,200d0)

      diag=diag+1
      binsize(diag) = 0.24d0
      call pwhgbookup(diag,'v eta'//cut,'LOG',binsize(diag),-6d0,6d0)

c     t-channel spin correlation
c     cth between charged lepton and spectator jet
      diag=diag+1
      binsize(diag) = 0.04d0
      call pwhgbookup(diag,'cth l-sj'//cut,'LIN',binsize(diag),-1d0,1d0)


c     total cross section sanity check
      diag=diag+1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'total'//cut,'LOG',binsize(diag),0d0,3d0)

      goto 111


 222  end

     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h' 
      include 'pwhg_math.h' 
      include  'LesHouches.h'
c     other common blocks
      integer numplots
      real * 8 binsize(100)
      common/pwhghistcommon/binsize,numplots
      real *8 cut_etabj,cut_ptbj,cut_etal,cut_ptl,
     $     cut_etasj,cut_ptsj,cut_ptv
      common/ccuts/cut_etabj,cut_ptbj,cut_etal,cut_ptl,
     $     cut_etasj,cut_ptsj,cut_ptv
      real *8 jet_ktRadius,jet_ktptmin
      common/cjetsparam/jet_ktRadius,jet_ktptmin

      integer i,ihep,mu,ist_top,ist_w,ist,id,diag,j,jjet
      logical cuts,ini,condition_b,condition_bbar,condition,skipjet
      data ini/.true./
      save ini

      integer maxnum
      parameter (maxnum=10)
      real *8 p_top(0:3,maxnum),p_w(0:3,maxnum),
     $     p_l(0:3,maxnum),p_v(0:3,maxnum),
     $     p_b(0:3,maxnum),p_bbar(0:3,maxnum)
      integer nt,nw,nl,nv,nb,nbbar,nbjets,njets
      integer jb,jl,jv,bflavjet,specjet,ij1,ij2

      real *8 pj1(0:3),pj2(0:3),pspecj(0:3),pbj(0:3)

      real *8 pt_t,y_t,eta_t,minv_t,pt_w,y_w,eta_w,minv_w,pt_tj1
      real *8 pttemp,pt_b,y_b,pt_bbar,y_bbar,pt_j1,eta_j1,pt_j2,eta_j2
      real *8 tmp,dphi_tj1,pt_l,eta_l,pt_v,eta_v,cth,beta(3)
      real *8 eta_bjet,pt_bjet,pt_sj,eta_sj,pttemp_bjet,pttemp_spec
      real *8 ptemp(0:3),ptemp2(0:3)

      integer maxtrack,maxjet
      parameter (maxtrack=2048, maxjet=2048)
      real *8 pjet(4,maxtrack),ptrack(4,maxtrack)
      integer jpart,ntracks,jetvec(maxtrack),ihep_of_track(maxtrack),
     $     jetvec_of_bjet(maxtrack)


c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/

      real *8 azi

      if (ini) then
         write(*,*) '*****************************'
         if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
            write(*,*) '   NLO analysis not implemented        '
            write(*,*) 'No analysis will be run'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'

         ini=.false.
      endif

      if (whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         return
      elseif (WHCPRG.eq.'HERWIG'.or.WHCPRG.eq.'PYTHIA') then
         if(WHCPRG.eq.'HERWIG') then
            ist_top=155
            ist_w  =195 !: for undecayed events
            ist_w  =155 !: for decayed events
         elseif(WHCPRG.eq.'PYTHIA') then
            ist_top=3
            ist_w  =3
         endif

         nt=0
         nw=0
         nl=0
         nv=0
         nb=0
         nbbar=0

         do ihep=1,nhep
            ist=isthep(ihep)
            id=idhep(ihep)
            condition_b=.false.
            condition_bbar=.false.
            condition_b=(ist.eq.1).and.
     $           (((id.gt.-600).and.(id.lt.-500)).or.
     $           ((id.gt.5000).and.(id.lt.6000)))
            condition_bbar=(ist.eq.1).and.
     $           (((id.gt.500).and.(id.lt.600)).or.
     $           ((id.gt.-6000).and.(id.lt.-5000)))
            if(condition_b.and.condition_bbar) then
               write(*,*) 'both b and bbar true at the same time'
               call exit(1)
            endif

c     top
            if(ist.eq.ist_top.and.abs(id).eq.6) then
c     !: no ttype used
               nt=nt+1
               p_top(0,nt)=phep(4,ihep)
               do mu=1,3
                  p_top(mu,nt)=phep(mu,ihep)
               enddo
c     W
            elseif(ist.eq.ist_w.and.abs(id).eq.24) then
c     !: no ttype used
               nw=nw+1
               p_w(0,nw)=phep(4,ihep)
               do mu=1,3
                  p_w(mu,nw)=phep(mu,ihep)
               enddo
c     electron
            elseif(ist.eq.1.and.abs(id).eq.11) then
c     !: no ttype, or decaytype, used
               nl=nl+1
               p_l(0,nl)=phep(4,ihep)
               do mu=1,3
                  p_l(mu,nl)=phep(mu,ihep)
               enddo
c     electron neutrino
            elseif(ist.eq.1.and.abs(id).eq.12) then
c     !: no ttype, or decaytype, used
               nv=nv+1
               p_v(0,nv)=phep(4,ihep)
               do mu=1,3
                  p_v(mu,nv)=phep(mu,ihep)
               enddo
c     b hadrons
            elseif(condition_b) then
               nb=nb+1
               p_b(0,nb)=phep(4,ihep)
               do mu=1,3
                  p_b(mu,nb)=phep(mu,ihep)
               enddo
c     bbar hadrons
            elseif(condition_bbar) then
               nbbar=nbbar+1
               p_bbar(0,nbbar)=phep(4,ihep)
               do mu=1,3
                  p_bbar(mu,nbbar)=phep(mu,ihep)
               enddo
            endif
         enddo
         
c     set up arrays for jet finding
         ntracks=0
         njets=0
         do jpart=1,maxtrack
            do mu=1,4
               ptrack(mu,jpart)=0d0
            enddo
            jetvec(jpart)=0
         enddo      
         do jjet=1,maxjet
            do mu=1,4
               pjet(mu,jjet)=0d0
            enddo
         enddo

c     loop over final state particles, to find particles that
c     can be clustered
         do ihep=1,nhep
c     exclude leptons, gauge and Higgs bosons, but include gluons
            if ((isthep(ihep).eq.1).and.(
     $           (abs(idhep(ihep)).le.10).or.
     $           (abs(idhep(ihep)).ge.40).or.
     $           (abs(idhep(ihep)).eq.21))) then
c     copy momenta to construct jets 
               ntracks=ntracks+1
               ihep_of_track(ntracks)=ihep
               do mu=1,4
                  ptrack(mu,ntracks)=phep(mu,ihep)
               enddo
               if(ntracks.eq.maxtrack) then
                  write(*,*)
     $                 'hwanal: too many particles, increase maxtrack'
                  call exit(1)
               endif
            endif
         enddo
         if (ntracks.eq.0) then
            write(*,*) 'No tracks found, drop analysis of this event'
            goto 999
         endif

      else
         write(*,*) 'Invalid WHCPRG'
         call exit(1)
      endif

**********************************************************************
*     kt algorithm, fastjet implementation
      njets=0
      call fastjetktwhich(ptrack,ntracks,jet_ktptmin,jet_ktRadius,
     $     pjet,njets,jetvec) 
**********************************************************************

      skipjet=.false.
c     tag b-flavoured jets
      nbjets=0
      do i=1,ntracks
         ihep=ihep_of_track(i)
         condition_b=(
     $  ((abs(idhep(ihep)).gt.500).and.(abs(idhep(ihep)).lt.600)).or.
     $  ((abs(idhep(ihep)).gt.5000).and.(abs(idhep(ihep)).lt.6000)))
         if(condition_b) then
            nbjets=nbjets+1
            jetvec_of_bjet(nbjets)=jetvec(i)
         endif
      enddo

c     find the hardest non-b-flavoured jet (spectator jet)
c     and the hardest b-flavoured jet
      pttemp_spec=-1d0
      pttemp_bjet=-1d0
      specjet=-1
      bflavjet=-1
      do jjet=1,njets
c     loop on all jets
         condition=.true.
         do j=1,nbjets
c     loop on all b-jets: condition becomes false if current
c     jet (jjet) corresponds to one of the b-jets.
            condition=condition.and.(jjet.ne.jetvec_of_bjet(j))
         enddo
         if(condition) then
c     non-b-jet
            if((pjet(1,jjet)**2+pjet(2,jjet)**2).gt.pttemp_spec) then
               pttemp_spec=(pjet(1,jjet)**2+pjet(2,jjet)**2)
               specjet=jjet
            endif
         else
c     b-jet
            if((pjet(1,jjet)**2+pjet(2,jjet)**2).gt.pttemp_bjet) then
               pttemp_bjet=(pjet(1,jjet)**2+pjet(2,jjet)**2)
               bflavjet=jjet
            endif
         endif
      enddo

c     find the 2 hardest jets
c     We want to exclude the jet arising from the top decay.
c     In this analysis, we assume that this jet is the hardest
c     b-flavoured one, i.e. the one tagged above as bflavjet.
c     Notice that this is different from what we did in the
c     POWHEG-st code, where an explicit tjet search was present.

      if(njets.ge.1) then
c     Sort jets according to decreasing pt
c     Find the hardest jet
         pttemp=-1d0
         ij1=0
         tmp=0d0
         do jjet=1,njets
            if(jjet.ne.bflavjet) then
               tmp=sqrt(pjet(1,jjet)**2+pjet(2,jjet)**2)
               if(tmp.gt.pttemp) then
                  ij1=jjet
                  pttemp=tmp
               endif
            endif
         enddo
c     Now ij1 is the hardest jet
         if(ij1.ne.0) then
            pj1(0)=pjet(4,ij1)
            do mu=1,3
               pj1(mu)=pjet(mu,ij1)
            enddo
         endif

         if(njets.ge.2) then
c     Find the next-to-hardest jet
            pttemp=-1d0
            ij2=0
            tmp=0d0
            do jjet=1,njets
               if(jjet.ne.bflavjet) then
                  tmp=sqrt(pjet(1,jjet)**2+pjet(2,jjet)**2)
                  if((tmp.gt.pttemp).and.(jjet.ne.ij1))then
                     ij2=jjet
                     pttemp=tmp
                  endif
               endif
            enddo
c     Now ij2 is the next-to-hardest jet
            if(ij2.ne.0) then
               pj2(0)=pjet(4,ij2)
               do mu=1,3
                  pj2(mu)=pjet(mu,ij2)
               enddo
            endif
            if(ij1.eq.ij2) then
               write(*,*) 'Suspicious event: ij1=ij2'
               skipjet=.true.
            endif
         endif
      else
         skipjet=.true.
      endif



ccccccccccccccccccccccccccccccc
c     Observables
ccccccccccccccccccccccccccccccc      

c     t
      pt_t      = sqrt(p_top(1,1)**2+p_top(2,1)**2)
      call getrapidity(p_top(0,1),y_t)
      call get_pseudorap(p_top(0,1),eta_t)
      call getinvmass(p_top(0,1),minv_t)

c     w
      pt_w      = sqrt(p_w(1,1)**2+p_w(2,1)**2)
      call getrapidity(p_w(0,1),y_w)
      call get_pseudorap(p_w(0,1),eta_w)
      call getinvmass(p_w(0,1),minv_w)

c     hardest b
      if(nb.gt.0) then
         pttemp=-1d0
         jb=0
         do i=1,nb
            pt_b = sqrt(p_b(1,i)**2+p_b(2,i)**2)
            if(pt_b.gt.pttemp) then
               pttemp=pt_b
               jb=i
            endif
         enddo
         call getrapidity(p_b(0,jb),y_b)
      endif

c     hardest bbar
      if(nbbar.gt.0) then
         pttemp=-1d0
         jb=0
         do i=1,nbbar
            pt_bbar = sqrt(p_bbar(1,i)**2+p_bbar(2,i)**2)
            if(pt_bbar.gt.pttemp) then
               pttemp=pt_bbar
               jb=i
            endif
         enddo
         call getrapidity(p_bbar(0,jb),y_bbar)
      endif

      if(.not.skipjet) then
c     j1
      if(ij1.gt.0.and.nt.gt.0) then
         pt_j1     = sqrt(pj1(1)**2+pj1(2)**2)
         call get_pseudorap(pj1,eta_j1)
         dphi_tj1  = dabs(atan2(p_top(2,1),p_top(1,1)) - 
     $     atan2(pj1(2),pj1(1)))
         dphi_tj1=min(dphi_tj1,2d0*pi-dphi_tj1)

         dphi_tj1 = abs(azi(p_top(0,1))-azi(pj1))
         if (dphi_tj1.gt.pi) then
            dphi_tj1 = 2*pi-dphi_tj1
         endif



      endif

c     j2
      if(ij2.gt.0) then
         pt_j2     = sqrt(pj2(1)**2+pj2(2)**2)
         call get_pseudorap(pj2,eta_j2)
      endif
      endif

c     hardest lepton
      if(nl.gt.0) then
         pttemp=-1d0
         jl=0
         do i=1,nl
            pt_l = sqrt(p_l(1,i)**2+p_l(2,i)**2)
            if(pt_l.gt.pttemp) then
               pttemp=pt_l
               jl=i
            endif
         enddo
         call get_pseudorap(p_l(0,jl),eta_l)
      endif

c     hardest neutrino
      if(nv.gt.0) then
         pttemp=-1d0
         jv=0
         do i=1,nv
            pt_v = sqrt(p_v(1,i)**2+p_v(2,i)**2)
            if(pt_v.gt.pttemp) then
               pttemp=pt_v
               jv=i
            endif
         enddo
         call get_pseudorap(p_v(0,jv),eta_v)
      endif

c     spin-correlation (theta asymmetry)
      if(specjet.gt.0.and.nt.gt.0) then
         pspecj(0)=pjet(4,specjet)
         do mu=1,3
            beta(mu)=-p_top(mu,1)/p_top(0,1)
            pspecj(mu)=pjet(mu,specjet)
         enddo
         call boost_pow(beta,pspecj,ptemp)
         call boost_pow(beta,p_l(0,1),ptemp2)
         cth=(ptemp(1)*ptemp2(1)+ptemp(2)*ptemp2(2)+
     #ptemp(3)*ptemp2(3))
         cth=cth/sqrt(ptemp(1)**2+ptemp(2)**2+ptemp(3)**2)
         cth=cth/sqrt(ptemp2(1)**2+ptemp2(2)**2+ptemp2(3)**2)
      endif
      

cccccccccccccccccccccccccccccccccc
c     Histograms filling, no cuts
cccccccccccccccccccccccccccccccccc
      diag=0


c-----top
c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     y_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,y_t,dsig/binsize(diag))

c     eta_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,eta_t,dsig/binsize(diag))

c     minv_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,minv_t,dsig/binsize(diag))

c-----w
c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     y_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,y_w,dsig/binsize(diag))

c     eta_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,eta_w,dsig/binsize(diag))

c     minv_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,minv_w,dsig/binsize(diag))

c-----b
c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     y_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,y_b,dsig/binsize(diag))

c-----bbar
c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     y_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,y_bbar,dsig/binsize(diag))

c-----J1
c     pt_j1
      diag=diag+1
      if(ij1.gt.0.and..not.skipjet) 
     $     call pwhgfill(diag,pt_j1,dsig/binsize(diag))

c     eta_j1
      diag=diag+1
      if(ij1.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,eta_j1,dsig/binsize(diag))

c-----J2
c     pt_j2
      diag=diag+1
      if(ij2.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,pt_j2,dsig/binsize(diag))

c     eta_j2
      diag=diag+1
      if(ij2.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,eta_j2,dsig/binsize(diag))

c-----pt(t,j1)
      diag=diag+1
      if(ij1.gt.0.and.nt.gt.0.and..not.skipjet) then
         pt_tj1=sqrt((p_top(1,1)+pj1(1))**2+(p_top(2,1)+pj1(2))**2 )
         call pwhgfill(diag,pt_tj1,dsig/binsize(diag))
      endif

c-----phi(t-j1)
c     dphi_tj1
      diag=diag+1
      if(ij1.gt.0.and.nt.gt.0.and..not.skipjet) 
     $     call pwhgfill(diag,dphi_tj1,dsig/binsize(diag))

c-----lept
c     pt_l
      diag=diag+1
      if(nl.gt.0) call pwhgfill(diag,pt_l,dsig/binsize(diag))

c     eta_l
      diag=diag+1
      if(nl.gt.0) call pwhgfill(diag,eta_l,dsig/binsize(diag))

c-----neutr
c     pt_l
      diag=diag+1
      if(nv.gt.0) call pwhgfill(diag,pt_v,dsig/binsize(diag))

c     eta_l
      diag=diag+1
      if(nv.gt.0) call pwhgfill(diag,eta_v,dsig/binsize(diag))

c-----spin corr
      diag=diag+1
      if(specjet.gt.0.and.nt.gt.0) 
     $     call pwhgfill(diag,cth,dsig/binsize(diag))

c-----total
      diag=diag+1
      call pwhgfill(diag,1.5d0,dsig/binsize(diag))


ccccccccccccccccccccccccccccc
c     Histogram filling, cuts
ccccccccccccccccccccccccccccc

c     Cuts, set 1
      cuts=.false.
      if(bflavjet.gt.0.and.specjet.gt.0) then
         pbj(0)=pjet(4,bflavjet)
         do mu=1,3
            pbj(mu)=pjet(mu,bflavjet)
         enddo
         pt_bjet=sqrt(pbj(1)**2+pbj(2)**2)
         call get_pseudorap(pbj,eta_bjet)

cccccccccccccccccccccccccccccccccccc
c$$$         pt_bjet=pt_b
c$$$         eta_bjet=y_b
cccccccccccccccccccccccccccccccccccc

         pt_sj=sqrt(pspecj(1)**2+pspecj(2)**2)
         call get_pseudorap(pspecj,eta_sj)

         cuts=(
     $        (dabs(eta_bjet) .lt.cut_etabj).and.
     $        (pt_bjet        .gt.cut_ptbj).and.
     $        (dabs(eta_l) .lt.cut_etal).and.
     $        (pt_l        .gt.cut_ptl).and.
     $        (pt_v        .gt.cut_ptv).and.
     $        (dabs(eta_sj).lt.cut_etasj).and.
     $        (pt_sj       .gt.cut_ptsj))
c$$$     $        (dabs(eta_sj).gt.cut_etasj).and.
c$$$     $        (pt_sj       .gt.cut_ptsj))


      endif

      if(cuts) then
c     reply same structure as above, without definition of variables
c-----top
c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     pt_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,pt_t,dsig/binsize(diag))

c     y_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,y_t,dsig/binsize(diag))

c     eta_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,eta_t,dsig/binsize(diag))

c     minv_top
      diag=diag+1
      if(nt.gt.0) call pwhgfill(diag,minv_t,dsig/binsize(diag))

c-----w
c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     pt_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,pt_w,dsig/binsize(diag))

c     y_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,y_w,dsig/binsize(diag))

c     eta_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,eta_w,dsig/binsize(diag))

c     minv_w
      diag=diag+1
      if(nw.gt.0) call pwhgfill(diag,minv_w,dsig/binsize(diag))

c-----b
c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     pt_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,pt_b,dsig/binsize(diag))

c     y_b
      diag=diag+1
      if(nb.gt.0) call pwhgfill(diag,y_b,dsig/binsize(diag))

c-----bbar
c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     pt_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,pt_bbar,dsig/binsize(diag))

c     y_bbar
      diag=diag+1
      if(nbbar.gt.0) call pwhgfill(diag,y_bbar,dsig/binsize(diag))

c-----J1
c     pt_j1
      diag=diag+1
      if(ij1.gt.0.and..not.skipjet) 
     $     call pwhgfill(diag,pt_j1,dsig/binsize(diag))

c     eta_j1
      diag=diag+1
      if(ij1.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,eta_j1,dsig/binsize(diag))

c-----J2
c     pt_j2
      diag=diag+1
      if(ij2.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,pt_j2,dsig/binsize(diag))

c     eta_j2
      diag=diag+1
      if(ij2.gt.0.and..not.skipjet)
     $     call pwhgfill(diag,eta_j2,dsig/binsize(diag))

c-----pt(t,j1)
      diag=diag+1
      if(ij1.gt.0.and.nt.gt.0.and..not.skipjet) then
         pt_tj1=sqrt((p_top(1,1)+pj1(1))**2+(p_top(2,1)+pj1(2))**2 )
         call pwhgfill(diag,pt_tj1,dsig/binsize(diag))
      endif

c-----phi(t-j1)
c     dphi_tj1
      diag=diag+1
      if(ij1.gt.0.and.nt.gt.0.and..not.skipjet) 
     $     call pwhgfill(diag,dphi_tj1,dsig/binsize(diag))

c-----lept
c     pt_l
      diag=diag+1
      if(nl.gt.0) call pwhgfill(diag,pt_l,dsig/binsize(diag))

c     eta_l
      diag=diag+1
      if(nl.gt.0) call pwhgfill(diag,eta_l,dsig/binsize(diag))

c-----neutr
c     pt_l
      diag=diag+1
      if(nv.gt.0) call pwhgfill(diag,pt_v,dsig/binsize(diag))

c     eta_l
      diag=diag+1
      if(nv.gt.0) call pwhgfill(diag,eta_v,dsig/binsize(diag))

c-----spin corr
      diag=diag+1
      if(specjet.gt.0.and.nt.gt.0) 
     $     call pwhgfill(diag,cth,dsig/binsize(diag))

c-----total
      diag=diag+1
      call pwhgfill(diag,1.5d0,dsig/binsize(diag))


      endif

 999  end
      

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(0:3),y
      real *8 tiny
      parameter (tiny=1.d-8)
c     !: protect for small p(0)-p(3) values
      if(dabs(p(0)-p(3)).lt.tiny) then
         y=sign(1.d0,p(3))*1.d8
      else
         y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      endif
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(dabs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end

      subroutine get_pseudorap(p,eta)
      implicit none
      real*8 p(0:3),eta,pt,th
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

      function azi(p)
      implicit none
      include 'pwhg_math.h'  
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end



