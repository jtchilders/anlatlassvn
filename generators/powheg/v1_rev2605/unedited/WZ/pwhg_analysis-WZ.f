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

      call inihists

c     total cross section sanity check
      call bookupeqbins('total',1d0,0d0,1d0)

      call  bookupeqbins('httot',20d0,0d0,500d0)

      call  bookupeqbins('Z_massa',5d0,0d0,200d0)
      call  bookupeqbins('Z_massb',2d0,70d0,110d0)
      call  bookupeqbins('Z_pt',5d0,0d0,200d0)
      call  bookupeqbins('Z_eta',0.2d0,-4d0,4d0)
      call  bookupeqbins('Z_y',0.2d0,-4d0,4d0)

      call  bookupeqbins('lW_eta',0.2d0,-4d0,4d0)
      call  bookupeqbins('lW_pt',5d0,0d0,150d0)
      call  bookupeqbins('lZm_eta',0.2d0,-4d0,4d0)
      call  bookupeqbins('lZm_pt',5d0,0d0,150d0)
      call  bookupeqbins('lZp_eta',0.2d0,-4d0,4d0)
      call  bookupeqbins('lZp_pt',5d0,0d0,150d0)

      call  bookupeqbins('miss_pt',5d0,0d0,200d0)

      call  bookupeqbins('3l_mass',10d0,0d0,500d0)
      call  bookupeqbins('Z_mta',5d0,0d0,150d0)
      call  bookupeqbins('Z_mtb',1d0,80d0,120d0)
      call  bookupeqbins('W_mta',5d0,0d0,150d0)
      call  bookupeqbins('W_mtb',1d0,60d0,100d0)

      call bookupeqbins('ch_y_asyml',0.2d0,-4d0,4d0)
      call bookupeqbins('ch_pt_asyml',5d0,0d0,500d0)

C     -- variables involving the jet 
      call bookupeqbins('j_ptinp', 5d0,0d0,150d0)
      call bookupeqbins('j_etainp',0.2d0,-4d0,4d0)
      call bookupeqbins('jZphiinp',0.314d0,0d0,3.14d0)
      call bookupeqbins('jZetainp',0.2d0,-6d0,6d0)
      call bookupeqbins('jlRinp',0.1d0,0d0,4d0)

      call bookupeqbins('j_ptA', 5d0,40d0,200d0)
      call bookupeqbins('j_etaA',0.2d0,-4d0,4d0)
      call bookupeqbins('jZphiA',0.314d0,0d0,3.14d0)
      call bookupeqbins('jZetaA',0.2d0,-6d0,6d0)
      call bookupeqbins('jlRA',0.1d0,0d0,4d0)

      call bookupeqbins('j_ptB', 5d0,10d0,400d0)
      call bookupeqbins('j_etaB',0.2d0,-4d0,4d0)
      call bookupeqbins('jZphiB',0.314d0,0d0,3.14d0)
      call bookupeqbins('jZetaB',0.2d0,-6d0,6d0)
      call bookupeqbins('jlRB',0.1d0,0d0,4d0)

      call bookupeqbins('lll_pt',10d0,0d0,500d0)
      call bookupeqbins('alllept_pt',10d0,0d0,500d0)
      call bookupeqbins('W_pt', 10d0,0d0,500d0)
      call bookupeqbins('W_phi', 0.314d0,0d0,3.14d0)
      call bookupeqbins('Z_phi', 0.314d0,0d0,3.14d0)
      
      call bookupeqbins('pt_WZ', 10d0,0d0,500d0)
      call bookupeqbins('m_WZ', 10d0,0d0,500d0)

      end 
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      include 'cvecbos.h'
      include 'nwz.f'
      include 'constants.f'
      logical passcuts, passjetcuts
      integer ihep,mu
      logical ini
      data ini/.true./
      save ini
c     binsize
c     we need to tell to this analysis file which program is running it
      integer   maxjet,mjets
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer vdecaytempW1,vdecaytempW2,vdecaytempZ1,vdecaytempZ2
      integer nlmZ, nlpZ,nlW,nnuW,ilmZ(100),ilpZ(100),ilW(100),inuW(100)
      real * 8 httot,y,eta,mass,maxptl,dphi,deta,dr,dy,pt
c      integer i3,i4,i5,i6,j,nu,ijc,iZmd,bestiZmd
      real * 8 powheginput
      external powheginput
c-----TM added
      double precision plZm(4),plZp(4),plW(4),pnu(4)
      double precision pt_lll,pt_alllept,pt_W
      double precision eta1,eta2,eta3,eta4,pt1,pt2,pt3,pt4
      double precision massZZ,massWZ
      double precision maxpt,minpt,maxeta
      double precision dphi_W, dphi_Z
      real * 8 pt_3456,mass_3456,eta_3456,y_3456
      integer j
      double precision mtrans,fbZm,fbZp,fbW
      double precision zmass
      parameter (zmass=91.1876d0)      
      double precision dr1,dr2,dr3,drmin

      if(dsig.eq.0) return

      if (ini) then
         vdecaymodeW = powheginput("vdecaymodeW")
         vdecaymodeZ = powheginput("vdecaymodeZ")
         if (vdecaymodeW.lt.0) then
            nwz = +1
            write(*,*)'POWHEG: W+ Z production and decay'
         elseif (vdecaymodeW.gt.0) then
            nwz = -1
            write(*,*)'POWHEG: W- Z production and decay'
         endif

         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis WZ'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'

         if ((abs(vdecaymodeZ).eq.11.or.abs(vdecaymodeZ).eq.13
     $        .or.abs(vdecaymodeZ).eq.15) .and.
     $        (abs(vdecaymodeW).eq.11 .or. abs(vdecaymodeW).eq. 13 .or.
     $        abs(vdecaymodeW).eq. 15)) then
            continue
         else
            write(*,*) '**************************************'
            write(*,*) ' template analysis works only for e, mu and tau'
            write(*,*) '                 STOP     '
            write(*,*) '**************************************'
            call exit(1)
         endif
         ini=.false.
      endif

      vdecaytempW1 = vdecaymodeW
      vdecaytempZ1 = abs(vdecaymodeZ)
      
      vdecaytempZ2 = -vdecaytempZ1  ! other l from Z
      vdecaytempW2 =-sign(1,vdecaytempW1)*(abs(vdecaytempW1)+1)   ! nu from W


c     find WZ decay products
         nlW=0                 ! l(+ or -) from W
         nnuW=0                ! nu(bar) from W
         nlmZ=0                ! l+ from Z
         nlpZ=0                ! l- from Z            
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.vdecaytempZ1) then
                  nlmZ = nlmZ+1
                  ilmZ(nlmZ) = ihep
               elseif(idhep(ihep).eq.vdecaytempZ2) then
                  nlpZ = nlpZ+1
                  ilpZ(nlpZ) = ihep
               elseif(idhep(ihep).eq.vdecaytempW1) then
                  nlW = nlW+1
                  ilW(nlW) = ihep
               elseif(idhep(ihep).eq.vdecaytempW2) then
                  nnuW = nnuW+1
                  inuW(nnuW) = ihep
               endif
            endif
         enddo

      if (nlmZ .ge. 100 .or. nlpZ .ge. 100 .or. nlW .ge. 100 .or.
     .     nnuW .ge. 100) then
         write(*,*) 'crazy event, too many leptons', nlmZ, nlpz,nlw, 
     .        nnuW
         return
      endif


      if (abs(vdecaymodeW).eq.abs(vdecaymodeZ)) then
c - all the 'W' leptons are in nlmZ for W- production
c - or in nlpZ for W+ production
         if (nwz.eq.1) then
            if(nlmZ.lt.1 .or. nlpZ.lt.2) then
               write(*,*)'crazy event, missing leptons',
     .              nlmz,nlpz,nlw,nnuw
            endif
         elseif (nwz.eq.-1) then
            if(nlmZ.lt.2 .or. nlpZ.lt.1) then
               write(*,*)'crazy event, missing leptons',
     .              nlmz,nlpz,nlw,nnuw
            endif
         endif
      else
         if (nlmZ.lt. 1 .or. nlpZ.lt.1 .or. nlW.lt.1 .or. nnuW.lt.1)then
            write(*,*) vdecaymodeW, vdecaymodeZ
            write(*,*) vdecaytempW1, vdecaytempW2
            write(*,*) vdecaytempZ1, vdecaytempZ2
            write(*,*) 'crazy event, missing leptons', nlmZ, nlpz,nlw, 
     .           nnuW
            write(*,*) idhep(:nhep) 
            stop 
            return
         endif
      endif

      
c-----sort the leptons by pt. If like flavour leptons
c-----are present then have to take second highest pt
c-----lepton from the Z pool of leptons, with the correct
c-----charge depending on whether we are doing W+ or W-
      if (abs(vdecaymodeW).eq.abs(vdecaymodeZ)) then
         call sortbypt(nlmZ,ilmZ(1:nlmZ))
         call sortbypt(nlpZ,ilpZ(1:nlpZ))
         call sortbypt(nnuW,inuW(1:nnuW))

         plZm(1:4)=phep(1:4,ilmZ(1))
         plZp(1:4)=phep(1:4,ilpZ(1))
         pnu(1:4) =phep(1:4,inuW(1))         
         if (nwz.eq.1) then
            plW(1:4) = phep(1:4,ilpZ(2))
         elseif (nwz.eq.-1) then
            plW(1:4) = phep(1:4,ilmZ(2))
         endif
      else
         call sortbypt(nlmZ,ilmZ(1:nlmZ))
         call sortbypt(nlpZ,ilpZ(1:nlpZ))
         call sortbypt(nlW,ilW(1:nlW))
         call sortbypt(nnuW,inuW(1:nnuW))

         plZm(1:4)=phep(1:4,ilmZ(1))
         plZp(1:4)=phep(1:4,ilpZ(1))
         plW(1:4)=phep(1:4,ilW(1))
         pnu(1:4) =phep(1:4,inuW(1))         
      endif
      rr=0.6d0
      call buildjets(1,mjets,rr,ktj,etaj,rapj,phij,pj)


c-----this is sorting out the order of the leptons of like flav,
c----- and opposite charge
      call getyetaptmass(plZm+plZp,y,eta,pt,massZZ)
      if (abs(vdecaymodeW).eq.abs(vdecaymodeZ)) then         
         if (nwz.eq.1) then
            call getyetaptmass(plZm+plW,y,eta,pt,massWZ)
         elseif (nwz.eq.-1) then
            call getyetaptmass(plZp+plW,y,eta,pt,massWZ)
         endif

         if (abs(massWZ-zmass).lt.abs(massZZ-zmass)) then
            if (nwz.eq.1) then
               plZp(1:4)=phep(1:4,ilpz(2))
               plW(1:4) =phep(1:4,ilpz(1))
            elseif (nwz.eq.-1) then
               plZm(1:4)=phep(1:4,ilmz(2))
               plW(1:4) =phep(1:4,ilmz(1))
            endif
         endif
      endif
c-----here we have plZp, plZm, plW, pnu ordered


c-----make cuts on the lepton pt,eta  and missing pt
      call getyetaptmass(plW,y,eta1,pt1,mass)
      call getyetaptmass(plZm,y,eta2,pt2,mass)
      call getyetaptmass(plZp,y,eta3,pt3,mass)
      call getyetaptmass(pnu,y,eta4,pt4,mass)

      maxpt = max(pt1,pt2,pt3)
      minpt = min(pt1,pt2,pt3)
      if (maxpt.lt.20d0) return
      if (minpt.lt.10d0) return
      maxeta = max(abs(eta1),abs(eta2),abs(eta3))
      if (maxeta.gt.2.5d0) return
      if (pt4.lt.20d0) return
c--------------------------------------------------


c-----start doing some plots


      httot=pt1+pt2+pt3+pt4
      do j=1,mjets
         httot=httot+ktj(j)
      enddo


      call filld('total',0d0,dsig)
      call filld('httot',httot,dsig)

c-----Zplots
      call getyetaptmass(plZm+plZp,y,eta,pt,mass)
      call filld('Z_massa', mass,dsig)
      call filld('Z_massb', mass,dsig)
      call filld('Z_pt', pt,dsig)
      call filld('Z_eta', eta,dsig)
      call filld('Z_y', y,dsig)


c-----individual leptons
      call filld('lW_eta',eta1,dsig)
      call filld('lW_pt',pt1,dsig)
      call filld('lZm_eta',eta2,dsig)
      call filld('lZm_pt',pt2,dsig)
      call filld('lZp_eta',eta3,dsig)
      call filld('lZp_pt',pt3,dsig)

c-----missing pt      
      call filld('miss_pt',pt4,dsig)

c-----massy plots
      call getyetaptmass(plZm+plZp+plW,y,eta,pt_lll,mass)
      call filld('3l_mass',mass,dsig)
      call getdydetadphidr(plZm,plZp,dy,deta,dphi_Z,dr)
      mtrans=dsqrt(2d0*pt2*pt3*(1d0-dcos(dphi)))
      call filld('Z_mta',mtrans,dsig)
      call filld('Z_mtb',mtrans,dsig)
      call getdydetadphidr(plW,pnu,dy,deta,dphi_W,dr)
      mtrans=dsqrt(2d0*pt1*pt4*(1d0-dcos(dphi_W)))
      call filld('W_mta',mtrans,dsig)
      call filld('W_mtb',mtrans,dsig)
      call getyetaptmass(plZm+plZp+plW+pnu,y_3456,eta_3456,
     &     pt_3456,mass_3456)
      call filld('m_WZ', mass_3456,dsig)
      call filld('pt_WZ', pt_3456,dsig)

c -- plots for ATGCs
      pt_alllept = pt4 + pt_lll
      call getyetaptmass(plW+pnu,y,eta,pt_W,mass)
      call filld('lll_pt', pt_lll, dsig)
      call filld('alllept_pt', pt_alllept,dsig)
      call filld('W_pt', pt_W, dsig)
      call filld('W_phi', dphi_W,dsig)
      call filld('Z_phi', dphi_Z,dsig)


c     y charge asymmetry 
      call filld('ch_y_asyml', eta3,  dsig/2d0)
      call filld('ch_y_asyml', eta2, -dsig/2d0)

c     pt charge asymmetry 
      call filld('ch_pt_asyml', pt3, dsig/2d0)
      call filld('ch_pt_asyml', pt2,-dsig/2d0)


c     --- jet plots 
      if(mjets.ge.1) then
      call getdydetadphidr(plZm,pj(:,1),dy,deta,dphi,dr1)
      call getdydetadphidr(plZp,pj(:,1),dy,deta,dphi,dr2)
      call getdydetadphidr(plW,pj(:,1),dy,deta,dphi,dr3)
      drmin = min(dr1,dr2,dr3)
      call getyetaptmass(pj(:,1),y,eta,pt,mass)
      call getdydetadphidr(plZm+plZp,pj(:,1),dy,deta,dphi,dr)
        
      call filld('j_ptinp',pt,dsig)
C      call filld('j_etainp',eta,dsig)
C      call filld('jZphiinp',dphi,dsig)
C      call filld('jZetainp',deta,dsig)
C      call filld('jlRinp',drmin,dsig)

      if (pt.gt.20) then
      call filld('j_ptA',pt,dsig)
      call filld('j_etaA',eta,dsig)
      call filld('jZphiA',dphi,dsig)
      call filld('jZetaA',deta,dsig)
      call filld('jlRA',drmin,dsig)
      endif

      if (pt.gt.100) then
      call filld('j_ptB',pt,dsig)
      call filld('j_etaB',eta,dsig)
      call filld('jZphiB',dphi,dsig)
      call filld('jZetaB',deta,dsig)
      call filld('jlRB',drmin,dsig)
      endif

      endif

      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(*),y,eta,pt,mass,pv
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      mass=sqrt(abs(p(4)**2-pv**2))
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(0:3),y
      y=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getinvmass(p,m)
      implicit none
      real * 8 p(0:3),m
      m=sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
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



      subroutine buildjets(iflag,mjets,rr,kt,eta,rap,phi,pjet)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer iflag,mjets
      real * 8  rr,kt(*),eta(*),rap(*),
     1     phi(*),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,jb
      real * 8 r,palg,ptmin,pp,tmp
      logical islept
      external islept
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         do j=1,nhep
            if (isthep(j).eq.1.and..not.islept(idhep(j))) then
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
      else
         do j=1,nup
            if (istup(j).eq.1.and..not.islept(idup(j))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C - Inclusive jet pT and Y spectra are to be compared to CDF data:    - C    
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
C     f = 0.75  overlapping fraction
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
      ptmin=1d0
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets
      if(njets.eq.0) return
c check consistency
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
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      end



      subroutine sortbypt(n,iarr)
      implicit none
      integer n,iarr(n)
      include 'hepevt.h'
      integer j,k
      real * 8 tmp,pt(nmxhep)
      logical touched
      do j=1,n
         pt(j)=sqrt(phep(1,iarr(j))**2+phep(2,iarr(j))**2)
      enddo
c bubble sort
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


      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end


      subroutine boostx(p_in,pt,ptt,p_out)
      implicit none
c--- Boost input vector p_in to output vector p_out using the same
c--- transformation as required to boost massive vector pt to ptt
      double precision p_in(4),pt(4),ptt(4),p_out(4),
     . p_tmp(4),beta(3),mass,gam,bdotp
      integer j
    
      mass=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2  
      if (mass .lt. 0d0) then
        write(6,*) 'mass**2 .lt. 0 in boostx.f, mass**2=',mass,pt
        stop
      endif
      mass=dsqrt(mass)

c--- boost to the rest frame of pt
      gam=pt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=-pt(j)/pt(4)
        bdotp=bdotp+beta(j)*p_in(j)
      enddo
      p_tmp(4)=gam*(p_in(4)+bdotp)
      do j=1,3
        p_tmp(j)=p_in(j)+gam*beta(j)/(1d0+gam)*(p_in(4)+p_tmp(4))
      enddo     

c--- boost from rest frame of pt to frame in which pt is identical
c--- with ptt, thus completing the transformation          
      gam=ptt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=+ptt(j)/ptt(4)
        bdotp=bdotp+beta(j)*p_tmp(j)
      enddo
      p_out(4)=gam*(p_tmp(4)+bdotp)
      do j=1,3
        p_out(j)=p_tmp(j)+gam*beta(j)/(1d0+gam)*(p_out(4)+p_tmp(4))
      enddo

      return
      end
      

      
