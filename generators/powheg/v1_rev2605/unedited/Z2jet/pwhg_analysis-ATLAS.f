C----------------------------------------------------------------------------
C     This analysis file is largely inspired by the ATLAS paper 1111.2690  
      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      character * 1 cnum(0:9)
      data cnum/'0','1','2','3','4','5','6','7','8','9'/
      character * 1 pr
      character * 6 suffix 
      common/pwhgprocess/pr
      real * 8 dy,dpt,dr,dptzoom
      integer j,i,jetcut

      call inihists

      pr='Z'
      dy=0.5d0
      dpt=5d0
      dr=0.4d0
      dptzoom=1d0

C     - Boson kinematics fully inclusive w.r.t radiation:
      call bookupeqbins(pr//'-tot-inc',1d0,0d0,1d0)
      call bookupeqbins(pr//'-y-inc',dy,-5d0,5d0)
      call bookupeqbins(pr//'-eta-inc',dy,-5d0,5d0)
      call bookupeqbins(pr//'-pt-inc',dpt,0d0,400d0)
      call bookupeqbins(pr//'-ptzoom-inc',dptzoom,1d0,101d0)
      call bookupeqbins(pr//'-m-inc',dpt,0d0,400d0)
      

         
C     -- loop over two sets of jet cuts (20 and 30 GeV) 
      do jetcut=2,2
         if (jetcut == 1) then 
            suffix = '-ptj20'
         else
            suffix = '-ptj30'
         endif

C     - Inclusive jet multiplicties
         call bookupeqbins('Njet-inc'//suffix,1d0,-0.5d0,5.5d0)
C     - Exclusive jet multiplicties
         call bookupeqbins('Njet-exc'//suffix,1d0,-0.5d0,5.5d0)
         
C     -- HT, mtw, m(jets), y_j(i), eta_j(i),ptj(i),m_j(i) for different 
C     jet multiplicities 
         do j=0,4
            if (j>0) then 
               do i=1,j
                  call bookupeqbins(
     $                 'j'//cnum(i)//'-y-ge'//cnum(j)//'j'//suffix,
     $                 dy,0d0,4d0)
                  call bookupeqbins(
     S                 'j'//cnum(i)//'-eta-ge'//cnum(j)//'j'//suffix,
     $                 dy,0d0,4d0)
                  call bookupeqbins(
     $                 'j'//cnum(i)//'-pt-ge'//cnum(j)//'j'//suffix,
     $                 dpt,30d0,180d0)
                  call bookupeqbins(
     $                 'j'//cnum(i)//'-ptzoom-ge'//
     $                 cnum(j)//'j'//suffix,dptzoom,1d0,101d0)
                  call bookupeqbins(
     $                 'j'//cnum(i)//'-m-ge'//cnum(j)//'j'//suffix,
     $                 dptzoom,1d0,101d0) 
               enddo
            endif
            
            call bookupeqbins('HT-ge'//cnum(j)//'j'//suffix,
     $           50d0,0d0,700d0)  
            call bookupeqbins('mtw-ge'//cnum(j)//'j'//suffix,
     $           5d0,0d0,120d0)  
            call bookupeqbins('massjets-ge'//cnum(j)//'j'//suffix,
     $           5d0,0d0,300d0)  
            
         enddo
         
         call bookupeqbins('yleptminusyj1-ge1j'//suffix,dy,-5d0,5d0)
         call bookupeqbins('yleptplusyj1-ge1j'//suffix,dy,-5d0,5d0)
         

         call bookupeqbins('drj1j2-ge2j'//suffix,dr,0d0,4.4d0)
         call bookupeqbins('dyj1j2-ge2j'//suffix,dy,0d0,3.5d0)
         call bookupeqbins('dphij1j2-ge2j'//suffix,pi/8d0,0d0,pi)

C     
         call bookupeqbins('mj1j2-ge2j'//suffix,dpt,0d0,300d0)

C     inclusive jet distribution 
         call bookupeqbins('inclptj'//suffix,dpt,30d0,180d0)
         call bookupeqbins('inclyj'//suffix,dy,0d0,4d0)

         call bookupeqbins(pr//'-ptzoom-ge2j'//suffix,dptzoom,1d0,101d0)

         call bookupeqbins('deltanjets'//suffix,1/30d0,0d0,1d0)
         call bookupeqbins('deltajets'//suffix,3d0,0d0,150d0)
         
      enddo

      call bookupeqbins('logweights',0.1d0,-3d0,10d0)    

      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h' 
      include 'pwhg_rad.h' 
      include 'pwhg_flg.h'
      include  'LesHouches.h'
      include 'pwhg_weights.h'
      real * 8 dsig(7),one7(7)
      integer nweights
      data one7/7*1d0/
      character * 1 pr
      common/pwhgprocess/pr
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer maxjet 
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),phi(maxjet), 
     1    phij(maxjet),pj(4,maxjet),rr,ptrel(4),pjet_in(4,maxjet) 
      character * 1 cnum(0:9)
      data cnum/'0','1','2','3','4','5','6','7','8','9'/
      real * 8 pz(4),p_el(4),p_pos(4)
      real * 8 y_el,eta_el,pt_el
      real * 8 y_pos,eta_pos,pt_pos
      real * 8 y,eta,pt,ptz,m,mll,mass_jets,psum_jets(4)
      real * 8 dy,deta,delphi,dr,yijs(4),ptmin,ht,deltajets,deltanjets
      real * 8 powheginput,dotp
      external powheginput,dotp
      integer j,k,i,jj,mjets,jetcut,ihep,njets
      character * 6 suffix 
      logical ini
      data ini/.true./
      save ini
      if(ini) then
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         ini=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
         nweights=1
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
          nweights=weights_num
      endif

      if(sum(abs(dsig)).eq.0) return

      pz = 0d0 
      p_pos = 0d0 
      p_el = 0d0 
C     select Z boson 
      if(WHCPRG.eq.'NLO') then
         pz=phep(1:4,3)+phep(1:4,4)
         p_pos = phep(1:4,3)
         p_el = phep(1:4,4)
      else
         do ihep=1,nhep
            if(isthep(ihep).eq.1.and.
     1           abs(idhep(ihep)).ge.11.and.abs(idhep(ihep)).le.16
     2           ) then
               if(abs(idhep(jmohep(1,ihep))).eq.23) then
                  pz=phep(1:4,jmohep(1,ihep))
               elseif(abs(idhep(jmohep(1,jmohep(1,ihep)))).eq.23) then
                  pz=phep(1:4,jmohep(1,jmohep(1,ihep)))
               endif
            endif
         enddo
         
C     select electron and positron 
         do ihep=1,nhep
            if(isthep(ihep).eq.1.and.
     1           idhep(ihep).eq.11 .or. idhep(ihep).eq.13 .or. 
     2           idhep(ihep).eq.15
     3           ) then
               if(abs(idhep(jmohep(1,ihep))).eq.23) then
                  p_el=phep(1:4,ihep)
               elseif(abs(idhep(jmohep(1,jmohep(1,ihep)))).eq.23) then
                  p_el=phep(1:4,ihep)
               endif
            endif
            if(isthep(ihep).eq.1.and.
     1           idhep(ihep).eq.-11 .or. idhep(ihep).eq.-13 .or. 
     2           idhep(ihep).eq.-15
     3           ) then
               if(abs(idhep(jmohep(1,ihep))).eq.23) then
                  p_pos=phep(1:4,ihep)
               elseif(abs(idhep(jmohep(1,jmohep(1,ihep)))).eq.23) then
                  p_pos=phep(1:4,ihep)
               endif
            endif
         enddo
      endif
      if (all(pz == 0d0)) then 
         write(*,*) 'Did not find Z boson'
      endif
      if (all(p_el == 0d0)) then 
         write(*,*) 'Did not find electron'
      endif
      if (all(p_pos == 0d0)) then 
         write(*,*) 'Did not find positron'
      endif



C     -- build jets in anti-kt algorithm 
      rr=0.4d0
      ptmin=20d0
      call buildjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj,yijs)
      
C     -- now loop over jet cuts (20 and 30 GeV) 
      do jetcut = 2,2 
         if (jetcut == 1) then 
            suffix = '-ptj20'
         else
            suffix = '-ptj30'
         endif         

C     -- select only jets in central rapidity and with pt_min 20 or 30 
C         + keep only jets isolated from electron and positron 
         call getyetaptmass(p_el,y_el,eta_el,pt_el,m)     
         call getyetaptmass(p_pos,y_pos,eta_pos,pt_pos,m)     
         njets = 0 
         do j=1,mjets
            if (abs(rapj(j)) < 4.4d0) then 
               if (jetcut == 1 .or. (jetcut == 2.and.ktj(j)>30d0)) then 
                  call getdydetadphidr(p_el,pj(:,j),dy,deta,delphi,dr)
                  if (dr < 0.5d0) return 
                  call getdydetadphidr(p_pos,pj(:,j),dy,
     $                 deta,delphi,dr)
                  if (dr < 0.5d0)  return 
                  njets = njets + 1
                  pjet_in(:,njets) = pj(:,j)
               endif
            endif
         enddo
         
C     -- now rename back the jets 
         do j=1,njets
            call getyetaptmass(pjet_in(:,j),rapj(j),etaj(j),ktj(j),m)
            phi(j)=atan2(pjet_in(2,j),pjet_in(1,j))
         enddo      


         
C     -- lepton cuts as in ATLAS 1111.2690 (currently e/mu average cuts) 
         mll = sqrt((p_el(4)+p_pos(4))**2-(p_el(1)+p_pos(1))**2
     .        -(p_el(2)+p_pos(2))**2-(p_el(3)+p_pos(3))**2)
         if (mll < 66d0 .or. mll > 116d0) return 

         if (pt_el < 20d0 .or. abs(eta_el) > 2.5d0) return 
         if (pt_pos < 20d0 .or. abs(eta_pos) > 2.5d0) return 

         
C     --- all cuts applied, now start filling histograms 

         if (njets.ge.2 .and. jetcut .eq.2) then 
            call filld('logweights',log10(abs(dsig(1:nweights))),one7)
         endif

         if (jetcut.eq.2) then 
C     - Boson kinematics fully inclusive w.r.t radiation:
            call getyetaptmass(pz,y,eta,ptz,m)
            call filld(pr//'-tot-inc', 0.5d0, dsig)
            call filld(pr//'-y-inc',    y, dsig)
            call filld(pr//'-eta-inc',eta, dsig)
            call filld(pr//'-pt-inc',  ptz, dsig)
            call filld(pr//'-ptzoom-inc',  ptz, dsig)
            call filld(pr//'-m-inc', m, dsig)
         endif


         if (njets .ge. 2) 
     1        call filld(pr//'-ptzoom-ge2j'//suffix,ptz,dsig)                 
      
C     -- Plot inclusive jet multiplicties at current jet pT threshold:
         if(njets.ge.0) call filld('Njet-inc'//suffix,0d0,dsig)
         if(njets.ge.1) call filld('Njet-inc'//suffix,1d0,dsig)
         if(njets.ge.2) call filld('Njet-inc'//suffix,2d0,dsig)
         if(njets.ge.3) call filld('Njet-inc'//suffix,3d0,dsig)
         if(njets.ge.4) call filld('Njet-inc'//suffix,4d0,dsig)
         if(njets.ge.5) call filld('Njet-inc'//suffix,5d0,dsig)
         
C     -- Plot exclusive jet multiplicties at current jet pT threshold:
         call filld('Njet-exc'//suffix,dble(njets),dsig)
         
C     -- compute mass(all jets) 
         psum_jets = sum(pjet_in(:,1:2),dim=2)
         mass_jets = sqrt(abs(psum_jets(4)**2
     .        -psum_jets(1)**2-psum_jets(2)**2-psum_jets(3)**2))
         
C     -- total transverse energy 
         ht = pt_el + pt_pos  
         do i=1,njets 
            ht = ht + ktj(i)
         enddo
         
         do j=0,min(4,njets)
C     -- fill inclusive pt and y distributions 
            if (j >0) then 
               call getyetaptmass(pjet_in(:,j),y,eta,pt,m)
               call filld('inclptj'//suffix,pt,dsig)
               call filld('inclyj'//suffix,abs(y),dsig)
               do i=1,j
C     -- plot variables (pt,y,eta,m) related to jet i for the j-th 
C     inclusive jet sample
C     if (njets .ge. j) then 
                  call getyetaptmass(pjet_in(:,i),y,eta,pt,m)
                  call filld('j'//cnum(i)//'-y-ge'//cnum(j)//'j'
     $//suffix,
     $                 abs(y),dsig)
                  call filld('j'//cnum(i)//'-eta-ge'//cnum(j)//'j'
     $//suffix,
     $                 abs(eta),dsig)
                  call filld('j'//cnum(i)//'-pt-ge'//cnum(j)//'j'
     $//suffix,
     $                 pt,dsig)
                  call filld('j'//cnum(i)//'-ptzoom-ge'//
     $                 cnum(j)//'j'//suffix,pt,dsig)
                  call filld('j'//cnum(i)//'-m-ge'//cnum(j)//'j'
     $//suffix,
     $                 m,dsig) 
C     endif

                  if (i == 1) then 
!                     call getyetaptmass(pjet_in(:,1),y,eta,pt,m)
                     call filld('yleptminusyj1-ge1j'//suffix,y_el-y,
     $                    dsig)
                     call filld('yleptplusyj1-ge1j'//suffix,y_el+y,dsig)
                  endif
               enddo 
            endif
            
C     -- plot HT, mll and m(jets) for the j-th inclusive jet sample 
            if (njets >= j) then 
               call filld('HT-ge'//cnum(j)//'j'//suffix,ht,dsig)
               call filld('mtw-ge'//cnum(j)//'j'//suffix,mll,dsig)
               call filld('massjets-ge'//cnum(j)//'j'//suffix,
     $              mass_jets,dsig)
            endif
            
         enddo   
         
         
         if (njets .ge. 2) then 
            call getdydetadphidr(pjet_in(:,1),pjet_in(:,2),dy,deta,
     $           delphi,dr)
            if (dr < 0.4) write(*,*) 'drj1j2', dr
            call filld('drj1j2-ge2j'//suffix,dr,dsig)
            call filld('dyj1j2-ge2j'//suffix,abs(dy),dsig)
            call filld('dphij1j2-ge2j'//suffix,delphi,dsig)
            
C     -- compute mass(j1j2) 
            psum_jets = sum(pjet_in(:,1:2),dim=2)
            mass_jets = sqrt(psum_jets(4)**2
     .           -psum_jets(1)**2-psum_jets(2)**2-psum_jets(3)**2)
            call filld('mj1j2-ge2j'//suffix,mass_jets,dsig)
         endif
         

         if(njets.eq.2) then
            deltajets=sqrt((pjet_in(1,1)+pjet_in(1,2))**2 +
     1           (pjet_in(2,1)+pjet_in(2,2))**2)
            deltanjets=deltajets/(sqrt(pjet_in(1,1)**2+pjet_in(2,1)**2)+
     1           sqrt(pjet_in(1,2)**2+pjet_in(2,2)**2))
            call filld('deltajets'//suffix,deltajets,dsig)
            call filld('deltanjets'//suffix,deltanjets,dsig)
         endif


         
      enddo
      
      
      end


      subroutine deltaplot(p1,p2,dsig,prefix,postfix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix,postfix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'-dy'//postfix,dy,dsig)
      call filld(prefix//'-deta'//postfix,deta,dsig)
      call filld(prefix//'-delphi'//postfix,delphi,dsig)
      call filld(prefix//'-dr'//postfix,dr,dsig)
      end


      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
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
      dr=sqrt(dy**2+dphi**2)
      end


      subroutine buildjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet,yijs)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet),pjet_in(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,jb,i
      real * 8 r,palg,pp,tmp
      logical islept
      external islept
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
      real * 8 yijs(4)
      character * 1 pr
      common/pwhgprocess/pr
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
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
c all but the Boson
            if (isthep(j).eq.1) then
               if(pr.eq.'H') then
                  if(idhep(j).eq.25) cycle
               elseif(pr.eq.'Z') then
                  if(WHCPRG.eq.'NLO') then
                     if(j.eq.3.or.j.eq.4) cycle
                  else
                     if(idhep(jmohep(1,j)).eq.23) cycle
                     if(jmohep(1,jmohep(1,j)).ne.0) then
                        if(idhep(jmohep(1,jmohep(1,j))).eq.23) cycle
                     endif
                  endif
               elseif(pr.eq.'W') then
                  if(WHCPRG.eq.'NLO') then
                     if(j.eq.3.or.j.eq.4) cycle
                  else
                     if(abs(idhep(jmohep(1,j))).eq.24) cycle
                     if(jmohep(1,jmohep(1,j)).ne.0) then
                        if(abs(idhep(jmohep(1,jmohep(1,j)))).eq.24)cycle
                     endif
                  endif
               endif
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
C     R = 0.7   radius parameter
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c      ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec,yijs)
      mjets=njets
      if(njets.eq.0) return
Cc check consistency
C      do k=1,ntracks
C         if(jetvec(k).gt.0) then
C            do mu=1,4
C               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
C            enddo
C         endif
C      enddo
C      tmp=0
C      do j=1,mjets
C         do mu=1,4
C            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
C         enddo
C      enddo
C      if(tmp.gt.1d-4) then
C         write(*,*) ' bug!'
C      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- 
      do j=1,mjets
         call getyetaptmass(pjet(:,j),rap(j),eta(j),kt(j),tmp)
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

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      subroutine pwhgfinalopshist
      end
