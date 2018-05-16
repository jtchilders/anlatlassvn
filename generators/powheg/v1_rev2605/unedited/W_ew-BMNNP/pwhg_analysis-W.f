c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,k
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      integer maxjet
      parameter (maxjet=1)

      call inihists

c     total cross section sanity check
      call bookupeqbins('Njet',1d0,-0.5d0,5.5d0)

      call bookupeqbins('V_y',0.2d0,-4d0,4d0)
      call bookupeqbins('V_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('V_pt',2d0,0d0,400d0)
      call bookupeqbins('V_m',2d0,0d0,400d0)

      do j=1,maxjet
         call bookupeqbins('j'//cnum(j)//'_y',0.2d0,-4d0,4d0)
         call bookupeqbins('j'//cnum(j)//'_eta',0.2d0,-4d0,4d0)
         call bookupeqbins('j'//cnum(j)//'_pt',2d0,0d0,400d0)
         call bookupeqbins('j'//cnum(j)//'_m',2d0,0d0,400d0)   
      enddo

      do j=1,2
      do k=j+1,maxjet
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_y',0.2d0,-4d0,4d0)  
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_eta',0.2d0,-4d0,4d0)
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_pt',2d0,0d0,400d0)  
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_m',2d0,0d0,400d0)   
      enddo
      enddo

      do j=1,maxjet
         call bookupeqbins('Vj'//cnum(j)//'_dy',0.2d0,-4d0,4d0)
         call bookupeqbins('Vj'//cnum(j)//'_deta',0.2d0,-4d0,4d0)
         call bookupeqbins('Vj'//cnum(j)//'_delphi',pi/20,0d0,pi)
         call bookupeqbins('Vj'//cnum(j)//'_dr',0.2d0,0d0,20d0)  
      enddo

      do j=1,2
      do k=j+1,maxjet
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_dy',0.2d0,-4d0,4d0)  
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_deta',0.2d0,-4d0,4d0)
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_delphi',pi/20,0d0,pi)
         call bookupeqbins('j'//cnum(j)//'j'//cnum(k)//
     1        '_dr',0.2d0,0d0,20d0)  
      enddo
      enddo
      
      do j=1,2
      do k=j+1,maxjet
         call bookupeqbins('Vj'//cnum(j)//'_j'//cnum(k)//
     1        '_dy',0.2d0,-4d0,4d0)  
         call bookupeqbins('Vj'//cnum(j)//'_j'//cnum(k)//
     1        '_deta',0.2d0,-4d0,4d0)
         call bookupeqbins('Vj'//cnum(j)//'_j'//cnum(k)//
     1        '_delphi',pi/20,0d0,pi)
         call bookupeqbins('Vj'//cnum(j)//'_j'//cnum(k)//
     1        '_dr',0.2d0,0d0,20d0)  
      enddo
      enddo

      if(maxjet.ge.3) then
         call bookupeqbins('Vj1j2_j3_dy',0.2d0,-4d0,4d0)  
         call bookupeqbins('Vj1j2_j3_deta',0.2d0,-4d0,4d0)
         call bookupeqbins('Vj1j2_j3_delphi',pi/20,0d0,pi)
         call bookupeqbins('Vj1j2_j3_dr',0.2d0,0d0,20d0)
      endif
      
      end
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
      integer   maxjet,mjets
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr
      character * 1 cnum(9)
      data cnum/'1','2','3','4','5','6','7','8','9'/
      save cnum
      integer j,k
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      real * 8 ph(4)
      real * 8 httot,y,eta,pt,m
      real * 8 dy,deta,delphi,dr
      integer ihep
      real * 8 powheginput,dotp
      external powheginput,dotp

      if(dsig.eq.0) return

      do ihep=1,nhep
         if(idhep(ihep).eq.25) then
            ph=phep(1:4,ihep)
         endif
      enddo

      rr=0.5d0 
      call buildjets(1,mjets,rr,ktj,etaj,rapj,phij,pj)

      if(mjets.eq.0) then
         call filld('Njet',0d0,dsig)
      elseif(mjets.eq.1) then
         call filld('Njet',1d0,dsig)
      elseif(mjets.eq.2) then
         call filld('Njet',2d0,dsig)
      elseif(mjets.eq.3) then
         call filld('Njet',3d0,dsig)
      elseif(mjets.eq.4) then
         call filld('Njet',4d0,dsig)
      elseif(mjets.eq.5) then
         call filld('Njet',5d0,dsig)
      else
c         write(*,*) ' Njet?',mjets
      endif


      mjets=min(1,mjets)

c Higgs
      call getyetaptmass(ph,y,eta,pt,m)
      call filld('V_y',    y, dsig)
      call filld('V_eta',eta, dsig)
      call filld('V_pt',  pt, dsig)
      call filld('V_m', m, dsig)
c jets
      do j=1,mjets
         call getyetaptmass(pj(:,j),y,eta,pt,m)
         call filld('j'//cnum(j)//'_y',     y, dsig)
         call filld('j'//cnum(j)//'_eta', eta, dsig)
         call filld('j'//cnum(j)//'_pt',   pt, dsig)
         call filld('j'//cnum(j)//'_m',     m, dsig)
      enddo
      do j=1,mjets
         do k=j+1,mjets
            call getyetaptmass(pj(:,j)+pj(:,k),y,eta,pt,m)
            call filld('j'//cnum(j)//'j'//cnum(k)//'_y',    y, dsig)
            call filld('j'//cnum(j)//'j'//cnum(k)//'_eta',eta, dsig)
            call filld('j'//cnum(j)//'j'//cnum(k)//'_pt',  pt, dsig)
            call filld('j'//cnum(j)//'j'//cnum(k)//'_m', m, dsig)
         enddo
      enddo

      do j=1,mjets
         call deltaplot(ph,pj(:,j),dsig,'Vj'//cnum(j))
      enddo

      do j=1,mjets
         do k=j+1,mjets
            call deltaplot(pj(:,j),pj(:,k),dsig,
     1           'j'//cnum(j)//'j'//cnum(k))
         enddo
      enddo

      do j=1,mjets
         do k=j+1,mjets
            call deltaplot(ph+pj(:,j),pj(:,k),dsig,
     1           'Vj'//cnum(j)//'_j'//cnum(k))
         enddo
      enddo
      if(mjets.ge.3) then
         call deltaplot(ph+pj(:,1)+pj(:,2),pj(:,3),dsig,'hj1j2_j3')
      endif

      end

      subroutine yetaptmassplot(p,dsig,prefix)
      implicit none
      real * 8 p(4),dsig
      character *(*) prefix
      real * 8 y,eta,pt,m
      call getyetaptmass(p,y,eta,pt,m)
      call filld(prefix//'_y',y,dsig)
      call filld(prefix//'_eta',eta,dsig)
      call filld(prefix//'_pt',pt,dsig)
      call filld(prefix//'_m',m,dsig)
      end

      subroutine deltaplot(p1,p2,dsig,prefix)
      implicit none
      real * 8 p1(4),p2(4),dsig
      character *(*) prefix
      real * 8 dy,deta,delphi,dr
      call getdydetadphidr(p1,p2,dy,deta,delphi,dr)
      call filld(prefix//'_dy',dy,dsig)
      call filld(prefix//'_deta',deta,dsig)
      call filld(prefix//'_delphi',delphi,dsig)
      call filld(prefix//'_dr',dr,dsig)
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
      dr=sqrt(deta**2+dphi**2)
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
      integer  n_v_leptons
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
         n_v_leptons=0
         do j=1,nhep
c all but the Higgs
            if (isthep(j).eq.1) then
               if(islept(idhep(j)).and.(idhep(jmohep(1,j)).eq.23
     1              .or.abs(idhep(jmohep(1,j))).eq.24)) then
                  n_v_leptons=n_v_leptons+1
                  if(n_v_leptons.gt.2) then
                     write(*,*) 'Fatal error found other than TWO'
                     write(*,*) 'leptons claiming to have the W/Z'
                     write(*,*) 'as their father. Illegitimate child'
                     write(*,*) 'detected! Quitting!'
                     call exit(-1)
                  endif
                  cycle
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
CXXX      ptmin=5d0 
      ptmin=20d0 
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
