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

      call inihists

c     total cross section sanity check
      call bookupeqbins('total',1d0,0d0,1d0)

      call bookupeqbins('lep1_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1_pt',2d0,0d0,400d0)
      call bookupeqbins('lep1_m',2d0,0d0,400d0)

      call bookupeqbins('lep2_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2_pt',2d0,0d0,400d0)
      call bookupeqbins('lep2_m',2d0,0d0,400d0)

      call bookupeqbins('alp1_y',0.2d0,-4d0,4d0)
      call bookupeqbins('alp1_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('alp1_pt',2d0,0d0,400d0)
      call bookupeqbins('alp1_m',2d0,0d0,400d0)

      call bookupeqbins('alp2_y',0.2d0,-4d0,4d0)
      call bookupeqbins('alp2_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('alp2_pt',2d0,0d0,400d0)
      call bookupeqbins('alp2_m',2d0,0d0,400d0)

      call bookupeqbins('lep1alp1_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp1_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp1_pt',2d0,0d0,400d0)
      call bookupeqbins('lep1alp1_m',2d0,0d0,400d0)

      call bookupeqbins('lep1alp2_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp2_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp2_pt',2d0,0d0,400d0)
      call bookupeqbins('lep1alp2_m',2d0,0d0,400d0)

      call bookupeqbins('lep2alp1_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp1_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp1_pt',2d0,0d0,400d0)
      call bookupeqbins('lep2alp1_m',2d0,0d0,400d0)

      call bookupeqbins('lep2alp2_y',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp2_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp2_pt',2d0,0d0,400d0)
      call bookupeqbins('lep2alp2_m',2d0,0d0,400d0)

      call bookupeqbins('4l_y',0.2d0,-4d0,4d0)
      call bookupeqbins('4l_eta',0.2d0,-4d0,4d0)
      call bookupeqbins('4l_pt',2d0,0d0,400d0)
      call bookupeqbins('4l_m',2d0,0d0,400d0)

      call bookupeqbins('j_y_20cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_y_40cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_y_60cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_y_100cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_y_200cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_y_400cut',0.2d0,-4d0,4d0)
      
      call bookupeqbins('j_dy_20cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_dy_40cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_dy_60cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_dy_100cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_dy_200cut',0.2d0,-4d0,4d0)
      call bookupeqbins('j_dy_400cut',0.2d0,-4d0,4d0)

      call bookupeqbins('lep1lep2_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1lep2_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1lep2_delphi',pi/20,0d0,pi)
      call bookupeqbins('lep1lep2_dr',0.2d0,-8d0,8d0)
      
      call bookupeqbins('alp1alp2_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('alp1alp2_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('alp1alp2_delphi',pi/20,0d0,pi)
      call bookupeqbins('alp1alp2_dr',0.2d0,-8d0,8d0)
      
      call bookupeqbins('lep1alp1_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp1_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp1_delphi',pi/20,0d0,pi)
      call bookupeqbins('lep1alp1_dr',0.2d0,-8d0,8d0)
      
      call bookupeqbins('lep1alp2_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp2_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep1alp2_delphi',pi/20,0d0,pi)
      call bookupeqbins('lep1alp2_dr',0.2d0,-8d0,8d0)
      
      call bookupeqbins('lep2alp1_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp1_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp1_delphi',pi/20,0d0,pi)
      call bookupeqbins('lep2alp1_dr',0.2d0,-8d0,8d0)
      
      call bookupeqbins('lep2alp2_dy',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp2_deta',0.2d0,-4d0,4d0)
      call bookupeqbins('lep2alp2_delphi',pi/20,0d0,pi)
      call bookupeqbins('lep2alp2_dr',0.2d0,-8d0,8d0)
      
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
c     binsize
c     we need to tell to this analysis file which program is running it
      integer   maxjet,mjets
      parameter (maxjet=2048)
      real * 8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1    phij(maxjet),pj(4,maxjet),rr
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer vdecaytemp1,vdecaytemp2
      save vdecaytemp1,vdecaytemp2

      integer lep1, lep2, alp1,alp2
      integer nlep1,nlep2,nalp1,nalp2
      integer ilep1(100),ilep2(100),ialp1(100),ialp2(100)
      real * 8 plep1(4),plep2(4),palp1(4),palp2(4)
      real * 8 httot,y,eta,pt,m
      real * 8 dy,deta,delphi,dr
      integer ihep
      real * 8 powheginput
      external powheginput
      real * 8 mllmin
      save mllmin, lep1, lep2, alp1, alp2

      if(dsig.eq.0) return

      if (ini) then
         if(powheginput("#vdecaymodeZ1").gt.0) then
            vdecaytemp1=powheginput("#vdecaymodeZ1")
            vdecaytemp2=powheginput("#vdecaymodeZ2")
            lep1= vdecaytemp1
            alp1=-vdecaytemp1
            lep2= vdecaytemp2
            alp2=-vdecaytemp2
         elseif(powheginput("#vdecaymodeZ").gt.0) then
            vdecaytemp1=powheginput("#vdecaymodeZ")
            vdecaytemp2=powheginput("#vdecaymodeW")
            lep1= vdecaytemp1
            alp1=-vdecaytemp1
            if(vdecaytemp2.gt.0) then
               lep2= vdecaytemp2
               alp2=-(vdecaytemp2+1)
            else
               lep2=-vdecaytemp2+1
               alp2= vdecaytemp2
            endif
         elseif(powheginput("#vdecaymodeWm").gt.0) then
            vdecaytemp1=powheginput("#vdecaymodeWm")
            vdecaytemp2=powheginput("#vdecaymodeWp")
            lep1= vdecaytemp1
            alp1=-(vdecaytemp1+1)
            lep2= -vdecaytemp2+1
            alp2= vdecaytemp2
         else
            write(*,*) ' which analysis?'
            call exit(-1)
         endif

         
            
         mllmin=powheginput('#mllmin')
         write (*,*)
         write (*,*) '********************************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif

         ini=.false.
      endif


c     find Z decay products
      nlep1=0                   ! number of lepton 1
      nalp1=0                   ! number of antilepton 1
      nlep2=0                   ! number of l-
      nalp2=0                   ! number of l+
      do ihep=1,nhep
         if(isthep(ihep).eq.1) then
            if(idhep(ihep).eq.lep1) then
               nlep1=nlep1+1
               ilep1(nlep1) = ihep
            elseif(idhep(ihep).eq.alp1) then
               nalp1=nalp1+1
               ialp1(nalp1) = ihep
            elseif(idhep(ihep).eq.lep2) then
               nlep2=nlep2+1
               ilep2(nlep2) = ihep
            elseif(idhep(ihep).eq.alp2) then
               nalp2=nalp2+1
               ialp2(nalp2) = ihep
            endif
         endif
      enddo

      if (nlep1 .ge. 100 .or. nalp1 .ge. 100 .or. nlep2 .ge. 100 .or.
     .     nalp2 .ge. 100) then
         write(*,*) 'crazy event, too many leptons'
         return
      endif
      
      call sortbypt(nlep1,ilep1)
      call sortbypt(nalp1,ialp1)
      call sortbypt(nlep2,ilep2)
      call sortbypt(nalp2,ialp2)

      if(lep1.eq.lep2) then
         if(nlep1.lt.2) then
            write(*,*) 'crazy event, not enough leptons'
            call exit(-1)
         endif
         nlep2=1
         ilep2(1)=ilep1(2)
         nlep1=1
      endif
      if(alp1.eq.alp2) then
         if(nalp1.lt.2) then
            write(*,*) 'crazy event, not enough leptons'
            call exit(-1)
         endif
         nalp2=1
         ialp2(1)=ialp1(2)
         nalp1=1
      endif

      call filld('total',0.5d0,dsig)

      rr=0.6d0 
      call buildjets(1,mjets,rr,ktj,etaj,rapj,phij,pj)

      plep1=phep(1:4,ilep1(1))
      palp1=phep(1:4,ialp1(1))
      plep2=phep(1:4,ilep2(1))
      palp2=phep(1:4,ialp2(1))

      call yetaptmassplot(plep1,dsig,'lep1')
      if(lep1.eq.lep2) then
         call yetaptmassplot(plep2,dsig,'lep1')
      else
         call yetaptmassplot(plep2,dsig,'lep2')
      endif
      call yetaptmassplot(palp1,dsig,'alp1')
      if(alp1.eq.alp2) then
         call yetaptmassplot(palp2,dsig,'alp1')
      else
         call yetaptmassplot(palp2,dsig,'alp2')
      endif

      if(lep2.eq.lep1.and.alp2.eq.alp1) then
         call yetaptmassplot(plep1+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep1+palp2,dsig,'lep1alp1')
         call yetaptmassplot(plep2+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep2+palp2,dsig,'lep1alp1')
      elseif(lep1.eq.lep2) then
         call yetaptmassplot(plep1+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep1+palp2,dsig,'lep1alp2')
         call yetaptmassplot(plep2+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep2+palp2,dsig,'lep1alp2')
      elseif(alp1.eq.alp2) then
         call yetaptmassplot(plep1+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep1+palp2,dsig,'lep1alp1')
         call yetaptmassplot(plep2+palp1,dsig,'lep2alp1')
         call yetaptmassplot(plep2+palp2,dsig,'lep2alp1')
      else
         call yetaptmassplot(plep1+palp1,dsig,'lep1alp1')
         call yetaptmassplot(plep1+palp2,dsig,'lep1alp2')
         call yetaptmassplot(plep2+palp1,dsig,'lep2alp1')
         call yetaptmassplot(plep2+palp2,dsig,'lep2alp2')
      endif

         
      call yetaptmassplot(plep1+plep2+palp1+palp2,dsig,'4l')

c jets
      if(mjets.gt.0) then
         if(ktj(1).gt.20) call filld('j_y_20cut',rapj(1),dsig)
         if(ktj(1).gt.40) call filld('j_y_40cut',rapj(1),dsig)
         if(ktj(1).gt.60) call filld('j_y_60cut',rapj(1),dsig)
         if(ktj(1).gt.100) call filld('j_y_100cut',rapj(1),dsig)
         if(ktj(1).gt.200) call filld('j_y_200cut',rapj(1),dsig)
         if(ktj(1).gt.400) call filld('j_y_400cut',rapj(1),dsig)
         call getdydetadphidr(plep1+plep2+palp1+palp2,pj(:,1),
     1        dy,deta,delphi,dr)
         if(ktj(1).gt.20) call filld('j_dy_20cut',dy,dsig)
         if(ktj(1).gt.40) call filld('j_dy_40cut',dy,dsig)
         if(ktj(1).gt.60) call filld('j_dy_60cut',dy,dsig)
         if(ktj(1).gt.100) call filld('j_dy_100cut',dy,dsig)
         if(ktj(1).gt.200) call filld('j_dy_200cut',dy,dsig)
        if(ktj(1).gt.400) call filld('j_dy_400cut',dy,dsig)
      endif

      call deltaplot(plep1,plep2,dsig,'lep1lep2')

      call deltaplot(plep1,palp1,dsig,'alp1alp2')

      if(lep2.eq.lep1.and.alp1.eq.alp2) then
         call deltaplot(plep1,palp1,dsig,'lep1alp1')
         call deltaplot(plep1,palp2,dsig,'lep1alp1')
         call deltaplot(plep2,palp1,dsig,'lep1alp1')
         call deltaplot(plep2,palp2,dsig,'lep1alp1')
      elseif(alp2.eq.alp1) then
         call deltaplot(plep1,palp1,dsig,'lep1alp1')
         call deltaplot(plep1,palp2,dsig,'lep1alp1')
         call deltaplot(plep2,palp1,dsig,'lep2alp1')
         call deltaplot(plep2,palp2,dsig,'lep2alp1')
      elseif(lep1.eq.lep2) then
         call deltaplot(plep1,palp1,dsig,'lep1alp1')
         call deltaplot(plep1,palp2,dsig,'lep1alp2')
         call deltaplot(plep2,palp1,dsig,'lep1alp1')
         call deltaplot(plep2,palp2,dsig,'lep1alp2')
      else
         call deltaplot(plep1,palp1,dsig,'lep1alp1')
         call deltaplot(plep1,palp2,dsig,'lep1alp2')
         call deltaplot(plep2,palp1,dsig,'lep2alp1')
         call deltaplot(plep2,palp2,dsig,'lep2alp2')
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
