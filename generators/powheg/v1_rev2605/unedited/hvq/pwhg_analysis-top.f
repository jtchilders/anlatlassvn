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
      integer diag
      real * 8 binsize(700)
      common/pwhghistcommon/binsize
      integer icut
      character * 3 cut

      call pwhginihist

      diag=0

c mass of the pair
      diag=diag+1
      binsize(diag)=10
      call pwhgbookup(diag,'mass-t-tbar','LOG',binsize(diag),0d0,1000d0)
c pt of top
      diag=diag+1
      binsize(diag)=10
      call pwhgbookup(diag,'t pt','LOG',binsize(diag),0d0,1000d0)
c pt asymmetry
      diag=diag+1
      binsize(diag)=10
      call pwhgbookup(diag,'t - tbar,pt','LOG',binsize(diag),0d0,1000d0)
c rapidity
      diag=diag+1
      binsize(diag)=0.2
      call pwhgbookup(diag,'t y','LOG',binsize(diag),-3d0,3d0)
c rapidity asymmetry
      diag=diag+1
      binsize(diag)=0.2
      call pwhgbookup(diag,'t - tbar,y','LOG',binsize(diag),-3d0,3d0)
c (Non b) jet pt
      diag=diag+1
      binsize(diag)=2
      call pwhgbookup(diag,'jet pt','LOG',binsize(diag),0d0,200d0)
c (Non b) jet y
      do icut=20,100,20
         write(cut,'(i3)') icut
         diag=diag+1
         binsize(diag)=0.2
         call pwhgbookup(diag,'yjet,pt>'//cut,
     1        'LOG',binsize(diag),-4d0,4d0)
      enddo
c (Non b) jet y-yttbar
      do icut=20,100,20
         write(cut,'(i3)') icut
         diag=diag+1
         binsize(diag)=0.2
         call pwhgbookup(diag,'yjet-yttbar,pt>'//cut,
     1        'LOG',binsize(diag),-4d0,4d0)
      enddo
c Lepton+
c Lepton energy in top rest frame 
      diag=diag+1
      binsize(diag)=1
      call pwhgbookup(diag,'eem','LOG',binsize(diag),0d0,100d0)
c Lepton cos theta in top rest frame 
      diag=diag+1
      binsize(diag)=0.1
      call pwhgbookup(diag,'cos thdec','LOG',binsize(diag),-1d0,1d0)
c Lepton azimuth in top rest frame. The origin is the top
c transverse direction
      diag=diag+1
      binsize(diag)=pi/20
      call pwhgbookup(diag,'phi dec','LOG',binsize(diag),-pi,pi)
      
      end



      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0,dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      logical ini
      data ini/.true./
      save ini
      integer diag,icut,j,mjets
      real * 8 pjet(4),ppairttbar(4),ptt,yt,etat,pttbar,ytbar,etatbar,
     1     ptttbar,mttbar,yttbar,etattbar,ptjet,yjet,etajet,eep,cthep,
     1     phiep,eem,cthem,phiem,dphi
      real * 8 binsize(700)
      common/pwhghistcommon/binsize
      integer ihep,it,itbar,iem,iep,ijet,mu
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/
      integer   maxjet
      parameter (maxjet=2048)
      real * 8  kt(maxjet),eta(maxjet),rap(maxjet),
     c    phi(maxjet),pj(4,maxjet)
      integer ibtag(maxjet),itags
      logical sonoftop
      external sonoftop
      dsig=dsig0
      it=0
      itbar=0
      iem=0
      iep=0
      ijet=0
      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         do ihep=3,nhep
            if(idhep(ihep).eq.6) then
               it=ihep
            elseif(idhep(ihep).eq.-6) then
               itbar=ihep
            elseif(idhep(ihep).eq.11) then
               iem=ihep
            elseif(idhep(ihep).eq.-11) then
               iep=ihep
            elseif(abs(idhep(ihep)).lt.6.or.idhep(ihep).eq.21) then
               ijet=ihep
            endif
         enddo
         do mu=1,4
            pjet(mu)=phep(mu,ijet)
         enddo
      else
         do ihep=1,nhep
            if(isthep(ihep).eq.1) then
               if(idhep(ihep).eq.11.and.sonoftop(ihep,j)) then
                  iem=ihep
                  itbar=j
               elseif(idhep(ihep).eq.-11.and.sonoftop(ihep,j)) then
                  iep=ihep
                  it=j
               endif
            endif
         enddo
         mjets=10
         call buildjets(mjets,kt,eta,rap,phi,pj,ibtag)
         do j=1,mjets
            if(ibtag(j).eq.0) then
               do mu=1,4
                  pjet(mu)=pj(mu,j)
               enddo
               goto 10
            endif
         enddo
      endif
 10   continue
      call ptyeta(phep(1,it),ptt,yt,etat)
      call ptyeta(phep(1,itbar),pttbar,ytbar,etatbar)
      if(ijet.eq.0) then
         ptjet=0
      else
         call ptyeta(phep(1,ijet),ptjet,yjet,etajet)
      endif
      diag=0
c mass of the pair
      diag=diag+1
      do mu=1,4
         ppairttbar(mu)=phep(mu,it)+phep(mu,itbar)
      enddo
      mttbar=sqrt(ppairttbar(4)**2
     1     -ppairttbar(1)**2-ppairttbar(2)**2-ppairttbar(3)**2)
      call ptyeta(ppairttbar,ptttbar,yttbar,etattbar)
      call pwhgfill(diag,mttbar,dsig/binsize(diag))
c pt of top
      diag=diag+1
      call pwhgfill(diag,ptt,dsig/binsize(diag))
c pt asymmetry
      diag=diag+1
      call pwhgfill(diag,ptt,dsig/binsize(diag))
      call pwhgfill(diag,pttbar,-dsig/binsize(diag))
c rapidity
      diag=diag+1
      call pwhgfill(diag,yt,dsig/binsize(diag))
c rapidity asymmetry
      diag=diag+1
      call pwhgfill(diag,yt,dsig/binsize(diag))
      call pwhgfill(diag,ytbar,-dsig/binsize(diag))
c (Non b) jet pt
      diag=diag+1
      call pwhgfill(diag,ptjet,dsig/binsize(diag))
c (Non b) jet y
      do icut=20,100,20
         diag=diag+1
         if(ptjet.gt.icut) then
            call pwhgfill(diag,yjet,dsig/binsize(diag))
         endif
      enddo
c (Non b) jet y-yttbar
      do icut=20,100,20
         diag=diag+1
         if(ptjet.gt.icut) then
            call pwhgfill(diag,yjet-yttbar,dsig/binsize(diag))
         endif
      enddo
c leptons may not be there, in case we are called by the NLO
c section
      if(iep.eq.0.or.iem.eq.0) return         
c Lepton+
c Compute lepton variables in top rest frame
      call decvariables(phep(1,iep),phep(1,it),eep,cthep,phiep)
      call decvariables(phep(1,iem),phep(1,itbar),eem,cthem,phiem)
c Lepton energy in top rest frame 
      diag=diag+1
      call pwhgfill(diag,eep,dsig/binsize(diag))
c Lepton cos theta in top rest frame 
      diag=diag+1
      call pwhgfill(diag,cthep*cthem,dsig/binsize(diag))
c Lepton azimuth in top rest frame. The origin is the top
c transverse direction
      diag=diag+1
      dphi=(phiep-phiem)-nint((phiep-phiem)/(2*pi))*2*pi
      call pwhgfill(diag,dphi,dsig/binsize(diag))
      end


      subroutine decvariables(pdec0,ppart0,edec,cthdec,phidec)
      implicit none
      include 'pwhg_math.h'
      real * 8 pdec0(4),ppart0(4),edec,cthdec,phidec
      real * 8 pdec(0:3),ppart(0:3)
      real * 8 vec(3),beta,pt
      integer mu
      do mu=1,3
         pdec(mu)=pdec0(mu)
         ppart(mu)=ppart0(mu)
      enddo
      pdec(0)=pdec0(4)
      ppart(0)=ppart0(4)
      vec(1)=0
      vec(2)=0
      vec(3)=-1
      beta=ppart(3)/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      pt=sqrt(ppart(1)**2+ppart(2)**2)
      vec(1)=ppart(1)/pt
      vec(2)=ppart(2)/pt
      vec(3)=0
      beta=-pt/ppart(0)
      call mboost(1,vec,beta,pdec,pdec)
      call mboost(1,vec,beta,ppart,ppart)
      edec=pdec(0)
      cthdec=pdec(3)/sqrt(pdec(1)**2+pdec(2)**2+pdec(3)**2)
      phidec=atan2(pdec(2),pdec(1))-atan2(vec(2),vec(1))
c bring it back between -pi and pi
      phidec=phidec-nint(phidec/(2*pi))*2*pi
      end

      subroutine ptyeta(p,pt,y,eta)
      implicit none
      real * 8 p(4),pt,y,eta
      real * 8 pp,tiny
      parameter (tiny=1d-12)
      pt=sqrt(p(1)**2+p(2)**2)
      y=log((p(4)+p(3))/(p(4)-p(3)))/2
      pp=sqrt(pt**2+p(3)**2)*(1+tiny)
      eta=log((pp+p(3))/(pp-p(3)))/2
      end





      subroutine buildjets(mjets,kt,eta,rap,phi,pjet,ibtag)
c     arrays to reconstruct jets
      implicit none
      integer mjets
      real * 8  kt(mjets),eta(mjets),rap(mjets),phi(mjets),pjet(4,mjets)
      integer ibtag(mjets)
      include   'hepevt.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=20)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8 r,palg,ptmin,pp,tmp
      logical bson
C - Initialize arrays and counters for output jets
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
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if (isthep(j).eq.1) then
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
         return
      endif
C --------------------------------------------------------------------- C
C - Inclusive jet pT and Y spectra are to be compared to CDF data:    - C    
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
C     f = 0.75  overlapping fraction
      palg=-1
      r=0.5d0
      ptmin=15
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c Find b decay products among tracks
      do j=1,njets
         ibtag(j)=0
      enddo
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
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            if(bson(itrackhep(k))) then
               ibtag(jetvec(k))=ibtag(jetvec(k))+1
            endif
         endif
      enddo
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

      logical function bson(j)
      implicit none
      integer J
      include   'hepevt.h'
      integer jcurr
      logical bhadr
      jcurr=j
c     This only happens in parton level analysis
      if(abs(idhep(jcurr)).eq.5) then
         bson=.true.
         return
      endif
 1    continue
      bson=.false.
      if(bhadr(idhep(jcurr))) then
         bson=.true.
         return
      endif
      jcurr=jmohep(1,jcurr)
      if(idhep(jcurr).eq.0) then
         bson=.false.
         return
      endif
      goto 1
      end

      logical function bhadr(idhep)
      implicit none
      integer idhep
      integer i1,i2,idigit
      i1=idigit(1,idhep)
      if(i1.eq.1) then
c         is a bottomed meson
         i2=idigit(3,idhep)
      elseif(i1.eq.2) then
c is a bottomed barion
         i2=idigit(5,idhep)
      endif
      if(i2.eq.5) then
         bhadr=.true.
      else
         bhadr=.false.
      endif
      end

      
      logical function sonoftop(j,jtop)
      implicit none
      integer j,jtop
      include   'hepevt.h'
      integer jcurr
      logical bhadr
      jcurr=j
c     This only happens in parton level analysis
 1    continue
      if(abs(idhep(jcurr)).eq.6) then
         sonoftop=.true.
         jtop=jcurr
         return
      endif
      jcurr=jmohep(1,jcurr)
      if(idhep(jcurr).eq.0) then
         sonoftop=.false.
         jtop=0
         return
      endif
      goto 1
      end

      
      function idigit(k,l)
      implicit none
      integer idigit,k,l
      idigit=abs(mod(l,10**k)/10**(k-1))
      end

