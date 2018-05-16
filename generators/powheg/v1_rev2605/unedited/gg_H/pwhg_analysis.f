c     The next subroutines, open some histograms and prepare them 
c     to receive data 
c     You can substitute these  with your favourite ones
c     init   :  opens the histograms
c     topout :  closes them
c     pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include '../pwhg_book.h'
      include 'PhysPars.h'
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
      character * 20 cuts
      real *8 emmin,emmax
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      real *8 powheginput
      external powheginput

      call pwhginihist

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         EMMIN=sqrt(ph_Hmass2low)
         EMMAX=sqrt(ph_Hmass2high)
      elseif((WHCPRG.eq.'HERWIG').or.(WHCPRG.eq.'PYTHIA')
     $.or.(WHCPRG.eq.'LHE   ')) then
         EMMIN=powheginput('hmass')-10*powheginput('hwidth')
         EMMAX=powheginput('hmass')+10*powheginput('hwidth')
      endif
      
      diag=1
      binsize(diag) = 1d0
      call pwhgbookup(diag,'total','LOG',binsize(diag),0d0,3d0)

      diag=diag+1
      binsize(diag) = (emmax-emmin)/50d0
      call pwhgbookup(diag,'m(H)','LOG',binsize(diag),emmin,emmax)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y(H)','LOG',binsize(diag),-3d0,3d0)

      diag=diag+1
      binsize(diag) = 6d0
      call pwhgbookup(diag,'pt(H)','LOG',binsize(diag),0d0,300d0)

      cuts=' Pt > 10 GeV '
      diag=diag+1
      binsize(diag) = 1d-1
      call pwhgbookup(diag,'log10(pi-phi)'//cuts,'LOG',binsize(diag),
     $     -4d0 ,1d0)

      diag=diag+1
      binsize(diag) = 6d0
      call pwhgbookup(diag,'pt jet'//cuts,'LOG',binsize(diag),0d0,300d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y jet'//cuts,'LOG',binsize(diag),-3d0,3d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Dy H-jet'//cuts,'LOG',binsize(diag),-3d0
     $     ,3d0)

      cuts=' Pt > 40 GeV '

      diag=diag+1
      binsize(diag) = 6d0
      call pwhgbookup(diag,'pt jet'//cuts,'LOG',binsize(diag),0d0,300d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y jet'//cuts,'LOG',binsize(diag),-3d0,3d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Dy H-jet'//cuts,'LOG',binsize(diag),-3d0
     $     ,3d0)

      cuts=' Pt > 80 GeV '

      diag=diag+1
      binsize(diag) = 6d0
      call pwhgbookup(diag,'pt jet'//cuts,'LOG',binsize(diag),0d0,300d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'y jet'//cuts,'LOG',binsize(diag),-3d0,3d0)

      diag=diag+1
      binsize(diag) = 0.2d0
      call pwhgbookup(diag,'Dy H-jet'//cuts,'LOG',binsize(diag),-3d0
     $     ,3d0)

      end

      
     
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'PhysPars.h'
      real*8 m1,pt1,y1,delphi
      real*8 e1,px1,py1,pz1,p1
      logical ini
      data ini/.true./
      save ini
c     binsize
      integer diag
      real * 8 binsize(100)
      common/pwhghistcommon/binsize
c     we need to tell to this analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer i
c arrays to reconstruct jets
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 pjet(4,maxjet) 
      real *8 ptjet,yjet,tmp
      real * 8 R,ptmin_fastkt
      integer jetvec(maxtrack),j1
      integer mu,jpart,jjet,found,njets,
     1     ihep,ntracks,ijet
      logical buildjets
      parameter (buildjets=.true.)
      real * 8 cut(3)
      data cut/10d0,40d0,80d0/
      
 
      if (ini) then
         write(*,*) '*****************************'
         if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
            write(*,*) '       NLO ANALYSIS'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE ANALYSIS'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         write(*,*) '*****************************'
         ini=.false.
      endif

      found=0

c     Loop over final state particles to find Higgs
      do ihep=1,nhep
         if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     #.or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     #.and.(idhep(ihep).eq.25)) then
            j1=ihep
            found=found+1
         endif
      enddo
         
      if(found.lt.1) then
         write(*,*) 'ERROR: Higgs not found'
         call exit(1)
      elseif(found.gt.1) then
         write(*,*) 'ERROR: more Higgs-like particles found'
         call exit(1)
      endif

      
c     Higgs
      
      e1=phep(4,j1)
      px1=phep(1,j1)
      py1=phep(2,j1)
      pz1=phep(3,j1)
      p1=sqrt(px1**2+py1**2+pz1**2)
      pt1=sqrt(px1**2+py1**2)
      y1=log((e1+pz1)/(e1-pz1))/2d0
      m1=sqrt((e1)**2-(px1)**2-(py1)**2-(pz1)**2)
      if (pt1.gt.10d0) then 
         delphi=abs(atan2(py1,px1))
         delphi=min(delphi,2d0*pi-delphi)
      endif 


c     total sigma
      diag=1
      call pwhgfill(diag,1.5d0,dsig/binsize(diag))

c     invariant mass of the Higgs
      diag=diag+1
      call pwhgfill(diag,m1,dsig/binsize(diag))

c     y(H)
      diag=diag+1
      call pwhgfill(diag,y1,dsig/binsize(diag))

c     pt H
      diag=diag+1
      call pwhgfill(diag,pt1,dsig/binsize(diag))

c     pi-delphi
      diag=diag+1
      if ((pi-delphi).gt.1d-13) then
      call pwhgfill(diag,log10(pi-delphi),dsig/binsize(diag))
      endif

c     set up arrays for jet finding
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
      j1=0
      ntracks=0
      njets=0
c     Loop over final state particles to find jets 
      do ihep=1,nhep
         if ((isthep(ihep).eq.1).and.
     1    (((abs(idhep(ihep)).le.10).or.(abs(idhep(ihep)).ge.40))
c     exclude leptons, gauge and higgs bosons
     2    .or.(abs(idhep(ihep)).eq.21)))
c     but  include gluons 
     3           then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'hwanal: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo    

      

      if(buildjets.and.ntracks.gt.0) then
************************************************************************
*     fastkt algorithm
**********************************************************************
c     R = 0.7  Radius parameter
c.....run the clustering 
      R = 0.7d0          
      ptmin_fastkt = 0d0
      call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
     #pjet,njets,jetvec) 
c     
c     ... now we have the jets      


      if (njets.gt.0)then
         ptjet=0d0
         j1=0
         do ijet=1,njets
c............find the hardest jet
            tmp=sqrt(pjet(1,ijet)**2 + pjet(2,ijet)**2)
            if (tmp.gt.ptjet) then
               ptjet=tmp
               j1=ijet
            endif
         enddo
         yjet=log((pjet(4,j1)+pjet(3,j1))
     #/(pjet(4,j1)-pjet(3,j1)))/2d0
         
         do i=1,3
c     pt jet
            diag=diag+1
            if(ptjet.gt.cut(i)) then
               call pwhgfill(diag,ptjet,dsig/binsize(diag))
            endif
c     y jet
            diag=diag+1
            if(ptjet.gt.cut(i)) then
               call pwhgfill(diag,yjet,dsig/binsize(diag))
            endif
c     Dy H-jet
            diag=diag+1
            if(ptjet.gt.cut(i)) then
               call pwhgfill(diag,y1-yjet,dsig/binsize(diag))
            endif
         enddo
      endif
      endif
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
