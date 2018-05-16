c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none
      end 

      
      subroutine analysis(dsig)
      implicit none
      real * 8 dsig
      include 'hepevt.h'
      end



      subroutine buildjets(process,njets,pjet,jetvec)
c     arrays to reconstruct jets
      implicit none
      include 'hepevt.h'
      integer njets
      real * 8 pjet(4,*)
      character * 20 process
      integer maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real *8 ptrack(4,maxtrack)
      real *8 r,ptmin_fastkt
      integer jetvec(maxtrack),ihep_of_track(maxtrack)
      integer ihep,j1,ntracks,jpart,jjet,mu
      real * 8 found
      real * 8 f,caf,sf
      integer seed
      data seed/1/
      save seed
      common/cptrack/ihep_of_track,ptrack,ntracks
      common/cptmin_fastkt/ptmin_fastkt
      logical ini
      save ini
      data ini/.true./

      if (ini) then
         write(*,*) '*******************************'
         write(*,*) '   ptmin_fastkt = ',ptmin_fastkt
         write(*,*) '*******************************'
         ini=.false.
      endif
c     get valid tracks
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
      found=0
      ntracks=0
      njets=0
c     loop over final state particles to find jets 
      do ihep=1,nhep
c     "stable" particles (particles that we detect)
         if (isthep(ihep).eq.1) then
            if(ntracks.eq.maxtrack) then
               write(*,*)
     #              'analyze: too many particles, increase maxtrack'
               stop
            endif
c     copy momenta to construct jets 
            ntracks=ntracks+1
            ihep_of_track(ntracks)=ihep
            do mu=1,4
               ptrack(mu,ntracks)=phep(mu,ihep)
            enddo
         endif
      enddo
      if (ntracks.eq.0) then
         njets=0
         return
      endif
c     siscone algorithm
c*********************************************************************
c      R = 0.7  radius parameter
c      f = 0.5  overlapping fraction
c.....run the clustering
c      call fastjetsiscone(ptrack,ntracks,0.7d0,0.5d0,pjet,njets) 
c*********************************************************************
c     fastkt algorithm
c*********************************************************************
c      R = 0.7  Radius parameter
c.....run the clustering 
c      R = 0.5d0          
c      ptmin_fastkt = 0d0
c      call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
c     #     pjet,njets,jetvec)
c     now we have the jets

c      call fastjetd0runiicone(ptrack,ntracks,0.7d0,6d0,0.5d0,pjet,njets)
c      return

c     kt algo
c      R=0.4
c      f=0.75
c      sf=1
c      caf=1
c      call fastjetcdfmidpoint(ptrack,ntracks,r,f,sf,caf,pjet,njets)

c      if     (process.eq."antikt R04") then
c         R=0.4d0 
c      elseif (process.eq."antikt R06") then
c         R=0.6d0
c      elseif (process.eq."antikt R02") then
c         R=0.2d0
c      elseif (process.eq."antikt R08") then
c         R=0.8d0
c      else
c         write(*,*) 'JET ANALYSIS TO USE UNKNOWN:',process
c         call exit(1)
c      endif
      R=0.4d0
      if (process.eq."antikt") then
         call fastjetantikt(ptrack,ntracks,ptmin_fastkt,R,
     $        pjet,njets,jetvec)
      elseif (process.eq."kt") then
         call fastjetktwhich(ptrack,ntracks,ptmin_fastkt,R,
     $        pjet,njets,jetvec)
      else
         write(*,*) 'JET ANALYSIS TO USE UNKNOWN:',process
         call exit(1)
      endif
      end





      function is_B_hadron(id)
      implicit none
      logical is_B_hadron
      integer id
      is_B_hadron=((id.gt.-600).and.(id.lt.-500)).or.
     $     ((id.gt.5000).and.(id.lt.6000)).or.(id.eq.5)
      end

      function is_BBAR_hadron(id)
      implicit none
      logical is_BBAR_hadron
      integer id
      is_BBAR_hadron=((id.gt.500).and.(id.lt.600)).or.
     $     ((id.gt.-6000).and.(id.lt.-5000)).or.(id.eq.-5)
      end


      
