      subroutine analysis_driver(dsig0,ikin)
      implicit none
      real * 8 dsig0
      integer ikin
      integer jpart,mu
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'hepevt.h'
      real * 8 powheginput
      logical testplots,ini
      data testplots/.true./
      data ini/.true./
      save testplots,ini
      if(ini) then
         if (powheginput('#testplots').eq.1d0) then
            testplots=.true.
         else
            testplots=.false.
         endif
         ini=.false.
      endif
      if(.not.testplots) return
      if(ikin.eq.0) then
         nhep=nlegborn
         do jpart=1,nhep
            do mu=1,3
               phep(mu,jpart)=kn_pborn(mu,jpart)
            enddo
            phep(4,jpart)=kn_pborn(0,jpart)
            phep(5,jpart)=sqrt(abs(phep(4,jpart)**2-
     #           phep(1,jpart)**2-phep(2,jpart)**2-phep(3,jpart)**2))
            if(jpart.le.2) then
               isthep(jpart)=-1
            else
               isthep(jpart)=1
            endif
c we may not know the flavour af all partons
            idhep(jpart)=0
         enddo
      else
         nhep=nlegreal
         do jpart=1,nhep
            do mu=1,3
               phep(mu,jpart)=kn_preal(mu,jpart)
            enddo
            phep(4,jpart)=kn_preal(0,jpart)
            phep(5,jpart)=sqrt(abs(phep(4,jpart)**2-
     #           phep(1,jpart)**2-phep(2,jpart)**2-phep(3,jpart)**2))
            if(jpart.le.2) then
               isthep(jpart)=-1
            else
               isthep(jpart)=1
            endif
c we may not know the flavour af all partons
            idhep(jpart)=0
         enddo
      endif
c     assign colorless particles'ID and massive partons ID
c     NB: all regions/flavstruct have the same white and colored particles
c     We can then use the flst_born(jpart,1) of the FIRST Born
      do jpart=3,flst_lightpart-1
         idhep(jpart)=flst_born(jpart,1)
      enddo
c     call analysis routine
      call analysis(dsig0)
      end


     
      subroutine lhtohep
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      integer j
      nhep=nup
      do j=1,nup
         idhep(j)=idup(j)
         isthep(j)=istup(j)
         phep(:,j)=pup(:,j)
         jmohep(:,j)=mothup(:,j)
      enddo
      end

