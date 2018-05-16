
      function suppfact2e(pin,flav)
      implicit none
      include 'nlegborn.h'
      include 'vvsettings.f'
      real * 8  suppfact2e
      real * 8 pin(0:3,nlegborn)
      integer flav(nlegborn)
      include 'PhysPars.h'
      real * 8 pcombz1(0:3),pcombz2(0:3),mz1,mz2
      real * 8 pcombw1(0:3),pcombw2(0:3),mw1,mw2
      real * 8 f1,f2
c     identical lepton 1, identical lepton 2, different lepton
      integer nu,l1,l2,dl
      real * 8 dotp
      external dotp
      if(.not.interference) then
         suppfact2e=1
         return
      endif
      if(2*(flav(3)/2).eq.flav(3)) then
         nu=3
         l1=4
      else
         nu=4
         l1=3
      endif
      if(abs(flav(nu)+flav(l1)).ne.1) then
         write(*,*) 'suppfact2e: error: wrong flavour structure'
         call exit(-1)
      endif
      if(flav(5).eq.flav(l1)) then
         l2=5
         dl=6
      elseif(flav(6).eq.flav(l1)) then
         l2=6
         dl=5
      else
         suppfact2e=1
         return
      endif
      if(flav(l1).ne.flav(l2).or.flav(dl)+flav(l1).ne.0) then
         write(*,*) 'suppfact2e: error: wrong flavour structure'
      endif
      pcombz1=pin(:,l1)+pin(:,dl)
      pcombw1=pin(:,nu)+pin(:,l2)
      pcombz2=pin(:,l2)+pin(:,dl)
      pcombw2=pin(:,nu)+pin(:,l1)
      mz1=dotp(pcombz1,pcombz1)
      mz2=dotp(pcombz2,pcombz2)
      mw1=dotp(pcombw1,pcombw1)
      mw2=dotp(pcombw2,pcombw2)
      f1=mz1*((mz1-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
     1      *((mw1-ph_Wmass2)**2+ph_Wwidth**2*ph_Wmass2)
      f2=mz2*((mz2-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
     1      *((mw2-ph_Wmass2)**2+ph_Wwidth**2*ph_Wmass2)
      suppfact2e=2*f1**2/(f1**2+f2**2)
      end
