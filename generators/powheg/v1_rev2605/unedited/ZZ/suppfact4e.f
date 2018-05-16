
      function suppfact4e(pin,flav)
      implicit none
      real * 8  suppfact4e
      real * 8 pin(0:3,*)
      integer flav(*)
      include 'PhysPars.h'
      include 'vvsettings.f'
      real * 8 p34(0:3),p56(0:3),p36(0:3),p45(0:3),
     1         m34,m56,m36,m45,f34,f56,f36,f45
      integer i,j,k,mu,nu,ipow,id,io
      real * 8 dotp,powheginput
      external dotp,powheginput
      logical ini
      data ini/.true./
      save ini
      if(.not.interference) then
         suppfact4e=1
         return
      endif
      if(ini) then
         ipow=powheginput("#suppfactpow")
         if(ipow.lt.0) ipow=2
         write(*,*) ' suppfact4e: ipow=',ipow
         ini=.false.
      endif
      if(flav(3).eq.flav(5)) then
         id=5
         io=6
      elseif(flav(3).eq.flav(6)) then
         id=6
         io=5
      else
         id=0
      endif
      if(id.ne.0) then
         p34=pin(:,3)+pin(:,4)
         p56=pin(:,id)+pin(:,io)
         p36=pin(:,3)+pin(:,io)
         p45=pin(:,4)+pin(:,id)
         m34=dotp(p34,p34)
         m56=dotp(p56,p56)
         m36=dotp(p36,p36)
         m45=dotp(p45,p45)
         
         f34=m34*((m34-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
         f56=m56*((m56-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
         f36=m36*((m36-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
         f45=m45*((m45-ph_Zmass2)**2+ph_Zwidth**2*ph_Zmass2)
         suppfact4e=2*(f36*f45)**ipow/((f36*f45)**ipow+(f34*f56)**ipow)
      else
         suppfact4e=1
      endif
      end
