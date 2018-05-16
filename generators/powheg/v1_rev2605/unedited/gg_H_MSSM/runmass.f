c$$$      double precision function asf(x,tc,asc,zm)
c$$$
c$$$c     compute the running (SM, MSbar) alpha_s
c$$$      
c$$$      implicit double precision (a-h,o-z)
c$$$
c$$$      pi = 4d0*atan(1d0)
c$$$
c$$$      fn=5d0
c$$$      
c$$$      b0=11d0-2*fn/3d0
c$$$      b1=102d0-38d0*fn/3d0
c$$$      vvv=1-b0*asc/(2d0*pi)*log(zm/x)
c$$$      
c$$$      if(x.le.tc) then          ! five flavors
c$$$
c$$$         asf=asc/vvv*(1-b1/b0*asc/(4d0*pi*vvv)*log(vvv))
c$$$         
c$$$      else                      
c$$$
c$$$         vvv=1-b0*asc/(2d0*pi)*log(zm/tc) ! first evolve up to q=mt
c$$$         
c$$$         ast=asc/vvv*(1-b1/b0*asc/(4*pi*vvv)*log(vvv))
c$$$         
c$$$         b0t=b0-2d0/3d0         ! six flavours
c$$$         b1t=b1-38d0/3d0
c$$$         vvv=1-b0t*ast/(2*pi)*log(tc/x) !     now evolve up to the scale >mt
c$$$
c$$$         asf=ast/vvv*(1-b1t/b0t*ast/(4*pi*vvv)*log(vvv))
c$$$
c$$$      endif
c$$$      
c$$$      return
c$$$      end
c$$$
c$$$c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      double precision function runt(x,tc,asc,zm)

c     compute the running (SM, MSbar) top mass
 
      implicit double precision (a-h,o-z)
      
      pi=4*atan(1d0)
      
      if(x.le.tc)then
         fn=5d0
      else
         fn=6d0
      endif

      b0=11-2*fn/3d0
      b1=102-38*fn/3d0
      g0=8d0
      g1=404d0/3d0-40*fn/9d0
      
      asx=asf(x,tc,asc,zm)
      ast=asf(tc,tc,asc,zm)
      rrr=tc*(asx/ast)**(g0/(2*b0))
      
c     this is the relation between mpole/mtrun(mtrun)
c      pol1=1+4*ast/(3*pi)+8.243d0*(ast/pi)**2

c     this is the relation between mpole/mt(mpole)
      pol2= 1+4*ast/(3*pi) +10.9d0*(ast/pi)**2
      
      corr=1+ast*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/ast-1)
      
      runt=rrr*corr/pol2

      runt = runt*(1d0-asx/pi/3d0-(asx/pi)**2*29d0/72d0) ! shift to DRbar
      
      return
      end

c     
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c
c$$$
c$$$      double precision function runb(x,mbmb,mt,asmz,mz)
c$$$
c$$$c     compute the running (SM, DRbar) bottom mass
c$$$ 
c$$$      implicit none
c$$$      
c$$$      double precision x,mbmb,mt,asmz,mz
c$$$      double precision pi,b0,b1,g0,g1,asx,asb,asf,ast,rrr,corr,mbmt
c$$$      integer nf
c$$$
c$$$      pi=4*atan(1d0)
c$$$      
c$$$      if(x.le.mt) then
c$$$
c$$$         nf = 5                 ! evolve from mb to x with nf=5
c$$$         b0 = 11-2*nf/3d0
c$$$         b1 = 102-38*nf/3d0
c$$$
c$$$         g0 = 8d0
c$$$         g1 = 404d0/3d0-40*nf/9d0
c$$$      
c$$$         asx=asf(x,mt,asmz,mz)
c$$$         asb=asf(mbmb,mt,asmz,mz)
c$$$
c$$$         rrr=mbmb*(asx/asb)**(g0/(2*b0))      
c$$$         corr=1+asb*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/asb-1)
c$$$      
c$$$         runb=rrr*corr
c$$$
c$$$      else
c$$$
c$$$         nf = 5                 ! first evolve from mb to mt with nf=5
c$$$         b0 = 11-2*nf/3d0
c$$$         b1 = 102-38*nf/3d0
c$$$
c$$$         g0 = 8d0
c$$$         g1 = 404d0/3d0-40*nf/9d0
c$$$      
c$$$         ast=asf(mt,mt,asmz,mz)
c$$$         asb=asf(mbmb,mt,asmz,mz)
c$$$
c$$$         rrr=mbmb*(ast/asb)**(g0/(2*b0))      
c$$$         corr=1+asb*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(ast/asb-1)
c$$$      
c$$$         mbmt=rrr*corr
c$$$         
c$$$         nf = 6                 ! then evolve from mt to x with nf=6
c$$$         b0 = 11-2*nf/3d0
c$$$         b1 = 102-38*nf/3d0
c$$$
c$$$         g0 = 8d0
c$$$         g1 = 404d0/3d0-40*nf/9d0
c$$$      
c$$$         asx=asf(x,mt,asmz,mz)
c$$$
c$$$         rrr=mbmt*(asx/ast)**(g0/(2*b0))      
c$$$         corr=1+ast*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/ast-1)
c$$$      
c$$$         runb=rrr*corr
c$$$
c$$$      endif
c$$$
c$$$      runb = runb*(1d0-asx/pi/3d0-(asx/pi)**2*29d0/72d0) ! shift to DRbar
c$$$
c$$$      return
c$$$      end
c$$$
