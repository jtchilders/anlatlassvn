      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real*8 p(0:3,nlegreal),amp2
      integer rflav(nlegreal)
      real*8 pp(0:3,nlegreal)
      real*8 dotp
      external dotp

      amp2=0d0
      pp=p

c     if it's W- PRODUCTION then pp=p, no exchange p4 <--> p5!!
      if(rflav(4).lt.0)then
c     it's  W+ PRODUCTION!!
         call swap_momenta(pp(:,4),pp(:,5))
      endif

      if (rflav(1).gt.0.and.rflav(2).eq.0) then
c     quark emission, rflav(1)>0
         call swap_momenta_minus(pp(:,2),pp(:,6))

      elseif (rflav(1).eq.0.and.rflav(2).lt.0) then
c     antiquark emission, rflav(1)>0
         call swap_momenta_minus(pp(:,1),pp(:,6))

      elseif (rflav(1).lt.0.and.rflav(2).gt.0) then
c     gluon emission, rflav(1)<0
         call swap_momenta(pp(:,1),pp(:,2))

      elseif (rflav(1).lt.0.and.rflav(2).eq.0) then
c     antiquark emission, rflav(1)<0
         call swap_momenta_minus(pp(:,2),pp(:,6))
         call swap_momenta(pp(:,1),pp(:,2))

      elseif (rflav(1).eq.0.and.rflav(2).gt.0) then
c     quark emission, rflav(1)<0
         call swap_momenta_minus(pp(:,1),pp(:,6))
         call swap_momenta(pp(:,1),pp(:,2))
      endif      

      call amplitude(pp,rflav,amp2)
      end

      subroutine swap_momenta(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),tmp(0:3)
      tmp=p1
      p1=p2
      p2=tmp
      end


      subroutine swap_momenta_minus(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3)
      call swap_momenta(p1,p2)
      p1=-p1
      p2=-p2
      end


c     computes the total squared amplitude of the real diagram,
c     with all the coefficients
      subroutine amplitude(pp,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 amp2,pp(0:3,nlegreal)
      integer rflav(nlegreal)
      integer i,j,n
      real * 8 p12,p13,p14,p16,p23,p24,p26,p33,p34,p36,p46,t0,props,t1
      real * 8 gw,mw,colfac
      external dotp
      real * 8 dotp

      n=3d0
      gw=ph_unit_e/ph_sthw
      mw=ph_Wmass

      p12=dotp(pp(0,1),pp(0,2))
      p13=dotp(pp(0,1),pp(0,3))
      p14=dotp(pp(0,1),pp(0,4))
      p23=dotp(pp(0,2),pp(0,3))
      p24=dotp(pp(0,2),pp(0,4))
      p33=dotp(pp(0,3),pp(0,3))
      p34=dotp(pp(0,3),pp(0,4))
      p16=dotp(pp(0,1),pp(0,6))
      p26=dotp(pp(0,2),pp(0,6))
      p36=dotp(pp(0,3),pp(0,6))
      p46=dotp(pp(0,4),pp(0,6))

      t1 = (-4*p24+2*p46)/p26*p16-4*p24+2*p14+(-4*p46*p12+2*p24*p34+2*p2
     #4*p46+2*p24*p36-4*p24*p13+2*p14*p46+8*p24*p12+p24*p33-4*p24*p14+2*
     #p46*p13)/p26-2*p24/p16*p26+(-2*p24*p14+4*p24*p12+2*p14*p13+2*p24*p
     #46-2*p24*p13+2*p14**2+2*p24*p36-2*p14*p12)/p16+(-2*p12*p46*p13-2*p
     #12*p24*p36+4*p12*p24*p14-2*p12*p14*p46-2*p12*p24*p46+2*p46*p12**2+
     #4*p12*p24*p13-4*p24*p12**2)/p16/p26

      t0 = -128*t1

c     COLOUR FACTORS, QUARK-GLUON SWITCH AND CKM MATRIX
c     colour factors: CF*n from sum over initial colours, 1/4 from
c     average over initial spins, 1/n from average over quark colours
c     and 1/(n^2-1) from average over gluon colours
      if(.not.rflav(6).eq.0) then
         amp2=-1d0              !quark-gluon switch
         colfac=CF*n/4d0/n/(n**2-1d0)
         j=rflav(6)
         if (rflav(1).eq.0) then
            i=rflav(2)
         else
            i=rflav(1)
         endif
      else
         amp2=1d0               !no quark-gluon switch
         colfac=CF*n/4d0/n**2
         i=rflav(1)
         j=rflav(2)
      endif
      if (mod(i,2).eq.0) then
         amp2=amp2*ph_CKM(abs(i)/2,(abs(j)+1)/2)**2
      else
         amp2=amp2*ph_CKM((abs(i)+1)/2,abs(j)/2)**2
      endif

      amp2=amp2*colfac*t0

c     coupling constants and W mass
c     factor gw^2/8 from each weak vertex: two vertices, gw^4/64
c     factor 4mw^4/v^2 from Higgs vertex; but v^2 = 4mw^2/gw^2
c     from Higgs vertex: mw^2 gw^2
c     factor 4*pi*alpha_s from strong factor: strip off alpha_s/2pi
c     from strong vertex: 8 pi^2
      amp2=amp2 * (gw**2/8d0) * (gw**2/8d0) * (mw**2*gw**2) *(8d0*Pi**2)
c
c     W propagators
      props = 1/((2*p12-2*p16-2*p26-mw**2)**2+ph_Wwidth**2*mw**2)/((-2*p
     #46-2*p34+2*p14+2*p24-mw**2)**2+ph_Wwidth**2*mw**2)
      amp2=amp2*props
      end


c$$$
c$$$
c$$$      subroutine realcolour_lh
c$$$c Wrapper subroutine to call the MadGraph code to associate
c$$$c a (leading) color structure to an event.
c$$$      implicit none
c$$$      include 'nlegborn.h'
c$$$      include 'LesHouches.h'
c$$$      integer rflav(nlegreal),color(2,nlegreal)
c$$$      integer i,j
c$$$      do i=1,nlegreal
c$$$         rflav(i)=idup(i)
c$$$         if (rflav(i).eq.21) rflav(i)=0
c$$$      enddo
c$$$      call real_color(rflav,color)
c$$$      do i=1,2
c$$$         do j=1,nlegreal
c$$$            icolup(i,j)=color(i,j)
c$$$         enddo
c$$$      enddo
c$$$      end
