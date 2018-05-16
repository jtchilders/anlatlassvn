      subroutine setreal(p,rflav,amp2)
      implicit none
      real*8 dotp
      external dotp
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real *8 p(0:3,nlegreal),amp2
      real *8 pp(0:3,nlegreal)
      integer rflav(nlegreal)

      amp2=0d0
      pp=p

      if(rflav(4).lt.0)then
         call swap_momenta(pp(:,4),pp(:,5))
      endif

      if (rflav(1).gt.0.and.rflav(2).eq.0) then
         call swap_momenta_minus(pp(:,2),pp(:,6))

      elseif (rflav(1).eq.0.and.rflav(2).lt.0) then
         call swap_momenta_minus(pp(:,1),pp(:,6))

      elseif (rflav(1).lt.0.and.rflav(2).gt.0) then
         call swap_momenta(pp(:,1),pp(:,2))

      elseif (rflav(1).lt.0.and.rflav(2).eq.0) then
         call swap_momenta_minus(pp(:,2),pp(:,6))
         call swap_momenta(pp(:,1),pp(:,2))

      elseif (rflav(1).eq.0.and.rflav(2).gt.0) then
         call swap_momenta_minus(pp(:,1),pp(:,6))
         call swap_momenta(pp(:,1),pp(:,2))
      endif

      call real_sqamp(pp,rflav,amp2)
      end


      subroutine swap_momenta(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3),tmp(0:3)

      tmp=p1
      p1=p2
      p2=tmp
      end
c
      subroutine swap_momenta_minus(p1,p2)
      implicit none
      real * 8 p1(0:3),p2(0:3)

      call swap_momenta(p1,p2)
      p1=-p1
      p2=-p2
      end
c


c     computes the total squared amplitude of the real diagram,
c     with all the coefficients
      subroutine real_sqamp(pp,rflav,amp2)
      implicit none
      external dotp
      real* 8 dotp
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 amp2,pp(0:3,nlegreal)
      integer rflav(nlegreal),i,j,n
      real * 8 p12,p13,p14,p16,p23,p24,p26,p33,p34,p36,p46,t0,props
      real * 8 gw,colfac,couplz,AL,VL,AQ,VQ,chargeQ,chargeL,T3L,T3Q
      real * 8 s1,s2,s3,s4,s5,s6,s7,s8,s9,s10

      n=3
      gw=ph_unit_e/ph_sthw

      couplz = ph_unit_e/(2*ph_sthw*ph_cthw)
c
c     vectorial and axial couplings to Z boson
      if (mod(abs(rflav(4)),2).eq.1) then
c     LEPTON
         chargeL = -1
         T3L = -1d0/2d0
      elseif (mod(abs(rflav(4)),2).eq.0) then
c     NEUTRINO
         chargeL = 0
         T3L = 1d0/2d0
      endif

      if (rflav(6).eq.0) then
         if (mod(abs(rflav(1)),2).eq.0) then
c     UP TYPE QUARK
            chargeQ = 2d0/3d0
            T3Q = 1d0/2d0
         elseif (mod(abs(rflav(1)),2).eq.1) then
c     DOWN TYPE QUARK
            chargeQ = -1d0/3d0
            T3Q = -1d0/2d0
         endif
      else
         if (mod(abs(rflav(6)),2).eq.0) then
c     UP TYPE QUARK
            chargeQ = 2d0/3d0
            T3Q = 1d0/2d0
         elseif (mod(abs(rflav(6)),2).eq.1) then
c     DOWN TYPE QUARK
            chargeQ = -1d0/3d0
            T3Q = -1d0/2d0
         endif
      endif

      VL = T3L - 2*chargeL*ph_sthw**2
      AL = -T3L
      VQ = T3Q - 2*chargeQ*ph_sthw**2
      AQ = -T3Q

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

      s3 = -16/p26/p16
      s6 = 6*p26*p14**2-4*p12**2*p24+2*p12*p14*p46-6*p24*p16**2+2*p12**2
     #*p46-2*p24*p26*p14-p14*p16*p33+4*p24*p26*p12-p26*p14*p33-2*p26*p14
     #*p36-2*p26*p14*p34-2*p14*p16*p34+2*p14*p16*p13-2*p12*p24*p36-2*p12
     #*p24*p46+2*p12*p14*p36+2*p46*p16*p14-2*p12*p14*p16+p26*p46*p33+2*p
     #26*p46*p36+2*p26*p46*p34-2*p26*p46*p13-4*p26*p46*p14-2*p26*p46*p16
     #+2*p24*p16*p33+4*p24*p16*p36+4*p24*p16*p46
      s5 = 4*p24*p16*p34-6*p24*p16*p13-6*p24*p16*p14+8*p24*p16*p12+6*p26
     #*p14*p13+4*p26*p14*p16+2*p14**2*p16+2*p26*p46**2-2*p24*p26**2+s6+2
     #*p46*p16**2-4*p12*p14**2-2*p12*p46**2-2*p12*p14*p26+2*p46*p16*p13-
     #2*p46*p16*p12+4*p12*p24*p13+4*p12*p24*p14+2*p24*p26*p36+2*p24*p26*
     #p46-4*p24*p26*p16+2*p12*p14*p33+4*p12*p14*p34-4*p12*p14*p13-p12*p4
     #6*p33-2*p12*p46*p36-2*p12*p46*p34-2*p24*p26*p13
      s6 = VQ**2*AL**2
      s4 = s5*s6
      s2 = s3*s4
      s4 = -16/p26/p16
      s7 = 6*p26*p14**2-4*p12**2*p24+2*p12*p14*p46-6*p24*p16**2+2*p12**2
     #*p46-2*p24*p26*p14-p14*p16*p33+4*p24*p26*p12-p26*p14*p33-2*p26*p14
     #*p36-2*p26*p14*p34-2*p14*p16*p34+2*p14*p16*p13-2*p12*p24*p36-2*p12
     #*p24*p46+2*p12*p14*p36+2*p46*p16*p14-2*p12*p14*p16+p26*p46*p33+2*p
     #26*p46*p36+2*p26*p46*p34-2*p26*p46*p13-4*p26*p46*p14-2*p26*p46*p16
     #+2*p24*p16*p33+4*p24*p16*p36+4*p24*p16*p46
      s6 = 4*p24*p16*p34-6*p24*p16*p13-6*p24*p16*p14+8*p24*p16*p12+6*p26
     #*p14*p13+4*p26*p14*p16+2*p14**2*p16+2*p26*p46**2-2*p24*p26**2+s7+2
     #*p46*p16**2-4*p12*p14**2-2*p12*p46**2-2*p12*p14*p26+2*p46*p16*p13-
     #2*p46*p16*p12+4*p12*p24*p13+4*p12*p24*p14+2*p24*p26*p36+2*p24*p26*
     #p46-4*p24*p26*p16+2*p12*p14*p33+4*p12*p14*p34-4*p12*p14*p13-p12*p4
     #6*p33-2*p12*p46*p36-2*p12*p46*p34-2*p24*p26*p13
      s7 = AQ**2*VL**2
      s5 = s6*s7
      s3 = s4*s5
      s1 = s2+s3
      s3 = s1
      s6 = 64/p26
      s8 = 1/p16
      s10 = 2*p26*p14**2+4*p12*p46*p13+4*p12**2*p24+6*p12*p14*p46+2*p24*
     #p16**2-2*p12**2*p46+2*p24*p26*p14-p14*p16*p33-4*p24*p26*p12-p26*p1
     #4*p33-2*p26*p14*p36-2*p26*p14*p34-2*p14*p16*p34+2*p14*p16*p13+2*p1
     #2*p24*p36+2*p12*p24*p46+2*p12*p14*p36-2*p46*p16*p14-2*p12*p14*p16+
     #p26*p46*p33+2*p26*p46*p36+2*p26*p46*p34-2*p26*p46*p13-4*p26*p46*p1
     #4-2*p26*p46*p16
      s9 = 2*p24*p16*p13+2*p24*p16*p14-8*p24*p16*p12+2*p26*p14*p13+2*p14
     #**2*p16+2*p26*p46**2+2*p24*p26**2-2*p46*p16**2-4*p12*p14**2-2*p12*
     #p46**2+s10+2*p12*p14*p26-2*p46*p16*p13+6*p46*p16*p12-4*p12*p24*p13
     #-4*p12*p24*p14-2*p24*p26*p36-2*p24*p26*p46+4*p24*p26*p16+2*p12*p14
     #*p33+4*p12*p14*p34-4*p12*p14*p13-p12*p46*p33-2*p12*p46*p36-2*p12*p
     #46*p34+2*p24*p26*p13
      s7 = s8*s9
      s5 = s6*s7
      s6 = VL*AL*VQ*AQ
      s4 = s5*s6
      s2 = s3+s4
      s3 = s2
      s6 = -16/p26/p16
      s9 = 6*p26*p14**2-4*p12**2*p24+2*p12*p14*p46-6*p24*p16**2+2*p12**2
     #*p46-2*p24*p26*p14-p14*p16*p33+4*p24*p26*p12-p26*p14*p33-2*p26*p14
     #*p36-2*p26*p14*p34-2*p14*p16*p34+2*p14*p16*p13-2*p12*p24*p36-2*p12
     #*p24*p46+2*p12*p14*p36+2*p46*p16*p14-2*p12*p14*p16+p26*p46*p33+2*p
     #26*p46*p36+2*p26*p46*p34-2*p26*p46*p13-4*p26*p46*p14-2*p26*p46*p16
     #+2*p24*p16*p33+4*p24*p16*p36+4*p24*p16*p46
      s8 = 4*p24*p16*p34-6*p24*p16*p13-6*p24*p16*p14+8*p24*p16*p12+6*p26
     #*p14*p13+4*p26*p14*p16+2*p14**2*p16+2*p26*p46**2-2*p24*p26**2+2*p4
     #6*p16**2-4*p12*p14**2-2*p12*p46**2+s9-2*p12*p14*p26+2*p46*p16*p13-
     #2*p46*p16*p12+4*p12*p24*p13+4*p12*p24*p14+2*p24*p26*p36+2*p24*p26*
     #p46-4*p24*p26*p16+2*p12*p14*p33+4*p12*p14*p34-4*p12*p14*p13-p12*p4
     #6*p33-2*p12*p46*p36-2*p12*p46*p34-2*p24*p26*p13
      s9 = VQ**2*VL**2
      s7 = s8*s9
      s5 = s6*s7
      s7 = -16/p26/p16
      s10 = 6*p26*p14**2-4*p12**2*p24+2*p12*p14*p46-6*p24*p16**2+2*p12**
     #2*p46-2*p24*p26*p14-p14*p16*p33+4*p24*p26*p12-p26*p14*p33-2*p26*p1
     #4*p36-2*p26*p14*p34-2*p14*p16*p34+2*p14*p16*p13-2*p12*p24*p36-2*p1
     #2*p24*p46+2*p12*p14*p36+2*p46*p16*p14-2*p12*p14*p16+p26*p46*p33+2*
     #p26*p46*p36+2*p26*p46*p34-2*p26*p46*p13-4*p26*p46*p14-2*p26*p46*p1
     #6+2*p24*p16*p33+4*p24*p16*p36+4*p24*p16*p46
      s9 = 4*p24*p16*p34-6*p24*p16*p13-6*p24*p16*p14+8*p24*p16*p12+6*p26
     #*p14*p13+4*p26*p14*p16+2*p14**2*p16+2*p26*p46**2-2*p24*p26**2+2*p4
     #6*p16**2-4*p12*p14**2-2*p12*p46**2+s10-2*p12*p14*p26+2*p46*p16*p13
     #-2*p46*p16*p12+4*p12*p24*p13+4*p12*p24*p14+2*p24*p26*p36+2*p24*p26
     #*p46-4*p24*p26*p16+2*p12*p14*p33+4*p12*p14*p34-4*p12*p14*p13-p12*p
     #46*p33-2*p12*p46*p36-2*p12*p46*p34-2*p24*p26*p13
      s10 = AQ**2*AL**2
      s8 = s9*s10
      s6 = s7*s8
      s4 = s5+s6
      t0 = s3+s4

c     COLOUR FACTORS AND QUARK-GLUON SWITCH
c     colour factors: CF*n from sum over initial colours, 1/4 from
c     average over initial spins, 1/n from average over quark colours
c     and 1/(n^2-1) from average over gluon colours
      if(.not.rflav(6).eq.0) then
         amp2=-1d0   !quark-gluon switch
         colfac=CF*n/4d0/n/(n**2-1d0)
      else
         amp2=1d0   !no quark-gluon switch
         colfac=CF*n/4d0/n**2
      endif

      amp2=amp2*colfac*t0

c     coupling constants and Z mass
c     factor couplz^2 from each weak vertex: two vertices, couplz^4
c     factor 4mz^4/v^2 from Higgs vertex; but v^2 = 4mw^2/gw^2
c     from Higgs vertex: mz^4 gw^2/mw^2
c     from strong vertex: 4*pi*st_alpha, strip off st_alpha/(2*Pi)
      amp2=amp2*(couplz)**2*(couplz)**2*(ph_Zmass2*gw)**2/ph_Wmass2
     $     *(8*Pi**2)

c     Z propagators

      props = 1/((-2*p26-2*p16+2*p12-ph_Zmass2)**2+ph_ZmZw**2)/((-2*p46-
     #2*p34+2*p14+2*p24-ph_Zmass2)**2+ph_ZmZw**2)

      amp2=amp2*props

      end
