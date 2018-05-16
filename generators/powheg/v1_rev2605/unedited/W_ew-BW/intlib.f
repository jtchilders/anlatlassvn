c---------------------------------------------------------------------
c     checked version of "integrallib.f"
c---------------------------------------------------------------------
c====================================================================
c     checked with hollik's routines ; status 16.2.93 ; 
c====================================================================
c--------------------------------------------------------------------
      subroutine fint(s,m1,m2,f)
c--------------------------------------------------------------------
c     the f-function as specified in : boehm et.al., 
c     fortschr.phys.34 ; and in : hollik, fortschr.phys. 38 
c--------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-6)
      
      complex*16 f,sc
      include 'pwhg_wzgrad.h'
      sc=s+ieps
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      p = s-(dabs(m1)+dabs(m2))**2
      m = s-(dabs(m1)-dabs(m2))**2
      imagf = 0d0
      rest = 0d0
      if (abs(s).lt.eps) then
         realf = 0d0
         goto 999
      endif
c
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x2+eps) then
               realf = 1d0+(1d0-x2/s)*cdlog(1d0/(sc/x2-1d0))
               imagf = pi*(1d0-x2/s)
            else if (s.lt.x2-eps) then
               realf = 1d0+(1d0-x2/s)*dlog(1d0/(1d0-s/x2))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      else if (x2.lt.eps) then
         if (x1.lt.eps) then
            realf = 0d0
         else
            if (s.gt.x1+eps) then
               realf = 1d0+(1d0-x1/s)*cdlog(1d0/(sc/x1-1d0))
               imagf = pi*(1d0-x1/s)
            else if (s.lt.x1-eps) then
               realf = 1d0+(1d0-x1/s)*dlog(1d0/(1d0-s/x1))
            else
               realf = 1d0
            endif 
         endif
         goto 999
      endif
c
      if ((abs(s).lt.x1/1d3).or.(abs(s).lt.x2/1d3)) then
         if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
            realf = s/6d0/x1
         else
            realf = s/(x1-x2)**2*((x1+x2)/2d0
     $           -x1*x2/(x1-x2)*dlog(x1/x2))
         endif
         goto 999
      endif
c
      if (dabs(dabs(m1)-dabs(m2)).lt.eps) then
         rest = 2d0
      else
         rest = 1d0-((x1-x2)/s-(x1+x2)/(x1-x2))*dlog(x1/x2)/2d0
      endif
c
      if (m.lt.0d0) then
         realf = dsqrt(p*m)*dlog((dsqrt(-p)+dsqrt(-m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
      else if (p.lt.0d0) then
         realf = -2d0*dsqrt(-p*m)*datan(dsqrt(-m/p))/s
      else
         realf = -dsqrt(p*m)*dlog((dsqrt(p)+dsqrt(m))**2
     $        /(4d0*dsqrt(x1*x2)))/s
         imagf = dsqrt(p*m)/s*pi
      endif
 999  continue
      f = dcmplx(rest+realf,imagf)
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine afunc(m,a0)
      implicit real*8(a-z)
      parameter (eps=1d-4)
      common/renorm/mue,mue2,mue4
      muedr=mue
      if (m**2.lt.eps) then 
         a0 = 0d0
      else
         a0 = m**2*(1d0-dlog(m**2/muedr**2))
      endif
      end
c---------------------------------------------------------------------
      subroutine bfunc(s,m1,m2,b0,b1,b20)
c---------------------------------------------------------------------
c     b0 and b1 are the b0"quer" and b1"quer"-quantities defined 
c     in : hollik, fortschr.phys. 38 (but with log-terms incl.)
c     b20 is the g_{mu}.{nu} coefficient of the 
c     tensor-2-point-integral (log-terms incl, too)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      parameter (eps=1d-4)
      complex*16 b0,b1,b20,ff,sc,ieps
      common/cc/ieps
      common/renorm/mue,mue2,mue4
      muedr=mue
      sc=s+ieps
      pi = 4d0*datan(1d0)
      x1 = m1**2
      x2 = m2**2
      if (dabs(s).lt.eps) then
         call bnull(m1,m2,realb0)
         call beins(m1,m2,realb1)
         call bzwnull(m1,m2,realb20)
         b0 = dcmplx(realb0,0d0)
         b1 = dcmplx(realb1,0d0)
         b20 = dcmplx(realb20,0d0)
         goto 999
      endif
      if ((x1.lt.eps).and.(x2.lt.eps)) goto 47 
      call fint(s,m1,m2,ff)
 47   continue
      if (x1.lt.eps) then
         if (x2.lt.eps) then
            b0 = 2d0-cdlog(sc/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x2/muedr**2)
            b1 = -(1d0+(1d0-x2/s)*ff-dlog(x2/muedr**2))/2d0
         endif
         goto 666
      endif
      if (x2.lt.eps) then
         if (x1.lt.eps) then
            b0 = 2d0-cdlog(sc/muedr**2)+dcmplx(0d0,1d0)*pi
            b1 = -b0/2d0
         else
            b0 = 1d0+ff-dlog(x1/muedr**2)
            b1 = -(1d0+(1d0+x1/s)*ff-dlog(x1/muedr**2))/2d0
         endif
         goto 666
      endif
      if (dabs(x1-x2).lt.1d3*eps) then
         b0 = ff
         b1 = -ff/2d0
      else
         b0 = 1d0-(x1+x2)/(x1-x2)*dlog(x1/x2)/2d0+ff
         b1 = -1d0/2d0+x1/(x1-x2)*dlog(x1/x2)/2d0-(s+x1-x2)/2d0/s*ff
      endif
      b0 = b0 - dlog(x1*x2/muedr**4)/2d0
      b1 = b1 + dlog(x2/muedr**2)/2d0 
c     
 666  continue
      call afunc(m2,a0)
      b20 = a0/6d0+x1/3d0*b0
     $     +(s+x1-x2)/6d0*b1+(x1+x2-s/3d0)/6d0
 999  continue
      end
c---------------------------------------------------------------------
      subroutine bderiv(x,m1,m2,db0,db1,db20)
c---------------------------------------------------------------------
c     real parts of d(b0(s,m1,m2))/ds , d(b1(s,m1,m2))/ds
c     and d(b20(s,m1,m2))/ds
c     (s -> x in the code)
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 cf,b0,b1,b20
      parameter (eps=1d-1)
*
      pi = 4d0*datan(1d0)
      xm1 = m1**2
      xm2 = m2**2
*
      if (x.lt.eps) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm2
            db1 = -1d0/12d0/xm2
         else if (xm1.lt.eps) then
            db0 = 1d0/2d0/xm2
            db1 = -1d0/6d0/xm2
         else if (xm2.lt.eps) then
            db0 = 1d0/2d0/xm1
            db1 = -1d0/6d0/xm1
         else
            db0 = (xm1**2-xm2**2-2d0*xm1*xm2*dlog(xm1/xm2))
     $           /2d0/(xm1-xm2)**3
            fss0 = (xm1**3-xm2**3+9d0*xm1*xm2*(xm1-xm2)
     $           -6d0*xm1*xm2*dlog(xm1/xm2)*(xm1+xm2))
     $           /3d0/(xm1-xm2)**5
            db1 = -db0/2d0-(xm1-xm2)/4d0*fss0
         endif
         goto 11
      endif
 13   continue
      if ((x.lt.xm1/1d3).or.(x.lt.xm2/1d3)) then
         if (dabs(xm1-xm2).lt.eps) then
            db0 = 1d0/6d0/xm1
         else
            db0 = ((xm1+xm2)/2d0-xm1*xm2/(xm1-xm2)*dlog(xm1/xm2))
     $           /(xm1-xm2)**2
         endif
         db1 = -db0/2d0
         goto 11
      endif
c
      if ((xm1.lt.eps).and.(xm2.lt.eps)) then
         db0 = -1d0/x
         db1 = -db0/2d0
         goto 11
      endif
c
      if (xm2.lt.eps) then
         if (x.gt.xm1+eps) then
            deriv = -(1d0+xm1/x*dlog(x/xm1-1d0))/x
         else if (x.lt.xm1-eps) then
            deriv = -(1d0+xm1/x*dlog(1d0-x/xm1))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      if (xm1.lt.eps) then
         if (x.gt.xm2+eps) then
            deriv = -(1d0+xm2/x*dlog(x/xm2-1d0))/x
         else if (x.lt.xm2-eps) then
            deriv = -(1d0+xm2/x*dlog(1d0-x/xm2))/x
         else
            deriv = 0d0
         endif 
         goto 10
      endif
c
      sm = xm1+xm2
      dm = xm2-xm1
      sm12 = (dabs(m1)+dabs(m2))**2
      dm12 = (dabs(m1)-dabs(m2))**2
      lm = dlog(xm2/xm1)/2d0
*
      if (x.lt.dm12) then
         s = dsqrt(sm12-x)
         d = dsqrt(dm12-x)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2+2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      else if (x.lt.sm12) then
         if (dabs(xm1-xm2).lt.eps) then
            if (x.lt.eps) then
               deriv = 1d0/6d0/xm1
            else
               sx = 4d0*xm1/x
               deriv = (sx/dsqrt(sx-1d0)
     $              *datan(1d0/dsqrt(sx-1d0))-1d0)/x
            endif
         else
*
c Achtung: fuer x=sm12 geht deriv -> infty! 
*
            s = dsqrt(sm12/x-1d0)
            d = dsqrt(1d0-dm12/x)
            fact = datan(d/s)
            deriv = (dm*lm/x+((d**2-s**2+2d0*d**2*s**2)
     $           /(s*d))*fact-1d0)/x
         endif
      else
         s = dsqrt(x-sm12)
         d = dsqrt(x-dm12)
         fact = dlog((d+s)**2/(4d0*dabs(m1)*dabs(m2)))
         deriv = (dm*lm/x-((d**2+s**2-2d0*d**2*s**2/x)
     $        /(2d0*s*d))*fact-1d0)/x
      endif
*
 10   continue
c--------------------------------------------------------------------
      db0 = deriv
      call fint(x,m1,m2,cf)
      f = dreal(cf)
      db1 = -(x+xm1-xm2)/2d0/x*deriv+(xm1-xm2)/2d0/x**2*f
c     
 11   continue
c
      call bfunc(x,m1,m2,b0,b1,b20)
      db20 = xm1/3d0*db0+dreal(b1)/6d0+(x+xm1-xm2)/6d0*db1
     $     -1d0/18d0
c
 99   continue
      end
c---------------------------------------------------------------------
      subroutine cfuncnew(m,s,m1,m2,m3,c0,c1p,c1m,c20,c2p,c2m)
c---------------------------------------------------------------------
c     contains the scalar coefficient functions of the 
c     tensor three point integral
c---------------------------------------------------------------------
      implicit real*8(a-z)
      complex*16 c0,c1p,c1m,c20,c2p,c2m,C0_
     $     ,b0x32,b0x31,b0s12,b1x32,b1x31,adummy,bdummy
      x = m**2
      x1 = m1**2
      x2 = m2**2
      x3 = m3**2
c
      call bfunc(s,m1,m2,b0s12,adummy,bdummy)
      call bfunc(x,m3,m2,b0x32,b1x32,adummy)
      call bfunc(x,m3,m1,b0x31,b1x31,adummy)
c new C0 function added on Aug.21 2004
c      call c0int(m,s,m1,m2,m3,c0)
      c0=C0_(m**2,s,m**2,m3,m1,m2,0)
c
      c1m = ((b0x32+b0x31)/2d0-b0s12
     $     -(x+x3-(x1+x2)/2d0)*c0)/(4d0*x-s)
      c1p = ((b0x32-b0x31)/2d0+(x1-x2)/2d0*c0)/s
      c20 = b0s12/4d0+x3/2d0*c0-(x1-x2)/4d0*c1p
     $     -(x1+x2-2d0*(x+x3))/4d0*c1m+1d0/4d0
      c2m = (-c20+(b1x32+b1x31)/4d0+b0s12/2d0
     $     -(x+x3-(x1+x2)/2d0)*c1m)/(4d0*x-s)
      c2p = (-c20-(b1x32+b1x31)/4d0+(x1-x2)/2d0*c1p)/s
      end
c---------------------------------------------------------------------
      subroutine c0int(mf,s,m1,m2,m3,c0)
c**************************************************************
c                                                             *
c  the scalar vertex integral with equal external masses mf   *
c                                                             *
c-------------------------------------------------------------------
c     s = momentum transfer; m1,m2,m3  are the internal masses
c
      implicit real*8 (a-y)
      complex*16 z1,z2,z11,z12,z21,z22,cl1,cl2,cl3,cspen,spence,
     &     int,c0
      xmf=mf*mf
      if (mf.lt.1d-1) then
         mfstrich = 1d-1
         xmf = mfstrich**2
      endif
c     xm's : are fermion and boson masses squared
      xm1=m1*m1
      xm2=m2*m2
      xm3=m3*m3
c     the t'hooft-veltman parameters
      a=1.d0
      b=xmf/s
      c=-1.d0
      d=xm1-xm2-s
      e=xm3-xm1-xmf+s
      f=xm2/s
      d=d/s
      e=e/s
c     discriminante for alpha-equation
      disc=c*c-4.d0*a*b
      if (disc .lt. 0.d0) goto 500
      al=(-c-dsqrt(disc))/2.d0/b
      nenner=c+2.d0*al*b
c..the first integral.............................................
      y0=-(d+e*al+2.d0*a+c*al)/nenner
      y01=y0-1.d0
      d1=(c+e)**2-4.d0*b*(a+d+f)
      x1=-(c+e)/2.d0/b
      if (d1.gt.0.d0) goto 10
c.......complex zeroes of logarithms
      sq1=dsqrt(-d1)
      x2=sq1/2.d0/b
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl1=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 15
10    continue
c........real zeroes
      sq1=dsqrt(d1)
      x2=sq1/2.d0/b
      y1=x1+x2
      y2=x1-x2
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      cl1=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
15    continue
c..the second integral............................................
      y0=-(d+e*al)/nenner/(1.d0-al)
      y01=y0-1.d0
      d2=(e+d)**2-4.d0*f*(a+b+c)
      x1=-(e+d)/2.d0/(a+b+c)
      if(d2.gt.0.d0) goto 20
c.......complex zeroes of logarithms
      sq2=dsqrt(-d2)
      x2=sq2/2.d0/(a+b+c)
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl2=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 25
20    continue
c........real zeroes
      x2=dsqrt(d2)/2.d0/(a+b+c)
      y1=x1+x2
      y2=x1-x2
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl2=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
25    continue
c..the third integral............................................
      y0=(d+e*al)/nenner/al
      y01=y0-1.d0
      d3=d*d-4.d0*a*f
      x1=-d/2.d0/a
      if (d3.gt.0.d0) goto 30
c........complex zeroes of logarithms
      sq3=dsqrt(-d3)
      x2=sq3/2.d0/a
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl3=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 35
30    continue
c........real zeroes
      x2=dsqrt(d3)/2.d0/a
      y1=x1+x2
      y2=x1-x2
 31   format(1h ,3e12.4)
      y11=y0 /(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl3=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
35    continue
c..summation of the 3 integrals ....................................
      int = -cl1+cl2-cl3
      c0 = int/nenner/s
      goto 501
500   continue
c..error message for complex alpha................................
      write(6,21)
21    format(1h ,'  i cannot handle a complex alpha')
501   return
      end
c--------------------------------------------------------------------
c     the x"null" subroutines calculate the b0-, b1-, c0-,
c     d0- and d20-integrals at zero external momenta  
c--------------------------------------------------------------------
      subroutine bnull(a,b,b0)
      implicit real*8(a-z)
      parameter (eps=1d-8)      
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b0-function !'
         b0 = 0d0
         goto 2
      endif
      if (xa*xb.eq.0d0) then
         zw = 1d0-dlog(x/muedr**2)
      else
         zw = -dlog(xa*xb/muedr**4)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw+1d0-(xa+xb)/(xa-xb)*dlog(xa/xb)/2d0
         endif
      endif
      b0 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine beins(a,b,b1)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      x = xa+xb
      if (x.lt.eps) then
c        write(6,*)'all args zero is not allowed for b1-function !'
         b1 = 0d0
         goto 2
      endif      
      if (xa.eq.0d0) then
         zw = -(1d0/2d0-dlog(xb/muedr**2))/2d0
      else if (xb.eq.0d0) then
         zw = -(3d0/2d0-dlog(xa/muedr**2))/2d0
      else
         zw = dlog(xb/muedr**2)/2d0
         if (dabs(dabs(a)-dabs(b)).gt.eps) then
            zw = zw -(xa+xb)/(xa-xb)/4d0-1d0/2d0
     $           +xa/(xa-xb)*dlog(xa/xb)/2d0
     $           +xa*xb/(xa-xb)**2*dlog(xa/xb)/2d0
         endif
      endif
      b1 = zw
 2    continue
      end
c-----------------------------------------------------------------     
      subroutine bzwnull(a,b,b20)
      implicit real*8(a-z)
      common/renorm/mue,mue2,mue4
      muedr=mue
      xa = a**2
      xb = b**2
      call bnull(a,b,b0)
      call beins(a,b,b1)
      zw = xa*b0/3d0+(xa-xb)*b1/6d0+(xa+xb)/6d0
      if (xb.ne.0d0) then 
         zw = zw+xb*(1d0-dlog(xb/muedr**2))/6d0
      endif
      b20 = zw
      end
c-----------------------------------------------------------------
      subroutine cnull(a,b,c,c0)
      implicit real*8(a-z)
      parameter (eps=1d-8)
      xa = a**2
      xb = b**2
      xc = c**2
      if (dabs(dabs(a)-dabs(b)).lt.eps) then
         if (dabs(dabs(a)-dabs(c)).lt.eps) then
            zw = -1d0/2d0/xa
         else 
            zw = (-1d0+xc/(xa-xc)*dlog(xa/xc))/(xa-xc)
         endif
      else 
         call bnull(a,c,bac)
         call bnull(b,c,bbc)
         zw = (bac-bbc)/(xa-xb)
      endif
      c0 = zw
      end
c-----------------------------------------------------------------     
      subroutine dnull(a,b,c,d,d0)
      implicit real*8(a-z)
      integer i,j
      real*8 m(4)
      parameter (eps=1d-8)
      m(1) = a
      m(2) = b
      m(3) = c
      m(4) = d
      do 30 i = 1,4
         do 40 j = i+1,4
            if (dabs(dabs(m(i))-dabs(m(j))).lt.eps) then
               m(i) = m(i)+eps
               m(j) = m(j)-eps
            endif
 40      continue
 30   continue
      call cnull(a,b,c,c0abc)
      call cnull(a,b,d,c0abd)
      zw = (c0abc-c0abd)/(c**2-d**2)
      d0 = zw
      end
c-----------------------------------------------------------------
      subroutine dzwnull(a,b,c,d,d20)
      implicit real*8(a-z)
      call cnull(b,c,d,c0bcd)
      call dnull(a,b,c,d,d0abcd)
      zw = c0bcd+a**2*d0abcd
      d20 = zw
      end
c-----------------------------------------------------------------     
************************************************************************
        function C0_(q01,q12,q20,m0,m1,m2,ext)                            
************************************************************************
*	scalar 3-point function
*-----------------------------------------------------------------------
*	General result from A.Denner, Fortschr. Phys. 41 (1993) 307
* d0=q^2-m0^2, d1=(q+p1)^2-m1^2, d2=(q+p2)^2-m2^2
* q01=p1^2, q20=p2^2, q12=(p1-p2)^2
*-----------------------------------------------------------------------
************************************************************************
        implicit real*8 (a-z)                                         
	integer i,j,ext
	complex*16 c0_,cspens,etass,ieps,alpha,alp(0:2),sc,xs
	complex*16 y0(0:2),y(0:2,-1:1),x(0:2,-1:1)
	real*8 thp(0:2),thy(0:2)

	kappa2(a,b,c) = a**2+b**2+c**2-2d0*a*b-2d0*a*c-2d0*b*c
        lambda2=lambda*lambda
        pi   = 4d0*datan(1d0)                                               
	eps  = 1d-17
	ieps = dcmplx(0d0,eps)
	m02 = m0**2
	m12 = m1**2
	m22 = m2**2
	p01 = q01
	p12 = q12
	p20 = q20
C*** Check for IR divergence
c    Permutate until m0=0
10	continue
	if ((m02*m12*m22.eq.0d0).and.(m02.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 20
	endif
	goto 10  
20	continue
	if ((p01.eq.m12).and.(p20.eq.m22).and.(m02.eq.0d0)) goto 500

C****** Regular C0 function
c    Permutate until p01=0
11	continue
	if ((p01*p12*p20.eq.0d0).and.(p01.ne.0d0)) then
	  m   = m02
	  m02 = m12
	  m12 = m22
	  m22 = m
	  p   = p01
	  p01 = p12
	  p12 = p20
	  p20 = p
	else
	  goto 21
	endif
	goto 11  
21	continue
	if (p01.eq.0d0) goto 600

C*** Regular C0 function with p01,p12,p20 =/= 0
	alpha  = sqrt( abs(kappa2(p01,p12,p20)) )
	alp(0) = sqrt( kappa2(p12,m12,m22)*(1d0+ieps*sign(1d0,p12)) )
	alp(1) = sqrt( kappa2(p20,m22,m02)*(1d0+ieps*sign(1d0,p20)) )
	alp(2) = sqrt( kappa2(p01,m02,m12)*(1d0+ieps*sign(1d0,p01)) )

	do 99 i=0,2
	  if (alp(i).eq.dcmplx(0d0,0d0)) alp(i) = ieps*abs(alpha)
99	continue

	y0(0)  = ( p12*(p12-p01-p20+2d0*m02-m12-m22)
     &	  - (p20-p01)*(m12-m22)+alpha*(p12-m12+m22) )/2d0/alpha/p12
	y0(1)  = ( p20*(p20-p12-p01+2d0*m12-m22-m02)
     &	  - (p01-p12)*(m22-m02)+alpha*(p20-m22+m02) )/2d0/alpha/p20
	y0(2)  = ( p01*(p01-p20-p12+2d0*m22-m02-m12)
     &	  - (p12-p20)*(m02-m12)+alpha*(p01-m02+m12) )/2d0/alpha/p01

	do 100 j=-1,1,2
	  x(0,j) = (p12-m12+m22+j*alp(0))/2d0/p12
	  x(1,j) = (p20-m22+m02+j*alp(1))/2d0/p20
	  x(2,j) = (p01-m02+m12+j*alp(2))/2d0/p01
	do 100 i=0,2
	  y(i,j) = y0(i)-x(i,j)
100	continue

	do 200 i=0,2
	  thp(i) = 0d0
	  thy(i) = 0d0
	  if (dimag(y(i,+1)*y(i,-1)).le.0d0) thy(i) = 1d0
200	continue
	if (p12.le.0d0) thp(0) = 1d0
	if (p20.le.0d0) thp(1) = 1d0
	if (p01.le.0d0) thp(2) = 1d0

	c0_ = 0d0
	do 300 i=0,2
	do 400 j=-1,1,2
	  c0_ = c0_ + cspens((y0(i)-1d0)/y(i,j)) - cspens(y0(i)/y(i,j))
     &	         + etass(1d0-x(i,j),1d0/y(i,j))*log((y0(i)-1d0)/y(i,j))
     &	            - etass(   -x(i,j),1d0/y(i,j))*log(y0(i)/y(i,j))
400	continue
	  c0_ = c0_ - log((y0(i)-1d0)/y0(i))*(
     &		        etass(-x(i,+1),-x(i,-1))-etass(y(i,+1),y(i,-1))
     &		      - dcmplx(0d0,2d0*pi)*thp(i)*thy(i) )
300	continue
	c0_ = c0_/alpha

	return
600	continue

	if (m02*m12.eq.0d0) goto 700
C*** Regular C0 function with p01=0 and m0,m1=/=0
	alpha   = 1d0+(p12-p20)/(m02-m12-ieps*(m02+m12))
	alp(0)  = sqrt( kappa2(p20,m02,m22)+ieps*p20*(p20-m02-m22) )
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,+1) = (p20-m02-m22+alp(0))/2d0/m02
	x(0,-1) = (p20-m02-m22-alp(0))/2d0/m02
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 610 i=-1,1,2
	do 620 j=0,1
	c0_ = c0_ + (1d0-2d0*j)*( 
     &		cspens(1d0+x(j,i)/alpha) - cspens(1d0+x(j,i)) 
     &		+ etass(1d0/alpha,-x(j,i))*log(1d0+x(j,i)/alpha) )
620	continue
610	continue
	c0_ = c0_ + log(alpha)*log(m02/m12)
	c0_ = c0_/(p20-p12) 

	return
700	continue

C*** Regular C0 function with p01=0 and m0=0
	if (m12.eq.0d0) then
	  m12 = m02
	  m02 = 0d0
	  p   = p12
	  p12 = p20
	  p20 = p
	endif
	alpha   = 1d0+(p12-p20)/(-m12-ieps*m12)
	alp(1)  = sqrt( kappa2(p12,m12,m22)+ieps*p12*(p12-m12-m22) )
	x(0,-1) = m22/(p20-m22+ieps*m12)
	x(1,+1) = (p12-m12-m22+alp(1))/2d0/m12
	x(1,-1) = (p12-m12-m22-alp(1))/2d0/m12
	c0_ = 0d0
	do 710 i=-1,1,2
	c0_ = c0_ - cspens(1d0+x(1,i)/alpha) + cspens(1d0+x(1,i)) 
     &		  - etass(1d0/alpha,-x(1,i))*log(1d0+x(1,i)/alpha) 
710	continue
	c0_ = c0_ + cspens(1d0+x(0,-1)/alpha) - cspens(1d0+x(0,-1)) 
     &		  + etass(1d0/alpha,-x(0,-1))*log(1d0+x(0,-1)/alpha) 
	c0_ = c0_ + log(alpha)*log((m22-ieps*m12-p20)/m12)
     &		  - log(alpha)**2/2d0
	c0_ = c0_/(p20-p12) 

	return
500	continue

C*** IR-divergent C0 function
        write(6,*)'warning: C0 is IR divergent'
c        stop
	sc  = p12+abs(p12)*ieps
	mm1 = sqrt(m12)
	mm4 = sqrt(m22)
	xs = -4d0*mm1*mm4/(sc-(mm1-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm1*mm4/( sc-(mm1-mm4)**2))+1d0 )**2 
	c0_ = xs/mm1/mm4/(1d0-xs**2)*(
     &	    log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)
     &		    -log(lambda2/mm1/mm4) )
     &	  - pi**2/6d0+cspens(xs**2)+log(mm1/mm4)**2/2d0
     &	  + cspens(1d0-xs*mm1/mm4) + cspens(1d0-xs*mm4/mm1) )
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPENS(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPENS,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)/                                                   
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=854513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPENS=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPENS=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPENS=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPENS=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPENS=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPENS=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
      END                                                            
***********************************************************************
        FUNCTION ETASS(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETASS,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          

	if (((IM1.eq.0d0).and.(DREAL(C1).lt.0d0)).or.
     &	    ((IM2.eq.0d0).and.(DREAL(C2).lt.0d0)).or.
     &	    ((IM12.eq.0d0).and.(DREAL(C1*C2).lt.0d0))) then
	  write(*,*) 'etass function on cut !!!'
	  write(*,*) 'C1    = ',C1
	  write(*,*) 'C2    = ',C2
	  write(*,*) 'C1*C2 = ',C1*C2
	  stop
	endif
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETASS = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETASS = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETASS = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            

************************************************************************
        FUNCTION D0_(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext) 
************************************************************************
*  SCALAR 4-POINT FUNCTION                                             *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*	General result from                                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*	IR-divergent case from                                         *
*        W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349         *
*----------------------------------------------------------------------*
************************************************************************
        IMPLICIT REAL*8 (A-Z) 
	REAL*8 M(4),P(4,4),K(4,4)                                     
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D0_,D0_ext,CSPENS,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	COMPLEX*16 SC,TC,XS,X2,X3,q2c,q3c,Y
	INTEGER I,J,ext


	if (ext.ne.0) then
	  D0_ = D0_ext(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext)
	  return
	endif

        lambda2=lambda*lambda
        PI = 4D0*DATAN(1D0)                                               

        MM1=M1    
        MM2=M2    
        MM3=M3    
        MM4=M4    
        M12=M1*M1 
        M22=M2*M2 
        M32=M3*M3 
        M42=M4*M4 
        Q1=P1 
        Q2=P2   
        Q3=P3   
	Q4=P4
        Q12=P12   
        Q23=P23
C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C--->	****** IR-divergent CASE ******
	IF ( ((Q1.EQ.M12).AND.(Q2.EQ.M32).AND.(M22.EQ.0D0)).OR.
     *	     ((Q2.EQ.M22).AND.(Q3.EQ.M42).AND.(M32.EQ.0D0)).OR.
     *	     ((Q3.EQ.M32).AND.(Q4.EQ.M12).AND.(M42.EQ.0D0)).OR.
     *	     ((Q4.EQ.M42).AND.(Q1.EQ.M22).AND.(M12.EQ.0D0)) ) goto 50

C--->	****** REGULAR CASE with at least one mass zero ******

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...   
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
C	ONLY MM1 AND MM3 ARE ZERO ==> m1 <-> m2
	IF ((MM2*MM4.NE.0D0).AND.(MM1.EQ.0D0)) then
	  mm0 = mm1
	  mm1 = mm2
	  mm2 = mm0
	  q0  = q2
	  q2  = q12
	  q12 = q0
	  q0  = q4
	  q4  = q23
	  q23 = q0
	  goto 40
	endif
C 	check whether all masses are zero
	if ((mm1.eq.0d0).and.(mm2.eq.0d0).and.(mm4.eq.0d0)) then
	  WRITE(*,*) 'D0 case with all mi=0 not implemented !'
	  stop
	endif
C 	permutate until mm1 is non-zero
	if (mm1.eq.0d0) goto 10
	goto 70

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)

c	XX(1)=XX0(1)-IEPS*DD/(XX0(1)-XX0(2))
c	XX(2)=XX0(2)+IEPS*DD/(XX0(1)-XX0(2))

c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D0_ = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D0_ = D0_ + (-1D0)**(I+J)*(
     *		CSPENS(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPENS(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) ) 
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) ) 

	   D0_ = D0_ + (3D0-2D0*I)*(
     *		 ETAS(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA(RS(1,3),1D0/RS(2,4))
     *		   *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
	   D0_ = D0_ + (3D0-2D0*I)*(
     *		+ETA(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA(-XX(I),R(1,3)/R(2,4))
     *		  +ETA(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA(R(1,3),1D0/R(2,4))
     *		   *ETA(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       -CSPENS(1D0+SS(1)*X(I,1))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPENS(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPENS(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D0_ = D0_ + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D0_ = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *		CSPENS(1D0+SS(4)*X(I,4))
     *	       +CSPENS(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	return

C	****** ONLY MM1 IS NON-ZERO ******
70	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=10D0
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	k(1,2) = (m(1)**2-p(1,2))/m(1)/m(2)
	k(1,3) = (m(1)**2-p(1,3))/m(1)/m(3)
	k(1,4) = (m(1)**2-p(1,4))/m(1)/m(4)
	k(2,3) = -p(2,3)/m(2)/m(3)
	k(2,4) = -p(2,4)/m(2)/m(4)
	k(3,4) = -p(3,4)/m(3)/m(4)
	DO 74 J=2,4
	DO 74 I=1,J-1
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
74	CONTINUE
	AA=K(2,4)*K(3,4)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	D0_ = DCMPLX(0D0,0D0)
	DO 76 I=1,2
	D0_ = D0_ + (2D0*I-3D0)*(
     *	        CSPENS(1D0+(K(1,4)-IEPS)*xx(i))
     *	       +eta(-xx(i),K(1,4)-IEPS)*log(1D0+(K(1,4)-IEPS)*xx(i))
     *	       -CSPENS(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPENS(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS)) 
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
76	CONTINUE
	D0_ = D0_/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	return

C	****** general IR-divergent D0-function ******
50	CONTINUE
        write(6,*)'warning: D0 is IR divergent'
c        stop
C	PERMUTATE UNTIL MM1 IS THE PHOTON
	GOTO 52
51	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
52	IF ((mm1.ne.0d0).or.(q1.ne.m22).or.(q4.ne.m42)) GOTO 51

	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	sc  = q23+abs(q23)*ieps
	tc  = q12+abs(q12)*ieps
	q2c = q2+abs(q2)*ieps
	q3c = q3+abs(q3)*ieps
	xs = -4d0*mm2*mm4/(sc-(mm2-mm4)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm4/( sc-(mm2-mm4)**2))+1d0 )**2 

	if (mm3.eq.0d0) goto 60

C	*** general case ***
	if (q2.ne.(mm2-mm3)**2) then
	  x2 = -4d0*mm2*mm3/(q2c-(mm2-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm2*mm3/( q2c-(mm2-mm3)**2))+1d0 )**2 
	else
	  x2 = 1d0
	endif
	if (q3.ne.(mm4-mm3)**2) then
	  x3 = -4d0*mm4*mm3/(q3c-(mm4-mm3)**2) /
     &	     ( sqrt(1d0-4d0*mm4*mm3/( q3c-(mm4-mm3)**2))+1d0 )**2 
	else
	  x3 = 1d0
	endif

	d0_ = xs/mm2/mm4/(q12-m32)/(1d0-xs**2)*(
     &	 2d0*cdlog(xs)*(cdlog(1d0-xs**2)-cdlog(lambda*mm3/(m32-tc)))
     &	+pi**2/2d0+cspens(xs**2)+cdlog(x2)**2+cdlog(x3)**2
     &	-cspens(xs*x2*x3)-(cdlog(xs)+cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs*x2*x3)
     &	-cspens(xs*x2/x3)-(cdlog(xs)+cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs*x2/x3)
     &	-cspens(xs/x2*x3)-(cdlog(xs)-cdlog(x2)+cdlog(x3))
     $       *cdlog(1d0-xs/x2*x3)
     &	-cspens(xs/x2/x3)-(cdlog(xs)-cdlog(x2)-cdlog(x3))
     $       *cdlog(1d0-xs/x2/x3) )
	return

60	continue
C	*** special case: mass mm3 opposite to photon is 0 ***
	if ((q2.eq.m22).or.(q3.eq.m42)) goto 61
	Y = mm2/mm4*(q3c-m42)/(q2c-m22)
	d0_ = xs/mm2/mm4/q12/(1d0-xs**2)*(
     &	 log(xs)*( -log(xs)/2d0+2d0*log(1d0-xs**2)-log(lambda2/mm2/mm4)
     &		   -log((q2-m22)/tc)-log((q3-m42)/tc) )
     &	+pi**2/6d0+cspens(xs**2)+log(y)**2/2d0
     &	-cspens(xs*y)-(log(xs)+log(y))*log(1d0-xs*y)
     &	-cspens(xs/y)-(log(xs)-log(y))*log(1d0-xs/y) )
	return

61	continue
C	*** special case: doubly IR-divergent D0 ***
	if ((q2.eq.m22).and.(q3.eq.m42)) then
	d0_ = -xs/mm2/mm4/q12/(1d0-xs**2)*2d0*log(xs)*log(-lambda2/tc)
	else
	  write(*,*) 'Special case of IR-divergent D0 not implemented!'
          stop
	endif

	END
************************************************************************
        function D0_ext(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,ext)
************************************************************************
*       scalar 4-point function
*-----------------------------------------------------------------------
*	special cases specified by "ext" 
*	default = general result (ext=0)
*-----------------------------------------------------------------------
************************************************************************
        implicit real*8 (a-z)
        complex*16 D0_,D0_ext,ieps,z1,z2,cspens
        integer ext

        common/param/pi,el,alpha,alpha0,alphaz,alphas,GF,cw,cw2,sw,sw2,
     &               mw,mw2,gw,mz,mz2,gz,mh,mh2,
     &               mxe,mxe2,ml(3),ml2(3),mqp(3),mqp2(3),mqm(3),mqm2(3)

	ieps = dcmplx(0d0,1d-20)
        lambda2=lambda*lambda
        
	if (ext.eq.401) then
c***    --------------------
c***	D0(me^2,kp^2,me^2,km^2,u,t,me,0,MW,0)  with me -> 0
	  if ((p1.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.mxe).and.
     &	      (m2.eq.0d0).and.(m3.eq.MW).and.(m4.eq.0d0)) then
    	    if ((p2.eq.MW2).and.(p4.eq.MW2)) then
	      u = p12
	      t = p23
	      goto 401
	    endif
	  endif
c***	D0(me^2,kp^2,me^2,km^2,u,t,MW,0,me,0)  with me -> 0
	  if ((p1.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p2.eq.MW2).and.(p4.eq.MW2)) then
	      u = p12
	      t = p23
	      goto 401
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,401) !'
	  stop
401	  continue
	  D0_ext = (-log(ml(1)/mw)**2
     &		    +log((mw2-u)/mw/ml(1)-ieps)*log(mw2/lambda2)
     &	            +log((mw2-u)/mw2-ieps)**2
     &	            +2d0*log((mw2-u)/ml(1)/mw-ieps)*log(t/mw2+ieps)
     &		    +cspens(1d0-(u-mw2)/mw2-ieps) )/t/(u-mw2)
	  return

	elseif (ext.eq.402) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,0,me,0,me)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.0d0).and.
     &	      (m2.eq.mxe).and.(m3.eq.0d0).and.(m4.eq.mxe)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2)) then
	      t = p12
	      s = p23
	      goto 402
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,402) !'
	  stop
402	  continue
	  D0_ext = ( -log(-s/ml2(1)-ieps)**2/2d0
     &		     +log(-s/ml2(1)-ieps)*( 
     &		        +log(-s/lambda2-ieps)
     &			+2d0*log(t/mw2+ieps) )
     &		    -pi**2/6d0 )/s/t
	  return

	elseif (ext.eq.403) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,0,me,0)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2)) then
	      t = p12
	      s = p23
	      goto 403
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,402) !'
	  stop
403	  continue
	  D0_ext = 2d0*log((mw2-t)/ml(1)/mw-ieps)
     &		   *log(-s/lambda2-ieps)/s/(t-mw2)
	  return

	elseif (ext.eq.404) then
c***    --------------------
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,0,me,MZ)  with me -> 0
	  if ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	      (m2.eq.0d0).and.(m3.eq.mxe).and.(m4.eq.MZ)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2).and.(p12.lt.0d0)) then
	      t = p12
	      s = p23
	      goto 404
	    endif
c***	D0(kp^2,me^2,me^2,km^2,t,s,MW,MZ,me,0)  with me -> 0
	  elseif ((p2.eq.mxe2).and.(p3.eq.mxe2).and.(m1.eq.MW).and.
     &	          (m2.eq.MZ).and.(m3.eq.mxe).and.(m4.eq.0d0)) then
    	    if ((p1.eq.MW2).and.(p4.eq.MW2).and.(p12.lt.0d0)) then
	      t = p12
	      s = p23
	      goto 404
	    endif
	  endif
	  write(*,*) 'Inconsistent call of D0_ext(...,404) !'
	  stop
404	  continue
	  z1     = dcmplx(1d0,+sqrt(4d0*mw2/mz2-1d0))/2d0
	  z2     = dcmplx(1d0,-sqrt(4d0*mw2/mz2-1d0))/2d0
	  D0_ext = (-log((mw2-t)/ml(1)/mw-ieps)**2
     &		    +log((mw2-t)/ml(1)/mw-ieps)
     &		     *( log(lambda2/ml2(1))-2d0*log(1d0-s/mz2-ieps) )
     &		    -cspens(1d0-(mw2-t)/mw2*z1)+pi**2/6d0
     &		    -cspens(1d0-(mw2-t)/mw2*z2)	    )/(t-mw2)/(mz2-s)
	  return

	endif

	D0_ext = D0_(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4,0)

	end
c******************************************************************
      subroutine dmn(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4,d1,d2,d3,ch1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       vierpunktfkt. mit bis zu 4 integrationsimpulsen in zaehler
c       realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       05.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-y)
      implicit complex*16(z)
      real*8 c1123(0:2),c2123(0:3),c1134(0:2),c2134(0:3),d2(0:6)
      real*8 c1234(0:2),c2234(0:3),c1124(0:2),c2124(0:3),d1(0:3)
      real*8 d3(0:3,0:3),ch1(0:2)
      complex*16 d0gen,d0reg
*
      q7 = q1+q2+q3+q6-q4-q5
      call chutmn(q2,q3,q7,m2,m3,m4,c1234,c2234)
      call chutmn(q1,q2,q4,m1,m2,m3,c1123,c2123)
      call chutmn(q1,q7,q6,m1,m2,m4,c1124,c2124)
      call chutmn(q4,q3,q6,m1,m3,m4,c1134,c2134)
      ch1(0) = c1234(0)
      ch1(1) = -c1234(0)-c1234(1)-c1234(2)
      ch1(2) = c1234(2)
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      p12 = (q1+q4-q2)/2d0
      p13 = (q1+q6-q7)/2d0
      p23 = (q4+q6-q3)/2d0
      det = q1*q4*q6+2d0*p12*p13*p23-q1*p23*p23-q4*p13*p13-q6*p12*p12
      mat11 = q4*q6-p23*p23
      mat12 = p13*p23-q6*p12
      mat13 = p12*p23-q4*p13
      mat21 = mat12
      mat22 = q1*q6-p13*p13
      mat23 = p12*p13-q1*p23
      mat31 = mat13
      mat32 = mat23
      mat33 = q1*q4-p12*p12
      cf1 = q1+m12-m22
      cf2 = q4+m12-m32
      cf3 = q6+m12-m42
c      if ((m22.lt.q1).and.(m32.lt.q1).and.(m42.lt.q1)) goto 50
c         d1(0) = ggttf(q7,q4,q1,m12,m22)
c      if (m12.ge.4d0*m22) then
c         d1(0) = dreal(d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4))
c      end if
c      goto 100
c--->   mass singular 3 fermion box
c 50   continue
c      if (dabs(m1-dsqrt(q1)).lt.1d-2) then
c         m1 = m1+1d-2
c         m12 = m1**2
c      end if
c      d1(0) = ggttb(q7,q4,q1,m12,m22)
c      goto 100
c       where is q4 coming from?
 100  d0new = dreal(d0reg(q1,q2,q3,q6,q7,q4,m1,m2,m3,m4))
      d1(0) = d0new
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+m12*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
c      write(6,*)q1,q2,q3,q6,q7,q4,q7,m1,m2,m3,m4
c      write(6,*)d0new,d2(0)
c      write(6,*)d1
c      write(6,*)d2
c      write(6,*)d3
      return
      end
***********************************************************************
      complex*16 function fs(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral, doppeltlang, 'regulaerer anteil'
*       f(q2,rm,rn)=b0(q2,rm,rn)-b0(0d0,rm,rn)  'subtrahiertes f'
*       q2=quadrat des die schleife durchlaufenden impulses
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       19.10.83
************************************************************************
      real*8 m,n,pi,a,s,t,b,c,d,q2,rm,rn,u,v,w,m2,n2
      data pi/3.1415926535897932384626438d0/
      m=rm
      n=rn
      m2 = m**2
      n2 = n**2
      if (dabs(m) .eq. dabs(n) ) goto 30
      if (n2 .eq. 0d0) goto 310
      if (m2 .eq. 0d0) goto 300
c---- >  allgemeiner fall
      if (q2 .ne. 0d0) goto 520
      b=0d0
      a=0d0
      goto 560
 520  u=m*m+n*n
      v=m*m-n*n
      w=m*n
      if (dabs(q2/v).le.1d-4) then
         b = u/v/v/2d0+2d0*w*w/v/v/v*dlog(n2/m2)/2d0
         b = b*q2
         a = 0d0
         goto 570
      end if
      s=dabs(m)+dabs(n)
      t=dabs(m)-dabs(n)
      c=dsqrt(dabs(s*s-q2))
      d=dsqrt(dabs(t*t-q2))
      b=1d0+(v/q2-u/v)*dlog(n2/m2)/2d0
      if (2d0*w .le. dabs(q2-u)) goto 550
      b=b-2d0*c*d/q2*datan(d/c)
      a=0d0
      goto 560
 550  a=c*d/q2
      b=b-dsign(1d0,q2-u)*a*dlog((c+d)*(c+d)/(4d0*w))
      a=pi*a
      if (q2 .ge. u) goto 560
      a=0d0
 560  continue
 570  fs=dcmplx(b,a)
      return
c---- >  gleiche massen
 30   if (q2 .ne. 0d0) goto 40
      b=0d0
      a=0d0
      goto 560
 40   u=4d0*m*m
      v=dsqrt(dabs(1d0-u/q2))
      if ((q2 .ge. 0d0) .and. (q2 .lt. u)) goto 50
      b=2d0-v*dlog((v+1d0)*(v+1d0)/u*dabs(q2))
      a=pi*v
      if (q2 .ge. u) goto 560
      a=0d0
      goto 560
 50   b=2d0-2d0*v*datan(1d0/v)
      a=0d0
      goto 560
c---- >  eine masse null
 300  m=n
 310  if (q2 .ne. m*m) goto 320
      a=0d0
      b=1d0
      goto 560
 320  b=1d0
      if(q2 .eq. 0d0) b=0d0
      a=b*(1d0-m*m/(q2+(1d0-b)))
      b=b-a*dlog(dabs(1d0-q2/m/m))
      a=pi*a
      if (q2 .gt. m*m) goto 560
      a=0d0
      goto 560
      end
***********************************************************************
      complex*16 function b0(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil (only (delta-log mue**2) subtracted)
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs
*      
      b0 = dcmplx(0d0)
      r12 = rm*rm
      r22 = rn*rn
      if (dabs(rm).eq.dabs(rn)) goto 110
      if ((r22.eq.0d0).or.(r12.eq.0d0)) goto 120
      b0 = dcmplx(-dlog(r12)/2d0-dlog(r22)/2d0+1d0-(r12+r22)/(r12-r22)
     1     *dlog(r12/r22)/2d0,0d0)
      goto 200
c---  >   beide massen gleich
 110  if (r12.eq.0d0) goto 300
      b0 = dcmplx(-dlog(r12))
      goto 200
c---  >   eine masse null
 120  b0 = 1d0-dcmplx(dlog(r12+r22))
 200  b0 = b0+fs(q2,rm,rn)
      return
c---  >   beide massen gleich null
 300  b0 = 2d0-cdlog(dcmplx(q2,1d-20))
      return
      end
***********************************************************************
      complex*16 function b1(q2,rm,rn)
***********************************************************************
*     skalares einschleifenintegral b1 minus
*     divergenter anteil (only -1/2(delta-log(mue**2)) subtracted)
*     rm,rn: massen der teilchen auf beiden armen
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,r)
      complex*16 fs,b0,xb0,xb1
*
      r12= rm*rm
      r22= rn*rn
      b1 = dcmplx(0d0)
      if (dabs(rm).eq.dabs(rn)) goto 200
      if (r12.eq.0d0) goto 120
      if (r22.eq.0d0) goto 130
      if (q2.eq.0d0) goto 140
      b1 = (r22-r12)*b0(q2,rm,rn)+r12*(1d0-dlog(r12))
     1     -r22*(1d0-dlog(r22))
      b1 = b1/2d0/q2
      goto 200
 120  b1 = r22*fs(q2,rm,rn)/2d0/q2
      goto 200
 130  b1 = -r12*fs(q2,rm,rn)/2d0/q2
 200  b1 = b1-b0(q2,rm,rn)/2d0
      goto 300
c Achtung: original code has been changed here!
 140  call bquer(q2,rm,rn,xb0,xb1)
      b1 = xb1+dlog(r22)/2d0
 300  continue
      return
      end
***********************************************************************
      complex*16 function b20(q2,rm,rn)
***********************************************************************
*       skalares einschleifenintegral b0 minus
*       divergenter anteil
*       rm,rn: massen der teilchen auf beiden armen
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       04.02.87 sa
************************************************************************
      implicit real*8(r,m,q)
      complex*16 b0,b1,b21
*      
      r12= rm*rm
      r22= rn*rn
      b20 = (q2+r12-r22)*b1(q2,rm,rn)+2d0*r12*b0(q2,rm,rn)
     1     +r12+2d0*r22-q2/3d0
      if (r22.ne.0d0) b20 = b20-r22*dlog(r22)
      b20 = b20/6d0
      return
      entry b21(q2,rm,rn)
      r12= rm*rm
      if ((dabs(rm).eq.dabs(rn)).and.(q2.eq.0d0)) goto 100
      r22= rn*rn
      b21 = 2d0*(r12-r22+q2)*b1(q2,rm,rn)+r12*b0(q2,rm,rn)
     1     +(r12-r22-q2/3d0)/2d0
      if (r22.ne.0d0) b21 = b21+r22*dlog(r22)
      b21 = -b21/3d0/q2
      return
 100  b21 = -dlog(r12)/3d0
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c0(cq1,cq2,cq3,mm1,mm2,mm3)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       q1,q2,q3=quadrate der impulse
*       m1,m2,m3: massenquadrate der teilchen auf inneren linien
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(m,q,c)
      implicit complex*16(z)
      complex*16 y(3,3),spence,cywur
      real*8 cj(3)
*
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
      m1 = mm1*mm1
      m2 = mm2*mm2
      m3 = mm3*mm3
      q1 = cq1
      q2 = cq2
      q3 = cq3
      if (q2 .eq. 0d0) then
         q2 = cq1
         q1 = cq2
         m3 = mm1*mm1
         m1 = mm3*mm3
      end if
      if (q1 .eq. 0d0) then
         j = 2
      end if
      c1 = 1d0+(q1-q3)/q2
      if (c1*c1-4d0*q1/q2 .lt. 0d0) then
         write(6,*) ' neg. determinante !'
      end if
      cwurz = dsqrt(c1*c1-4d0*q1/q2)
      calphs = (c1-cwurz)/2d0
      c0 = 0d0
      ca = q1
      cb = q2
      cc = q3-q2-q1
      cd = m2-m1-q1
      ce = q1-q3+m3-m2
      cf = m1
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = dcmplx((cnum-2d0*ca-cc*calphs)/cdenom,0d0)
      cy0 = -cc-ce
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cb*(ca+cd+cf),1d-20))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      cy0 = -ce-cd
      cywur = cdsqrt(dcmplx(cy0*cy0-4d0*cf*(ca+cb+cc),1d-20))
      y(2,3) = dcmplx(cnum/cdenom/(1d0-calphs),0d0)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      cywur = cdsqrt(dcmplx(cd*cd-4d0*ca*cf,1d-20))
      if (j .eq. 3) then
         y(3,3) = dcmplx(-cnum/calphs/cdenom,0d0)
         y(3,1) = (-cd+cywur)/2d0/ca
         y(3,2) = (-cd-cywur)/2d0/ca
      end if
      do 100 i=1,j
         c0 = c0+cj(i)*dreal(spence(y(i,3)/(y(i,3)-y(i,1)))
     1        -spence((y(i,3)-1d0)/(y(i,3)-y(i,1)))
     2        +spence(y(i,3)/(y(i,3)-y(i,2)))
     3        -spence((y(i,3)-1d0)/(y(i,3)-y(i,2))))
 100  continue
      c0 = c0/cdenom
      return
      end
***********************************************************************
      subroutine cmue(q1,q2,q3,m1,m2,m3,c1)
***********************************************************************
*     nicht-skalare dreipunktfunktion, doppeltlang,
*     q1,q2,q3=quadrate der impulse
*     m1,m2,m3: massen der teilchen auf inneren linien
*     nicht ir divergent
*     d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       20.01.87 sa
************************************************************************
      implicit real*8(a,c,g,q,m,s,d)
      implicit complex*16(z)
      complex*16 b0,spence,cscal,C0_
      real*8 c1(0:2)
*
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
c--->   two gluons on shell, 3 equal internal fermions
      if ((q1.eq.0d0).and.(q2.eq.0d0).and.(q3.gt.0d0)) goto 230
      goto 240
 230  if (4d0*m12.ge.q3) then
         awur = dsqrt(4d0*m12/q3-1d0)
         z1 = dcmplx(1d0,awur)/2d0
         z2 = dcmplx(1d0,-awur)/2d0
      else
         awur = dsqrt(1d0-4d0*m12/q3)
         al1 = 1d0+awur
         al2 = 1d0-awur
         z1 = dcmplx(al1,al2*1d-10)/2d0
         z2 = dcmplx(al2,-al2*1d-10)/2d0
      end if
      c1(0) = -spence(1d0/z1)-spence(1d0/z2)
      c1(0) = c1(0)/q3
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c     write(6,*)'Check this point in ''boxlib.f'' at line 408'
c      write(6,*)'                                       ---'
c     wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c new C0 function added on Aug.21 2004
c      c1(0) = dreal(cscal(q3,0d0,m1,m3,m2))
      c1(0)=dreal(C0_(0d0,q3,0d0,m2,m1,m3,0))
      goto 100
 240  qq2 = q2
      qq1 = q1
      qq3 = q3
      mm1 = m1
      mm2 = m2
      mm3 = m3
      if ((dabs(q1).gt.0d0).and.(q3.eq.0d0)) then
         qq3 = q1
         qq1 = q3
         mm2 = m3
         mm3 = m2
      end if
 300  c1(0) = c0(qq1,qq2,qq3,mm1,mm2,mm3)
      goto 100
 100  matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cb3 = dreal(b0(q3,m1,m3))
      cb2 = dreal(b0(q2,m2,m3))
      cb1 = dreal(b0(q1,m1,m2))
      cmp1 = (cb3-cb2-(m12-m22+q1)*c1(0))/2d0
      cmp2 = (cb1-cb3-(m22-m32+q3-q1)*c1(0))/2d0
      c1(1) = (q2*cmp1-matr*cmp2)/det
      c1(2) = (q1*cmp2-matr*cmp1)/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmn(q1,q2,q3,m1,m2,m3,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     04.02.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      complex*16 b0,b1
      real*8 det,c1(0:2),c2(0:3)
*      
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      matr = (q3-q2-q1)/2d0
      det = q1*q2-matr*matr
      cf1 = m12-m22+q1
      cf2 = m22-m32+q3-q1
      cb0 = dreal(b0(q2,m2,m3))
      call cmue(q1,q2,q3,m1,m2,m3,c1)
      c2(0) = m12/2d0*c1(0)+(cb0+cf1*c1(1)+cf2*c1(2)+1d0)/4d0
      cb1 = dreal(b1(q3,m1,m3))
      cb2 = dreal(b1(q1,m1,m2))
      cmp1 = cb1+cb0-cf1*c1(1)-2d0*c2(0)
      cmp2 = cb2-cb1-cf2*c1(1)
      c2(1) = q2*cmp1-matr*cmp2
      c2(1) = c2(1)/2d0/det
      c2(3) =-matr*cmp1+q1*cmp2
      c2(3) = c2(3)/2d0/det
      cb3 = b1(q2,m2,m3)
      cmp1 = cb1-cb3-cf1*c1(2)
      cmp2 = -cb1-cf2*c1(2)-2d0*c2(0)
      c2(2) = -matr*cmp1+q1*cmp2
      c2(2) = c2(2)/2d0/det
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chutmn(q1,q2,q3,m1,m2,m3,ch1,ch2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dreipunktfkt. mit 2 integrationsimpulsen in zaehler
c     realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     06.05.87 sa
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(c,q,m)
      real*8 c1(0:2),c2(0:3),ch1(0:2),ch2(0:3)
*
      call cmn(q1,q2,q3,m1,m2,m3,c1,c2)
      ch1(0) = c1(0)
      ch1(1) = c1(1)-c1(2)
      ch1(2) = c1(2)
      ch2(0) = c2(0)
      ch2(1) = c2(1)-2d0*c2(3)+c2(2)
      ch2(2) = c2(2)
      ch2(3) = c2(3)-c2(2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function d0gen(q1,q2,q3,q4,q5,q6,m1,m2,m3,m4)
***********************************************************************
*       skalare vierpunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q,p)
      implicit complex*16(c,z)
      complex*16 cm2(4),cmm2(4),cl(4,4),cq(4,4),ca(4),ca2(4)
      real*8 p(4,4)
*
      cm2(1) = dcmplx(m1*m1,-1d-20)
      cm2(2) = dcmplx(m2*m2,-1d-20)
      cm2(3) = dcmplx(m3*m3,-1d-20)
      cm2(4) = dcmplx(m4*m4,-1d-20)
      p(1,2) = q1
      p(1,3) = q4
      p(1,4) = q6
      p(2,3) = q2
      p(2,4) = q1+q2+q3+q6-q4-q5
      p(3,4) = q3
      do 110 i1=1,4
         do 100 i2=i1+1,4
            cl(i1,i2) = p(i1,i2)-cm2(i1)-cm2(i2)
 100     continue
 110  continue
      ca(2) = dcmplx(1d0/m1/m1)
      ca2(2) = ca(2)*ca(2)
      ca(1) = cl(1,2)-cdsqrt(cl(1,2)*cl(1,2)-4d0*cm2(1)*cm2(2))
      ca(1) = -ca(1)/2d0/cm2(1)*ca(2)
      ca2(1) = ca(1)*ca(1)
      ca(3) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,3)*ca(2)-cl(1,3)*ca(1))
      ca2(3) = ca(3)*ca(3)
      ca(4) = (-cm2(2)*ca2(2)+cm2(1)*ca2(1))
     1     /(cl(2,4)*ca(2)-cl(1,4)*ca(1))
      ca2(4) = ca(4)*ca(4)
      do 210 i1=1,4
         cmm2(i1) = -cm2(i1)*ca2(i1)
         do 200 i2=i1+1,4
            cq(i1,i2) = cl(i1,i2)*ca(i1)*ca(i2)+cm2(i1)*ca2(i1)
     1           +cm2(i2)*ca2(i2)
 200     continue
 210  continue
      caa = -cq(3,4)
      cb = -cq(2,3)
      cg = -cq(1,2)
      cc = -cq(2,4)+cq(2,3)+cq(3,4)
      ch = -cq(1,4)-cq(2,3)+cq(1,3)+cq(2,4)
      cj = -cq(1,3)+cq(1,2)+cq(2,3)
      cd = cmm2(3)-cmm2(4)+cq(3,4)
      ce = cmm2(2)-cmm2(3)+cq(2,4)-cq(3,4)
      ck = cmm2(1)-cmm2(2)+cq(1,4)-cq(2,4)
      cf = cmm2(4)
      d0gen = c0gen(caa,cb,cc,cd,ce,cf)
     1     -c0gen(caa,cb,cc,cd,(ce+ck),cf)
      d0gen = -ca(1)*ca(2)*ca(3)*ca(4)*d0gen/ck
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function c0gen(za,zb,zc,zd,ze,zf)
***********************************************************************
*       skalare dreipunktfunktion, doppeltlang,
*       za-zf: parameters a-f from 't hooft veltman
*       d o u b l e p r e c i s i o n
*-----------------------------------------------------------------------
*       07.08.90 sa
************************************************************************
      implicit real*8(m,q)
      implicit complex*16 (c,z)
      complex*16 y(3,3),c1(3),c2(3),yd1,yd2
      complex*16 spence
      real*8 cj(3)
      j = 3
      cj(1) = -1d0
      cj(2) = 1d0
      cj(3) = -1d0
*      
      cwurz = cdsqrt(zc*zc-4d0*za*zb)
      calphs = (-zc-cwurz)/2d0/zb
      c0 = dcmplx(0d0)
      ca = za
      cb = zb
      cc = zc
      cd = zd
      ce = ze
      cf = zf
      cnum = -(cd+ce*calphs)
      cdenom = cc+2d0*calphs*cb
      y(1,3) = (cnum-2d0*ca-cc*calphs)/cdenom
      cy0 = -cc-ce
      cywur = cdsqrt(cy0*cy0-4d0*cb*(ca+cd+cf))
      y(1,1) = (cy0+cywur)/2d0/cb
      y(1,2) = (cy0-cywur)/2d0/cb
      c1(1) = cb
      c2(1) = ca+cd+cf
      cy0 = -ce-cd
      cywur = cdsqrt(cy0*cy0-4d0*cf*(ca+cb+cc))
      y(2,3) = cnum/cdenom/(1d0-calphs)
      y(2,1) = (cy0+cywur)/2d0/(ca+cb+cc)
      y(2,2) = (cy0-cywur)/2d0/(ca+cb+cc)
      c1(2) = ca+cb+cc
      c2(1) = cf
      cywur = cdsqrt(cd*cd-4d0*ca*cf)
      y(3,3) = -cnum/calphs/cdenom
      y(3,1) = (-cd+cywur)/2d0/ca
      y(3,2) = (-cd-cywur)/2d0/ca
      c1(3) = ca
      c2(3) = cf
      do 100 i=1,3
         yd1 = y(i,3)-y(i,1)
         yd2 = y(i,3)-y(i,2)
         c0 = c0+cj(i)*(spence(y(i,3)/yd1)
     1        -spence((y(i,3)-1d0)/yd1)
     2        +spence(y(i,3)/yd2)
     3        -spence((y(i,3)-1d0)/yd2)
     5        )
 100  continue
      c0gen = c0/cdenom
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttb(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       mass singular scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       21.06.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      pi = 4d0*datan(1d0)
      cs = dcmplx(s,s*1d-10)
      cmh2 = dcmplx(mib2,-mib2*1d-10)
      tmt = t-mex2
      cmht = cmh2-t
      cmhmt = cmh2-mex2
      mht = mib2-t
      mhmt = mib2-mex2
      c1 = -tmt/cmht
      c2 = -s*cmh2/cmhmt/cmhmt
      zlogi = dcmplx(dlog(mif2/s),pi)
      zlogh = cdlog(cmhmt)
      zlogs = dcmplx(dlog(mib2/s),pi)
      log1 = dlog(mht)
      log2 = dlog(mib2)
      log3 = dlog((s*mib2+mhmt*mhmt)/mhmt/mhmt)
      cd0 = 4d0*spence(c1)+spence(c2)
     1     -zlogi*(2d0*zlogh-2d0*log1+zlogi/2d0)
     2     +log3*(2d0*log2-2d0*zlogh-zlogs)
      ggttb = dreal(cd0)/s/mht
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function ggttf(s,t,mex2,mib2,mif2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       scalar box integral gg-->tt
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       11.07.90 sa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-b,d-y)
      implicit complex*16(c,z)
      complex*16 spence
*
      eps1 = 1d-25
      eps = mib2*eps1
      cs = dcmplx(s,eps1*s)
      ct = dcmplx(t,-eps1*t)
      cmh2 = dcmplx(mib2,-eps)
      ctm = ct-mex2
      tm = dreal(ctm)
      tmth = t+mex2-mib2
      tmth2 = tmth*tmth
      st = s*t+tm*tm
      wur1 = dsqrt(tmth2-4d0*mex2/s*st)
      w12 = wur1*s
      y1 = (tmth+wur1)/2d0*s/st
      y2 = (tmth-wur1)/2d0*s/st
      cy1 = dcmplx(y1,eps1)
      cy2 = dcmplx(y2,-eps1)
      y11 = (y1-1d0)/y1
      y21 = (y2-1d0)/y2
      wur3 = dsqrt(tmth2-4d0*mex2*t)
      cy3 = dcmplx(tmth+wur3,eps)/2d0/t
      cy4 = dcmplx(tmth-wur3,-eps)/2d0/t
      cy31 = (cy3-1d0)/cy3
      cy41 = (cy4-1d0)/cy4
      cy13 = y1-cy3
      cy14 = y1-cy4
      cy23 = y2-cy3
      cy24 = y2-cy4
      mth = 2d0*mex2-mib2
      if (4d0*mex2.le.mib2) then
         wur5 = mib2*dsqrt(1d0-4d0*mex2/mib2)
         cy5 = dcmplx(mth+wur5,eps)/2d0/mex2
         cy6 = dcmplx(mth-wur5,-eps)/2d0/mex2
      else
         wur5 = mib2*dsqrt(4d0*mex2/mib2-1d0)
         cy5 = dcmplx(mth,wur5)/2d0/mex2
         cy6 = dcmplx(mth,-wur5)/2d0/mex2
      end if
      cy51 = (cy5-1d0)/cy5
      cy61 = (cy6-1d0)/cy6
      cy15 = y1-cy5
      cy16 = y1-cy6
      cy25 = y2-cy5
      cy26 = y2-cy6
      cy7 = s/(s+ctm)
      cy17 = y1-cy7
      cy27 = y2-cy7
      cd0 =
     1     -spence((cy1-1d0)/cy1)
     2     -spence(y1/cy17)+spence((y1-1d0)/cy17)
     3     -spence(y1/cy13)+spence((y1-1d0)/cy13)
     4     -spence(y1/cy14)+spence((y1-1d0)/cy14)
     5     +spence(y1/cy15)-spence((y1-1d0)/cy15)
     6     +spence(y1/cy16)-spence((y1-1d0)/cy16)
     7     +dlog((y1-1d0)/y1)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy1)
     8     -cdlog(cy17)-cdlog(cy13*cy14)+cdlog(cy15*cy16) )
      cd0 = cd0-(
     1     -spence((cy2-1d0)/cy2)
     2     -spence(y2/cy27)+spence((y2-1d0)/cy27)
     3     -spence(y2/cy23)+spence((y2-1d0)/cy23)
     4     -spence(y2/cy24)+spence((y2-1d0)/cy24)
     5     +spence(y2/cy25)-spence((y2-1d0)/cy25)
     6     +spence(y2/cy26)-spence((y2-1d0)/cy26)
     7     +dlog((y2-1d0)/y2)*( dlog(tm*mex2/t/(s+tm))+cdlog(cy2)
     8     -cdlog(cy27)-cdlog(cy23*cy24)+cdlog(cy25*cy26) )
     9     )
*     
      beta = dsqrt(1d0-4d0*mex2/s)
      al1 = (1d0+beta)/2d0
      al2 = (1d0-beta)/2d0
      sat = al1*s+tm
      rz1 = (-beta*tm+mib2+wur1)/2d0/sat
      rz2 = (-beta*tm+mib2-wur1)/2d0/sat
      z1 = dcmplx(rz1,eps1)
      z2 = dcmplx(rz2,-eps1)
      z3 = cmplx(al1,eps1)
      z4 = cmplx(al2,-eps1)
      z13 = rz1-z3
      z14 = rz1-z4
      z23 = rz2-z3
      z24 = rz2-z4
      z5 = -cy5+1d0
      z6 = -cy6+1d0
      z15 = rz1-al2*z5
      z16 = rz1-al2*z6
      z25 = rz2-al2*z5
      z26 = rz2-al2*z6
      cd0 = cd0+(
     1     -spence((rz1-al2)/z15)+spence(rz1/z15)
     2     -spence((rz1-al2)/z16)+spence(rz1/z16)
     3     -spence(z1/(z1-al2))
     7     -cdlog(z1/(z1-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z1)
     8     +cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     -spence((rz2-al2)/z25)+spence(rz2/z25)
     2     -spence((rz2-al2)/z26)+spence(rz2/z26)
     3     -spence(z2/(z2-al2))
     7     -cdlog(z2/(z2-al2))*( dlog(-mex2/tm/al2)-cdlog(al2-z2)
     8     +cdlog(z25*z26) )
     7     )
      z15 = rz1+al1*z5
      z16 = rz1+al1*z6
      z25 = rz2+al1*z5
      z26 = rz2+al1*z6
      z17 = rz1+al1*(1d0-cy7)
      z27 = rz2+al1*(1d0-cy7)
      zst = dcmplx(s+tm,eps)
      cd0 = cd0+(
     1     +spence((rz1+al1)/z15)-spence(rz1/z15)
     2     +spence((rz1+al1)/z16)-spence(rz1/z16)
     3     -spence((rz1+al1)/z17)+spence(rz1/z17)
     7     +cdlog(z1/(z1+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z17)+cdlog(z15*z16) )
     7     )
      cd0 = cd0-(
     1     +spence((rz2+al1)/z25)-spence(rz2/z25)
     2     +spence((rz2+al1)/z26)-spence(rz2/z26)
     3     -spence((rz2+al1)/z27)+spence(rz2/z27)
     7     +cdlog(z2/(z2+al1))*( cdlog(mex2/zst/al1)
     8     -cdlog(-z27)+cdlog(z25*z26) )
     7     )
      z15 = rz1+z3
      z25 = rz2+z3
      z16 = z15-1d0
      z26 = z25-1d0
*     
      cd0 = cd0+(
     1     +spence(z15/z16)
     2     +spence(z16/(z16+al1))-spence(z15/(z16+al1))
     3     +spence(z16/z1)-spence(z15/z1)
     7     -cdlog(z15/z16)*( cdlog(-z16)-cdlog(-z1)
     8     -cdlog(-z16-al1) )
     7     )
      cd0 = cd0-(
     1     +spence(z25/z26)
     2     +spence(z26/(z26+al1))-spence(z25/(z26+al1))
     3     +spence(z26/z2)-spence(z25/z2)
     7     -cdlog(z25/z26)*( cdlog(-z26)-cdlog(-z2)
     8     -cdlog(-z26-al1) )
     7     )
*     
 999  ggttf = dreal(cd0)/w12
      return
      end
************************************************************************
*                                                                      *
*    routines for the calculation of d0 for arbitrary real masses      *
*    no ir singularities are allowed                                   *
*    q13 < 0 required                                                  *
*                                                                      *
*    ansgar denner    18.4.96                                          *
************************************************************************
*   formulae are from:                                                 *
*     a. denner, u. nierste, r. scharf, nucl. phys. b367 (1991) 637    *
*   conventions:                                                       *
*     denominators: (q**2 - mm1**2)*[(q+k2)**2 - mm2**2]*              *
*                   [(q+k3)**2 - mm3**2]*[(q+k4)**2 - mm4**2]          *
*     invariants:   q12=k2**2, q23=(k2-k3)**2, q34=(k3-k4)**2,         *
*                   q14=k4**2, q24=(k2-k4)**2, q13=k3**2               *
************************************************************************
*                                                                      *
*   warning: d0m0,  d0m00, d0m000,d00000 go wrong for b*b = 4*a*c !!!  *
*                                                                      *
************************************************************************
* functions:                                                           *
* d0reg, d0m0,  d0m00, d0m000,d00000                                   *
* cdln,  cspenc,cspenh,cspen_new, cspcon,cspcoe, eta,ettile,etae       *
************************************************************************
      function d0reg(q12,q23,q34,q14,q24,q13,m1,m2,m3,m4)
************************************************************************
*  scalar 4-point function  for  q13 < 0 (one r real positive)         *
*  regular case                                                        *
*  imaginary part sometimes wrong                                      *
*  qij = (pi-pj)**2   and   pi is momentum of propagator with mass mi  *
*----------------------------------------------------------------------*
*  07.01.94 ansgar denner       last changed 23.10.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m3,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m12,m22,m32,m42
      real*8     ir12,ir14,ir23,ir24,ir34
      real*8     ix(2,4),is(4)
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 a,b,c
      complex*16 x(2,4),s(4)
      complex*16 d0reg,cspcoe,etae,cdln,d0m0
      integer    k,j
 
      if (m1.eq.0d0) then
        d0reg = d0m0(q23,q12,q14,q34,q24,q13,m3,m2,m4)
      else if (m2.eq.0d0) then
        d0reg = d0m0(q14,q12,q23,q34,q13,q24,m4,m1,m3)
      else if (m3.eq.0d0) then
        d0reg = d0m0(q12,q23,q34,q14,q24,q13,m1,m2,m4)
      else if (m4.eq.0d0) then
        d0reg = d0m0(q23,q34,q14,q12,q13,q24,m2,m3,m1)
      else
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12+m32-q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22+m32-q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (m32+m42-q34)/m3/m4
      if (k13.lt.2d0) then
       write(*,*) ' d0reg: case not implemented'
       write(*,*) ' k13 = ',k13
       write(*,*) ' q13 = ',q13
       write(*,*) ' m1  = ',m1
       write(*,*) ' m3  = ',m3
      end if
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2,0d0)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2,0d0)))
      r13 = 1d0/r13
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2,0d0)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2,0d0)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2,0d0)))
      r24 = 1d0/r24
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2,0d0)))
      a   =  k34/r24-k23 + (k12-k14/r24)*r13
      b   =  (1d0/r13-r13)*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k34*r24-k23 + (k12-k14*r24)/r13
      x(1,4) = (-b+sqrt(b*b-4d0*a*c))/2d0/a
      x(2,4) = (-b-sqrt(b*b-4d0*a*c))/2d0/a
      if(abs(x(1,4)).gt.abs(x(2,4))) then
        x(2,4) = c/a/x(1,4)
      else
        x(1,4) = c/a/x(2,4)
      end if

      if(k12.lt.-2d0) then
        ir12 = sign(1d1,1d0-abs(r12))
      else
        ir12 = 0d0
      end if
      if(k14.lt.-2d0) then
        ir14 = sign(1d1,1d0-abs(r14))
      else
        ir14 = 0d0
      end if
      if(k23.lt.-2d0) then
        ir23 = sign(1d1,1d0-abs(r23))
      else
        ir23 = 0d0
      end if
      if(k24.lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      if(k34.lt.-2d0) then
        ir34 = sign(1d1,1d0-abs(r34))
      else
        ir34 = 0d0
      end if
      if(dreal(x(1,4)).gt.0d0) then
         ix(1,4) = 1d0
      else
         ix(1,4) = 0d0
      end if
      if(dreal(x(2,4)).gt.0d0) then
         ix(2,4) = -1d0
      else
         ix(2,4) = 0d0
      end if

      x(1,1) = x(1,4)/r24
      x(2,1) = x(2,4)/r24
      x(1,2) = x(1,4)/r24*r13
      x(2,2) = x(2,4)/r24*r13
      x(1,3) = x(1,4)*r13
      x(2,3) = x(2,4)*r13
      s(1)  = r12
      s(2)  = r23
      s(3)  = r34
      s(4)  = r14

      is(1)  = ir12
      is(2)  = ir23
      is(3)  = ir34
      is(4)  = ir14
c --> changed 23.10.95
      if(dreal(x(1,1)).gt.0d0) then
         ix(1,1) = ix(1,4) + ir24
      else
c        ix(1,1) = 0d0
         ix(1,1) = -ix(1,4) - ir24
      end if
      if(dreal(x(2,1)).gt.0d0) then
         ix(2,1) = ix(2,4) + ir24
      else
c        ix(2,1) = 0d0
         ix(2,1) = -ix(2,4) - ir24
      end if
c --> changed 23.10.95
      ix(1,3) = ix(1,4)
      ix(2,3) = ix(2,4)
      ix(1,2) = ix(1,1) 
      ix(2,2) = ix(2,1) 
 
      d0reg = dcmplx(0d0,0d0)
      do 20 k=1,2
         do 10 j=1,4
            d0reg = d0reg + (-1d0)**(j+k) * (
     &           cspcoe(-x(k,j),s(j),-ix(k,j),is(j))
     &           + cspcoe(-x(k,j),1d0/s(j),-ix(k,j),-is(j)) )
 10      continue
         d0reg = d0reg - (-1d0)**k* 
     &        etae(-x(k,4),1d0/r24,-ix(k,4),-ir24,-ix(k,1))*
     &        cdln((1d0+k14*x(k,4)+x(k,4)**2)/(1d0+k34*x(k,3)
     $        +x(k,3)**2),
     &        -dreal(1d0+k34*x(k,3)+x(k,3)**2))
 20   continue
      d0reg = d0reg/m1/m2/m3/m4/sqrt(b*b-4d0*a*c)
      end if
      end
************************************************************************
      function d0m0(q12,q23,q34,q14,q24,q13,m1,m2,m4)
************************************************************************
*  scalar 4-point function  for m3 = 0                                 *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner       last changed 08.07.94 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     ir12,ir14,ir24
      real*8     ix1(2),ix4(2)
      complex*16 r12,r13,r14,r24
      complex*16 a,b,c,d,det
      complex*16 x1(2),x4(2)
      complex*16 d0m0,cspcoe,ettile,d0m00,cdln
      integer    i
 
      if (m1.eq.0d0) then
        d0m0 = d0m00(q12,q13,q34,q24,q14,q23,m2,m4)
      else if (m2.eq.0d0) then
        d0m0 = d0m00(q12,q23,q34,q14,q24,q13,m1,m4)
      else if (m4.eq.0d0) then
        d0m0 = d0m00(q14,q34,q23,q12,q24,q13,m1,m2)
      else
      m3  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22    -q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      a   =  k34/r24-k23
      b   =  k13*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k13*(k12-r24*k14)+r24*k34-k23
      d   = -k34*r24+k23
      det =  k12*k12*k34*k34 + k14*k14*k23*k23 + k24*k24*k13*k13
     &     - 2d0*(k12*k23*k34*k14 + k12*k24*k34*k13 + k14*k24*k23*k13)
     &     + 4d0*(k12*k23*k13 + k23*k34*k24 + k14*k34*k13)
     &     - 4d0*(k23*k23 + k34*k34 + k13*k13)
      x4(1) = (-b+sqrt(det))/2d0/a 
      x4(2) = (-b-sqrt(det))/2d0/a 
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = c/a/x4(1)
      else
        x4(1) = c/a/x4(2)
      end if
      x1(1) = x4(1)/r24
      x1(2) = x4(2)/r24

      if(k12.lt.-2d0) then
        ir12 = sign(1d1,1d0-abs(r12))
      else
        ir12 = 0d0
      end if
      if(k14.lt.-2d0) then
        ir14 = sign(1d1,1d0-abs(r14))
      else
        ir14 = 0d0
      end if
      if(k24.lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      ix4(1) = -sign(1d0,dreal(d)) 
      ix4(2) =  sign(1d0,dreal(d))
c  choice of impart avoids  i*pi*log(1-(1+eps)) terms:
c     ix4(2) = -sign(1d0,dreal(d))
c  but yields wrong imaginary part!  16.03.95
      ix1(1) =  sign(1d0,ix4(1)*dreal(r24))
      ix1(2) =  sign(1d0,ix4(2)*dreal(r24))
 
      d0m0 = dcmplx(0d0)
      do 10 i=1,2
      d0m0 = d0m0 + (2*i-3) * (
     &       cspcoe(-x4(i),r14,-ix4(i),ir14)
     &     + cspcoe(-x4(i),1d0/r14,-ix4(i),-ir14)
     &     - cspcoe(-x1(i),r12,-ix1(i),ir12) 
     &     - cspcoe(-x1(i),1d0/r12,-ix1(i),-ir12)
     &     - cspcoe(-x4(i),dcmplx(k34/k13),-ix4(i),-k13) 
     &     + cspcoe(-x1(i),dcmplx(k23/k13),-ix1(i),-k13)
     &     - ettile(-x4(i),1d0/r24,-ix4(i),-ir24) *
     &       (cdln((k12-r24*k14-(r24-1d0/r24)*x4(i))/d,
     &            dreal(-(r24-1d0/r24)*ix4(i))/d)
     &       +cdln(dcmplx(k13),-1d0)) ) 
10    continue
      d0m0 = d0m0/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d0m00(q12,q23,q34,q14,q24,q13,m1,m4)
************************************************************************
*  scalar 4-point function  for m2 = m3 = 0                            *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  10.04.92 ansgar denner       last changed 13.07.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     eps
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 k12e,k13e,r14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d0m00,cspcon,d0m000
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      if (m1.eq.0d0) then
        d0m00  = d0m000(q24,q34,q13,q12,q14,q23,m4)
      else if (m4.eq.0d0) then
        d0m00  = d0m000(q12,q23,q34,q14,q24,q13,m1)
      else
c     eps = 1d-15
      m3  = m1
      m2  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12    -q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (       -q23)/m2/m3
      k24 = (    m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
c --> changed 13.07.94
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
c --> changed 13.07.94
      a   =  k34*k24-k23
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12-k23
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      r14e  = r14*dcmplx(1d0,sign(eps,dreal(1d0/r14-r14)))
 
      d0m00 = dcmplx(0d0)
      do 10 i=1,2
      d0m00 = d0m00 + (2*i-3) * (
     &       cspcon(-x4(i),r14e) + cspcon(-x4(i),1d0/r14e)
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k23e))  )
10    continue
      d0m00 = d0m00/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d0m000(q12,q23,q34,q14,q24,q13,m1)
************************************************************************
*  scalar 4-point function  for m2 = m3 = m4 = 0                       *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  23.04.92 ansgar denner       last changed 17.03.95 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24,m1
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m2,m3,m4,m12,m22,m32,m42
      real*8     eps
      complex*16 k12e,k13e,k14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d0m000,d00000,cspcon
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      if (m1.eq.0d0) then
        d0m000 = d00000(q12,q23,q34,q14,q24,q13)
      else
c     eps = 1d-15
      m4  = m1
      m3  = m1
      m2  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12    -q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12    -q14)/m1/m4
      k23 = (       -q23)/m2/m3
      k24 = (       -q24)/m2/m4
      k34 = (       -q34)/m3/m4
      a   =  k34*k24
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12-k23
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      k14e  = k14*dcmplx(1d0,-sign(eps,k14))
 
      d0m000 = dcmplx(0d0)
      do 10 i=1,2
      d0m000 = d0m000 + (2*i-3) * (
     &       cspcon(-x4(i),k14e)
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k23e))  )
10    continue
      d0m000 = d0m000/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end if
      end
************************************************************************
      function d00000(q12,q23,q34,q14,q24,q13)
************************************************************************
*  scalar 4-point function  for m1 = m2 = m3 = m4 = 0                  *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  20.01.95 ansgar denner       last changed 16.05.93 ansgar denner    *
************************************************************************
      implicit   none
      real*8     q12,q23,q34,q14,q13,q24
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m2
      real*8     eps
      complex*16 k12e,k13e,k14e,k23e,k24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x4(2)
      complex*16 d00000,cspcon
      integer    i
c     common /eps/    eps
 
      eps = 1d-20
      m2  = abs(q24)
      k12 = (       -q12)/m2
      k13 = (       -q13)/m2
      k14 = (       -q14)/m2
      k23 = (       -q23)/m2
      k24 = (       -q24)/m2
      k34 = (       -q34)/m2
      a   =  k34*k24
      b   =  k13*k24+k12*k34-k14*k23
      c   =  k13*k12
      d   =  k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      k12e  = k12*dcmplx(1d0,-sign(eps,k12))
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k24e  = k24*dcmplx(1d0,-sign(eps,k24))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      k14e  = k14*dcmplx(1d0,-sign(eps,k14))
 
      d00000 = dcmplx(0d0)
      do 10 i=1,2
      d00000 = d00000 + (2*i-3) * (
     &     - log(-x4(i))**2/2d0
     &     - cspcon(-x4(i),k34e/k13e) - cspcon(-x4(i),k24e/k12e)
     &     + log(-x4(i))*(log(k12e)+log(k13e)-log(k14e)-log(k23e))  )
10    continue
      d00000 = d00000/m2/m2/a/(x4(1)-x4(2))
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function spen(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       spence function
c       realteil
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       02.10.89 ansgar denner
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit logical(a-z)
        complex*16 cspen_new,cx
        real*8     spen,x
 
        cx=dcmplx(x)
        spen=dreal(cspen_new(cx))
        end
************************************************************************
      function d0m0ne(q12,q23,q34,q14,q24,q13,m1,m2,m4)
************************************************************************
*  scalar 4-point function  for m3 = 0                                 *
*  regular case                                                        *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner
************************************************************************
      implicit logical (a-z)
      real*8     q12,q23,q34,q14,q13,q24,m1,m2,m4
      real*8     k12,k13,k14,k23,k24,k34
      real*8     m3,m12,m22,m32,m42
      real*8     eps
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r12e,k13e,r14e,k23e,r24e,k34e
      complex*16 a,b,c,d,cd
      complex*16 x1(2),x4(2)
      complex*16 d0m0ne,cspcon,ettile
      integer    i
 
      eps = 1d-15
      m3  = m1
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      k12 = (m12+m22-q12)/m1/m2
      k13 = (m12    -q13)/m1/m3
      k14 = (m12+m42-q14)/m1/m4
      k23 = (m22    -q23)/m2/m3
      k24 = (m22+m42-q24)/m2/m4
      k34 = (    m42-q34)/m3/m4
      r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
      a   =  k34/r24-k23
      b   =  k13*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k13*(k12-r24*k14)+r24*k34-k23
      d   = -k34*r24+k23
      cd  = c+dcmplx(0d0,eps)*d
      x4(1) = (-b+sqrt(b*b-4d0*a*cd))/2d0/a
      x4(2) = (-b-sqrt(b*b-4d0*a*cd))/2d0/a
      if(abs(x4(1)).gt.abs(x4(2))) then
        x4(2) = cd/a/x4(1)
      else
        x4(1) = cd/a/x4(2)
      end if
      x1(1) = x4(1)/r24
      x1(2) = x4(2)/r24
      k13e  = k13*dcmplx(1d0,-sign(eps,k13))
      k23e  = k23*dcmplx(1d0,-sign(eps,k23))
      k34e  = k34*dcmplx(1d0,-sign(eps,k34))
      r12e  = r12*dcmplx(1d0,sign(eps,dreal(1d0/r12-r12)))
      r14e  = r14*dcmplx(1d0,sign(eps,dreal(1d0/r14-r14)))
      r24e  = r24*dcmplx(1d0,sign(eps,dreal(1d0/r24-r24)))
 
      d0m0ne = dcmplx(0d0)
      do 10 i=1,2
      d0m0ne = d0m0ne + (2*i-3) * (
     &       cspcon(-x4(i),r14e) + cspcon(-x4(i),1d0/r14e)
     &     - cspcon(-x1(i),r12e) - cspcon(-x1(i),1d0/r12e)
     &     - cspcon(-x4(i),k34e/k13e) + cspcon(-x1(i),k23e/k13e)
     &     - ettile(-x4(i),1d0/r24,1d0/r24e) *
     &       (log((k12-r24*k14-(r24-1d0/r24)*x4(i))/d)+log(k13e))  )
10    continue
      d0m0ne = d0m0ne/m1/m2/m3/m4/a/(x4(1)-x4(2))
      end
************************************************************************
      function cdln(cz,eps)
************************************************************************
*     complex logarithm of cz + i*eps                                  *
*----------------------------------------------------------------------*
*     09.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cdln,cz,cdlog
      real*8     pi

      pi = 4d0*datan(1d0)
      if(dimag(cz).eq.0d0.and.dreal(cz).le.0.d0)then
        if (eps.eq.0d0) then
          write(80,*) 'cdln:  argument on cut '
          write(80,*) 'cdln:  eps = 0'
          write(80,*) 'cdln:  cz  = ',cz
        end if
        cdln=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
      else
        cdln=cdlog(cz)
      end if
      end
************************************************************************
      function cspenc(cz,eps)
************************************************************************
*       complex spence function  of cz + i*eps                         *
*       calculated by mapping on the area where there is a quickly     *
*       convergent series                                              *
*----------------------------------------------------------------------*
*     08.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cspenc,cz
      real*8     pi,pi26
      real*8     az,rz,az1
      complex*16 cz1,cspenh,cdln

      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      cz1     = 1d0-cz
      az1     = cdabs(cz1)
      az      = cdabs(cz)
      rz      = dreal(cz)

      if (az1.lt.1d-15) then
         cspenc = pi26
      else if (rz.lt.0.5d0) then
         if (az.lt.1d0) then
            cspenc = cspenh(cz,eps)
         else
            cspenc = -pi26 - .5d0*cdln(-cz,-eps)**2
     &           - cspenh(1d0/cz,-eps)
         end if
      else
         if (az1.lt.1d0) then
            cspenc =  pi26 - cdln(cz,eps)*cdln(cz1,-eps)
     &           - cspenh(cz1,-eps)
         else
            cspenc = 2d0*pi26 + .5d0*cdln(-cz1,-eps)**2
     &           - cdln(cz,eps)*cdln(cz1,-eps)
     &           + cspenh(1d0/cz1,eps)
         end if
      end if
      end
************************************************************************
      function cspenh(cz,eps)
************************************************************************
*       complex spence function of cz + i*eps                          *
*       in convergence region                                          *
*       calculation of bernoulli series                                *
*----------------------------------------------------------------------*
*     09.01.90 ansgar denner                                           *
************************************************************************
      implicit   none
      complex*16 cspenh,cdln,cz,x,x2
      real*8     eps
c     real*8     b(11),eps
      real*8 b(11)/
     1   0.1666666666666666666666666667d0,
     2  -0.0333333333333333333333333333d0,
     3   0.0238095238095238095238095238d0,
     4  -0.0333333333333333333333333333d0,
     5   0.0757575757575757575757575758d0,
     6  -0.2531135531135531135531135531d0,
     7   1.1666666666666666666666666667d0,
     8  -7.0921568627450980392156862745d0,
     9  54.97117794486215538847117794486d0,
     +  -529.124242424242424242424242424242d0,
     1  6192.123188405797101449275362318d0  /
c     beachte:                 b(n)=b2n
c     b(1)=1./6.
c     b(2)=-1./30.
c     b(3)=1./42.
c     b(4)=-1./30.
c     b(5)=5./66.
c     b(6)=-691./2730.
c     b(7)=7./6.
c     b(8)=-3617./510.
c     b(9)=43867./798.
c     b(10)=-174611./330.
c     b(11)=854513./138.
c     pi=3.1415926535897932384
c     pi*pi/6.=1.6449..., pi*pi/3=3.28986...
c
      integer    j
      real*8     factor
      complex*16 power,term,csp

      b(11)  =    854513d0/ 138d0
      b(10)  =  - 174611d0/ 330d0
      b(9)   =     43867d0/ 798d0
      b(8)   =  -   3617d0/ 510d0
      b(7)   =         7d0/   6d0
      b(6)   =  -    691d0/2730d0
      b(5)   =         5d0/  66d0
      b(4)   =  -      1d0/  30d0
      b(3)   =         1d0/  42d0
      b(2)   =  -      1d0/  30d0
      b(1)   =         1d0/   6d0
      x      =  -cdln(1d0-cz,-eps)
c     write(80,*)  'cspenh'
      x2     =  x*x
      power  =  x
      factor =  1d0
      cspenh =  x - x2/4d0
      do 10 j=2,22,2
         factor = factor / j / (j+1)
         power  = power * x2
         term   = b(j/2) * factor * power
         csp    = cspenh + term
         if (csp.eq.cspenh) return
         cspenh = csp
10    continue
      if (cdabs(term/csp).gt.1d-15) then
        write(80,*) 'cspenh converges badly  ',cz,x
        write(80,*) 'cspenh converges badly  ',csp-term,cspenh,term
      end if 
      end
************************************************************************
      function cspen_new(cz)
************************************************************************
*       complex spence function                                        *
*----------------------------------------------------------------------*
*     08.07.94 ansgar denner                                           *
************************************************************************
      implicit   none
      complex*16 cspen_new,cspenc,cz

      if((dimag(cz).eq.0d0).and.(dreal(cz-1d0).gt.5d-13))then
        write(80,*) 'cspen_new:  argument on cut '
        write(80,*) 'cspen_new:  cz  = ',cz
      end if
      cspen_new = cspenc(cz,0d0)
      end
************************************************************************
      function cspcon(z1,z2)
************************************************************************
*  complex spence function plus continuation terms                     *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner                                              *
************************************************************************
      implicit   none
      complex*16 cspcon,cspen_new,eta,z1,z2,cdlog
      real*8     pi,pi26
 
      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      if(dreal(z1*z2).gt.0d0) then
        cspcon = cspen_new(1d0-z1*z2) + eta(z1,z2)*cdlog(1d0-z1*z2)
      else
        cspcon = pi26-cspen_new(z1*z2)
     &           - (cdlog(z1)+cdlog(z2))*cdlog(1d0-z1*z2)
      end if
      end
************************************************************************
      function cspcoe(z1,z2,i1,i2)
************************************************************************

*  complex spence function plus continuation terms                     *
*  i@ is assumed to dominate i1                                        *
*----------------------------------------------------------------------*
*  08.07.94 ansgar denner      last changed  17.03.95                  *
************************************************************************
      implicit   none
      complex*16 cspcoe,cspenc,etae,cdln,z1,z2,z12,etaa
      real*8     i1,i2,i12
      real*8     pi,pi26

      pi = 4d0*datan(1d0)
      pi26 = pi*pi/6d0
      z12 = z1*z2
      i12 = i2*sign(1d0,dreal(z1))
      if(dreal(z12).gt.0.5d0) then
        cspcoe = cspenc(1d0-z12,0d0) 
        etaa   = etae(z1,z2,i1,i2,i12)
        if(etaa.ne.0d0) then
          cspcoe = cspcoe + etaa*cdln(1d0-z12,-i12)
        end if
      else if(abs(z12).lt.1d-4) then
         cspcoe = pi26-cspenc(z12,0d0)
     &        + (cdln(z1,i1)+cdln(z2,i2))
     &        *z12*(1d0+z12/2d0+z12*z12/3d0+z12*z12*z12/4d0)
      else
         cspcoe = pi26-cspenc(z12,0d0)
     &        - (cdln(z1,i1)+cdln(z2,i2))
     &        *cdln(1d0-z12,-0d0)
      end if

      end
************************************************************************
      function eta(c1,c2)
************************************************************************
*     complex eta-function
*----------------------------------------------------------------------*
*     8.06.90    ansgar denner       last changed   11.07.94
************************************************************************
      implicit     none
      complex*16 eta,c1,c2
      real*8     im1,im2,im12,re1,re2
      real*8     pi

      pi = 4d0*datan(1d0)
      re1    = dreal(c1)
      re2    = dreal(c2)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          eta = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          eta = dcmplx(0d0,-2d0*pi)
      else
          eta = dcmplx(0d0,0d0)
          if(.not.(im2.eq.0d0.and.re2.gt.0d0.or.
     &             im1.eq.0d0.and.re1.gt.0d0).and.
     &       (im1.eq.0d0.and.re1.lt.0d0 .or.
     &        im2.eq.0d0.and.re2.lt.0d0 .or.
     &        im12.eq.0d0.and.dreal(c1*c2).lt.0d0)) then
             write(80,*) ' eta not defined '
             write(80,*) ' eta:  c1  = ',c1
             write(80,*) ' eta:  c2  = ',c2
             write(80,*) ' eta:  c12 = ',c1*c2
          end if
      end if
      end
************************************************************************
      function etae(c1,c2,i1,i2,i12)
************************************************************************
*     complex eta-function
*----------------------------------------------------------------------*
*     8.06.90    ansgar denner       last changed   11.07.94
************************************************************************
      implicit     none
      complex*16 etae,c1,c2
      real*8     i1,i2,i12
      real*8     im1,im2,im12
      real*8     pi

      pi = 4d0*datan(1d0)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
      if(im1 .eq.0d0) im1  = i1
      if(im2 .eq.0d0) im2  = i2
      if(im12.eq.0d0) im12 = i12
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
         etae = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
         etae = dcmplx(0d0,-2d0*pi)
      else
         etae = dcmplx(0d0,0d0)
         if(im1.eq.0d0.and.dreal(c1).lt.0d0 .or.
     &        im2.eq.0d0.and.dreal(c2).lt.0d0 .or.
     &        im12.eq.0d0.and.dreal(c1*c2).lt.0d0) then
            write(80,*) ' eta not defined '
            write(80,*) ' eta:  c1  = ',c1
            write(80,*) ' eta:  c2  = ',c2
            write(80,*) ' eta:  c12 = ',c1*c2
         end if
      end if
      end
************************************************************************
      function ettile(c1,c2,i1,i2)
************************************************************************
*  complex eta-tilde   for 16-spence-d0                                *
*----------------------------------------------------------------------*
*  29.03.92 ansgar denner        last changed 14.07.94                 *
************************************************************************
      implicit     none
      complex*16 ettile,etae,c1,c2
      real*8     i1,i2
      real*8     im1,im2,re2
      real*8     pi

      pi = 4d0*datan(1d0)
      im1    = dimag(c1)
      if(im1.eq.0d0) im1 = i1
      im2    = dimag(c2)
      re2    = dreal(c2)
      if(im2.ne.0d0) then
         ettile = etae(c1,c2,i1,0d0,0d0)
      else if (re2.gt.0d0) then
         ettile = dcmplx(0d0,0d0)
      else if (im1.gt.0d0.and.i2.gt.0d0) then
         ettile = dcmplx(0d0,-2d0*pi)
      else if (im1.lt.0d0.and.i2.lt.0d0) then
         ettile = dcmplx(0d0, 2d0*pi)
      else
         ettile = dcmplx(0d0,0d0)
         if(im1.eq.0d0.and.dreal(c1).lt.0d0 .or.
     &        im2.eq.0d0.and.i2.eq.0d0.and.dreal(c1*c2).lt.0d0) then
            write(80,*) ' ettile not defined '
            write(80,*) ' ettile:  c1  = ',c1
            write(80,*) ' ettile:  c2  = ',c2
            write(80,*) ' ettile:  i1  = ',i1
            write(80,*) ' ettile:  i2  = ',i2
         end if
      end if
      end

      double precision function df(xs,ma,mb)
*
c derivative of the real part of the function f(s,ma,mb)
*
      implicit real*8(a-z)
      parameter(eps=1.d-6)
      s=xs
      if(s.lt.(ma-mb)**2) then
         rplus = dsqrt((ma+mb)**2-s)
         rmin = dsqrt((ma-mb)**2-s)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s-
     .   ((rmin**2+2.d0*rmin**2*rplus**2/s+rplus**2)/(2.d0*rmin*rplus))*
     .        dlog((rmin+rplus)**2/(4.d0*ma*mb))-1.d0 )/s
      elseif(s.lt.(ma+mb)**2) then
         rplus = dsqrt((ma+mb)**2/s-1.d0)
         rmin = dsqrt(1.d0-(ma-mb)**2/s)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s+
     .        ((rmin**2+2.d0*rmin**2*rplus**2-rplus**2)/(rmin*rplus))*
     .        datan(rmin/rplus)-1.d0 )/s
      else
         rplus = dsqrt(s-(ma+mb)**2)
         rmin = dsqrt(s-(ma-mb)**2)
         df = ((mb**2-ma**2)*dlog(mb**2/ma**2)/2d0/s-
     .   ((rmin**2-2.d0*rmin**2*rplus**2/s+rplus**2)/(2.d0*rmin*rplus))*
     .        dlog((rmin+rplus)**2/(4.d0*ma*mb))-1.d0 )/s
*        y = dsqrt(s-4.d0*ma**2)
*        z = dsqrt(s)
*        df = ((y**2-s)*dlog((z+y)/(z-y))/(2.d0*z*y)-1.d0)/s
      endif
      end
*
      double precision function g(xs,ma,mb)
*
c imaginary part of the function f(s,ma,mb)
*
      implicit double precision (a-z)
      s=xs
      pi = 4d0*datan(1d0)
      g = 0.d0
      pm = (dabs(ma)+dabs(mb))**2
      mm = (dabs(ma)-dabs(mb))**2
      if(s.gt.pm) g = pi*dsqrt((s-pm)*(s-mm))/s
      end
*
      double precision function p(xs,m)
***************************************************************
* real part of the 1-loop qed vacuumpolarisation contribution *
* from a fermion with mass m                                  *
*                                                             *
*   relation with the function f:                             *
*     p(s,m) = 1/3 - (1 + 2m**2/s)*f(s,m,m)                   *
***************************************************************
      implicit double precision (a-z)
      s=xs
      if(s.eq.0) then
        p = 0.d0
      else if(s.lt.0) then
        x = dsqrt(1.d0-4.d0*m**2/s)
        p = -8.d0/3.d0+x**2-
     .       x*(3.d0-x**2)*dlog(-4.d0*m**2/(s*(1.d0+x)**2))/2.d0
      else if(s.lt.4.d0*m**2) then
        x = dsqrt(4.d0*m**2/s-1.d0)
        p = -8.d0/3.d0-x**2+x*(3.d0+x**2)*datan(1.d0/x)
      else
        x = dsqrt(1.d0-4.d0*m**2/s)
        p = -8.d0/3.d0+x**2-
     .       x*(3.d0-x**2)*dlog(4.d0*m**2/(s*(1.d0+x)**2))/2.d0
      endif
      end
*
      complex*16 function l2(x)
      implicit double precision (a-z)
      parameter(z2=1.64493406684823d0)
* statement function (sp(u)=li2(x), u=-log(1-x) ):
      sp(u) = ((((((-1.d0/10886400.d0)*u**2+1.d0/211680.d0)*u**2-
     .     1.d0/3600.d0)*u**2+1.d0/36.d0)*u-1.d0/4.d0)*u+1.d0)*u
*
      pi = 4d0*datan(1d0)
      rel2 = -7.d0/2.d0-2.d0*x-(x+3.d0/2.d0)*dlog(x**2)
      if(x.gt.1) then
        u  = -dlog(1.d0+1.d0/x)
        iml2 = -pi*(2.d0*(1.d0+x)**2*u+3.d0+2.d0*x)
        rel2 = rel2+2.d0*(1.d0+x)**2*(-sp(u)-u*dlog(x))
      elseif(x.gt.0.d0) then
        u    = -dlog(1.d0+x)
        u1   =  dlog(x)
        iml2 = -pi*(2.d0*(1.d0+x)**2*(u+log(x))+3.d0+2.d0*x)
        rel2 =  rel2+2.d0*(1.d0+x)**2*(sp(u)+z2-u1**2/2.d0-u*u1)
      elseif(x.gt.-1.d0/2.d0) then
        u    = -log(1.d0+x)
        u1   =  log(-x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(sp(u)-2.d0*z2-u1**2/2.d0-u*u1)
      elseif(x.gt.-2.d0) then
        u    = -log(-x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(-sp(u)-z2-u**2/2.d0)
      else
        u    = -log(1.d0+1.d0/x)
        iml2 =  0.d0
        rel2 =  rel2+2.d0*(1.d0+x)**2*(-sp(u)-u*log(-x))
      endif
      l2=dcmplx(rel2,iml2)
      end
*
      complex*16 function l3(x)
      implicit double precision (a-z)
      pi = 4d0*datan(1d0)
      if(x.gt.1.d0/4.d0) then
        sq  =sqrt(4.d0*x-1.d0)
        f   =atan(1.d0/sq)
        iml3=0.d0
        rel3=(2.5d0-2.d0*x+(4.d0*x+2.d0)*sq*f-
     $       8.d0*x*(x+2.d0)*f**2)/3.d0
      elseif(x.gt.0.d0) then
        sq  =sqrt(1.d0-4.d0*x)
        f   =log((1.d0+sq)/(1.d0-sq))
        iml3=-pi*((2.d0*x+1.d0)*sq+2.d0*x*(x+2.d0)*f)/3.d0
        rel3=(2.5d0-2.d0*x+(2.d0*x+1.d0)*sq*f+
     .       2.d0*x*(x+2.d0)*(f**2-pi**2))/3.d0
      else
        sq  =sqrt(1.d0-4.d0*x)
        f   =log((sq+1.d0)/(sq-1.d0))
        iml3=0.d0
        rel3=(2.5d0-2.d0*x+(2.d0*x+1.d0)*sq*f+
     $       2.d0*x*(x+2.d0)*f**2)/3.d0
      endif
      l3=dcmplx(rel3,iml3)
      end
*
      double precision function f(xs,ma,mb)
*
c real part of the function f(s,ma,mb)
*
      implicit real*8(a-z)
      parameter(eps=1.d-6)
      s=xs
      ma2=ma**2
      mb2=mb**2
      if(abs(s).lt.eps) then
         f=0.0d0
      elseif(dabs(ma).lt.eps) then
         if(s.gt.mb2+eps) then
            f=1.d0+(1.d0-mb2/s)*dlog(1./(s/mb2-1.d0))
         elseif(s.lt.mb2-eps) then
            f=1.d0+(1.d0-mb2/s)*dlog(1./(1.d0-s/mb2))
         else
            f=1.d0
         endif
      elseif(dabs(mb).lt.eps) then
         if(s.gt.ma2+eps) then
            f=1.d0+(1.d0-ma2/s)*dlog(1./(s/ma2-1.d0))
         elseif(s.lt.ma2-eps) then
            f=1.d0+(1.d0-ma2/s)*dlog(1./(1.d0-s/ma2))
         else
            f=1.d0
         endif
      else
         if(dabs(dabs(mb)-dabs(ma)).lt.eps) then
            f=2.d0
         else
            f=1.d0+((ma2-mb2)/s-(ma2+mb2)/(ma2-mb2))*dlog(mb2/ma2)/2d0
         endif
         if(s.ge.(dabs(ma)+dabs(mb))**2) then
            rplus=dsqrt(s-(dabs(ma)+dabs(mb))**2)
            rmin =dsqrt(s-(dabs(ma)-dabs(mb))**2)
            f=f-rplus*rmin*dlog((rplus+rmin)**4/(16.d0*ma2*mb2))/2d0/s
         elseif(s.lt.(dabs(ma)-dabs(mb))**2) then
            rplus=dsqrt((dabs(ma)+dabs(mb))**2-s)
            rmin =dsqrt((dabs(ma)-dabs(mb))**2-s)
            f=f+rplus*rmin*dlog((rplus+rmin)**4/(16.d0*ma2*mb2))/2d0/s
         else
            rplus=dsqrt((dabs(ma)+dabs(mb))**2-s)
            rmin =dsqrt(s-(dabs(ma)-dabs(mb))**2)
            f=f-2.d0*rplus*rmin*datan(rmin/rplus)/s
         endif
      endif
      end
c*********************************************************
*
      subroutine cfunc(xs,mf,m1,m2,m3,c0,c1p,c1m,c20,c2p,c2m)
*
c  definition of the invariant functions in the 3-point integrals
c  with equal external masses mf.
c  s = momentum transfer; m1,m2,m3 are the internal masses
c
c                                   p2 ( = pf)
c                        m2   .
c
c               s .     .       m3
c
c                        m1   .
c                                  p1 ( = -pf)
c
c
c  c0 = scalar integral, c1p/m = c1+/-,  c2p/m= c2+/-
c
      implicit real*8(a-z)
      complex*16 c0, c1p,c1m,c20,c2p,c2m,cscal,b0s12,b031,b032,
     &           b132,b131,b1s12,c1,c2
      s=xs
      xmf=mf*mf
      c0=cscal(s,mf,m1,m2,m3)
      call bquer(s,m1,m2,b0s12,b1s12)
      call bquer(xmf,m3,m1,b031,b131)
      call bquer(xmf,m3,m2,b032,b132)
c
c  b0jk := b0(xmf,mj,mk),  b1jk:= b1(xmf,mf,mk)
c  b0s12 := b0(s,m1,m2)
c
c  the c1 functions:
c
      ms=4.d0*xmf-s
      xm1=m1**2
      xm2=m2**2
      xm3=m3**2
      c1=.25d0*dlog(xm3**2/xm1/xm2)+b0s12-.5d0*(b031+b032)
     &   +(xmf+xm3-xm1/2.d0-xm2/2.d0)*c0
      c1p=c1/ms
      c1m=(dlog(xm2/xm1)/2d0+b031-b032+(xm2-xm1)*c0)/2.d0/s
c
c  the c2 functions:
c
      c2=1.d0+b0s12+(xm1+xm2-2.d0*xm3-2.d0*xmf)*c1p+(xm1-xm2)*c1m
     &   +2.d0*xm3*c0
      c20=c2/4.d0
      c2=                     (b131+b132+2.d0*b0s12-.5d0)/4.d0
     &   +(2.d0*xm3-xm1-xm2+2.d0*xmf)/2.d0*c1p-c20
      c2p=c2/ms
      c2=                    -(b131+b132-.5d0)/4.d0
     &   -(xm1-xm2)/2.d0*c1m-c20
      c2m=c2/s
c
      return
      end
c
c********************************************************
c
      subroutine bquer(xs,m1,m2,b0,b1)
c
c  b0 and b1 are the (finite) invariant functions in the
c  2-point integrals.
c  s = q**2;  m1,m2 are the internal masses
c
      implicit real*8(a-z)
      external f,g
      complex*16 b0,b1,cf
      s=xs
      xm1=m1**2
      xm2=m2**2
      lm=dlog(xm2/xm1)/2d0
      if ((s.eq.0d0).and.(dabs(dabs(m1)-dabs(m2)).gt.1d-8)
     $     .and.(dabs(m1*m2).gt.0d0)) then
         b0 = dcmplx(1d0+(xm1+xm2)/(xm1-xm2)*lm,0d0)
         b1 = dcmplx(-1d0/4d0-(xm1+xm2)/(xm1-xm2)/4d0-
     $        xm1/(xm1-xm2)*lm-xm1*xm2/(xm1-xm2)**2*lm,0d0)
         goto 20
      endif
      cf=dcmplx(f(s,m1,m2),g(s,m1,m2))
      if (dabs(m1).eq.dabs(m2)) goto 10
      b0=1.d0-(xm2+xm1)/(xm2-xm1)*lm+cf
      b1=-.25d0+xm1/(xm2-xm1)*lm+(xm2-xm1-s)/2.d0/s*cf
      goto 20
10    continue
      b0=cf
c achtung: Aenderung am 14.8.94: factor 1/4 dazugezaehlt!
      b1=-.5d0*cf+1d0/4d0
20    return
      end
c
c********************************************************
c
      subroutine bquer1(x,m1,m2,b0,b1,p0,p1)
c
c  b0 and b1 are the (finite) inariant functions in the
c  2-point integrals, p0 and p1 their derivatives.
c  real parts only, needed for fermion renormalization.
c  x = q**2;  m1,m2 are the internal masses
c
      implicit real*8(a-z)
      external f
      cf=f(x,m1,m2)
      xm1=m1**2
      xm2=m2**2
      lm=dlog(xm2/xm1)/2d0
      if (m1.eq.m2) goto 10
      b0=1.d0-(xm2+xm1)/(xm2-xm1)*lm+cf
      b1=-.25d0+xm1/(xm2-xm1)*lm+(xm2-xm1-x)/2.d0/x*cf
      goto 20
10    b0=cf
      b1=-b0/2.d0+1d0/4d0
20    continue
c
c   calculation of the derivatives:
c
      sm=xm1+xm2
      dm=xm2-xm1
      sm12=(m1+m2)**2
      dm12=(m1-m2)**2
      s=dsqrt(dabs(sm12-x))
      d=dsqrt(dabs(dm12-x))
      klam=(dm*dm/(x*x)-sm/x)/s/d
      anf=-1.d0/x+dm/(x*x)*lm
      if (x.lt.dm12) goto 30
      if (x.gt.sm12) goto 40
      fact=2.d0*datan(d/s)
      goto 41
30    eps=1.d0
      fact=dlog(dabs((s+d)/(s-d)))
      goto 41
40    eps=-1.d0
      fact=-dlog(dabs((s+d)/(s-d)))
41    continue
      deriv=anf-klam*fact
      p0=deriv
      b1p=.5d0-lm-2.d0*b1-b0+(xm2-xm1-x)*deriv
      p1=b1p/2.d0/x
      return
      end
c
c**************************************************************
c                                                             *
c  the scalar vertex integral with equal external masses mf   *
c                                                             *
c**************************************************************
      complex*16 function cscal(xs,mf,m1,m2,m3)
c
c  s = momentum transfer; m1,m2,m3  are the internal masses
c
      implicit real*8 (a-y)
      complex*16 z1,z2,z11,z12,z21,z22,cl1,cl2,cl3,cspen,spence,int
      s=xs
      xmf=mf*mf
      if (mf.lt.1d-1) then
         mfstrich = 1d-1
         xmf=mfstrich**2
      end if
c.........xmf etc.   are fermion and boson masses squared
      xm1=m1*m1
      xm2=m2*m2
      xm3=m3*m3
c
c..t'hooft-veltman parameters
      a=1.d0
      b=xmf/s
      c=-1.d0
      d=xm1-xm2-s
      e=xm3-xm1-xmf+s
      f=xm2/s
      d=d/s
      e=e/s
c..discriminante for alpha-equation
      disc=c*c-4.d0*a*b
      if (disc .lt. 0.d0) goto 500
      al=(-c-dsqrt(disc))/2.d0/b
      nenner=c+2.d0*al*b
c..the first integral.............................................
      y0=-(d+e*al+2.d0*a+c*al)/nenner
      y01=y0-1.d0
      d1=(c+e)**2-4.d0*b*(a+d+f)
      x1=-(c+e)/2.d0/b
      if (d1.gt.0.d0) goto 10
c.......complex zeroes of logarithms
      sq1=dsqrt(-d1)
      x2=sq1/2.d0/b
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl1=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 15
10    continue
c........real zeroes
      sq1=dsqrt(d1)
      x2=sq1/2.d0/b
      y1=x1+x2
      y2=x1-x2
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      cl1=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
15    continue
c..the second integral............................................
      y0=-(d+e*al)/nenner/(1.d0-al)
      y01=y0-1.d0
      d2=(e+d)**2-4.d0*f*(a+b+c)
      x1=-(e+d)/2.d0/(a+b+c)
      if(d2.gt.0.d0) goto 20
c.......complex zeroes of logarithms
      sq2=dsqrt(-d2)
      x2=sq2/2.d0/(a+b+c)
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl2=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 25
20    continue
c........real zeroes
      x2=dsqrt(d2)/2.d0/(a+b+c)
      y1=x1+x2
      y2=x1-x2
      y11=y0/(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl2=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
25    continue
c..the third integral............................................
      y0=(d+e*al)/nenner/al
      y01=y0-1.d0
      d3=d*d-4.d0*a*f
      x1=-d/2.d0/a
      if (d3.gt.0.d0) goto 30
c........complex zeroes of logarithms
      sq3=dsqrt(-d3)
      x2=sq3/2.d0/a
      z1=dcmplx(x1,x2)
      z2=dcmplx(x1,-x2)
      z11=y0/(y0-z1)
      z12=y01/(y0-z1)
      z21=y0/(y0-z2)
      z22=y01/(y0-z2)
      cl3=spence(z11)-spence(z12)+spence(z21)-spence(z22)
      goto 35
30    continue
c........real zeroes
      x2=dsqrt(d3)/2.d0/a
      y1=x1+x2
      y2=x1-x2
 31   format(1h ,3e12.4)
      y11=y0 /(y0-y1)
      y12=y01/(y0-y1)
      y21=y0/(y0-y2)
      y22=y01/(y0-y2)
      sig1= y0/dabs(y0)
      sig2= y01/dabs(y01)
      cl3=cspen(y11,sig1)-cspen(y12,sig2)+cspen(y21,-sig1)
     &   -cspen(y22,-sig2)
35    continue
c..summation of the 3 integrals ....................................
      int=-cl1+cl2-cl3
      cscal=int/nenner/s
      goto 501
500   continue
c..error message for complex alpha................................
      write(6,21)
21    format(1h ,'  i cannot handle a complex alpha!')
      stop
501   return
      end
c
c*******************************************************************
c
      complex*16 function cspen(x,sig)
      implicit real*8(a-y)
      complex*16 z,cpi,spence,zx
      pi=4d0*datan(1d0)
      pi6=pi*pi/6.d0
      cpi=dcmplx(0.d0,pi)
      if (x.lt.1.d0) goto 10
      if(x.eq.1.d0) goto 11
      lx=dlog(x)
      x1=1.d0-x
      lx1=dlog(-x1)
      z=dcmplx(x1,0.d0)
      if (sig.gt.0.d0) goto 5
      cspen=-spence(z)+pi6-lx*(lx1+cpi)
      goto 20
5     cspen=-spence(z)+pi6-lx*(lx1-cpi)
      goto 20
10    zx=dcmplx(x,0.d0)
      cspen=spence(zx)
      goto 20
11    cspen=dcmplx(pi6,0.d0)
20    return
      end
c
c**************************************************************
c
      complex*16 function cln(x,sig)
      implicit real*8(a-z)
      complex*16 cpi
      pi=4d0*datan(1d0)
      cpi=dcmplx(0.d0,pi)
      if (x.gt.0.d0) goto 10
      x1=-x
      if (sig.gt.0.d0) goto 5
      cln=dlog(x1)-cpi
      goto 20
5     cln=dlog(x1)+cpi
      goto 20
10    y=dlog(x)
      cln=dcmplx(y,0.d0)
20    return
      end
      complex*16 function spence(xx)
c  the dilogarithm for general complex argument.
c  not allowed: real(xx) gt 1 with imag(xx)=0.
      implicit real*8(a-z)
      integer n
      complex*16 xx,x,z,d,p,cdlog
      dimension a(19)
      pi=4d0*datan(1d0)
      x=xx
      xr=dreal(x)
      xi=dimag(x)
      if(xr.ne.1.) goto 111
      if(xi.eq.0.) goto 20
111   continue
c    projection into the convergence radius
      vor=1.d0
      p=dcmplx(0.d0,0.d0)
      r=dreal(x)
      if (r .le. 0.5d0) goto 1
      p=pi*pi/6.d0- cdlog(x)*cdlog(1.d0-x)
      vor=-1.d0
      x=1.d0-x
    1 continue
      b=cdabs(x)
      if (b .lt. 1.d0) goto 2
      p=p - (pi*pi/6.d0+ cdlog(-x)*cdlog(-x)/2.d0)*vor
      vor=vor*(-1.d0)
      x=1.d0/x
    2 continue
c    calculation of the spence function
      a(1)=1.d0
      a(2)=-0.5d0
      a(3)=1.d0/6.d0
      a(5)=-1.d0/30.d0
      a(7)=1.d0/42.d0
      a(9)=-1.d0/30.d0
      a(11)=5.d0/66.d0
      a(13)=-691.d0/2730.d0
      a(15)=7.d0/6.d0
      a(17)=-3617.d0/510.d0
      a(19)=43867.d0/798.d0
      do 5 n=2,9,1
      a(2*n)=0.d0
    5 continue
      z=(-1.d0)*cdlog(1.d0-x)
      d=dcmplx(a(19),0.d0)
      do 10 n=1,18,1
      d=d*z/(20.d0-n) + a(19-n)
   10 continue
      d=d*z
      spence=d*vor + p
      goto 30
   20 continue
      spence=pi*pi/6.d0
   30 continue
      return
      end
c---------------------------------------------------------------------
*
c---------------------------------------------------------------------
      subroutine ddfunc(xpi1,xpi2,xpi3,pi12,pi13,pi23,
     $     m0,mi1,mi2,mi3,d1,d2,d3,c1,c2)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle vierpkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2, d3=(k+pi3)^2-m3^2
c xpi1=pi1^2,xpi2=pi2^2,xpi3=pi3^2,pi12=pi1.pi2,pi13=pi1.pi3,pi23=pi2.pi3
c Up to rank 3 !
c Only IR finite D-functions ! 
*---------------------------------------------------------------------
      implicit none
      integer i
      real*8 xpi1,xpi2,xpi3,pi12,pi13,pi23,q12,q13,q23,
     $     mi1,mi2,mi3,m0,xmi1,xmi2,xmi3,xm0,cf1,cf2,cf3,det
      real*8 mat11,mat12,mat13,mat21,mat22,mat23,mat31,mat32,mat33
      complex*16 d1(0:3),d2(0:6),d3(0:3,0:3),c1(0:2,4),c2(0:3,4)
      complex*16 D0_
      complex*16 c1234(0:2),c2234(0:3),c1134(0:2),c2134(0:3),
     $     c1124(0:2),c2124(0:3),c1123(0:2),c2123(0:3),ccc
      complex*16 s11,s12,s13,s20,s211,s212,s213,s221,s222,s223,s231,
     $     s232,s233,s311,s312,s313,s321,s322,s323,s331,s332,s333,
     $     s3113,s3213,s3313,s310,s320,s330
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      xmi3 = mi3*mi3
      q12=xpi1+xpi2-2d0*pi12
      q13=xpi1+xpi3-2d0*pi13
      q23=xpi2+xpi3-2d0*pi23
      cf1=xpi1-xmi1+xm0
      cf2=xpi2-xmi2+xm0
      cf3=xpi3-xmi3+xm0
c Gram determinant
      det=xpi1*xpi2*xpi3-xpi1*pi23**2-xpi2*pi13**2-xpi3*pi12**2+
     $     2d0*pi12*pi13*pi23
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ddfunc'
c         write(6,*)xpi1,xpi2,xpi3,pi23,pi13,pi12
c         stop
c      endif
c matrix elements of the inverse Gram matrix x det
      mat11=xpi2*xpi3-pi23**2
      mat12=pi13*pi23-pi12*xpi3
      mat13=pi12*pi23-pi13*xpi2
      mat21=mat12
      mat22=xpi1*xpi3-pi13**2
      mat23=pi12*pi13-pi23*xpi1
      mat31=mat13
      mat32=mat23
      mat33=xpi1*xpi2-pi12**2
*
      d1(0)=D0_(xpi1,q12,q23,xpi3,xpi2,q13,
     $     m0,mi1,mi2,mi3,0)
      call ccfunc2(q12,q13,pi23-pi12-pi13+xpi1,mi1,mi2,mi3,
     $     c1234(0),c1234(1),c1234(2),
     $     c2234(1),c2234(2),c2234(0),c2234(3))
      call ccfunc2(xpi2,xpi3,pi23,m0,mi2,mi3,
     $     c1134(0),c1134(1),c1134(2),
     $     c2134(1),c2134(2),c2134(0),c2134(3))
      call ccfunc2(xpi1,xpi3,pi13,m0,mi1,mi3,
     $     c1124(0),c1124(1),c1124(2),
     $     c2124(1),c2124(2),c2124(0),c2124(3))
      call ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c1123(0),c1123(1),c1123(2),
     $     c2123(1),c2123(2),c2123(0),c2123(3))
      do i=0,2
         c1(i,1)=c1234(i)
         c1(i,2)=c1134(i)
         c1(i,3)=c1124(i)
         c1(i,4)=c1123(i)
      end do
      do i=0,3
         c2(i,1)=c2234(i)
         c2(i,2)=c2134(i)
         c2(i,3)=c2124(i)
         c2(i,4)=c2123(i)
      end do
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      s11 = (c1134(0)-c1234(0)-cf1*d1(0))/2d0
      s12 = (c1124(0)-c1234(0)-cf2*d1(0))/2d0
      s13 = (c1123(0)-c1234(0)-cf3*d1(0))/2d0
      d1(1) = (mat11*s11+mat12*s12+mat13*s13)/det
      d1(2) = (mat21*s11+mat22*s12+mat23*s13)/det
      d1(3) = (mat31*s11+mat32*s12+mat33*s13)/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 2):-------------------
      ccc = c1234(1)+c1234(2)+c1234(0)
      s20 = c1234(0)+xm0*d1(0)
      s211 = (ccc-cf1*d1(1))/2d0
      s212 = (c1134(1)-c1234(1)-cf1*d1(2))/2d0
      s213 = (c1134(2)-c1234(2)-cf1*d1(3))/2d0
      s221 = (c1124(1)+ccc-cf2*d1(1))/2d0
      s222 = -(c1234(1)+cf2*d1(2))/2d0
      s223 = (c1124(2)-c1234(2)-cf2*d1(3))/2d0
      s231 = (c1123(1)+ccc-cf3*d1(1))/2d0
      s232 = (c1123(2)-c1234(1)-cf3*d1(2))/2d0
      s233 = -(c1234(2)+cf3*d1(3))/2d0
      d2(0) = s20-s211-s222-s233
      d2(1) = (mat11*(s211-d2(0))+mat12*s221+mat13*s231)/det
      d2(2) = (mat21*s212+mat22*(s222-d2(0))+mat23*s232)/det
      d2(3) = (mat31*s213+mat32*s223+mat33*(s233-d2(0)))/det
c---  >   d2(1,2)
      d2(4) = (mat11*s212+mat12*(s222-d2(0))+mat13*s232)/det
c---  >   d2(1,3)
      d2(5) = (mat11*s213+mat12*s223+mat13*(s233-d2(0)))/det
c---  >   d2(2,3)=d2(2,1) with m2<-->m4 for square boxes
      d2(6) = (mat21*s213+mat22*s223+mat23*(s233-d2(0)))/det
*---------------------------------------------------------------------
*--entwicklungskoeff. des tensorintegrals (Rank 3):-------------------
      s310 = (c2134(0)-c2234(0)-cf1*d2(0))/2d0
      s320 = (c2124(0)-c2234(0)-cf2*d2(0))/2d0
      s330 = (c2123(0)-c2234(0)-cf3*d2(0))/2d0
      ccc = ccc+(c2234(1)+c2234(2)-c1234(0))/2d0+c2234(3)
      s311 = -ccc-cf1*d2(1)/2d0
      s321 = (c2124(1)-cf2*d2(1))/2d0-ccc
      s331 = (c2123(1)-cf3*d2(1))/2d0-ccc
      s312 = (c2134(1)-c2234(1)-cf1*d2(2))/2d0
      s322 = -(c2234(1)+cf2*d2(2))/2d0
      s332 = (c2123(2)-c2234(1)-cf3*d2(2))/2d0
      s313 = (c2134(2)-c2234(2)-cf1*d2(3))/2d0
      s323 = (c2124(2)-c2234(2)-cf2*d2(3))/2d0
      s333 = -(c2234(2)+cf3*d2(3))/2d0
      s3113 = (c2234(3)+c2234(2)+c1234(2)-cf1*d2(5))/2d0
      s3213 = (c2124(3)+c2234(3)+c2234(2)+c1234(2)-cf2*d2(5))/2d0
      s3313 = (c2234(3)+c2234(2)+c1234(2)-cf3*d2(5))/2d0
      d3(0,1) = (mat11*s310+mat12*s320+mat13*s330)/det
      d3(0,2) = (mat21*s310+mat22*s320+mat23*s330)/det
      d3(0,3) = (mat31*s310+mat32*s320+mat33*s330)/det
      d3(1,1) = (mat11*(s311-2d0*d3(0,1))+mat12*s321+mat13*s331)/det
      d3(3,3) = (mat31*s313+mat32*s323+mat33*(s333-2d0*d3(0,3)))/det
      d3(1,2) = (mat21*(s311-2d0*d3(0,1))+mat22*s321+mat23*s331)/det
      d3(1,3) = (mat31*(s311-2d0*d3(0,1))+mat32*s321+mat33*s331)/det
      d3(2,1) = (mat11*s312+mat12*(s322-2d0*d3(0,2))+mat13*s332)/det
      d3(2,2) = (mat21*s312+mat22*(s322-2d0*d3(0,2))+mat23*s332)/det
      d3(2,3) = (mat31*s312+mat32*(s322-2d0*d3(0,2))+mat33*s332)/det
      d3(3,1) = (mat11*s313+mat12*s323+mat13*(s333-2d0*d3(0,3)))/det
      d3(3,2) = (mat21*s313+mat22*s323+mat23*(s333-2d0*d3(0,3)))/det
      d3(0,0) = (mat21*(s3113-d3(0,3))+mat22*s3213+mat23*
     1     (s3313-d3(0,1)))/det
 99   continue
      return
      end
c---------------------------------------------------------------------
      subroutine ccfunc2(xpi1,xpi2,pi12,m0,mi1,mi2,
     $     c0,c11,c12,c21,c22,c20,c23)
*---------------------------------------------------------------------
c entwicklungskoeff. fuer das vektorielle und tensorielle dreipkt.-integral
c d0=k^2-m0^2, d1=(k+pi1)^2-m1^2, d2=(k+pi2)^2-m2^2
c xpi1=pi1^2,xpi2=pi2^2,pi12=pi1.pi2
*---------------------------------------------------------------------
      implicit none
      real*8 xpi1,xpi2,pi12,mi1,mi2,m0,xmi1,xmi2,xm0,f1,f2,q12,det
      complex*16 c0,c11,c12,c21,c22,c20,c23,
     &     b01,b02,b03,b11,b12,b13,b20,C0_
      xm0 = m0*m0
      xmi1 = mi1*mi1
      xmi2 = mi2*mi2
      q12=xpi1+xpi2-2d0*pi12
      f1=xpi1-xmi1+xm0
      f2=xpi2-xmi2+xm0
      det=xpi1*xpi2-pi12**2
c      if(dabs(det).lt.1d-1)then
c         write(6,*)'warning: potential instability in ccfunc2'
c         stop
c      endif
      if(dabs(q12).ne.0d0.and.dabs(q12).lt.1d-6)q12=0d0
      if(dabs(xpi1).ne.0d0.and.dabs(xpi1).lt.1d-6)xpi1=0d0
      if(dabs(xpi2).ne.0d0.and.dabs(xpi2).lt.1d-6)xpi2=0d0
      c0=C0_(xpi1,q12,xpi2,m0,mi1,mi2,0)
c
c UV divergent B0 and B1 functions
c UV divergences have been subtracted: N_eps=2/eps-gamma_E+ln(4pi)
c
      call bfunc(q12,mi1,mi2,b01,b11,b20)
      call bfunc(xpi2,m0,mi2,b02,b12,b20)
      call bfunc(xpi1,m0,mi1,b03,b13,b20)
*---------------------------------------------------------------------
*--entwicklungskoeff. des vektorintegrals:----------------------------
      c11=1d0/2d0/det*(pi12*(b01-b03)+xpi2*(b02-b01)
     $     +(-xpi2*f1+pi12*f2)*c0)

      c12=1d0/2d0/det*(pi12*(b01-b02)+xpi1*(b03-b01)
     $     +(-xpi1*f2+pi12*f1)*c0)
*---------------------------------------------------------------------
*------------entwicklungskoeff. des tensorintegrals--------------
      c20 = (b01+1d0)/4d0+xm0/2d0*c0+(f1*c11+f2*c12)/4d0

      c21 = 1d0/2d0/det*(pi12*(-b01-b11-b13+f2*c11)+
     $     xpi2*(b01+b11-2d0*c20-f1*c11))

      c22 = 1d0/2d0/det*(pi12*(b11-b12+f1*c12)+
     $     xpi1*(-b11-2d0*c20-f2*c12))

      c23 = 1d0/2d0/det*(pi12*(b11+2d0*c20+f2*c12)+
     $     xpi2*(-b11+b12-f1*c12))
      return
      end

