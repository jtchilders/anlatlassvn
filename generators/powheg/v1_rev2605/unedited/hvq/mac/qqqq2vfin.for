Starts dribbling to qqqq2vfin.for (2010/8/18, 21:35:37).
NIL
      function qqqq2(s,t,m2,mu2)
      implicit double precision (a-z)
      integer nl
      character * 2 scheme
      common/nl/nl
      common/betfac/betfac
      common/scheme/scheme
      data pi/3.141 592 653 589 793/
      ro = 4*m2/s
      t1 = -t/s
      zg = 1
      nlf = nl
      t2 = 1-t1
      b = dsqrt(1-ro)
      lp = (b+1)/2.0E+0
      lm = (1-b)/2.0E+0
      at = s*t1
      aw = s*t2
      vlm2 = dlog(m2/mu2)
      vltm = dlog(at/m2)
      vlpm = dlog(lp/lm)
      vlsm = dlog(s/m2)
      vlsmu = dlog(s/mu2)
      vlwm = dlog(aw/m2)
      vlbl = dlog(b/lm)
      vdw = ddilog((aw-m2)/aw)-vlwm**2/2.0E+0
      vdt = ddilog((at-m2)/at)-vltm**2/2.0E+0
      vdmp = ddilog(-lm/lp)
      vdmb = vlbl**2/2.0E+0+ddilog(-lm/b)
      auinv = 1/(m2-aw)
      atinv = 1/(m2-at)
      ss = (-3.2E+1)*pi*(4*t1**2-4*t1+ro+2)*vlsmu*vlwm*zg**6/(2.7E+1*s)
      ss = 3.2E+1*pi*(4*t1**2-4*t1+ro+2)*vlsm*vlwm*zg**6/(2.7E+1*s)+ss
      ss = 6.4E+1*pi*vlwm*zg**6/(2.7E+1*s)+ss
      ss = (-1.12E+2)*pi*(4*t1**2-4*t1+ro+2)*vlsmu*vltm*zg**6/(2.7E+1*s)
     1   +ss
      ss = 1.12E+2*pi*(4*t1**2-4*t1+ro+2)*vlsm*vltm*zg**6/(2.7E+1*s)+ss
      ss = 2.24E+2*pi*vltm*zg**6/(2.7E+1*s)+ss
      ss = (-3.2E+1)*pi*(4*t1**2-4*t1+ro+2)*vlsmu**2*zg**6/(2.7E+1*s)+ss
      ss = 8.0E+0*pi*(4*t1**2-4*t1+ro+2)*vlsm*vlsmu*zg**6/(3.0E+0*s)+ss
      ss = (-4.0E+0)*pi*(ro-2)*(4*t1**2-4*t1+ro+2)*vlpm*vlsmu*zg**6/(2.7
     1   E+1*b*s)+ss
      ss = 1.6E+2*pi*(4*t1**2-4*t1+ro+2)*vlsmu*zg**6/(2.7E+1*s)+ss
      ss = (-4.0E+1)*pi*(4*t1**2-4*t1+ro+2)*vlsm**2*zg**6/(2.7E+1*s)+ss
      ss = 4.0E+0*pi*(ro-2)*(4*t1**2-4*t1+ro+2)*vlpm*vlsm*zg**6/(2.7E+1*
     1   b*s)+ss
      ss = (-1.6E+1)*pi*(40*t1**2-40*t1+10*ro+21)*vlsm*zg**6/(2.7E+1*s)+
     1   ss
      ss = 8.0E+0*pi*(ro-2)*vlpm*zg**6/(2.7E+1*b*s)+ss
      ss = (-3.2E+2)*pi*zg**6/(2.7E+1*s)+ss
      dd = (-1.6E+1)*pi*(8*t1**2-12*t1+ro+6)*vlsm*vlwm*zg**6/(2.7E+1*s)
      dd = 1.6E+1*auinv*pi*(2*t1**2-2*t1+ro)*vlwm*zg**6/2.7E+1+dd
      dd = (-6.4E+1)*pi*vlwm*zg**6/(2.7E+1*s)+dd
      dd = (-5.6E+1)*pi*(8*t1**2-4*t1+ro+2)*vlsm*vltm*zg**6/(2.7E+1*s)+d
     1   d
      dd = 5.6E+1*atinv*pi*(2*t1**2-2*t1+ro)*vltm*zg**6/2.7E+1+dd
      dd = (-2.24E+2)*pi*vltm*zg**6/(2.7E+1*s)+dd
      dd = dd-pi*((-1.6E+2)*t1**2/3.0E+0+40*t1+(-2.2E+1)*ro/3.0E+0-20)*v
     1   lsm**2*zg**6/s/9.0E+0
      dd = (-4.0E+0)*pi*(2*t1-1)*(12*ro*t1-36*t1-6*ro+8.0E+1/3.0E+0)*vls
     1   m*zg**6/(9.0E+0*b**2*s)+dd
      dd = (-4.0E+0)*pi*(12*t1**2+8*ro*t1-20*t1-3*ro+6)*vlsm*zg**6/(3.0E
     1   +0*b**4*s)+dd
      dd = 8.0E+0*pi*(8*nlf*t1**2-84*t1**2-8*nlf*t1+84*t1+2*nlf*ro-12*ro
     1   +4*nlf-40)*vlsm*zg**6/(2.7E+1*s)+dd
      dd = dd-pi*((-8.0E+0)*ro*t1**2/3.0E+0+(-2.72E+2)*t1**2/3.0E+0+(-8.
     1   0E+0)*ro*t1/3.0E+0+3.2E+2*t1/3.0E+0+(-2.0E+0)*ro**2/3.0E+0+(-1.
     2   0E+1)*ro/3.0E+0-44)*vlpm**2*zg**6/(b*s)/9.0E+0
      dd = pi*(2*t1+ro-2)*(12*ro*t1-36*t1-6*ro+8.0E+1/3.0E+0)*vlpm**2*zg
     1   **6/(b**3*s)/9.0E+0+dd
      dd = pi*(4*ro*t1**2+8*t1**2+4*ro*t1-16*t1+3*ro**2-8*ro+8)*vlpm**2*
     1   zg**6/(b**5*s)/3.0E+0+dd
      dd = 8.0E+0*b*pi*(ro+2)*(4*t1**2-4*t1+ro+1)*vlpm*zg**6/(2.7E+1*s)+
     1   dd
      dd = 4.0E+0*pi*(8*ro*t1**2-12*t1**2-8*ro*t1+12*t1-ro+2)*vlpm*zg**6
     1   /(2.7E+1*b*s)+dd
      dd = 8.0E+0*(2*nlf-33)*pi*(4*t1**2-4*t1+ro+2)*vlm2*zg**6/(2.7E+1*s
     1   )+dd
      dd = 1.6E+1*pi*(4*t1+ro-2)*vdw*zg**6/(2.7E+1*s)+dd
      dd = (-5.6E+1)*pi*(4*t1-ro-2)*vdt*zg**6/(2.7E+1*s)+dd
      dd = 8.0E+0*pi*(48*t1**2+8.0E+0*ro*t1/3.0E+0-56*t1+5.0E+0*ro/3.0E+
     1   0+7.0E+1/3.0E+0)*vdmp*zg**6/(9.0E+0*b*s)+dd
      dd = 4.0E+0*pi*(2*t1+ro-2)*(12*ro*t1-36*t1-6*ro+8.0E+1/3.0E+0)*vdm
     1   p*zg**6/(9.0E+0*b**3*s)+dd
      dd = 4.0E+0*pi*(4*ro*t1**2+8*t1**2+4*ro*t1-16*t1+3*ro**2-8*ro+8)*v
     1   dmp*zg**6/(3.0E+0*b**5*s)+dd
      dd = (-8.0E+0)*pi*(ro-2)*(4*t1**2-4*t1+ro+2)*vdmb*zg**6/(2.7E+1*b*
     1   s)+dd
      dd = (-2.0E+0)*pi**3*((-1.6E+1)*ro*t1**2/3.0E+0+(-1.12E+2)*t1**2/3
     1   .0E+0+8.0E+0*ro*t1/3.0E+0+1.36E+2*t1/3.0E+0+(-4.0E+0)*ro**2/3.0
     2   E+0+(-5.0E+0)*ro/3.0E+0-18)*zg**6/(2.7E+1*b*s)+dd
      dd = 8.0E+0*pi*(ro-1)*zg**6/(2.7E+1*b**2*s)+dd
      dd = pi**3*(2*t1+ro-2)*(12*ro*t1-36*t1-6*ro+8.0E+1/3.0E+0)*zg**6/(
     1   b**3*s)/2.7E+1+dd
      dd = (-8.0E+0)*pi*(ro-1)*(2*t1-1)**2*zg**6/(3.0E+0*b**4*s)+dd
      dd = pi**3*(4*ro*t1**2+8*t1**2+4*ro*t1-16*t1+3*ro**2-8*ro+8)*zg**6
     1   /(b**5*s)/9.0E+0+dd
      dd = -pi*(192*ro*t1**2-80*pi**2*t1**2+320*nlf*t1**2-2368*t1**2-192
     1   *ro*t1+40*pi**2*t1-320*nlf*t1+2368*t1+48*ro**2-2*pi**2*ro+80*nl
     2   f*ro-496*ro-20*pi**2+160*nlf-2384)*zg**6/s/8.1E+1+(-3.2E+1)*pi*
     3   *3*(4*t1**2-4*t1+ro+2)*zg**6/(8.1E+1*s)+dd
      qqqq2 = (ss+dd)/pi/4.0E+0
      return
      end
