Starts dribbling to ggqqvfin.for (2010/8/18, 15:46:35).
NIL
      function ggqq(s,t,m2,mu2)
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
      vca = 3
      vcf = 4
      vcf = vcf/3.0E+0
      vtf = 1
      vtf = vtf/2.0E+0
      vda = 8
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
      lt1 = log(t1)
      lt2 = log(t2)
      ss = (8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*vca
     1   *(vcf+t1**2*vca-t1*vca)*vlsmu**2*vtf*zg**6/(s*(t1-1)**2*t1**2*v
     2   da)/2.0E+0
      ss = (8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*vca
     1   *(4*vcf+2*t1**2*vca-2*t1*vca-vca)*vlsm*vlsmu*vtf*zg**6/(s*(t1-1
     2   )**2*t1**2*vda)/4.0E+0+ss
      tmp0 = 2*vcf+2*t1**2*vca-2*t1*vca-vca
      tmp0 = -(ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+
     1   ro**2)*tmp0*(2*vcf-vca)*vlpm*vlsmu*vtf*zg**6/(b*s*(t1-1)**2*t1*
     2   *2*vda)/8.0E+0
      ss = tmp0+ss
      ss = lt2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)
     1   *vca*(2*vcf+t1**2*vca-vca)*vlsmu*vtf*zg**6/(s*(t1-1)**2*t1**2*v
     2   da)/2.0E+0+ss
      ss = lt1*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)
     1   *vca*(2*vcf+t1**2*vca-2*t1*vca)*vlsmu*vtf*zg**6/(s*(t1-1)**2*t1
     2   **2*vda)/2.0E+0+ss
      tmp0 = 4*nlf*vtf-6*vcf-11*vca
      tmp0 = (8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*t
     1   mp0*(vcf+t1**2*vca-t1*vca)*vlsmu*vtf*zg**6/(s*(t1-1)**2*t1**2*v
     2   da)/6.0E+0
      ss = tmp0+ss
      ss = ss-(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*
     1   vca*(6*vcf+4*t1**2*vca-4*t1*vca-vca)*vlsm**2*vtf*zg**6/(s*(t1-1
     2   )**2*t1**2*vda)/4.0E+0
      tmp0 = 2*vcf+2*t1**2*vca-2*t1*vca-vca
      tmp0 = (ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+r
     1   o**2)*tmp0*(2*vcf-vca)*vlpm*vlsm*vtf*zg**6/(b*s*(t1-1)**2*t1**2
     2   *vda)/8.0E+0
      ss = tmp0+ss
      ss = ss-lt2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro*
     1   *2)*vca*(2*vcf+t1**2*vca-vca)*vlsm*vtf*zg**6/(s*(t1-1)**2*t1**2
     2   *vda)/2.0E+0
      ss = ss-lt1*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro*
     1   *2)*vca*(2*vcf+t1**2*vca-2*t1*vca)*vlsm*vtf*zg**6/(s*(t1-1)**2*
     2   t1**2*vda)/2.0E+0
      tmp0 = 32*nlf*t1**4*vcf*vtf-64*nlf*t1**3*vcf*vtf+16*nlf*ro*t1**2*v
     1   cf*vtf+48*nlf*t1**2*vcf*vtf-16*nlf*ro*t1*vcf*vtf-16*nlf*t1*vcf*
     2   vtf+4*nlf*ro**2*vcf*vtf+32*nlf*t1**6*vca*vtf-96*nlf*t1**5*vca*v
     3   tf+16*nlf*ro*t1**4*vca*vtf
      tmp0 = 112*nlf*t1**4*vca*vtf-32*nlf*ro*t1**3*vca*vtf-64*nlf*t1**3*
     1   vca*vtf+4*nlf*ro**2*t1**2*vca*vtf+16*nlf*ro*t1**2*vca*vtf+16*nl
     2   f*t1**2*vca*vtf-4*nlf*ro**2*t1*vca*vtf-48*t1**4*vcf**2+96*t1**3
     3   *vcf**2-24*ro*t1**2*vcf**2-72*t1**2*vcf**2+24*ro*t1*vcf**2+24*t
     4   1*vcf**2-6*ro**2*vcf**2-48*t1**6*vca*vcf+144*t1**5*vca*vcf-24*r
     5   o*t1**4*vca*vcf-352*t1**4*vca*vcf+48*ro*t1**3*vca*vcf+464*t1**3
     6   *vca*vcf-6*ro**2*t1**2*vca*vcf-164*ro*t1**2*vca*vcf-252*t1**2*v
     7   ca*vcf+6*ro**2*t1*vca*vcf+140*ro*t1*vca*vcf+44*t1*vca*vcf-35*ro
     8   **2*vca*vcf-160*t1**6*vca**2+480*t1**5*vca**2-116*ro*t1**4*vca*
     9   *2-512*t1**4*vca**2+232*ro*t1**3*vca**2+224*t1**3*vca**2-29*ro*
     :   *2*t1**2*vca**2-104*ro*t1**2*vca**2-32*t1**2*vca**2+29*ro**2*t1
     ;   *vca**2-12*ro*t1*vca**2+3*ro**2*vca**2+tmp0
      tmp0 = -tmp0*vlsm*vtf*zg**6/(s*(t1-1)**2*t1**2*vda)/6.0E+0
      ss = tmp0+ss
      tmp0 = 2*vcf+2*t1**2*vca-2*t1*vca-vca
      tmp0 = -(ro-2)*(2*t1**2-2*t1+ro)**2*tmp0*(2*vcf-vca)*vlpm*vtf*zg**
     1   6/(b*s*(t1-1)**2*t1**2*vda)/4.0E+0
      ss = tmp0+ss
      ss = lt2*(2*t1**2-2*t1+ro)**2*vca*(2*vcf+t1**2*vca-vca)*vtf*zg**6/
     1   (s*(t1-1)**2*t1**2*vda)+ss
      ss = lt1*(2*t1**2-2*t1+ro)**2*vca*(2*vcf+t1**2*vca-2*t1*vca)*vtf*z
     1   g**6/(s*(t1-1)**2*t1**2*vda)+ss
      tmp0 = 8*nlf*t1**4*vtf-16*nlf*t1**3*vtf+8*nlf*ro*t1**2*vtf+8*nlf*t
     1   1**2*vtf-8*nlf*ro*t1*vtf+2*nlf*ro**2*vtf-12*t1**4*vcf+24*t1**3*
     2   vcf-12*ro*t1**2*vcf-12*t1**2*vcf+12*ro*t1*vcf-3*ro**2*vcf-34*t1
     3   **4*vca+68*t1**3*vca-40*ro*t1**2*vca-34*t1**2*vca+40*ro*t1*vca-
     4   10*ro**2*vca
      tmp0 = 2.0E+0*tmp0*(vcf+t1**2*vca-t1*vca)*vtf*zg**6/(3.0E+0*s*(t1-
     1   1)**2*t1**2*vda)
      ss = tmp0+ss
      dd = 2*vcf+t1*vca-vca
      dd = -auinv*dd*pi*(4*t1+ro-4)*(8*t1**4-16*t1**3+12*ro*t1**2+8*t1**
     1   2-12*ro*t1+3*ro**2)*(2*vcf-vca)*vlwm**2*vtf*zg**6/((t1-1)**3*t1
     2   *vda)/8.0E+0
      tmp0 = 8*ro*t1**3*vcf+16*t1**3*vcf+4*ro**2*t1**2*vcf-8*ro*t1**2*vc
     1   f-48*t1**2*vcf-8*ro**2*t1*vcf+4*ro*t1*vcf+48*t1*vcf+6*ro**2*vcf
     2   -4*ro*vcf-16*vcf+4*ro*t1**4*vca+4*t1**4*vca+2*ro**2*t1**3*vca-1
     3   2*ro*t1**3*vca-24*t1**3*vca-8*ro**2*t1**2*vca+18*ro*t1**2*vca+4
     4   8*t1**2*vca+11*ro**2*t1*vca-16*ro*t1*vca-40*t1*vca-5*ro**2*vca+
     5   6*ro*vca+12*vca
      tmp0 = pi*tmp0*(2*vcf-vca)*vlwm**2*vtf*zg**6/(s*(t1-1)**3*t1*vda)
      dd = tmp0+dd
      dd = 4*pi*(4*t1**4-8*t1**3+ro**2*t1**2+2*ro*t1**2+6*t1**2-ro**2*t1
     1   -2*ro*t1-2*t1+ro**2)*vca*(2*vcf-vca)*vltm*vlwm*vtf*zg**6/(s*(t1
     2   -1)**2*t1**2*vda)+dd
      dd = 2*pi*(8*t1**4-24*t1**3+2*ro*t1**2+30*t1**2-2*ro*t1-18*t1+ro**
     1   2+4)*vca**2*vlsm*vlwm*vtf*zg**6/(s*(t1-1)**2*vda)+dd
      tmp0 = 2*vcf+t1*vca-2*vca
      tmp0 = pi*(4*t1**3+2*ro**2*t1**2+2*ro*t1**2-12*t1**2-6*ro**2*t1-2*
     1   ro*t1+16*t1-ro**3+6*ro**2-8)*tmp0*(2*vcf-vca)*vlpm*vlwm*vtf*zg*
     2   *6/(b*s*(t1-1)**2*t1*vda)
      dd = tmp0+dd
      tmp0 = 96*t1**5*vcf+16*ro*t1**4*vcf-256*t1**4*vcf+96*ro*t1**3*vcf+
     1   256*t1**3*vcf+24*ro**2*t1**2*vcf-192*ro*t1**2*vcf-192*t1**2*vcf
     2   +12*ro**2*t1*vcf+32*ro*t1*vcf+160*t1*vcf+6*ro**3*vcf-36*ro**2*v
     3   cf+48*ro*vcf-64*vcf-32*t1**5*vca-8*ro*t1**4*vca-56*ro*t1**3*vca
     4   +192*t1**3*vca-12*ro**2*t1**2*vca+104*ro*t1**2*vca-256*t1**2*vc
     5   a-8*ro**2*t1*vca-8*ro*t1*vca+96*t1*vca-3*ro**3*vca+20*ro**2*vca
     6   -32*ro*vca
      tmp0 = auinv**2*pi*s*tmp0*(2*vcf+t1*vca-vca)*vlwm*vtf*zg**6/((t1-1
     1   )**2*t1*vda)/8.0E+0
      dd = tmp0+dd
      tmp0 = 16*t1**4*vcf+8*ro*t1**3*vcf-20*t1**3*vcf+2*ro**2*t1**2*vcf+
     1   2*ro*t1**2*vcf+2*ro**2*t1*vcf-10*ro*t1*vcf-4*t1*vcf+ro**3*vcf-2
     2   *ro**2*vcf+8*vcf+4*t1**4*vca-28*t1**3*vca-2*ro*t1**2*vca+44*t1*
     3   *2*vca-4*ro*t1*vca-20*t1*vca-ro**2*vca+6*ro*vca
      tmp0 = -auinv*pi*tmp0*(2*vcf+t1*vca-vca)*vlwm*vtf*zg**6/((t1-1)**2
     1   *t1*vda)
      dd = tmp0+dd
      dd = dd-4*pi*(2*t1**2-2*t1+ro)**2*vca*(2*vcf+t1**2*vca-vca)*vlwm*v
     1   tf*zg**6/(s*(t1-1)**2*t1**2*vda)
      tmp0 = 2*vcf-t1*vca
      tmp0 = atinv*pi*(4*t1-ro)*(8*t1**4-16*t1**3+12*ro*t1**2+8*t1**2-12
     1   *ro*t1+3*ro**2)*tmp0*(2*vcf-vca)*vltm**2*vtf*zg**6/((t1-1)*t1**
     2   3*vda)/8.0E+0
      dd = tmp0+dd
      dd = dd-pi*(2*vcf-vca)*(8*ro*t1**3*vcf+16*t1**3*vcf-4*ro**2*t1**2*
     1   vcf-16*ro*t1**2*vcf+12*ro*t1*vcf-2*ro**2*vcf-4*ro*t1**4*vca-4*t
     2   1**4*vca+2*ro**2*t1**3*vca+4*ro*t1**3*vca-8*t1**3*vca+2*ro**2*t
     3   1**2*vca-6*ro*t1**2*vca+ro**2*t1*vca)*vltm**2*vtf*zg**6/(s*(t1-
     4   1)*t1**3*vda)
      dd = 2*pi*(8*t1**4-8*t1**3+2*ro*t1**2+6*t1**2-2*ro*t1-2*t1+ro**2)*
     1   vca**2*vlsm*vltm*vtf*zg**6/(s*t1**2*vda)+dd
      tmp0 = 2*vcf-t1*vca-vca
      tmp0 = pi*(4*t1**3-2*ro**2*t1**2-2*ro*t1**2-2*ro**2*t1+2*ro*t1+4*t
     1   1+ro**3-2*ro**2)*tmp0*(2*vcf-vca)*vlpm*vltm*vtf*zg**6/(b*s*(t1-
     2   1)*t1**2*vda)
      dd = tmp0+dd
      dd = atinv**2*pi*s*(2*vcf-t1*vca)*(96*t1**5*vcf-16*ro*t1**4*vcf-22
     1   4*t1**4*vcf+160*ro*t1**3*vcf+192*t1**3*vcf-24*ro**2*t1**2*vcf-1
     2   92*ro*t1**2*vcf+60*ro**2*t1*vcf-6*ro**3*vcf-32*t1**5*vca+8*ro*t
     3   1**4*vca+160*t1**4*vca-88*ro*t1**3*vca-128*t1**3*vca+12*ro**2*t
     4   1**2*vca+112*ro*t1**2*vca-32*ro**2*t1*vca+3*ro**3*vca)*vltm*vtf
     5   *zg**6/((t1-1)*t1**2*vda)/8.0E+0+dd
      dd = atinv*pi*(2*vcf-t1*vca)*(16*t1**4*vcf-8*ro*t1**3*vcf-44*t1**3
     1   *vcf+2*ro**2*t1**2*vcf+26*ro*t1**2*vcf+36*t1**2*vcf-6*ro**2*t1*
     2   vcf-18*ro*t1*vcf+ro**3*vcf+2*ro**2*vcf+4*t1**4*vca+12*t1**3*vca
     3   -2*ro*t1**2*vca-16*t1**2*vca+8*ro*t1*vca-ro**2*vca)*vltm*vtf*zg
     4   **6/((t1-1)*t1**2*vda)+dd
      dd = dd-4*pi*(2*t1**2-2*t1+ro)**2*vca*(2*vcf+t1**2*vca-2*t1*vca)*v
     1   ltm*vtf*zg**6/(s*(t1-1)**2*t1**2*vda)
      dd = pi*(ro+1)*(2*t1**2-2*t1+1)*vca**2*vlsm**2*vtf*zg**6/(s*(t1-1)
     1   *t1*vda)+dd
      tmp0 = 16*t1**4*vcf-32*t1**3*vcf-12*ro**2*t1**2*vcf+4*ro*t1**2*vcf
     1   +40*t1**2*vcf+12*ro**2*t1*vcf-4*ro*t1*vcf-24*t1*vcf-2*ro**3*vcf
     2   +4*ro**2*vcf+4*ro**2*t1**4*vca+4*ro*t1**4*vca-20*t1**4*vca-8*ro
     3   **2*t1**3*vca-8*ro*t1**3*vca+40*t1**3*vca-2*ro**3*t1**2*vca+18*
     4   ro**2*t1**2*vca+2*ro*t1**2*vca-40*t1**2*vca+2*ro**3*t1*vca-14*r
     5   o**2*t1*vca+2*ro*t1*vca+20*t1*vca+ro**3*vca-2*ro**2*vca
      tmp0 = -pi*tmp0*(2*vcf-vca)*vlpm*vlsm*vtf*zg**6/(b*s*(t1-1)**2*t1*
     1   *2*vda)/2.0E+0
      dd = tmp0+dd
      dd = pi*(2*t1-1)**2*(4*ro*t1**2-12*t1**2-4*ro*t1+12*t1-ro**2+ro-4)
     1   *vca**2*vlsm*vtf*zg**6/(b**2*s*(t1-1)*t1*vda)+dd
      dd = 3.0E+0*pi*(2*t1-1)**2*(4*t1**2-4*t1-ro+2)*vca**2*vlsm*vtf*zg*
     1   *6/(2.0E+0*b**4*s*(t1-1)*t1*vda)+dd
      dd = dd-pi*(16*t1**6-48*t1**5+24*ro*t1**4+48*t1**4-48*ro*t1**3-16*
     1   t1**3+4*ro**2*t1**2+34*ro*t1**2-t1**2-4*ro**2*t1-10*ro*t1+t1+2*
     2   ro**2)*vca**2*vlsm*vtf*zg**6/(s*(t1-1)**2*t1**2*vda)
      dd = pi*(32*ro*t1**2*vcf**2-32*t1**2*vcf**2-32*ro*t1*vcf**2+32*t1*
     1   vcf**2-8*ro**2*vcf**2-8*ro*vcf**2+16*vcf**2+32*ro*t1**4*vca*vcf
     2   -64*t1**4*vca*vcf-64*ro*t1**3*vca*vcf+128*t1**3*vca*vcf+24*ro**
     3   2*t1**2*vca*vcf-8*ro*t1**2*vca*vcf-88*t1**2*vca*vcf-24*ro**2*t1
     4   *vca*vcf+40*ro*t1*vca*vcf+24*t1*vca*vcf+16*ro**2*vca*vcf+8*ro*v
     5   ca*vcf-32*vca*vcf-16*ro*t1**4*vca**2-32*t1**4*vca**2+32*ro*t1**
     6   3*vca**2+64*t1**3*vca**2-12*ro**2*t1**2*vca**2-5*ro*t1**2*vca**
     7   2-66*t1**2*vca**2+12*ro**2*t1*vca**2-11*ro*t1*vca**2+34*t1*vca*
     8   *2-6*ro**2*vca**2-2*ro*vca**2+4*vca**2)*vlpm**2*vtf*zg**6/(b*s*
     9   (t1-1)*t1*vda)/4.0E+0+dd
      dd = dd-pi*(16*ro*t1**4-48*t1**4-32*ro*t1**3+96*t1**3-7*ro**2*t1**
     1   2+49*ro*t1**2-98*t1**2+7*ro**2*t1-33*ro*t1+50*t1+ro**3-3*ro**2+
     2   6*ro-8)*vca**2*vlpm**2*vtf*zg**6/(b**3*s*(t1-1)*t1*vda)/4.0E+0
      dd = dd-pi*(2*t1-1)**2*(4*ro*t1**2+8*t1**2-4*ro*t1-8*t1+3*ro**2-8*
     1   ro+8)*vca**2*vlpm**2*vtf*zg**6/(b**5*s*(t1-1)*t1*vda)/8.0E+0
      dd = dd-pi*vlpm**2*vtf*(8*ro**2*t1**2*vca*vtf-8*ro**2*t1*vca*vtf+2
     1   *ro**2*vca*vtf+32*t1**2*vcf**2-32*t1*vcf**2-8*ro**2*vcf**2+48*v
     2   cf**2-16*ro*t1**2*vca*vcf-56*t1**2*vca*vcf+16*ro*t1*vca*vcf+56*
     3   t1*vca*vcf+10*ro**2*vca*vcf-4*ro*vca*vcf-64*vca*vcf+8*ro*t1**2*
     4   vca**2+20*t1**2*vca**2-8*ro*t1*vca**2-20*t1*vca**2-3*ro**2*vca*
     5   *2+2*ro*vca**2+20*vca**2)*zg**6/(s*(t1-1)*t1*vda)/4.0E+0
      dd = b*pi*vlpm*vtf*(8*ro**2*t1**2*vca*vtf-8*ro**2*t1*vca*vtf+2*ro*
     1   *2*vca*vtf+8*vcf**2-16*ro*t1**2*vca*vcf-16*t1**2*vca*vcf+16*ro*
     2   t1*vca*vcf+16*t1*vca*vcf-4*ro*vca*vcf-14*vca*vcf+8*ro*t1**2*vca
     3   **2+8*t1**2*vca**2-8*ro*t1*vca**2-8*t1*vca**2+2*ro*vca**2+5*vca
     4   **2)*zg**6/(s*(t1-1)*t1*vda)+dd
      tmp0 = 16*ro*t1**4*vcf-64*t1**4*vcf-32*ro*t1**3*vcf+128*t1**3*vcf+
     1   16*ro**2*t1**2*vcf-16*ro*t1**2*vcf-72*t1**2*vcf-16*ro**2*t1*vcf
     2   +32*ro*t1*vcf+8*t1*vcf+4*ro**3*vcf-8*ro**2*vcf+48*ro*t1**6*vca-
     3   80*t1**6*vca-144*ro*t1**5*vca+240*t1**5*vca+8*ro**2*t1**4*vca+1
     4   00*ro*t1**4*vca-196*t1**4*vca-16*ro**2*t1**3*vca+40*ro*t1**3*vc
     5   a-8*t1**3*vca+4*ro**3*t1**2*vca-10*ro**2*t1**2*vca-31*ro*t1**2*
     6   vca+54*t1**2*vca-4*ro**3*t1*vca+18*ro**2*t1*vca-13*ro*t1*vca-10
     7   *t1*vca-2*ro**3*vca+4*ro**2*vca
      tmp0 = pi*tmp0*(2*vcf-vca)*vlpm*vtf*zg**6/(b*s*(t1-1)**2*t1**2*vda
     1   )/2.0E+0
      dd = tmp0+dd
      tmp0 = 4*nlf*vtf-11*vca
      tmp0 = (-2.0E+0)*pi*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-
     1   4*t1+ro**2)*tmp0*(vcf+t1**2*vca-t1*vca)*vlm2*vtf*zg**6/(3.0E+0*
     2   s*(t1-1)**2*t1**2*vda)
      dd = tmp0+dd
      tmp0 = 16*t1**4*vcf-32*t1**3*vcf+24*ro*t1**2*vcf+16*t1**2*vcf-24*r
     1   o*t1*vcf+6*ro**2*vcf-8*t1**4*vca+16*t1**3*vca-12*ro*t1**2*vca-4
     2   *t1**2*vca+12*ro*t1*vca-8*t1*vca-3*ro**2*vca+4*vca
      tmp0 = -auinv*pi*(4*t1+ro-4)*tmp0*(2*vcf+t1*vca-vca)*vdw*vtf*zg**6
     1   /((t1-1)**3*t1*vda)/8.0E+0
      dd = tmp0+dd
      tmp0 = 16*ro*t1**3*vcf**2
      tmp0 = 32*t1**3*vcf**2+8*ro**2*t1**2*vcf**2-16*ro*t1**2*vcf**2-96*
     1   t1**2*vcf**2-16*ro**2*t1*vcf**2+8*ro*t1*vcf**2+96*t1*vcf**2+12*
     2   ro**2*vcf**2-8*ro*vcf**2-32*vcf**2+8*ro*t1**4*vca*vcf-8*t1**4*v
     3   ca*vcf+4*ro**2*t1**3*vca*vcf-32*ro*t1**3*vca*vcf-16*t1**3*vca*v
     4   cf-16*ro**2*t1**2*vca*vcf+36*ro*t1**2*vca*vcf+92*t1**2*vca*vcf+
     5   22*ro**2*t1*vca*vcf-20*ro*t1*vca*vcf-104*t1*vca*vcf-12*ro**2*vc
     6   a*vcf+8*ro*vca*vcf+36*vca*vcf-16*t1**5*vca**2-8*ro*t1**4*vca**2
     7   +56*t1**4*vca**2-2*ro**2*t1**3*vca**2+20*ro*t1**3*vca**2-62*t1*
     8   *3*vca**2+6*ro**2*t1**2*vca**2-18*ro*t1**2*vca**2+10*t1**2*vca*
     9   *2-7*ro**2*t1*vca**2+8*ro*t1*vca**2+22*t1*vca**2+3*ro**2*vca**2
     :   -2*ro*vca**2-10*vca**2+tmp0
      tmp0 = pi*tmp0*vdw*vtf*zg**6/(s*(t1-1)**3*t1*vda)
      dd = tmp0+dd
      tmp0 = 16*t1**4*vcf-32*t1**3*vcf+24*ro*t1**2*vcf+16*t1**2*vcf-24*r
     1   o*t1*vcf+6*ro**2*vcf-8*t1**4*vca+16*t1**3*vca-12*ro*t1**2*vca-4
     2   *t1**2*vca+12*ro*t1*vca-3*ro**2*vca
      tmp0 = atinv*pi*(4*t1-ro)*tmp0*(2*vcf-t1*vca)*vdt*vtf*zg**6/((t1-1
     1   )*t1**3*vda)/8.0E+0
      dd = tmp0+dd
      dd = dd-pi*(16*ro*t1**3*vcf**2+32*t1**3*vcf**2-8*ro**2*t1**2*vcf**
     1   2-32*ro*t1**2*vcf**2+24*ro*t1*vcf**2-4*ro**2*vcf**2-8*ro*t1**4*
     2   vca*vcf+8*t1**4*vca*vcf+4*ro**2*t1**3*vca*vcf-48*t1**3*vca*vcf+
     3   4*ro**2*t1**2*vca*vcf+12*ro*t1**2*vca*vcf+4*t1**2*vca*vcf+2*ro*
     4   *2*t1*vca*vcf-12*ro*t1*vca*vcf+2*ro**2*vca*vcf-16*t1**5*vca**2+
     5   8*ro*t1**4*vca**2+24*t1**4*vca**2-2*ro**2*t1**3*vca**2-12*ro*t1
     6   **3*vca**2+2*t1**3*vca**2+6*ro*t1**2*vca**2-ro**2*t1*vca**2)*vd
     7   t*vtf*zg**6/(s*(t1-1)*t1**3*vda)
      tmp0 = 32*t1**4*vcf**2-64*t1**3*vcf**2-24*ro**2*t1**2*vcf**2+8*ro*
     1   t1**2*vcf**2
      tmp0 = 80*t1**2*vcf**2+24*ro**2*t1*vcf**2-8*ro*t1*vcf**2-48*t1*vcf
     1   **2-4*ro**3*vcf**2+8*ro**2*vcf**2+8*ro**2*t1**4*vca*vcf+8*ro*t1
     2   **4*vca*vcf-56*t1**4*vca*vcf-16*ro**2*t1**3*vca*vcf-16*ro*t1**3
     3   *vca*vcf+112*t1**3*vca*vcf-4*ro**3*t1**2*vca*vcf+48*ro**2*t1**2
     4   *vca*vcf-120*t1**2*vca*vcf+4*ro**3*t1*vca*vcf-40*ro**2*t1*vca*v
     5   cf+8*ro*t1*vca*vcf+64*t1*vca*vcf+4*ro**3*vca*vcf-8*ro**2*vca*vc
     6   f-64*t1**6*vca**2+192*t1**5*vca**2-4*ro**2*t1**4*vca**2-5*ro*t1
     7   **4*vca**2-226*t1**4*vca**2+8*ro**2*t1**3*vca**2+10*ro*t1**3*vc
     8   a**2+132*t1**3*vca**2+2*ro**3*t1**2*vca**2-18*ro**2*t1**2*vca**
     9   2-3*ro*t1**2*vca**2-22*t1**2*vca**2-2*ro**3*t1*vca**2+14*ro**2*
     :   t1*vca**2-2*ro*t1*vca**2-12*t1*vca**2-ro**3*vca**2+2*ro**2*vca*
     ;   *2+tmp0
      tmp0 = pi*tmp0*vdmp*vtf*zg**6/(b*s*(t1-1)**2*t1**2*vda)
      dd = tmp0+dd
      dd = dd-pi*(16*ro*t1**4-48*t1**4-32*ro*t1**3+96*t1**3-7*ro**2*t1**
     1   2+49*ro*t1**2-98*t1**2+7*ro**2*t1-33*ro*t1+50*t1+ro**3-3*ro**2+
     2   6*ro-8)*vca**2*vdmp*vtf*zg**6/(b**3*s*(t1-1)*t1*vda)
      dd = dd-pi*(2*t1-1)**2*(4*ro*t1**2+8*t1**2-4*ro*t1-8*t1+3*ro**2-8*
     1   ro+8)*vca**2*vdmp*vtf*zg**6/(b**5*s*(t1-1)*t1*vda)/2.0E+0
      tmp0 = 2*vcf+2*t1**2*vca-2*t1*vca-vca
      tmp0 = -pi*(ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*
     1   t1+ro**2)*tmp0*(2*vcf-vca)*vdmb*vtf*zg**6/(b*s*(t1-1)**2*t1**2*
     2   vda)
      dd = tmp0+dd
      tmp0 = 128*ro*t1**4*vcf**2-224*t1**4*vcf**2-256*ro*t1**3*vcf**2+44
     1   8*t1**3*vcf**2+40*ro**2*t1**2*vcf**2+72*ro*t1**2*vcf**2-304*t1*
     2   *2*vcf**2-40*ro**2*t1*vcf**2+56*ro*t1*vcf**2+80*t1*vcf**2+12*ro
     3   **3*vcf**2-24*ro**2*vcf**2+128*ro*t1**6*vca*vcf
      tmp0 = -256*t1**6*vca*vcf-384*ro*t1**5*vca*vcf+768*t1**5*vca*vcf+7
     1   2*ro**2*t1**4*vca*vcf+200*ro*t1**4*vca*vcf-696*t1**4*vca*vcf-14
     2   4*ro**2*t1**3*vca*vcf+240*ro*t1**3*vca*vcf+112*t1**3*vca*vcf+12
     3   *ro**3*t1**2*vca*vcf+16*ro**2*t1**2*vca*vcf-128*ro*t1**2*vca*vc
     4   f+136*t1**2*vca*vcf-12*ro**3*t1*vca*vcf+56*ro**2*t1*vca*vcf-56*
     5   ro*t1*vca*vcf-64*t1*vca*vcf-12*ro**3*vca*vcf+24*ro**2*vca*vcf-6
     6   4*ro*t1**6*vca**2+64*t1**6*vca**2+192*ro*t1**5*vca**2-192*t1**5
     7   *vca**2-36*ro**2*t1**4*vca**2-133*ro*t1**4*vca**2+158*t1**4*vca
     8   **2+72*ro**2*t1**3*vca**2-54*ro*t1**3*vca**2+4*t1**3*vca**2-6*r
     9   o**3*t1**2*vca**2-18*ro**2*t1**2*vca**2+45*ro*t1**2*vca**2-54*t
     :   1**2*vca**2+6*ro**3*t1*vca**2-18*ro**2*t1*vca**2+14*ro*t1*vca**
     ;   2+20*t1*vca**2+3*ro**3*vca**2-6*ro**2*vca**2+tmp0
      tmp0 = pi**3*tmp0*vtf*zg**6/(b*s*(t1-1)**2*t1**2*vda)/1.2E+1
      dd = tmp0+dd
      dd = pi*(ro-1)*(2*t1-1)**2*vca*(2*vcf-vca)*vtf*zg**6/(b**2*s*(t1-1
     1   )*t1*vda)+dd
      dd = dd-pi**3*(16*ro*t1**4-48*t1**4-32*ro*t1**3+96*t1**3-7*ro**2*t
     1   1**2+49*ro*t1**2-98*t1**2+7*ro**2*t1-33*ro*t1+50*t1+ro**3-3*ro*
     2   *2+6*ro-8)*vca**2*vtf*zg**6/(b**3*s*(t1-1)*t1*vda)/1.2E+1
      dd = pi*(ro-1)*(2*t1-1)**4*vca**2*vtf*zg**6/(b**4*s*(t1-1)*t1*vda)
     1   +dd
      dd = dd-pi**3*(2*t1-1)**2*(4*ro*t1**2+8*t1**2-4*ro*t1-8*t1+3*ro**2
     1   -8*ro+8)*vca**2*vtf*zg**6/(b**5*s*(t1-1)*t1*vda)/2.4E+1
      tmp0 = 64*pi**2*t1**5*vcf+192*t1**5*vcf+16*pi**2*ro*t1**4*vcf-192*
     1   pi**2*t1**4*vcf+192*t1**4*vcf+64*pi**2*ro*t1**3*vcf
      tmp0 = 240*ro*t1**3*vcf+192*pi**2*t1**3*vcf-960*t1**3*vcf+24*pi**2
     1   *ro**2*t1**2*vcf-24*ro**2*t1**2*vcf-176*pi**2*ro*t1**2*vcf+288*
     2   ro*t1**2*vcf-64*pi**2*t1**2*vcf-192*t1**2*vcf+192*ro**2*t1*vcf+
     3   96*pi**2*ro*t1*vcf-1296*ro*t1*vcf+1536*t1*vcf+6*pi**2*ro**3*vcf
     4   -24*pi**2*ro**2*vcf-168*ro**2*vcf+768*ro*vcf-768*vcf-32*pi**2*t
     5   1**5*vca-8*pi**2*ro*t1**4*vca+96*pi**2*t1**4*vca-1152*t1**4*vca
     6   -32*pi**2*ro*t1**3*vca-288*ro*t1**3*vca-96*pi**2*t1**3*vca+3456
     7   *t1**3*vca-12*pi**2*ro**2*t1**2*vca+88*pi**2*ro*t1**2*vca+288*r
     8   o*t1**2*vca+32*pi**2*t1**2*vca-3456*t1**2*vca-72*ro**2*t1*vca-4
     9   8*pi**2*ro*t1*vca+288*ro*t1*vca+1152*t1*vca-3*pi**2*ro**3*vca+1
     :   2*pi**2*ro**2*vca+72*ro**2*vca-288*ro*vca+tmp0
      tmp0 = -auinv*pi*tmp0*(2*vcf+t1*vca-vca)*vtf*zg**6/((t1-1)**3*t1*v
     1   da)/4.8E+1
      dd = tmp0+dd
      tmp0 = 64*pi**2*t1**5*vcf+192*t1**5*vcf-16*pi**2*ro*t1**4*vcf-128*
     1   pi**2*t1**4*vcf-1152*t1**4*vcf+128*pi**2*ro*t1**3*vcf+240*ro*t1
     2   **3*vcf+64*pi**2*t1**3*vcf+1728*t1**3*vcf-24*pi**2*ro**2*t1**2*
     3   vcf+24*ro**2*t1**2*vcf-112*pi**2*ro*t1**2*vcf-1008*ro*t1**2*vcf
     4   +48*pi**2*ro**2*t1*vcf+144*ro**2*t1*vcf-6*pi**2*ro**3*vcf-32*pi
     5   **2*t1**5*vca+8*pi**2*ro*t1**4*vca+64*pi**2*t1**4*vca+1152*t1**
     6   4*vca-64*pi**2*ro*t1**3*vca-288*ro*t1**3*vca-32*pi**2*t1**3*vca
     7   -1152*t1**3*vca+12*pi**2*ro**2*t1**2*vca+56*pi**2*ro*t1**2*vca+
     8   576*ro*t1**2*vca-24*pi**2*ro**2*t1*vca-72*ro**2*t1*vca+3*pi**2*
     9   ro**3*vca
      tmp0 = atinv*pi*tmp0*(2*vcf-t1*vca)*vtf*zg**6/((t1-1)*t1**3*vda)/4
     1   .8E+1
      dd = tmp0+dd
      tmp0 = 256*nlf*t1**6*vcf*vtf-768*nlf*t1**5*vcf*vtf+256*nlf*ro*t1**
     1   4*vcf*vtf+768*nlf*t1**4*vcf*vtf-512*nlf*ro*t1**3*vcf*vtf-256*nl
     2   f*t1**3*vcf*vtf+64*nlf*ro**2*t1**2*vcf*vtf+256*nlf*ro*t1**2*vcf
     3   *vtf
      tmp0 = -64*nlf*ro**2*t1*vcf*vtf+256*nlf*t1**8*vca*vtf-1024*nlf*t1*
     1   *7*vca*vtf-24*pi**2*ro**2*t1**6*vca*vtf+192*ro**2*t1**6*vca*vtf
     2   +288*nlf*ro*t1**6*vca*vtf+32*ro*t1**6*vca*vtf+1536*nlf*t1**6*vc
     3   a*vtf+72*pi**2*ro**2*t1**5*vca*vtf-576*ro**2*t1**5*vca*vtf-864*
     4   nlf*ro*t1**5*vca*vtf-96*ro*t1**5*vca*vtf-1024*nlf*t1**5*vca*vtf
     5   -78*pi**2*ro**2*t1**4*vca*vtf+64*nlf*ro**2*t1**4*vca*vtf+624*ro
     6   **2*t1**4*vca*vtf+872*nlf*ro*t1**4*vca*vtf+104*ro*t1**4*vca*vtf
     7   +256*nlf*t1**4*vca*vtf+36*pi**2*ro**2*t1**3*vca*vtf-128*nlf*ro*
     8   *2*t1**3*vca*vtf-288*ro**2*t1**3*vca*vtf-304*nlf*ro*t1**3*vca*v
     9   tf-48*ro*t1**3*vca*vtf-6*pi**2*ro**2*t1**2*vca*vtf+64*nlf*ro**2
     :   *t1**2*vca*vtf+48*ro**2*t1**2*vca*vtf+8*nlf*ro*t1**2*vca*vtf+8*
     ;   ro*t1**2*vca*vtf-96*pi**2*t1**6*vcf**2-1152*t1**6*vcf**2+288*pi
     <   **2*t1**5*vcf**2+3456*t1**5*vcf**2-8*pi**2*ro**2*t1**4*vcf**2-9
     =   6*pi**2*ro*t1**4*vcf**2-576*ro*t1**4*vcf**2-368*pi**2*t1**4*vcf
     >   **2-4992*t1**4*vcf**2+16*pi**2*ro**2*t1**3*vcf**2+tmp0
      tmp0 = 192*pi**2*ro*t1**3*vcf**2+1152*ro*t1**3*vcf**2+256*pi**2*t1
     1   **3*vcf**2+4224*t1**3*vcf**2-24*pi**2*ro**2*t1**2*vcf**2-192*ro
     2   **2*t1**2*vcf**2-144*pi**2*ro*t1**2*vcf**2-768*ro*t1**2*vcf**2-
     3   80*pi**2*t1**2*vcf**2-1536*t1**2*vcf**2+16*pi**2*ro**2*t1*vcf**
     4   2+192*ro**2*t1*vcf**2+48*pi**2*ro*t1*vcf**2+192*ro*t1*vcf**2-8*
     5   pi**2*ro**2*vcf**2-1152*t1**8*vca*vcf+4608*t1**7*vca*vcf+16*pi*
     6   *2*ro*t1**6*vca*vcf-768*ro*t1**6*vca*vcf+264*pi**2*t1**6*vca*vc
     7   f-8864*t1**6*vca*vcf-48*pi**2*ro*t1**5*vca*vcf+2304*ro*t1**5*vc
     8   a*vcf-792*pi**2*t1**5*vca*vcf+10464*t1**5*vca*vcf+42*pi**2*ro**
     9   2*t1**4*vca*vcf-192*ro**2*t1**4*vca*vcf+124*pi**2*ro*t1**4*vca*
     :   vcf-3920*ro*t1**4*vca*vcf+968*pi**2*t1**4*vca*vcf-6264*t1**4*vc
     ;   a*vcf-84*pi**2*ro**2*t1**3*vca*vcf+384*ro**2*t1**3*vca*vcf-168*
     <   pi**2*ro*t1**3*vca*vcf+4000*ro*t1**3*vca*vcf-616*pi**2*t1**3*vc
     =   a*vcf+464*t1**3*vca*vcf+78*pi**2*ro**2*t1**2*vca*vcf-512*ro**2*
     >   t1**2*vca*vcf+tmp0
      tmp0 = 100*pi**2*ro*t1**2*vca*vcf-1520*ro*t1**2*vca*vcf+176*pi**2*
     1   t1**2*vca*vcf+744*t1**2*vca*vcf-36*pi**2*ro**2*t1*vca*vcf+320*r
     2   o**2*t1*vca*vcf-24*pi**2*ro*t1*vca*vcf-96*ro*t1*vca*vcf+4*pi**2
     3   *ro**2*vca*vcf+160*pi**2*t1**8*vca**2-1088*t1**8*vca**2-640*pi*
     4   *2*t1**7*vca**2+4352*t1**7*vca**2+64*pi**2*ro*t1**6*vca**2-1296
     5   *ro*t1**6*vca**2+1004*pi**2*t1**6*vca**2-5952*t1**6*vca**2-192*
     6   pi**2*ro*t1**5*vca**2+3888*ro*t1**5*vca**2-772*pi**2*t1**5*vca*
     7   *2+2624*t1**5*vca**2+pi**2*ro**2*t1**4*vca**2-320*ro**2*t1**4*v
     8   ca**2+214*pi**2*ro*t1**4*vca**2-3772*ro*t1**4*vca**2+260*pi**2*
     9   t1**4*vca**2+640*t1**4*vca**2-2*pi**2*ro**2*t1**3*vca**2+640*ro
     :   **2*t1**3*vca**2-108*pi**2*ro*t1**3*vca**2+1064*ro*t1**3*vca**2
     ;   +20*pi**2*t1**3*vca**2-576*t1**3*vca**2-3*pi**2*ro**2*t1**2*vca
     <   **2-320*ro**2*t1**2*vca**2+22*pi**2*ro*t1**2*vca**2+116*ro*t1**
     =   2*vca**2-32*pi**2*t1**2*vca**2+4*pi**2*ro**2*t1*vca**2+tmp0
      tmp0 = -pi*tmp0*vtf*zg**6/(s*(t1-1)**3*t1**3*vda)/1.2E+1
      dd = tmp0+dd
      ggqq = ss+dd/pi/4.0E+0
      return
      end
