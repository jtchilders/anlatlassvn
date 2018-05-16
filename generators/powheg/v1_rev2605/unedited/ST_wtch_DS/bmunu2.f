c     list of query replaces:
c      mup - nu
c     LorentzIndex(mu) - mu
c     LorentzIndex(nu) - nu
c     Pair(mu,nu) - metric(mu,nu)
c     Momentum(k1) - k1
c     Momentum(k2) - k2
c     Momentum(p1) - p1
c     Momentum(p2) - p2


      function born_bmunu(mu,nu,p1,p2,k1,k2,mt,mw,s,t,u)
c     it returns, as a number, the value of bmunu for singletop wt processes
c     p1 is incoming quark momentum
c     p2 is gluon momentum 
c     k1 is W momentum
c     k2 is top momentum
      implicit none
      integer mu,nu
      real *8 p1(0:3),p2(0:3),k1(0:3),k2(0:3),mt,mw,s,t,u
      real *8 born_bmunu
      real *8 bmunutmp

      real *8 Pair
      external Pair

      real *8 metric(0:3,0:3)
      data metric/1d0, 0d0, 0d0, 0d0,
     #            0d0,-1d0, 0d0, 0d0,
     #            0d0, 0d0,-1d0, 0d0,
     #            0d0, 0d0, 0d0,-1d0/
      save metric


      bmunutmp=(5*mt**4*metric(mu,nu))/(mt**2 - t)**2 + 
     -  (mt**2*mw**2*metric(mu,nu))/(mt**2 - t)**2 - 
     -  (2*mw**4*metric(mu,nu))/(mt**2 - t)**2 - 
     -  (mt**4*mw**2*metric(mu,nu))/(s*(mt**2 - t)**2) - 
     -  (mt**2*s*metric(mu,nu))/(mt**2 - t)**2 + 
     -  (2*mw**2*s*metric(mu,nu))/(mt**2 - t)**2 - 
     -  (6*mt**2*t*metric(mu,nu))/(mt**2 - t)**2 - 
     -  (3*mt**4*t*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (2*t**2*metric(mu,nu))/(mt**2 - t)**2 + 
     -  (7*mt**2*t**2*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (mt**4*t**2*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) + 
     -  (mw**2*t**2*metric(mu,nu))/(s*(mt**2 - t)**2) - 
     -  (s*t**2*metric(mu,nu))/(mw**2*(mt**2 - t)**2) - 
     -  (4*t**3*metric(mu,nu))/(s*(mt**2 - t)**2) - 
     -  (2*mt**2*t**3*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) + 
     -  (t**4*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) - 
     -  (5*mt**2*u*metric(mu,nu))/(mt**2 - t)**2 - 
     -  (mt**4*u*metric(mu,nu))/(mw**2*(mt**2 - t)**2) + 
     -  (2*mw**2*u*metric(mu,nu))/(mt**2 - t)**2 + 
     -  (3*mt**4*u*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (mt**6*u*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*mt**2*mw**2*u*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (4*t*u*metric(mu,nu))/(mt**2 - t)**2 + 
     -  (2*mt**2*t*u*metric(mu,nu))/(mw**2*(mt**2 - t)**2) - 
     -  (5*mt**2*t*u*metric(mu,nu))/(s*(mt**2 - t)**2) - 
     -  (3*mt**4*t*u*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) - 
     -  (2*mw**2*t*u*metric(mu,nu))/(s*(mt**2 - t)**2) - 
     -  (2*t**2*u*metric(mu,nu))/(mw**2*(mt**2 - t)**2) + 
     -  (2*t**2*u*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (2*mt**2*t**2*u*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) - 
     -  (2*mt**2*u**2*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (2*t*u**2*metric(mu,nu))/(s*(mt**2 - t)**2) + 
     -  (mt**2*t*u**2*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) - 
     -  (t**2*u**2*metric(mu,nu))/(mw**2*s*(mt**2 - t)**2) - 
     -  (4*mt**4*Pair(mu,k1)*Pair(nu,k1))/
     -   (mw**2*(mt**2 - t)**2) + (4*mt**2*t*Pair(mu,k1)*
     -     Pair(nu,k1))/(mw**2*(mt**2 - t)**2) + 
     -  (2*mt**2*Pair(mu,k2)*Pair(nu,k1))/
     -   (mt**2 - t)**2 - (4*mw**2*Pair(mu,k2)*
     -     Pair(nu,k1))/(mt**2 - t)**2 + 
     -  (2*mt**2*t*Pair(mu,k2)*Pair(nu,k1))/
     -   (mw**2*(mt**2 - t)**2) + (2*mt**4*Pair(mu,p1)*
     -     Pair(nu,k1))/(mw**2*(mt**2 - t)**2) - 
     -  (4*mt**4*Pair(mu,p1)*Pair(nu,k1))/
     -   (s*(mt**2 - t)**2) - (2*mt**6*Pair(mu,p1)*
     -     Pair(nu,k1))/(mw**2*s*(mt**2 - t)**2) - 
     -  (4*mt**2*t*Pair(mu,p1)*Pair(nu,k1))/
     -   (mw**2*(mt**2 - t)**2) + (10*mt**2*t*Pair(mu,p1)*
     -     Pair(nu,k1))/(s*(mt**2 - t)**2) + 
     -  (4*mt**4*t*Pair(mu,p1)*Pair(nu,k1))/
     -   (mw**2*s*(mt**2 - t)**2) + (2*t**2*Pair(mu,p1)*
     -     Pair(nu,k1))/(mw**2*(mt**2 - t)**2) - 
     -  (6*t**2*Pair(mu,p1)*Pair(nu,k1))/
     -   (s*(mt**2 - t)**2) - (4*mt**2*t**2*Pair(mu,p1)*
     -     Pair(nu,k1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*t**3*Pair(mu,p1)*Pair(nu,k1))/
     -   (mw**2*s*(mt**2 - t)**2) + (4*mt**2*u*Pair(mu,p1)*
     -     Pair(nu,k1))/(s*(mt**2 - t)**2) - 
     -  (2*mt**4*u*Pair(mu,p1)*Pair(nu,k1))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*t*u*Pair(mu,p1)*
     -     Pair(nu,k1))/(s*(mt**2 - t)**2) + 
     -  (2*t**2*u*Pair(mu,p1)*Pair(nu,k1))/
     -   (mw**2*s*(mt**2 - t)**2) + (2*mt**4*Pair(mu,p2)*
     -     Pair(nu,k1))/(s*(mt**2 - t)**2) - 
     -  (2*mt**2*t*Pair(mu,p2)*Pair(nu,k1))/
     -   (s*(mt**2 - t)**2) - (2*mt**4*t*Pair(mu,p2)*
     -     Pair(nu,k1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*mt**2*t**2*Pair(mu,p2)*Pair(nu,k1))/
     -   (mw**2*s*(mt**2 - t)**2) + (2*mt**2*Pair(mu,k1)*
     -     Pair(nu,k2))/(mt**2 - t)**2 - 
     -  (4*mw**2*Pair(mu,k1)*Pair(nu,k2))/
     -   (mt**2 - t)**2 + (2*mt**2*t*Pair(mu,k1)*
     -     Pair(nu,k2))/(mw**2*(mt**2 - t)**2) - 
     -  (6*mt**2*Pair(mu,p1)*Pair(nu,k2))/
     -   (mt**2 - t)**2 + (4*mw**2*Pair(mu,p1)*
     -     Pair(nu,k2))/(mt**2 - t)**2 + 
     -  (2*mt**4*Pair(mu,p1)*Pair(nu,k2))/
     -   (s*(mt**2 - t)**2) - (4*mt**2*mw**2*Pair(mu,p1)*
     -     Pair(nu,k2))/(s*(mt**2 - t)**2) + 
     -  (4*t*Pair(mu,p1)*Pair(nu,k2))/(mt**2 - t)**2 - 
     -  (4*mt**2*t*Pair(mu,p1)*Pair(nu,k2))/
     -   (s*(mt**2 - t)**2) + (4*mw**2*t*Pair(mu,p1)*
     -     Pair(nu,k2))/(s*(mt**2 - t)**2) - 
     -  (2*t**2*Pair(mu,p1)*Pair(nu,k2))/
     -   (mw**2*(mt**2 - t)**2) + (2*t**2*Pair(mu,p1)*
     -     Pair(nu,k2))/(s*(mt**2 - t)**2) + 
     -  (2*mt**2*t**2*Pair(mu,p1)*Pair(nu,k2))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*t**3*Pair(mu,p1)*
     -     Pair(nu,k2))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*mt**2*t*u*Pair(mu,p1)*Pair(nu,k2))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*t**2*u*Pair(mu,p1)*
     -     Pair(nu,k2))/(mw**2*s*(mt**2 - t)**2) - 
     -  (4*mt**2*mw**2*Pair(mu,p2)*Pair(nu,k2))/
     -   (s*(mt**2 - t)**2) + (4*mt**2*t*Pair(mu,p2)*
     -     Pair(nu,k2))/(s*(mt**2 - t)**2) + 
     -  (4*mw**2*t*Pair(mu,p2)*Pair(nu,k2))/
     -   (s*(mt**2 - t)**2) - (4*t**2*Pair(mu,p2)*
     -     Pair(nu,k2))/(s*(mt**2 - t)**2) + 
     -  (2*mt**4*Pair(mu,k1)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) - (4*mt**4*Pair(mu,k1)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (2*mt**6*Pair(mu,k1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*mt**2*t*Pair(mu,k1)*
     -     Pair(nu,p1))/(mw**2*(mt**2 - t)**2) + 
     -  (10*mt**2*t*Pair(mu,k1)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) + (4*mt**4*t*Pair(mu,k1)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*t**2*Pair(mu,k1)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) - (6*t**2*Pair(mu,k1)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (4*mt**2*t**2*Pair(mu,k1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (2*t**3*Pair(mu,k1)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (4*mt**2*u*Pair(mu,k1)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) - (2*mt**4*u*Pair(mu,k1)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) - 
     -  (4*t*u*Pair(mu,k1)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) + (2*t**2*u*Pair(mu,k1)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) - 
     -  (6*mt**2*Pair(mu,k2)*Pair(nu,p1))/
     -   (mt**2 - t)**2 + (4*mw**2*Pair(mu,k2)*
     -     Pair(nu,p1))/(mt**2 - t)**2 + 
     -  (2*mt**4*Pair(mu,k2)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) - (4*mt**2*mw**2*Pair(mu,k2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (4*t*Pair(mu,k2)*Pair(nu,p1))/(mt**2 - t)**2 - 
     -  (4*mt**2*t*Pair(mu,k2)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) + (4*mw**2*t*Pair(mu,k2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (2*t**2*Pair(mu,k2)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) + (2*t**2*Pair(mu,k2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (2*mt**2*t**2*Pair(mu,k2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*t**3*Pair(mu,k2)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*mt**2*t*u*Pair(mu,k2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*t**2*u*Pair(mu,k2)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (8*mt**2*Pair(mu,p1)*Pair(nu,p1))/
     -   (mt**2 - t)**2 + (8*mt**4*mw**2*Pair(mu,p1)*
     -     Pair(nu,p1))/(s**2*(mt**2 - t)**2) - 
     -  (8*mt**2*mw**2*Pair(mu,p1)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) - (8*t*Pair(mu,p1)*
     -     Pair(nu,p1))/(mt**2 - t)**2 - 
     -  (4*mt**2*t*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) - (4*mt**6*t*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) - 
     -  (16*mt**2*mw**2*t*Pair(mu,p1)*Pair(nu,p1))/
     -   (s**2*(mt**2 - t)**2) + (4*mt**2*t*Pair(mu,p1)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (8*mt**4*t*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (8*mw**2*t*Pair(mu,p1)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (4*t**2*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) + (8*mt**4*t**2*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (8*mw**2*t**2*Pair(mu,p1)*Pair(nu,p1))/
     -   (s**2*(mt**2 - t)**2) - (4*t**2*Pair(mu,p1)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (12*mt**2*t**2*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*mt**2*t**3*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*t**3*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*mt**6*u*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*mt**4*u*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (8*mt**4*t*u*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) - 
     -  (8*mt**2*t*u*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*mt**2*t**2*u*Pair(mu,p1)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*t**2*u*Pair(mu,p1)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (4*mt**2*Pair(mu,p2)*
     -     Pair(nu,p1))/(mt**2 - t)**2 + 
     -  (4*mt**4*mw**2*Pair(mu,p2)*Pair(nu,p1))/
     -   (s**2*(mt**2 - t)**2) - (8*mt**4*Pair(mu,p2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (4*mt**2*mw**2*Pair(mu,p2)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) - (4*t*Pair(mu,p2)*
     -     Pair(nu,p1))/(mt**2 - t)**2 - 
     -  (2*mt**2*t*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) - (2*mt**6*t*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) - 
     -  (8*mt**2*mw**2*t*Pair(mu,p2)*Pair(nu,p1))/
     -   (s**2*(mt**2 - t)**2) + (14*mt**2*t*Pair(mu,p2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (4*mt**4*t*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (4*mw**2*t*Pair(mu,p2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) + 
     -  (2*t**2*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*(mt**2 - t)**2) + (4*mt**4*t**2*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*mw**2*t**2*Pair(mu,p2)*Pair(nu,p1))/
     -   (s**2*(mt**2 - t)**2) - (6*t**2*Pair(mu,p2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (6*mt**2*t**2*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**2*t**3*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (2*t**3*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**6*u*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*mt**2*u*Pair(mu,p2)*Pair(nu,p1))/
     -   (s*(mt**2 - t)**2) + (2*mt**4*u*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s*(mt**2 - t)**2) + 
     -  (4*mt**4*t*u*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s**2*(mt**2 - t)**2) - (4*t*u*Pair(mu,p2)*
     -     Pair(nu,p1))/(s*(mt**2 - t)**2) - 
     -  (4*mt**2*t*u*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**2*t**2*u*Pair(mu,p2)*
     -     Pair(nu,p1))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (2*t**2*u*Pair(mu,p2)*Pair(nu,p1))/
     -   (mw**2*s*(mt**2 - t)**2) + (2*mt**4*Pair(mu,k1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) - 
     -  (2*mt**2*t*Pair(mu,k1)*Pair(nu,p2))/
     -   (s*(mt**2 - t)**2) - (2*mt**4*t*Pair(mu,k1)*
     -     Pair(nu,p2))/(mw**2*s*(mt**2 - t)**2) + 
     -  (2*mt**2*t**2*Pair(mu,k1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2) - (4*mt**2*mw**2*Pair(mu,k2)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) + 
     -  (4*mt**2*t*Pair(mu,k2)*Pair(nu,p2))/
     -   (s*(mt**2 - t)**2) + (4*mw**2*t*Pair(mu,k2)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) - 
     -  (4*t**2*Pair(mu,k2)*Pair(nu,p2))/
     -   (s*(mt**2 - t)**2) + (4*mt**2*Pair(mu,p1)*
     -     Pair(nu,p2))/(mt**2 - t)**2 + 
     -  (4*mt**4*mw**2*Pair(mu,p1)*Pair(nu,p2))/
     -   (s**2*(mt**2 - t)**2) - (8*mt**4*Pair(mu,p1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) - 
     -  (4*mt**2*mw**2*Pair(mu,p1)*Pair(nu,p2))/
     -   (s*(mt**2 - t)**2) - (4*t*Pair(mu,p1)*
     -     Pair(nu,p2))/(mt**2 - t)**2 - 
     -  (2*mt**2*t*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*(mt**2 - t)**2) - (2*mt**6*t*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s**2*(mt**2 - t)**2) - 
     -  (8*mt**2*mw**2*t*Pair(mu,p1)*Pair(nu,p2))/
     -   (s**2*(mt**2 - t)**2) + (14*mt**2*t*Pair(mu,p1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) + 
     -  (4*mt**4*t*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2) + (4*mw**2*t*Pair(mu,p1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) + 
     -  (2*t**2*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*(mt**2 - t)**2) + (4*mt**4*t**2*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*mw**2*t**2*Pair(mu,p1)*Pair(nu,p2))/
     -   (s**2*(mt**2 - t)**2) - (6*t**2*Pair(mu,p1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) - 
     -  (6*mt**2*t**2*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**2*t**3*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (2*t**3*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**6*u*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (4*mt**2*u*Pair(mu,p1)*Pair(nu,p2))/
     -   (s*(mt**2 - t)**2) + (2*mt**4*u*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s*(mt**2 - t)**2) + 
     -  (4*mt**4*t*u*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s**2*(mt**2 - t)**2) - (4*t*u*Pair(mu,p1)*
     -     Pair(nu,p2))/(s*(mt**2 - t)**2) - 
     -  (4*mt**2*t*u*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2) - (2*mt**2*t**2*u*Pair(mu,p1)*
     -     Pair(nu,p2))/(mw**2*s**2*(mt**2 - t)**2) + 
     -  (2*t**2*u*Pair(mu,p1)*Pair(nu,p2))/
     -   (mw**2*s*(mt**2 - t)**2)






      born_bmunu=bmunutmp
      end


      function Pair(index,lorentzvec)
      implicit none
      integer index
      real *8 lorentzvec(0:3)
      real *8 Pair
      
      Pair=lorentzvec(index)
      return
      end
