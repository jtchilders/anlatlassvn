      subroutine agg(s,t,u,agg2)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,agg2,aijgg
      double precision mbot,mtop,l1,l2
      double complex bdva2, bdva4

      agg2=0d0
c     gg->h/Hg
      if((ih.eq.1).or.(ih.eq.2)) then
        agg2 = abs(bdva2(s,t,u))**2+abs(bdva2(u,s,t))**2+
     &       abs(bdva2(t,u,s))**2+abs(bdva4(s,t,u))**2
c     gg->Ag
      else
          mtop = mfer(1)
          mbot = mfer(2)
          l1 = lambdafer(1)
          l2 = lambdafer(2)

ccccccccccccc    top top
          agg2 = agg2 + aijgg(mtop,mtop,l1,l1,mh,u,s,t)
          agg2 = agg2 + aijgg(mtop,mtop,l1,l1,mh,t,u,s)
          agg2 = agg2 + aijgg(mtop,mtop,l1,l1,mh,s,t,u)

ccccccccccccc    bottom bottom
          agg2 = agg2 + aijgg(mbot,mbot,l2,l2,mh,u,s,t)
          agg2 = agg2 + aijgg(mbot,mbot,l2,l2,mh,t,u,s)
          agg2 = agg2 + aijgg(mbot,mbot,l2,l2,mh,s,t,u)

ccccccccccccc    top bottom
          agg2 = agg2 + 2d0*aijgg(mtop,mbot,l1,l2,mh,u,s,t)
          agg2 = agg2 + 2d0*aijgg(mtop,mbot,l1,l2,mh,t,u,s)
          agg2 = agg2 + 2d0*aijgg(mtop,mbot,l1,l2,mh,s,t,u)

ccccccccccccc tf**2
          agg2 = agg2*trfer(1)**2

      endif

      end

      subroutine aqqbar(s,aqqbar2)
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      double precision aqqbar2
      double complex ampl
      double precision s,t,u,realfac
      double precision m0,m12,y0,y12
      double precision s0,t0,u0,s12,t12,u12
      double precision mtop,mbot,lbot,ltop
      double precision yt,yb,st,sb
      double complex bdvd12,bdvd0,fermion,scalar,bdsvh2
      external bdsvh2
      integer i

      ampl = dcmplx(0d0)

      if((ih.eq.1).or.(ih.eq.2)) then
      do i=1,afer
         m12 = mfer(i)
         s12 = s/m12**2
         y12 = m12**2/mh2
         fermion = lambdafer(i)*trfer(i)*y12*bdvd12(s12,y12)
         ampl = ampl + fermion
      end do

      if (flg_lhscalars.ne.1) then
       do i=1,asca
        m0 = msca(i)
        s0 = s/m0**2
        y0 = m0**2/mh2
        scalar = lambdasca(i)*(mmaa**2/m0**2)*trsca(i)*y0*bdvd0(s0,y0)
        ampl = ampl + scalar
       end do
      else
         realfac = -(mh2-s)/(2*(mh2))
         do i=1,asca
            m0 = msca(i)
            y0 = m0**2/mh2
            scalar = realfac*lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &           (-1d0/3d0)
            ampl = ampl + scalar
         end do
      endif
c     Pseudoscalar
      else
          mtop = mfer(1)
          mbot = mfer(2)
          ltop = lambdafer(1)
          lbot = lambdafer(2)
          yt=mtop**2/mh**2
          yb=mbot**2/mh**2
          st=s/mtop**2
          sb=s/mbot**2
          ampl = ltop*yt*bdsvh2(st,yt) + lbot*yb*bdsvh2(sb,yb)

          ampl = trfer(1) * ampl * 2
      endif

      aqqbar2 = abs(ampl)**2

      end

      double complex function bdva2(s,t,u)
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      double precision s,t,u,m0,m12,realfac
      double precision y0,y12
      double precision s0,t0,u0,s12,t12,u12
      double complex bdvb12,bdvb0,fermion,scalar
      external bdvb12,bdvb0
      integer i

      bdva2 = dcmplx(0d0)
      do i=1,afer
         m12 = mfer(i)
         s12 = s/m12**2
         t12 = t/m12**2
         u12 = u/m12**2
         y12 = m12**2/mh2

         fermion = lambdafer(i)*trfer(i)*y12**2*
     &        (bdvb12(s12,t12,u12,y12)+bdvb12(s12,u12,t12,y12))

         bdva2 = bdva2 + fermion
      end do

      if (flg_lhscalars.ne.1) then
         do i=1,asca
            m0 = msca(i)
            s0 = s/m0**2
            t0 = t/m0**2
            u0 = u/m0**2
            y0 = m0**2/mh2
            scalar = lambdasca(i)*(mmaa**2/m0**2)*trsca(i)*y0**2*
     &           (bdvb0(s0,t0,u0,y0)+bdvb0(s0,u0,t0,y0))

            bdva2 = bdva2 + scalar
         end do
      else
         realfac = s**2/(4*(s+t+u)**2)
         do i=1,asca
            m0 = msca(i)
            y0 = m0**2/mh2
            scalar = realfac*lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &           (-1d0/3d0)
            bdva2 = bdva2 + scalar
         end do
      endif

      return
      end

      double complex function bdva4(s,t,u)
      implicit none
      include 'PhysPars.h'
      include 'Flags.h'
      double precision s,t,u,m0,m12,y0,y12,realfac
      double precision s0,t0,u0,s12,t12,u12
      double complex bdvc12,bdvc0,fermion,scalar
      external bdvc12,bdvc0
      integer i

      bdva4 = dcmplx(0d0)

      do i=1,afer
      m12 = mfer(i)
      s12 = s/m12**2
      t12 = t/m12**2
      u12 = u/m12**2
      y12 = m12**2/mh2
      fermion = lambdafer(i)*trfer(i)*y12**2*
     &     (bdvc12(s12,t12,u12,y12)+bdvc12(t12,u12,s12,y12)+
     &      bdvc12(u12,s12,t12,y12))

      bdva4 = bdva4 + fermion
      end do


      if (flg_lhscalars.ne.1) then
         do i=1,asca
            m0 = msca(i)
            s0 = s/m0**2
            t0 = t/m0**2
            u0 = u/m0**2
            y0 = m0**2/mh2
            scalar = lambdasca(i)*(mmaa**2/m0**2)*trsca(i)*y0**2*
     &      (bdvc0(s0,t0,u0,y0)+bdvc0(t0,u0,s0,y0)+bdvc0(u0,s0,t0,y0))

            bdva4 = bdva4 + scalar
         end do
      else
         realfac = 1d0/4d0
         do i=1,asca
            m0 = msca(i)
            y0 = m0**2/mh2
            scalar = realfac*lambdasca(i)*trsca(i)*(mmaa/m0)**2*
     &           (-1d0/3d0)
            bdva4 = bdva4 + scalar
         end do
      endif

      return
      end



      double complex function bdvb12(s,t,u,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,y12,m2
      double complex xs,xt,x12
      double complex bdvbigb12, bdvh3, reduced
      external bdvbigb12, bdvh3, reduced

      x12 = reduced(1d0/y12)
      xs = reduced(s)
      xt = reduced(t)
      m2 = y12*mh2

      bdvb12= bdvbigb12(s,t,u,y12) +
     &       s/4d0*(0.5d0*log(x12)**2-0.5d0*log(xs)**2)
     &       -(s/2d0-s**2/(s+u))*(0.5d0*log(x12)**2-0.5d0*log(xt)**2)
     &       -s/8d0*bdvh3(s,u,t,m2)+s/4d0*bdvh3(t,s,u,m2)


      return
      end

      double complex function bdvb0(s,t,u,y0)
      implicit none
      double precision s,t,u,y0
      double complex bdvbigb0
      external bdvbigb0

      bdvb0= -0.5d0*bdvbigb0(s,t,u,y0) 

      return
      end


      double complex function bdvbigb12(s,t,u,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,m12, y12, m2
      double complex xs,xt,x12,r1,r2
      double complex bdvh3, reduced
      external bdvh3, reduced

      x12 = reduced(1d0/y12)
      xs = reduced(s)
      xt = reduced(t)
      m2 = y12*mh2

      if (y12.lt.0.25d0) then
         r1 = dsqrt(1d0-4d0*y12)
      else
         r1 = cdsqrt(dcmplx(1d0-4d0*y12))
      end if
      if ((t.lt.0.d0).or.(t.gt.4d0)) then
         r2 = dsqrt(1d0-4d0/t)
      else
         r2 = cdsqrt(dcmplx(1d0-4d0/t))
      end if


      bdvbigb12= s*(t-s)/(s+t) +
     &     2d0*(t*u**2+2d0*s*t*u)/(s+u)**2*(r1*log(x12)-r2*log(xt))
     &     -(1d0+t*u/s)*0.5d0*log(x12)**2+0.5d0*log(xs)**2
     &     -2d0*(2d0*s**2/(s+u)**2 -1d0-t*u/s)*
     &                0.5d0*(log(x12)**2-log(xt)**2)
     &     +0.5d0*(t*u/s+3d0)*bdvh3(s,u,t,m2)-bdvh3(t,s,u,m2)


      return
      end

      double complex function bdvbigb0(s,t,u,y0)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u, m0,y0,m2
      double complex xs,xt,x0, r1, r2
      double complex bdvh3, reduced
      external bdvh3, reduced

      x0 = reduced(1d0/y0)
      xs = reduced(s)
      xt = reduced(t)
      m2 = y0*mh2

      if (y0.lt.0.25d0) then
         r1 = dsqrt(1d0-4d0*y0)
      else
         r1 = cdsqrt(dcmplx(1d0-4d0*y0))
      end if
      if ((t.lt.0.d0).or.(t.gt.4d0)) then
         r2 = dsqrt(1d0-4d0/t)
      else
         r2 = cdsqrt(dcmplx(1d0-4d0/t))
      end if


      bdvbigb0= s*(t-s)/(s+t) 
     &    +2d0*(t*u**2+2d0*s*t*u)/(s+u)**2*(r1*log(x0)-r2*log(xt))
     &     -(1d0+t*u/s)*0.5d0*log(x0)**2+0.5d0*log(xs)**2
     &     -2d0*(2d0*s**2/(s+u)**2 -1d0-t*u/s)*
     &                0.5d0*(log(x0)**2-log(xt)**2)
     &     +0.5d0*(t*u/s+3d0)*bdvh3(s,u,t,m2)-bdvh3(t,s,u,m2)

      return
      end


      double complex function bdvc12(s,t,u,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,m12, y12, m2
      double complex xs,x12
      double complex bdvbigc12, bdvh3, reduced
      external bdvbigc12, bdvh3, reduced

      x12 = reduced(1d0/y12)
      xs = reduced(s)
      m2 = y12*mh2

      bdvc12= bdvbigc12(s,t,u,y12) +
     &    0.5d0/y12*0.5d0*(log(x12)**2-log(xs)**2)+
     &    0.25d0/y12*bdvh3(s,u,t,1d0)

      return
      end

      double complex function bdvc0(s,t,u, y0)
      implicit none
      double precision s,t,u, y0
      double complex bdvbigc0
      external bdvbigc0

      bdvc0 = -0.5d0*bdvbigc0(s,t,u, y0)

      return
      end



      double complex function bdvbigc12(s,t,u,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,y12,m2
      double complex xs,x12
      double complex bdvh3,reduced
      external bdvh3, reduced

      x12 = reduced(1d0/y12)
      xs = reduced(s)
      m2 = y12*mh2

      bdvbigc12 = -2d0*s-2d0*0.5d0*(log(x12)**2-log(xs)**2)-
     &            bdvh3(u,s,t,1d0)

      return
      end

      double complex function bdvbigc0(s,t,u, y0)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,y0,m2
      double complex xs,x0
      double complex bdvh3, reduced
      external bdvh3, reduced

      x0 = reduced(1d0/y0)
      xs = reduced(s)
      m2 = y0*mh2
      bdvbigc0 = -2d0*s-2d0*0.5d0*(log(x0)**2-log(xs)**2)-
     &            bdvh3(u,s,t,1d0)
      
      return
      end


      double complex function bdvd12(s,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,m12, y12, m2
      double complex xs,x12
      double complex bdvbigd12, reduced
      external bdvbigd12, bdvh3, reduced
      
      x12 = reduced(1d0/y12)
      xs = reduced(s)
      m2 = y12*mh2

      bdvd12= bdvbigd12(s,y12) -
     &        2d0*(0.5d0*(log(x12)**2-log(xs)**2))


      return
      end

      double complex function bdvd0(s, y0)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u, y0
      double complex bdvbigd0
      external bdvbigd0

      bdvd0 = -0.5d0*bdvbigd0(s, y0)

      return
      end



      double complex function bdvbigd12(s,y12)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,y12,m2
      double complex xs,x12,r1,r2
      double complex reduced
      external reduced

      x12 = reduced(1d0/y12)
      xs = reduced(s)
      m2 = y12*mh2

      if (y12.lt.0.25d0) then
         r1 = dsqrt(1d0-4d0*y12)
      else
         r1 = cdsqrt(dcmplx(1d0-4d0*y12))
      end if
      if ((s.lt.0.d0).or.(s.gt.4d0)) then
         r2 = dsqrt(1d0-4d0/s)
      else
         r2 = cdsqrt(dcmplx(1d0-4d0/s))
      end if

      bdvbigd12 = 4d0 + 4d0*s/(1d0/y12-s)*(r1*log(x12)-r2*log(xs))+
     &            8d0/(1d0/y12-s)*0.5d0*(log(x12)**2-log(xs)**2)

      return
      end

      double complex function bdvbigd0(s, y0)
      implicit none
      include 'PhysPars.h'
      double precision s,t,u,y0,m2
      double complex xs,x0,r1,r2
      double complex bdvh3, reduced
      external bdvh3, reduced

      x0 = reduced(1d0/y0)
      xs = reduced(s)
      m2 = y0*mh2

      if (y0.lt.0.25d0) then
         r1 = dsqrt(1d0-4d0*y0)
      else
         r1 = cdsqrt(dcmplx(1d0-4d0*y0))
      end if
      if ((s.lt.0.d0).or.(s.gt.4d0)) then
         r2 = dsqrt(1d0-4d0/s)
      else
         r2 = cdsqrt(dcmplx(1d0-4d0/s))
      end if

      bdvbigd0  = 4d0 + 4d0*s/(1d0/y0-s)*(r1*log(x0)-r2*log(xs))+
     &            8d0/(1d0/y0-s)*0.5d0*(log(x0)**2-log(xs)**2)

      return
      end



      double complex function bdvh3(s,t,u,m2)
      implicit none
      double precision s,t,u,m2
      double complex bdvi3,w3
      external bdvi3,w3

      bdvh3 = -w3(t,s,u,(s+t+u), 1d0)

c      if (1>0) print*,'h3=',bdvh3,s,t,u
c      bdvh3 = dcmplx(0d0)

c      bdvh3 = bdvi3(s,t,u,t) + bdvi3(s,t,u,u) - bdvi3(s,t,u,s+t+u)

      return
      end

      double complex function w3(s,t,u,v,m2)
      implicit none
      double precision s,t,u,v,xf,xh,tau,m2
      double complex ei3
      external ei3

      w3 = ei3(s,t,u,v,m2) -ei3(s,t,u,s,m2) - ei3(s,t,u,u,m2)
c      w3 = dcmplx(0d0)

      return
      end


c      subroutine eei3(s,t,u,varg,mbsq,ei3)
      double complex function ei3(s,t,u,varg,mbsq)
C     ehsv:EqnA.21
      implicit none
c      double complex ei3
      double precision pi
      double precision s,t,u,varg,rat,al,be,ga,r,theta,phi,cosphi
      double precision arg1,arg2,arg3,arg4,ddilog,mbsq
      double complex eli2,zth,zph,CLI2

      external cli2,eli2

      pi = 3.141592653589793d0
      rat=4d0/varg
      if (rat .lt. 0d0) then

           be=0.5d0*(1d0+dsqrt(1d0+4d0*t/(u*s)))
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=(be-ga)/be
           arg4=(be-ga)/(be-1d0)
           ei3=2d0/(2d0*be-1d0)
     .     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     .     +0.5d0*(dlog(be)**2-dlog(be-1d0)**2)
     .     +dlog(ga)*dlog((ga+be-1d0)/be)
     .     +dlog(ga-1d0)*dlog((be-1d0)/(ga+be-1d0)))
      elseif (rat .gt. 1d0) then
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t/(u*s)))
           al=dsqrt(rat-1d0)
           r=dsqrt((al**2+1d0)/(al**2+(2d0*be-1d0)**2))
           cosphi=r*(al**2+2d0*be-1d0)/(1d0+al**2)
           if (cosphi.gt.1) then
              if (cosphi.lt.1.000001) then
                 write(18,*) "Detected rounding error - cosphi: "
                 write(18,*) cosphi
                 cosphi = 1.d0
              else
                 write(18,*) "Rounding error seems too large! Aborting!"
                 write(18,*) cosphi
                 stop
              endif
           endif
           phi=acos(cosphi)
           theta=acos(r*(al**2-2d0*be+1d0)/(1d0+al**2))
           zth=r*dcmplx(cos(theta),sin(theta))
           zph=r*dcmplx(cos(phi),sin(phi))
           ei3=2d0/(2d0*be-1d0)
     .     *(2d0*dble(eli2(zth))-2d0*dble(eli2(zph))
     .     +(phi-theta)*(phi+theta-pi))
      else
           be=0.5d0*(1d0+dsqrt(1d0+4d0*t/(u*s)))
           ga=0.5d0*(1d0+dsqrt(1d0-rat))
           arg1=ga/(ga+be-1d0)
           arg2=(ga-1d0)/(ga+be-1d0)
           arg3=ga/(ga-be)
           arg4=(ga-1d0)/(ga-be)
      
           ei3=2d0/(2d0*be-1d0)
     .     *(-ddilog(arg1)+ddilog(arg2)+ddilog(arg3)-ddilog(arg4)
     .     +dlog(ga/(1d0-ga))*dlog((ga+be-1d0)/(be-ga))
     .     -dcmplx(0d0,1d0)*pi*dlog((ga+be-1d0)/(be-ga)))
      endif

      return
      end

      double precision function aijgg(mi,mj,l1,l2,mh,s,t,u)
      implicit none
      double precision mi,mj,l1,l2,s,t,u,mh, m4
      double precision si,ti,ui,sj,tj,uj,yi,yj
      double precision b1,b2,b3,b4,b5,b6,b7
      double complex h2i,h2j,h3istu,h3jstu,h3itus,h3jtus,h3iust,h3just
      double complex h2ti,h2tj
      double complex atmp,bdsvh2,bdvh3,atmp1,atmp2
      external bdsvh2,bdvh3

      m4=mh**4

      si = s/mi**2
      ti = t/mi**2
      ui = u/mi**2
      yi = mi**2/mh**2

      sj = s/mj**2
      tj = t/mj**2
      uj = u/mj**2
      yj = mj**2/mh**2

      b1=1d0/2d0*( 4d0*t**2*u**2/(t+u)**2 +s**2 -3d0*t*u+
     -             s*(t+u) +(t+u)**2 )

      b2=s**2 +t**2 +u**2 +s*t +2d0*s**2*t*u/(s-t)/(s+u) -
     -   2d0*s*t**2*u/(s-t)/(t+u)

      b3=1d0/8d0*( s**2 +t**2 +u**2 +t*u +s*(t+u) )

      b4=1d0/4d0*(s+t)*(t+u)

      b5=-1d0/2d0*(t**2 +u**2 +s*(t+u) )

      b6=-1d0/2d0*(s**2 +(t+u)*(s+u) +(t-u)*t*u/(t+u) )

      b7=-1d0/2d0*(s**2 +(t+u)*(s+t) -(t-u)*t*u/(t+u) )

      h2i = bdsvh2(si,yi)
      h2j = bdsvh2(sj,yj)
      h2ti = bdsvh2(ti,yi)
      h2tj = bdsvh2(tj,yj)

c      h2i = dconjg(bdsvh2(si,yi))
c      h2j = dconjg(bdsvh2(sj,yj))
c      h2ti = dconjg(bdsvh2(ti,yi))
c      h2tj = dconjg(bdsvh2(tj,yj))


      h3istu = bdvh3(si,ti,ui,mi**2)
      h3jstu = bdvh3(sj,tj,uj,mj**2)
      h3iust = bdvh3(ui,si,ti,mi**2)
      h3just = bdvh3(uj,sj,tj,mj**2)
      h3itus = bdvh3(ti,ui,si,mi**2)
      h3jtus = bdvh3(tj,uj,sj,mj**2)


c      print*,'invar0=',s,t,u
c      print*,'invar=',si,ti,ui,yi,sj,tj,uj,yj
c      print*,'b=',b1/m4,b2/m4,b3/m4,b4/m4,b5/m4,b6/m4,b7/m4
c      print*,h2i,h2j,h2ti,h2tj,h3istu,h3jstu,h3itus,h3jtus,h3iust,h3just

      atmp1 = b1*h2i*dconjg(h2j)+
     -       b2*h2i*dconjg(h2tj)+
     -       b3*h3istu*dconjg(h3jstu)+
     -       b4*h3istu*dconjg(h3just)+
     -       b5*h2i*dconjg(h3jstu)+
     -       b6*h2i*dconjg(h3jtus)+
     -       b7*h2i*dconjg(h3just)

      atmp2 =
     -       b1*h2j*dconjg(h2i)+
     -       b2*h2j*dconjg(h2ti)+
     -       b3*h3jstu*dconjg(h3istu)+
     -       b4*h3jstu*dconjg(h3iust)+
     -       b5*h2j*dconjg(h3istu)+
     -       b6*h2j*dconjg(h3itus)+
     -       b7*h2j*dconjg(h3iust)

c$$$      print*,'b11=',h2i*dconjg(h2j)
c$$$      print*,'b21=',      h2i*dconjg(h2tj)
c$$$      print*,'b31=',      h3istu*dconjg(h3jstu)
c$$$      print*,'b41=',      h3istu*dconjg(h3just)
c$$$      print*,'b51=',      h2i*dconjg(h3jstu)
c$$$      print*,'b61=',      h2i*dconjg(h3jtus)
c$$$      print*,'b71=',      h2i*dconjg(h3just)
c$$$
c$$$
c$$$      print*,'b12=',      b1*h2j*dconjg(h2i)
c$$$      print*,'b22=',      b2*h2j*dconjg(h2ti)
c$$$      print*,'b32=',      b3*h3jstu*dconjg(h3istu)
c$$$      print*,'b42=',      b4*h3jstu*dconjg(h3iust)
c$$$      print*,'b52=',      b5*h2j*dconjg(h3istu)
c$$$      print*,'b62=',      b6*h2j*dconjg(h3itus)
c$$$      print*,'b72=',      b7*h2j*dconjg(h3iust)


      atmp = atmp1+atmp2


      aijgg = 2d0*dble(atmp)

      aijgg=yi*yj/mh**4/4d0*aijgg*l1*l2


c      print*,'aijgg=',aijgg,atmp


      end

      double complex function bdsvh2(s,y)
      implicit none
      double complex mylog,arg1,arg2,arlog1,arlog2
      double precision s,y
      external mylog

      if (abs(s).eq.0d0) then
         arlog1=dcmplx(1d0)
      else
         arg1=dcmplx(1d0-4d0/s)
         arlog1=(sqrt(arg1)-1)/(sqrt(arg1)+1)
      endif

      arg2=dcmplx(1d0-4d0*y)
      arlog2=(sqrt(arg2)-1)/(sqrt(arg2)+1)

c      print*,s,y,arlog1,arlog2

      bdsvh2=1d0/2d0*(log(arlog1)**2-log(arlog2)**2)
c      print*,'---------------------',bdsvh2,s,y,arg1,arlog1,arg2,arlog2
c      bdsvh2=dcmplx(0d0,0d0)

      end

