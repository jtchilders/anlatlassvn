      subroutine gen_real_phsp_remn(xx,xjac)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_rad.h'
      include 'PhysPars.h'
      real * 8 xx(ndiminteg)
      real * 8 m2,xjac,tau,y,beta,cth,s,rs,
     # z,zhigh,zlow,m2inv

*
      real*8 mlep2
      common/leptmass/mlep2
*
      real*8 random
      external random

      real*8 dotp
      external dotp
*
      real*8 pl(0:3),pn(0:3),pg(0:3)
      real*8 betal,betag
      real*8 phi,sph,cph,sth
      real*8 phiph,sphg,cphg,sthg,cthg
      real*8 v(0:3)
      real*8 odv,abig,bbig,cbig,arg
      real*8 p1mod,p1mod_1,p1mod_2,den
      real*8 eps,emax
      integer mu
      integer mode
      real*8 vec(3)
      data vec/0d0,0d0,1d0/
*
      real*8 mpart2,ppart2,ppart,anorm,bnorm,f12,f123
      real*8 ombetal
      real*8 cth12,phi12
      real*8 q1s(0:3),q2s(0:3)
      real*8 q1(0:3),q2(0:3)
*
      mode = 1
*
c 1 /(16 pi S) d m^2 d cth d y
      xjac=1d0/kn_sbeams
      zlow = atan((ph_Wmass2low  - ph_Wmass2)/ph_WmWw)
      zhigh= atan((min(10*ph_Wmass2high,kn_sbeams) - ph_Wmass2)/ph_WmWw)
      z =zlow+(zhigh-zlow)*xx(1)
      xjac=xjac*(zhigh-zlow)
      m2=ph_WmWw*tan(z)+ph_Wmass2
c d m^2 jacobian
      xjac=xjac*ph_WmWw/cos(z)**2
c d x1 d x2 = d tau d y;
      tau=m2/kn_sbeams
      kn_sborn=m2
      rs=sqrt(m2)
c ymax=|log(tau)|/2
      y=-(1-2*xx(2))*log(tau)/2
      xjac=-xjac*log(tau)

c initial state particles
      kn_cmpreal(0,1)= rs/2
      kn_cmpreal(1,1)= 0
      kn_cmpreal(2,1)= 0
      kn_cmpreal(3,1)= kn_cmpreal(0,1)

      kn_cmpreal(0,2)= kn_cmpreal(0,1)
      kn_cmpreal(1,2)= 0
      kn_cmpreal(2,2)= 0
      kn_cmpreal(3,2)=-kn_cmpreal(0,2)

c Build kinematics
      kn_x1=sqrt(tau)*exp(y)
      kn_x2=tau/kn_x1

      z=1-2*xx(3)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)*2d0*pi

      sth= sqrt( (1d0-cth)*(1d0+cth) )
      phi= 0d0
      cph= cos(phi)
      sph= sin(phi)

      if (mode.eq.1) then
*
* 2 body kinematics -->
*
* -> FS    
          v(0) = 1d0
          v(1) = sth*cph
          v(2) = sth*sph
          v(3) = cth
*
* <-- || 
*     || 3 bodies kinematics
*     || sceglie energia / momento fotone -->
* -> E
          eps = 1d-7*rs
          emax = ( m2 - mlep2 )/( 2d0*rs )
          pg(0) = emax**xx(4)*eps**(1d0-xx(4))
          xjac=xjac*log(emax/eps) * pg(0)
* -> p
          phiph = 2d0*pi*xx(5)
          xjac = xjac*2d0*pi
          cphg  = cos(phiph)
          sphg  = sin(phiph)
*
          betal = sqrt( (1d0-kn_masses(3)/rs*2)
     +                 *(1d0+kn_masses(3)/rs*2) )
          if (kn_masses(3)/rs*2d0.lt.1d-4) then

              ombetal =  + mlep2/(2.*(rs/2d0)**2) 
     +                   + mlep2**2/(8.*(rs/2d0)**4) 
     +                   + mlep2**3/(16.*(rs/2d0)**6) 
     +                   + (5*mlep2**4)/(128.*(rs/2d0)**8)

          else

              ombetal = 1d0 - betal

          endif

          anorm = log( (1d0+betal)/ombetal )/betal
          cthg = ( 1d0 - ( 1d0 + betal)
     +                   *((1d0+betal)/ombetal)**(-xx(6)) )
     +           /betal
          xjac = xjac*anorm * (1d0 - betal*cthg)
          sthg = sqrt( (1.d0-cthg)*(1d0+cthg) )
*
          pg(1) = pg(0)*sthg*cphg
          pg(2) = pg(0)*sthg*sphg
          pg(3) = pg(0)*cthg
          call rot(-1,v,pg,pg)
* <-- || 
          odv  = pg(1)*v(1) + pg(2)*v(2) + pg(3)*v(3)
          abig = m2 - 2.d0 * rs * pg(0) + mlep2 
          bbig = -2.d0 * (rs - pg(0))
          cbig = -2.d0 * odv
          arg =   bbig*bbig*( abig*abig
     +                       - mlep2*(bbig*bbig-cbig*cbig) )
*
          p1mod_1 = ( abig*cbig + sqrt(arg) )
     +                         / (bbig*bbig - cbig*cbig)
          p1mod_2 = ( abig*cbig - sqrt(arg) )
     +                         / (bbig*bbig - cbig*cbig)

          if (p1mod_1.gt.0.d0.and.p1mod_2.gt.0.d0) then
             if (random().lt.0.5) then
                p1mod = p1mod_1
             else 
                p1mod = p1mod_2
             endif
             xjac = xjac*2d0
          endif
          if (p1mod_1.gt.0.d0.and.p1mod_2.le.0.d0) then
             p1mod = p1mod_1
          endif
          if (p1mod_1.le.0.d0.and.p1mod_2.gt.0.d0) then
             p1mod = p1mod_2
          endif

          pl(0) = sqrt(p1mod*p1mod + mlep2)

          pn(0) = rs - pl(0) - pg(0)

          do mu=1,3
              pl(mu) =   p1mod * v(mu)
              pn(mu) = - pg(mu) - pl(mu)
          enddo

          den = abs( pn(0)*p1mod + pl(0)*(p1mod + odv) )
          xjac = xjac *pg(0) * p1mod*p1mod / den
     +           /(2d0*pi)**5 / 8d0

          do mu=0,3
              kn_cmpreal(mu,3) = pl(mu) 
              kn_cmpreal(mu,4) = pn(mu)
              kn_cmpreal(mu,5) = pg(mu)
          enddo

          if (abs(dotp(pn,pn)).gt.1d-10) then
              xjac = 0
          endif

      elseif (mode.eq.2) then

c decay products in their rest frame
          anorm = log( (m2-kn_masses(3)**2)/mlep2 )
          mpart2 = exp(xx(4)*anorm)*mlep2
          

          call dphi2(m2,mpart2,0d0,cth,phi,q1,q2,f123)

          phi12 = xx(5)*2d0*pi
          cth12 = xx(6)*2d0 - 1d0

          call dphi2(mpart2,mlep2,0d0,cth12,phi12,q1s,q2s,f12)

          xjac = xjac * f12 * f123
     +                  *2d0*(2*pi)
     +                  *mpart2 * anorm
     +                  *(2*pi)**4
     +                  /16d0/(2*pi)**9

          do mu=0,3
              kn_cmpreal(mu,4) = q2(mu)
          enddo

          call boost(q1s,q1,kn_cmpreal(0,3))
          call boost(q2s,q1,kn_cmpreal(0,5))

      endif

      beta=(kn_x1-kn_x2)/(kn_x1+kn_x2)
      call mboost(nlegreal-2,vec,beta,kn_cmpreal(0,3),kn_preal(0,3))
      do mu=0,3
         kn_preal(mu,1)=kn_x1*kn_beams(mu,1)
         kn_preal(mu,2)=kn_x2*kn_beams(mu,2)
      enddo

      kn_csi = kn_cmpreal(0,5)/sqrt(kn_sborn)*2d0
      kn_y   = kn_cmpreal(3,5)/kn_cmpreal(0,5)
      kn_azi = atan(kn_cmpreal(2,5)/kn_cmpreal(1,5))

      kn_xb1 = kn_x1
      kn_xb2 = kn_x2

      m2inv = dotp( kn_cmpreal(0,3)+kn_cmpreal(0,4),
     +              kn_cmpreal(0,3)+kn_cmpreal(0,4) )

      if ( m2inv.lt.ph_Wmass2low.or.m2inv.gt.ph_Wmass2high ) xjac = 0d0

      call compdij
      call setsoftvecisr
      call compdijsoft
      call compcsimax

      end
*
**
*
      subroutine dphi2(qi2,qj2,qk2,cthj,phij,qj,qk,flux)
*
      implicit none
      real*8 qj(0:3),qk(0:3)
      real*8 qi2,qj2,qk2,cthj,sthj,phij,flux
      real*8 cphij,sphij,qjm,qjm2
*
      sthj= sqrt( (1.d0-cthj)*(1d0+cthj) )
      cphij= cos(phij)
      sphij= sin(phij)
*
      qjm2= (qi2+qk2-qj2)*(qi2+qk2-qj2)/4.d0/qi2 - qk2

      qjm= sqrt(qjm2)
*
* momenta are calculated
*
      qj(0)= sqrt(qj2+qjm2)
      qj(1)= qjm*sthj*cphij
      qj(2)= qjm*sthj*sphij
      qj(3)= qjm*cthj
*
      qk(0)= sqrt(qk2+qj(1)*qj(1)+qj(2)*qj(2)+qj(3)*qj(3))
      qk(1)= -qj(1)
      qk(2)= -qj(2)
      qk(3)= -qj(3)
*
* flux is calculated
*
      flux= qjm/(qj(0)+qk(0))
*
      return
      end

*
* input:  qjs_0 ,qjs_x, qjs_y, qjs_z, qi_0, qi_1, qi_2, qi_3
* output: qj_0, qj_x, qj_y, qj_z, dj2
*
      subroutine boost(qjs,qi,qj)
*
      implicit real*8 (a-h,o-z)
      dimension qjs(0:3),qi(0:3),qj(0:3)
      dimension qjt(3),an(3),qjsl(3)
*
      qil= sqrt(qi(1)*qi(1)+qi(2)*qi(2)+qi(3)*qi(3))
*
      zero = 1.d-30
*
      if(qi(0).gt.zero.and.qil.gt.zero) then
         bi= qil/qi(0)
         if(bi.ge.1.d0) then
            gml = zero
         else
            gml= 1.d0/sqrt(1.d0-bi*bi)
         endif
      else
         bi = 0.d0
         gml = 1.d0
      endif
*
      if(abs(qil).lt.zero) then
         an(1)= 0.d0
         an(2)= 0.d0
         an(3)= 0.d0
      else
         an(1)= -qi(1)/qil
         an(2)= -qi(2)/qil
         an(3)= -qi(3)/qil
      endif
*
      qjslm= qjs(1)*an(1)+qjs(2)*an(2)+qjs(3)*an(3)
      qjsl(1)= qjslm*an(1)
      qjsl(2)= qjslm*an(2)
      qjsl(3)= qjslm*an(3)
      qjt(1)= qjs(1)-qjsl(1)
      qjt(2)= qjs(2)-qjsl(2)
      qjt(3)= qjs(3)-qjsl(3)
*
      qj(0)= gml*(qjs(0)-bi*qjslm)
      qjlm= gml*(-bi*qjs(0)+qjslm)
*
      qj(1)= qjlm*an(1)+qjt(1)
      qj(2)= qjlm*an(2)+qjt(2)
      qj(3)= qjlm*an(3)+qjt(3)
*
      return
      end
*
**
*
      subroutine rot(idir,vect,pin,pout)
! written by CMCC, last modified 9/10/2005
      implicit double precision (a-h,o-z)       
      double precision pin(0:3),pout(0:3),pp(0:3),r(3,3),
     >     vers(3),vect(0:3)
* This subroutine rotates the 4-vector pin in the frame where the z-axis is
* directed along the 4-vector vect(0,1,2,3). The rotated vector is stored
* in pout
* idir =  1 ---> direct rotation matrix
* idir = -1 ---> inverse rotation matrix
      pp(0) = pin(0)
      pp(1) = pin(1)
      pp(2) = pin(2)
      pp(3) = pin(3)

      vmo = 1.d0/sqrt(vect(1)**2+vect(2)**2+vect(3)**2)
      vers(1) = vect(1)*vmo
      vers(2) = vect(2)*vmo
      vers(3) = vect(3)*vmo
      vt = sqrt(vers(1)**2+vers(2)**2)

!   BUG - pointed out by CLEO people
!      v1ovt = vers(1)/vt
!      if (vt.eq.0.d0) v1ovt = 0.d0
!      v2ovt = vers(2)/vt
!      if (vt.eq.0.d0) v2ovt = 1.d0
 
      v1ovt = 0.d0
      v2ovt = 1.d0
      if (vt.gt.0.d0) then
         v1ovt = vers(1)/vt
         v2ovt = vers(2)/vt
      endif
      
      if (idir.eq.(-1)) then    !! INVERSE rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(1,2) = -v2ovt
         r(1,3) =  vers(1)      
         r(2,1) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(2,3) =  vers(2)
         r(3,1) = -vt
         r(3,2) =  0.d0
         r(3,3) =  vers(3)
      else  ! if (idir.eq.1) !! DIRECT rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(2,1) = -v2ovt
         r(3,1) =  vers(1)
         r(1,2) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(3,2) =  vers(2)
         r(1,3) = -vt
         r(2,3) =  0.d0
         r(3,3) =  vers(3)
      endif
      pout(0) = pp(0)
      pout(1) = r(1,1)*pp(1) + r(1,2)*pp(2) + r(1,3)*pp(3)
      pout(2) = r(2,1)*pp(1) + r(2,2)*pp(2) + r(2,3)*pp(3)
      pout(3) = r(3,1)*pp(1) + r(3,2)*pp(2) + r(3,3)*pp(3)
      return
      end
*
