c     fills the four-momenta of the Born-process particles,
c     both in the laboratory and in the center-of-mass frame,
c     the Bjorken x of the two incoming partons,
c     the partonic Mandelstam variable s,
c     the masses for Born and real emission,
c     and the Jacobi factor for the phase space volume
      subroutine born_phsp(xborn)
      implicit none
#include "PhysPars.h"
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      real *8 xborn(ndiminteg-3),vec(3),beta
      integer mu

      kn_jacborn=1

      ! set masses for Born and real emission (-> phasespace)
      kn_masses(1)=0
      kn_masses(2)=0
      kn_masses(3)=ph_slepton3mass !slepton_mass
      kn_masses(4)=ph_slepton4mass !slepton_mass
      kn_masses(5)=0
      kn_minmass =  kn_masses(3)+kn_masses(4)

      call x1x2phspace(kn_sbeams,kn_minmass,xborn(1),kn_xb1,kn_xb2,
     &     kn_sborn,kn_jacborn)
      kn_cmpborn(0,1)=sqrt(kn_sborn)/2
      kn_cmpborn(3,1)=kn_cmpborn(0,1)
      kn_cmpborn(0,2)=kn_cmpborn(0,1)
      kn_cmpborn(3,2)=-kn_cmpborn(0,1)

      call twobody0(xborn(3),kn_sborn,kn_masses(3),kn_masses(4),
     &     kn_cmpborn(0,3),kn_cmpborn(0,4),kn_jacborn)

      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn-2,vec,beta,kn_cmpborn(0,3),kn_pborn(0,3))
      do mu=0,3
         kn_pborn(mu,1)=kn_xb1*kn_beams(mu,1)
         kn_pborn(mu,2)=kn_xb2*kn_beams(mu,2)
      enddo

      end



c     uses the two random numbers provided in xx
c     to generate the fraction of the beam momentum
c     for the two partons entering the Born process,
c     calculates the partonic Mandelstam variable s
c     and the Jacobi factor for the phase space volume
      subroutine x1x2phspace(sbeams,minmass,xx,x1,x2,s,jac)
      implicit none
      real *8 sbeams,minmass,xx(2),x1,x2,s,jac
      real *8 taumin,tau,y

      taumin=minmass**2/sbeams
      tau=exp(log(taumin)*(1-xx(1)**2))
      jac=jac*tau*abs(log(taumin))*2*xx(1)
      s=sbeams*tau
      y=-(1-2*xx(2))*log(tau)/2
      jac=-jac*log(tau)
      x1=sqrt(tau)*exp(y)
      x2=tau/x1

      end



c     uses the random number provided in xx
c     to generate four-momenta in the cms frame
c     for the particles produced in the Born process,
c     calculates the Jacobi factor for the phase space volume
      subroutine twobody0(xx,s,m1,m2,p1,p2,jac)
      implicit none
      include "../include/pwhg_math.h"
      real *8 xx(1),s,m1,m2,p1(0:3),p2(0:3),jac
      real *8 k,m,ctheta,stheta

      m=sqrt(s)
      if(m.le.m1+m2) then
         jac=0
         return
      endif
      k=sqrt((m**2-(m1+m2)**2)*(m**2-(m1-m2)**2))/(2*m)
      ctheta=1-2*xx(1)
      stheta=sqrt(abs(1-ctheta**2))
      p1(3)=ctheta*k
      p1(1)=stheta*k
      p1(2)=0
      p2(1)=-p1(1)
      p2(2)=-p1(2)
      p2(3)=-p1(3)
      p1(0)=sqrt(m1**2+k**2)
      p2(0)=sqrt(m2**2+k**2)
      jac=jac/(8*pi)*k/p1(0)

      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real *8 fact,pt
      real *8 powheginput
      external powheginput

      if (ini) then
         pt = powheginput("#ptsupp")
         if(pt.gt.0) then
            write(*,*) ' ******** WARNING: ptsupp is deprecated'
            write(*,*) ' ******** Replace it with bornsuppfact'
         else
            pt = powheginput("#bornsuppfact")
         endif
         if(pt.ge.0) then
            write(*,*) '**************************'
            write(*,*) 'No Born suppression factor'
            write(*,*) '**************************'
         endif
         ini=.false.
      endif
      fact=1d0

      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
#include "PhysPars.h"
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      real *8 muf,mur
      real *8 muref
      real *8 kn34(0:3)
      real *8 dotp
      external dotp
      logical ini
      data ini/.true./
      logical runningscales
      save runningscales
      real *8 powheginput
      external powheginput

      if(ini) then
         if(powheginput('#runningscale').ge.1) then
            if (powheginput('#runningscale').ne.1) then
               write(*,*) 'you set runningscale to ', runningscales
               write(*,*) 'this runningscale value is not allowed'
               call exit(1)
            endif
            runningscales=.true.
         else
            runningscales=.false.
         endif
      endif

      if (runningscales) then
         kn34(:) = kn_pborn(:,3)+kn_pborn(:,4)
         muref = dsqrt(abs(kn34(0)**2-kn34(1)**2
     &                    -kn34(2)**2-kn34(3)**2))
      else
         muref=ph_slepton3mass+ph_slepton4mass
      endif

      if (ini) then
        write(*,*) '***********************************************'
        write(*,*) '  Factorization and renormalization scales:'
        if (runningscales) then
          write(*,*) '  using invariant mass of sleptons for muref'
        else
          write(*,*) '  using sum of slepton masses for muref'
          write(*,*) '  muref = ', muref
        endif
        ! NOTE: scales will be multiplied (outside of this function) by
        ! an additional factor if renscfact or facscfact is set > 0
        ! in powheg.input
        write(*,*) '  muren = muref * ', st_renfact
        write(*,*) '  mufac = muref * ', st_facfact
        write(*,*) '***********************************************'
        ini = .false.
      endif

      muref = max(muref,dsqrt(2d0))  ! avoid too small values of the scales

      muf=muref
      mur=muref

      end
