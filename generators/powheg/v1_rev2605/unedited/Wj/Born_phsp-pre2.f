      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,taumin,tau,y,beta,vec(3),cth,cthdec,phidec,s,
     # z,zhigh,zlow,dir(3),kzed,cthmax
      integer mu,k,j
      logical ini
      data ini/.true./
      save ini
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         ini=.false.
      endif
c Phase space:
c d m2/(2 pi) d omega /(8*(2 pi)^2)  (s-m2)/2
c d omegadec/(8*(2 pi)^2)
c omega: 3d angle in CM system
c omegadec: 3d angle in CM system of W decay products
      zlow=atan((ph_Wmass2low  - ph_Wmass2)/ph_WmWw)
      zhigh=atan((min(ph_Wmass2high,kn_sbeams)  - ph_Wmass2)/ph_WmWw)
      z=zlow+(zhigh-zlow)*xborn(1)
      xjac=zhigh-zlow
      m2=ph_WmWw*tan(z)+ph_Wmass2
c d m^2/(2pi) jacobian
      xjac=xjac*ph_WmWw/cos(z)**2/(2*pi)
c d x1 d x2 = d tau d y;
      taumin=( sqrt(m2+kn_ktmin**2) + kn_ktmin )**2/kn_sbeams
      tau=exp(log(taumin)*(1-xborn(2)**2))
      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
      s=kn_sbeams*tau
      kn_sborn=s
c compute W momentum in partonic cm
      kzed=(s-m2)/(2*sqrt(s))
c ymax=|log(tau)|/2
      y=-(1-2*xborn(3))*log(tau)/2
      xjac=-xjac*log(tau)
c abs to protect from tiny negative values
      cthmax=sqrt(abs(1d0-(kn_ktmin/kzed)**2))
      z=1-2*xborn(4)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)
      cth=cth*cthmax
      kn_born_pt2=(1-cth**2)*kzed**2
      xjac=xjac*cthmax
c supply 2 pi for azimuthal integration (not performed)
      xjac=xjac*2*pi
      xjac=xjac*(s-m2)/s/(8*(2*pi)**2)
c
      cthdec=1-2*xborn(5)
      kn_cthdec=cthdec
      xjac=xjac*2
      phidec=2*pi*xborn(6)
      xjac=xjac*2*pi
      xjac=xjac/(8*(2*pi)**2)
      kn_jacborn=xjac
c Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1
c decay products in their rest frame
      kn_cmpborn(0,3)=sqrt(m2)/2
      kn_cmpborn(0,4)=kn_cmpborn(0,3)
      kn_cmpborn(3,3)=cthdec*kn_cmpborn(0,3)
      kn_cmpborn(1,3)=sqrt(1-cthdec**2)*sin(phidec)*kn_cmpborn(0,3)
      kn_cmpborn(2,3)=sqrt(1-cthdec**2)*cos(phidec)*kn_cmpborn(0,3)
      kn_cmpborn(1,4)=-kn_cmpborn(1,3)
      kn_cmpborn(2,4)=-kn_cmpborn(2,3)
      kn_cmpborn(3,4)=-kn_cmpborn(3,3)
c velocity of W in partonic CM
      beta=(s-m2)/(s+m2)
      vec(1)=sqrt(1-cth**2)
      vec(2)=0
      vec(3)=cth
      call mboost(2,vec,beta,kn_cmpborn(0,3),kn_cmpborn(0,3))
      kn_cmpborn(1,5)=-kn_cmpborn(1,3)-kn_cmpborn(1,4)
      kn_cmpborn(2,5)=-kn_cmpborn(2,3)-kn_cmpborn(2,4)
      kn_cmpborn(3,5)=-kn_cmpborn(3,3)-kn_cmpborn(3,4)
      kn_cmpborn(0,5)=sqrt(s)-kn_cmpborn(0,3)-kn_cmpborn(0,4)
      kn_cmpborn(0,1)=sqrt(s)/2
      kn_cmpborn(0,2)=kn_cmpborn(0,1)
      kn_cmpborn(3,1)=kn_cmpborn(0,1)
      kn_cmpborn(3,2)=-kn_cmpborn(0,2)
      kn_cmpborn(1,1)=0
      kn_cmpborn(1,2)=0
      kn_cmpborn(2,1)=0
      kn_cmpborn(2,2)=0      
c now boost everything along 3
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn-2,vec,beta,kn_cmpborn(0,3),kn_pborn(0,3))
      do mu=0,3
         kn_pborn(mu,1)=kn_xb1*kn_beams(mu,1)
         kn_pborn(mu,2)=kn_xb2*kn_beams(mu,2)
      enddo
c      call checkmomzero(nlegborn,kn_pborn)
c      call checkmass(2,kn_pborn(0,3))

c minimal final state mass 
      kn_minmass=sqrt(ph_Wmass2low)

      end


      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real * 8 fact,pt2,pt2supp,powheginput,pt
      save ini,pt2supp,pt     
c CAVEAT!!!  process dependent subroutine
      if (ini) then
         pt = powheginput("#ptsupp")
         if(pt.gt.0) then
            write(*,*) ' ******** WARNING: ptsupp is deprecated'
            write(*,*) ' ******** Replace it with bornsuppfact'
            call flush
            call exit(-1)
         else
            pt = powheginput("#bornsuppfact")
         endif
         ini = .false.
         pt2supp = pt**2
      endif
      if (pt.lt.0) then
         fact=1d0
      else         
         pt2=kn_pborn(1,5)**2+kn_pborn(2,5)**2
         fact=pt2/(pt2+pt2supp)         
      endif
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 muf,mur,powheginput
      logical ini
      data ini/.true./
      logical runningscales
      real * 8 pt2
      save ini,runningscales
      if(ini) then
         if(powheginput("#runningscales").eq.1) then
            runningscales = .true.
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=pt  used for Bbar function   **'
            write(*,*) '**   muf=pt  used for Bbar function   **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         else
c fixed scale is the default
            runningscales = .false.
            write(*,*) '********************************************'
            write(*,*) '********************************************'
            write(*,*) '**   mur=W mass  used for Bbar function   **'
            write(*,*) '**   muf=W mass  used for Bbar function   **'
            write(*,*) '********************************************'
            write(*,*) '********************************************'
         endif
         ini=.false.
      endif
      if (runningscales) then
         pt2=kn_pborn(1,5)**2+kn_pborn(2,5)**2
         mur=sqrt(pt2)
         muf=mur
      else
         muf=ph_Wmass
         mur=ph_Wmass
      endif
      end
