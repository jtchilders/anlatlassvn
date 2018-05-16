      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,tau,y,beta,vec(3),cth,cthdec,phidec,s,
     # z,zhigh,zlow
      integer mu,k
      logical ini
      data ini/.true./
      save ini
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      real*8 rkallen
      external rkallen
      real*8 pmod,sqla,betal
      real * 8 masswindow_low,masswindow_high
      real * 8 phsp_Wm,phsp_Ww,phsp_Wmass2,phsp_WmWw
      save phsp_Wm,phsp_Ww,phsp_Wmass2,phsp_WmWw
      real * 8 powheginput
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegreal
            kn_masses(k)=0
         enddo
         kn_masses(3)=decmass
         ini=.false.
         phsp_Wm=powheginput('#phsp_Wm')
         phsp_Ww=powheginput('#phsp_Ww')
         if(phsp_Wm.lt.0) phsp_Wm=ph_Wmass
         if(phsp_Ww.lt.0) phsp_Ww=ph_Wwidth
         phsp_Wmass2=phsp_Wm**2
         phsp_WmWw=phsp_Wm*phsp_Ww

c     mass window
         masswindow_low = powheginput("#masswindow_low")
         if (masswindow_low.le.0d0) masswindow_low=30d0
         masswindow_high = powheginput("#masswindow_high")
         if (masswindow_high.le.0d0) masswindow_high=30d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     set mass window around W-mass peak in unit of ph_Wwidth
c     It is used in the generation of the Born phase space
         ph_Wmass2low=max(decmass*10d0,phsp_Wm-masswindow_low*phsp_Ww)
         ph_Wmass2low=ph_Wmass2low**2
         ph_Wmass2high=phsp_Wm+masswindow_high*phsp_Ww
         ph_Wmass2high=min(kn_sbeams,ph_Wmass2high**2)
      endif

c 1 /(16 pi S) d m^2 d cth d y
      xjac=1d0/kn_sbeams/(16*pi)
      zlow=atan((ph_Wmass2low - phsp_Wmass2)/phsp_WmWw)
      zhigh=atan((ph_Wmass2high - phsp_Wmass2)/phsp_WmWw)
      z=zlow+(zhigh-zlow)*xborn(1)
      m2=phsp_WmWw*tan(z)+phsp_Wmass2
c d m^2 jacobian
      xjac=xjac*(zhigh-zlow)*(phsp_WmWw/cos(z)**2)
c d x1 d x2 = d tau d y;
      tau=m2/kn_sbeams
      s=kn_sbeams*tau
      kn_sborn=s
c ymax=|log(tau)|/2
      y=-(1-2*xborn(2))*log(tau)/2
      xjac=-xjac*log(tau)
      z=1-2*xborn(3)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)
      kn_born_pt2=0d0
      
c
      sqla = rkallen(m2,kn_masses(3)**2,0d0)
      xjac = xjac * sqla/m2

      cthdec=cth
      phidec=0d0
      kn_cthdec=cthdec
      kn_jacborn=xjac
c Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1
c decay products in their rest frame
      pmod = rkallen(m2,kn_masses(3)**2,0d0)/sqrt(m2) / 2d0

      kn_cmpborn(0,3)=sqrt( pmod**2 + kn_masses(3)**2 )
      kn_cmpborn(1,3)=pmod * sqrt(1-kn_cthdec**2)*cos(phidec)
      kn_cmpborn(2,3)=pmod * sqrt(1-kn_cthdec**2)*sin(phidec)
      kn_cmpborn(3,3)=pmod * kn_cthdec

      kn_cmpborn(0,4)= sqrt(m2) - kn_cmpborn(0,3)
      kn_cmpborn(1,4)= -kn_cmpborn(1,3)
      kn_cmpborn(2,4)= -kn_cmpborn(2,3)
      kn_cmpborn(3,4)= -kn_cmpborn(3,3)
c initial state particles
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
      real * 8 fact,pt
      real * 8 powheginput
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
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'mathx.h'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      real *8 muref
      real *8 dotp
      external dotp
      logical runningscales
      save runningscales
      real * 8 pt2
      real * 8 powheginput
      external powheginput
      if(ini) then
         if(powheginput('#runningscale').ge.1) then
            runningscales=.true.
         else
            runningscales=.false.
         endif
      endif
      if (runningscales) then
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            if (powheginput('#runningscale').eq.1) then
               write(*,*) '    scales set to the W virtuality ' 
            elseif (powheginput('#runningscale').eq.2) then
               write(*,*) '    scales set to the W transverse mass ' 
            else 
               write(*,*) "runningscale value not allowed"
               call exit(1)
            endif
            write(*,*) '*************************************'
            ini=.false.
         endif
         if(powheginput('#runningscale').eq.2) then
            pt2=(kn_pborn(1,3)+kn_pborn(1,4))**2+(kn_pborn(2,3)
     $           +kn_pborn(2,4))**2
            muref=sqrt(pt2+ph_Wmass*ph_Wmass)
         else
            muref=sqrt(2d0*dotp(kn_pborn(0,3),kn_pborn(0,4)))
         endif
      else
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to the W mass '
            write(*,*) '*************************************'
            ini=.false.
         endif
         muref=ph_Wmass
      endif
      muf=muref
      mur=muref
      mudim2=1d0
      end

