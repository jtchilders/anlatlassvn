      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,tau,y,beta,vec(3),s,
     #     z,zhigh,zlow
      integer mu,k
      logical ini
      data ini/.true./
      save ini
      real * 8 Hmass2,Hmass2low,Hmass2high,HmHw  
      real * 8 mass,bwfactor,powheginput      
      external mass,powheginput      
      logical check
      parameter(check=.false.)
      logical BW
      integer BWshape
      save BW,BWshape
      integer BWflag
      real * 8 weight,m
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         if (powheginput("#zerowidth").eq.1) then
            BW=.false.
         else
            BW=.true.
         endif
         if (BW) then 
            BWshape=powheginput("#bwshape") 
            if (BWshape.lt.0) BWshape=1
            if (BWshape.gt.3) then 
               write(*,*) 'Unsupported value for bwshape in the '//
     $              'powheg.input file: ',powheginput("#bwshape")
               write(*,*) 'The POWHEG BOX aborts'
               call exit(1)
            endif
         endif
         ini=.false.
      endif

      Hmass2 = ph_Hmass2
      Hmass2low = ph_Hmass2low
      Hmass2high = ph_Hmass2high
      HmHw = ph_HmHw

c     Factor due to physical phase space
c     2*pi/S dM^2 dY delta(s- M^2)
      xjac=2d0*pi/kn_sbeams

      if (BW) then
         zlow=atan((Hmass2low  - Hmass2)/HmHw)
         zhigh=atan((Hmass2high - Hmass2)/HmHw)
         z=zlow+(zhigh-zlow)*xborn(1)
         xjac=xjac*(zhigh-zlow)
         m2=HmHw*tan(z)+Hmass2
c     d m^2 jacobian
         xjac=xjac*HmHw/cos(z)**2
      else
         m2 = Hmass2
         Hmass2low=m2
         Hmass2high = m2
      endif
c     Assign the Higgs boson mass
      kn_masses(3)=ph_Hmass

      if(BW) then
C     Evaluate BW factor due to the presence of an 
c     off-shell Higgs boson (see eq. 2.6 of 0812.0578)
c     In the case of Higgs running width,
c     the cross section is obtained with the formula
c     (32) of  hep-ph/9504378 
         if (bwshape.eq.1) then
            bwfactor=m2*HmHw/Hmass2/pi/
     $           ((m2- Hmass2)**2+(m2*HmHw/Hmass2)**2)
C     MC@NLO and MCFM use instead:
c     BW(M,M0,Ga)=M0 Ga/pi * 1/((M^2-M0^2)^2+M0^2 Ga^2
         elseif (bwshape.eq.2) then   
            bwfactor=HmHw/pi/
     $        ((m2-Hmass2)**2+HmHw**2)    
         elseif (bwshape.eq.3) then   
c     Passarino complex-pole scheme
c     first: attach the fixed-width BW. the 1/pi factor is attached together to 
c     the weight
            bwfactor=HmHw/
     $        ((m2-Hmass2)**2+HmHw**2)          
c     second: compute the rescale factor
            m=sqrt(m2)
            BWflag=0
            call pwhg_cpHTO_reweight(ph_Hmass,ph_Hwidth,ph_topmass,
     $           BWflag,m,weight)
            xjac = xjac * weight/pi
         else 
            write(*,*) 'ERROR in bwshape'
            call exit(1)
         endif
c     BW factor 
         xjac=xjac*bwfactor
      endif

c     d x1 d x2 = d tau d y;
      tau=m2/kn_sbeams
      s=kn_sbeams*tau
      kn_sborn=s
c ymax=|log(tau)|/2
      y=-(1-2*xborn(2))*log(tau)/2
      xjac=-xjac*log(tau)
      kn_born_pt2=0d0

      kn_jacborn = xjac

c     Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1

c initial state particles
      kn_cmpborn(0,1)=sqrt(s)/2
      kn_cmpborn(0,2)=kn_cmpborn(0,1)
      kn_cmpborn(3,1)=kn_cmpborn(0,1)
      kn_cmpborn(3,2)=-kn_cmpborn(0,2)
      kn_cmpborn(1,1)=0
      kn_cmpborn(1,2)=0
      kn_cmpborn(2,1)=0
      kn_cmpborn(2,2)=0      
      do mu=0,3
         kn_cmpborn(mu,3)=kn_cmpborn(mu,1)+kn_cmpborn(mu,2)
      enddo
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
      if(check) then
         call checkmomzero(nlegborn,kn_pborn)
         call checkmass(2,kn_pborn(0,3))
      endif
c minimal value of sqrt(s)
      kn_minmass=sqrt(Hmass2low)
      end


      function mass(p)
      implicit none
      real * 8 p(0:3),mass
      mass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
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
               write(*,*) '    scales set to the H virtuality '            
            else 
               write(*,*) "runningscale value not allowed"
               call exit(1)
            endif
            write(*,*) '*************************************'
            ini=.false.
         endif
         muref=sqrt(dotp(kn_pborn(0,3),kn_pborn(0,3)))
      else
         if (ini) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to the H mass '
            write(*,*) '*************************************'
            ini=.false.
         endif
         muref=ph_Hmass
      endif
      muf=muref
      mur=muref
      end
