      subroutine born_phsp2(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'

      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow
      integer mu,k,j
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw
      real * 8 m2jj,pv(4),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2)
      real * 8 mass      
      external mass      
      logical check
      parameter(check=.false.)
      logical BW
      parameter (BW=.true.)
      real * 8 epsilon
      parameter (epsilon=1d-10)

c generation cuts:
      real*8 mllmin
      real * 8 pt1cut,pt2cut,pt1,pt2
      real * 8 m2ll,m2llmin,m2min,q21,q43
      logical generation_cuts

      real * 8 BW_fixed
      real * 8 powheginput
      external powheginput
c
c extra stuff for Zjj:
      real*8 z1
      real*8 wbw,vmass,vwidth
      real*8 phiv,costhv,sinthv,cphiv,sphiv
      real*8 phil,costhl,sinthl,cphil,sphil
      real*8 plmod,pl1_cm(4),pl(0:3,2),pl1(4),pl2(4)

c------------------------------------
c
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         ini=.false.

         select case (procid)
         case (Zjj_l)
           if(powheginput('#mll_gencut').lt.15d0) then
              m2llmin=225d0
           else
              m2llmin=powheginput('#mll_gencut')   
              m2llmin=m2llmin**2
           endif
           Vmass2 = ph_Zmass2
           VmVw = ph_ZmZw
           Vmass = ph_Zmass
           Vwidth = ph_Zwidth         
         
         case default
           m2llmin=0d0        
           Vmass2 = ph_wmass2
           VmVw = ph_wmww
           Vmass = ph_wmass
           Vwidth = ph_wwidth         
         end select
         
        
c code allows for additional generation cuts 
c     (to be used instead of Born_suppression factor)
c     (disabled by default):
         pt1cut = powheginput('#ptj_gencut') 
         pt2cut = powheginput('#ptj_gencut')
         if(pt1cut.lt.0d0) then
            pt1cut=0d0
            pt2cut=0d0
         endif

         m2min=0.1d0**2d0
         if ((pt1cut.ne.0d0).or.(pt2cut.ne.0d0)
     &       .or.(m2llmin.ne.0d0)) then
            generation_cuts = .true.
            write(*,*) '*************************************'
            write(*,*) '****    PS-CUTS IN PLACE!!!      ****' 
            write(*,*) '*************************************'
            write(*,*) ''
            write(*,*) 'jet cuts:'
            write(*,*) 'ptj1_min = ',pt1cut
            write(*,*) 'ptj2_min = ',pt2cut
            write(*,*) ''
            write(*,*) 'lepton cuts:'
            write(*,*) 'mll_min    = ',dsqrt(m2llmin)
            write(*,*) ''
         endif

      endif

      Vmass2low = max(m2llmin,m2min)
      Vmass2high = kn_sbeams/2d0 
      
      ! init:
      xjac = 1d0


      if (BW) then
c variable trafo:
         z1 = xborn(1)**6
         xjac = xjac*6d0*xborn(1)**5
         call breitw(z1,Vmass2low,Vmass2high,Vmass,Vwidth,m2,wbw)
         xjac = xjac*wbw
      else
         xjac = 1d0
         m2 = Vmass2
      endif


c     d x1 d x2 = d tau d y;
      taumin=m2/kn_sbeams
      tau=exp(log(taumin)*(1-xborn(2)**2))
      xjac=xjac*tau*abs(log(taumin))*2*xborn(2)
      s=kn_sbeams*tau
      kn_sborn=s
c     ymax=|log(tau)|/2
      y=-(1-2*xborn(3))*log(tau)/2
      xjac=-xjac*log(tau)

c     generate dijet squared mass m2jj
      m2jj = xborn(4)*(sqrt(s)-sqrt(m2))**2
      xjac = xjac * (sqrt(s)-sqrt(m2))**2
      
      pV(4)=(s+m2-m2jj)/(2*sqrt(s))
      pVmod2 = pV(4)**2-m2
      
      if (pVmod2.lt.epsilon) then
         pVmod = 0d0
c         write(*,*) '============================>',pVmod2
      else
         pVmod = sqrt(pVmod2)
      endif

      costhv=2d0*xborn(8)-1d0      
      phiv=2d0*pi*xborn(9)
      sinthv=dsqrt(1d0-costhv**2)
      cphiv=dcos(phiv)
      sphiv=dsin(phiv)
      xjac = xjac*4d0*pi

      pV(1) = pVmod*sinthv*sphiv
      pV(2) = pVmod*sinthv*cphiv
      pV(3) = pVmod*costhv

      if (check) then
         print*,'pV_mass=',dsqrt(pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2)
      endif    

c lepton1:
      costhl=2d0*xborn(5)-1d0      
      phil=2d0*pi*xborn(10)
      sinthl=dsqrt(1d0-costhl**2)
      cphil=dcos(phil)
      sphil=dsin(phil)
      xjac = xjac*4d0*pi

      plmod = dsqrt(m2)/2d0
      pl1_cm(1) = plmod*sinthl*sphil
      pl1_cm(2) = plmod*sinthl*cphil
      pl1_cm(3) = plmod*costhl      
      pl1_cm(4) = plmod

c boost into lab frame:
      call boost(dsqrt(m2),pV,pl1_cm,pl1)
      pl2(:) = pV(:)-pl1(:)

      pl(1:3,1) = pl1(1:3)
      pl(1:3,2) = pl2(1:3)
      pl(0,1) = pl1(4)
      pl(0,2) = pl2(4)

c     supply the other factors to the jacobian
c     factor for the two-body jet phase space
      xjac=xjac/(8*(2*pi)**2)
c
C     anothor factor for 1->2 body phase space 
      xjac=xjac/(8*(2*pi)**2)
c
cc     factor for V production
      xjac=xjac/(4*(2*pi)**4)*pVmod/sqrt(s)

c     Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1
c     boost back in the lab frame
c     now boost everything along 3rd axis
      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1

      call mboost(2,vec,betaCM,pl,kn_pborn(0,3))      

c     build jet momenta in the jet CM frame
      pJ(0,1) = sqrt(m2jj)/2
      pJ(0,2) = pJ(0,1)
c     azimuth and polar angle of a jet
      z=1-2*xborn(6)**2
      xjac=xjac*4*xborn(6)
      cthj=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)

      phij = 2*pi*xborn(7)
      xjac=xjac*2*pi

      kn_jacborn = xjac

      pJ(1,1) = pJ(0,1)*sqrt(1-cthj**2)*sin(phij)
      pJ(2,1) = pJ(0,1)*sqrt(1-cthj**2)*cos(phij)
      pJ(3,1) = pJ(0,1)*cthj
      do mu=1,3
         pJ(mu,2)=-pJ(mu,1)
      enddo

      do mu=0,3
         kn_pborn(mu,1)=kn_xb1*kn_beams(mu,1)
         kn_pborn(mu,2)=kn_xb2*kn_beams(mu,2)
      enddo

      if (check) then
         call mboost(2,vec,-betaCM,kn_pborn(0,1),ptmp(0,1))      
         write(*,*) 'CM vec1',(ptmp(mu,1),mu=0,3)
         write(*,*) 'CM vec2',(ptmp(mu,2),mu=0,3)
      endif

c     boost in the lab frame
c     compute first p_plus+p_minus-pV
      do mu=0,3
         pcmjj(mu)= kn_pborn(mu,1) + kn_pborn(mu,2) 
     &            - kn_pborn(mu,3) - kn_pborn(mu,4)
      enddo
      pcmjjmod = sqrt(pcmjj(1)**2+pcmjj(2)**2+pcmjj(3)**2)
c     recompute pcmjj(0) from m2jj, otherwise there are points where 
c     beta > 1 or beta < 0
      pcmjj(0) = sqrt(m2jj+pcmjjmod**2)

c      write(*,*) '1 ===========> ',(pcmjj(0)**2-pcmjjmod**2)/m2jj

      beta=pcmjjmod/pcmjj(0)

      do mu=1,3
         vec(mu)=pcmjj(mu)/pcmjjmod
      enddo

      call mboost(2,vec,beta,pJ(0,1),kn_pborn(0,5))


      if (generation_cuts) then
c     jet cuts:
      pt1 = sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)
      pt2 = sqrt(kn_pborn(1,6)**2+kn_pborn(2,6)**2) 
      if ((pt1.lt.pt1cut).or.(pt2.lt.pt2cut)) then      
         kn_jacborn=0d0
      endif
c
c     lepton cuts:
      m2ll = (kn_pborn(0,3)+kn_pborn(0,4))**2-
     &       (kn_pborn(1,3)+kn_pborn(1,4))**2-
     &       (kn_pborn(2,3)+kn_pborn(2,4))**2-
     &       (kn_pborn(3,3)+kn_pborn(3,4))**2
      if (m2ll.lt.m2llmin) then  
         kn_jacborn=0d0
      endif


      endif !generation_cuts

      q21=   (kn_pborn(0,1)-kn_pborn(0,5))**2-

     &       (kn_pborn(1,1)-kn_pborn(1,5))**2-

     &       (kn_pborn(2,1)-kn_pborn(2,5))**2-

     &       (kn_pborn(3,1)-kn_pborn(3,5))**2



      q43=   (kn_pborn(0,2)-kn_pborn(0,6))**2-

     &       (kn_pborn(1,2)-kn_pborn(1,6))**2-

     &       (kn_pborn(2,2)-kn_pborn(2,6))**2-

     &       (kn_pborn(3,2)-kn_pborn(3,6))**2

      if((abs(q21).lt.m2min).or.(abs(q43).lt.m2min)) then

         kn_jacborn=0d0

      endif      
      

      if (check) then
         call mboost(1,vec,-beta,pcmjj(0),ptmp(0,1))
         write(*,*) 'only time component ==> ',(ptmp(mu,1),mu=0,3)
         write(*,*) ''
         write(*,*) 'new set'
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_pborn(mu,j),mu=0,3)
            write(*,*) 'mass ',j,mass(kn_pborn(0,j))
         enddo
         call checkmomzero(nlegborn,kn_pborn)
      endif
      
c     boost in the CM frame
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,-betaCM,kn_pborn(0,1),kn_cmpborn(0,1))      
      if (check) then
         write(*,*) ''
         write(*,*) 'new set'     
         do j=1,nlegborn
            write(*,*) 'mom ',j,(kn_cmpborn(mu,j),mu=0,3)
         enddo
      endif

      select case (procID)
      case(Zjj_l)
         kn_minmass=sqrt(m2llmin)
      case default
         kn_minmass=sqrt(m2min)
      end select

      end

c================

      function mass(p)
      implicit none
      real * 8 p(0:3),mass
      mass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end

c================
    
      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'process.inc'
      include 'vbfnlo-files/global.inc'      
      real * 8 muf,mur
      logical ini
      data ini/.true./
      real *8 muref
      real *8 muref_sq
      real*8 q12(0:4),q34(0:4)
      integer mu
      real *8 dotp,powheginput
      external dotp,powheginput
      real * 8 fixedscale,renscfact,facscfact
      logical runningscales
      save runningscales
      save renscfact,facscfact
      logical htotscale,qscale
      save htotscale,qscale
      real * 8 ptj1,ptj2,ptj3
      real * 8 htot,etp,etm

         if (ini) then
            select case (procID) 
            case (Zjj_l)
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to the Z mass '
            write(*,*) '*************************************'
            ini=.false.
            case default
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to the W mass '
            write(*,*) '*************************************'
            ini=.false.   
            end select
         endif
         select case (procID)
         case (Zjj_l)
           muref=ph_Zmass
         case default
           muref=ph_wmass
         end select

        muf=muref
        mur=muref

      end

cccccc

      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'

      real * 8 fact,ptmin
      real * 8 pt6,pt5
      real * 8 kp
      real * 8 p56(0:3)
      logical, save :: ini = .true.  
      real*8 powheginput
      external powheginput

      ptmin=powheginput('#ptsuppfact') 
      kp = 2d0

      if (ini) then 
         if (flg_weightedev) then 
            write(*,*) 'Using Born suppression' 
            write(*,*) 'with ptmin[GeV]=',ptmin 
            write(*,*) 'exponents: kp = ',kp 
         else
            write(*,*) 'Using no Born suppression' 
         endif
         ini = .false. 
      endif

      if(flg_weightedev) then
         pt6=kn_cmpborn(1,6)**2+kn_cmpborn(2,6)**2
         pt5=kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2
         fact=(pt5/(pt5+ptmin**2))**kp*(pt6/(pt6+ptmin**2))**kp
      else
         fact=1
      endif
      end

