      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      real * 8 xborn(ndiminteg-3)
      call born_phsp_2(xborn)
      end

      subroutine born_suppression(fact)
      implicit none
      double precision fact
      call born_suppression_2(fact)
      end


      subroutine born_phsp_2(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,taumin,tau,y,beta,betaCM,vec(3),cth,s,
     #     z,zhigh,zlow
      integer mu,k,j,ileg
      logical ini
      data ini/.true./
      save ini
      real * 8 Vmass2,Vmass2low,Vmass2high,VmVw
      real * 8 m2jj,pV(0:3),pVmod,pVmod2,pJ(0:3,2),cthj,phij,pcmjj(0:3),
     #     pcmjjmod,ptmp(0:3,2)
      real * 8 mass
      external mass
      logical check
      parameter(check=.true.)
      logical BW
      parameter (BW=.true.)
      real * 8 epsilon
      parameter (epsilon=1d-15)
      real * 8 pt1cut,pt2cut,pt1,pt2,m2jjmin
      real *8 cthdec,phidec
      integer psgen
      parameter (psgen=2)
c     1=no imp sampling at all on mjj^2
c     2=imp sampling on mjj^2
cccccccccccccccccc
      double precision pt2j1,pt2j2,pt2Z,m2inv56
      double precision dotp,supp
      double precision p1(4),p2(4),dR12,rsepn_p

      double precision PScut_ptj,PScut_mjj,PScut_dRjj
      save PScut_ptj,PScut_mjj,PScut_dRjj
      double precision powheginput
      external powheginput
cccccccccccccccccc

      if(ini) then
c     set initial- and final-state masses for Born and real
         do ileg=1,nlegborn
            kn_masses(ileg)=0.
         enddo
         kn_masses(nlegreal)=0.


         if(powheginput('#bornktmin').ne.-1d6) then
            write(*,*) 'Error: bornktmin not allowed in Zjj'
            write(*,*) 'Suppression factors on Born phase space'//
     &           ' should be used instead'
            write(*,*) 'please set supp_m2inv56 and supp_pt2j'
            call exit(-1)
         endif

         write(*,*) 'PSGEN= ',psgen

ccccccccccccccccccccccccccccccc
c     ISR generation cuts
ccccccccccccccccccccccccccccccc
         PScut_ptj=powheginput('#PScut_ptj')
         if(PScut_ptj.lt.0) PScut_ptj=1d0

c$$$         write(*,*) '************************************'
c$$$         write(*,*) '    GENERATION CUTS'
c$$$         write(*,*) 'Generation cut on UB pt_j =  ',PScut_ptj

cccccccccccccccccccccccccccccc
c     FSR generation cuts
cccccccccccccccccccccccccccccc
         PScut_mjj=-1d0
         PScut_mjj=powheginput('#PScut_mjj')
         if(PScut_mjj.lt.0) PScut_mjj=1d0
         PScut_dRjj=-1d0
         PScut_dRjj=powheginput('#PScut_dRjj')
         if((PScut_mjj.lt.0d0).and.(PScut_dRjj.lt.0d0)) then
c$$$            write(*,*) '**********************'
c$$$            write(*,*) '    GENERATION CUTS'
            write(*,*) 'Error in PS: no generation cuts for FSR set'
            write(*,*) 'Set either PScut_mjj or PScut_dRjj'
            call exit(-1)
         elseif((PScut_mjj.gt.0d0).and.(PScut_dRjj.lt.0d0)) then
c$$$            write(*,*) '**********************'
c$$$            write(*,*) '    GENERATION CUTS'
c$$$            write(*,*) 'Generation cut on mjj = ',PScut_mjj
         elseif((PScut_mjj.lt.0d0).and.(PScut_dRjj.gt.0d0)) then
c$$$            write(*,*) '**********************'
c$$$            write(*,*) '    GENERATION CUTS'
            write(*,*) 'Generation cut on dRjj = ',PScut_dRjj
            write(*,*) 'Error: cut on dRjj is an untested feature'
            call exit(-1)
         else
            write(*,*) '**********************'
            write(*,*) '    GENERATION CUTS'
            write(*,*) 'Error in PS: 2 generation cuts for FSR set'
            call exit(-1)
         endif
         write(*,*) '************************************'
         ini=.false.
      endif
      
      Vmass2 = ph_Zmass2
      Vmass2low = ph_Zmass2low
      Vmass2high = ph_Zmass2high
      VmVw = ph_ZmZw

c     xborn(1): s_ll (m2) distributed between
c     Vmass2low and Vmass2high, according to a BW shape
      zlow=atan((Vmass2low  - Vmass2)/VmVw)
      zhigh=atan((min(Vmass2high,kn_sbeams)  - Vmass2)/VmVw)
      z=zlow+(zhigh-zlow)*xborn(1)
      xjac=zhigh-zlow
      m2=VmVw*tan(z)+Vmass2
c     Jac=dm2/2/pi
      xjac=xjac*VmVw/cos(z)**2/(2*pi)

ccccccccccccccccccc
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
      if(psgen.eq.1) then
         m2jj = xborn(4)*(sqrt(s)-sqrt(m2))**2
         xjac = xjac * (sqrt(s)-sqrt(m2))**2
      elseif(psgen.eq.2) then
         m2jj = (xborn(4))**2 *(sqrt(s)-sqrt(m2))**2
         xjac = xjac * 2 * xborn(4) * (sqrt(s)-sqrt(m2))**2
      else
         write(*,*) 'Error: wrong psgen'
         call exit(-1)
      endif

c     pZ in the CM frame
      pV(0)=(s+m2-m2jj)/(2*sqrt(s))
      pVmod2 = pV(0)**2-m2
      
      if (pVmod2.lt.epsilon) then
         pVmod = 0d0
      else
         pVmod = sqrt(pVmod2)
      endif
      
      z=1-2*xborn(5)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)
      
      pV(1) = pVmod*sqrt(1-cth**2)
      pV(2) = 0d0
      pV(3) = pVmod*cth
c     now pV is the Z momentum in the CM frame

c     supply 2 pi for azimuthal integration (not performed)
      xjac=xjac*2*pi
c     supply the other factors to the jacobian
c     factor for the two-body jet phase space
      xjac=xjac/(8*(2*pi)**2)
c     factor for V production
      xjac=xjac/(4*(2*pi)**3)*pVmod/sqrt(s)
ccccccccccccccccccc
c     decay variables
      cthdec=1-2*xborn(8)
      kn_cthdec=cthdec
      xjac=xjac*2
      phidec=2*pi*xborn(9)
      xjac=xjac*2*pi
      xjac=xjac/(8*(2*pi)**2)
ccccccccccccccccc
c     Build kinematics
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
c     boost decay products in the CM frame, using
c     informations already used to compute pV
      beta=pVmod/pV(0)
      vec(1)=sqrt(1-cth**2)
      vec(2)=0d0
      vec(3)=cth
      call mboost(2,vec,beta,kn_cmpborn(0,3),kn_cmpborn(0,3))
      
c     boost back in the lab frame
c     now boost everything along 3rd axis
      betaCM=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(2,vec,betaCM,kn_cmpborn(0,3),kn_pborn(0,3))      

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

c     boost in the lab frame
c     compute first p_plus+p_minus-pV
      do mu=0,3
         pcmjj(mu)= kn_pborn(mu,1) + kn_pborn(mu,2) 
     $        - (kn_pborn(mu,3)+kn_pborn(mu,4))
      enddo
      pcmjjmod = sqrt(pcmjj(1)**2+pcmjj(2)**2+pcmjj(3)**2)
c     recompute pcmjj(0) from m2jj, otherwise there are points where 
c     beta > 1 or beta < 0
      pcmjj(0) = sqrt(m2jj+pcmjjmod**2)

      beta=pcmjjmod/pcmjj(0)

      do mu=1,3
         vec(mu)=pcmjj(mu)/pcmjjmod
      enddo

      call mboost(2,vec,beta,pJ(0,1),kn_pborn(0,5))

      if (check) then
         call checkmomzero(nlegborn,kn_pborn)
      endif
      
c     boost in the CM frame
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,-betaCM,kn_pborn(0,1),kn_cmpborn(0,1))      

      kn_minmass=sqrt(Vmass2low)


cccccccccccccccccccccccccccccccccccccccccccccccc
c     Beginning of PS cuts

      pt2j1=kn_pborn(1,5)**2+kn_pborn(2,5)**2
      pt2j2=kn_pborn(1,6)**2+kn_pborn(2,6)**2
      pt2Z=(kn_pborn(1,3)+kn_pborn(1,4))**2+
     $     (kn_pborn(2,3)+kn_pborn(2,4))**2
      m2inv56=2*dotp(kn_pborn(0,5),kn_pborn(0,6))
      p1(4)=kn_pborn(0,5)
      p2(4)=kn_pborn(0,6)
      do mu=1,3
         p1(mu)=kn_pborn(mu,5)
         p2(mu)=kn_pborn(mu,6)
      enddo
c$$$      dR12 =rsepn_p(p1,p2)


ccccccccccccccccccccccccccccccccccccccccccccc
c     ISR singularities at the Born level
ccccccccccccccccccccccccccccccccccccccccccccc
      if(
     $     (sqrt(pt2j1).lt.PScut_ptj).or.
     $     (sqrt(pt2j2).lt.PScut_ptj) ) then
        kn_jacborn=0.
      endif

ccccccccccccccccccccccccccccccccccccccccccccc
c     FSR singularity at the Born level
ccccccccccccccccccccccccccccccccccccccccccccc
c     mjj cut
      if(PScut_mjj.gt.0d0) then
         if(sqrt(m2inv56).lt.PScut_mjj) then
            kn_jacborn=0.
         endif
c     dRjj
      elseif(PScut_dRjj.gt.0d0) then
c     Experimental feature. It should never enter here
         write(*,*) 'Error in phase space,PScut_dRjj'
         call exit(-1)
         if(dR12.lt.PScut_dRjj) then
            kn_jacborn=0.
         endif
      endif

c     End of PS cuts
cccccccccccccccccccccccccccccccccccccccccccccccc
         
      end


      function mass(p)
      implicit none
      real * 8 p(0:3),mass
      mass = sqrt(abs(p(0)**2-p(1)**2-p(2)**2-p(3)**2))
      end


      subroutine born_suppression_2(fact)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      include '../include/pwhg_kn.h'
      double precision fact
      logical ini,suppression_enabled
      data ini/.true./
      double precision pt2j1,pt2j2,pt2Z,m2inv56
      double precision cpt2j,cpt2Z,cm2inv56,powpt2j,powpt2Z,powm2inv56
      save cpt2j,cpt2Z,cm2inv56,powpt2j,powpt2Z,powm2inv56
      save ini,suppression_enabled
      double precision dotp,supp
      double precision powheginput
      external powheginput
      integer mu
      double precision p1(4),p2(4),dR12,rsepn_p

      suppression_enabled=powheginput("#bornsuppfact").gt.0d0

      if (ini) then
         if(suppression_enabled) then
            write(*,*) '*****************************'
            write(*,*) '    SUPPRESSION FACTOR ON'
            write(*,*) '    supp = (Q2 / (Q2+Lam2) )**pow'
            
            cpt2j = powheginput('#supp_pt2j')
            if(cpt2j.lt.0d0) cpt2j = ( 10. )**2
            cm2inv56 = powheginput('#supp_m2inv56')
            if(cm2inv56.lt.0d0) cm2inv56 = ( 5.  )**2
            powpt2j = powheginput('#supp_powpt2j')
            if(powpt2j.lt.0d0) powpt2j = 2.
            powm2inv56 = powheginput('#supp_powm2inv56')
            if(powm2inv56.lt.0d0) powm2inv56 = 2.
            
            write(*,*) '      Lam2_ptj= ',cpt2j
            write(*,*) '      pow_ptj= ',powpt2j
            write(*,*) '      Lam2_mjj= ',cm2inv56
            write(*,*) '      pow_mjj= ',powm2inv56
            
            write(*,*) '*****************************'
            
         else
            write(*,*) '*****************************'
            write(*,*) '    NO BORN SUPPRESSION FACTOR'
            write(*,*) ' Error'
            call exit(-1)
         endif
         ini=.false.
      endif

      
      if(.not.suppression_enabled) then
         fact=1d0 !:
         return
      else
         pt2j1=kn_pborn(1,5)**2+kn_pborn(2,5)**2
         pt2j2=kn_pborn(1,6)**2+kn_pborn(2,6)**2
         pt2Z=(kn_pborn(1,3)+kn_pborn(1,4))**2+
     $        (kn_pborn(2,3)+kn_pborn(2,4))**2
         m2inv56=2*dotp(kn_pborn(0,5),kn_pborn(0,6))
      
         p1(4)=kn_pborn(0,5)
         p2(4)=kn_pborn(0,6)
         do mu=1,3
            p1(mu)=kn_pborn(mu,5)
            p2(mu)=kn_pborn(mu,6)
         enddo
         
c$$$         dR12 =rsepn_p(p1,p2)

         fact=
     $        supp(pt2j1,cpt2j,powpt2j)*
     $        supp(pt2j2,cpt2j,powpt2j)*
     $        supp(m2inv56,cm2inv56,powm2inv56)
      endif


      end

      function supp(x,fact,pow)
      implicit none
      double precision x,fact,pow,supp
      supp=(x/(x+fact))**pow
      end


c$$$      subroutine born_suppression_NO(fact)
c$$$      implicit none
c$$$      include 'PhysPars.h'
c$$$      include 'nlegborn.h'
c$$$      include 'pwhg_flst.h'
c$$$      include 'pwhg_kn.h'
c$$$c      include 'coupl.inc'
c$$$      real * 8 fact,zmass
c$$$      real * 8 pt5,pt6,pt56,dotp
c$$$      integer il,ilbar
c$$$      il=3
c$$$      ilbar=4
c$$$      zmass=ph_Zmass
c$$$      pt5=kn_cmpborn(1,il)**2+kn_cmpborn(2,il)**2
c$$$      pt6=kn_cmpborn(1,ilbar)**2+kn_cmpborn(2,ilbar)**2
c$$$      pt56=dotp(kn_cmpborn(0,il),kn_cmpborn(0,ilbar))
c$$$      fact=((1/zmass**2)/(1/pt5+1/pt6+1/pt56+1/zmass**2))**2
c$$$      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      double precision muf,mur
      integer mu
      double precision pZ(0:3),ptZ,pt5,pt6,Ht

      logical ini
      data ini/.true./

      double precision powheginput
      external powheginput
      integer choice
      save choice
c      parameter (choice=2)


      if (ini) then
         choice = powheginput('scalechoice')
         if(choice.eq.1) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to Z mass         '
            write(*,*) '*************************************'
         elseif(choice.eq.2) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to Ht/2'
            write(*,*) '    (defined in set_fac_ren_scales)'
            write(*,*) '*************************************'
         elseif(choice.eq.3) then
            write(*,*) '*************************************'
            write(*,*) '    Factorization and renormalization '
            write(*,*) '    scales set to MtZ'
            write(*,*) '    (defined in set_fac_ren_scales)'
            write(*,*) '*************************************'
         else
            write(*,*) 'Unsupported scale choice'
            call exit(-1)
         endif
         ini=.false.
      endif
      
      if(choice.eq.1) then 
         muf=ph_Zmass
         mur=ph_Zmass
      elseif(choice.eq.2) then
         do mu=0,3
            pZ(mu)=kn_cmpborn(mu,3)+kn_cmpborn(mu,4)
         enddo
         ptZ=sqrt(pZ(1)**2+pZ(2)**2)
         pt5=sqrt(kn_cmpborn(1,5)**2+kn_cmpborn(2,5)**2)
         pt6=sqrt(kn_cmpborn(1,6)**2+kn_cmpborn(2,6)**2)

         Ht=sqrt(ph_Zmass**2+ptZ**2)+pt5+pt6

         muf=Ht/2d0
         mur=Ht/2d0
      elseif(choice.eq.3) then
         do mu=0,3
            pZ(mu)=kn_cmpborn(mu,3)+kn_cmpborn(mu,4)
         enddo
         ptZ=sqrt(pZ(1)**2+pZ(2)**2)

         muf=sqrt(ph_Zmass**2+ptZ**2)
         mur=sqrt(ph_Zmass**2+ptZ**2)
      endif
      end

