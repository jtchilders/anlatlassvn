      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'pwhg_flg.h'
      real * 8 xborn(ndiminteg-3)
      integer k
      logical ini
      data ini/.true./
      save ini
      real * 8 xjac,smin,smax,z,s,wt,sqrts,
     1     m3,m45,taumin,lntaum,
     2     tau,ymax,ycm,xx(2),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     3     p12(4),p45(4),p345(4),beta,vec(3),s345min,s345max,s345,
     4     sratio,lnsratio,expon
      integer mu

      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         kn_masses(3)=ph_Hmass
         ph_HmHw=ph_Hmass*ph_Hwidth
         ini = .false.
      endif      

      sqrts = sqrt(kn_sbeams)
      xjac=1
c     First determine virtuality of the Higgs
      smin=ph_Hmass2low
      smax=ph_Hmass2high
      call breitw(xborn(1),smin,smax,ph_hmass,ph_hwidth,s,wt)
c breitw includes in wt a factor
c   ((s-ph_hmass)**2+(ph_hmass*ph_hwidth)**2)/ph_hmass*ph_hwidth
c Take it off
      xjac=xjac*wt/(pi)*ph_hmass*ph_hwidth/
     1     ((s-ph_hmass**2)**2+(ph_hmass*ph_hwidth)**2)

c If you want Passarino's shape, put it here

      m3=sqrt(s)
c      write(*,*)"--> mH mass in PWHG: ",m3

      smin=ph_Zmass2low
      smax=ph_Zmass2high
c the following better for Z/gamma
c      z=xborn(2)**4
c      xjac=xjac*4*xborn(2)**3
      z=xborn(2)

      call breitw(z,smin,smax,ph_zmass,ph_zwidth,s,wt)
      xjac=xjac*wt/(2*pi)
      m45=sqrt(s)
      
      s345min=(m3+m45)**2 
      s345max = (sqrts-kn_ktmin)**2 - kn_ktmin**2

      if (s345max.lt.s345min) then
         write(*,*) 'Too high values for H/W virtualities'
         write(*,*) 'Return Jacobian=0'
         kn_jacborn = 0
         return
      endif

c      s345 = (s345max - s345min)*xborn(5)**3 + s345min      
c      xjac = xjac*(s345max - s345min)*3*xborn(5)**2

      sratio = s345min/s345max
      lnsratio = log(sratio)
      expon=1d0/4
      s345 = exp(lnsratio*(xborn(5)**expon))*s345max
      xjac = xjac*(-lnsratio*s345)*expon*xborn(5)**(expon-1)

      expon=1d0/5
      taumin = (kn_ktmin + sqrt(s345 + kn_ktmin**2))**2/kn_sbeams
      lntaum = log(taumin)
      tau = exp(lntaum*(1d0-xborn(3))**expon)
      xjac = xjac*(-lntaum*tau)*expon*(1d0-xborn(3))**(expon-1)

      kn_sborn = kn_sbeams*tau

      ymax=-0.5d0*log(tau)
      ycm=xborn(4)*2*ymax-ymax
      xjac = xjac*2*ymax
      xx(1)=sqrt(tau)*exp(ycm)
      xx(2)=tau/xx(1)

c---if x's out of normal range abort
      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. 1d-8)
     & .or. (xx(2) .lt. 1d-8)) then
         write(*,*) ' error in Born phase space!, x1,x2 our of range'
         write(*,*) xx(1),xx(2)
         kn_jacborn = 0
         return
      endif

      p1(4)=xx(1)*sqrts*0.5d0
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=xx(1)*sqrts*0.5d0

      p2(4)=xx(2)*sqrts*0.5d0
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-xx(2)*sqrts*0.5d0

C     total incoming momentum 
      p12 = p1+p2 

c      s345min=(m3+m45)**2

      call phi1_2m_nobw_ktmin(p12,0d0,sqrt(s345),kn_ktmin,
     $     xborn(6),xborn(7),
     $     p6,p345,wt)
      xjac=xjac*wt


      if (sqrt(p6(1)**2+p6(2)**2).lt.kn_ktmin) then
         write(*,*) '** ERROR in phi1_2m_nobw_ktmin: p6_T< kn_ktmin **'
         write(*,*) sqrt(p6(1)**2+p6(2)**2),' < ',kn_ktmin         
         write(*,*) 'The POWHEG BOX continues'
      endif


      call phi1_2(xborn(8),xborn(9),p345,p3,p45,m3,m45,wt)
      xjac=xjac*wt

      call phi3m0(xborn(10),xborn(11),p45,p4,p5,wt)
      xjac=xjac*wt

      kn_pborn(0,1) = p1(4)
      kn_pborn(0,2) = p2(4)
      kn_pborn(0,3) = p3(4)
      kn_pborn(0,4) = p4(4)
      kn_pborn(0,5) = p5(4)
      kn_pborn(0,6) = p6(4)

      kn_pborn(1:3,1) = p1(1:3)
      kn_pborn(1:3,2) = p2(1:3)
      kn_pborn(1:3,3) = p3(1:3)
      kn_pborn(1:3,4) = p4(1:3)
      kn_pborn(1:3,5) = p5(1:3)
      kn_pborn(1:3,6) = p6(1:3)

      kn_jacborn = xjac/(2d0*pi)

c     now boost everything BACK along z-axis 
      kn_xb1 = xx(1)
      kn_xb2 = xx(2)
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=-1
      call mboost(nlegborn-2,vec,beta,kn_pborn(:,3:),
     1     kn_cmpborn(:,3:))
      do mu=0,3
         kn_cmpborn(mu,1)=sqrt(kn_xb1*kn_xb2)*kn_beams(mu,1)
         kn_cmpborn(mu,2)=sqrt(kn_xb1*kn_xb2)*kn_beams(mu,2)
      enddo
      kn_minmass = sqrt(ph_Hmass2low+ph_Zmass2low)

      end

      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      real * 8 fact,ptmin,ptminZ
      real * 8 pt2,pt2Z
      logical ini
      data ini/.true./
      real * 8 powheginput
      save ini,ptmin,ptminZ    
      if (ini) then
         ptmin=powheginput("#bornsuppfact")      
         ptminZ=powheginput("#bornsuppfactZ")      
         if (ptmin.lt.0d0) then
            ptmin=0d0
         endif
         if (ptminZ.lt.0d0) then
            ptminZ=0d0
         endif
         ini=.false.
      endif
      if(flg_weightedev) then
         pt2=kn_cmpborn(1,6)**2+kn_cmpborn(2,6)**2
         fact = pt2/(ptmin**2+pt2)
         pt2Z=(kn_cmpborn(1,4)+kn_cmpborn(1,5))**2+
     $        (kn_cmpborn(2,4)+kn_cmpborn(2,5))**2
         fact=fact*(pt2Z+1d0)/(pt2Z+1d0+ptminZ**2)
      else
         fact=1
      endif
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      real * 8 muf,mur
      logical ini
      data ini/.true./
      integer runningscales
      real * 8 pt1,pt2,ptHsq,Ht
      real * 8 powheginput
      external powheginput      
      save ini,runningscales

      if (ini) then
         runningscales=powheginput("#runningscales")

         if(powheginput("#minlo").eq.1) then
            write(*,*) '****************************************'
            write(*,*) '*******          MINLO ACTIVE    *******'
            write(*,*) '****************************************'
            write(*,*) '*******     FIXED SCALES!          *****'
            runningscales=0
         endif

         if (runningscales.eq.1) then
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '** mur=sqrt(MH^2+pT_H^2)+pT_Z         **'
            write(*,*) '** muf=mur   used for Bbar function   **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         elseif (runningscales.eq.2) then
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=muf=sqrt(pT_lep1*pT_lep2)    **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         else
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=muf=MH+MZ                    **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'            
         endif
         ini=.false.
      endif
      
      if (runningscales.eq.1) then
         pt1=sqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)
         pt2=sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)
         ptHsq=kn_pborn(1,3)**2+kn_pborn(2,3)**2
         Ht=sqrt(ph_hmass**2+ptHsq)+pt1+pt2
         mur=Ht
         muf=mur
      elseif (runningscales.eq.2) then
         pt1=sqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)
         pt2=sqrt(kn_pborn(1,5)**2+kn_pborn(2,5)**2)
         mur=sqrt(pt1*pt2)
         if(mur.lt.2) mur=2
         muf=mur         
      else
         muf=ph_hmass+ph_zmass
         mur=ph_hmass+ph_zmass
      endif
      end



