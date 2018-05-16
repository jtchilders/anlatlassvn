      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      real * 8 xborn(ndiminteg-3),tmpvec(0:3)
      integer k
      logical ini,fullphsp
      data ini/.true./
      save ini,fullphsp
      real * 8 powheginput
      real * 8 xjac,smin,smax,z,s,wt,sqrts,
     1         mllminsq,mllmaxsq,m3,m45,taumin,lntaum,
     2         tau,ymax,ycm,xx(2),p1(4),p2(4),p3(4),p4(4),p5(4),
     3         p12(4),p45(4),beta,vec(3),expon
      integer mu

      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0
         kn_masses(3)=ph_Hmass
         ini = .false.
      endif      

      sqrts = sqrt(kn_sbeams)
      xjac=1
c     First determine virtualities of the Higgs
      smin=ph_Hmass2low
      smax=ph_Hmass2high
      mllminsq=ph_Wmass2low
      mllmaxsq=ph_Wmass2high
c
      z=xborn(1)
c
      call breitw(z,smin,smax,ph_hmass,ph_hwidth,s,wt)
c breitw includes in wt a factor
c   ((s-ph_hmass)**2+(ph_hmass*ph_hwidth)**2)/ph_hmass*ph_hwidth
c Take it off
      xjac=xjac*wt/(pi)*ph_hmass*ph_hwidth/
     1     ((s-ph_hmass**2)**2+(ph_hmass*ph_hwidth)**2)

c If you want Passarino's shape, put it here

      m3=sqrt(s)
      smin=mllminsq
      smax=mllmaxsq
c the following better for Z/gamma
c      z=xborn(2)**4
c      xjac=xjac*4*xborn(2)**3
      z=xborn(2)

      call breitw(z,smin,smax,ph_Wmass,ph_Wwidth,s,wt)
      xjac=xjac*wt/(2*pi)
      m45=sqrt(s)  
      expon=1d0/5
      taumin = ((m3+m45)/sqrts)**2
      lntaum = dlog(taumin)      
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

      call phi1_2(xborn(5),xborn(6),p12,p3,p45,m3,m45,wt)
      xjac=xjac*wt

      call phi3m0(xborn(7),xborn(8),p45,p4,p5,wt)
      xjac=xjac*wt

      kn_pborn(0,1) = p1(4)
      kn_pborn(0,2) = p2(4)
      kn_pborn(0,3) = p3(4)
      kn_pborn(0,4) = p4(4)
      kn_pborn(0,5) = p5(4)


      kn_pborn(1:3,1) = p1(1:3)
      kn_pborn(1:3,2) = p2(1:3)
      kn_pborn(1:3,3) = p3(1:3)
      kn_pborn(1:3,4) = p4(1:3)
      kn_pborn(1:3,5) = p5(1:3)

      kn_jacborn = xjac

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

c      kn_minmass = sqrt(ph_Hmass2low+ph_Wmass2low)
      
      end

      subroutine born_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
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
         if (runningscales.eq.1) then
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=Ht  used for Bbar function   **'
            write(*,*) '**   muf=Ht  used for Bbar function   **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         elseif (runningscales.eq.2) then
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=muf=sqrt(pt1*pt2)            **'
            write(*,*) '****************************************'
            write(*,*) '****************************************'
         else
            write(*,*) '****************************************'
            write(*,*) '****************************************'
            write(*,*) '**   mur=muf=MH+MW                    **'
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
         muf=ph_hmass+ph_wmass
         mur=ph_hmass+ph_wmass
      endif
      end



