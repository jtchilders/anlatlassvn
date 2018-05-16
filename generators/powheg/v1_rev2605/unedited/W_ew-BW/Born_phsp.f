      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
c     lepton masses
      real *8 lepmass(3),decmass
      common/clepmass/lepmass,decmass
      real*8 alam,kn_cmpborn_vec
      real*8 powheginput
      external powheginput
      real * 8 xborn(ndiminteg-3)
      real * 8 m2,xjac,tau,y,beta,vec(3),cth,cthdec,phidec,s,
     # z,zhigh,zlow
      integer mu,k
      logical ini
      data ini/.true./
      save ini
      if(ini) then
c     set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0d0
         enddo
c WZGRAD edit start
      if(powheginput('idvecbos').eq.24)then
         kn_masses(3)=decmass
      elseif(powheginput('idvecbos').eq.-24)then
         kn_masses(4)=decmass
      endif
c WZGRAD edit end
         kn_masses(nlegreal)=0d0
         ini=.false.
      endif
c Phase space:
c 1 /(16 pi S) d m^2 d cth d y
      xjac=1d0/kn_sbeams/(16*pi)
      zlow=atan((ph_Wmass2low  - ph_Wmass2)/ph_WmWw)
      zhigh=atan((ph_Wmass2high  - ph_Wmass2)/ph_WmWw)
      z=zlow+(zhigh-zlow)*xborn(1)
      xjac=xjac*(zhigh-zlow)
      m2=ph_WmWw*tan(z)+ph_Wmass2
c d m^2 jacobian
      xjac=xjac*ph_WmWw/cos(z)**2
c d x1 d x2 = d tau d y;
      tau=m2/kn_sbeams
      s=kn_sbeams*tau
      kn_sborn=s
c ymax=|log(tau)|/2
      y=-(1d0-2d0*xborn(2))*log(tau)/2d0
      xjac=-xjac*log(tau)
c generation of cth
      z=1d0-2d0*xborn(3)
      xjac=xjac*2
      cth=1.5d0*(z-z**3/3)
      xjac=xjac*1.5d0*(1-z**2)
      kn_born_pt2=0d0
      cthdec=cth
      phidec=0d0
      kn_cthdec=cthdec
      kn_jacborn=xjac
c WZGRAD edit
      kn_cmpborn_vec=alam(sqrt(m2),kn_masses(3),kn_masses(4))
     $     /2d0/sqrt(m2)
      kn_jacborn=kn_jacborn*
     $     2d0/sqrt(m2)*kn_cmpborn_vec
c Build kinematics
      kn_xb1=sqrt(tau)*exp(y)
      kn_xb2=tau/kn_xb1
c decay products in their rest frame
c      kn_cmpborn(0,3)=sqrt(m2)/2
c      kn_cmpborn(0,4)=kn_cmpborn(0,3)
c      kn_cmpborn(3,3)=kn_cthdec*kn_cmpborn(0,3)
c      kn_cmpborn(1,3)=sqrt(1-kn_cthdec**2)*cos(phidec)*kn_cmpborn(0,3)
c      kn_cmpborn(2,3)=sqrt(1-kn_cthdec**2)*sin(phidec)*kn_cmpborn(0,3)
c      kn_cmpborn(1,4)=-kn_cmpborn(1,3)
c      kn_cmpborn(2,4)=-kn_cmpborn(2,3)
c      kn_cmpborn(3,4)=-kn_cmpborn(3,3)
c include lepton mass
c WZGRAD edit start
      kn_cmpborn(0,3)=sqrt(m2)/2+
     $     (kn_masses(3)**2-kn_masses(4)**2)/2d0/sqrt(m2)
      kn_cmpborn(0,4)=sqrt(m2)/2+
     $     (kn_masses(4)**2-kn_masses(3)**2)/2d0/sqrt(m2)
      kn_cmpborn(3,3)=kn_cthdec*kn_cmpborn_vec
      kn_cmpborn(1,3)=sqrt(1-kn_cthdec**2)*cos(phidec)*kn_cmpborn_vec
      kn_cmpborn(2,3)=sqrt(1-kn_cthdec**2)*sin(phidec)*kn_cmpborn_vec
      kn_cmpborn(1,4)=-kn_cmpborn(1,3)
      kn_cmpborn(2,4)=-kn_cmpborn(2,3)
      kn_cmpborn(3,4)=-kn_cmpborn(3,3)
c WZGRAD edit end
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
c
c      call checkmomzero(nlegborn,kn_pborn)
c      call checkmass(2,kn_pborn(0,3))
c      do mu=1,4
c         write(6,*)dsqrt(dabs(kn_pborn(0,mu)**2-kn_pborn(1,mu)**2-
c     $        kn_pborn(2,mu)**2-kn_pborn(3,mu)**2)),decmass
c      enddo
c      stop
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
         if(powheginput('#runningscale').eq.1) then
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
               write(*,*) ' scales set to the W virtuality '
            else 
               write(*,*) "runningscale value not allowed"
               call exit(1)
            endif
            write(*,*) '*************************************'
            ini=.false.
         endif
         muref=sqrt(2d0*dotp(kn_pborn(0,3),kn_pborn(0,4)))
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
      end
c****************************************************************
      subroutine LO_observables
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'

      real*8 thetal,thetan,etan
      real*8 ptl,ptn,etal
      common/observsLO/ptl,ptn,etal

      integer calocuts
      real*8 p1(4),p2(4),p3(4)
      real*8 sp1(4),sp2(4),sp3(4)
      real*8 powheginput
      external powheginput
c smearing and recombination
      calocuts = powheginput('calo')
      if(calocuts.eq.1)then
         call momentumprep_LO(p1,p2,p3)
         call smear(p1,p2,p3,sp1,sp2,sp3)
         ptl = dsqrt(sp1(1)**2+sp1(2)**2)
         ptn = dsqrt(sp2(1)**2+sp2(2)**2)

         thetal=atan2(ptl,sp1(3))
         thetan=atan2(ptn,sp2(3))

         etal=-dlog(dtan(thetal/2d0))
         etan=-dlog(dtan(thetan/2d0))
      else
         ptl=dsqrt(kn_pborn(1,3)**2+kn_pborn(2,3)**2)
         ptn=dsqrt(kn_pborn(1,4)**2+kn_pborn(2,4)**2)

         thetal=atan2(ptl,kn_pborn(3,3))
         thetan=atan2(ptn,kn_pborn(3,4))
         
         etal=-dlog(dtan(thetal/2d0))
         etan=-dlog(dtan(thetan/2d0))
      endif
      end
c****************************************************************
c***********************************************************************
      subroutine momentumprep_LO(p1,p2,p3)
c.....
c     this subroutine prepares the momenta for smearing - only for LO
c     input: kn_born from common/pwhg_kn 
c     output: p1,p2 of common/smomenta
c.....
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      real*8 p1(4),p2(4),p3(4)
      real*8 powheginput
      external powheginput
!     =========================================
!     assign FS momenta     
!     =========================================
      if(powheginput('idvecbos').eq.24)then
!     lepton
      p1(1) = kn_pborn(1,3)       
      p1(2) = kn_pborn(2,3)       
      p1(3) = kn_pborn(3,3)       
      p1(4) = kn_pborn(0,3)
!     neutrino
      p2(1) = kn_pborn(1,4)       
      p2(2) = kn_pborn(2,4)       
      p2(3) = kn_pborn(3,4)       
      p2(4) = kn_pborn(0,4)       
      elseif(powheginput('idvecbos').eq.-24)then
!     lepton
      p1(1) = kn_pborn(1,4)       
      p1(2) = kn_pborn(2,4)       
      p1(3) = kn_pborn(3,4)       
      p1(4) = kn_pborn(0,4)       
!     neutrino
      p2(1) = kn_pborn(1,3)       
      p2(2) = kn_pborn(2,3)       
      p2(3) = kn_pborn(3,3)       
      p2(4) = kn_pborn(0,3)       
      endif
      p3(1) = 0d0
      p3(2) = 0d0
      p3(3) = 0d0
      p3(4) = 0d0
      end
c***********************************************************************
      real*8 function alam(m,m1,m2)
      implicit none
      real*8 m,m1,m2
      real*8 s,s1,s2
      real*8 aux
      s=m**2
      s1=m1**2
      s2=m2**2
      aux=s**2+s1**2+s2**2-2*s*s1-2*s*s2-2*s1*s2
      if(aux.lt.0.d0)aux=0.d0	
      alam=dsqrt(aux)
      return
      end	
