      subroutine real_EW(flavor,res)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_wzgrad.h'
      integer flavor(5),i,j
      real*8 res,sig(2),sigh(3,2)
      real*8 epmin,ephoton,smin,mes
      integer flag1,flag2
      common/leptonidf/flag1,flag2
      real*8 powheginput
      external powheginput

      if(kn_emitter.eq.0)then
         flag1=0 !flag for hard-noncollinear PS
         flag2=0 !flag for dR cut
!        =================================
!        call in all pieces: sigh(1:3,1:2)
!        =================================
!        first calculate the invariants from input momenta
!        they will be stored in pwhg_wzgrad.h
         call dotprods
!        separately call IS, FS and INT contributions to MES 
         do j=1,3
         call mathard_cc(j,sig)
         sigh(j,1)=sig(1)
         sigh(j,2)=sig(2)
         enddo
!        =================================
!        do the kinematic cuts 
!        =================================
!        hard non-collinear cut parameters
         epmin=deltas*dsqrt(rl_sinv(1,2))/2d0
         if(test(3).eq.2.and.test(4).eq.2)epmin=deltas*mw/2d0
         ephoton=kn_cmpreal(0,5) !photon energy in COM frame!
         smin=deltac*dsqrt(rl_sinv(1,2))*ephoton
!        cut 1 ==> is photon soft?
!        if photon is soft, all contributions must vanish (IS, FS, INT)
         if(ephoton.lt.epmin)then
            do j=1,3
               sigh(j,1)=0d0
               sigh(j,2)=0d0
            enddo
            flag1=1
            goto 10
         endif 
!        cut2 ==> quark-photon angle too small?
!        this affects IS and INT contributions
!        and possibility of del_R cuts
         if(dabs(rl_sinv(1,5)).lt.smin.or.
     1      dabs(rl_sinv(2,5)).lt.smin)then
         sigh(1,1)=0d0
         sigh(1,2)=0d0
         sigh(3,1)=0d0
         sigh(3,2)=0d0
         flag1=1 
         endif
!        cut 3 ==> lepton-photon angle too small?
!        this affects only FS and INT contributions
!        and possibility of del_R cuts
         if(collcut.eq.1)then
         smin=smin*(1d0-2d0*ephoton/dsqrt(rl_sinv(1,2)))
            if(dabs(rl_sinv(3,5)).lt.smin)then
            sigh(2,1)=0d0
            sigh(2,2)=0d0
            sigh(3,1)=0d0
            sigh(3,2)=0d0
            flag1=1 
            endif
         endif
 10      continue
c dR cut
c smearing and application of dR cut is done in NLO_observables
         if(powheginput('calo').eq.1)then 
            call NLO_observables
            if(flag2.eq.1)then
               res=0d0
               return
            endif
         endif

         if(flag1.eq.1.and.test(3).eq.4)then
            res=0d0
            return
         endif
!        ============================================
!        put result together wrt QED subset from user
!        ============================================
         do i=1,2
            if(flavor(i).eq.2.or.flavor(i).eq.1)then
               if(test(3).eq.1)then
               mes=sigh(1,i)
               elseif(test(3).eq.2)then
               mes=sigh(2,i)
               elseif(test(3).eq.3)then
               mes=sigh(3,i)
               elseif(test(3).eq.4)then
               mes=sigh(1,i)+sigh(2,i)+sigh(3,i)
               endif
            endif
         enddo

         if(powheginput('idvecbos').eq.24)then
         res=ph_CKM(1,3)**2*mes
         elseif(powheginput('idvecbos').eq.-24)then
         res=ph_CKM(2,1)**2*mes
         endif
      else
      res=0d0
      endif
      end
!***********************************************************************
      subroutine dotprods  
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_wzgrad.h'
      real*8 dotp,powheginput
      external dotp,powheginput

!     kn_preal(0:3,1)==>incoming +rapidity
!     kn_preal(0:3,2)==>incoming -rapidity
!     kn_preal(0:3,3)==>lepton
!     kn_preal(0:3,4)==>neutrino
!     kn_preal(0:3,5)==>photon

      rl_sinv(1,1)=2d0*dotp(kn_preal(0,1),kn_preal(0,1))
      rl_sinv(1,2)=2d0*dotp(kn_preal(0,1),kn_preal(0,2))

      if(powheginput('idvecbos').eq.24)then
      rl_sinv(1,3)=2d0*dotp(kn_preal(0,1),kn_preal(0,3))
      rl_sinv(1,4)=2d0*dotp(kn_preal(0,1),kn_preal(0,4))
      rl_sinv(2,3)=2d0*dotp(kn_preal(0,2),kn_preal(0,3))
      rl_sinv(2,4)=2d0*dotp(kn_preal(0,2),kn_preal(0,4))
      rl_sinv(3,3)=2d0*dotp(kn_preal(0,3),kn_preal(0,3))
      rl_sinv(3,4)=2d0*dotp(kn_preal(0,3),kn_preal(0,4))
      rl_sinv(4,4)=2d0*dotp(kn_preal(0,4),kn_preal(0,4))
      rl_sinv(3,5)=2d0*dotp(kn_preal(0,3),kn_preal(0,5))
      rl_sinv(4,5)=2d0*dotp(kn_preal(0,4),kn_preal(0,5))
      elseif(powheginput('idvecbos').eq.-24)then
      rl_sinv(1,3)=2d0*dotp(kn_preal(0,1),kn_preal(0,4))
      rl_sinv(1,4)=2d0*dotp(kn_preal(0,1),kn_preal(0,3))
      rl_sinv(2,3)=2d0*dotp(kn_preal(0,2),kn_preal(0,4))
      rl_sinv(2,4)=2d0*dotp(kn_preal(0,2),kn_preal(0,3))
      rl_sinv(3,3)=2d0*dotp(kn_preal(0,4),kn_preal(0,4))
      rl_sinv(3,4)=2d0*dotp(kn_preal(0,4),kn_preal(0,3))
      rl_sinv(4,4)=2d0*dotp(kn_preal(0,3),kn_preal(0,3))
      rl_sinv(3,5)=2d0*dotp(kn_preal(0,4),kn_preal(0,5))
      rl_sinv(4,5)=2d0*dotp(kn_preal(0,3),kn_preal(0,5))
      endif

      rl_sinv(1,5)=2d0*dotp(kn_preal(0,1),kn_preal(0,5))

      rl_sinv(2,1)=rl_sinv(1,2)
      rl_sinv(2,2)=2d0*dotp(kn_preal(0,2),kn_preal(0,2))
      rl_sinv(2,5)=2d0*dotp(kn_preal(0,2),kn_preal(0,5))

      rl_sinv(3,1)=rl_sinv(1,3)
      rl_sinv(3,2)=rl_sinv(2,3)

      rl_sinv(4,1)=rl_sinv(1,4)
      rl_sinv(4,2)=rl_sinv(2,4)
      rl_sinv(4,3)=rl_sinv(3,4)

      rl_sinv(5,1)=rl_sinv(1,5)
      rl_sinv(5,2)=rl_sinv(2,5)
      rl_sinv(5,3)=rl_sinv(3,5)
      rl_sinv(5,4)=rl_sinv(4,5)
      rl_sinv(5,5)=2d0*dotp(kn_preal(0,5),kn_preal(0,5))

      end
c***********************************************************************
      subroutine momentumprep_MES(p1,p2,p3)
c.....
c     this subroutine prepares the momenta for smearing - only for MES
c     input: kn_preal from common/pwhg_kn 
c     output: p1,p2,p3 of common/smomenta_NLO
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
      p1(1) = kn_preal(1,3)       
      p1(2) = kn_preal(2,3)       
      p1(3) = kn_preal(3,3)       
      p1(4) = kn_preal(0,3)
!     neutrino
      p2(1) = kn_preal(1,4)       
      p2(2) = kn_preal(2,4)       
      p2(3) = kn_preal(3,4)       
      p2(4) = kn_preal(0,4)       
      elseif(powheginput('idvecbos').eq.-24)then
!     lepton
      p1(1) = kn_preal(1,4)       
      p1(2) = kn_preal(2,4)       
      p1(3) = kn_preal(3,4)       
      p1(4) = kn_preal(0,4)       
!     neutrino
      p2(1) = kn_preal(1,3)       
      p2(2) = kn_preal(2,3)       
      p2(3) = kn_preal(3,3)       
      p2(4) = kn_preal(0,3)       
      endif
!     photon
      p3(1) = kn_preal(1,5)
      p3(2) = kn_preal(2,5)
      p3(3) = kn_preal(3,5)
      p3(4) = kn_preal(0,5)
      end
c***********************************************************************
