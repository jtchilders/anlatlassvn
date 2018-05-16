c ----------------------------------------------
c WWA Vertices 
c ---------------------------------------------


      subroutine wwa_anomal(q1,q2,qa,wwa)
c Calculates anomalous WWA couplings
c  q1: momentum of  W+, q1(mu)
c  q2: momentum of the W-,q2(nu)
c  qa: momentum of the photon,qa(rho)
c  wwa: leptonic tensor      
      implicit none
      double precision q1(0:3),q2(0:3),qa(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu,epsrmunurho
      integer mu,nu,rho,i
      complex*16 wwa(0:3,0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu,epsrmunurho

      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"

      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)      
      im=(0.,1.)
      do mu=0,3
         do nu=0,3
            do rho=0,3
               wwa(mu,nu,rho)=(0.d0,0.d0)
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               if (fbw .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fbw*gw*wmass**2*sinw*
     -                 (qa(mu)*gmunu(nu,rho)-qa(nu)*gmunu(mu,rho))
               endif
               if (fwww .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+3./2.*fwww*gw**3*sinw*
     &                 (q2(mu)*qa(nu)*q1(rho)-qa(mu)*q1(nu)*q2(rho)+
     -                 gmunu(nu,rho)*(qa(mu)*dotrr(q1,q2)-q2(mu)*
     &                 dotrr(q1,qa))-gmunu(mu,rho)*(qa(nu)*dotrr(q1,q2)-
     &                 q1(nu)*dotrr(q2,qa))+gmunu(mu,nu)*(q2(rho)*
     &                 dotrr(q1,qa)-q1(rho)*dotrr(q2,qa)))
               endif
               if (fw .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fw/2.*wmass**2*gw*sinw*
     &                 (qa(nu)*gmunu(mu,rho)-qa(mu)*gmunu(nu,rho))
               endif
               if (fb .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fb/2.*wmass**2*gw*sinw*
     &                 (qa(nu)*gmunu(mu,rho)-qa(mu)*gmunu(nu,rho))      
               endif
               if (fbtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)-fbtilde*wmass**2*gw*
     &                 sinw/2.*epsrmunurho(qa,mu,nu,rho)
               endif
               if (fbwtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fbwtilde*wmass**2*gw*
     &                 sinw*epsrmunurho(qa,mu,nu,rho)
               endif
               if (fwwtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)-fwwtilde*2.*gw*sinw*
     &                 wmass**2*(epsrmunurho(q2,mu,nu,rho)+
     &                 epsrmunurho(q1,mu,nu,rho)+
     &                 epsrmunurho(qa,mu,nu,rho))
               endif
               if (fwtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fwtilde*gw*sinw*
     &                 wmass**2*(epsrmunurho(q2,mu,nu,rho)+
     &                 epsrmunurho(q1,mu,nu,rho)+1./2.*
     &                 epsrmunurho(qa,mu,nu,rho))
               endif
               if (fwwwtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)-fwwwtilde/2.*gw**3*sinw*
     &                 (-epsrrrmu(q1,q2,qa,rho)*gmunu(mu,nu)-
     &                 epsrrrmu(q1,q2,qa,nu)*gmunu(mu,rho)-
     &                 epsrrmunu(q1,qa,nu,rho)*q2(mu)
     &                 -epsrrmunu(q1,q2,nu,rho)*qa(mu)-
     &                 epsrrrmu(q1,q2,qa,mu)*gmunu(nu,rho)
     &                 +epsrrmunu(q2,qa,mu,rho)*q1(nu)-
     &                 epsrrmunu(q1,q2,mu,rho)*qa(nu)+
     &                 epsrrmunu(q1,qa,mu,nu)*q2(rho)
     &                 +epsrrmunu(q2,qa,mu,nu)*q1(rho)-
     &                 epsrmunurho(qa,mu,nu,rho)*dotrr(q1,q2)-
     &                 epsrmunurho(q1,mu,nu,rho)*dotrr(q2,qa)-
     &                 epsrmunurho(q2,mu,nu,rho)*dotrr(q1,qa))
               endif
               if (fdwtilde .ne. 0) then
                  wwa(mu,nu,rho)=wwa(mu,nu,rho)+fdwtilde*2.*gw**3*sinw*
     &                 (epsrrmunu(qa,q2,nu,rho)*(qa(mu)-q2(mu))
     &                 +epsrrmunu(qa,q1,mu,rho)*(q1(nu)-qa(nu))+
     &                 epsrrmunu(q2,q1,mu,nu)*(q2(rho)-q1(rho))
     &                 +epsrmunurho(qa,mu,nu,rho)*(dotrr(qa,q2)+
     &                 dotrr(qa,q1))+epsrmunurho(q2,mu,nu,rho)*
     &                 (dotrr(qa,q2)+dotrr(q2,q1))+
     &                 epsrmunurho(q1,mu,nu,rho)*
     &                 (dotrr(qa,q1)+dotrr(q2,q1)))
               endif
            enddo
         enddo
      enddo

      end
      

            

      subroutine wwa_anomalwmin(q1,q2,j,r,mu,nu,wwa)
c Calculates anomalous WWA couplings when a current j is attached to a W-
c   q1: Momentum of th W+, mu
c   q2: Momentum of the A, nu
c   j: leptonic current
c   r: -q1-q2
c   mu,nu: Lorentz Indices of the leptonic tensor
c   wwa: component of the leptonic tensor  
    
      implicit none

      double precision q1(0:3), q2(0:3),r(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu
      double complex j(0:3,-1:1),wwa, dotrc,epscrmunu,epscrrmu,epscrrr,im
      integer mu, nu
      external dotrr, dotrc,epscrmunu,epsrrmunu,epsrrrmu,epscrrmu,epscrrr
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"

      im=(0.,1.)    
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)      
      wwa=(0.d0,0.d0)
      if (fbw .ne. 0) then
         wwa=wwa + fbw*gw*wmass**2*sinw*(-dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)+q2(mu)*j(nu,-1))
      endif
      if (fwww .ne. 0) then
         wwa= wwa + 3./2.*gw**3*sinw*fwww*(dotrc(q2,j(0,-1))*q1(nu)*
     -        r(mu)-q2(mu)*dotrc(q1,j(0,-1))
     &        *r(nu)+dotrr(q2,r)*(gmunu(mu,nu)*dotrc(q1,j(0,-1))
     &        -q1(nu)*j(mu,-1))+dotrr(r,q1)*(q2(mu)*j(nu,-1)-
     -        gmunu(mu,nu)*dotrc(q2,j(0,-1)))+
     &        dotrr(q1,q2)*(r(nu)*j(mu,-1)-r(mu)*j(nu,-1)))
      endif
      if (fw .ne. 0) then
         wwa = wwa + fw*gw*wmass**2*sinw/2.*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fb .ne. 0)then
         wwa = wwa +fb*gw*wmass**2*sinw/2.*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fbtilde .ne. 0) then
         wwa = wwa-fbtilde*wmass**2*gw*sinw/2.*
     -        epscrmunu(j(0,-1),q2,mu,nu)
      endif
      if (fbwtilde .ne. 0) then
         wwa = wwa+fbwtilde*wmass**2*gw*sinw*epscrmunu(j(0,-1),q2,mu,nu)
      endif      
      if (fwwtilde .ne. 0)then
         wwa =wwa +2.*fwwtilde*gw*sinw*wmass**2*
     -        (-epscrmunu(j(0,-1),r,mu,nu)
     &        -epscrmunu(j(0,-1),q1,mu,nu)-epscrmunu(j(0,-1),q2,mu,nu))
      endif 
      if (fwtilde .ne. 0) then
         wwa = wwa-fwtilde*gw*sinw*wmass**2*(-epscrmunu(j(0,-1),r,mu,nu)
     &        -epscrmunu(j(0,-1),q1,mu,nu)-1./2.*
     -        epscrmunu(j(0,-1),q2,mu,nu))
      endif
      if (fwwwtilde .ne. 0)then
         wwa = wwa + fwwwtilde/2.*sinw*gw**3*
     &        (-epsrrrmu(r,q1,q2,nu)*j(mu,-1)+epscrrr(j(0,-1),r,q1,q2)*
     -        gmunu(mu,nu)+epscrrmu(j(0,-1),q1,q2,nu)
     &        *r(mu)-epscrrmu(j(0,-1),r,q1,nu)*q2(mu)-
     -        epsrrrmu(r,q1,q2,mu)*j(nu,-1)
     &        -epsrrmunu(r,q2,mu,nu)*dotrc(q1,j(0,-1))-
     -        epsrrmunu(r,q1,mu,nu)*dotrc(q2,j(0,-1))+
     &        epscrrmu(j(0,-1),q1,q2,mu)*r(nu)+
     -        epscrrmu(j(0,-1),r,q2,mu)*q1(nu)+
     -        epscrmunu(j(0,-1),q2,mu,nu)*
     &        dotrr(r,q1)+epscrmunu(j(0,-1),q1,mu,nu)*dotrr(r,q2)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q1,q2))     
      endif       
      if (fdwtilde .ne. 0) then
         wwa = wwa +fdwtilde*2*gw**3*sinw*(-epscrrmu(j(0,-1),q2,r,nu)*
     -        (r(mu)-q2(mu))-epscrrmu(j(0,-1),q1,r,mu)
     &        *(q1(nu)-r(nu))-epsrrmunu(q1,q2,mu,nu)*dotrc(q1,j(0,-1))
     &        +epsrrmunu(q1,q2,mu,nu)*dotrc(q2,j(0,-1))+
     -        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,q2)+
     -        epscrmunu(j(0,-1),q2,mu,nu)
     &        *dotrr(q1,q2)+epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)
     &        *dotrr(q1,r)+epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q2,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q2,r))
      endif

      end


      
      subroutine wwa_anomalwplus(q1,q2,j,r,mu,nu,wwa)
c Calculates anomalous WWA couplings when a current j is attached to a W+
c   q1: Momentum of th W-, mu
c   q2: Momentum of the A, nu
c   j: leptonic current
c   r: -q1-q2
c   mu,nu: Lorentz Indices of the leptonic tensor
c   wwa: component of the leptonic tensor      
      implicit none
      double precision q1(0:3), q2(0:3),r(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu
      double complex j(0:3,-1:1),wwa, dotrc,epscrmunu,epscrrmu,epscrrr
      integer mu, nu
      external dotrr, dotrc,epscrmunu,epsrrmunu,epsrrrmu,epscrrmu,epscrrr
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"

      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)      
      wwa=(0.d0,0.d0)
      if (fbw .ne. 0) then
         wwa=wwa+fbw*gw*wmass**2*sinw*(dotrc(q2,j(0,-1))*gmunu(mu,nu)-
     -        q2(mu)*j(nu,-1))
      endif
      if (fwww .ne. 0) then
         wwa=wwa+3./2.*fwww*gw**3*sinw*(-j(mu,-1)*r(nu)*dotrr(q1,q2)+
     -        r(mu)*(j(nu,-1)*dotrr(q1,q2)-q1(nu)*dotrc(q2,j(0,-1)))+
     &        gmunu(mu,nu)*dotrc(q2,j(0,-1))*dotrr(q1,r)+q2(mu)*(r(nu)*
     &        dotrc(q1,j(0,-1))-j(nu,-1)*dotrr(q1,r))+j(mu,-1)*q1(nu)*
     &        dotrr(q2,r)-gmunu(mu,nu)*dotrc(q1,j(0,-1))*dotrr(q2,r))    
      endif
      if (fw .ne. 0) then
         wwa =wwa +fw*gw*wmass**2*sinw/2.*(-gmunu(mu,nu)*
     &        dotrc(q2,j(0,-1))+j(nu,-1)*q2(mu))
      endif
      if (fb .ne. 0) then
         wwa = wwa +fb*gw*wmass**2*sinw/2.*(-gmunu(mu,nu)*
     &        dotrc(q2,j(0,-1))+j(nu,-1)*q2(mu))
      endif
      if (fbtilde .ne. 0) then
         wwa = wwa+fbtilde*wmass**2*gw*sinw/2.*
     &        epscrmunu(j(0,-1),q2,mu,nu)
      endif
      if (fbwtilde .ne. 0) then
         wwa = wwa-fbwtilde*wmass**2*gw*sinw*epscrmunu(j(0,-1),q2,mu,nu)
      endif
      if (fwwtilde .ne. 0) then
         wwa = wwa +2.*fwwtilde*gw*sinw*wmass**2*
     &        (epscrmunu(j(0,-1),q1,mu,nu)
     &        +epscrmunu(j(0,-1),q2,mu,nu)+epscrmunu(j(0,-1),r,mu,nu))
      endif
      if (fwtilde .ne. 0) then
         wwa = wwa-fwtilde*gw*sinw*wmass**2*(epscrmunu(j(0,-1),q1,mu,nu)
     &        +1./2.*epscrmunu(j(0,-1),q2,mu,nu)+
     &        epscrmunu(j(0,-1),r,mu,nu))
      endif
      if (fwwwtilde .ne. 0) then
         wwa = wwa-fwwwtilde/2.*gw**3*sinw*(epscrrr(j(0,-1),q1,q2,r)*
     &        gmunu(mu,nu)-epsrrrmu(q1,q2,r,nu)
     &        *j(mu,-1)+epscrrmu(j(0,-1),q1,r,nu)*q2(mu)+
     &        epscrrmu(j(0,-1),q1,q2,nu)*r(mu)
     &        -epsrrrmu(q1,q2,r,mu)*j(nu,-1)-epscrrmu(j(0,-1),q2,r,mu)*
     &        q1(nu)+epscrrmu(j(0,-1),q1,q2,mu)*r(nu)+
     &        epsrrmunu(q2,r,mu,nu)*dotrc(q1,j(0,-1))+
     &        epsrrmunu(q1,r,mu,nu)*dotrc(q2,j(0,-1))
     &        +epscrmunu(j(0,-1),r,mu,nu)*dotrr(q1,q2)+
     &        epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q1,r)+
     &        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q2,r))
      endif
      if (fdwtilde .ne. 0) then
         wwa = wwa -fdwtilde*2*gw**3*sinw*(-epscrrmu(j(0,-1),q2,r,nu)*
     &        (r(mu)-q2(mu))
     &        -epscrrmu(j(0,-1),q1,r,mu)*(q1(nu)-r(nu))-
     &        epsrrmunu(q1,q2,mu,nu)*dotrc(q1,j(0,-1))
     &        +epsrrmunu(q1,q2,mu,nu)*dotrc(q2,j(0,-1))+
     &        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,q2)
     &        +epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q1,q2)+
     &        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,r)+
     &        epscrmunu(j(0,-1),r,mu,nu)
     &        *dotrr(q1,r)+epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q2,r)+
     &        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q2,r))
      endif      
      end      
      

      subroutine wwa_anomalww(qwp,qwm,qa,jp,jm,mu,wwa)
calculates anomalous WWA couplings when two currents are attached to both W"s
c qwp : momentum of the W+
c qwm: momentum of the W-
c qa : momentum of the photon
c jp: current which is attached to the W+
c jm: current which is attached to the W-
c mu: Lorentz index
c wwa: component of the leptonic tensor

      IMPLICIT NONE
      double precision qwp(0:3), qwm(0:3),qa(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu
      double complex jp(0:3,-1:1),jm(0:3,-1:1),wwa, dotrc
      double complex epscrmunu,epscrrmu,epscrrr
      integer mu
      external dotrr, dotrc,epscrmunu,epsrrmunu,epsrrrmu,epscrrmu,epscrrr

      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"      

      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      wwa=0.d0
      
      end
      

      
c ------------------------------------------------------------------------------
c WWZ VERTICES
c ------------------------------------------------------------------------------
      
       subroutine wwz_anomal(q1,q2,qz,wwz)
c Calculates anomalous WWA couplings
c  q1: momentum of W+, q1(mu)
c  q2: momentum of W-,q2(nu)
c  qz: momentum of the Z,qz(rho)
c  wwz: leptonic tensor      
      implicit none
      double precision q1(0:3),q2(0:3),qz(0:3),dotrr,sinw,cosw,epsrrrmu,epsrrmunu,epsrmunurho
      integer mu,nu,rho
      complex*16 wwz(0:3,0:3,0:3),im
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"  

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)

      do mu=0,3
         do nu=0,3
            do rho=0,3
               wwz(mu,nu,rho)=(0.d0,0.d0)
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               if (fbw .ne. 0) then  
                  wwz(mu,nu,rho)=wwz(mu,nu,rho) +fbw*gw*wmass*zmass*
     -                 sin2w*(qz(nu)*gmunu(mu,rho)-qz(mu)*gmunu(nu,rho))
               endif
               if (fwww .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)+3./2.*gw**3*fwww*cosw*
     -                 (qz(nu)*q1(rho)
     &                 *q2(mu)-qz(mu)*q1(nu)*q2(rho)+dotrr(qz,q2)*
     -                 (q1(nu)*gmunu(mu,rho)-q1(rho)*gmunu(mu,nu))
     &                 +dotrr(q2,q1)*(qz(mu)*gmunu(nu,rho)-qz(nu)*
     -                 gmunu(mu,rho))+dotrr(q1,qz)*(q2(rho)*
     -                 gmunu(mu,nu)-q2(mu)*gmunu(nu,rho)))
               endif
               if (fw .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)+fw*gw/2.*wmass*zmass*
     -                 ((q2(mu)-qz(mu))
     &                 *gmunu(nu,rho)+(qz(nu)-q1(nu))*gmunu(rho,mu)+
     -                 (q1(rho)-q2(rho))*gmunu(mu,nu)-sin2w*(qz(nu)*
     &                 gmunu(mu,rho)-qz(mu)*gmunu(nu,rho)))
               endif
               if (fb .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)-fb*gw/2.*wmass*zmass*
     &                 sin2w*(qz(nu)*gmunu(mu,rho)-qz(mu)*gmunu(nu,rho))
               endif
               if (fwwtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)-fwwtilde*2.*gw*cosw*
     &                 wmass**2*
     &                 (+epsrmunurho(q2,mu,nu,rho)+
     &                 epsrmunurho(q1,mu,nu,rho)+
     &                 epsrmunurho(qz,mu,nu,rho))
               endif
               if (fbtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)+fbtilde*wmass**2*sin2w/
     &                 cosw*gw/2.*epsrmunurho(qz,mu,nu,rho)
               endif
               if (fbwtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)-fbwtilde*wmass**2*sin2w/
     &                 cosw*gw*epsrmunurho(qz,mu,nu,rho)
               endif
               if (fwtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)+fwtilde*gw*cosw*
     &                 wmass**2*(epsrmunurho(q2,mu,nu,rho)
     &                 +epsrmunurho(q1,mu,nu,rho)+(1.+1./2.*sin2w/
     &                 cosw**2)*epsrmunurho(qz,mu,nu,rho))
               endif
               if (fwwwtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)-fwwwtilde/2.*gw**3*cosw*
     &                 (-epsrrrmu(q1,q2,qz,rho)
     &                 *gmunu(mu,nu)-epsrrrmu(q1,q2,qz,nu)*
     &                 gmunu(mu,rho)-epsrrmunu(q1,qz,nu,rho)*q2(mu)
     &                 -epsrrmunu(q1,q2,nu,rho)*qz(mu)-
     &                 epsrrrmu(q1,q2,qz,mu)*gmunu(nu,rho)+
     &                 epsrrmunu(q2,qz,mu,rho)*q1(nu)
     &                 -epsrrmunu(q1,q2,mu,rho)*qz(nu)+
     &                 epsrrmunu(q1,qz,mu,nu)*q2(rho)
     &                 +epsrrmunu(q2,qz,mu,nu)*q1(rho)-
     &                 epsrmunurho(qz,mu,nu,rho)*dotrr(q1,q2)
     &                 -epsrmunurho(q1,mu,nu,rho)*dotrr(q2,qz)-
     &                 epsrmunurho(q2,mu,nu,rho)*dotrr(q1,qz))
               endif 
               if (fdwtilde .ne. 0) then
                  wwz(mu,nu,rho)=wwz(mu,nu,rho)+fdwtilde*2.*cosw*gw**3*
     &                 (epsrrmunu(q2,qz,nu,rho)*(q2(mu)-qz(mu))
     &                 +epsrrmunu(q1,qz,mu,rho)*(qz(nu)-q1(nu))+
     &                 epsrrmunu(q2,q1,mu,nu)*(q2(rho)-q1(rho))
     &                 +epsrmunurho(q2,mu,nu,rho)*(dotrr(q2,q1)+
     &                 dotrr(q2,qz))+epsrmunurho(q1,mu,nu,rho)
     &                 *(dotrr(q2,q1)+dotrr(q1,qz))+
     &                 epsrmunurho(qz,mu,nu,rho)*(dotrr(q2,qz)+
     &                 dotrr(q1,qz)))
               endif   
            enddo
         enddo
      enddo
      end
      

      subroutine wwz_anomalwmin(q1,q2,j,r,mu,nu,wwz)
c Calculates anomalous WWZ couplings when a current j is attached to a W-
c   q1: Momentum of the W+, mu
c   q2: Momentum of the Z, nu
c   j: leptonic current
c   r: -q1-q2
c   mu,nu: Lorentz Indices of the leptonic tensor
c   wwz: component of the leptonic tensor      
      implicit none
      double precision q1(0:3), q2(0:3),r(0:3),dotrr,sinw,cosw,epsrrrmu,epsrrmunu
      double complex j(0:3,-1:1),wwz, dotrc,epscrmunu,epscrrmu,epscrrr,im
      integer mu, nu
      external dotrr, dotrc,epscrmunu,epsrrmunu,epsrrrmu,epscrrmu,epscrrr
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"               

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      wwz=(0.d0,0.d0)
      if (fbw .ne. 0) then
         wwz = wwz + fbw*gw*wmass*zmass*sinw**2*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fwww .ne. 0) then 
         wwz = wwz+3./2.*gw**3*cosw*fwww*(dotrc(q2,j(0,-1))*q1(nu)*r(mu)
     &        -q2(mu)*dotrc(q1,j(0,-1))*r(nu)+dotrr(q2,r)*(gmunu(mu,nu)*
     -        dotrc(q1,j(0,-1))
     &        -q1(nu)*j(mu,-1))+dotrr(r,q1)*(q2(mu)*j(nu,-1)-
     -        gmunu(mu,nu)*dotrc(q2,j(0,-1)))
     &        +dotrr(q1,q2)*(r(nu)*j(mu,-1)-r(mu)*j(nu,-1)))
      endif
      if(fw .ne. 0) then
         wwz = wwz +fw*gw*wmass*zmass/2.*((-q2(mu)+r(mu))*j(nu,-1)+
     -        (-dotrc(q1,j(0,-1))
     &        +dotrc(q2,j(0,-1)))*gmunu(mu,nu)+(-r(nu)+q1(nu))*j(mu,-1)
     &        -sin2w*(dotrc(q2,j(0,-1))*gmunu(mu,nu)-q2(mu)*j(nu,-1)))
      endif
      if (fb .ne. 0) then
         wwz = wwz -fb*gw*wmass*zmass*sin2w/2.*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fwwtilde .ne. 0) then
         wwz = wwz -2.*fwwtilde*gw*cosw*wmass**2*
     -        (epscrmunu(j(0,-1),r,mu,nu)
     &        +epscrmunu(j(0,-1),q1,mu,nu)+epscrmunu(j(0,-1),q2,mu,nu))
      endif
      if (fbtilde .ne. 0) then
         wwz = wwz + fbtilde*wmass**2*sin2w/cosw*gw/2.*
     -        epscrmunu(j(0,-1),q2,mu,nu)
      endif      
      if (fbwtilde .ne. 0) then
         wwz = wwz-fbwtilde*wmass**2*sin2w/cosw*gw*
     -        epscrmunu(j(0,-1),q2,mu,nu)
      endif
      if (fwtilde .ne. 0) then
         wwz = wwz +fwtilde*gw*cosw*wmass**2*(epscrmunu(j(0,-1),r,mu,nu)
     &        +epscrmunu(j(0,-1),q1,mu,nu)+(1.+1./2.*sin2w/cosw**2)*
     -        epscrmunu(j(0,-1),q2,mu,nu))
      endif
      if (fwwwtilde .ne. 0) then
         wwz = wwz +fwwwtilde/2.*cosw*gw**3*
     &        (-epsrrrmu(r,q1,q2,nu)*j(mu,-1)+epscrrr(j(0,-1),r,q1,q2)*
     -        gmunu(mu,nu)
     &        +epscrrmu(j(0,-1),q1,q2,nu)*r(mu)-
     -        epscrrmu(j(0,-1),r,q1,nu)*q2(mu)-epsrrrmu(r,q1,q2,mu)*j(nu,-1)
     &        -epsrrmunu(r,q2,mu,nu)*dotrc(q1,j(0,-1))-
     -        epsrrmunu(r,q1,mu,nu)*dotrc(q2,j(0,-1))
     &        +epscrrmu(j(0,-1),q1,q2,mu)*r(nu)+
     -        epscrrmu(j(0,-1),r,q2,mu)*q1(nu)+
     -        epscrmunu(j(0,-1),q2,mu,nu)*
     &        dotrr(r,q1)+epscrmunu(j(0,-1),q1,mu,nu)*dotrr(r,q2)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q1,q2))
      endif
      if (fdwtilde .ne. 0) then
         wwz = wwz +fdwtilde*2*gw**3*cosw*(-epscrrmu(j(0,-1),q2,r,nu)*
     -        (r(mu)-q2(mu))
     &        -epscrrmu(j(0,-1),q1,r,mu)*(q1(nu)-r(nu))-
     -        epsrrmunu(q1,q2,mu,nu)*dotrc(q1,j(0,-1))
     &        +epsrrmunu(q1,q2,mu,nu)*dotrc(q2,j(0,-1))+
     -        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,q2)
     &        +epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q1,q2)+
     -        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)
     &        *dotrr(q1,r)+epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q2,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q2,r))
      endif                     
      end   

      
      subroutine wwz_anomalwplus(q1,q2,j,r,mu,nu,wwz)
c Calculates anomalous WWZ couplings when a current j is attached to a W+
c   q1: Momentum of the W-, mu
c   q2: Momentum of the Z, nu
c   j: leptonic current
c   r: -q1-q2
c   mu,nu: Lorentz Indices of the leptonic tensor
c   wwz: component of the leptonic tensor      
      implicit none
      double precision q1(0:3), q2(0:3),r(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu
      double complex j(0:3,-1:1),wwz,dotrc,epscrmunu,epscrrmu,epscrrr,im
      integer mu, nu
      external dotrr,dotrc,epscrmunu,epsrrmunu,epsrrrmu,epscrrmu,epscrrr
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"            

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      wwz=(0.d0,0.d0)
      if (fbw .ne. 0) then
         wwz = wwz - fbw*gw*wmass*zmass*sinw**2*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fwww .ne. 0) then 
         wwz = wwz-3./2.*gw**3*cosw*fwww*(dotrc(q2,j(0,-1))*q1(nu)*r(mu)
     &        -q2(mu)*dotrc(q1,j(0,-1))*r(nu)+dotrr(q2,r)*
     -        (gmunu(mu,nu)*dotrc(q1,j(0,-1))
     &        -q1(nu)*j(mu,-1))+dotrr(r,q1)*(q2(mu)*j(nu,-1)-gmunu(mu,nu)
     &        *dotrc(q2,j(0,-1)))+dotrr(q1,q2)*(r(nu)*j(mu,-1)-r(mu)*
     -        j(nu,-1)))
      endif
      if(fw .ne. 0) then
         wwz = wwz -fw*gw*wmass*zmass/2.*((-q2(mu)+r(mu))*j(nu,-1)+
     -        (-dotrc(q1,j(0,-1))
     &        +dotrc(q2,j(0,-1)))*gmunu(mu,nu)+(-r(nu)+q1(nu))*j(mu,-1)
     &        -sin2w*(dotrc(q2,j(0,-1))*gmunu(mu,nu)-q2(mu)*j(nu,-1)))
      endif
      if (fb .ne. 0) then
         wwz = wwz +fb*gw*wmass*zmass*sin2w/2.*(dotrc(q2,j(0,-1))*
     -        gmunu(mu,nu)-q2(mu)*j(nu,-1))
      endif
      if (fwwtilde .ne. 0) then
         wwz = wwz +2.*fwwtilde*gw*cosw*wmass**2*
     -        (epscrmunu(j(0,-1),r,mu,nu)
     &        +epscrmunu(j(0,-1),q1,mu,nu)+epscrmunu(j(0,-1),q2,mu,nu))
      endif
      if (fbtilde .ne. 0) then
         wwz = wwz - fbtilde*wmass**2*sin2w/cosw*gw/2.*
     -        epscrmunu(j(0,-1),q2,mu,nu)
      endif      
      if (fbwtilde .ne. 0) then
         wwz = wwz+fbwtilde*wmass**2*sin2w/cosw*gw*
     -        epscrmunu(j(0,-1),q2,mu,nu)
      endif
      if (fwtilde .ne. 0) then
         wwz = wwz -fwtilde*gw*cosw*wmass**2*(epscrmunu(j(0,-1),r,mu,nu)
     &        +epscrmunu(j(0,-1),q1,mu,nu)+(1.+1./2.*sin2w/cosw**2)*
     -        epscrmunu(j(0,-1),q2,mu,nu))
      endif
      if (fwwwtilde .ne. 0) then
         wwz = wwz -fwwwtilde/2.*cosw*gw**3*
     &        (-epsrrrmu(r,q1,q2,nu)*j(mu,-1)+epscrrr(j(0,-1),r,q1,q2)*
     -        gmunu(mu,nu)
     &        +epscrrmu(j(0,-1),q1,q2,nu)*r(mu)-
     -        epscrrmu(j(0,-1),r,q1,nu)*q2(mu)-epsrrrmu(r,q1,q2,mu)*
     -        j(nu,-1)
     &        -epsrrmunu(r,q2,mu,nu)*dotrc(q1,j(0,-1))-
     -        epsrrmunu(r,q1,mu,nu)*dotrc(q2,j(0,-1))
     &        +epscrrmu(j(0,-1),q1,q2,mu)*r(nu)+
     -        epscrrmu(j(0,-1),r,q2,mu)*q1(nu)+
     -        epscrmunu(j(0,-1),q2,mu,nu)*
     &        dotrr(r,q1)+epscrmunu(j(0,-1),q1,mu,nu)*dotrr(r,q2)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q1,q2))
      endif 
      if (fdwtilde .ne. 0) then
         wwz = wwz -fdwtilde*2*gw**3*cosw*(-epscrrmu(j(0,-1),q2,r,nu)*
     -        (r(mu)-q2(mu))
     &        -epscrrmu(j(0,-1),q1,r,mu)*(q1(nu)-r(nu))-
     -        epsrrmunu(q1,q2,mu,nu)*dotrc(q1,j(0,-1))
     &        +epsrrmunu(q1,q2,mu,nu)*dotrc(q2,j(0,-1))+
     -        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,q2)
     &        +epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q1,q2)+
     -        epscrmunu(j(0,-1),q1,mu,nu)*dotrr(q1,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)
     &        *dotrr(q1,r)+epscrmunu(j(0,-1),q2,mu,nu)*dotrr(q2,r)+
     -        epscrmunu(j(0,-1),r,mu,nu)*dotrr(q2,r))
      endif                    
      end   
      
                             
c -----------------------------------------------------------------------
c     WWH VERTICES
c -----------------------------------------------------------------------
      
      subroutine wwh_anomal(q1,q2,qh,wwh)
c calculates anomalous WWH couplings
c q1: momentum of  W+, q1(mu)
c q2: momentum of the W-, q2(nu)
c qh: momentum of the higgs boson
c wwh: leptonic tensor
      implicit none
      double precision q1(0:3),q2(0:3),qh(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu
      integer mu,nu
      complex*16 wwh(0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc" 

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            wwh(mu,nu)=(0.d0,0.d0)
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            if (fww .ne. 0) then
               wwh(mu,nu)=wwh(mu,nu)-fww*2.*gw*wmass*(q2(mu)*q1(nu)-
     -              dotrr(q1,q2)*gmunu(mu,nu))
            endif
            if (fw .ne. 0) then
               wwh(mu,nu)=wwh(mu,nu)+fw/2.*gw*wmass*((dotrr(q1,qh)+
     -              dotrr(q2,qh))*gmunu(mu,nu)-qh(mu)*q1(nu)-qh(nu)*
     -              q2(mu))
            endif
            if (fwwtilde .ne. 0) then
               wwh(mu,nu) = wwh(mu,nu)+fwwtilde*2.*gw*wmass*
     -              epsrrmunu(q2,q1,mu,nu)
            endif
            if (fwtilde .ne. 0) then
               wwh(mu,nu) = wwh(mu,nu)-fwtilde*gw*wmass*
     -              epsrrmunu(q2,q1,mu,nu)
            endif
         enddo
      enddo
      end
      
      

c ------------------------------------------------------------------
c      ZZH-Vertices
c ------------------------------------------------------------------

      subroutine zzh_anomal(q1,q2,qh,zzh)
c calculates anomalous zzh couplings
c q1,q2: momenta of the Z"s: q1(mu),q2(nu)
c qh: momentum of th Higgs boson
c zzh: leptonic tensor
      IMPLICIT NONE
      double precision q1(0:3),q2(0:3),qh(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu 
      integer mu,nu
      complex*16  zzh(0:3,0:3),im      
      external dotrr,epsrrrmu,epsrrmunu   
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc" 

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            zzh(mu,nu)=(0.d0,0.d0)
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            if (fbw .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+2.*fbw*wmass*gw*sin2w*
     -              (dotrr(q1,q2)*gmunu(mu,nu)-q1(nu)*q2(mu))
            endif
            if (fww .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu) +fww*wmass*gw*cosw**2*2.*
     -              (dotrr(q1,q2)*gmunu(mu,nu)-q1(nu)*q2(mu))
            endif
            if (fbb .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu) +fbb*sin2w**2/cosw*gw*zmass*2.*
     -              (dotrr(q1,q2)*gmunu(mu,nu)-q1(nu)*q2(mu))
            endif
            if (fw .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+fw*gw*wmass/2.*((dotrr(q1,qh)+
     -              dotrr(q2,qh))*gmunu(mu,nu)-qh(mu)*q1(nu)-qh(nu)*q2(mu))
            endif
            if (fb .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+fb*sin2w/cosw*gw*zmass/2.*
     -              ((dotrr(q1,qh)+dotrr(q2,qh))*gmunu(mu,nu)-qh(mu)*
     -              q1(nu)-qh(nu)*q2(mu))
            endif
            if (fwwtilde .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)-fwwtilde*2.*gw*wmass*cosw**2*
     -              epsrrmunu(q1,q2,mu,nu)
            endif
            if (fbtilde .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+fbtilde*sin2w/cosw*gw*zmass/2.*
     -              (epsrrmunu(qh,q1,mu,nu)-epsrrmunu(qh,q2,mu,nu))
            endif
            if (fbbtilde .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)-fbbtilde*2.*sin2w**2/cosw*gw*zmass*
     -              epsrrmunu(q1,q2,mu,nu)
            endif
            if (fbwtilde .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+fbwtilde*(-sin2w/cosw*gw*zmass*
     -              (epsrrmunu(qh,q1,mu,nu)
     &              -epsrrmunu(qh,q2,mu,nu))+2.*sin2w**2/cosw*gw*zmass*
     -              epsrrmunu(q1,q2,mu,nu))
            endif
            if (fwtilde .ne. 0) then
               zzh(mu,nu)=zzh(mu,nu)+fwtilde*(sin2w/cosw*gw*zmass/2.*
     -              (epsrrmunu(qh,q1,mu,nu)
     &              -epsrrmunu(qh,q2,mu,nu))+gw*wmass*cosw**2*
     -              epsrrmunu(q1,q2,mu,nu)
     &              -sin2w**2/cosw*gw*zmass*epsrrmunu(q1,q2,mu,nu))
            endif      
         enddo
      enddo
      end
      
c ------------------------------------------------------------------
c     AZH-Vertices
c ------------------------------------------------------------------

      subroutine azh_anomal(qa,qz,qh,azh)
c calculates anomalous AZH couplings
c qa: momentum of the photon, qa(mu)
c qz: momentum of the Z, qz(nu)
c qh: momentum of the Higgs boson
c azh: azh tensor of rank 2            
      IMPLICIT NONE
      double precision qa(0:3),qz(0:3),qh(0:3),dotrr,sinw,cosw,epsrrrmu,epsrrmunu 
      integer mu,nu
      complex*16 azh(0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu   
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"  

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            azh(mu,nu)=(0.d0,0.d0)
      enddo
      enddo
      do mu=0,3
         do nu=0,3
            if (fbw .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fbw*gw*sinw/cosw*wmass*(cosw**2-
     -              sin2w)*(-gmunu(mu,nu)*dotrr(qa,qz)+qz(mu)*qa(nu))
            endif
            if (fww .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fww*2.*gw*wmass*cosw*sinw*
     -              (gmunu(mu,nu)*dotrr(qa,qz)-qz(mu)*qa(nu))
            endif
            if (fbb .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fbb*2.*gw*sinw/cosw*wmass*sin2w*
     -              (-gmunu(mu,nu)*dotrr(qa,qz)+qz(mu)*qa(nu))
            endif
            if (fw .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fw*gw/2.*wmass*sinw/cosw*
     -              (gmunu(mu,nu)*dotrr(qa,qh)-qh(mu)*qa(nu))
            endif
            if (fb .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fb*gw/2.*wmass*sinw/cosw*(-gmunu(mu,nu)*
     -              dotrr(qa,qh)+qh(mu)*qa(nu))
            endif
            if (fwwtilde .ne. 0) then
               azh(mu,nu)=azh(mu,nu)-fwwtilde*2.*gw*wmass*cosw*sinw*
     -              epsrrmunu(qa,qz,mu,nu)
            endif
            if (fbtilde .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fbtilde/2.*zmass*gw*sinw*
     -              epsrrmunu(qa,qh,mu,nu)
            endif
            if (fbbtilde .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fbbtilde*2.*zmass*gw*sinw**3*
     -              epsrrmunu(qa,qz,mu,nu)
            endif
            if (fbwtilde .ne. 0) then
               azh(mu,nu)=azh(mu,nu)+fbwtilde*(-zmass*gw*sinw*
     -              epsrrmunu(qa,qh,mu,nu)-2.*zmass*gw*sinw**3*
     -              epsrrmunu(qa,qz,mu,nu))
            endif      
            if (fwtilde .ne. 0) then
               azh(mu,nu)=azh(mu,nu)-fwtilde/2.*gw*wmass*sinw/cosw*epsrrmunu(qa,qh,mu,nu)
            endif            
         enddo
      enddo
      end
      
c -----------------------------------------------------------------
c          AAH - Vertices
c -----------------------------------------------------------------

      subroutine aah_anomal(qa1,qa2,qh,aah)
c calculates anomalous AAH couplings
c qa1,qa2: momenta of the photons: qa1(mu), qa2(nu)
c qH: momentum of the Higgs boson
c aah: leptonic tensor of rank 2  
      IMPLICIT NONE
      double precision qa1(0:3),qa2(0:3),qh(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu 
      integer mu,nu
      complex*16 aah(0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu   
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc" 

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            aah(mu,nu)=(0.d0,0.d0)
         enddo
      enddo
      do mu=0,3
         do nu=0,3          
            if (fbw .ne. 0) then
               aah(mu,nu)=aah(mu,nu)+fbw*gw*wmass*sin2w*2.*
     -              (-gmunu(mu,nu)*dotrr(qa1,qa2)+qa1(nu)*qa2(mu))
            endif
            if (fww .ne. 0) then
               aah(mu,nu)=aah(mu,nu)+fww*gw*wmass*sin2w*2.*
     -              (gmunu(mu,nu)*dotrr(qa1,qa2)-qa1(nu)*qa2(mu))
            endif
            if (fbb .ne. 0) then
               aah(mu,nu)=aah(mu,nu)+fbb*gw*wmass*sin2w*2.*
     -              (gmunu(mu,nu)*dotrr(qa1,qa2)-qa1(nu)*qa2(mu))
            endif
            if (fwwtilde .ne. 0) then
               aah(mu,nu)=aah(mu,nu)-fwwtilde*2.*gw*wmass*sin2w*
     -              epsrrmunu(qa1,qa2,mu,nu)
            endif
            if (fbbtilde .ne. 0) then
               aah(mu,nu)=aah(mu,nu)-fbbtilde*2.*zmass*gw*sin2w*cosw*
     -              epsrrmunu(qa1,qa2,mu,nu)
            endif
            if (fbwtilde .ne. 0) then
               aah(mu,nu)=aah(mu,nu)+fbwtilde*2.*zmass*gw*sin2w*cosw*
     -              epsrrmunu(qa1,qa2,mu,nu)
            endif
            if (fwtilde .ne. 0) then
               aah(mu,nu)=aah(mu,nu)+fwtilde*(gw*wmass*sin2w*
     -              epsrrmunu(qa1,qa2,mu,nu)-zmass*gw*sin2w*cosw*
     -              epsrrmunu(qa1,qa2,mu,nu))
            endif
         enddo
      enddo
      end
     
c -----------------------------------------------------------------
c          WWWW - Vertices
c ------------------------------------------------------------------

      subroutine wwww_anomal(q1,q2,q3,q4,wwww)
c Calculates anomalous wwww couplings
c  q1-q4 momenta of the W"s: q1(mu) W+ ,q2(nu) W- ,q3(rho) W+ ,q4(sigma) W-
c  wwww: tensor of the 4 W vertex 
      IMPLICIT NONE     
      double precision q1(0:3),q2(0:3),q3(0:3),q4(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      integer mu,nu,rho,sigma
      complex*16 wwww(0:3,0:3,0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc" 

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  wwww(mu,nu,rho,sigma)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  if (fwww .ne. 0) then
                     wwww(mu,nu,rho,sigma)=wwww(mu,nu,rho,sigma)-fwww*
     &                    3./2.*gw**4
     &                    *((q4(mu)*(q1(nu)-q3(nu))+q2(mu)*q3(nu)+
     &                    q3(mu)*q4(nu))*gmunu(rho,sigma)
     &                    -gmunu(mu,nu)*(dotrr(q1,q4)+dotrr(q2,q3))*
     &                    gmunu(rho,sigma)
     &                    +gmunu(mu,sigma)*(q4(nu)*q1(rho)+q1(nu)*
     &                    (q2(rho)-q4(rho))+q3(nu)*q4(rho))
     &                    -gmunu(nu,sigma)*((q2(mu)+q4(mu))*q1(rho)+
     &                    q3(mu)*(q2(rho)+q4(rho)))
     &                    +gmunu(nu,rho)*(q3(mu)*q2(sigma)+q2(mu)*
     &                    (q1(sigma)-q3(sigma))+q4(mu)*q3(sigma))
     &                    +gmunu(mu,nu)*(q4(rho)*q1(sigma)+q1(rho)*
     &                    q2(sigma)+q2(rho)*(q3(sigma)-q1(sigma)))
     &                    -gmunu(mu,rho)*((q1(nu)+q3(nu))*q2(sigma)+
     &                    q4(nu)*(q1(sigma)+q3(sigma)))
     &                    -gmunu(mu,sigma)*gmunu(nu,rho)*(dotrr(q1,q2)+
     &                    dotrr(q3,q4))
     &                    +gmunu(mu,rho)*gmunu(nu,sigma)*(dotrr(q1,q2)+
     &                    dotrr(q1,q4)+dotrr(q2,q3)+dotrr(q3,q4)))
                  endif
                  if (fw .ne. 0) then
                     wwww(mu,nu,rho,sigma)=wwww(mu,nu,rho,sigma)+fw*
     &                    gw**2*wmass**2*(2.*gmunu(mu,rho)*
     &                    gmunu(nu,sigma)
     &                    -gmunu(mu,sigma)*gmunu(nu,rho)-gmunu(mu,nu)*
     &                    gmunu(sigma,rho))
                  endif
                  if (fwwwtilde .ne. 0) then
                     wwww(mu,nu,rho,sigma)=wwww(mu,nu,rho,sigma)-
     &                    fwwwtilde/2.*gw**4*((epsrrmunu(q1,q2,rho,sigma)
     &                    +epsrrmunu(q1,q4,rho,sigma)-
     &                    epsrrmunu(q2,q3,rho,sigma))*gmunu(mu,nu)
     &                    +epsrrmunu(q1,q2,nu,sigma)*gmunu(mu,rho)-
     &                    (epsrrmunu(q1,q4,nu,sigma)
     &                    +epsrrmunu(q2,q3,nu,sigma)+
     &                    epsrrmunu(q3,q4,nu,sigma))*gmunu(mu,rho)
     &                    -(epsrrmunu(q1,q2,nu,rho)+
     &                    epsrrmunu(q1,q4,nu,rho)+
     &                    epsrrmunu(q3,q4,nu,rho))*gmunu(mu,sigma)
     &                    +(-epsrmunurho(q2,nu,rho,sigma)+
     &                    epsrmunurho(q4,nu,rho,sigma))*q3(mu)
     &                    +(-epsrmunurho(q1,nu,rho,sigma)-
     &                    epsrmunurho(q3,nu,rho,sigma))*(q2(mu)-q4(mu))
     &                    +(epsrrmunu(q1,q2,mu,sigma)-
     &                    epsrrmunu(q2,q3,mu,sigma)+
     &                    epsrrmunu(q3,q4,mu,sigma))*gmunu(nu,rho)
     &                    -(epsrrmunu(q1,q2,mu,rho)+
     &                    epsrrmunu(q1,q4,mu,rho)+
     &                    epsrrmunu(q2,q3,mu,rho)
     &                    -epsrrmunu(q3,q4,mu,rho))*gmunu(nu,sigma)-
     &                    (-epsrmunurho(q2,mu,rho,sigma)-
     &                    epsrmunurho(q4,mu,rho,sigma))*(q1(nu)-q3(nu))
     &                    +(-epsrmunurho(q3,mu,rho,sigma)+
     &                    epsrmunurho(q1,mu,rho,sigma))*q4(nu)
     &                    +(epsrrmunu(q1,q4,mu,nu)-
     &                    epsrrmunu(q2,q3,mu,nu)+
     &                    epsrrmunu(q3,q4,mu,nu))*gmunu(rho,sigma)
     &                    +(-epsrmunurho(q4,mu,nu,sigma)+
     &                    epsrmunurho(q2,mu,nu,sigma))*q1(rho)
     &                    -(-epsrmunurho(q1,mu,nu,sigma)-
     &                    epsrmunurho(q3,mu,nu,sigma))*(q2(rho)-q4(rho))
     &                    +(-epsrmunurho(q1,mu,nu,rho)+
     &                    epsrmunurho(q3,mu,nu,rho))*q2(sigma)
     &                    +(-epsrmunurho(q2,mu,nu,rho)-
     &                    epsrmunurho(q4,mu,nu,rho))*(q1(sigma)-
     &                    q3(sigma))+epsmunurhosigma(mu,nu,rho,sigma)
     &                    *(dotrr(q1,q2)-dotrr(q1,q4)-dotrr(q2,q3)+
     &                    dotrr(q3,q4)))
                  endif
                  if (fdwtilde .ne. 0) then
                     wwww(mu,nu,rho,sigma)=wwww(mu,nu,rho,sigma)-
     &                    fdwtilde*2.*gw**4
     &                    *(epsrrmunu(q3,q4,rho,sigma)*gmunu(mu,nu)-2.*
     &                    epsrrmunu(q2,q4,nu,sigma)*gmunu(mu,rho)
     &                    +epsrrmunu(q2,q3,nu,rho)*gmunu(mu,sigma)-
     &                    epsrmunurho(q2,nu,rho,sigma)*(q2(mu)-q3(mu)-
     &                    q4(mu))
     &                    -epsrmunurho(q4,nu,rho,sigma)*(q2(mu)+q3(mu)-
     &                    q4(mu))+epsrrmunu(q1,q4,mu,sigma)*gmunu(nu,rho)
     &                    -2.*epsrrmunu(q1,q3,mu,rho)*gmunu(nu,sigma)+
     &                    epsrmunurho(q3,mu,rho,sigma)*(q1(nu)-q3(nu)+
     &                    q4(nu))
     &                    -epsrmunurho(q1,mu,rho,sigma)*(-q1(nu)+q3(nu)+
     &                    q4(nu))+epsrrmunu(q1,q2,mu,nu)
     &                    *gmunu(rho,sigma)+epsrmunurho(q4,mu,nu,sigma)*
     &                    (q1(rho)+q2(rho)-q4(rho))
     &                    -epsrmunurho(q2,mu,nu,sigma)*(q1(rho)-q2(rho)+
     &                    q4(rho))-epsrmunurho(q1,mu,nu,rho)*(q1(sigma)-
     &                    q2(sigma)-q3(sigma))
     &                    -epsrmunurho(q3,mu,nu,rho)*(q1(sigma)+
     &                    q2(sigma)-q3(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (-dotrr(q1,q2)+dotrr(q1,q4)+dotrr(q2,q3)-
     &                    dotrr(q3,q4)))
                  endif 
               enddo
            enddo
         enddo
      enddo     
      end
      

c -----------------------------------------------------------
c     WWZZ-Vertices
c -----------------------------------------------------------

      subroutine wwzz_anomal(qwp,qwm,qz1,qz2,wwzz)
c calculates anomalous wwzz couplings
c qwp: momentum of the W+, qwp(rho)
c qwm: momentum of the W-, qwm(sigma)
c qz1,qz2: momenta of the Z"s: qz1(mu),qz2(nu)
c wwzz: wwzz-tensor
      IMPLICIT NONE

      double precision qwp(0:3),qwm(0:3),qz1(0:3),qz2(0:3)
      double precision dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      integer mu,nu,rho,sigma
      complex*16 wwzz(0:3,0:3,0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc" 

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  wwzz(mu,nu,rho,sigma)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  if (fwww .ne. 0) then
                     wwzz(mu,nu,rho,sigma)=wwzz(mu,nu,rho,sigma)-
     &                    fwww*3./2.*cosw**2*gw**4
     &                    *((qz2(mu)*(qwm(nu)+qwp(nu))+(qwm(mu)+qwp(mu))
     &                    *qz1(nu))*gmunu(rho,sigma)
     &                    -gmunu(mu,nu)*(dotrr(qwm,qz1)+dotrr(qwm,qz2)+
     &                    dotrr(qwp,qz1)+dotrr(qwp,qz2))*gmunu(rho,sigma)
     &                    -gmunu(nu,sigma)*(qz2(mu)*qwm(rho)+qwm(mu)*
     &                    (qz1(rho)-qz2(rho))
     &                    +qwp(mu)*qz2(rho))-gmunu(mu,sigma)*(qz1(nu)*
     &                    qwm(rho)+qwp(nu)*qz1(rho)
     &                    +qwm(nu)*(qz2(rho)-qz1(rho)))-gmunu(nu,rho)*
     &                    (qz2(mu)*qwp(sigma)+qwp(mu)*(qz1(sigma)
     &                    -qz2(sigma))+qwm(mu)*qz2(sigma))-gmunu(mu,rho)
     &                    *(qz1(nu)*qwp(sigma)+qwm(nu)*qz1(sigma)
     &                    +qwp(nu)*(qz2(sigma)-qz1(sigma)))+gmunu(mu,nu)
     &                    *((qz1(rho)+qz2(rho))*qwp(sigma)
     &                    +qwm(rho)*(qz1(sigma)+qz2(sigma)))+
     &                    gmunu(mu,sigma)*gmunu(nu,rho)
     &                    *(dotrr(qwm,qz2)+dotrr(qwp,qz1))+gmunu(mu,rho)
     &                    *gmunu(nu,sigma)*(dotrr(qwm,qz1)+dotrr(qwp,qz2)))
                  endif      
                  if (fw .ne. 0) then
                     wwzz(mu,nu,rho,sigma)=wwzz(mu,nu,rho,sigma)+fw*
     &                    gw**2*wmass**2*(-2.*gmunu(mu,nu)*gmunu(rho,sigma)
     &                    +gmunu(mu,rho)*gmunu(nu,sigma)+gmunu(mu,sigma)
     &                    *gmunu(nu,rho))
                  endif
                  if (fwwwtilde .ne. 0) then
                     wwzz(mu,nu,rho,sigma)=wwzz(mu,nu,rho,sigma)+
     &                    fwwwtilde/2.*cosw**2*gw**4
     &                    *((epsrrmunu(qwm,qz1,rho,sigma)+
     &                    epsrrmunu(qwm,qz2,rho,sigma)-
     &                    epsrrmunu(qwp,qz1,rho,sigma)
     &                    -epsrrmunu(qwp,qz2,rho,sigma))*gmunu(mu,nu)-
     &                    (epsrrmunu(qwm,qz1,nu,sigma)
     &                    +epsrrmunu(qwp,qz1,nu,sigma)+
     &                    epsrrmunu(qwp,qz2,nu,sigma))*gmunu(mu,rho)
     &                    -(epsrrmunu(qwm,qz1,nu,rho)+
     &                    epsrrmunu(qwm,qz2,nu,rho)+
     &                    epsrrmunu(qwp,qz1,nu,rho))*gmunu(mu,sigma)
     &                    +(-epsrmunurho(qz1,nu,rho,sigma)-
     &                    epsrmunurho(qz2,nu,rho,sigma))*(qwm(mu)-qwp(mu))
     &                    +(-epsrmunurho(qwm,nu,rho,sigma)+
     &                    epsrmunurho(qwp,nu,rho,sigma))*qz2(mu)
     &                    -(epsrrmunu(qwm,qz2,mu,sigma)+
     &                    epsrrmunu(qwp,qz1,mu,sigma)+
     &                    epsrrmunu(qwp,qz2,mu,sigma))*gmunu(nu,rho)
     &                    -(epsrrmunu(qwm,qz1,mu,rho)+
     &                    epsrrmunu(qwm,qz2,mu,rho)+
     &                    epsrrmunu(qwp,qz2,mu,rho))*gmunu(nu,sigma)
     &                    +(-epsrmunurho(qz1,mu,rho,sigma)-
     &                    epsrmunurho(qz2,mu,rho,sigma))*
     &                    (qwm(nu)-qwp(nu))
     &                    +(-epsrmunurho(qwm,mu,rho,sigma)+
     &                    epsrmunurho(qwp,mu,rho,sigma))*qz1(nu)
     &                    +(epsrrmunu(qwm,qz1,mu,nu)-
     &                    epsrrmunu(qwm,qz2,mu,nu)+
     &                    epsrrmunu(qwp,qz1,mu,nu)-
     &                    epsrrmunu(qwp,qz2,mu,nu))*gmunu(rho,sigma)
     &                    +(-epsrmunurho(qz2,mu,nu,sigma)+
     &                    epsrmunurho(qz1,mu,nu,sigma))*qwm(rho)
     &                    -(-epsrmunurho(qwm,mu,nu,sigma)-
     &                    epsrmunurho(qwp,mu,nu,sigma))*(qz1(rho)-qz2(rho))
     &                    +(-epsrmunurho(qz2,mu,nu,rho)+
     &                    epsrmunurho(qz1,mu,nu,rho))*qwp(sigma)
     &                    -(-epsrmunurho(qwm,mu,nu,rho)-
     &                    epsrmunurho(qwp,mu,nu,rho))*(qz1(sigma)-
     &                    qz2(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (dotrr(qwm,qz1)-dotrr(qwm,qz2)-dotrr(qwp,qz1)+
     &                    dotrr(qwp,qz2)))
                  endif
                  if (fdwtilde .ne. 0) then
                     wwzz(mu,nu,rho,sigma)=wwzz(mu,nu,rho,sigma)+
     &                    fdwtilde*2.*gw**4*cosw**2
     &                    *(2.*epsrrmunu(qwm,qwp,rho,sigma)*gmunu(mu,nu)
     &                    +epsrrmunu(qz2,qwm,nu,sigma)*gmunu(mu,rho)
     &                    +epsrrmunu(qz2,qwp,nu,rho)*gmunu(mu,sigma)-
     &                    epsrmunurho(qwp,nu,rho,sigma)*(qz2(mu)+qwm(mu)
     &                    -qwp(mu))
     &                    +epsrmunurho(qwm,nu,rho,sigma)*(qz2(mu)-
     &                    qwm(mu)+qwp(mu))+epsrrmunu(qz1,qwm,mu,sigma)
     &                    *gmunu(nu,rho)+epsrrmunu(qz1,qwp,mu,rho)*
     &                    gmunu(nu,sigma)-epsrmunurho(qwp,mu,rho,sigma)*
     &                    (qz1(nu)+qwm(nu)-qwp(nu))
     &                    +epsrmunurho(qwm,mu,rho,sigma)*(qz1(nu)-
     &                    qwm(nu)+qwp(nu))
     &                    -2.*epsrrmunu(qz1,qz2,mu,nu)*gmunu(rho,sigma)+
     &                    epsrmunurho(qz2,mu,nu,sigma)*(qz1(rho)-
     &                    qz2(rho)+qwm(rho))
     &                    -epsrmunurho(qz1,mu,nu,sigma)*(-qz1(rho)+
     &                    qz2(rho)+qwm(rho))+epsrmunurho(qz2,mu,nu,rho)
     &                    *(qz1(sigma)-qz2(sigma)+qwp(sigma))-
     &                    epsrmunurho(qz1,mu,nu,rho)*(-qz1(sigma)+
     &                    qz2(sigma)+qwp(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (-dotrr(qz1,qwm)+dotrr(qz1,qwp)+dotrr(qz2,qwm)
     &                    -dotrr(qz2,qwp)))
                  endif
               enddo
            enddo
         enddo
      enddo
      end


c -------------------------------------------------------------------------
c    WWAZ-Vertices
c -------------------------------------------------------------------------

      subroutine wwaz_anomal(qwp,qwm,qa,qz,wwaz)
c calculates anomalous WWAZ couplings
c qwp: momentum of the W+, qwp(rho)
c qwm: momentum of the W-, qwm(sigma)
c qa: momentum of the photon, qa(mu)
c qz: momentum of the Z, qz(nu)
c wwaz: wwaz-tensor of rank 4  
      IMPLICIT NONE          
      double precision qwp(0:3),qwm(0:3),qa(0:3),qz(0:3),dotrr,sinw,cosw
      double precision epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      integer mu,nu,rho,sigma
      complex*16 wwaz(0:3,0:3,0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"
  
      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  wwaz(mu,nu,rho,sigma)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  if (fwww .ne. 0) then
                     wwaz(mu,nu,rho,sigma)=wwaz(mu,nu,rho,sigma)-
     &                    fwww*3./2.*cosw*sinw*gw**4*((qz(mu)
     &                    *(qwm(nu)+qwp(nu))+(qwm(mu)+qwp(mu))*qa(nu))*
     &                    gmunu(rho,sigma)
     &                    -gmunu(mu,nu)*(dotrr(qwm,qa)+dotrr(qwm,qz)+
     &                    dotrr(qwp,qa)+dotrr(qwp,qz))*gmunu(rho,sigma)
     &                    -gmunu(nu,sigma)*(qz(mu)*qwm(rho)+qwm(mu)*
     &                    (qa(rho)-qz(rho))
     &                    +qwp(mu)*qz(rho))-gmunu(mu,sigma)*(qa(nu)*
     &                    qwm(rho)+qwp(nu)*qa(rho)+qwm(nu)*(qz(rho)-
     &                    qa(rho)))
     &                    -gmunu(nu,rho)*(qz(mu)*qwp(sigma)+qwp(mu)*
     &                    (qa(sigma)
     &                    -qz(sigma))+qwm(mu)*qz(sigma))-gmunu(mu,rho)*
     &                    (qa(nu)*qwp(sigma)+qwm(nu)*qa(sigma)+qwp(nu)*
     &                    (qz(sigma)
     &                    -qa(sigma)))+gmunu(mu,nu)*((qa(rho)+qz(rho))*
     &                    qwp(sigma)
     &                    +qwm(rho)*(qa(sigma)+qz(sigma)))+
     &                    gmunu(mu,sigma)*gmunu(nu,rho)*(dotrr(qwm,qz)+
     &                    dotrr(qwp,qa))
     &                    +gmunu(mu,rho)*gmunu(nu,sigma)*(dotrr(qwm,qa)+
     &                    dotrr(qwp,qz)))
                  endif
                  if (fw .ne. 0) then
                     wwaz(mu,nu,rho,sigma)=wwaz(mu,nu,rho,sigma)+fw/2.*
     &                    gw**2*wmass**2
     &                    *sinw/cosw*(-2.*gmunu(mu,nu)*gmunu(rho,sigma)+
     &                    gmunu(mu,rho)*gmunu(nu,sigma)
     &                    +gmunu(mu,sigma)*gmunu(nu,rho))
                  endif 
                  if (fwwwtilde .ne. 0) then
                     wwaz(mu,nu,rho,sigma)=wwaz(mu,nu,rho,sigma)+
     &                    fwwwtilde/2.*cosw*sinw*gw**4*
     &                    ((epsrrmunu(qwm,qa,rho,sigma)+
     &                    epsrrmunu(qwm,qz,rho,sigma)-
     &                    epsrrmunu(qwp,qa,rho,sigma)
     &                    -epsrrmunu(qwp,qz,rho,sigma))*gmunu(mu,nu)-
     &                    (epsrrmunu(qwm,qa,nu,sigma)
     &                    +epsrrmunu(qwp,qa,nu,sigma)+
     &                    epsrrmunu(qwp,qz,nu,sigma))*gmunu(mu,rho)
     &                    -(epsrrmunu(qwm,qa,nu,rho)+
     &                    epsrrmunu(qwm,qz,nu,rho)+
     &                    epsrrmunu(qwp,qa,nu,rho))*gmunu(mu,sigma)
     &                    +(-epsrmunurho(qa,nu,rho,sigma)-
     &                    epsrmunurho(qz,nu,rho,sigma))*(qwm(mu)-qwp(mu))
     &                    +(-epsrmunurho(qwm,nu,rho,sigma)+
     &                    epsrmunurho(qwp,nu,rho,sigma))*qz(mu)-
     &                    (epsrrmunu(qwm,qz,mu,sigma)
     &                    +epsrrmunu(qwp,qa,mu,sigma)+
     &                    epsrrmunu(qwp,qz,mu,sigma))*gmunu(nu,rho)
     &                    -(epsrrmunu(qwm,qa,mu,rho)+
     &                    epsrrmunu(qwm,qz,mu,rho)+
     &                    epsrrmunu(qwp,qz,mu,rho))*gmunu(nu,sigma)
     &                    +(-epsrmunurho(qa,mu,rho,sigma)-
     &                    epsrmunurho(qz,mu,rho,sigma))*(qwm(nu)-qwp(nu))
     &                    +(-epsrmunurho(qwm,mu,rho,sigma)+
     &                    epsrmunurho(qwp,mu,rho,sigma))*qa(nu)+
     &                    (epsrrmunu(qwm,qa,mu,nu)
     &                    -epsrrmunu(qwm,qz,mu,nu)+
     &                    epsrrmunu(qwp,qa,mu,nu)-
     &                    epsrrmunu(qwp,qz,mu,nu))*gmunu(rho,sigma)
     &                    +(-epsrmunurho(qz,mu,nu,sigma)+
     &                    epsrmunurho(qa,mu,nu,sigma))*qwm(rho)
     &                    -(-epsrmunurho(qwm,mu,nu,sigma)-
     &                    epsrmunurho(qwp,mu,nu,sigma))*(qa(rho)-qz(rho))
     &                    +(-epsrmunurho(qz,mu,nu,rho)+
     &                    epsrmunurho(qa,mu,nu,rho))*qwp(sigma)
     &                    -(-epsrmunurho(qwm,mu,nu,rho)-
     &                    epsrmunurho(qwp,mu,nu,rho))*(qa(sigma)-qz(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (dotrr(qwm,qa)-dotrr(qwm,qz)-dotrr(qwp,qa)+
     &                    dotrr(qwp,qz)))
                  endif
                  if (fdwtilde .ne. 0) then
                     wwaz(mu,nu,rho,sigma)=wwaz(mu,nu,rho,sigma)+
     &                    fdwtilde*2.*gw**4*cosw*sinw*(2.*
     &                    epsrrmunu(qwm,qwp,rho,sigma)
     &                    *gmunu(mu,nu)+epsrrmunu(qz,qwm,nu,sigma)*
     &                    gmunu(mu,rho)
     &                    +epsrrmunu(qz,qwp,nu,rho)*gmunu(mu,sigma)-
     &                    epsrmunurho(qwp,nu,rho,sigma)
     &                    *(qz(mu)+qwm(mu)-qwp(mu))+
     &                    epsrmunurho(qwm,nu,rho,sigma)*(qz(mu)-qwm(mu)+
     &                    qwp(mu))+epsrrmunu(qa,qwm,mu,sigma)
     &                    *gmunu(nu,rho)+epsrrmunu(qa,qwp,mu,rho)*
     &                    gmunu(nu,sigma)-epsrmunurho(qwp,mu,rho,sigma)*
     &                    (qa(nu)+qwm(nu)-qwp(nu))
     &                    +epsrmunurho(qwm,mu,rho,sigma)*(qa(nu)-qwm(nu)
     &                    +qwp(nu))
     &                    -2.*epsrrmunu(qa,qz,mu,nu)*gmunu(rho,sigma)+
     &                    epsrmunurho(qz,mu,nu,sigma)*(qa(rho)-qz(rho)+
     &                    qwm(rho))
     &                    -epsrmunurho(qa,mu,nu,sigma)*(-qa(rho)+qz(rho)
     &                    +qwm(rho))+epsrmunurho(qz,mu,nu,rho)
     &                    *(qa(sigma)-qz(sigma)+qwp(sigma))-
     &                    epsrmunurho(qa,mu,nu,rho)*(-qa(sigma)+
     &                    qz(sigma)+qwp(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (-dotrr(qa,qwm)+dotrr(qa,qwp)+dotrr(qz,qwm)
     &                    -dotrr(qz,qwp)))
                  endif     
               enddo
            enddo
         enddo
      enddo      
      end
      
c --------------------------------------------------------------
c      WWAA-Vertices
c --------------------------------------------------------------

      subroutine wwaa_anomal(qwp,qwm,qa1,qa2,wwaa)
c calculates anomalous wwaa couplings
c qwp: momentum of the W+, qwp(rho)
c qwm: momentum of the W-, qwm(sigma)
c qa1,qza: momenta of the photons: qa1(mu),qa2(nu)
c wwaa: wwa-tensor            
      IMPLICIT NONE
      double precision qwp(0:3),qwm(0:3),qa1(0:3),qa2(0:3),dotrr
      double precision sinw,cosw
      double precision epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      integer mu,nu,rho,sigma
      complex*16 wwaa(0:3,0:3,0:3,0:3),im
      external dotrr,epsrrrmu,epsrrmunu,epsrmunurho,epsmunurhosigma
      include "gmunu.inc"
      include "coupl.inc"
      include "an_couplings.inc"   

      im=(0.,1.)
      sinw=sqrt(sin2w)
      cosw=sqrt(1-sin2w)
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3
                  wwaa(mu,nu,rho,sigma)=(0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      do mu=0,3
         do nu=0,3
            do rho=0,3
               do sigma=0,3      
                  if (fwww .ne. 0) then
                     wwaa(mu,nu,rho,sigma)=wwaa(mu,nu,rho,sigma)-
     &                    fwww*3./2.*sin2w*gw**4
     &                    *((qa2(mu)*(qwm(nu)+qwp(nu))+(qwm(mu)+qwp(mu))
     &                    *qa1(nu))*gmunu(rho,sigma)
     &                    -gmunu(mu,nu)*(dotrr(qwm,qa1)+dotrr(qwm,qa2)+
     &                    dotrr(qwp,qa1)+dotrr(qwp,qa2))*gmunu(rho,sigma)
     &                    -gmunu(nu,sigma)*(qa2(mu)*qwm(rho)+qwm(mu)*
     &                    (qa1(rho)-qa2(rho))
     &                    +qwp(mu)*qa2(rho))-gmunu(mu,sigma)*(qa1(nu)*
     &                    qwm(rho)+qwp(nu)*qa1(rho)
     &                    +qwm(nu)*(qa2(rho)-qa1(rho)))-gmunu(nu,rho)*
     &                    (qa2(mu)*qwp(sigma)+qwp(mu)*(qa1(sigma)
     &                    -qa2(sigma))+qwm(mu)*qa2(sigma))-gmunu(mu,rho)
     &                    *(qa1(nu)*qwp(sigma)+qwm(nu)*qa1(sigma)
     &                    +qwp(nu)*(qa2(sigma)-qa1(sigma)))+gmunu(mu,nu)
     &                    *((qa1(rho)+qa2(rho))*qwp(sigma)
     &                    +qwm(rho)*(qa1(sigma)+qa2(sigma)))+
     &                    gmunu(mu,sigma)*gmunu(nu,rho)*(dotrr(qwm,qa2)
     &                    +dotrr(qwp,qa1))+gmunu(mu,rho)*gmunu(nu,sigma)
     &                    *(dotrr(qwm,qa1)+dotrr(qwp,qa2)))
                  endif      
                  if (fwwwtilde .ne. 0) then
                     wwaa(mu,nu,rho,sigma)=wwaa(mu,nu,rho,sigma)+
     &                    fwwwtilde/2.*sin2w*gw**4*
     &                    ((epsrrmunu(qwm,qa1,rho,sigma)
     &                    +epsrrmunu(qwm,qa2,rho,sigma)-
     &                    epsrrmunu(qwp,qa1,rho,sigma)
     &                    -epsrrmunu(qwp,qa2,rho,sigma))*gmunu(mu,nu)-
     &                    (epsrrmunu(qwm,qa1,nu,sigma)+
     &                    epsrrmunu(qwp,qa1,nu,sigma)
     &                    +epsrrmunu(qwp,qa2,nu,sigma))*gmunu(mu,rho)
     &                    -(epsrrmunu(qwm,qa1,nu,rho)+
     &                    epsrrmunu(qwm,qa2,nu,rho)+
     &                    epsrrmunu(qwp,qa1,nu,rho))*gmunu(mu,sigma)
     &                    +(-epsrmunurho(qa1,nu,rho,sigma)-
     &                    epsrmunurho(qa2,nu,rho,sigma))*(qwm(mu)-qwp(mu))
     &                    +(-epsrmunurho(qwm,nu,rho,sigma)+
     &                    epsrmunurho(qwp,nu,rho,sigma))*qa2(mu)-
     &                    (epsrrmunu(qwm,qa2,mu,sigma)
     &                    +epsrrmunu(qwp,qa1,mu,sigma)+
     &                    epsrrmunu(qwp,qa2,mu,sigma))*gmunu(nu,rho)
     &                    -(epsrrmunu(qwm,qa1,mu,rho)+
     &                    epsrrmunu(qwm,qa2,mu,rho)+
     &                    epsrrmunu(qwp,qa2,mu,rho))*gmunu(nu,sigma)
     &                    +(-epsrmunurho(qa1,mu,rho,sigma)-
     &                    epsrmunurho(qa2,mu,rho,sigma))*(qwm(nu)-qwp(nu))
     &                    +(-epsrmunurho(qwm,mu,rho,sigma)+
     &                    epsrmunurho(qwp,mu,rho,sigma))*qa1(nu)+
     &                    (epsrrmunu(qwm,qa1,mu,nu)
     &                    -epsrrmunu(qwm,qa2,mu,nu)+
     &                    epsrrmunu(qwp,qa1,mu,nu)-
     &                    epsrrmunu(qwp,qa2,mu,nu))*gmunu(rho,sigma)
     &                    +(-epsrmunurho(qa2,mu,nu,sigma)+
     &                    epsrmunurho(qa1,mu,nu,sigma))*qwm(rho)
     &                    -(-epsrmunurho(qwm,mu,nu,sigma)-
     &                    epsrmunurho(qwp,mu,nu,sigma))*(qa1(rho)-qa2(rho))
     &                    +(-epsrmunurho(qa2,mu,nu,rho)+
     &                    epsrmunurho(qa1,mu,nu,rho))*qwp(sigma)
     &                    -(-epsrmunurho(qwm,mu,nu,rho)-
     &                    epsrmunurho(qwp,mu,nu,rho))*(qa1(sigma)-
     &                    qa2(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (dotrr(qwm,qa1)-dotrr(qwm,qa2)-dotrr(qwp,qa1)+
     &                    dotrr(qwp,qa2)))
                  endif
                  if (fdwtilde .ne. 0) then
                     wwaa(mu,nu,rho,sigma)=wwaa(mu,nu,rho,sigma)+
     &                    fdwtilde*2.*gw**4*sin2w*(2.*
     &                    epsrrmunu(qwm,qwp,rho,sigma)*gmunu(mu,nu)
     &                    +epsrrmunu(qa2,qwm,nu,sigma)*gmunu(mu,rho)
     &                    +epsrrmunu(qa2,qwp,nu,rho)*gmunu(mu,sigma)-
     &                    epsrmunurho(qwp,nu,rho,sigma)*(qa2(mu)+qwm(mu)
     &                    -qwp(mu))
     &                    +epsrmunurho(qwm,nu,rho,sigma)*(qa2(mu)-
     &                    qwm(mu)+qwp(mu))+epsrrmunu(qa1,qwm,mu,sigma)
     &                    *gmunu(nu,rho)+epsrrmunu(qa1,qwp,mu,rho)*
     &                    gmunu(nu,sigma)-epsrmunurho(qwp,mu,rho,sigma)*
     &                    (qa1(nu)
     &                    +qwm(nu)-qwp(nu))+
     &                    epsrmunurho(qwm,mu,rho,sigma)*(qa1(nu)-
     &                    qwm(nu)+qwp(nu))
     &                    -2.*epsrrmunu(qa1,qa2,mu,nu)*gmunu(rho,sigma)+
     &                    epsrmunurho(qa2,mu,nu,sigma)*(qa1(rho)-
     &                    qa2(rho)+qwm(rho))
     &                    -epsrmunurho(qa1,mu,nu,sigma)*(-qa1(rho)+
     &                    qa2(rho)+qwm(rho))+epsrmunurho(qa2,mu,nu,rho)
     &                    *(qa1(sigma)-qa2(sigma)+qwp(sigma))-
     &                    epsrmunurho(qa1,mu,nu,rho)*(-qa1(sigma)+
     &                    qa2(sigma)+qwp(sigma))
     &                    +epsmunurhosigma(mu,nu,rho,sigma)*
     &                    (-dotrr(qa1,qwm)+dotrr(qa1,qwp)+dotrr(qa2,qwm)
     &                    -dotrr(qa2,qwp)))
                  endif
               enddo
            enddo
         enddo
      enddo
      end
      

********************************************************************************
********************************************************************************

** This subroutine applies the formfactor for the dim-6 operators for anomalous
** gauge boson couplings.
** Input via anomV.dat.
** The coefficients FB, FW and FWWW (which parametrise the WWZ and WWgamma 
** couplings) can have individual formfactors
      
      subroutine anomal_formfactor(p1,p2)

      implicit none
      include "mssm.inc"
      include "an_couplings.inc"


      real*8 p1(0:3),p2(0:3),s
      double precision ff, ff_fwww, ff_fw, ff_fb, ff_KZ, ff_KA



      if (formfact) then

** energy scale:
         s=(p1(0)+p2(0))**2-(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-
     -        (p1(3)+p2(3))**2
         s=abs(s)   ! newly added

** universal formfactor:
         ff=(1./(1.+s/ffmassscale2))**DBLE(ffexponent)

** the individual formfactors:
         ff_fwww = (1./(1.+s/massscale2FWWW))**DBLE(ffexpFWWW)
         ff_fw = (1./(1.+s/massscale2FW))**DBLE(ffexpFW)
         if (trianom .eq. 1) then
            ff_fb = (1./(1.+s/massscale2FB))**DBLE(ffexpFB)
         else
c$$$            ff_fb = CW2*MZ2*fw_0*ff_fw - 
c$$$     &           2d0*zDkappa0_0/((1D0 + s/massscale2KZ)**ffexpKZ)
c$$$            ff_fb = ff_fb/(CW2*MZ2*fw_0 - 2d0*zDkappa0_0)
            ff_fb = 2d0*aDkappa0_0/(1d0 + s/massscale2KA)**dble(ffexpKA)
     &           - MW2*fw_0*ff_fw
            ff_fb = ff_fb/(2d0*aDkappa0_0 - MW2*fw_0)
         end if

c      print*, Sqrt(s),ff!lambda,ff,p1(0)+p2(0)
c      write(11,*) dsig(1),ff

         if (fbw .ne. 0) then
            fbw=fbw_0
            fbw=fbw*ff
         endif
         if (fdw .ne. 0) then
            fdw=fdw_0
            fdw=fdw*ff
         endif
         if (fwww .ne. 0) then
            fwww=fwww_0
            fwww = fwww*ff_fwww
         endif
         if (fww .ne. 0) then
            fww=fww_0
            fww=fww*ff
         endif
         if (fbb .ne. 0) then
            fbb=fbb_0
            fbb=fbb*ff
         endif
         if (fw .ne. 0) then
            fw=fw_0
            fw = fw*ff_fw
         endif
         if (fb .ne. 0) then
            fb=fb_0
            fb = fb*ff_fb
         endif
         if (fwwtilde .ne. 0) then
            fwwtilde=fwwtilde_0
            fwwtilde=fwwtilde*ff
         endif
         if (fbwtilde .ne. 0) then
            fbwtilde=fbwtilde_0
            fbwtilde=fbwtilde*ff
         endif
         if (fbbtilde .ne. 0) then
            fbbtilde=fbbtilde_0
            fbbtilde=fbbtilde*ff
         endif
         if (fwtilde .ne. 0) then
            fwtilde=fwtilde_0
            fwtilde=fwtilde*ff
         endif
         if (fbtilde .ne. 0) then
            fbtilde=fbtilde_0
            fbtilde=fbtilde*ff
         endif
         if (fwwwtilde .ne. 0) then
            fwwwtilde=fwwwtilde_0
            fwwwtilde=fwwwtilde*ff
         endif
         if (fdwtilde .ne. 0) then
            fdwtilde=fdwtilde_0
            fdwtilde=fdwtilde*ff
         endif


** Now moving to the trianom=2 parametrisation
         if (trianom .eq. 2) then
            ff_KA = (1./(1.+s/massscale2KA))**DBLE(ffexpKA)
            ff_KZ = (1./(1.+s/massscale2KZ))**DBLE(ffexpKZ)
         else
! numerator is multiplied by appropriate formfactors
            ff_KA = (fw + fb)/(fw_0 + fb_0)
            ff_KZ = (CW2*fw - SW2*fb)/(fw_0*CW2 - SW2*fb_0)  
         end if
         if (lambda0 .ne. 0D0) then
            lambda0 = lambda0_0
            lambda0 = lambda0*ff_fwww
         end if
         if (zDg0 .ne. 0D0) then
            zDg0 = zDg0_0
            zDg0 = zDg0*ff_fw
         end if
         if ((aDkappa0 .ne. 0d0) .or. (zDkappa0 .ne. 0d0)) then
            aDkappa0 = aDkappa0_0
            aDkappa0 = aDkappa0*ff_KA
            zDkappa0 = zDkappa0_0
            zDkappa0 = zDkappa0*ff_KZ
         end if


      endif


      end



********************************************************************************
                        
