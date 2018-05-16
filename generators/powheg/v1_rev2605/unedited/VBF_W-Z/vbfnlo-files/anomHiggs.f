c==================  subroutine read_anomHVV ================================
c
c       this subroutine reads in anomalous couplings for the higgs 
c       boson in different parametrisation and transforms them 
c       into the Parametrisation of Phys. Lett. B591, 297
c
c       Note that in some cases (WW) both anomV.dat and anom_HVV.dat are
c       read-in.  In this case, the shared values fb,fw,dkgam,dg1Z are 
c       taken from anomV.dat.  If switch1=true in anom_HVV.dat, the HVV
c       couplings are taken directly from these inputs, BUT there is a
c       consistency check with the existing shared values from anomV.dat
c
c========================================================================
c
      subroutine read_anomHVV
      implicit none

      include "an_couplings.inc"
      include "mssm.inc"
      include "../../include/pwhg_math.h"
      include "global.inc"

      real*8 g5hvv(2,1:4,1:4), lambda5, 
     &     treefacW, treefacZ, treefac,loopfac, 
     &     d, dodd, dB,dBo, dg1Z, dkgam, kgamo,
     &     fWWe, fWWo, fBBe, fBBo, fWe, fBe, fBo
      common /hcoupl/ g5hvv, lambda5
      common /lhcoup/ treefacW, treefacZ, loopfac
      integer ffac
      real*8 m2ff,mff
      logical lff, switch1, switch2, switch3
      common/ formfacmass/ m2ff,mff,ffac,lff
      integer i,j,k, hvv1, hvv2
      character*50 gtext

      real*8 vbfnloinput
      external vbfnloinput

c     elektroweak parameters used here:
      double precision s2, al, g, lff1, switch11
!       real*8 pi
!       parameter(pi =3.141592653589793d0) 

      double precision  alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf
      common /bkopin/   alfas, xmt, alfa, xmz, xmw, sin2w, xmh, gf

*      s2 = 0.2312d0
*      mW = 80.425d0
*      mZ = 91.189d0
*      al = 0.00729927d0

      s2 = sin2w
      al = Alfa

      g = gf

c     initialize couplings to zero
      do i=1,4
         do j =1,4
            g5hvv(1,i,j) = 0.0d0
            g5hvv(2,i,j) = 0.0d0
         enddo
      enddo
      d = 0.0d0
      dB = 0.0d0
      dg1Z = 0.0d0
      dkgam = 0.0d0
      dodd = 0.0d0
      dBo = 0.0d0
      kgamo = 0.0d0
      fWWe = 0.0d0
      fBBe = 0.0d0
      fWe = 0.0d0
      fBe = 0.0d0
      fWWo = 0.0d0
      fBBo = 0.0d0
      fBo = 0.0d0
      hvv1 = 0
      hvv2 = 0
      switch1 = .false.
      switch2 = .false.
      switch3 = .false.

      print *," "
      print *,"  Information on anomalous Higgs coupling parameters  "
      print *,"------------------------------------------------------"
!       call loadfile("anom_HVV.dat",.true.)

c     Form Factor
       lff=.false.
       lff1=vbfnloinput("FORMFACTOR")
       if (lff1.eq.1) lff=.true.
       
       if(lff) then
          mff=vbfnloinput("MASS_SCALE")
          ffac=vbfnloinput("FFAC")
          m2ff = mff**2
       endif


c     Parametrisation: Phys. Lett. B591, 297
       switch1=.false.
       switch11=vbfnloinput("PARAMETR1")
       if (switch11.eq.1) switch1=.true.
       if(switch1) then
          lambda5=vbfnloinput("LAMBDA5")
c     read cp even couplings
          g5hvv(1,3,4)=vbfnloinput("G5E_HWW")
          g5hvv(1,4,3) = g5hvv(1,3,4)
          g5hvv(1,2,2)=vbfnloinput("G5E_HZZ")
          g5hvv(1,1,1)=vbfnloinput("G5E_HGG")
          g5hvv(1,1,2)= vbfnloinput("G5E_HGZ")
          g5hvv(1,2,1) = g5hvv(1,1,2)
c     read cp odd couplings
          g5hvv(2,3,4)=vbfnloinput("G5O_HWW")
          g5hvv(2,4,3) = g5hvv(2,3,4)
          g5hvv(2,2,2)=vbfnloinput("G5O_HZZ")
          g5hvv(2,1,1)=vbfnloinput("G5O_HGG")
          g5hvv(2,1,2)=vbfnloinput("G5O_HGZ")
          g5hvv(2,2,1) = g5hvv(2,1,2)
       endif

c     Parametrisation of the L3-Collaboration
       switch2=.false.
       switch11=vbfnloinput("PARAMETR2")
       if (switch11.eq.1) switch2=.true.

       if(switch2) then
c     read cp even couplings
          d=vbfnloinput("D_EVEN")
          dB=vbfnloinput("DB_EVEN")
          dg1Z=vbfnloinput("DG1Z_EVEN")
          dkgam=vbfnloinput("DKGAM_EVEN")
c     read cp odd couplings
          dodd=vbfnloinput("D_ODD")
          dBo=vbfnloinput("DB_ODD")
          kgamo=vbfnloinput("KGAM_ODD")
          hvv1=vbfnloinput("HVV1")
          if (switch1) then
             write(*,*)'WARNING! You have set both PARAMETR2 and'
             write(*,*)'PARAMETR1 to .true. in anom_HVV.dat.'
             write(*,*)'We will use parametrisation 1.'
          end if
       endif


c     Parametrisation: Phys. Lett. B318, 155
       switch3=.false.
       switch11=vbfnloinput("PARAMETR3")
       if (switch11.eq.1) switch3=.true.
       if(switch3) then
c     read cp even couplings
          fWWe=vbfnloinput("FWW_EVEN")
          fBBe=vbfnloinput("FBB_EVEN")
          fWe=vbfnloinput("FW_EVEN")
          fBe=vbfnloinput("FB_EVEN")
c     read cp odd couplings
          fWWo=vbfnloinput("FWW_ODD")
          fBBo=vbfnloinput("FBB_ODD")
          fBo=vbfnloinput("FB_ODD")
          hvv2=vbfnloinput("HVV2")
          if (switch1) then
             write(*,*)'WARNING! You have set both PARAMETR3 and'
             write(*,*)'PARAMETR1 to .true. in anom_HVV.dat.'
             write(*,*)'We will use parametrisation 1.'
          end if
          if (switch2) then
             write(*,*)'WARNING! You have set both PARAMETR3 and'
             write(*,*)'PARAMETR2 to .true. in anom_HVV.dat.'
             write(*,*)'We will use parametrisation 2.'
             switch3 = .false.
          end if
       endif

c     Loopfactor and Treefactor
      treefacZ=vbfnloinput("TREEFACZ")
      treefacW=vbfnloinput("TREEFACW")
      if ((treefacZ .eq. -9999d0) .or. (treefacW .eq. -9999d0)) then
         treefac=vbfnloinput("TREEFAC")
         if (treefacZ .eq. -9999d0) treefacZ = treefac
         if (treefacW .eq. -9999d0) treefacW = treefac
      end if
      if(treefacZ.ne.1d0) then
         write(6,*) " SM-type Higgs-ZZ coupling multiplied by ",treefacZ 
      endif
      if(treefacW.ne.1d0) then
         write(6,*) " SM-type Higgs-WW coupling multiplied by ",treefacW
      endif
      loopfac=vbfnloinput("LOOPFAC")
      if(loopfac.ne.1d0) then
         write(6,*) " SM contributions to h-gamma-gamma and"
         write(6,*) " h-Z-gamma loop induced Higgs couplings"
         write(6,*) " multiplied by ",loopfac 
      endif
      print *," "

      if((.not.switch1).and.(.not.switch2).and.(.not.switch3)) then
         lambda5 = 480.0d0
         write(*,*)'WARNING! You have chosen to include anomalous '
         write(*,*)'couplings, but have not chosen a parametrisation.'
         write(*,*)'Anomalous couplings will be set to 0.'
      endif

 800  format(a)

** Converting from alternative parameterisation
      if (.not. switch1) then
         do i=1,4
            do j=1,4
               g5hvv(2,i,j) = 0.0d0
               g5hvv(1,i,j) = 0.0d0
            enddo
         enddo
         call anomH_convert(hvv1,hvv2,switch2,switch3,
     &        d, dodd, dB, dBo, dg1Z, dkgam, kgamo,
     &        fWWe, fWWo, fBBe, fBBo, fWe, fBe, fBo)
      end if


c$$$** Debugging output:
c$$$      do i = 1, 2
c$$$         do j = 1, 4
c$$$            do k = 1, 4
c$$$               write(*,*)'ijk g5hvv =', i, j, k, g5hvv(i,j,k)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      stop
      

      end


c===================subroutine fillanomhcoupl============================
c
c     this subroutine fills the array ahvv with g5hvv, lambda5
c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
c     
c     v1, v2 = 1: photon
c     v1, v2 = 2: Z0 boson
c     v1, v2 = 3,4: W boson
c
c     init = 1: only initialization to zero
c     else:     fill the array ahvv with g5hvv, lambda
c
c========================================================================
c
      subroutine fillanomhcoupl(init,dum_hzg,dum_hgg,loopfac)
      
      implicit none

! #include "VBFNLO/utilities/global.inc"

      real*8  g5hvv(2,1:4,1:4), lambda5, loopfac
      integer ffac
      real*8 m2ff,mff
      logical lff
      common/ formfacmass/ m2ff,mff,ffac,lff
      integer i,j,k,init
      complex*16 dum_hzg,dum_hgg !sm 1-loop contributions 
      common /hcoupl/ g5hvv, lambda5
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL


      if(init.eq.1) then   ! initialze ahvv only
         do i = 1,4
            do j = 1,4
               ahvv(1,i,j) = (0.0d0,0.0d0)
               ahvv(2,i,j) = (0.0d0,0.0d0)
               ahvv(3,i,j) = (0.0d0,0.0d0)
            enddo
         enddo

c         if (with_anomHiggs) then
c            write(*,*)
c            write(6,*) " Anomalous HVV couplings are being used."
c            write(*,*)
c            write(*,"(a)") " G5HVV: "
c            write(*,"(a)") " CP-EVEN COUPLINGS: "
c            write(*,"(a6,g12.4)") " HGAMGAM: ",g5hvv(1,1,1)
c            write(*,"(a6,g12.4,a2,g12.4)") " HZGAM: ",g5hvv(1,1,2),", ",
c     &                                   g5hvv(1,2,1)
c            write(*,"(a6,g12.4)") " HZZ: ",g5hvv(1,2,2)
c            write(*,"(a6,g12.4,a2,g12.4)") " HWW: ",g5hvv(1,3,4),
c     &              ", ",g5hvv(1,4,3)
c            write(*,"(a)") " CP-ODD COUPLINGS: "
c            write(*,"(a6,g12.4)") " HGAMMA: ",g5hvv(2,1,1)
c            write(*,"(a6,g12.4,a2,g12.4)") " HZGAM: ",g5hvv(2,1,2),", ",
c     &                                   g5hvv(2,2,1)
c            write(*,"(a6,g12.4)") " HZZ: ",g5hvv(2,2,2)
c            write(*,"(a6,g12.4,a2,g12.4)") " HWW: ",g5hvv(2,3,4),
c     &              ", ",g5hvv(2,4,3)
c            write(*,*)
c            write(*,"(a11,g12.4,a3)") " LAMBDA5 = ",lambda5,"GeV"
c            write(*,*)

c            write(*,"(a11,g12.4)") " TREEFAC = ",treefac
c            write(*,"(a11,g12.4)") " LOOPFAC = ",loopfac
c            write(*,*)
c            write(*,"(a7,l4)") " LFF = ",lff  
c            write(*,"(a8,g12.4)") " M2FF = ",m2ff
c            write(*,"(a8,i4)") " FFAC = ",ffac
c            write(*,*)
c         else
c            write(*,*)
c            lff = .false. 
c            mff = 1.0d0
c            lambda5 = 1.0d0
c            do i =1,4
c               do j = 1,4
c                  g5hvv(2,i,j) =0.0d0        
c                  g5hvv(1,i,j) =0.0d0
c               enddo 
c            enddo
c            treefac = 1.0d0
c            loopfac = 1.0d0
c         endif

         return
      endif
        
c     fill ahvv coefficients with lambda5, g5hvv
c        h-w-w anomalous couplings
      do i = 3,4
         do j = 3,4
            ahvv(2,i,j) = -g5hvv(1,i,j)/lambda5*(2.0d0,0.0d0)
            ahvv(3,i,j) = g5hvv(2,i,j)/lambda5*(2.0d0,0.0d0)
         enddo
      enddo
c
c        h-z-z anomalous couplings
      ahvv(2,2,2) = -2.0d0*g5hvv(1,2,2)/lambda5*(1.0d0,0.0d0)
      ahvv(3,2,2) = 2.0d0*g5hvv(2,2,2)/lambda5*(1.0d0,0.0d0)
c
c        h-gamma-gamma
c        SM + CPeven
      ahvv(2,1,1) = loopfac*dum_hgg - 2.0d0*g5hvv(1,1,1)/lambda5
     $     *(1.0d0,0.0d0) 
c        CP odd
      ahvv(3,1,1) = 2.0d0*g5hvv(2,1,1)/lambda5 *(1.0d0,0.0d0)
c
c        h-z-gamma
c        SM + CP-even 
      ahvv(2,1,2) = -(loopfac*dum_hzg + 2.0d0*g5hvv(1,1,2)/lambda5
     $     *(1.0d0,0.0d0))
      ahvv(2,2,1) = -(loopfac*dum_hzg + 2.0d0*g5hvv(1,2,1)/lambda5
     $     *(1.0d0,0.0d0))
c     CP-odd
      ahvv(3,1,2) = 2.0d0*g5hvv(2,1,2)/lambda5*(1.0d0,0.0d0)
      ahvv(3,2,1) = 2.0d0*g5hvv(2,2,1)/lambda5*(1.0d0,0.0d0)

      return

      end



********************************************************************************
********************************************************************************

c     This subroutine computes the decay width of 
c     h -> Z/gamma Z/gamma -> f1 f1bar f3 f3bar for 
c     flav1,flav3
c     Final state fermions are treated as non-identical particles.
c
      subroutine flavor_sum(width)
c     only works for h->Z*Z* decay

      implicit none

      include "koppln_ew.inc"
      include "mssm.inc"

      real*8 width,dum
      real*8 matrix(11,11),mass(11),Nc(4),gm(4),matrix1(4,4)
      real*8 Nf(4),factor
c     nu_e,e,nu_mu,mu,nu_tau,tau,u,d,c, s, b
c     1    2   3   4    5    6   7 8 9 10 11
      real*8 flavor(11)
      integer i,j,k,l
      real*8 pi
      parameter(pi=3.141592653589793d0)
      logical lcomputeALL
      data flavor/1,2,1,2,1,2,3,4,3,4,4/
      parameter(lcomputeALL=.false.)

c     color factors 
      Nc(1) = 1.d0
      Nc(2) = 1.d0
      Nc(3) = 3.d0*(1.d0 + ALFAS/PI)   
      Nc(4) = 3.d0*(1.d0 + ALFAS/PI)

c     Number of fermion in generation
      Nf(1) = 3.0d0
      Nf(2) = 3.0d0 
      Nf(3) = 2.0d0             !no top
      Nf(4) = 3.0d0

c     compute geometric means of masses 
      gm(1) = 0d0     ! NEUTRINO MASS
      gm(2) = (Mf(2,1)*Mf(2,2)*Mf(2,3))**(1.0d0/3.d0)
      gm(3) = (Mf(3,1)*Mf(3,2))**(1.0d0/2.d0)
      gm(4) = (Mf(4,1)*Mf(4,2)*Mf(4,3))**(1.0d0/3.d0)
c      print*,"geom. mean of masses=",gm
c      stop

c    setting mass
      mass(1) = Mf(1,1)
      mass(2) = Mf(2,1)
      mass(3) = Mf(1,2)
      mass(4) = Mf(2,2)
      mass(5) = Mf(1,3)
      mass(6) = Mf(2,3)
      mass(7) = Mf(3,1)
      mass(8) = Mf(4,1)
      mass(9) = Mf(3,2)
      mass(10) = Mf(4,2)
      mass(11) = Mf(4,3)

      if(lcomputeALL) then
         dum = 0.0d0
         do k=1,11
            do l = k,11
               if (l.eq.k) then
                  factor = 2D0
               else
                  factor = 1D0
               end if
               i = flavor(k)
               j = flavor(l)
               call hzgammawidth(mass(k),mass(l),i,j,matrix(k,l))
               dum = matrix(k,l)*Nc(i)*Nc(j)/factor + dum
            enddo
         enddo
      else
         dum =0.0d0
         factor = 2.0d0
         do i=1,4
            do j = 1,4         
               call hzgammawidth(gm(i),gm(j),i,j,matrix1(i,j))
               dum = matrix1(i,j)*Nc(i)*Nc(j)*Nf(i)*Nf(j)/factor + dum
            enddo
         enddo
      endif
      width = dum

      return
      end


********************************************************************************
c
*** Subroutine called in the calculation of h -> ZZ -> 4f, from flavor_sum
*** mf1 and mf3 and the fermion masses 
*** iflav1 and iflav3 are the fermion flavours 
*** flavor_sum runs over all of the end-state fermions

      subroutine hzgammawidth(mf1,mf3,iflav1,iflav3,width)

      implicit none

      real*8 width,mf1,mf3,mh,qsq1,qsq2,ss1,ss2
      integer flav1,flav3
      integer iflav1,iflav3      
      common/flavors/flav1,flav3
      real*8 mu1,mu3,x1,x2
      common/masses/mu1,mu3
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id 
         
      id = 2                    !always -- i.e. Z
      flav1 = iflav1
      flav3 = iflav3
      
      mu1 = mf1
      mu3 = mf3
      

      x1 = 0.0d0
      x2 = 1.0d0

      call quad2dim(x1,x2,ss2)
      width = ss2 

      return
      end
c    
c
********************************************************************************
*** This is the function that's (eventually) called to calculate the 
*** H -> ZZ -> 4f width
*** It's (I think) the equivalent of the function 'fun' in koppln.F that's 
*** used to calculated H -> ZZ in the SM
      
      real*8 function dgam(xx,yy)

      implicit none

      integer flav1,flav3,i,j
      common/flavors/flav1,flav3
      real*8 qsq1,qsq2,xx,yy
      double complex F1,F2,F3,F1c,F2c,F3c
      real*8 mu1,mu3
      common/masses/mu1,mu3
      real*8 dum1
      real*8 ints(3,3),lamb1,lamb2,factor,pi
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      parameter(pi=3.141592653589793d0)
      real*8 Jac2,Jac1
      real*8 lowq2(2),hiq2(2)
 
      
c     convention here is 2 is yy
c     1 is xx 
c      qsq1 = yy
c      qsq2 = xx
      lowq2(1) = 4.0d0*mu1**2
      hiq2(1) = xm2(6)
c      print*,"xm2(6)=",xm2(6)
         
      call com_Jq2(1,lowq2,hiq2,xx,qsq1,Jac1) !outer integral
      
c     compute bound on inner integration dy or 2
      lowq2(2) = 4.0d0*mu3**2
      hiq2(2) = (sqrt(xm2(6))-sqrt(qsq1))**2

      call com_Jq2(2,lowq2,hiq2,yy,qsq2,Jac2) !inner integral

c      print*,"hiq2 in unit g(x,y)=",hiq2
c      print*,"lowq2 in unit g(x,y)=",lowq2
      dum1 = 0.0d0
      call computeIs(qsq1,qsq2,ints)

      do i=-1,1,2
         do j = -1,1,2

** this subroutine is in HtoWW.F.  It calculates the combinations of H-V and
** V-f couplings
            call computeFs(qsq1,qsq2,i,j,flav1,flav3,2,2,
     $           F1,F2,F3)
  
c            print*,"flav11,flav3=",flav1,flav3

            F1c = conjg(F1)
            F2c = conjg(F2)
            F3c = conjg(F3)

            dum1 = (dble(f1)**2+dimag(f1)**2)*ints(1,1)+
     $           (dble(f2)**2+dimag(f2)**2)*ints(2,2)+
     $           (dble(f3)**2+dimag(f3)**2)*ints(3,3)+
     $           2.0d0*dble(f1*f2c)*ints(2,1) + dum1
         enddo
      enddo

      lamb1 = qsq1/xm2(6)
      lamb2 = qsq2/xm2(6)
      call lambda1(1.0d0,lamb1,lamb2,factor)
c      print*,"dum=",dum
      dgam = dum1*factor*(1.0d0/(8.0d0*pi))*
     $     (1.0d0/(32.d0*pi**2))**2 
     $     *(1.d0/(2.0d0*pi))**2 * 1.d0/(2.0d0*sqrt(xm2(6)))

      dgam = dgam*Jac1*Jac2     !mult by jacobians

      return
      end
c      


********************************************************************************
c
*** Function used in calculation of H -> ZZ -> 4f

      real*8 function ddwidth(yy)
      implicit none
      real*8 yy,x,y,dgam
      common/xy/ x,y
      external dgam

      y = yy
      ddwidth = dgam(x,y)

      return
      end


********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      real*8 function dwidth(xx)
      implicit none
      real*8 eps
      parameter(eps = 1.0d-5)
      real*8 xx,y1,y2,x,y
      common /xy/ x,y
      real*8 ss,gaus,ddwidth      
      
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
     
      integer flav1,flav3      
      common/flavors/flav1,flav3
      
     
      external ddwidth
      external gaus     
c
      x = xx
      
      y1 = 0.0d0
      y2 = 1.0d0
      ss = gaus(ddwidth,y1,y2,eps)
        
      dwidth = ss
      return
      end


********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      subroutine quad2dim(x1,x2,ss)
      implicit none
      integer i
      real*8 eps
      parameter(eps = -1.0d-2)
      real*8 ss,x1,x2,dwidth,gaus2
      external dwidth
      external gaus2
     
      ss = gaus2(dwidth,x1,x2,eps)  

      return
      end



********************************************************************************

*** Function used in calculation of H -> ZZ -> 4f

      subroutine computeIs(qsq1,qsq2,arrayOFints)

      implicit none

      integer i,j
      real*8 qsq1,qsq2
      real*8 arrayOFints(3,3)
      real*8 q1q2,pi               !q1.q2
      parameter(pi =3.141592653589793d0)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c     initialize
      do i=1,3
         do j=1,3
            arrayOFints(i,j) = 0.0d0
         enddo
      enddo
c     compute q1.q2
      q1q2 = (xm2(6) - qsq1 - qsq2)/2.d0
c
c      Fill arrayOFints
      arrayOFints(1,1) = (64.d0/9.d0)*pi**2 *
     1     (q1q2**2 + 2.0d0*qsq1*qsq2)
      arrayOFints(1,2) = (64.d0/3.d0)*pi**2* qsq1*qsq2* q1q2
      arrayOFints(3,3) = (128.d0/9.d0)*pi**2* 
     $     qsq1*qsq2*(q1q2**2 - qsq1*qsq2)
      arrayOFints(2,1) = arrayOFints(1,2)
      arrayOFints(2,2) = (64.d0/9.d0)*pi**2 *qsq1*qsq2 *
     1     (2.0d0*q1q2**2 + qsq1*qsq2)
      return
      end


********************************************************************************
c
c     compute jacobian and q2.  USed in H -> ZZ -> 4f and H -> Zgamma -> 2f gam

      subroutine com_Jq2(ich,lowq2,hiq2,r,qsq,J)

      implicit none

      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      real*8 r,qsq,J             ! J is the jacobian and r is the variable 
                                 ! being integrated over
      real*8 rmid,rmidl,a,c,b         !c is the midpoint on [a,b]
      real*8 x,xmin,q2cut,xmax
      integer i,ich
 
      real*8 lowq2(2),hiq2(2)
     
      logical lmap,logmapon
      parameter(logmapon=.true.)
c
      q2cut = 60.d0**2
      rmid = 0.2d0
      
      i = ich                   !1 or 2 = inner and outer bounds
c     if hiq2(i) is less than q2cut then let q2cut = hiq2(i)
      if(hiq2(i).lt.q2cut) then ! whole is integrated using ln map
         a = lowq2(i)
         b = hiq2(i)
         q2cut = hiq2(i)
         c = q2cut
         rmidl = 1.0d0
      else 
         a = lowq2(i)
         b = hiq2(i)
         c = q2cut
         rmidl = rmid
      endif
         
c      print*,"rmid=",rmid
c      print*,"imap=",imap
      lmap = (r.lt.rmidl).and.(a.gt.0.0d0).and.logmapon
      if(lmap) then        !ln map

         qsq = a * exp(R*(dlog(c/a))/rmidl)
         J = qsq *(dlog(c/a))/rmidl
c         print*,J,q2,a
c         print*,"imap=",imap
c         print*,"masses=",mu1,mu3
      else
         if(a.eq.0.0d0) then ! for a =0 always use bw map
            c = a
            rmidl =0.0d0
         else 
            c= q2cut
            rmidl = rmid
         endif
c         print*,"rmid=",rmid
         call calZ(c,xmin)
         call calZ(b,xmax)
         x = xmin +(r - rmidl)/(1.0d0 - rmidl)*(xmax - xmin)
         call calQ2(x,qsq)
         
         J = (xmax - xmin)/(1.0d0 - rmidl) 
     $        * xmgs(id)*(dtan(x)**2 + 1.0d0)
      endif
c      print*,J,q2
      return
      end
c


********************************************************************************

*** This is the subroutine that calculates the width 
***          h -> gamma Z -> gamma 2f 
*** It's called from width_hgamff (which sums over fermion type)
*** idf is fermion type, mf is fermion mass

      subroutine widthhgz(idf,mf,q2,width)
      implicit none
      integer idf,i,j               !1 for neutrino,2 for electron,3 for up, 4
                                    !for down
      real*8 cc(3),width(4),ratio(3) !1 is gg* 2 is gz* and 3 is the 
                                     ! interference term
      common/constants/cc
      real*8 q2,FF2(3),psfactor ! q squared of decaying boson
      real*8 delta,mg,mz,mh,mf,Pi,Lamb1,c1,c2,c3
      parameter(Pi= 3.141592653589793d0)
     
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      real*8 mf1
*      parameter(mf1=0d0)
      mf1 = mf

      mh = Sqrt(xm2(6))
c      write(6,*) mh
      mz = Sqrt(xm2(2))
      mg = xmg(2)
c      mg =0.001d0
c     Put coupling factors here:
      do i =1,3
         FF2(i) = 0.0d0
         cc(i) = 0.0d0
      enddo

c      print*,"lambda =" ,lambda5e,lambda5o
* H -> gamma gamma -> ff gamma
      cc(1) = 2.0d0*v(idf,2)**2 * 
     1     ( dble(ahvv(2,1,1))**2+dimag(ahvv(2,1,1))**2 +
     2       dble(ahvv(3,1,1))**2+dimag(ahvv(3,1,1))**2 ) 
* H -> Z gamma -> ff gamma
      cc(2) = (CLR(idf,2,-1)**2 + CLR(idf,2,+1)**2) * 
     1     ( dble(ahvv(2,2,1))**2+dimag(ahvv(2,2,1))**2 +
     2       dble(ahvv(3,2,1))**2 + dimag(ahvv(3,2,1))**2 ) 
* interference: H -> Z gamma -> ff gamma WITH H -> gamma gamma -> ff gamma
      cc(3) = 2.d0*(CLR(idf,2,-1) + CLR(idf,2,+1))*v(idf,1) * 
     1     ( dble(ahvv(2,1,1)*conjg(ahvv(2,1,2)))+
     2       dble(ahvv(3,1,1)*conjg(ahvv(3,1,2))) )
c
c      write(6,*) cc
c      print*,"mz",mz
c      print*,"mg",mg
c
      FF2(1) = cc(1)/q2**2 
      FF2(2) = cc(2)/((q2 - mz**2)**2 + mg**2)
      FF2(3) = cc(3)*(q2-mz**2)/((q2 - mz**2)**2 + mg**2)/q2
c
    
c      write(6,*) "ff2=",ff2
       psfactor =(1.0d0/(mh**3))*(1.0d0/(32.0d0*Pi**2))**2 * 
     1           (mh**2 - q2)*sqrt((q2 - 4.0d0*mf1**2)/q2) 

c     We use a fermion mass to regulate photon pole.
      do i = 1,3

         width(i) = FF2(i) * psfactor * (4.0d0/3.0d0)* Pi *
     1             (q2-mf1**2)*(mh**2 - q2)**2 ! m2s is computed with 
                                               ! massive fermions
c         write(6,*)"q2 =",q2
c         write(6,*)"dw/dq2 =",width(i),"for",i
      enddo
      
      width(4) = width(1)+width(2)+width(3)

      end


********************************************************************************

*** Function used in calculation of H -> gamma Z -> gamma ff

      function dw(xx)  ! interface with subroutine dwidth

      implicit none

      integer idf,iflav
      common/flav/ idf
      integer id                !id = 3,4 W"s 2 is Z
      common /partid/ id
      real*8 dw,xx,q2
      integer jj
      common/choice/jj
      real*8 mass(2),mf
* mass1 is set by the arguments of quad1d, mass(1) = mfermion, mass(2) = mhiggs
      common/mass1/mass
      real*8 ss,gaus,dpw(4)
      real*8 jac,lowq2(2),hiq2(2)
      external gaus
      
      id = 2
c     
      lowq2(1) = 4.0d0*mass(1)**2
      hiq2(1) = mass(2)**2
      call com_Jq2(1,lowq2,hiq2,xx,q2,Jac)
      iflav = idf
      mf = mass(1)
      call widthhgz(iflav,mf,q2,dpw)

c      dw  = dpw(2)+dpw(3)+dpw(1)
      dw = dpw(jj)*jac          !mult by jacobian

      return
      end


********************************************************************************
c
** USed in calculation of H -> Z gamma -> ff gamma

      subroutine quad1d(ii,iflav,mf,hmass,ss)
 
      implicit none

      integer idf,iflav
      common/flav/ idf
      integer ii,jj
      common/choice/jj
      real*8 mass(2),mf,hmass
      common/mass1/mass
      real*8 eps
      parameter(eps = 1.0d-7)
      real*8 ss,x1,x2,dw,gaus
      external dw
      external gaus

      jj = ii
      idf = iflav
      mass(1) = mf
      mass(2) = hmass
      x1 = 0.0d0
      x2 = 1.0d0
      ss = gaus(dw,x1,x2,eps)  
c      write(6,*) ss

      return

      end


********************************************************************************

*** Subroutine called from koppln.F to calculate 
***          h-> gamma/Z* gamma -> f fbar gamma
*** for anomalous higgs-V couplings

      subroutine width_hgamff(width)

      implicit none

      include "koppln_ew.inc"
      include "mssm.inc"

      real*8 width,dum
      real*8 Nc(4),gm(4),matrix1(4),matrix(11)
      real*8 Nf(4)
c     nu_e,e,nu_mu,mu,nu_tau,tau,u,d,c, s, b
c     1    2   3   4    5    6   7 8 9 10 11
      integer i,j,k,l,flavor(11)
       real*8 pi,q2min,q2max
      parameter(pi=3.141592653589793d0)
      logical ldoALL
      data flavor/1,2,1,2,1,2,3,4,3,4,4/
      data ldoALL/.false./

c
c     color factors
      Nc(1) = 1.d0
      Nc(2) = 1.d0
      Nc(3) = 3.d0*(1.d0 + ALFAS/PI)
      Nc(4) = 3.d0*(1.d0 + ALFAS/PI)

c     Number of fermion in generation
      Nf(1) = 3.0d0
      Nf(2) = 3.0d0
      Nf(3) = 2.0d0             !no top
      Nf(4) = 3.0d0

c     compute geometric means of masses
      gm(1) = 0d0  ! NEUTRINO MASS
      gm(2) = (Mf(2,1)*Mf(2,2)*Mf(2,3))**(1.0d0/3.d0)
      gm(3) = (Mf(3,1)*Mf(3,2))**(1.0d0/2.d0)
      gm(4) = (Mf(4,1)*Mf(4,2)*Mf(4,3))**(1.0d0/3.d0)

      dum =0.0d0
      if(ldoALL) then
          do k=1,11
             j = flavor(k)
             call quad1d(4,j,gm(j),xmh,matrix(k))
             dum = matrix(k)*Nc(j) + dum
          enddo
       else
          do j = 1,4
             call quad1d(4,j,gm(j),xmh,matrix1(j))
             dum = matrix1(j)*Nc(j)*Nf(j) + dum
          enddo
       endif
       
      width = dum

      return
      end



********************************************************************************
********************************************************************************

** This subroutine fills the Higgs-V-V couplings g5hvv for different 
** parametrisations

      subroutine anomH_convert(hvv1,hvv2,switch2,switch3,
     &     d, dodd, dB, dBo, dg1Z, dkgam, kgamo,
     &     fWW, fWWo, fBB, fBBo, fW, fB, fBo)

      implicit none

c     elektroweak parameters used here:
      double precision g, pi
      parameter(pi =3.141592653589793d0) 
      include "mssm.inc"

      double precision g5hvv(2,1:4,1:4), lambda5
      common /hcoupl/ g5hvv, lambda5

** Flags determining which anomalous HVV couplings are used.
**        = 0 : only HZA
**        = 1 : only HAA
**        = 2 : only HZZ
**        = 3 : only HWW
**        = 4 : all 
      integer hvv1, hvv2

** Switches for anomalous Higgs parameterisations
      logical switch2, switch3

** input for anomalous couplings
      double precision treefacW, treefacZ, loopfac
      common /lhcoup/ treefacW, treefacZ, loopfac
      double precision d, dodd, dB, dBo, dg1Z, dkgam, kgamo
      double precision fWW, fWWo, fBB, fBBo, fW, fB, fBo


      integer i, j, k


      g = Sqrt(4*pi*AlfaQED)/SW


c     fill g5hvv with d, dB, dg1Z, dkgam for L3 input
      if(switch2) then
c     
c        h-w-w anomalous couplings
         if((hvv1.eq.3).or.(hvv1.eq.4).or.(hvv1.eq.6)) then
            lambda5 = 480.0d0
            g5hvv(1,3,4) = lambda5*g/mW*(d+(MW**2)/(MZ**2)*dg1Z)
            g5hvv(1,4,3) = g5hvv(1,3,4)
            g5hvv(2,3,4) = lambda5*g/mW*dodd
            g5hvv(2,4,3) = g5hvv(2,3,4)
c        contribution to a1
            treefacW = treefacW-2*CW2*dg1Z
         endif
c
c        h-z-z anomalous couplings
         if((hvv1.eq.2).or.(hvv1.eq.4).or.(hvv1.eq.6).or.(hvv1.eq.7)) then
            lambda5 = 480.0d0
            g5hvv(1,2,2) = lambda5*g/mW*(CW2*d+SW2*dB+
     &           (1d0-2d0*SW2)*dg1Z+SW2/CW2*dkgam)
            g5hvv(2,2,2) = lambda5*g/mW*(CW2*
     &           dodd+SW2*dBo+SW2/CW2*kgamo)
c        contribution to a1
            treefacZ = treefacZ-2*((1d0-2d0*SW2)*dg1Z+SW2/CW2*dkgam) 
         endif
c
c        h-gamma-gamma
         if(hvv1.eq.1.or.hvv1.eq.4.or.hvv1.eq.5) then
c           SM + CPeven
            lambda5 = 480.0d0
            g5hvv(1,1,1) = lambda5*g/mW*(SW2*d+CW2*dB)
c           CP odd
            g5hvv(2,1,1) = lambda5*g/mW*(SW2*dodd+CW2*dBo)
         endif
c
c        h-z-gamma
         if(hvv1.eq.0.or.hvv1.eq.4.or.hvv1.eq.5.or.(hvv1.eq.7)) then
c           SM + CP-even
            lambda5 = 480.0d0
            g5hvv(1,1,2) = lambda5*g/(2d0*mW)*(2d0*CW*SW*
     &      (d-dB)+2d0*CW*SW*dg1Z-SW/CW*dkgam)
            g5hvv(1,2,1) = g5hvv(1,1,2)
c           CP-odd
            g5hvv(2,1,2) = lambda5*g/(2d0*mW)*(2d0*CW*SW*
     &      (dodd-dBo)-SW/CW*kgamo)
            g5hvv(2,2,1) = g5hvv(2,1,2)
         endif

c     fill ahvv with fWW, fBB, fW, fB:  converting from PLB 318, 515
      elseif(switch3) then
c
c        h-w-w anomalous couplings
         if(hvv2.eq.3.or.hvv2.eq.4) then
            lambda5 = 480.0d0
            g5hvv(1,3,4) = lambda5*g*mW*(-fWW+1d0/2d0*fW)
            g5hvv(1,4,3) = g5hvv(1,3,4)
            g5hvv(2,3,4) = -lambda5*g*mW*fWWo
            g5hvv(2,4,3) = g5hvv(2,3,4)
c        contribution to a1
            treefacW = treefacW - MW2*fW
         endif
c
c        h-z-z anomalous couplings
         if(hvv2.eq.2.or.hvv2.eq.4) then
            lambda5 = 480.0d0
            g5hvv(1,2,2) = lambda5*g*mW/CW2*(-SW2**2*fBB-
     &      CW2**2*fWW+1d0/2d0*(CW2*fW+SW2*fB))
            g5hvv(2,2,2) = lambda5*g*mW/CW2*(-SW2**2*
     &      fBBo-CW2**2*fWWo+SW2*fBo)
c        contribution to a1
            treefacZ = treefacZ-mZ**2*(CW2*fW+SW2*fB)
         endif
c
c        h-gamma-gamma
         if(hvv2.eq.1.or.hvv2.eq.4) then
            lambda5 = 480.0d0
c           SM + CPeven
            g5hvv(1,1,1) = -lambda5*g*mW*SW2*(fBB+fWW)
c           CP odd
            g5hvv(2,1,1) = -lambda5*g*mW*SW2*(fBBo+fWWo)
         endif
c
c        h-z-gamma
         if(hvv2.eq.0.or.hvv2.eq.4) then
            lambda5 = 480.0d0
c           SM + CP-even
            g5hvv(1,1,2) = lambda5*g*mW*SW/CW*(SW2*fBB-
     &           CW2*fWW+1d0/4d0*(fW-fB))
            g5hvv(1,2,1) = g5hvv(1,1,2)
c           CP-odd
            g5hvv(2,1,2) = lambda5*g*mW*SW/CW*(SW2*fBBo-
     &           CW2*fWWo-1d0/4d0*fBo)
            g5hvv(2,2,1) = g5hvv(2,1,2)
         endif
      endif      

c$$$** Debugging output
c$$$      do i = 1, 2
c$$$         do j = 1, 4
c$$$            do k = 1, 4
c$$$               write(*,*)'ijk g5hvv =', i, j, k, g5hvv(i,j,k)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      stop


      end 
