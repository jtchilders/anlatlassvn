

!---- This serves as a driver/wrapper for the MCFM pieces of H+j in POWHEG 
!----- C. Williams August 2011 

      subroutine i2MCFM_2_POWHEG(p,vflav,res) 
      implicit none 
      include 'nlegborn.h'
      include 'MCFM_Include/constants.f' 
      include 'MCFM_Include/scale.f'
      include 'MCFM_Include/interface_settings.f' 
      include 'pwhg_st.h' 
      include 'PhysPars.h' 
      integer vflav(nlegborn) 
      double precision p(0:3,nlegborn) 
      double precision res(0:3)
      double precision pmcfm(mxpart,4)

      call i2mcfm_powheg_momswap(p,pmcfm)
!----- set st_alpha and scale 
      call i2mcfm_fill_alphas(st_alpha)
      scale=dsqrt(st_muren2) !-- CHECK This
      musq=st_muren2

      ret1=vflav(1) 
      ret2=vflav(2) 
      ret3=vflav(4) 
      
  !    call pgen(pin,pmcfm)
      call gg_hg_vi(pmcfm,res)

      return 
      end subroutine 


      subroutine i2MCFM_2_POWHEG_IP
      implicit none 
      include 'nlegborn.h'
      include 'MCFM_Include/constants.f' 
      include 'pwhg_st.h'
      include 'PhysPars.h' 
        
!------ setup interface settings as needed for POWHEG 
      call set_interface_MCFM
!------ Setup common blocks (ew parameters) 
!----- Need to sort this out! 
c      ph_alphaem = 0.0949563/(4*pi)
c      ph_gfermi=0.116639D-04
c      ph_sthw2=0.2226459d0

c     in the following, set ewscheme = 3, USER choice, i.e. the couplings are set exactly at
c     the user's values
      call i2mcfm_fill_ew(3,ph_gfermi,ph_alphaem,ph_sthw2,5)
!------- Setup initial value of alpha_s 
c      call i2mcfm_set_alphas(st_alpha) 

      return 
      end subroutine 


      subroutine i2mcfm_powheg_momswap(p,pmcfm) 
      implicit none 
      include 'nlegborn.h' 
      include 'MCFM_Include/constants.f'
      integer nu,i
      double precision p(0:3,nlegborn),pmcfm(mxpart,4),swap 
!---- Initalize pmcfm 
      do i=1,mxpart 
         do nu=1,4 
            pmcfm(i,nu)=0d0 
         enddo
      enddo
!---- swap ingoing/outgoing 
      do i=1,2 
         do nu=0,3 
            pmcfm(i,nu+1)=-p(nu,i) 
         enddo 
      enddo 

!------ MCFM needs jet = p5 so set p3 = 0 
      do i=4,nlegborn+1 
         do nu=0,3
            pmcfm(i,nu+1)=p(nu,i-1) 
         enddo 
      enddo


!------- POWHEG notation 0 = E ,1 =x 2=y 3=z 
!------- MCFM notation   1=x,2=y,3=z,4=E 

!----- swap POWHEG => MCFM 
      do i=1,nlegborn+1 
         swap=pmcfm(i,1)
         pmcfm(i,1)=pmcfm(i,2)
         pmcfm(i,2)=pmcfm(i,3)
         pmcfm(i,3)=pmcfm(i,4) 
         pmcfm(i,4)=swap
      enddo
!---- Check that MCFM momentum conservation is implmented properly, i.e. p1+p2+p3+p4+p5=0 

!      write(6,*) 'PMCFM' 
!      do i=1,5 
!         do nu=1,4 
!            write(6,*) 'p(',i,nu,') = ',pmcfm(i,nu) 
!         enddo
!      enddo

      
!      write(6,*) 'POWHEG' 
!      do i=1,nlegborn 
!         do nu=0,3 
!            write(6,*) 'p(',nu,i,') = ',p(nu,i) 
!         enddo
!      enddo

      do i=1,4 
         if(pmcfm(1,i)+pmcfm(2,i)+pmcfm(3,i)+pmcfm(4,i)+pmcfm(5,i)
     &        .gt.1d-10) then 
            write(6,*) 'Momentum problem in i2mcfm_powheg_momswap' 
      
         endif
      enddo

      return 
      end subroutine 
      
      
