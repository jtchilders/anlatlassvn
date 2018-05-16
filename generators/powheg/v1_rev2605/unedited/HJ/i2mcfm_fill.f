!---- A Collection of Filler routines to avoid file exchange if desired by host MC

!      subroutine i2mcfm_fill_higgs(mh,sp) 
!----- not needed with Madgraph input
!      implicit none 
!      include 'masses.f'
!      double precision mh
!      logical sp,spira
!      common/spira/spira
!      hmass=mh 
!     spira=sp
!      return 
!      end subroutine 
      

      subroutine i2mcfm_fill_alphas(alphas) 
      implicit none
      double precision alphas
      include 'MCFM_Include/qcdcouple.f'
      include 'pwhg_math.h'
      as=alphas
      ason2pi=alphas/(2*pi)
      ason4pi=alphas/(4*pi)
      gsq=4*pi*as
      end 


c      subroutine i2mcfm_set_alphas(as) 
c      implicit none 
c!---- This routine is a little different, will be called 
c!---- at run time and will set the value of inh_as
c      include 'MCFM_Include/interface_settings.f'
c      double precision as 
c      inh_as = as 
c      end 
      

      subroutine i2mcfm_fill_ew(ews_in,Gf_in,aemz_in,xw_in,nflav_in) 
      implicit none 
      include 'MCFM_Include/ewcouple.f' 
      include 'MCFM_Include/ewinput.f' 
      include 'MCFM_Include/nflav.f' 
      integer ews_in,nflav_in
      double precision Gf_in,aemz_in,xw_in
      ewscheme=ews_in
      Gf_inp=Gf_in
      aemmz_inp=aemz_in
      xw_inp=xw_in 
      nflav=nflav_in 
      return 
      end subroutine 

!      subroutine i2mcfm_fill_mv(wm_in,ww_in,zm_in,zw_in) 
!      implicit none 
!      include 'masses.f' 
!      include 'ewinput.f'
!      double precision wm_in,ww_in,zm_in,zw_in
     
!      wmass_inp=wm_in
!      wmass=wm_in
!      zmass=zm_in 
!      zmass_inp=zm_in
!      wwidth=ww_in
!      zwidth=zw_in 
    
!      return 
!      end subroutine 

!      subroutine i2mcfm_fill_mq(mqs)
!      implicit none 
!      include 'masses.f' 
!      double precision mqs(6) 
!      
!      md=mqs(1) 
!      mu=mqs(2) 
!      ms=mqs(3)
!      mc=mqs(4) 
!      mb=mqs(5) 
!      mt=mqs(6) 
!      mcsq=mc**2
!      mbsq=mb**2

!      return 
!      end subroutine 

     


!      subroutine i2mcfm_fill_ml(ml,wtau_in) 
!      implicit none 
!      include 'masses.f' 
!      double precision ml(3),wtau_in 

!      mel=ml(1)
!      mmu=ml(2) 
!      mtau=ml(3) 
!      tauwidth=wtau_in 

!      return 
!      end subroutine 


      subroutine i2mcfm_fill_ckm(Vckm)
      implicit none 
      double precision Vtd,Vts,Vtb
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      double precision Vckm(9) 
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      Vud=Vckm(1) 
      Vus=Vckm(2) 
      Vub=Vckm(3) 
      Vcd=Vckm(4) 
      Vcs=Vckm(5)
      Vcb=Vckm(6) 

      return 
      end subroutine 

      subroutine i2mcfm_fill_scheme(instring) 
      implicit none 
      include 'MCFM_Include/interface_settings.f'
      character *4 instring 
      
      inscheme=instring 
      return 
      end subroutine
