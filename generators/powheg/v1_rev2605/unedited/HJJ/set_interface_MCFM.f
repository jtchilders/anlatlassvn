      subroutine set_interface_MCFM
      implicit none 
      include 'nlegborn.h' 
      include 'MCFM_Include/interface_settings.f' 
c      integer vflav(nlegborn) 
      real * 8 powheginput
      external powheginput
!---- PROCESS NUMBER SHOULD BE CHANGED TO FIT NEED 
!     203 '  f(p1)+f(p2) --> H(-->b(p3)+b~(p4)) + f(p5)' 'N'
!     204 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4)) + f(p5)' 'N'
!     272 '  f(p1)+f(p2) --> H(tau^-(p3)+tau^+(p4))+f(p5)+f(p6) [in heavy top limit]' 'N'
!     pid NEVER USED
      pid = 272
!-------- In POWHEG Higgs does not decay, H=p3+f(p4) 
!----- Set up vflav 
c      ret1=vflav(1) 
c      ret2=vflav(2) 
c      ret3=vflav(4) 
c      ret4=vflav(5)
!------inscheme 
      inscheme='dred' 
!------ Number of particles 
      n_particles=nlegborn
!------ Inherit alpha_s from host 
      inherit_as =.true. 
!------ return poles in epsilon
      if(powheginput("#polecheck").gt.0) then
         ret_poles=.true. 
      else
         ret_poles=.false.
      endif
!------ efficient (should always = true) 
      efficient=.true. 
!----- alphas_eq1,alphaew_eq1
      alphas_eq1=.false. 
      alphaew_eq1=.false. 
      
      return 
      end subroutine
