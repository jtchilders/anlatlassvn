module     p6_ubbar_hepneg_abbrevd34h0
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(14), public :: abb34
   complex(ki), public :: R2d34
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p6_ubbar_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p6_ubbar_hepneg_kinematics
      use p6_ubbar_hepneg_model
      use p6_ubbar_hepneg_color, only: TR
      use p6_ubbar_hepneg_globalsl1, only: epspow
      implicit none
      abb34(1)=1.0_ki/(-es61-es12+es345)
      abb34(2)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb34(3)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb34(4)=NC**(-1)
      abb34(5)=abb34(4)-NC
      abb34(5)=abb34(5)*gHWW*spbk6e6*spak2e6*c1*gW**2*i_*TR*CVBU*abb34(3)*abb34&
      &(2)*abb34(1)
      abb34(6)=abb34(5)*spak5k6
      abb34(7)=spbk4k1*abb34(6)
      abb34(8)=abb34(5)*spbk4k1
      abb34(9)=abb34(8)*spak1k6
      abb34(10)=spbk2k1*spak2k5*abb34(9)
      abb34(11)=-es61*abb34(7)
      abb34(10)=abb34(10)+abb34(11)
      abb34(10)=2.0_ki*abb34(10)
      abb34(11)=-4.0_ki*abb34(7)
      abb34(12)=2.0_ki*abb34(7)
      abb34(8)=-2.0_ki*abb34(8)*spak2k5
      abb34(6)=2.0_ki*abb34(6)
      abb34(13)=-spbk6k4*abb34(6)
      abb34(14)=2.0_ki*abb34(5)
      abb34(5)=spak2k6*spbk4k2*abb34(5)
      abb34(5)=-abb34(9)+abb34(5)
      abb34(5)=2.0_ki*abb34(5)
      abb34(6)=-spbk4k2*abb34(6)
      R2d34=-abb34(7)
      rat2 = rat2 + R2d34
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='34' value='", &
          & R2d34, "'/>"
      end if
   end subroutine
end module p6_ubbar_hepneg_abbrevd34h0
