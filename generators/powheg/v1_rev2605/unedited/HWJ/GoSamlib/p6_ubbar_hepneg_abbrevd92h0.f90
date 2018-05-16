module     p6_ubbar_hepneg_abbrevd92h0
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(19), public :: abb92
   complex(ki), public :: R2d92
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
      abb92(1)=1.0_ki/(-es61-es12+es345)
      abb92(2)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb92(3)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb92(4)=NC**(-1)
      abb92(5)=gHWW*c1*gW**2*i_*TR*CVBU*abb92(3)*abb92(2)*abb92(1)
      abb92(6)=abb92(4)*abb92(5)*spbk4k1
      abb92(7)=abb92(5)*NC
      abb92(8)=spbk4k1*abb92(7)
      abb92(9)=abb92(6)-abb92(8)
      abb92(10)=spak2e6*spbk6e6
      abb92(11)=abb92(10)*spak5k6
      abb92(12)=abb92(9)*abb92(11)
      abb92(13)=-es61+es345-es12
      abb92(14)=2.0_ki*abb92(11)*abb92(8)*abb92(13)
      abb92(8)=abb92(6)+abb92(8)
      abb92(15)=abb92(8)*spak2k5
      abb92(13)=4.0_ki*abb92(13)*abb92(15)
      abb92(9)=spbk6e6*abb92(9)*spak5e6*spak2k6
      abb92(16)=-4.0_ki*abb92(9)
      abb92(17)=-16.0_ki*abb92(15)
      abb92(11)=-abb92(6)*abb92(11)
      abb92(9)=abb92(11)+abb92(9)
      abb92(9)=4.0_ki*abb92(9)
      abb92(11)=8.0_ki*abb92(15)
      abb92(15)=2.0_ki*abb92(12)
      abb92(6)=2.0_ki*spbk6k2*spae6k6*spak2k5*abb92(6)
      abb92(7)=spak2k6*abb92(10)*abb92(7)
      abb92(7)=3.0_ki*abb92(7)
      abb92(10)=spbk2k1*abb92(7)
      abb92(18)=abb92(4)+NC
      abb92(19)=8.0_ki*spak2k6
      abb92(5)=abb92(19)*abb92(5)*abb92(18)
      abb92(18)=spbk6k1*abb92(5)
      abb92(7)=spbk4k2*abb92(7)
      abb92(5)=-spbk6k4*abb92(5)
      abb92(8)=4.0_ki*abb92(8)
      abb92(8)=-spak5k6*abb92(8)
      R2d92=-abb92(12)
      rat2 = rat2 + R2d92
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='92' value='", &
          & R2d92, "'/>"
      end if
   end subroutine
end module p6_ubbar_hepneg_abbrevd92h0
