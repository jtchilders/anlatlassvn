module     p4_ubaru_hepemg_abbrevd541h1
   use p4_ubaru_hepemg_config, only: ki
   use p4_ubaru_hepemg_globalsh1
   implicit none
   private
   complex(ki), dimension(30), public :: abb541
   complex(ki), public :: R2d541
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p4_ubaru_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_color, only: TR
      use p4_ubaru_hepemg_globalsl1, only: epspow
      implicit none
      abb541(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb541(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb541(3)=spak1e6*spbe6k1
      abb541(4)=gHZZ*c1*NC*i_*TR*gUr*gel*abb541(2)*abb541(1)
      abb541(5)=abb541(4)*spbk4k1
      abb541(6)=abb541(5)*spak2k5
      abb541(7)=abb541(6)*abb541(3)
      abb541(8)=-es61+es345-es12
      abb541(9)=abb541(7)*abb541(8)
      abb541(10)=abb541(5)*spak5k6
      abb541(11)=spak2e6*abb541(10)*spbk6e6
      abb541(12)=2.0_ki*abb541(11)
      abb541(13)=es61*abb541(12)
      abb541(9)=abb541(13)+abb541(9)
      abb541(13)=-es12+es61
      abb541(13)=abb541(6)*abb541(13)
      abb541(14)=spbk6k4*spak2k6
      abb541(15)=abb541(4)*spbk6k1
      abb541(16)=spak5k6*abb541(15)
      abb541(17)=-abb541(16)*abb541(14)
      abb541(18)=spbk6k1*spak5k6
      abb541(19)=-spak1k2*abb541(5)*abb541(18)
      abb541(13)=abb541(19)+abb541(17)+abb541(13)
      abb541(13)=4.0_ki*abb541(13)
      abb541(17)=2.0_ki*abb541(4)
      abb541(18)=spak2e6*spbe6k4*abb541(18)*abb541(17)
      abb541(19)=abb541(7)+abb541(18)
      abb541(19)=2.0_ki*abb541(19)
      abb541(20)=8.0_ki*abb541(6)
      abb541(7)=abb541(18)-abb541(7)
      abb541(18)=abb541(4)*spbe6k1
      abb541(21)=abb541(18)*spak5e6
      abb541(22)=abb541(14)*abb541(21)
      abb541(22)=abb541(22)-abb541(7)
      abb541(22)=2.0_ki*abb541(22)
      abb541(23)=4.0_ki*abb541(6)
      abb541(24)=2.0_ki*abb541(21)
      abb541(25)=-abb541(14)*abb541(24)
      abb541(7)=-abb541(11)+abb541(25)+abb541(7)
      abb541(6)=2.0_ki*abb541(6)
      abb541(11)=spae6k6*spbk6k4
      abb541(25)=-spbk4k2*spak2e6
      abb541(11)=abb541(25)+abb541(11)
      abb541(11)=abb541(16)*abb541(11)
      abb541(25)=abb541(5)*spak5e6
      abb541(26)=es61*abb541(25)
      abb541(11)=abb541(26)+abb541(11)
      abb541(11)=2.0_ki*abb541(11)
      abb541(26)=-2.0_ki*abb541(25)
      abb541(27)=abb541(18)*spak2e6
      abb541(8)=abb541(27)*abb541(8)
      abb541(15)=-8.0_ki*spak2k6*abb541(15)
      abb541(28)=2.0_ki*abb541(27)
      abb541(17)=spae6k6*spbk6k1*abb541(17)
      abb541(29)=4.0_ki*abb541(4)
      abb541(14)=abb541(14)*abb541(29)
      abb541(3)=abb541(10)*abb541(3)
      abb541(30)=spbk4k2*abb541(21)*spak2k6
      abb541(3)=2.0_ki*abb541(3)+abb541(30)
      abb541(10)=4.0_ki*abb541(10)
      abb541(18)=2.0_ki*abb541(18)
      abb541(18)=spae6k6*abb541(18)
      abb541(4)=-spbk2k1*abb541(4)*spak2k5
      abb541(4)=-abb541(16)+abb541(4)
      abb541(4)=4.0_ki*abb541(4)
      abb541(5)=-4.0_ki*spak1k5*abb541(5)
      R2d541=0.0_ki
      rat2 = rat2 + R2d541
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='541' value='", &
          & R2d541, "'/>"
      end if
   end subroutine
end module p4_ubaru_hepemg_abbrevd541h1
