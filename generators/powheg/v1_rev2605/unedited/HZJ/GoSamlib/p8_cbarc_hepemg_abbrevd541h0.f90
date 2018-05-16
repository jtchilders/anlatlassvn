module     p8_cbarc_hepemg_abbrevd541h0
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_globalsh0
   implicit none
   private
   complex(ki), dimension(29), public :: abb541
   complex(ki), public :: R2d541
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p8_cbarc_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_color, only: TR
      use p8_cbarc_hepemg_globalsl1, only: epspow
      implicit none
      abb541(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb541(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb541(3)=spbe6k1*spak1e6
      abb541(4)=gHZZ*c1*NC*i_*TR*gCl*gel*abb541(2)*abb541(1)
      abb541(5)=abb541(4)*spak1k5
      abb541(6)=abb541(5)*spbk4k2
      abb541(7)=abb541(6)*abb541(3)
      abb541(8)=-es61+es345-es12
      abb541(9)=abb541(7)*abb541(8)
      abb541(10)=abb541(5)*spbk6k4
      abb541(11)=spbe6k2*abb541(10)*spae6k6
      abb541(12)=2.0_ki*abb541(11)
      abb541(13)=es61*abb541(12)
      abb541(9)=abb541(13)+abb541(9)
      abb541(13)=-es12+es61
      abb541(13)=abb541(6)*abb541(13)
      abb541(14)=spak5k6*spbk6k2
      abb541(15)=abb541(4)*spak1k6
      abb541(16)=spbk6k4*abb541(15)
      abb541(17)=-abb541(16)*abb541(14)
      abb541(18)=spak1k6*spbk6k4
      abb541(19)=-spbk2k1*abb541(5)*abb541(18)
      abb541(13)=abb541(19)+abb541(17)+abb541(13)
      abb541(13)=4.0_ki*abb541(13)
      abb541(17)=2.0_ki*abb541(4)
      abb541(18)=spbe6k2*spak5e6*abb541(18)*abb541(17)
      abb541(19)=abb541(7)+abb541(18)
      abb541(19)=2.0_ki*abb541(19)
      abb541(20)=8.0_ki*abb541(6)
      abb541(7)=abb541(18)-abb541(7)
      abb541(18)=abb541(4)*spak1e6
      abb541(21)=abb541(18)*spbe6k4
      abb541(22)=abb541(14)*abb541(21)
      abb541(22)=abb541(22)-abb541(7)
      abb541(22)=2.0_ki*abb541(22)
      abb541(23)=4.0_ki*abb541(6)
      abb541(24)=2.0_ki*abb541(21)
      abb541(25)=-abb541(14)*abb541(24)
      abb541(7)=-abb541(11)+abb541(25)+abb541(7)
      abb541(6)=2.0_ki*abb541(6)
      abb541(11)=spbk6e6*spak5k6
      abb541(25)=-spak2k5*spbe6k2
      abb541(11)=abb541(25)+abb541(11)
      abb541(11)=abb541(16)*abb541(11)
      abb541(25)=abb541(5)*spbe6k4
      abb541(26)=es61*abb541(25)
      abb541(11)=abb541(26)+abb541(11)
      abb541(11)=2.0_ki*abb541(11)
      abb541(26)=-2.0_ki*abb541(25)
      abb541(3)=abb541(10)*abb541(3)
      abb541(27)=spak2k5*abb541(21)*spbk6k2
      abb541(3)=2.0_ki*abb541(3)+abb541(27)
      abb541(10)=4.0_ki*abb541(10)
      abb541(27)=abb541(18)*spbe6k2
      abb541(8)=abb541(27)*abb541(8)
      abb541(15)=-8.0_ki*spbk6k2*abb541(15)
      abb541(28)=2.0_ki*abb541(27)
      abb541(17)=spbk6e6*spak1k6*abb541(17)
      abb541(18)=2.0_ki*abb541(18)
      abb541(18)=spbk6e6*abb541(18)
      abb541(29)=-spak1k2*abb541(4)*spbk4k2
      abb541(16)=-abb541(16)+abb541(29)
      abb541(16)=4.0_ki*abb541(16)
      abb541(4)=4.0_ki*abb541(4)
      abb541(14)=abb541(14)*abb541(4)
      abb541(5)=-4.0_ki*spbk4k1*abb541(5)
      R2d541=0.0_ki
      rat2 = rat2 + R2d541
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='541' value='", &
          & R2d541, "'/>"
      end if
   end subroutine
end module p8_cbarc_hepemg_abbrevd541h0
