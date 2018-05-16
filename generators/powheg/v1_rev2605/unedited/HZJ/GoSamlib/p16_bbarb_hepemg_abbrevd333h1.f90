module     p16_bbarb_hepemg_abbrevd333h1
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_globalsh1
   implicit none
   private
   complex(ki), dimension(21), public :: abb333
   complex(ki), public :: R2d333
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p16_bbarb_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_color, only: TR
      use p16_bbarb_hepemg_globalsl1, only: epspow
      implicit none
      abb333(1)=1.0_ki/(-es61-es12+es345)
      abb333(2)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb333(3)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb333(4)=NC**(-1)
      abb333(5)=spbk6e6*spbk4k1
      abb333(6)=abb333(5)*spak5k6
      abb333(7)=NC-abb333(4)
      abb333(7)=abb333(7)*gBr*gHZZ*c1*i_*TR*gel*abb333(3)*abb333(2)*abb333(1)
      abb333(8)=spak2e6*abb333(7)
      abb333(9)=abb333(6)*abb333(8)
      abb333(10)=spbk2k1*spak2k5
      abb333(5)=-spak1k6*abb333(10)*abb333(5)
      abb333(6)=es61*abb333(6)
      abb333(5)=abb333(6)+abb333(5)
      abb333(5)=2.0_ki*abb333(8)*abb333(5)
      abb333(6)=-es345+es12+es61
      abb333(11)=spak2k5*spbk4k1
      abb333(6)=abb333(11)*abb333(6)
      abb333(12)=spak1k5*spbk4k1
      abb333(13)=spak2k6*spbk6k1
      abb333(14)=-abb333(12)*abb333(13)
      abb333(6)=abb333(14)+abb333(6)
      abb333(6)=-abb333(7)*abb333(6)
      abb333(14)=spbk6k1*spak5k6
      abb333(10)=abb333(14)+abb333(10)
      abb333(14)=-spbk6k4*abb333(10)*abb333(7)
      abb333(15)=-spak2k6*abb333(14)
      abb333(6)=abb333(15)+abb333(6)
      abb333(6)=4.0_ki*abb333(6)
      abb333(10)=abb333(10)*abb333(8)
      abb333(15)=abb333(10)*spbe6k4
      abb333(16)=abb333(8)*spbe6k1
      abb333(17)=abb333(16)*abb333(12)
      abb333(15)=abb333(15)+abb333(17)
      abb333(15)=4.0_ki*abb333(15)
      abb333(17)=4.0_ki*abb333(9)
      abb333(18)=-2.0_ki*abb333(9)
      abb333(19)=abb333(7)*spae6k6
      abb333(20)=-spbk6k1*abb333(19)
      abb333(21)=abb333(8)*spbk2k1
      abb333(20)=abb333(20)-abb333(21)
      abb333(12)=abb333(20)*abb333(12)
      abb333(14)=spae6k6*abb333(14)
      abb333(10)=spbk4k2*abb333(10)
      abb333(11)=abb333(11)*abb333(19)
      abb333(19)=-abb333(8)*spak5k6*spbk4k1
      abb333(11)=abb333(19)+abb333(11)
      abb333(11)=spbk6k2*abb333(11)
      abb333(10)=abb333(11)+abb333(10)+abb333(12)+abb333(14)
      abb333(10)=2.0_ki*abb333(10)
      abb333(11)=2.0_ki*spbk6e6
      abb333(12)=abb333(11)*spak2k6
      abb333(14)=abb333(12)*abb333(21)
      abb333(7)=-8.0_ki*abb333(7)*abb333(13)
      abb333(13)=8.0_ki*abb333(16)
      abb333(16)=-4.0_ki*abb333(20)
      abb333(11)=abb333(8)*abb333(11)
      abb333(8)=spbk4k2*abb333(12)*abb333(8)
      R2d333=abb333(9)
      rat2 = rat2 + R2d333
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='333' value='", &
          & R2d333, "'/>"
      end if
   end subroutine
end module p16_bbarb_hepemg_abbrevd333h1
