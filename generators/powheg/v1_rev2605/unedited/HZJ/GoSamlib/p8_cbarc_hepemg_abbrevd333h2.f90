module     p8_cbarc_hepemg_abbrevd333h2
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_globalsh2
   implicit none
   private
   complex(ki), dimension(21), public :: abb333
   complex(ki), public :: R2d333
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
      abb333(1)=1.0_ki/(-es61-es12+es345)
      abb333(2)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb333(3)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb333(4)=NC**(-1)
      abb333(5)=spae6k6*spak1k4
      abb333(6)=abb333(5)*spbk6k5
      abb333(7)=NC-abb333(4)
      abb333(7)=abb333(7)*gCl*gHZZ*c1*i_*TR*ger*abb333(3)*abb333(2)*abb333(1)
      abb333(8)=spbe6k2*abb333(7)
      abb333(9)=abb333(6)*abb333(8)
      abb333(10)=spak1k2*spbk5k2
      abb333(5)=-spbk6k1*abb333(10)*abb333(5)
      abb333(6)=es61*abb333(6)
      abb333(5)=abb333(6)+abb333(5)
      abb333(5)=2.0_ki*abb333(8)*abb333(5)
      abb333(6)=-es345+es12+es61
      abb333(11)=spbk5k2*spak1k4
      abb333(6)=abb333(11)*abb333(6)
      abb333(12)=spbk5k1*spak1k4
      abb333(13)=spbk6k2*spak1k6
      abb333(14)=-abb333(12)*abb333(13)
      abb333(6)=abb333(14)+abb333(6)
      abb333(6)=-abb333(7)*abb333(6)
      abb333(14)=spak1k6*spbk6k5
      abb333(10)=abb333(14)+abb333(10)
      abb333(14)=-spak4k6*abb333(10)*abb333(7)
      abb333(15)=-spbk6k2*abb333(14)
      abb333(6)=abb333(15)+abb333(6)
      abb333(6)=4.0_ki*abb333(6)
      abb333(10)=abb333(10)*abb333(8)
      abb333(15)=abb333(10)*spak4e6
      abb333(16)=abb333(8)*spak1e6
      abb333(17)=abb333(16)*abb333(12)
      abb333(15)=abb333(15)+abb333(17)
      abb333(15)=4.0_ki*abb333(15)
      abb333(17)=4.0_ki*abb333(9)
      abb333(18)=-2.0_ki*abb333(9)
      abb333(19)=abb333(7)*spbk6e6
      abb333(20)=-spak1k6*abb333(19)
      abb333(21)=abb333(8)*spak1k2
      abb333(20)=abb333(20)-abb333(21)
      abb333(12)=abb333(20)*abb333(12)
      abb333(14)=spbk6e6*abb333(14)
      abb333(10)=spak2k4*abb333(10)
      abb333(11)=abb333(11)*abb333(19)
      abb333(19)=-abb333(8)*spbk6k5*spak1k4
      abb333(11)=abb333(19)+abb333(11)
      abb333(11)=spak2k6*abb333(11)
      abb333(10)=abb333(11)+abb333(10)+abb333(12)+abb333(14)
      abb333(10)=2.0_ki*abb333(10)
      abb333(11)=2.0_ki*spae6k6
      abb333(12)=abb333(11)*spbk6k2
      abb333(14)=abb333(12)*abb333(21)
      abb333(7)=-8.0_ki*abb333(7)*abb333(13)
      abb333(13)=8.0_ki*abb333(16)
      abb333(16)=-4.0_ki*abb333(20)
      abb333(11)=abb333(8)*abb333(11)
      abb333(8)=spak2k4*abb333(12)*abb333(8)
      R2d333=abb333(9)
      rat2 = rat2 + R2d333
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='333' value='", &
          & R2d333, "'/>"
      end if
   end subroutine
end module p8_cbarc_hepemg_abbrevd333h2
