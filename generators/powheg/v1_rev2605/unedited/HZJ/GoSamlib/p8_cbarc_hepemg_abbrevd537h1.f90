module     p8_cbarc_hepemg_abbrevd537h1
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_globalsh1
   implicit none
   private
   complex(ki), dimension(29), public :: abb537
   complex(ki), public :: R2d537
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
      abb537(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb537(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb537(3)=NC**(-1)
      abb537(4)=spak2e6*spbk6e6
      abb537(5)=gCr*gel*abb537(3)*abb537(1)*abb537(2)*gHZZ*c1*TR*i_
      abb537(6)=2.0_ki*abb537(5)
      abb537(7)=abb537(4)*abb537(6)
      abb537(8)=spbk4k1*spak5k6
      abb537(9)=-es12*abb537(8)*abb537(7)
      abb537(10)=spak2k5*spbk2k1
      abb537(11)=abb537(10)*spbk6k4
      abb537(12)=spak2k6*abb537(11)
      abb537(13)=spbk4k1*spak2k5
      abb537(14)=-es12*abb537(13)
      abb537(12)=abb537(12)+abb537(14)
      abb537(14)=4.0_ki*abb537(5)
      abb537(12)=abb537(12)*abb537(14)
      abb537(15)=spak2e6*spbe6k4
      abb537(16)=abb537(15)*abb537(14)
      abb537(17)=-abb537(10)*abb537(16)
      abb537(18)=abb537(10)*abb537(15)
      abb537(19)=spbe6k1*spak1e6*abb537(13)
      abb537(18)=abb537(18)+abb537(19)
      abb537(18)=abb537(18)*abb537(14)
      abb537(5)=abb537(5)*abb537(13)
      abb537(5)=8.0_ki*abb537(5)
      abb537(13)=abb537(14)*spbk4k1
      abb537(19)=abb537(13)*spak5k6
      abb537(4)=abb537(19)*abb537(4)
      abb537(20)=-spbk2k1*spbe6k4
      abb537(21)=-spbe6k1*spbk4k2
      abb537(20)=abb537(20)+abb537(21)
      abb537(20)=spak2e6*abb537(20)
      abb537(21)=spae6k6*spbk6k4*spbe6k1
      abb537(20)=abb537(21)+abb537(20)
      abb537(20)=spak2k5*abb537(20)
      abb537(21)=spbk4k1*spak5e6
      abb537(22)=spbe6k1*spak1k2*abb537(21)
      abb537(20)=abb537(22)+abb537(20)
      abb537(20)=abb537(20)*abb537(6)
      abb537(22)=spak5k6*spbk6k1
      abb537(22)=abb537(22)+abb537(10)
      abb537(23)=-spak1e6*abb537(22)
      abb537(24)=-spak5e6*es12
      abb537(23)=abb537(24)+abb537(23)
      abb537(23)=spbk4k1*abb537(23)
      abb537(11)=-spae6k6*abb537(11)
      abb537(10)=spak2e6*spbk4k2*abb537(10)
      abb537(10)=abb537(23)+abb537(11)+abb537(10)
      abb537(10)=abb537(10)*abb537(6)
      abb537(11)=abb537(13)*spak5e6
      abb537(21)=-abb537(6)*abb537(21)
      abb537(23)=-spak2k6*spbk2k1*abb537(7)
      abb537(24)=abb537(14)*spbe6k1
      abb537(25)=spak2e6*abb537(24)
      abb537(26)=abb537(6)*spbe6k1
      abb537(27)=spak2e6*abb537(26)
      abb537(28)=-spae6k6*spbk6k1
      abb537(29)=2.0_ki*spak2e6
      abb537(29)=-spbk2k1*abb537(29)
      abb537(28)=abb537(28)+abb537(29)
      abb537(28)=abb537(28)*abb537(6)
      abb537(29)=abb537(14)*spbk6k4*spak2k6
      abb537(6)=-abb537(6)*abb537(15)
      abb537(8)=spak1e6*abb537(26)*abb537(8)
      abb537(15)=spae6k6*abb537(26)
      abb537(22)=-abb537(22)*abb537(14)
      abb537(24)=spak5e6*abb537(24)
      abb537(26)=spak5e6*abb537(26)
      abb537(13)=-spak1k5*abb537(13)
      R2d537=0.0_ki
      rat2 = rat2 + R2d537
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='537' value='", &
          & R2d537, "'/>"
      end if
   end subroutine
end module p8_cbarc_hepemg_abbrevd537h1
