module     p8_cbarc_hepemg_abbrevd537h0
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_globalsh0
   implicit none
   private
   complex(ki), dimension(30), public :: abb537
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
      abb537(4)=spbe6k2*spae6k6
      abb537(5)=gCl*gel*abb537(3)*abb537(1)*abb537(2)*gHZZ*c1*TR*i_
      abb537(6)=2.0_ki*abb537(5)
      abb537(7)=abb537(4)*abb537(6)
      abb537(8)=spak1k5*spbk6k4
      abb537(9)=-es12*abb537(8)*abb537(7)
      abb537(10)=spbk4k2*spak1k2
      abb537(11)=abb537(10)*spak5k6
      abb537(12)=spbk6k2*abb537(11)
      abb537(13)=spak1k5*spbk4k2
      abb537(14)=-es12*abb537(13)
      abb537(12)=abb537(12)+abb537(14)
      abb537(14)=4.0_ki*abb537(5)
      abb537(12)=abb537(12)*abb537(14)
      abb537(15)=spbe6k2*spak5e6
      abb537(16)=abb537(15)*abb537(14)
      abb537(17)=-abb537(10)*abb537(16)
      abb537(18)=abb537(10)*abb537(15)
      abb537(19)=spak1e6*spbe6k1*abb537(13)
      abb537(18)=abb537(18)+abb537(19)
      abb537(18)=abb537(18)*abb537(14)
      abb537(5)=abb537(5)*abb537(13)
      abb537(5)=8.0_ki*abb537(5)
      abb537(13)=abb537(14)*spak1k5
      abb537(19)=abb537(13)*spbk6k4
      abb537(4)=abb537(19)*abb537(4)
      abb537(20)=-spak1k2*spak5e6
      abb537(21)=-spak1e6*spak2k5
      abb537(20)=abb537(20)+abb537(21)
      abb537(20)=spbe6k2*abb537(20)
      abb537(21)=spbk6e6*spak5k6*spak1e6
      abb537(20)=abb537(21)+abb537(20)
      abb537(20)=spbk4k2*abb537(20)
      abb537(21)=spak1k5*spbe6k4
      abb537(22)=spak1e6*spbk2k1*abb537(21)
      abb537(20)=abb537(22)+abb537(20)
      abb537(20)=abb537(20)*abb537(6)
      abb537(22)=spbk6k4*spak1k6
      abb537(22)=abb537(22)+abb537(10)
      abb537(23)=-spbe6k1*abb537(22)
      abb537(24)=-spbe6k4*es12
      abb537(23)=abb537(24)+abb537(23)
      abb537(23)=spak1k5*abb537(23)
      abb537(11)=-spbk6e6*abb537(11)
      abb537(10)=spbe6k2*spak2k5*abb537(10)
      abb537(10)=abb537(23)+abb537(11)+abb537(10)
      abb537(10)=abb537(10)*abb537(6)
      abb537(11)=abb537(13)*spbe6k4
      abb537(21)=-abb537(6)*abb537(21)
      abb537(23)=abb537(6)*spak1e6
      abb537(8)=spbe6k1*abb537(23)*abb537(8)
      abb537(24)=-spbk6k2*spak1k2*abb537(7)
      abb537(25)=abb537(14)*spak1e6
      abb537(26)=spbe6k2*abb537(25)
      abb537(27)=spbe6k2*abb537(23)
      abb537(28)=-spbk6e6*spak1k6
      abb537(29)=2.0_ki*spbe6k2
      abb537(29)=-spak1k2*abb537(29)
      abb537(28)=abb537(28)+abb537(29)
      abb537(28)=abb537(28)*abb537(6)
      abb537(29)=spbk6e6*abb537(23)
      abb537(22)=-abb537(22)*abb537(14)
      abb537(25)=spbe6k4*abb537(25)
      abb537(23)=spbe6k4*abb537(23)
      abb537(30)=abb537(14)*spak5k6*spbk6k2
      abb537(6)=-abb537(6)*abb537(15)
      abb537(13)=-spbk4k1*abb537(13)
      R2d537=0.0_ki
      rat2 = rat2 + R2d537
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='537' value='", &
          & R2d537, "'/>"
      end if
   end subroutine
end module p8_cbarc_hepemg_abbrevd537h0
