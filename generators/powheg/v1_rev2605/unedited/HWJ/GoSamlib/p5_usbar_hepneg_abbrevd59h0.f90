module     p5_usbar_hepneg_abbrevd59h0
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(28), public :: abb59
   complex(ki), public :: R2d59
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p5_usbar_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_color, only: TR
      use p5_usbar_hepneg_globalsl1, only: epspow
      implicit none
      abb59(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb59(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb59(3)=NC**(-1)
      abb59(4)=spbk6k4*spae6k6
      abb59(5)=spbk4k1*spak1e6
      abb59(4)=abb59(4)-abb59(5)
      abb59(5)=CVSU*abb59(3)*gW**2*abb59(1)*abb59(2)*gHWW*c1*TR*i_
      abb59(6)=abb59(5)*spbe6k1
      abb59(6)=2.0_ki*abb59(6)
      abb59(7)=abb59(6)*abb59(4)*spak2k5*es12
      abb59(8)=spak5k6*spbk6k1
      abb59(9)=spbk4k1*spak1k2
      abb59(10)=abb59(9)*abb59(8)
      abb59(11)=spak2k5*spbk4k1
      abb59(12)=es12*abb59(11)
      abb59(10)=abb59(10)+abb59(12)
      abb59(12)=4.0_ki*abb59(5)
      abb59(10)=abb59(10)*abb59(12)
      abb59(13)=abb59(12)*spbe6k1
      abb59(14)=spak5e6*abb59(13)
      abb59(15)=abb59(14)*abb59(9)
      abb59(16)=abb59(12)*spak2k5
      abb59(4)=-spbe6k1*abb59(4)*abb59(16)
      abb59(11)=abb59(5)*abb59(11)
      abb59(11)=8.0_ki*abb59(11)
      abb59(17)=spak2e6*spak1k5
      abb59(18)=-spak1k2*spak5e6
      abb59(17)=abb59(17)+abb59(18)
      abb59(17)=spbe6k1*abb59(17)
      abb59(18)=spbk6e6*spak5k6
      abb59(19)=-spak2e6*abb59(18)
      abb59(17)=abb59(17)+abb59(19)
      abb59(17)=spbk4k1*abb59(17)
      abb59(19)=spak2k5*spbe6k4
      abb59(20)=spak2e6*spbk2k1*abb59(19)
      abb59(17)=abb59(20)+abb59(17)
      abb59(5)=2.0_ki*abb59(5)
      abb59(17)=abb59(17)*abb59(5)
      abb59(20)=spbe6k1*spak1k5
      abb59(18)=abb59(20)-abb59(18)
      abb59(18)=abb59(9)*abb59(18)
      abb59(20)=spbk6k4*spak2k6
      abb59(9)=abb59(20)-abb59(9)
      abb59(20)=spbe6k2*abb59(9)
      abb59(21)=spbe6k4*es12
      abb59(20)=abb59(21)+abb59(20)
      abb59(20)=spak2k5*abb59(20)
      abb59(18)=abb59(20)+abb59(18)
      abb59(18)=abb59(18)*abb59(5)
      abb59(20)=abb59(16)*spbe6k4
      abb59(19)=abb59(5)*abb59(19)
      abb59(21)=-spbk6k4*abb59(16)
      abb59(22)=abb59(6)*spae6k6
      abb59(23)=-spak1k2*spbk6k1*abb59(22)
      abb59(13)=-spak2e6*abb59(13)
      abb59(24)=-spak2e6*abb59(6)
      abb59(25)=spbk6e6*spak2k6
      abb59(26)=2.0_ki*spbe6k1
      abb59(26)=-spak1k2*abb59(26)
      abb59(25)=abb59(25)+abb59(26)
      abb59(25)=abb59(25)*abb59(5)
      abb59(5)=abb59(5)*spak2e6
      abb59(26)=-spbk6e6*abb59(5)
      abb59(9)=abb59(9)*abb59(12)
      abb59(27)=spbe6k4*abb59(12)*spak2e6
      abb59(5)=-spbe6k4*abb59(5)
      abb59(8)=-abb59(12)*abb59(8)
      abb59(28)=spak5e6*abb59(6)
      abb59(16)=spbk4k2*abb59(16)
      abb59(6)=spak1e6*abb59(6)
      R2d59=0.0_ki
      rat2 = rat2 + R2d59
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='59' value='", &
          & R2d59, "'/>"
      end if
   end subroutine
end module p5_usbar_hepneg_abbrevd59h0
