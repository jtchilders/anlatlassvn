module     p0_dbard_hepemg_abbrevd613h2
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_globalsh2
   implicit none
   private
   complex(ki), dimension(22), public :: abb613
   complex(ki), public :: R2d613
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p0_dbard_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_color, only: TR
      use p0_dbard_hepemg_globalsl1, only: epspow
      implicit none
      abb613(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb613(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb613(3)=NC**(-1)
      abb613(4)=es61**(-1)
      abb613(5)=spak1e6*abb613(3)
      abb613(6)=NC*spak1e6
      abb613(7)=abb613(6)-abb613(5)
      abb613(8)=spbk6e6*spak4k6
      abb613(9)=spak1k4*spbe6k1
      abb613(10)=abb613(8)-abb613(9)
      abb613(7)=abb613(10)*abb613(7)
      abb613(10)=gDl*ger*abb613(1)*abb613(2)*gHZZ*c1*TR*i_
      abb613(11)=abb613(10)*abb613(4)
      abb613(12)=spbk5k2*abb613(11)
      abb613(13)=-abb613(7)*abb613(12)
      abb613(14)=-2.0_ki*abb613(8)+abb613(9)
      abb613(14)=abb613(14)*abb613(6)*abb613(10)*spbk5k2
      abb613(15)=NC*spak1k4
      abb613(16)=abb613(3)*spak1k4
      abb613(15)=abb613(15)+abb613(16)
      abb613(10)=4.0_ki*abb613(10)
      abb613(10)=-spbk5k2*abb613(15)*abb613(10)
      abb613(17)=abb613(9)*spak1e6
      abb613(18)=spbk6e6*spak1k6
      abb613(19)=abb613(18)*spak4e6
      abb613(17)=abb613(17)+abb613(19)
      abb613(17)=abb613(17)*NC
      abb613(19)=abb613(19)*abb613(3)
      abb613(17)=abb613(17)-abb613(19)
      abb613(9)=-abb613(9)*abb613(5)
      abb613(9)=abb613(9)-abb613(17)
      abb613(19)=4.0_ki*abb613(12)
      abb613(9)=abb613(9)*abb613(19)
      abb613(15)=abb613(15)*abb613(12)
      abb613(20)=16.0_ki*abb613(15)
      abb613(5)=abb613(8)*abb613(5)
      abb613(5)=abb613(5)+abb613(17)
      abb613(5)=abb613(5)*abb613(19)
      abb613(8)=-8.0_ki*abb613(15)
      abb613(12)=2.0_ki*abb613(12)
      abb613(7)=abb613(7)*abb613(12)
      abb613(12)=-spae6k6*spbk6k1*abb613(16)*abb613(12)
      abb613(6)=abb613(11)*abb613(6)
      abb613(15)=spak1k6*spbk6k2
      abb613(16)=spbe6k1*abb613(15)
      abb613(17)=3.0_ki*abb613(18)
      abb613(18)=spbk2k1*abb613(17)
      abb613(16)=abb613(16)+abb613(18)
      abb613(16)=abb613(16)*abb613(6)
      abb613(11)=8.0_ki*abb613(11)
      abb613(18)=NC+abb613(3)
      abb613(11)=abb613(11)*abb613(18)
      abb613(15)=-abb613(15)*abb613(11)
      abb613(21)=spak1k6*spbk6k5
      abb613(22)=-spbe6k1*abb613(21)
      abb613(17)=-spbk5k1*abb613(17)
      abb613(17)=abb613(22)+abb613(17)
      abb613(6)=abb613(17)*abb613(6)
      abb613(11)=abb613(21)*abb613(11)
      abb613(17)=abb613(19)*spak4k6*abb613(18)
      R2d613=abb613(13)
      rat2 = rat2 + R2d613
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='613' value='", &
          & R2d613, "'/>"
      end if
   end subroutine
end module p0_dbard_hepemg_abbrevd613h2
