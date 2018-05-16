module     p0_dbard_hepemg_abbrevd279h3
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_globalsh3
   implicit none
   private
   complex(ki), dimension(61), public :: abb279
   complex(ki), public :: R2d279
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
      abb279(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb279(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb279(3)=es12**(-1)
      abb279(4)=dotproduct(k3,e6)
      abb279(5)=dotproduct(k3,spvae6k1)
      abb279(6)=dotproduct(k3,spvak2e6)
      abb279(7)=dotproduct(k3,spvak4e6)
      abb279(8)=dotproduct(k3,spvae6k5)
      abb279(9)=dotproduct(k3,spvak2k1)
      abb279(10)=dotproduct(k3,spvak2k5)
      abb279(11)=dotproduct(k3,spvak4k1)
      abb279(12)=1.0_ki/2.0_ki*spak2k4
      abb279(13)=spbk5k4*abb279(12)
      abb279(13)=abb279(13)-1.0_ki/4.0_ki*abb279(10)
      abb279(14)=spbe6k1*spak4e6
      abb279(13)=abb279(14)*abb279(13)
      abb279(15)=1.0_ki/2.0_ki*spbe6k5
      abb279(16)=spak4e6*spbk5k1
      abb279(17)=spak2k5*abb279(16)*abb279(15)
      abb279(18)=spbe6k1*spak2k4
      abb279(19)=abb279(8)*abb279(18)
      abb279(20)=spbe6k5*spak2k4
      abb279(21)=1.0_ki/2.0_ki*abb279(20)
      abb279(22)=-spak5e6*spbk5k1*abb279(21)
      abb279(13)=abb279(17)+abb279(22)+1.0_ki/4.0_ki*abb279(19)+abb279(13)
      abb279(17)=gHZZ*Nfrat*c1*i_*TR*ger*abb279(3)*abb279(2)*abb279(1)
      abb279(19)=abb279(17)*gBr
      abb279(17)=abb279(17)*gBl
      abb279(22)=abb279(19)+1.0_ki/3.0_ki*abb279(17)
      abb279(13)=abb279(22)*abb279(13)
      abb279(23)=-spak4k5*abb279(15)
      abb279(23)=abb279(23)+1.0_ki/4.0_ki*abb279(7)
      abb279(24)=spak2e6*spbk5k1
      abb279(23)=abb279(24)*abb279(23)
      abb279(25)=-spbe6k4*abb279(16)*abb279(12)
      abb279(26)=abb279(11)*spak2e6*spbe6k5
      abb279(21)=spbk4k1*spak4e6*abb279(21)
      abb279(21)=abb279(21)+abb279(25)-1.0_ki/4.0_ki*abb279(26)+abb279(23)
      abb279(23)=abb279(17)+1.0_ki/3.0_ki*abb279(19)
      abb279(21)=abb279(23)*abb279(21)
      abb279(25)=spbk6k5*spae6k6
      abb279(26)=abb279(22)*abb279(18)*abb279(25)
      abb279(27)=spak4k6*spbk6e6
      abb279(28)=abb279(24)*abb279(27)
      abb279(29)=abb279(23)*abb279(28)
      abb279(26)=abb279(26)+abb279(29)
      abb279(29)=abb279(19)-abb279(17)
      abb279(30)=abb279(29)*abb279(20)
      abb279(31)=spae6k6*abb279(30)
      abb279(23)=-spak4k6*abb279(23)*abb279(15)*spak2e6
      abb279(32)=abb279(19)+abb279(17)
      abb279(33)=spbe6k5*abb279(32)*spak4e6
      abb279(34)=spak2k6*abb279(33)
      abb279(23)=2.0_ki/3.0_ki*abb279(34)+1.0_ki/3.0_ki*abb279(31)+abb279(23)
      abb279(23)=spbk6k1*abb279(23)
      abb279(31)=abb279(5)*abb279(30)
      abb279(34)=abb279(29)*abb279(16)
      abb279(35)=abb279(6)*abb279(34)
      abb279(31)=abb279(31)-abb279(35)
      abb279(35)=spbk6e6*abb279(34)
      abb279(36)=1.0_ki/2.0_ki*spak4e6
      abb279(22)=-abb279(22)*abb279(36)*spbk6k5*spbe6k1
      abb279(22)=-1.0_ki/3.0_ki*abb279(35)+abb279(22)
      abb279(22)=spak2k6*abb279(22)
      abb279(35)=abb279(32)*spak2k4
      abb279(36)=abb279(35)*spbk5k1
      abb279(37)=abb279(4)*abb279(36)
      abb279(38)=abb279(9)*abb279(33)
      abb279(13)=1.0_ki/3.0_ki*abb279(38)-2.0_ki/3.0_ki*abb279(37)+abb279(23)+1&
      &.0_ki/2.0_ki*abb279(26)+abb279(22)+abb279(21)+abb279(13)+1.0_ki/6.0_ki*a&
      &bb279(31)
      abb279(21)=spbk2k1*spak2k4
      abb279(22)=abb279(19)*abb279(21)*spak2k6*spbk6k5
      abb279(23)=spak1k2*spbk5k1
      abb279(26)=-abb279(17)*abb279(23)*spbk6k1*spak4k6
      abb279(22)=abb279(22)+abb279(26)
      abb279(22)=4.0_ki*abb279(22)
      abb279(26)=spbk2k1*spak2e6
      abb279(31)=-abb279(19)*abb279(20)*abb279(26)
      abb279(37)=spak1k2*spbe6k1
      abb279(38)=abb279(16)*abb279(17)*abb279(37)
      abb279(31)=abb279(31)+abb279(38)
      abb279(31)=4.0_ki*abb279(31)
      abb279(28)=abb279(32)*abb279(28)
      abb279(38)=abb279(19)+3.0_ki*abb279(17)
      abb279(39)=abb279(38)*abb279(25)
      abb279(40)=-abb279(18)*abb279(39)
      abb279(28)=abb279(40)-2.0_ki*abb279(28)
      abb279(40)=8.0_ki*abb279(36)
      abb279(41)=abb279(25)*spbe6k1
      abb279(42)=abb279(35)*abb279(41)
      abb279(17)=abb279(17)+3.0_ki*abb279(19)
      abb279(19)=abb279(17)*abb279(27)
      abb279(43)=-abb279(24)*abb279(19)
      abb279(42)=-2.0_ki*abb279(42)+abb279(43)
      abb279(36)=-4.0_ki*abb279(36)
      abb279(18)=abb279(18)*abb279(29)
      abb279(43)=abb279(17)*abb279(12)
      abb279(44)=spbe6k1*abb279(43)
      abb279(45)=abb279(24)*abb279(29)
      abb279(46)=1.0_ki/2.0_ki*abb279(38)
      abb279(24)=abb279(24)*abb279(46)
      abb279(25)=spbk2k1*abb279(29)*abb279(12)*abb279(25)
      abb279(47)=1.0_ki/2.0_ki*abb279(29)
      abb279(48)=abb279(27)*abb279(47)*abb279(23)
      abb279(43)=-abb279(43)*abb279(26)*spbk6e6
      abb279(49)=abb279(29)*spbk6k5
      abb279(50)=-spak2k4*abb279(49)
      abb279(51)=spae6k6*spbe6k1
      abb279(52)=abb279(51)*abb279(23)*abb279(46)
      abb279(53)=abb279(26)*abb279(29)
      abb279(54)=-spak2k6*spbk6e6*abb279(53)
      abb279(55)=abb279(29)*spbk6k1
      abb279(56)=-spak1k2*abb279(51)*abb279(55)
      abb279(54)=abb279(54)+abb279(56)
      abb279(54)=1.0_ki/2.0_ki*abb279(54)
      abb279(56)=-spbk6k1*spae6k6*abb279(47)
      abb279(53)=abb279(53)+abb279(56)
      abb279(37)=abb279(29)*abb279(37)
      abb279(56)=abb279(47)*spbk6e6
      abb279(57)=spak2k6*abb279(56)
      abb279(37)=abb279(37)+abb279(57)
      abb279(56)=spak2e6*abb279(56)
      abb279(46)=-spak1k2*abb279(41)*abb279(46)
      abb279(23)=abb279(23)*abb279(38)
      abb279(49)=-spak2k6*abb279(49)
      abb279(23)=-2.0_ki*abb279(23)+abb279(49)
      abb279(49)=abb279(29)*spak2e6
      abb279(57)=spbe6k5*abb279(49)
      abb279(38)=abb279(38)*abb279(15)
      abb279(58)=-spak2e6*abb279(38)
      abb279(59)=abb279(29)*spbk5k1
      abb279(60)=spak4k6*abb279(59)
      abb279(51)=-abb279(51)*abb279(47)
      abb279(61)=1.0_ki/2.0_ki*abb279(19)
      abb279(26)=abb279(26)*abb279(61)
      abb279(21)=abb279(17)*abb279(21)
      abb279(55)=spak4k6*abb279(55)
      abb279(21)=2.0_ki*abb279(21)+abb279(55)
      abb279(55)=-abb279(29)*abb279(14)
      abb279(17)=1.0_ki/2.0_ki*abb279(17)
      abb279(14)=-abb279(14)*abb279(17)
      abb279(27)=spbk5k2*spak2e6*abb279(27)
      abb279(41)=spak1k4*abb279(41)
      abb279(27)=abb279(27)+abb279(41)
      abb279(27)=abb279(32)*abb279(27)
      abb279(12)=spbe6k2*abb279(12)*abb279(39)
      abb279(41)=abb279(61)*spak1e6*spbk5k1
      abb279(12)=abb279(41)+abb279(12)+abb279(27)
      abb279(27)=-spak4k6*spbk6k5
      abb279(41)=spak1k4*spbk5k1
      abb279(27)=abb279(41)+abb279(27)
      abb279(27)=abb279(32)*abb279(27)
      abb279(32)=spbk5k2*abb279(35)
      abb279(27)=abb279(32)+abb279(27)
      abb279(27)=4.0_ki*abb279(27)
      abb279(32)=4.0_ki*abb279(33)
      abb279(33)=2.0_ki*abb279(33)
      abb279(35)=-spak1k4*spbe6k1
      abb279(41)=-spbe6k2*spak2k4
      abb279(35)=abb279(41)+abb279(35)
      abb279(35)=abb279(29)*abb279(35)
      abb279(19)=abb279(19)+abb279(35)
      abb279(19)=1.0_ki/2.0_ki*abb279(19)
      abb279(35)=spbk5k2*abb279(49)
      abb279(41)=spak1e6*abb279(59)
      abb279(35)=abb279(41)+abb279(39)+abb279(35)
      abb279(35)=1.0_ki/2.0_ki*abb279(35)
      abb279(20)=-abb279(20)*abb279(47)
      abb279(39)=spbk5k2*spak4e6*abb279(47)
      abb279(15)=-spak1k4*abb279(29)*abb279(15)
      abb279(16)=abb279(16)*abb279(47)
      abb279(17)=spbk6e6*spak4e6*abb279(17)
      abb279(29)=spae6k6*abb279(38)
      R2d279=abb279(13)
      rat2 = rat2 + R2d279
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='279' value='", &
          & R2d279, "'/>"
      end if
   end subroutine
end module p0_dbard_hepemg_abbrevd279h3
