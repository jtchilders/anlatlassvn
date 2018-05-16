module     p0_dbard_hepemg_abbrevd279h0
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_globalsh0
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
      abb279(5)=dotproduct(k3,spvak1e6)
      abb279(6)=dotproduct(k3,spvae6k2)
      abb279(7)=dotproduct(k3,spvae6k4)
      abb279(8)=dotproduct(k3,spvak5e6)
      abb279(9)=dotproduct(k3,spvak1k2)
      abb279(10)=dotproduct(k3,spvak1k4)
      abb279(11)=dotproduct(k3,spvak5k2)
      abb279(12)=1.0_ki/2.0_ki*spak5e6
      abb279(13)=-spbk5k4*abb279(12)
      abb279(13)=abb279(13)+1.0_ki/4.0_ki*abb279(7)
      abb279(14)=spbe6k2*spak1k5
      abb279(13)=abb279(14)*abb279(13)
      abb279(15)=1.0_ki/2.0_ki*spbk4k2
      abb279(16)=spbe6k4*spak1k5
      abb279(17)=-spak4e6*abb279(16)*abb279(15)
      abb279(18)=abb279(10)*spbe6k2*spak5e6
      abb279(19)=spak5e6*spbk4k2
      abb279(20)=1.0_ki/2.0_ki*abb279(19)
      abb279(21)=spak1k4*spbe6k4*abb279(20)
      abb279(13)=abb279(21)+abb279(17)-1.0_ki/4.0_ki*abb279(18)+abb279(13)
      abb279(17)=gHZZ*Nfrat*c1*i_*TR*gel*abb279(3)*abb279(2)*abb279(1)
      abb279(18)=abb279(17)*gBr
      abb279(17)=abb279(17)*gBl
      abb279(21)=abb279(18)+1.0_ki/3.0_ki*abb279(17)
      abb279(13)=abb279(21)*abb279(13)
      abb279(22)=spak4k5*abb279(15)
      abb279(22)=abb279(22)-1.0_ki/4.0_ki*abb279(11)
      abb279(23)=spak1e6*spbe6k4
      abb279(22)=abb279(23)*abb279(22)
      abb279(24)=spbk5k2*abb279(16)*abb279(12)
      abb279(25)=spak1e6*spbk4k2
      abb279(26)=abb279(8)*abb279(25)
      abb279(20)=-spbe6k5*spak1k5*abb279(20)
      abb279(20)=abb279(24)+abb279(20)+1.0_ki/4.0_ki*abb279(26)+abb279(22)
      abb279(22)=abb279(17)+1.0_ki/3.0_ki*abb279(18)
      abb279(20)=abb279(22)*abb279(20)
      abb279(24)=spak5k6*spbk6e6
      abb279(26)=abb279(22)*abb279(25)*abb279(24)
      abb279(27)=spbk6k4*spae6k6
      abb279(28)=abb279(14)*abb279(27)
      abb279(29)=abb279(21)*abb279(28)
      abb279(26)=abb279(26)+abb279(29)
      abb279(29)=abb279(18)-abb279(17)
      abb279(30)=abb279(29)*abb279(19)
      abb279(31)=spbk6e6*abb279(30)
      abb279(21)=-spbk6k4*spbe6k2*abb279(21)*abb279(12)
      abb279(32)=abb279(18)+abb279(17)
      abb279(33)=spak5e6*abb279(32)*spbe6k4
      abb279(34)=spbk6k2*abb279(33)
      abb279(21)=2.0_ki/3.0_ki*abb279(34)-1.0_ki/3.0_ki*abb279(31)+abb279(21)
      abb279(21)=spak1k6*abb279(21)
      abb279(31)=abb279(5)*abb279(30)
      abb279(34)=abb279(29)*abb279(16)
      abb279(35)=abb279(6)*abb279(34)
      abb279(31)=abb279(31)-abb279(35)
      abb279(35)=spae6k6*abb279(34)
      abb279(36)=1.0_ki/2.0_ki*spbe6k4
      abb279(22)=-abb279(22)*abb279(36)*spak5k6*spak1e6
      abb279(22)=1.0_ki/3.0_ki*abb279(35)+abb279(22)
      abb279(22)=spbk6k2*abb279(22)
      abb279(35)=abb279(32)*spbk4k2
      abb279(36)=abb279(35)*spak1k5
      abb279(37)=abb279(4)*abb279(36)
      abb279(38)=abb279(9)*abb279(33)
      abb279(13)=1.0_ki/3.0_ki*abb279(38)-2.0_ki/3.0_ki*abb279(37)+abb279(21)+1&
      &.0_ki/2.0_ki*abb279(26)+abb279(22)+abb279(20)+abb279(13)-1.0_ki/6.0_ki*a&
      &bb279(31)
      abb279(20)=spak1k2*spbk4k2
      abb279(21)=abb279(17)*abb279(20)*spbk6k2*spak5k6
      abb279(22)=spbk2k1*spak1k5
      abb279(26)=-abb279(18)*abb279(22)*spak1k6*spbk6k4
      abb279(21)=abb279(21)+abb279(26)
      abb279(21)=4.0_ki*abb279(21)
      abb279(26)=spbk2k1*spak1e6
      abb279(31)=abb279(16)*abb279(18)*abb279(26)
      abb279(37)=spak1k2*spbe6k2
      abb279(38)=-abb279(17)*abb279(19)*abb279(37)
      abb279(31)=abb279(31)+abb279(38)
      abb279(31)=4.0_ki*abb279(31)
      abb279(28)=abb279(32)*abb279(28)
      abb279(38)=abb279(17)+3.0_ki*abb279(18)
      abb279(39)=abb279(38)*abb279(24)
      abb279(40)=-abb279(25)*abb279(39)
      abb279(28)=abb279(40)-2.0_ki*abb279(28)
      abb279(40)=8.0_ki*abb279(36)
      abb279(41)=abb279(24)*spak1e6
      abb279(42)=abb279(35)*abb279(41)
      abb279(17)=abb279(18)+3.0_ki*abb279(17)
      abb279(18)=abb279(17)*abb279(27)
      abb279(43)=-abb279(14)*abb279(18)
      abb279(42)=-2.0_ki*abb279(42)+abb279(43)
      abb279(36)=-4.0_ki*abb279(36)
      abb279(25)=abb279(25)*abb279(29)
      abb279(43)=abb279(17)*abb279(15)
      abb279(44)=spak1e6*abb279(43)
      abb279(45)=abb279(14)*abb279(29)
      abb279(46)=1.0_ki/2.0_ki*abb279(38)
      abb279(14)=abb279(14)*abb279(46)
      abb279(24)=-spak1k2*abb279(29)*abb279(15)*abb279(24)
      abb279(47)=1.0_ki/2.0_ki*abb279(29)
      abb279(48)=-abb279(27)*abb279(47)*abb279(22)
      abb279(49)=spbk6e6*spak1e6
      abb279(50)=abb279(49)*abb279(22)*abb279(46)
      abb279(51)=abb279(29)*spak1k5
      abb279(52)=-spbk6k4*abb279(51)
      abb279(43)=-abb279(43)*abb279(37)*spae6k6
      abb279(53)=abb279(37)*abb279(29)
      abb279(54)=spbk6k2*spae6k6*abb279(53)
      abb279(55)=abb279(29)*spak1k6
      abb279(56)=spbk2k1*abb279(49)*abb279(55)
      abb279(54)=abb279(54)+abb279(56)
      abb279(54)=1.0_ki/2.0_ki*abb279(54)
      abb279(56)=spak1k6*spbk6e6*abb279(47)
      abb279(53)=-abb279(53)+abb279(56)
      abb279(26)=-abb279(29)*abb279(26)
      abb279(56)=abb279(47)*spae6k6
      abb279(57)=-spbk6k2*abb279(56)
      abb279(26)=abb279(26)+abb279(57)
      abb279(49)=abb279(49)*abb279(47)
      abb279(57)=1.0_ki/2.0_ki*abb279(18)
      abb279(37)=abb279(37)*abb279(57)
      abb279(20)=abb279(17)*abb279(20)
      abb279(55)=-spbk6k4*abb279(55)
      abb279(20)=2.0_ki*abb279(20)+abb279(55)
      abb279(55)=abb279(29)*abb279(23)
      abb279(17)=1.0_ki/2.0_ki*abb279(17)
      abb279(23)=-abb279(23)*abb279(17)
      abb279(58)=abb279(29)*spak5k6
      abb279(59)=spbk4k2*abb279(58)
      abb279(56)=-spbe6k2*abb279(56)
      abb279(46)=-spbk2k1*abb279(41)*abb279(46)
      abb279(22)=abb279(22)*abb279(38)
      abb279(58)=spbk6k2*abb279(58)
      abb279(22)=-2.0_ki*abb279(22)+abb279(58)
      abb279(58)=abb279(29)*spbe6k2
      abb279(60)=-spak5e6*abb279(58)
      abb279(38)=abb279(38)*abb279(12)
      abb279(61)=-spbe6k2*abb279(38)
      abb279(27)=spak2k5*spbe6k2*abb279(27)
      abb279(41)=spbk4k1*abb279(41)
      abb279(27)=abb279(27)+abb279(41)
      abb279(27)=abb279(32)*abb279(27)
      abb279(15)=spak2e6*abb279(15)*abb279(39)
      abb279(41)=abb279(57)*spbe6k1*spak1k5
      abb279(15)=abb279(41)+abb279(15)+abb279(27)
      abb279(27)=-spbk6k4*spak5k6
      abb279(41)=spbk4k1*spak1k5
      abb279(27)=abb279(41)+abb279(27)
      abb279(27)=abb279(32)*abb279(27)
      abb279(32)=spak2k5*abb279(35)
      abb279(27)=abb279(32)+abb279(27)
      abb279(27)=4.0_ki*abb279(27)
      abb279(32)=4.0_ki*abb279(33)
      abb279(33)=2.0_ki*abb279(33)
      abb279(35)=spbk4k1*spak1e6
      abb279(41)=spak2e6*spbk4k2
      abb279(35)=abb279(41)+abb279(35)
      abb279(35)=abb279(29)*abb279(35)
      abb279(18)=abb279(18)+abb279(35)
      abb279(18)=1.0_ki/2.0_ki*abb279(18)
      abb279(35)=-spak2k5*abb279(58)
      abb279(41)=-spbe6k1*abb279(51)
      abb279(35)=abb279(41)+abb279(39)+abb279(35)
      abb279(35)=1.0_ki/2.0_ki*abb279(35)
      abb279(39)=-spak2k5*spbe6k4*abb279(47)
      abb279(19)=abb279(19)*abb279(47)
      abb279(16)=-abb279(16)*abb279(47)
      abb279(12)=spbk4k1*abb279(29)*abb279(12)
      abb279(29)=spbk6e6*abb279(38)
      abb279(17)=spae6k6*spbe6k4*abb279(17)
      R2d279=abb279(13)
      rat2 = rat2 + R2d279
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='279' value='", &
          & R2d279, "'/>"
      end if
   end subroutine
end module p0_dbard_hepemg_abbrevd279h0
