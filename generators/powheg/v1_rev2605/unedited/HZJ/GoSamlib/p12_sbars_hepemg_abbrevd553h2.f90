module     p12_sbars_hepemg_abbrevd553h2
   use p12_sbars_hepemg_config, only: ki
   use p12_sbars_hepemg_globalsh2
   implicit none
   private
   complex(ki), dimension(37), public :: abb553
   complex(ki), public :: R2d553
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p12_sbars_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_color, only: TR
      use p12_sbars_hepemg_globalsl1, only: epspow
      implicit none
      abb553(1)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb553(2)=1.0_ki/(mH**2+es45+es61-es23-es345)
      abb553(3)=1.0_ki/(mH**2+es45-es123+es12-es345)
      abb553(4)=sqrt(mT**2)
      abb553(5)=1.0_ki/(es45-es123-es61+es23)
      abb553(6)=spak2l3**(-1)
      abb553(7)=spbl3k2**(-1)
      abb553(8)=spak1k4*abb553(5)
      abb553(9)=i_*TR*c1*gHT*gSl*ger*abb553(3)*abb553(1)
      abb553(10)=abb553(4)*abb553(9)
      abb553(11)=abb553(10)*abb553(8)
      abb553(12)=abb553(11)*spbk5k1
      abb553(10)=abb553(10)*abb553(2)
      abb553(13)=spbk5k2*abb553(10)
      abb553(14)=abb553(13)*spak2k4
      abb553(12)=abb553(14)-abb553(12)
      abb553(14)=-spbe6k2*abb553(12)
      abb553(13)=abb553(13)*spak4k5
      abb553(15)=abb553(13)*spbe6k5
      abb553(16)=abb553(14)-abb553(15)
      abb553(17)=-spak1e6*abb553(16)
      abb553(11)=abb553(11)*spbk5k4
      abb553(18)=abb553(11)*spak4e6
      abb553(19)=abb553(18)*spbe6k2
      abb553(17)=abb553(17)+abb553(19)
      abb553(19)=4.0_ki*abb553(17)
      abb553(20)=spak1e6*abb553(9)*spbk5k2*abb553(2)
      abb553(21)=spak2k4*spbe6k2
      abb553(22)=abb553(21)*abb553(20)
      abb553(8)=abb553(9)*abb553(8)
      abb553(9)=-spbk5k1*abb553(8)*spak1e6*spbe6k2
      abb553(9)=abb553(22)+abb553(9)
      abb553(22)=2.0_ki*spbe6k2
      abb553(8)=spak4e6*spbk5k4*abb553(8)*abb553(22)
      abb553(20)=abb553(20)*spbe6k5*spak4k5
      abb553(8)=2.0_ki*abb553(9)+abb553(8)+2.0_ki*abb553(20)
      abb553(8)=abb553(8)*abb553(4)**3
      abb553(9)=abb553(6)*abb553(7)*mH**2
      abb553(20)=-abb553(21)*abb553(9)
      abb553(21)=-spal3k4*spbe6l3
      abb553(20)=abb553(21)+abb553(20)
      abb553(21)=spbk6k2*spae6k6
      abb553(23)=abb553(21)*abb553(11)
      abb553(20)=abb553(23)*abb553(20)
      abb553(10)=spbk5k2**2*abb553(9)*spak4k5*abb553(10)
      abb553(24)=abb553(10)*spak2e6
      abb553(25)=spal3e6*spbl3k2
      abb553(26)=abb553(25)*abb553(12)
      abb553(27)=abb553(13)*spbk5l3
      abb553(28)=abb553(27)*spal3e6
      abb553(24)=abb553(28)+abb553(24)-abb553(26)
      abb553(26)=spak1k6*spbk6e6
      abb553(28)=-abb553(26)*abb553(24)
      abb553(29)=spbe6l3*spak1l3
      abb553(30)=abb553(29)*abb553(12)
      abb553(31)=abb553(21)*abb553(30)
      abb553(32)=abb553(13)*spbk6k5
      abb553(33)=abb553(32)*spae6k6
      abb553(34)=abb553(29)*abb553(33)
      abb553(14)=-abb553(14)*abb553(21)
      abb553(35)=spbe6k2*abb553(33)
      abb553(14)=abb553(14)+abb553(35)
      abb553(9)=abb553(9)*spak1k2
      abb553(14)=abb553(14)*abb553(9)
      abb553(35)=abb553(11)*spak4k6
      abb553(36)=abb553(35)*spbk6e6
      abb553(37)=abb553(25)*abb553(36)
      abb553(8)=abb553(37)+abb553(14)+abb553(34)+abb553(31)+abb553(28)+abb553(2&
      &0)+abb553(8)
      abb553(8)=4.0_ki*abb553(8)
      abb553(14)=-spbk6k2*abb553(12)
      abb553(14)=abb553(14)-abb553(32)
      abb553(14)=spak1k6*abb553(14)
      abb553(20)=-spbk6k2*abb553(35)
      abb553(14)=abb553(20)+abb553(14)
      abb553(14)=8.0_ki*abb553(14)
      abb553(20)=16.0_ki*abb553(17)
      abb553(17)=-8.0_ki*abb553(17)
      abb553(28)=spak1k6*abb553(16)
      abb553(31)=-spbe6k2*abb553(35)
      abb553(28)=abb553(31)+abb553(28)
      abb553(28)=2.0_ki*abb553(28)
      abb553(31)=-spak1e6*abb553(12)
      abb553(18)=abb553(31)-abb553(18)
      abb553(31)=spbk6k2*abb553(18)
      abb553(32)=-spak1e6*abb553(32)
      abb553(31)=abb553(32)+abb553(31)
      abb553(31)=2.0_ki*abb553(31)
      abb553(32)=-spbe6k2*abb553(9)
      abb553(29)=abb553(32)-abb553(26)-abb553(29)
      abb553(29)=4.0_ki*abb553(13)*abb553(29)
      abb553(25)=-abb553(11)*abb553(25)
      abb553(23)=-abb553(23)+abb553(25)
      abb553(23)=4.0_ki*abb553(23)
      abb553(25)=-spak1l3*abb553(16)
      abb553(32)=abb553(11)*spal3k4
      abb553(34)=-spbe6k2*abb553(32)
      abb553(25)=abb553(34)+abb553(25)
      abb553(25)=4.0_ki*abb553(25)
      abb553(34)=-spbl3k2*abb553(18)
      abb553(27)=-spak1e6*abb553(27)
      abb553(27)=abb553(27)+abb553(34)
      abb553(27)=4.0_ki*abb553(27)
      abb553(26)=-abb553(12)*abb553(26)
      abb553(9)=abb553(15)*abb553(9)
      abb553(15)=spbe6l3*abb553(32)
      abb553(9)=abb553(15)-abb553(36)+abb553(9)+abb553(26)-abb553(30)
      abb553(9)=4.0_ki*abb553(9)
      abb553(10)=-4.0_ki*spak1e6*abb553(10)
      abb553(15)=-abb553(12)*abb553(21)
      abb553(15)=-abb553(33)+abb553(15)+abb553(24)
      abb553(15)=4.0_ki*abb553(15)
      abb553(21)=2.0_ki*spak1e6
      abb553(21)=spbk6e6*abb553(13)*abb553(21)
      abb553(18)=-2.0_ki*spbk6e6*abb553(18)
      abb553(22)=abb553(22)*abb553(11)*spae6k6
      abb553(11)=32.0_ki*abb553(11)
      abb553(16)=-2.0_ki*spae6k6*abb553(16)
      abb553(13)=32.0_ki*abb553(13)
      abb553(12)=32.0_ki*abb553(12)
      R2d553=abb553(19)
      rat2 = rat2 + R2d553
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='553' value='", &
          & R2d553, "'/>"
      end if
   end subroutine
end module p12_sbars_hepemg_abbrevd553h2
