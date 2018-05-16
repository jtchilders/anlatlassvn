module     p16_bbarb_hepemg_abbrevd553h3
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_globalsh3
   implicit none
   private
   complex(ki), dimension(37), public :: abb553
   complex(ki), public :: R2d553
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
      abb553(1)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb553(2)=1.0_ki/(mH**2+es45+es61-es23-es345)
      abb553(3)=1.0_ki/(mH**2+es45-es123+es12-es345)
      abb553(4)=sqrt(mT**2)
      abb553(5)=1.0_ki/(es45-es123-es61+es23)
      abb553(6)=spak2l3**(-1)
      abb553(7)=spbl3k2**(-1)
      abb553(8)=spbk5k1*abb553(5)
      abb553(9)=i_*TR*c1*gHT*gBr*ger*abb553(3)*abb553(1)
      abb553(10)=abb553(4)*abb553(9)
      abb553(11)=abb553(10)*abb553(8)
      abb553(12)=abb553(11)*spak1k4
      abb553(10)=abb553(10)*abb553(2)
      abb553(13)=spak2k4*abb553(10)
      abb553(14)=abb553(13)*spbk5k2
      abb553(12)=abb553(14)-abb553(12)
      abb553(14)=-spak2e6*abb553(12)
      abb553(13)=abb553(13)*spbk5k4
      abb553(15)=abb553(13)*spak4e6
      abb553(16)=abb553(14)+abb553(15)
      abb553(17)=-spbe6k1*abb553(16)
      abb553(11)=abb553(11)*spak4k5
      abb553(18)=abb553(11)*spbe6k5
      abb553(19)=abb553(18)*spak2e6
      abb553(17)=abb553(17)-abb553(19)
      abb553(19)=4.0_ki*abb553(17)
      abb553(20)=spbe6k1*abb553(9)*spak2k4*abb553(2)
      abb553(21)=spbk5k2*spak2e6
      abb553(22)=abb553(21)*abb553(20)
      abb553(8)=abb553(9)*abb553(8)
      abb553(9)=-spak1k4*abb553(8)*spbe6k1*spak2e6
      abb553(9)=abb553(22)+abb553(9)
      abb553(22)=2.0_ki*spak2e6
      abb553(8)=-spbe6k5*spak4k5*abb553(8)*abb553(22)
      abb553(20)=abb553(20)*spak4e6*spbk5k4
      abb553(8)=2.0_ki*abb553(9)+abb553(8)-2.0_ki*abb553(20)
      abb553(8)=abb553(8)*abb553(4)**3
      abb553(9)=abb553(6)*abb553(7)*mH**2
      abb553(20)=abb553(21)*abb553(9)
      abb553(21)=spbk5l3*spal3e6
      abb553(20)=abb553(21)+abb553(20)
      abb553(21)=spak2k6*spbk6e6
      abb553(23)=abb553(21)*abb553(11)
      abb553(20)=abb553(23)*abb553(20)
      abb553(10)=spak2k4**2*abb553(9)*spbk5k4*abb553(10)
      abb553(24)=abb553(10)*spbe6k2
      abb553(25)=spbe6l3*spak2l3
      abb553(26)=abb553(25)*abb553(12)
      abb553(27)=abb553(13)*spal3k4
      abb553(28)=abb553(27)*spbe6l3
      abb553(24)=abb553(28)+abb553(24)+abb553(26)
      abb553(26)=spbk6k1*spae6k6
      abb553(28)=abb553(26)*abb553(24)
      abb553(29)=spal3e6*spbl3k1
      abb553(30)=abb553(29)*abb553(12)
      abb553(31)=abb553(21)*abb553(30)
      abb553(32)=abb553(13)*spak4k6
      abb553(33)=abb553(32)*spbk6e6
      abb553(34)=-abb553(29)*abb553(33)
      abb553(14)=-abb553(14)*abb553(21)
      abb553(35)=-spak2e6*abb553(33)
      abb553(14)=abb553(14)+abb553(35)
      abb553(9)=abb553(9)*spbk2k1
      abb553(14)=abb553(14)*abb553(9)
      abb553(35)=abb553(11)*spbk6k5
      abb553(36)=abb553(35)*spae6k6
      abb553(37)=-abb553(25)*abb553(36)
      abb553(8)=abb553(37)+abb553(14)+abb553(34)+abb553(31)+abb553(28)+abb553(2&
      &0)+abb553(8)
      abb553(8)=4.0_ki*abb553(8)
      abb553(14)=-spak2k6*abb553(12)
      abb553(14)=abb553(14)+abb553(32)
      abb553(14)=spbk6k1*abb553(14)
      abb553(20)=spak2k6*abb553(35)
      abb553(14)=abb553(20)+abb553(14)
      abb553(14)=8.0_ki*abb553(14)
      abb553(20)=16.0_ki*abb553(17)
      abb553(17)=-8.0_ki*abb553(17)
      abb553(28)=-spbe6k1*abb553(12)
      abb553(18)=abb553(28)+abb553(18)
      abb553(28)=spak2k6*abb553(18)
      abb553(31)=spbe6k1*abb553(32)
      abb553(28)=abb553(31)+abb553(28)
      abb553(28)=2.0_ki*abb553(28)
      abb553(31)=spbk6k1*abb553(16)
      abb553(32)=spak2e6*abb553(35)
      abb553(31)=abb553(32)+abb553(31)
      abb553(31)=2.0_ki*abb553(31)
      abb553(25)=abb553(11)*abb553(25)
      abb553(23)=abb553(23)+abb553(25)
      abb553(23)=4.0_ki*abb553(23)
      abb553(25)=spak2e6*abb553(9)
      abb553(25)=abb553(25)+abb553(26)+abb553(29)
      abb553(25)=4.0_ki*abb553(13)*abb553(25)
      abb553(29)=-spak2l3*abb553(18)
      abb553(27)=spbe6k1*abb553(27)
      abb553(27)=abb553(27)+abb553(29)
      abb553(27)=4.0_ki*abb553(27)
      abb553(29)=-spbl3k1*abb553(16)
      abb553(32)=abb553(11)*spbk5l3
      abb553(34)=spak2e6*abb553(32)
      abb553(29)=abb553(34)+abb553(29)
      abb553(29)=4.0_ki*abb553(29)
      abb553(10)=4.0_ki*spbe6k1*abb553(10)
      abb553(26)=-abb553(12)*abb553(26)
      abb553(9)=-abb553(15)*abb553(9)
      abb553(15)=-spal3e6*abb553(32)
      abb553(9)=abb553(15)+abb553(36)+abb553(9)+abb553(26)-abb553(30)
      abb553(9)=4.0_ki*abb553(9)
      abb553(15)=-abb553(12)*abb553(21)
      abb553(15)=abb553(33)+abb553(15)-abb553(24)
      abb553(15)=4.0_ki*abb553(15)
      abb553(21)=-abb553(22)*abb553(11)*spbk6e6
      abb553(16)=-2.0_ki*spbk6e6*abb553(16)
      abb553(22)=2.0_ki*spbe6k1
      abb553(22)=-spae6k6*abb553(13)*abb553(22)
      abb553(13)=-32.0_ki*abb553(13)
      abb553(18)=-2.0_ki*spae6k6*abb553(18)
      abb553(11)=-32.0_ki*abb553(11)
      abb553(12)=32.0_ki*abb553(12)
      R2d553=abb553(19)
      rat2 = rat2 + R2d553
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='553' value='", &
          & R2d553, "'/>"
      end if
   end subroutine
end module p16_bbarb_hepemg_abbrevd553h3
