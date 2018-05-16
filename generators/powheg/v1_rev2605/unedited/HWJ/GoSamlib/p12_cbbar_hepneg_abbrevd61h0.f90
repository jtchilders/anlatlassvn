module     p12_cbbar_hepneg_abbrevd61h0
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(29), public :: abb61
   complex(ki), public :: R2d61
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p12_cbbar_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_color, only: TR
      use p12_cbbar_hepneg_globalsl1, only: epspow
      implicit none
      abb61(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb61(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb61(3)=spak1e6*spbe6k1
      abb61(4)=c1*NC*gW**2*i_*TR*gHWW*CVBC*abb61(2)*abb61(1)
      abb61(5)=abb61(4)*spak2k5
      abb61(6)=abb61(5)*spbk4k1
      abb61(7)=abb61(6)*abb61(3)
      abb61(8)=es345-es12
      abb61(9)=abb61(7)*abb61(8)
      abb61(10)=abb61(4)*spbk4k1
      abb61(11)=abb61(10)*spak5k6
      abb61(12)=spak2e6*abb61(11)*spbk6e6
      abb61(13)=2.0_ki*abb61(12)
      abb61(14)=abb61(13)-abb61(7)
      abb61(14)=es61*abb61(14)
      abb61(9)=abb61(14)+abb61(9)
      abb61(14)=-3.0_ki*es12+es345-es61
      abb61(14)=abb61(6)*abb61(14)
      abb61(15)=spak1k2*spbk6k1*abb61(11)
      abb61(14)=-2.0_ki*abb61(15)+abb61(14)
      abb61(14)=2.0_ki*abb61(14)
      abb61(13)=-abb61(13)-abb61(7)
      abb61(13)=2.0_ki*abb61(13)
      abb61(15)=4.0_ki*abb61(6)
      abb61(16)=2.0_ki*abb61(7)
      abb61(17)=2.0_ki*abb61(4)
      abb61(18)=abb61(17)*spak2e6
      abb61(19)=spbk6k1*spbe6k4
      abb61(20)=abb61(19)*spak5k6*abb61(18)
      abb61(17)=abb61(17)*spak2k6
      abb61(21)=spbe6k1*abb61(17)*spbk6k4
      abb61(22)=spak5e6*abb61(21)
      abb61(12)=-abb61(12)+abb61(20)-abb61(22)
      abb61(20)=2.0_ki*abb61(12)
      abb61(7)=-abb61(7)+abb61(12)
      abb61(6)=2.0_ki*abb61(6)
      abb61(8)=-es61+abb61(8)
      abb61(12)=abb61(5)*spbe6k4
      abb61(8)=abb61(12)*abb61(8)
      abb61(21)=-spak1k5*abb61(21)
      abb61(8)=abb61(21)+abb61(8)
      abb61(21)=2.0_ki*abb61(12)
      abb61(22)=-spbk6k4*spbk6e6*spak5k6
      abb61(19)=spak1k5*abb61(19)
      abb61(19)=abb61(22)+abb61(19)
      abb61(22)=abb61(4)*spak2e6
      abb61(19)=2.0_ki*abb61(22)*abb61(19)
      abb61(5)=4.0_ki*abb61(5)
      abb61(23)=spbk6k4*abb61(5)
      abb61(24)=spbk6k1*spak2k6
      abb61(3)=abb61(24)*abb61(4)*abb61(3)
      abb61(25)=abb61(22)*spbe6k1
      abb61(26)=es61*abb61(25)
      abb61(3)=abb61(3)+abb61(26)
      abb61(3)=2.0_ki*abb61(3)
      abb61(26)=4.0_ki*abb61(4)
      abb61(24)=abb61(26)*abb61(24)
      abb61(27)=-4.0_ki*abb61(25)
      abb61(28)=-spbe6k1*abb61(18)
      abb61(17)=-spbk6e6*abb61(17)
      abb61(29)=-spbk6e6*abb61(18)
      abb61(4)=-spbk6k4*abb61(4)*spak2k6
      abb61(10)=-spak1k2*abb61(10)
      abb61(4)=abb61(4)+abb61(10)
      abb61(4)=4.0_ki*abb61(4)
      abb61(10)=spbe6k4*abb61(18)
      abb61(18)=spbe6k4*abb61(22)
      abb61(11)=-4.0_ki*abb61(11)
      abb61(5)=spbk4k2*abb61(5)
      R2d61=0.0_ki
      rat2 = rat2 + R2d61
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='61' value='", &
          & R2d61, "'/>"
      end if
   end subroutine
end module p12_cbbar_hepneg_abbrevd61h0
