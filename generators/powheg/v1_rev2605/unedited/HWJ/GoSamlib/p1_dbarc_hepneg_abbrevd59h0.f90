module     p1_dbarc_hepneg_abbrevd59h0
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(32), public :: abb59
   complex(ki), public :: R2d59
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p1_dbarc_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_model
      use p1_dbarc_hepneg_color, only: TR
      use p1_dbarc_hepneg_globalsl1, only: epspow
      implicit none
      abb59(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb59(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb59(3)=NC**(-1)
      abb59(4)=spbk6e6*spak5k6
      abb59(5)=spak1k5*spbe6k1
      abb59(4)=abb59(4)-abb59(5)
      abb59(4)=abb59(4)*spbk4k2
      abb59(6)=CVDC*abb59(3)*gW**2*abb59(1)*abb59(2)*gHWW*c1*TR*i_
      abb59(7)=2.0_ki*abb59(6)
      abb59(8)=abb59(7)*spak1e6
      abb59(9)=abb59(8)*es12*abb59(4)
      abb59(10)=spak1k5*es12
      abb59(11)=-spak1k2*spak5k6*spbk6k2
      abb59(10)=abb59(11)+abb59(10)
      abb59(11)=4.0_ki*abb59(6)
      abb59(12)=abb59(11)*spbk4k2
      abb59(10)=abb59(10)*abb59(12)
      abb59(13)=spbe6k2*spak1k2
      abb59(14)=spbk4k2*spak5e6
      abb59(15)=abb59(13)*abb59(14)
      abb59(16)=abb59(11)*abb59(15)
      abb59(4)=spak1e6*abb59(4)
      abb59(4)=-abb59(15)+abb59(4)
      abb59(4)=abb59(4)*abb59(11)
      abb59(6)=8.0_ki*spak1k5*abb59(6)*spbk4k2
      abb59(17)=-spbe6k4*spbk2k1
      abb59(18)=spbe6k2*spbk4k1
      abb59(17)=abb59(17)+abb59(18)
      abb59(17)=abb59(17)*spak1e6
      abb59(18)=-spbe6k2*spae6k6*spbk6k4
      abb59(17)=abb59(17)+abb59(18)
      abb59(17)=spak1k5*abb59(17)
      abb59(15)=abb59(15)+abb59(17)
      abb59(15)=abb59(15)*abb59(7)
      abb59(17)=abb59(11)*spak1k5
      abb59(18)=-spbk4k2*abb59(17)
      abb59(19)=abb59(17)*spbe6k4
      abb59(20)=spbe6k4*abb59(7)*spak1k5
      abb59(21)=abb59(14)*abb59(11)
      abb59(14)=-abb59(7)*abb59(14)
      abb59(22)=spbk6e6*spbk4k1
      abb59(23)=abb59(22)*abb59(8)*spak1k5
      abb59(24)=-spbk6k4*abb59(17)
      abb59(25)=abb59(13)*spae6k6
      abb59(26)=spbk6k2*abb59(7)*abb59(25)
      abb59(27)=abb59(11)*spak1e6
      abb59(28)=spbe6k2*abb59(27)
      abb59(29)=spak1k6*abb59(7)*spbk6e6
      abb59(30)=-spae6k6*spbk6k2
      abb59(31)=spak1e6*spbk2k1
      abb59(30)=abb59(30)+abb59(31)
      abb59(30)=abb59(30)*abb59(7)
      abb59(31)=abb59(8)*spbk6e6
      abb59(22)=-spak1k6*abb59(22)
      abb59(32)=-spbe6k4*es12
      abb59(13)=spbk4k1*abb59(13)
      abb59(13)=abb59(13)+abb59(22)+abb59(32)
      abb59(13)=spak1e6*abb59(13)
      abb59(22)=-spbk6k4*abb59(25)
      abb59(13)=abb59(22)+abb59(13)
      abb59(7)=abb59(13)*abb59(7)
      abb59(13)=spbk6k4*spak1k6
      abb59(22)=-spbk4k2*spak1k2
      abb59(13)=abb59(13)+abb59(22)
      abb59(13)=abb59(13)*abb59(11)
      abb59(22)=abb59(27)*spbe6k4
      abb59(25)=-spbe6k4*abb59(8)
      abb59(12)=spak5k6*abb59(12)
      abb59(5)=abb59(5)*abb59(8)*spbk4k1
      abb59(8)=-spbk4k1*abb59(17)
      abb59(17)=-spbe6k1*abb59(27)
      R2d59=0.0_ki
      rat2 = rat2 + R2d59
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='59' value='", &
          & R2d59, "'/>"
      end if
   end subroutine
end module p1_dbarc_hepneg_abbrevd59h0
