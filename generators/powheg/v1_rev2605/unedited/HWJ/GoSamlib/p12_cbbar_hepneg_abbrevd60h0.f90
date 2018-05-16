module     p12_cbbar_hepneg_abbrevd60h0
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(30), public :: abb60
   complex(ki), public :: R2d60
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
      abb60(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb60(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb60(3)=NC**(-1)
      abb60(4)=spbk4k1*spak5k6
      abb60(5)=CVBC*abb60(3)*gW**2*abb60(1)*abb60(2)*gHWW*c1*TR*i_
      abb60(6)=2.0_ki*abb60(5)
      abb60(7)=abb60(6)*spak2e6
      abb60(8)=abb60(7)*spbk6e6
      abb60(9)=-es12*abb60(8)*abb60(4)
      abb60(10)=-spak1k2*spak5k6*spbk6k1
      abb60(11)=-spak2k5*es12
      abb60(10)=abb60(10)+abb60(11)
      abb60(11)=4.0_ki*abb60(5)
      abb60(12)=abb60(11)*spbk4k1
      abb60(10)=abb60(10)*abb60(12)
      abb60(13)=spbe6k1*spak1k2
      abb60(14)=spbk4k1*spak5e6
      abb60(15)=abb60(13)*abb60(14)
      abb60(16)=abb60(11)*abb60(15)
      abb60(4)=-abb60(4)*spak2e6*spbk6e6
      abb60(4)=-abb60(15)+abb60(4)
      abb60(4)=abb60(4)*abb60(11)
      abb60(5)=-8.0_ki*spak2k5*abb60(5)*spbk4k1
      abb60(17)=-spbe6k4*spbk2k1
      abb60(18)=-spbe6k1*spbk4k2
      abb60(17)=abb60(17)+abb60(18)
      abb60(17)=spak2e6*abb60(17)
      abb60(18)=spbe6k1*spae6k6*spbk6k4
      abb60(17)=abb60(17)+abb60(18)
      abb60(17)=spak2k5*abb60(17)
      abb60(15)=abb60(15)+abb60(17)
      abb60(15)=abb60(15)*abb60(6)
      abb60(17)=abb60(11)*spak2k5
      abb60(18)=spbk4k1*abb60(17)
      abb60(19)=abb60(14)*abb60(11)
      abb60(14)=abb60(6)*abb60(14)
      abb60(20)=abb60(17)*spbe6k4
      abb60(21)=-spak2k5*spbe6k4*abb60(6)
      abb60(22)=-spak2k5*spbk4k2*abb60(8)
      abb60(23)=spbk6k4*abb60(17)
      abb60(24)=abb60(13)*spae6k6
      abb60(25)=spbk6k1*abb60(6)*abb60(24)
      abb60(26)=abb60(11)*spak2e6
      abb60(27)=-spbe6k1*abb60(26)
      abb60(28)=spae6k6*spbk6k1
      abb60(29)=spak2e6*spbk2k1
      abb60(28)=abb60(28)+abb60(29)
      abb60(28)=abb60(28)*abb60(6)
      abb60(29)=spbk6e6*spak2k6
      abb60(30)=-abb60(6)*abb60(29)
      abb60(13)=abb60(13)+abb60(29)
      abb60(13)=spbk4k2*abb60(13)
      abb60(29)=spbe6k4*es12
      abb60(13)=abb60(29)+abb60(13)
      abb60(13)=spak2e6*abb60(13)
      abb60(24)=-spbk6k4*abb60(24)
      abb60(13)=abb60(24)+abb60(13)
      abb60(6)=abb60(13)*abb60(6)
      abb60(13)=-spbk6k4*spak2k6
      abb60(24)=-spbk4k1*spak1k2
      abb60(13)=abb60(13)+abb60(24)
      abb60(13)=abb60(13)*abb60(11)
      abb60(24)=spbe6k4*abb60(26)
      abb60(7)=spbe6k4*abb60(7)
      abb60(12)=-spak5k6*abb60(12)
      abb60(17)=spbk4k2*abb60(17)
      R2d60=0.0_ki
      rat2 = rat2 + R2d60
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='60' value='", &
          & R2d60, "'/>"
      end if
   end subroutine
end module p12_cbbar_hepneg_abbrevd60h0
