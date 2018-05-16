module     p1_dbarc_hepneg_abbrevd60h0
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(30), public :: abb60
   complex(ki), public :: R2d60
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
      abb60(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb60(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb60(3)=NC**(-1)
      abb60(4)=spbe6k2*spae6k6
      abb60(5)=CVDC*abb60(3)*gW**2*abb60(1)*abb60(2)*gHWW*c1*TR*i_
      abb60(6)=2.0_ki*abb60(5)
      abb60(7)=abb60(4)*abb60(6)
      abb60(8)=spak1k5*spbk6k4
      abb60(9)=-es12*abb60(8)*abb60(7)
      abb60(10)=spbk4k2*spak1k2
      abb60(11)=abb60(10)*spak5k6
      abb60(12)=spbk6k2*abb60(11)
      abb60(13)=spak1k5*spbk4k2
      abb60(14)=-es12*abb60(13)
      abb60(12)=abb60(12)+abb60(14)
      abb60(14)=4.0_ki*abb60(5)
      abb60(12)=abb60(12)*abb60(14)
      abb60(15)=spbe6k2*spak5e6
      abb60(16)=abb60(15)*abb60(14)
      abb60(17)=-abb60(10)*abb60(16)
      abb60(18)=abb60(10)*abb60(15)
      abb60(19)=spak1e6*spbe6k1*abb60(13)
      abb60(18)=abb60(18)+abb60(19)
      abb60(18)=abb60(18)*abb60(14)
      abb60(5)=abb60(5)*abb60(13)
      abb60(5)=8.0_ki*abb60(5)
      abb60(13)=abb60(14)*spak1k5
      abb60(19)=abb60(13)*spbk6k4
      abb60(4)=abb60(19)*abb60(4)
      abb60(20)=-spak1k2*spak5e6
      abb60(21)=-spak1e6*spak2k5
      abb60(20)=abb60(20)+abb60(21)
      abb60(20)=spbe6k2*abb60(20)
      abb60(21)=spbk6e6*spak5k6*spak1e6
      abb60(20)=abb60(21)+abb60(20)
      abb60(20)=spbk4k2*abb60(20)
      abb60(21)=spak1k5*spbe6k4
      abb60(22)=spak1e6*spbk2k1*abb60(21)
      abb60(20)=abb60(22)+abb60(20)
      abb60(20)=abb60(20)*abb60(6)
      abb60(22)=spbk6k4*spak1k6
      abb60(22)=abb60(22)+abb60(10)
      abb60(23)=-spbe6k1*abb60(22)
      abb60(24)=-spbe6k4*es12
      abb60(23)=abb60(24)+abb60(23)
      abb60(23)=spak1k5*abb60(23)
      abb60(11)=-spbk6e6*abb60(11)
      abb60(10)=spbe6k2*spak2k5*abb60(10)
      abb60(10)=abb60(23)+abb60(11)+abb60(10)
      abb60(10)=abb60(10)*abb60(6)
      abb60(11)=abb60(13)*spbe6k4
      abb60(21)=-abb60(6)*abb60(21)
      abb60(23)=abb60(6)*spak1e6
      abb60(8)=spbe6k1*abb60(23)*abb60(8)
      abb60(24)=-spbk6k2*spak1k2*abb60(7)
      abb60(25)=abb60(14)*spak1e6
      abb60(26)=spbe6k2*abb60(25)
      abb60(27)=spbe6k2*abb60(23)
      abb60(28)=-spbk6e6*spak1k6
      abb60(29)=2.0_ki*spbe6k2
      abb60(29)=-spak1k2*abb60(29)
      abb60(28)=abb60(28)+abb60(29)
      abb60(28)=abb60(28)*abb60(6)
      abb60(29)=spbk6e6*abb60(23)
      abb60(22)=-abb60(22)*abb60(14)
      abb60(25)=spbe6k4*abb60(25)
      abb60(23)=spbe6k4*abb60(23)
      abb60(30)=abb60(14)*spak5k6*spbk6k2
      abb60(6)=-abb60(6)*abb60(15)
      abb60(13)=-spbk4k1*abb60(13)
      R2d60=0.0_ki
      rat2 = rat2 + R2d60
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='60' value='", &
          & R2d60, "'/>"
      end if
   end subroutine
end module p1_dbarc_hepneg_abbrevd60h0
