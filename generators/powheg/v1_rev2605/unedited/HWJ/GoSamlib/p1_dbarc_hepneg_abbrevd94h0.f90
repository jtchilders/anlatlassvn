module     p1_dbarc_hepneg_abbrevd94h0
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(22), public :: abb94
   complex(ki), public :: R2d94
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
      abb94(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb94(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb94(3)=NC**(-1)
      abb94(4)=es61**(-1)
      abb94(5)=spak1e6*abb94(3)
      abb94(6)=NC*spak1e6
      abb94(7)=abb94(6)-abb94(5)
      abb94(8)=spbk6e6*spak5k6
      abb94(9)=spak1k5*spbe6k1
      abb94(10)=abb94(8)-abb94(9)
      abb94(7)=abb94(10)*abb94(7)
      abb94(10)=CVDC*gW**2*abb94(1)*abb94(2)*gHWW*c1*TR*i_
      abb94(11)=abb94(10)*abb94(4)
      abb94(12)=spbk4k2*abb94(11)
      abb94(13)=-abb94(7)*abb94(12)
      abb94(14)=-2.0_ki*abb94(8)+abb94(9)
      abb94(14)=abb94(14)*abb94(6)*abb94(10)*spbk4k2
      abb94(15)=NC*spak1k5
      abb94(16)=abb94(3)*spak1k5
      abb94(15)=abb94(15)+abb94(16)
      abb94(10)=4.0_ki*abb94(10)
      abb94(10)=-spbk4k2*abb94(15)*abb94(10)
      abb94(17)=abb94(9)*spak1e6
      abb94(18)=spbk6e6*spak1k6
      abb94(19)=abb94(18)*spak5e6
      abb94(17)=abb94(17)+abb94(19)
      abb94(17)=abb94(17)*NC
      abb94(19)=abb94(19)*abb94(3)
      abb94(17)=abb94(17)-abb94(19)
      abb94(9)=-abb94(9)*abb94(5)
      abb94(9)=abb94(9)-abb94(17)
      abb94(19)=4.0_ki*abb94(12)
      abb94(9)=abb94(9)*abb94(19)
      abb94(15)=abb94(15)*abb94(12)
      abb94(20)=16.0_ki*abb94(15)
      abb94(5)=abb94(8)*abb94(5)
      abb94(5)=abb94(5)+abb94(17)
      abb94(5)=abb94(5)*abb94(19)
      abb94(8)=-8.0_ki*abb94(15)
      abb94(12)=2.0_ki*abb94(12)
      abb94(7)=abb94(7)*abb94(12)
      abb94(12)=-spae6k6*spbk6k1*abb94(16)*abb94(12)
      abb94(6)=abb94(11)*abb94(6)
      abb94(15)=spak1k6*spbk6k2
      abb94(16)=spbe6k1*abb94(15)
      abb94(17)=3.0_ki*abb94(18)
      abb94(18)=spbk2k1*abb94(17)
      abb94(16)=abb94(16)+abb94(18)
      abb94(16)=abb94(16)*abb94(6)
      abb94(11)=8.0_ki*abb94(11)
      abb94(18)=NC+abb94(3)
      abb94(11)=abb94(11)*abb94(18)
      abb94(15)=-abb94(15)*abb94(11)
      abb94(21)=spak1k6*spbk6k4
      abb94(22)=-spbe6k1*abb94(21)
      abb94(17)=-spbk4k1*abb94(17)
      abb94(17)=abb94(22)+abb94(17)
      abb94(6)=abb94(17)*abb94(6)
      abb94(11)=abb94(21)*abb94(11)
      abb94(17)=abb94(19)*spak5k6*abb94(18)
      R2d94=abb94(13)
      rat2 = rat2 + R2d94
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='94' value='", &
          & R2d94, "'/>"
      end if
   end subroutine
end module p1_dbarc_hepneg_abbrevd94h0
