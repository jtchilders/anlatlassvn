module     p1_dbarc_hepneg_abbrevd34h0
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(21), public :: abb34
   complex(ki), public :: R2d34
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
      abb34(1)=1.0_ki/(-es61-es12+es345)
      abb34(2)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb34(3)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb34(4)=NC**(-1)
      abb34(5)=spae6k6*spak1k5
      abb34(6)=abb34(5)*spbk6k4
      abb34(7)=NC-abb34(4)
      abb34(7)=abb34(7)*gHWW*c1*gW**2*i_*TR*CVDC*abb34(3)*abb34(2)*abb34(1)
      abb34(8)=spbe6k2*abb34(7)
      abb34(9)=abb34(6)*abb34(8)
      abb34(10)=spak1k2*spbk4k2
      abb34(5)=-spbk6k1*abb34(10)*abb34(5)
      abb34(6)=es61*abb34(6)
      abb34(5)=abb34(6)+abb34(5)
      abb34(5)=2.0_ki*abb34(8)*abb34(5)
      abb34(6)=-es345+es12+es61
      abb34(11)=spbk4k2*spak1k5
      abb34(6)=abb34(11)*abb34(6)
      abb34(12)=spbk4k1*spak1k5
      abb34(13)=spbk6k2*spak1k6
      abb34(14)=-abb34(12)*abb34(13)
      abb34(6)=abb34(14)+abb34(6)
      abb34(6)=-abb34(7)*abb34(6)
      abb34(14)=spak1k6*spbk6k4
      abb34(10)=abb34(14)+abb34(10)
      abb34(14)=-spak5k6*abb34(10)*abb34(7)
      abb34(15)=-spbk6k2*abb34(14)
      abb34(6)=abb34(15)+abb34(6)
      abb34(6)=4.0_ki*abb34(6)
      abb34(10)=abb34(10)*abb34(8)
      abb34(15)=abb34(10)*spak5e6
      abb34(16)=abb34(8)*spak1e6
      abb34(17)=abb34(16)*abb34(12)
      abb34(15)=abb34(15)+abb34(17)
      abb34(15)=4.0_ki*abb34(15)
      abb34(17)=4.0_ki*abb34(9)
      abb34(18)=-2.0_ki*abb34(9)
      abb34(19)=abb34(7)*spbk6e6
      abb34(20)=-spak1k6*abb34(19)
      abb34(21)=abb34(8)*spak1k2
      abb34(20)=abb34(20)-abb34(21)
      abb34(12)=abb34(20)*abb34(12)
      abb34(14)=spbk6e6*abb34(14)
      abb34(10)=spak2k5*abb34(10)
      abb34(11)=abb34(11)*abb34(19)
      abb34(19)=-abb34(8)*spbk6k4*spak1k5
      abb34(11)=abb34(19)+abb34(11)
      abb34(11)=spak2k6*abb34(11)
      abb34(10)=abb34(11)+abb34(10)+abb34(12)+abb34(14)
      abb34(10)=2.0_ki*abb34(10)
      abb34(11)=2.0_ki*spae6k6
      abb34(12)=abb34(11)*spbk6k2
      abb34(14)=abb34(12)*abb34(21)
      abb34(7)=-8.0_ki*abb34(7)*abb34(13)
      abb34(13)=8.0_ki*abb34(16)
      abb34(16)=-4.0_ki*abb34(20)
      abb34(11)=abb34(8)*abb34(11)
      abb34(8)=spak2k5*abb34(12)*abb34(8)
      R2d34=abb34(9)
      rat2 = rat2 + R2d34
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='34' value='", &
          & R2d34, "'/>"
      end if
   end subroutine
end module p1_dbarc_hepneg_abbrevd34h0
