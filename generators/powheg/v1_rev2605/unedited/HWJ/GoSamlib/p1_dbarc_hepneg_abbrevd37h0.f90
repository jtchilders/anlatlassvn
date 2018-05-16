module     p1_dbarc_hepneg_abbrevd37h0
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(16), public :: abb37
   complex(ki), public :: R2d37
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
      abb37(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb37(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb37(3)=NC**(-1)
      abb37(4)=es61**(-1)
      abb37(5)=spak1e6*c1*i_*TR*gHWW*CVDC*abb37(2)*abb37(1)
      abb37(6)=abb37(5)*spbk6e6
      abb37(7)=abb37(6)*abb37(3)
      abb37(8)=abb37(5)*NC
      abb37(9)=-spbk6e6*abb37(8)
      abb37(7)=abb37(7)+abb37(9)
      abb37(7)=spak5k6*abb37(7)
      abb37(9)=abb37(5)*spbe6k1
      abb37(10)=-abb37(3)*abb37(9)
      abb37(11)=spbe6k1*abb37(8)
      abb37(10)=abb37(10)+abb37(11)
      abb37(10)=spak1k5*abb37(10)
      abb37(7)=abb37(7)+abb37(10)
      abb37(7)=2.0_ki*abb37(7)*abb37(4)*spbk4k2*gW**2
      abb37(10)=abb37(3)-NC
      abb37(11)=gW*abb37(4)
      abb37(11)=abb37(11)**2
      abb37(5)=abb37(10)*abb37(5)*abb37(11)*spbk4k2
      abb37(10)=-spbk6e6*abb37(5)
      abb37(12)=-spak5k6*abb37(10)
      abb37(5)=-spbe6k1*abb37(5)
      abb37(13)=2.0_ki*spak1k5
      abb37(14)=abb37(5)*abb37(13)
      abb37(12)=abb37(12)+abb37(14)
      abb37(12)=4.0_ki*abb37(12)
      abb37(14)=2.0_ki*spak5k6
      abb37(15)=abb37(10)*abb37(14)
      abb37(16)=-spak1k5*abb37(5)
      abb37(15)=abb37(15)+abb37(16)
      abb37(15)=4.0_ki*abb37(15)
      abb37(10)=abb37(10)*abb37(13)
      abb37(8)=abb37(8)*abb37(11)
      abb37(13)=abb37(8)*spbe6k1
      abb37(11)=abb37(11)*abb37(3)
      abb37(9)=abb37(11)*abb37(9)
      abb37(9)=abb37(13)-abb37(9)
      abb37(13)=-spbk6k2*abb37(9)
      abb37(8)=abb37(8)*spbk6e6
      abb37(6)=abb37(11)*abb37(6)
      abb37(6)=abb37(8)-abb37(6)
      abb37(8)=spbk2k1*abb37(6)
      abb37(8)=abb37(13)+abb37(8)
      abb37(11)=4.0_ki*spak1k6
      abb37(8)=abb37(8)*abb37(11)
      abb37(9)=spbk6k4*abb37(9)
      abb37(6)=-spbk4k1*abb37(6)
      abb37(6)=abb37(9)+abb37(6)
      abb37(6)=abb37(6)*abb37(11)
      abb37(5)=abb37(5)*abb37(14)
      R2d37=0.0_ki
      rat2 = rat2 + R2d37
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='37' value='", &
          & R2d37, "'/>"
      end if
   end subroutine
end module p1_dbarc_hepneg_abbrevd37h0
