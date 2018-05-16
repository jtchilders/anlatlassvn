module     p0_dbaru_hepneg_abbrevd40h0
   use p0_dbaru_hepneg_config, only: ki
   use p0_dbaru_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(11), public :: abb40
   complex(ki), public :: R2d40
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p0_dbaru_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_model
      use p0_dbaru_hepneg_color, only: TR
      use p0_dbaru_hepneg_globalsl1, only: epspow
      implicit none
      abb40(1)=1.0_ki/(-es61-es12+es345)
      abb40(2)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb40(3)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb40(4)=NC**(-1)
      abb40(5)=NC-abb40(4)
      abb40(6)=gW*abb40(1)
      abb40(6)=spak1k5*abb40(6)**2*c1*i_*TR*gHWW*CVDU*abb40(3)*abb40(2)
      abb40(7)=abb40(5)*abb40(6)*spbk6k4*spbe6k2
      abb40(8)=-spae6k6*abb40(7)
      abb40(9)=-es61+es345-es12
      abb40(10)=-2.0_ki*abb40(8)*abb40(9)
      abb40(5)=-spbk4k2*abb40(6)*abb40(5)
      abb40(6)=8.0_ki*abb40(9)*abb40(5)
      abb40(9)=-4.0_ki*abb40(8)
      abb40(8)=8.0_ki*abb40(8)
      abb40(11)=spbk6e6*abb40(5)
      abb40(7)=abb40(7)+abb40(11)
      abb40(7)=4.0_ki*spak2k6*abb40(7)
      abb40(5)=2.0_ki*spae6k6*spbe6k2*abb40(5)
      R2d40=0.0_ki
      rat2 = rat2 + R2d40
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='40' value='", &
          & R2d40, "'/>"
      end if
   end subroutine
end module p0_dbaru_hepneg_abbrevd40h0
