module     p6_ubbar_hepneg_abbrevd37h0
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(14), public :: abb37
   complex(ki), public :: R2d37
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p6_ubbar_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p6_ubbar_hepneg_kinematics
      use p6_ubbar_hepneg_model
      use p6_ubbar_hepneg_color, only: TR
      use p6_ubbar_hepneg_globalsl1, only: epspow
      implicit none
      abb37(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb37(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb37(3)=NC**(-1)
      abb37(4)=es61**(-1)
      abb37(5)=spae6k6*spbk6k4
      abb37(6)=spak2k5*c1*i_*TR*gHWW*CVBU*abb37(2)*abb37(1)
      abb37(7)=abb37(6)*spbe6k1
      abb37(8)=abb37(5)*abb37(7)
      abb37(9)=abb37(7)*spbk4k1
      abb37(10)=-spak1e6*abb37(9)
      abb37(8)=abb37(8)+abb37(10)
      abb37(10)=NC-abb37(3)
      abb37(11)=-abb37(10)*abb37(4)*gW**2
      abb37(8)=2.0_ki*abb37(11)*abb37(8)
      abb37(6)=abb37(6)*spbk4k1
      abb37(11)=-8.0_ki*abb37(11)*abb37(6)
      abb37(12)=gW*abb37(4)
      abb37(10)=abb37(10)*abb37(12)**2
      abb37(7)=abb37(7)*abb37(10)
      abb37(5)=abb37(5)*abb37(7)
      abb37(9)=abb37(9)*abb37(10)
      abb37(12)=2.0_ki*spak1e6
      abb37(13)=abb37(9)*abb37(12)
      abb37(13)=-abb37(5)+abb37(13)
      abb37(13)=4.0_ki*abb37(13)
      abb37(14)=-spak1e6*abb37(9)
      abb37(5)=2.0_ki*abb37(5)+abb37(14)
      abb37(5)=4.0_ki*abb37(5)
      abb37(7)=abb37(7)*spbk6k4
      abb37(6)=spbk6e6*abb37(10)*abb37(6)
      abb37(6)=-abb37(7)+abb37(6)
      abb37(6)=4.0_ki*spak1k6*abb37(6)
      abb37(7)=abb37(7)*abb37(12)
      abb37(9)=2.0_ki*spae6k6*abb37(9)
      R2d37=0.0_ki
      rat2 = rat2 + R2d37
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='37' value='", &
          & R2d37, "'/>"
      end if
   end subroutine
end module p6_ubbar_hepneg_abbrevd37h0
