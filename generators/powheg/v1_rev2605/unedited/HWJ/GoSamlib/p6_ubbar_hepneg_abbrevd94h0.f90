module     p6_ubbar_hepneg_abbrevd94h0
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(18), public :: abb94
   complex(ki), public :: R2d94
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
      abb94(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb94(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb94(3)=NC**(-1)
      abb94(4)=es61**(-1)
      abb94(5)=spbk6k4*spae6k6
      abb94(6)=spbk4k1*spak1e6
      abb94(7)=abb94(5)-abb94(6)
      abb94(7)=abb94(7)*NC
      abb94(8)=abb94(5)*abb94(3)
      abb94(8)=abb94(7)-abb94(8)
      abb94(9)=abb94(6)*abb94(3)
      abb94(10)=abb94(9)+abb94(8)
      abb94(11)=gHWW*CVBU*gW**2*abb94(1)*abb94(2)*spak2k5*c1*TR*i_
      abb94(12)=abb94(11)*abb94(4)
      abb94(13)=spbe6k1*abb94(12)
      abb94(14)=-abb94(10)*abb94(13)
      abb94(5)=-2.0_ki*abb94(5)+abb94(6)
      abb94(5)=abb94(5)*NC*abb94(11)*spbe6k1
      abb94(11)=2.0_ki*abb94(11)
      abb94(15)=NC*spbk4k1
      abb94(11)=abb94(15)*abb94(11)
      abb94(16)=4.0_ki*abb94(13)
      abb94(8)=-abb94(8)*abb94(16)
      abb94(7)=-abb94(7)*abb94(16)
      abb94(16)=spbk4k1*abb94(3)
      abb94(17)=-abb94(16)-abb94(15)
      abb94(17)=8.0_ki*abb94(17)*abb94(12)
      abb94(13)=2.0_ki*abb94(13)
      abb94(10)=abb94(10)*abb94(13)
      abb94(15)=abb94(15)+2.0_ki*abb94(16)
      abb94(15)=spbk6e6*abb94(15)
      abb94(18)=spbk6k4*spbe6k1*NC
      abb94(15)=abb94(18)+abb94(15)
      abb94(15)=abb94(12)*spak1k6*abb94(15)
      abb94(6)=NC*abb94(6)
      abb94(6)=-abb94(9)+abb94(6)
      abb94(6)=2.0_ki*abb94(12)*spbk6e6*abb94(6)
      abb94(9)=-abb94(3)-NC
      abb94(9)=4.0_ki*abb94(12)*spbk6k4*abb94(9)
      abb94(12)=spae6k6*abb94(16)*abb94(13)
      R2d94=abb94(14)
      rat2 = rat2 + R2d94
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='94' value='", &
          & R2d94, "'/>"
      end if
   end subroutine
end module p6_ubbar_hepneg_abbrevd94h0
