module     p16_bbarb_hepemg_abbrevd413h3
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_globalsh3
   implicit none
   private
   complex(ki), dimension(16), public :: abb413
   complex(ki), public :: R2d413
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p16_bbarb_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_color, only: TR
      use p16_bbarb_hepemg_globalsl1, only: epspow
      implicit none
      abb413(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb413(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb413(3)=NC**(-1)
      abb413(4)=es61**(-1)
      abb413(5)=gHZZ*spbe6k1*c1*i_*TR*gBr*ger*abb413(2)*abb413(1)
      abb413(6)=abb413(5)*spae6k6
      abb413(7)=abb413(6)*abb413(3)
      abb413(8)=abb413(5)*NC
      abb413(9)=-spae6k6*abb413(8)
      abb413(7)=abb413(7)+abb413(9)
      abb413(7)=spbk6k5*abb413(7)
      abb413(9)=abb413(8)*spak1e6
      abb413(10)=abb413(5)*spak1e6
      abb413(11)=abb413(10)*abb413(3)
      abb413(9)=-abb413(11)+abb413(9)
      abb413(9)=spbk5k1*abb413(9)
      abb413(7)=abb413(7)+abb413(9)
      abb413(7)=2.0_ki*abb413(7)*spak2k4*abb413(4)
      abb413(9)=abb413(3)-NC
      abb413(11)=abb413(4)**2
      abb413(5)=abb413(9)*abb413(5)*abb413(11)*spak2k4
      abb413(9)=-spae6k6*abb413(5)
      abb413(12)=-spbk6k5*abb413(9)
      abb413(5)=-spak1e6*abb413(5)
      abb413(13)=2.0_ki*spbk5k1
      abb413(14)=abb413(5)*abb413(13)
      abb413(12)=abb413(12)+abb413(14)
      abb413(12)=4.0_ki*abb413(12)
      abb413(14)=2.0_ki*spbk6k5
      abb413(15)=abb413(9)*abb413(14)
      abb413(16)=-spbk5k1*abb413(5)
      abb413(15)=abb413(15)+abb413(16)
      abb413(15)=4.0_ki*abb413(15)
      abb413(5)=abb413(5)*abb413(14)
      abb413(8)=abb413(8)*abb413(11)
      abb413(14)=abb413(8)*spak1e6
      abb413(11)=abb413(11)*abb413(3)
      abb413(10)=abb413(11)*abb413(10)
      abb413(10)=abb413(14)-abb413(10)
      abb413(14)=-spak2k6*abb413(10)
      abb413(8)=abb413(8)*spae6k6
      abb413(6)=abb413(11)*abb413(6)
      abb413(6)=abb413(8)-abb413(6)
      abb413(8)=spak1k2*abb413(6)
      abb413(8)=abb413(14)+abb413(8)
      abb413(11)=4.0_ki*spbk6k1
      abb413(8)=abb413(8)*abb413(11)
      abb413(10)=spak4k6*abb413(10)
      abb413(6)=-spak1k4*abb413(6)
      abb413(6)=abb413(10)+abb413(6)
      abb413(6)=abb413(6)*abb413(11)
      abb413(9)=abb413(9)*abb413(13)
      R2d413=0.0_ki
      rat2 = rat2 + R2d413
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='413' value='", &
          & R2d413, "'/>"
      end if
   end subroutine
end module p16_bbarb_hepemg_abbrevd413h3
