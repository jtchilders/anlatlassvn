module     p12_cbbar_hepneg_abbrevd40h0
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(10), public :: abb40
   complex(ki), public :: R2d40
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
      abb40(1)=1.0_ki/(-es61-es12+es345)
      abb40(2)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb40(3)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb40(4)=NC**(-1)
      abb40(5)=NC-abb40(4)
      abb40(6)=gW*abb40(1)
      abb40(5)=abb40(5)*spbk6e6*spak2e6*abb40(6)**2*c1*i_*TR*gHWW*CVBC*abb40(3)&
      &*abb40(2)
      abb40(6)=abb40(5)*spbk4k1
      abb40(7)=spak5k6*abb40(6)
      abb40(8)=es345-es61-es12
      abb40(8)=2.0_ki*abb40(7)*abb40(8)
      abb40(9)=4.0_ki*abb40(7)
      abb40(7)=-8.0_ki*abb40(7)
      abb40(6)=-2.0_ki*spak2k5*abb40(6)
      abb40(10)=4.0_ki*spak2k6
      abb40(5)=-abb40(5)*abb40(10)
      abb40(10)=-spbk2k1*abb40(5)
      abb40(5)=-spbk4k2*abb40(5)
      R2d40=0.0_ki
      rat2 = rat2 + R2d40
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='40' value='", &
          & R2d40, "'/>"
      end if
   end subroutine
end module p12_cbbar_hepneg_abbrevd40h0
