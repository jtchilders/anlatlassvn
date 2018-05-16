module     p5_usbar_hepneg_abbrevd30h0
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_globalsh0
   implicit none
   private
   complex(ki), dimension(21), public :: abb30
   complex(ki), public :: R2d30
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p5_usbar_hepneg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_color, only: TR
      use p5_usbar_hepneg_globalsl1, only: epspow
      implicit none
      abb30(1)=1.0_ki/(-mW**2+es345+i_*mW*wW)
      abb30(2)=1.0_ki/(-mW**2+es45+i_*mW*wW)
      abb30(3)=NC**(-1)
      abb30(4)=es61**(-1)
      abb30(5)=spbk4k1*spak2k5
      abb30(6)=abb30(5)*spak1e6
      abb30(7)=spbk6k4*spak2k5
      abb30(8)=abb30(7)*spae6k6
      abb30(9)=abb30(6)-abb30(8)
      abb30(10)=NC-abb30(3)
      abb30(10)=abb30(10)*c1*gW**2*i_*TR*gHWW*CVSU*abb30(2)*abb30(1)
      abb30(11)=abb30(10)*spbe6k1
      abb30(12)=-abb30(4)*abb30(11)
      abb30(9)=abb30(9)*abb30(12)
      abb30(6)=-abb30(8)-abb30(6)
      abb30(6)=es12*abb30(6)
      abb30(13)=spak1e6*spak2k6
      abb30(14)=spbk2k1*abb30(7)*abb30(13)
      abb30(15)=spak1k2*spbk4k1
      abb30(16)=spbk6k2*abb30(15)*spae6k6*spak2k5
      abb30(17)=es345*abb30(8)
      abb30(6)=abb30(6)+abb30(17)+abb30(14)+abb30(16)
      abb30(6)=abb30(12)*abb30(6)
      abb30(8)=abb30(8)*abb30(11)
      abb30(6)=abb30(8)+abb30(6)
      abb30(6)=2.0_ki*abb30(6)
      abb30(8)=-abb30(10)*abb30(5)
      abb30(11)=spak2k6*spbk6k4
      abb30(11)=abb30(11)-abb30(15)
      abb30(10)=abb30(10)*abb30(4)
      abb30(14)=abb30(11)*abb30(10)
      abb30(15)=spak5k6*spbk6k1
      abb30(16)=-abb30(14)*abb30(15)
      abb30(17)=spbk4k2*spak2k5
      abb30(18)=spbk6k1*abb30(10)*spak2k6
      abb30(19)=-abb30(17)*abb30(18)
      abb30(8)=abb30(19)+abb30(8)+abb30(16)
      abb30(8)=4.0_ki*abb30(8)
      abb30(11)=abb30(11)*abb30(12)
      abb30(16)=abb30(11)*spak5e6
      abb30(19)=abb30(12)*spak2e6
      abb30(20)=abb30(19)*abb30(17)
      abb30(16)=abb30(16)+abb30(20)
      abb30(16)=4.0_ki*abb30(16)
      abb30(20)=-4.0_ki*abb30(9)
      abb30(21)=2.0_ki*abb30(9)
      abb30(7)=-abb30(12)*abb30(7)
      abb30(10)=abb30(10)*spbk6e6
      abb30(5)=-abb30(5)*abb30(10)
      abb30(5)=abb30(7)+abb30(5)
      abb30(5)=spak1k6*abb30(5)
      abb30(7)=spak5k6*spbk6e6*abb30(14)
      abb30(14)=abb30(12)*spak1k2
      abb30(10)=spak2k6*abb30(10)
      abb30(10)=abb30(14)+abb30(10)
      abb30(14)=abb30(10)*abb30(17)
      abb30(11)=spak1k5*abb30(11)
      abb30(5)=abb30(5)+abb30(11)+abb30(7)+abb30(14)
      abb30(5)=2.0_ki*abb30(5)
      abb30(7)=-spak1k2*spae6k6
      abb30(7)=abb30(7)+abb30(13)
      abb30(11)=2.0_ki*abb30(12)
      abb30(7)=spbk6k1*abb30(7)*abb30(11)
      abb30(12)=8.0_ki*abb30(18)
      abb30(13)=8.0_ki*abb30(19)
      abb30(10)=-4.0_ki*abb30(10)
      abb30(14)=spae6k6*abb30(11)
      abb30(15)=-spak1e6*abb30(15)
      abb30(17)=spak1k5*spbk6k1*spae6k6
      abb30(15)=abb30(15)+abb30(17)
      abb30(15)=abb30(15)*abb30(11)
      abb30(11)=spak1e6*abb30(11)
      R2d30=-abb30(9)
      rat2 = rat2 + R2d30
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='30' value='", &
          & R2d30, "'/>"
      end if
   end subroutine
end module p5_usbar_hepneg_abbrevd30h0
