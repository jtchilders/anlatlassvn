module     p4_ubaru_hepemg_abbrevd421h0
   use p4_ubaru_hepemg_config, only: ki
   use p4_ubaru_hepemg_globalsh0
   implicit none
   private
   complex(ki), dimension(11), public :: abb421
   complex(ki), public :: R2d421
   public :: init_abbrev
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
contains
   subroutine     init_abbrev()
      use p4_ubaru_hepemg_config, only: deltaOS, &
     &    logfile, debug_nlo_diagrams
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_color, only: TR
      use p4_ubaru_hepemg_globalsl1, only: epspow
      implicit none
      abb421(1)=1.0_ki/(-es61-es12+es345)
      abb421(2)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb421(3)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb421(4)=NC**(-1)
      abb421(5)=NC-abb421(4)
      abb421(6)=gHZZ*spak1k5*c1*abb421(1)**2*i_*TR*gUl*gel*abb421(3)*abb421(2)
      abb421(7)=abb421(5)*abb421(6)*spbk6k4*spbe6k2
      abb421(8)=-spae6k6*abb421(7)
      abb421(9)=-es61+es345-es12
      abb421(10)=-2.0_ki*abb421(8)*abb421(9)
      abb421(5)=-spbk4k2*abb421(6)*abb421(5)
      abb421(6)=8.0_ki*abb421(9)*abb421(5)
      abb421(9)=-4.0_ki*abb421(8)
      abb421(8)=8.0_ki*abb421(8)
      abb421(11)=spbk6e6*abb421(5)
      abb421(7)=abb421(7)+abb421(11)
      abb421(7)=4.0_ki*spak2k6*abb421(7)
      abb421(5)=2.0_ki*spae6k6*spbe6k2*abb421(5)
      R2d421=0.0_ki
      rat2 = rat2 + R2d421
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='421' value='", &
          & R2d421, "'/>"
      end if
   end subroutine
end module p4_ubaru_hepemg_abbrevd421h0
