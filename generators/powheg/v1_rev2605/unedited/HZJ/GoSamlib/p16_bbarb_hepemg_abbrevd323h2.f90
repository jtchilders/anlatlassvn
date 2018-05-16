module     p16_bbarb_hepemg_abbrevd323h2
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_globalsh2
   implicit none
   private
   complex(ki), dimension(19), public :: abb323
   complex(ki), public :: R2d323
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
      abb323(1)=1.0_ki/(-mZ**2+es345+i_*mZ*wZ)
      abb323(2)=1.0_ki/(-mZ**2+es45+i_*mZ*wZ)
      abb323(3)=NC**(-1)
      abb323(4)=es61**(-1)
      abb323(5)=abb323(3)-NC
      abb323(6)=-spbk6e6*spak4k6*abb323(5)
      abb323(7)=-spbe6k1*abb323(5)
      abb323(8)=-spak1k4*abb323(7)
      abb323(9)=abb323(6)+abb323(8)
      abb323(10)=abb323(4)*spbk5k2
      abb323(11)=gHZZ*gBl*ger*abb323(1)*abb323(2)*spak1e6*c1*TR*i_
      abb323(12)=abb323(10)*abb323(11)
      abb323(13)=-abb323(9)*abb323(12)
      abb323(14)=spak1k4*es12
      abb323(15)=spbk6k2*spak4k6*spak1k2
      abb323(14)=abb323(14)-abb323(15)
      abb323(14)=abb323(14)*abb323(7)
      abb323(15)=es12-es345
      abb323(15)=abb323(15)*spak4k6
      abb323(16)=spbk2k1*spak1k4*spak2k6
      abb323(15)=abb323(15)-abb323(16)
      abb323(5)=abb323(5)*spbk6e6
      abb323(15)=-abb323(15)*abb323(5)
      abb323(14)=abb323(15)+abb323(14)
      abb323(10)=abb323(14)*abb323(10)
      abb323(14)=spbk5k2*abb323(6)
      abb323(10)=abb323(14)+abb323(10)
      abb323(10)=2.0_ki*abb323(10)*abb323(11)
      abb323(12)=4.0_ki*abb323(12)
      abb323(6)=-abb323(6)*abb323(12)
      abb323(8)=abb323(8)*abb323(12)
      abb323(12)=2.0_ki*abb323(4)
      abb323(11)=abb323(12)*abb323(11)
      abb323(12)=abb323(11)*spbk5k2
      abb323(14)=abb323(9)*abb323(12)
      abb323(15)=-abb323(11)*spbk6k5*abb323(9)
      abb323(16)=abb323(12)*spak1k4*abb323(5)
      abb323(17)=-abb323(11)*abb323(5)
      abb323(18)=-spbk6k5*abb323(7)
      abb323(19)=-spbk5k1*abb323(5)
      abb323(18)=abb323(19)+abb323(18)
      abb323(18)=spak1k6*abb323(18)
      abb323(19)=spak1k2*abb323(7)
      abb323(5)=spak2k6*abb323(5)
      abb323(5)=abb323(5)+abb323(19)
      abb323(5)=spbk5k2*abb323(5)
      abb323(5)=abb323(5)+abb323(18)
      abb323(5)=abb323(5)*abb323(11)
      abb323(12)=-abb323(12)*spak4k6*abb323(7)
      abb323(9)=-abb323(11)*spbk5k1*abb323(9)
      abb323(7)=abb323(7)*abb323(11)
      R2d323=abb323(13)
      rat2 = rat2 + R2d323
      if (debug_nlo_diagrams) then
          write (logfile,*) "<result name='r2' index='323' value='", &
          & R2d323, "'/>"
      end if
   end subroutine
end module p16_bbarb_hepemg_abbrevd323h2
