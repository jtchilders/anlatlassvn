module     p16_bbarb_hepemg_d553h2l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity2d553h2l1d.f90
   ! generator: buildfortran_d.py
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d553
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd553h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd553
      complex(ki) :: brack
      acd553(1)=abb553(8)
      brack=acd553(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd553h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd553
      complex(ki) :: brack
      acd553(1)=k6(iv1)
      acd553(2)=abb553(20)
      acd553(3)=e6(iv1)
      acd553(4)=abb553(14)
      acd553(5)=spvak1k6(iv1)
      acd553(6)=abb553(16)
      acd553(7)=spvak4k6(iv1)
      acd553(8)=abb553(22)
      acd553(9)=spvak6k2(iv1)
      acd553(10)=abb553(18)
      acd553(11)=spvak6k5(iv1)
      acd553(12)=abb553(21)
      acd553(13)=spvak1e6(iv1)
      acd553(14)=abb553(15)
      acd553(15)=spvak2e6(iv1)
      acd553(16)=abb553(10)
      acd553(17)=spvae6k2(iv1)
      acd553(18)=abb553(9)
      acd553(19)=spval3e6(iv1)
      acd553(20)=abb553(27)
      acd553(21)=spvae6l3(iv1)
      acd553(22)=abb553(25)
      acd553(23)=spvak4e6(iv1)
      acd553(24)=abb553(23)
      acd553(25)=spvae6k5(iv1)
      acd553(26)=abb553(29)
      acd553(27)=spvak6e6(iv1)
      acd553(28)=abb553(31)
      acd553(29)=spvae6k6(iv1)
      acd553(30)=abb553(28)
      acd553(31)=-acd553(2)*acd553(1)
      acd553(32)=-acd553(4)*acd553(3)
      acd553(33)=-acd553(6)*acd553(5)
      acd553(34)=-acd553(8)*acd553(7)
      acd553(35)=-acd553(10)*acd553(9)
      acd553(36)=-acd553(12)*acd553(11)
      acd553(37)=-acd553(14)*acd553(13)
      acd553(38)=-acd553(16)*acd553(15)
      acd553(39)=-acd553(18)*acd553(17)
      acd553(40)=-acd553(20)*acd553(19)
      acd553(41)=-acd553(22)*acd553(21)
      acd553(42)=-acd553(24)*acd553(23)
      acd553(43)=-acd553(26)*acd553(25)
      acd553(44)=-acd553(28)*acd553(27)
      acd553(45)=-acd553(30)*acd553(29)
      brack=acd553(31)+acd553(32)+acd553(33)+acd553(34)+acd553(35)+acd553(36)+a&
      &cd553(37)+acd553(38)+acd553(39)+acd553(40)+acd553(41)+acd553(42)+acd553(&
      &43)+acd553(44)+acd553(45)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd553h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(17) :: acd553
      complex(ki) :: brack
      acd553(1)=d(iv1,iv2)
      acd553(2)=abb553(17)
      acd553(3)=e6(iv1)
      acd553(4)=spvak1k2(iv2)
      acd553(5)=abb553(12)
      acd553(6)=spvak1k5(iv2)
      acd553(7)=abb553(13)
      acd553(8)=spvak4k2(iv2)
      acd553(9)=abb553(11)
      acd553(10)=e6(iv2)
      acd553(11)=spvak1k2(iv1)
      acd553(12)=spvak1k5(iv1)
      acd553(13)=spvak4k2(iv1)
      acd553(14)=acd553(9)*acd553(13)
      acd553(15)=acd553(7)*acd553(12)
      acd553(16)=acd553(5)*acd553(11)
      acd553(14)=acd553(16)+acd553(14)+acd553(15)
      acd553(14)=acd553(10)*acd553(14)
      acd553(15)=acd553(9)*acd553(8)
      acd553(16)=acd553(7)*acd553(6)
      acd553(17)=acd553(5)*acd553(4)
      acd553(15)=acd553(17)+acd553(15)+acd553(16)
      acd553(15)=acd553(3)*acd553(15)
      acd553(16)=acd553(1)*acd553(2)
      brack=acd553(14)+acd553(15)+2.0_ki*acd553(16)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd553h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd553
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd553h2
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = 0
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
      end if
      if(present(i2)) then
          iv2=i2
          deg=2
      else
          iv2=1
      end if
      if(present(i3)) then
          iv3=i3
          deg=3
      else
          iv3=1
      end if
      t1 = 0
      if(deg.eq.0) then
         numerator = cond(epspow.eq.t1,brack_1,Q,mu2)
         return
      end if
      if(deg.eq.1) then
         numerator = cond(epspow.eq.t1,brack_2,Q,mu2)
         return
      end if
      if(deg.eq.2) then
         numerator = cond(epspow.eq.t1,brack_3,Q,mu2)
         return
      end if
      if(deg.eq.3) then
         numerator = cond(epspow.eq.t1,brack_4,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
!---#[ subroutine reconstruct_d553:
   subroutine     reconstruct_d553(coeffs)
      use p16_bbarb_hepemg_groups, only: tensrec_info_group3
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group3), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_553:
      coeffs%coeffs_553%c0 = derivative(czip)
      coeffs%coeffs_553%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_553%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_553%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_553%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_553%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_553%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_553%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_553%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_553%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_553%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_553%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_553%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_553%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_553%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_553%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_553%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_553%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_553%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_553%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_553%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_553%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_553%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_553%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_553%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_553%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_553%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_553%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_553%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_553%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_553%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_553%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_553%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_553%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_553%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_553:
   end subroutine reconstruct_d553
!---#] subroutine reconstruct_d553:
end module     p16_bbarb_hepemg_d553h2l1d
