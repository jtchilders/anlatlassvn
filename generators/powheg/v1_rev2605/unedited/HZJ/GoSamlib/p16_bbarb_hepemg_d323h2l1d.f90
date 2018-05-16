module     p16_bbarb_hepemg_d323h2l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity2d323h2l1d.f90
   ! generator: buildfortran_d.py
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d323
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd323h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(28) :: acd323
      complex(ki) :: brack
      acd323(1)=dotproduct(k1,qshift)
      acd323(2)=abb323(8)
      acd323(3)=dotproduct(k6,qshift)
      acd323(4)=abb323(6)
      acd323(5)=dotproduct(qshift,qshift)
      acd323(6)=abb323(14)
      acd323(7)=dotproduct(qshift,spvak1k2)
      acd323(8)=dotproduct(qshift,spvak4k5)
      acd323(9)=abb323(7)
      acd323(10)=abb323(9)
      acd323(11)=dotproduct(qshift,spvak6k2)
      acd323(12)=abb323(17)
      acd323(13)=dotproduct(qshift,spvak1k6)
      acd323(14)=abb323(12)
      acd323(15)=dotproduct(qshift,spvak4k2)
      acd323(16)=abb323(5)
      acd323(17)=abb323(15)
      acd323(18)=dotproduct(qshift,spvak6k1)
      acd323(19)=abb323(16)
      acd323(20)=abb323(10)
      acd323(21)=acd323(8)*acd323(12)
      acd323(21)=acd323(21)-acd323(17)
      acd323(21)=acd323(11)*acd323(21)
      acd323(22)=-acd323(18)*acd323(19)
      acd323(23)=-acd323(15)*acd323(16)
      acd323(24)=-acd323(13)*acd323(14)
      acd323(25)=acd323(5)*acd323(6)
      acd323(26)=-acd323(3)*acd323(4)
      acd323(27)=-acd323(1)*acd323(2)
      acd323(28)=acd323(8)*acd323(9)
      acd323(28)=-acd323(10)+acd323(28)
      acd323(28)=acd323(7)*acd323(28)
      brack=acd323(20)+acd323(21)+acd323(22)+acd323(23)+acd323(24)+acd323(25)+a&
      &cd323(26)+acd323(27)+acd323(28)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd323h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(31) :: acd323
      complex(ki) :: brack
      acd323(1)=k1(iv1)
      acd323(2)=abb323(8)
      acd323(3)=k6(iv1)
      acd323(4)=abb323(6)
      acd323(5)=qshift(iv1)
      acd323(6)=abb323(14)
      acd323(7)=spvak1k2(iv1)
      acd323(8)=dotproduct(qshift,spvak4k5)
      acd323(9)=abb323(7)
      acd323(10)=abb323(9)
      acd323(11)=spvak4k5(iv1)
      acd323(12)=dotproduct(qshift,spvak1k2)
      acd323(13)=dotproduct(qshift,spvak6k2)
      acd323(14)=abb323(17)
      acd323(15)=spvak1k6(iv1)
      acd323(16)=abb323(12)
      acd323(17)=spvak4k2(iv1)
      acd323(18)=abb323(5)
      acd323(19)=spvak6k2(iv1)
      acd323(20)=abb323(15)
      acd323(21)=spvak6k1(iv1)
      acd323(22)=abb323(16)
      acd323(23)=acd323(12)*acd323(9)
      acd323(24)=acd323(13)*acd323(14)
      acd323(23)=acd323(24)+acd323(23)
      acd323(23)=acd323(11)*acd323(23)
      acd323(24)=acd323(9)*acd323(8)
      acd323(24)=-acd323(10)+acd323(24)
      acd323(24)=acd323(7)*acd323(24)
      acd323(25)=acd323(14)*acd323(8)
      acd323(25)=-acd323(20)+acd323(25)
      acd323(25)=acd323(19)*acd323(25)
      acd323(26)=-acd323(2)*acd323(1)
      acd323(27)=-acd323(4)*acd323(3)
      acd323(28)=acd323(6)*acd323(5)
      acd323(29)=-acd323(16)*acd323(15)
      acd323(30)=-acd323(18)*acd323(17)
      acd323(31)=-acd323(22)*acd323(21)
      brack=acd323(23)+acd323(24)+acd323(25)+acd323(26)+acd323(27)+2.0_ki*acd32&
      &3(28)+acd323(29)+acd323(30)+acd323(31)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd323h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(13) :: acd323
      complex(ki) :: brack
      acd323(1)=d(iv1,iv2)
      acd323(2)=abb323(14)
      acd323(3)=spvak1k2(iv1)
      acd323(4)=spvak4k5(iv2)
      acd323(5)=abb323(7)
      acd323(6)=spvak1k2(iv2)
      acd323(7)=spvak4k5(iv1)
      acd323(8)=spvak6k2(iv2)
      acd323(9)=abb323(17)
      acd323(10)=spvak6k2(iv1)
      acd323(11)=acd323(9)*acd323(8)
      acd323(12)=acd323(5)*acd323(6)
      acd323(11)=acd323(12)+acd323(11)
      acd323(11)=acd323(7)*acd323(11)
      acd323(12)=acd323(9)*acd323(10)
      acd323(13)=acd323(5)*acd323(3)
      acd323(12)=acd323(12)+acd323(13)
      acd323(12)=acd323(4)*acd323(12)
      acd323(13)=acd323(1)*acd323(2)
      brack=acd323(11)+acd323(12)+2.0_ki*acd323(13)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd323h2
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = -k2
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
   end function derivative
!---#] function derivative:
!---#[ subroutine reconstruct_d323:
   subroutine     reconstruct_d323(coeffs)
      use p16_bbarb_hepemg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_323:
      coeffs%coeffs_323%c0 = derivative(czip)
      coeffs%coeffs_323%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_323%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_323%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_323%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_323%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_323%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_323%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_323%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_323%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_323%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_323%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_323%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_323%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_323%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_323:
   end subroutine reconstruct_d323
!---#] subroutine reconstruct_d323:
end module     p16_bbarb_hepemg_d323h2l1d
