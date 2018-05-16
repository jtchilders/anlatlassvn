module     p8_cbarc_hepemg_d605h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p8_cbarc_hepemg/helicity3d605h3l1d.f90
   ! generator: buildfortran_d.py
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d605
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd605h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd605
      complex(ki) :: brack
      acd605(1)=abb605(13)
      brack=acd605(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd605h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(18) :: acd605
      complex(ki) :: brack
      acd605(1)=k2(iv1)
      acd605(2)=abb605(11)
      acd605(3)=k6(iv1)
      acd605(4)=abb605(14)
      acd605(5)=e6(iv1)
      acd605(6)=abb605(12)
      acd605(7)=spvak2k6(iv1)
      acd605(8)=abb605(7)
      acd605(9)=spvak6k2(iv1)
      acd605(10)=abb605(5)
      acd605(11)=spvak2e6(iv1)
      acd605(12)=abb605(18)
      acd605(13)=acd605(11)*acd605(12)
      acd605(14)=acd605(9)*acd605(10)
      acd605(15)=acd605(7)*acd605(8)
      acd605(16)=acd605(5)*acd605(6)
      acd605(17)=acd605(3)*acd605(4)
      acd605(18)=acd605(1)*acd605(2)
      brack=acd605(13)+acd605(14)+acd605(15)+acd605(16)+acd605(17)+acd605(18)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd605h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(13) :: acd605
      complex(ki) :: brack
      acd605(1)=d(iv1,iv2)
      acd605(2)=abb605(17)
      acd605(3)=k2(iv1)
      acd605(4)=e6(iv2)
      acd605(5)=abb605(16)
      acd605(6)=k2(iv2)
      acd605(7)=e6(iv1)
      acd605(8)=spvak2k6(iv2)
      acd605(9)=abb605(6)
      acd605(10)=spvak2k6(iv1)
      acd605(11)=acd605(9)*acd605(8)
      acd605(12)=acd605(5)*acd605(6)
      acd605(11)=acd605(12)+acd605(11)
      acd605(11)=acd605(7)*acd605(11)
      acd605(12)=acd605(9)*acd605(10)
      acd605(13)=acd605(5)*acd605(3)
      acd605(12)=acd605(12)+acd605(13)
      acd605(12)=acd605(4)*acd605(12)
      acd605(13)=acd605(1)*acd605(2)
      brack=acd605(11)+acd605(12)+2.0_ki*acd605(13)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd605h3
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
!---#[ subroutine reconstruct_d605:
   subroutine     reconstruct_d605(coeffs)
      use p8_cbarc_hepemg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_605:
      coeffs%coeffs_605%c0 = derivative(czip)
      coeffs%coeffs_605%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_605%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_605%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_605%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_605%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_605%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_605%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_605%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_605%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_605%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_605%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_605%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_605%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_605%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_605:
   end subroutine reconstruct_d605
!---#] subroutine reconstruct_d605:
end module     p8_cbarc_hepemg_d605h3l1d
