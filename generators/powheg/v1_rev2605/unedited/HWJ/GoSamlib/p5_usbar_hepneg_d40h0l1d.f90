module     p5_usbar_hepneg_d40h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p5_usbar_hepneg/helicity0d40h0l1d.f90
   ! generator: buildfortran_d.py
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d40
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(16) :: acd40
      complex(ki) :: brack
      acd40(1)=dotproduct(k2,qshift)
      acd40(2)=abb40(7)
      acd40(3)=dotproduct(k6,qshift)
      acd40(4)=abb40(9)
      acd40(5)=dotproduct(qshift,spvak5k1)
      acd40(6)=abb40(5)
      acd40(7)=dotproduct(qshift,spvak5k4)
      acd40(8)=abb40(10)
      acd40(9)=dotproduct(qshift,spvak6k2)
      acd40(10)=abb40(6)
      acd40(11)=abb40(8)
      acd40(12)=-acd40(9)*acd40(10)
      acd40(13)=-acd40(7)*acd40(8)
      acd40(14)=-acd40(5)*acd40(6)
      acd40(15)=-acd40(3)*acd40(4)
      acd40(16)=-acd40(1)*acd40(2)
      brack=acd40(11)+acd40(12)+acd40(13)+acd40(14)+acd40(15)+acd40(16)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(15) :: acd40
      complex(ki) :: brack
      acd40(1)=k2(iv1)
      acd40(2)=abb40(7)
      acd40(3)=k6(iv1)
      acd40(4)=abb40(9)
      acd40(5)=spvak5k1(iv1)
      acd40(6)=abb40(5)
      acd40(7)=spvak5k4(iv1)
      acd40(8)=abb40(10)
      acd40(9)=spvak6k2(iv1)
      acd40(10)=abb40(6)
      acd40(11)=acd40(9)*acd40(10)
      acd40(12)=acd40(7)*acd40(8)
      acd40(13)=acd40(5)*acd40(6)
      acd40(14)=acd40(3)*acd40(4)
      acd40(15)=acd40(1)*acd40(2)
      brack=acd40(11)+acd40(12)+acd40(13)+acd40(14)+acd40(15)
   end function brack_2
!---#] function brack_2:
!---#[ function derivative:
   function derivative(mu2,i1) result(numerator)
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = k6
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
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
   end function derivative
!---#] function derivative:
!---#[ subroutine reconstruct_d40:
   subroutine     reconstruct_d40(coeffs)
      use p5_usbar_hepneg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 1 case :
      !---[# reconstruct coeffs%coeffs_40:
      coeffs%coeffs_40%c0 = derivative(czip)
      coeffs%coeffs_40%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_40%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_40%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_40%c1(4,1) = -derivative(czip,4)
      !---#] reconstruct coeffs%coeffs_40:
   end subroutine reconstruct_d40
!---#] subroutine reconstruct_d40:
end module     p5_usbar_hepneg_d40h0l1d
