module     p0_dbard_hepemg_d413h1l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbard_hepemg/helicity1d413h1l1d.f90
   ! generator: buildfortran_d.py
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d413
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd413h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(19) :: acd413
      complex(ki) :: brack
      acd413(1)=dotproduct(k1,qshift)
      acd413(2)=abb413(15)
      acd413(3)=dotproduct(k6,qshift)
      acd413(4)=abb413(12)
      acd413(5)=dotproduct(qshift,spvak1k6)
      acd413(6)=abb413(9)
      acd413(7)=dotproduct(qshift,spvak2k4)
      acd413(8)=abb413(6)
      acd413(9)=dotproduct(qshift,spvak5k4)
      acd413(10)=abb413(8)
      acd413(11)=dotproduct(qshift,spvak6k1)
      acd413(12)=abb413(5)
      acd413(13)=abb413(7)
      acd413(14)=-acd413(2)*acd413(1)
      acd413(15)=-acd413(4)*acd413(3)
      acd413(16)=-acd413(6)*acd413(5)
      acd413(17)=-acd413(8)*acd413(7)
      acd413(18)=-acd413(10)*acd413(9)
      acd413(19)=-acd413(12)*acd413(11)
      brack=acd413(13)+acd413(14)+acd413(15)+acd413(16)+acd413(17)+acd413(18)+a&
      &cd413(19)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd413h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(18) :: acd413
      complex(ki) :: brack
      acd413(1)=k1(iv1)
      acd413(2)=abb413(15)
      acd413(3)=k6(iv1)
      acd413(4)=abb413(12)
      acd413(5)=spvak1k6(iv1)
      acd413(6)=abb413(9)
      acd413(7)=spvak2k4(iv1)
      acd413(8)=abb413(6)
      acd413(9)=spvak5k4(iv1)
      acd413(10)=abb413(8)
      acd413(11)=spvak6k1(iv1)
      acd413(12)=abb413(5)
      acd413(13)=-acd413(2)*acd413(1)
      acd413(14)=-acd413(4)*acd413(3)
      acd413(15)=-acd413(6)*acd413(5)
      acd413(16)=-acd413(8)*acd413(7)
      acd413(17)=-acd413(10)*acd413(9)
      acd413(18)=-acd413(12)*acd413(11)
      brack=acd413(13)+acd413(14)+acd413(15)+acd413(16)+acd413(17)+acd413(18)
   end function brack_2
!---#] function brack_2:
!---#[ function derivative:
   function derivative(mu2,i1) result(numerator)
      use p0_dbard_hepemg_globalsl1, only: epspow
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_abbrevd413h1
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = -k3-k5-k4
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
!---#[ subroutine reconstruct_d413:
   subroutine     reconstruct_d413(coeffs)
      use p0_dbard_hepemg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 1 case :
      !---[# reconstruct coeffs%coeffs_413:
      coeffs%coeffs_413%c0 = derivative(czip)
      coeffs%coeffs_413%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_413%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_413%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_413%c1(4,1) = -derivative(czip,4)
      !---#] reconstruct coeffs%coeffs_413:
   end subroutine reconstruct_d413
!---#] subroutine reconstruct_d413:
end module     p0_dbard_hepemg_d413h1l1d
