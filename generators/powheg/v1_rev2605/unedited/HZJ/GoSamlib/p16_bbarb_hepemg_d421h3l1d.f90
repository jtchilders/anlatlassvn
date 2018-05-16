module     p16_bbarb_hepemg_d421h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity3d421h3l1d.f90
   ! generator: buildfortran_d.py
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d421
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd421h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(16) :: acd421
      complex(ki) :: brack
      acd421(1)=dotproduct(k2,qshift)
      acd421(2)=abb421(8)
      acd421(3)=dotproduct(k6,qshift)
      acd421(4)=abb421(9)
      acd421(5)=dotproduct(e6,qshift)
      acd421(6)=abb421(6)
      acd421(7)=dotproduct(qshift,spvak6k2)
      acd421(8)=abb421(5)
      acd421(9)=dotproduct(qshift,spvak2e6)
      acd421(10)=abb421(7)
      acd421(11)=abb421(10)
      acd421(12)=-acd421(9)*acd421(10)
      acd421(13)=-acd421(7)*acd421(8)
      acd421(14)=-acd421(5)*acd421(6)
      acd421(15)=-acd421(3)*acd421(4)
      acd421(16)=-acd421(1)*acd421(2)
      brack=acd421(11)+acd421(12)+acd421(13)+acd421(14)+acd421(15)+acd421(16)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd421h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(15) :: acd421
      complex(ki) :: brack
      acd421(1)=k2(iv1)
      acd421(2)=abb421(8)
      acd421(3)=k6(iv1)
      acd421(4)=abb421(9)
      acd421(5)=e6(iv1)
      acd421(6)=abb421(6)
      acd421(7)=spvak6k2(iv1)
      acd421(8)=abb421(5)
      acd421(9)=spvak2e6(iv1)
      acd421(10)=abb421(7)
      acd421(11)=acd421(2)*acd421(1)
      acd421(12)=acd421(4)*acd421(3)
      acd421(13)=acd421(6)*acd421(5)
      acd421(14)=acd421(8)*acd421(7)
      acd421(15)=acd421(10)*acd421(9)
      brack=acd421(11)+acd421(12)+acd421(13)+acd421(14)+acd421(15)
   end function brack_2
!---#] function brack_2:
!---#[ function derivative:
   function derivative(mu2,i1) result(numerator)
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd421h3
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
!---#[ subroutine reconstruct_d421:
   subroutine     reconstruct_d421(coeffs)
      use p16_bbarb_hepemg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 1 case :
      !---[# reconstruct coeffs%coeffs_421:
      coeffs%coeffs_421%c0 = derivative(czip)
      coeffs%coeffs_421%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_421%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_421%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_421%c1(4,1) = -derivative(czip,4)
      !---#] reconstruct coeffs%coeffs_421:
   end subroutine reconstruct_d421
!---#] subroutine reconstruct_d421:
end module     p16_bbarb_hepemg_d421h3l1d
