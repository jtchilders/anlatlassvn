module     p12_cbbar_hepneg_d37h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d37h0l1d.f90
   ! generator: buildfortran_d.py
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d37
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd37h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(19) :: acd37
      complex(ki) :: brack
      acd37(1)=dotproduct(k1,qshift)
      acd37(2)=abb37(5)
      acd37(3)=dotproduct(k6,qshift)
      acd37(4)=abb37(13)
      acd37(5)=dotproduct(e6,qshift)
      acd37(6)=abb37(11)
      acd37(7)=dotproduct(qshift,spvak1k6)
      acd37(8)=abb37(9)
      acd37(9)=dotproduct(qshift,spvak6k1)
      acd37(10)=abb37(7)
      acd37(11)=dotproduct(qshift,spvae6k1)
      acd37(12)=abb37(6)
      acd37(13)=abb37(8)
      acd37(14)=-acd37(11)*acd37(12)
      acd37(15)=-acd37(9)*acd37(10)
      acd37(16)=-acd37(7)*acd37(8)
      acd37(17)=-acd37(5)*acd37(6)
      acd37(18)=-acd37(3)*acd37(4)
      acd37(19)=-acd37(1)*acd37(2)
      brack=acd37(13)+acd37(14)+acd37(15)+acd37(16)+acd37(17)+acd37(18)+acd37(1&
      &9)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd37h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(18) :: acd37
      complex(ki) :: brack
      acd37(1)=k1(iv1)
      acd37(2)=abb37(5)
      acd37(3)=k6(iv1)
      acd37(4)=abb37(13)
      acd37(5)=e6(iv1)
      acd37(6)=abb37(11)
      acd37(7)=spvak1k6(iv1)
      acd37(8)=abb37(9)
      acd37(9)=spvak6k1(iv1)
      acd37(10)=abb37(7)
      acd37(11)=spvae6k1(iv1)
      acd37(12)=abb37(6)
      acd37(13)=-acd37(11)*acd37(12)
      acd37(14)=-acd37(9)*acd37(10)
      acd37(15)=-acd37(7)*acd37(8)
      acd37(16)=-acd37(5)*acd37(6)
      acd37(17)=-acd37(3)*acd37(4)
      acd37(18)=-acd37(1)*acd37(2)
      brack=acd37(13)+acd37(14)+acd37(15)+acd37(16)+acd37(17)+acd37(18)
   end function brack_2
!---#] function brack_2:
!---#[ function derivative:
   function derivative(mu2,i1) result(numerator)
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd37h0
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
!---#[ subroutine reconstruct_d37:
   subroutine     reconstruct_d37(coeffs)
      use p12_cbbar_hepneg_groups, only: tensrec_info_group3
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group3), intent(out) :: coeffs
      ! rank 1 case :
      !---[# reconstruct coeffs%coeffs_37:
      coeffs%coeffs_37%c0 = derivative(czip)
      coeffs%coeffs_37%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_37%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_37%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_37%c1(4,1) = -derivative(czip,4)
      !---#] reconstruct coeffs%coeffs_37:
   end subroutine reconstruct_d37
!---#] subroutine reconstruct_d37:
end module     p12_cbbar_hepneg_d37h0l1d
