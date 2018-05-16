module     p12_cbbar_hepneg_d34h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d34h0l1d.f90
   ! generator: buildfortran_d.py
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d34
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(21) :: acd34
      complex(ki) :: brack
      acd34(1)=dotproduct(k6,qshift)
      acd34(2)=abb34(11)
      acd34(3)=dotproduct(qshift,qshift)
      acd34(4)=abb34(12)
      acd34(5)=dotproduct(qshift,spvak2k1)
      acd34(6)=abb34(6)
      acd34(7)=dotproduct(qshift,spvak5k1)
      acd34(8)=abb34(5)
      acd34(9)=dotproduct(qshift,spvak5k4)
      acd34(10)=dotproduct(qshift,spvak6k1)
      acd34(11)=abb34(14)
      acd34(12)=abb34(13)
      acd34(13)=dotproduct(qshift,spvak6k2)
      acd34(14)=abb34(8)
      acd34(15)=abb34(10)
      acd34(16)=-acd34(13)*acd34(14)
      acd34(17)=-acd34(7)*acd34(8)
      acd34(18)=-acd34(5)*acd34(6)
      acd34(19)=acd34(3)*acd34(4)
      acd34(20)=-acd34(1)*acd34(2)
      acd34(21)=acd34(9)*acd34(11)
      acd34(21)=-acd34(12)+acd34(21)
      acd34(21)=acd34(10)*acd34(21)
      brack=acd34(15)+acd34(16)+acd34(17)+acd34(18)+acd34(19)+acd34(20)+acd34(2&
      &1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(23) :: acd34
      complex(ki) :: brack
      acd34(1)=k6(iv1)
      acd34(2)=abb34(11)
      acd34(3)=qshift(iv1)
      acd34(4)=abb34(12)
      acd34(5)=spvak2k1(iv1)
      acd34(6)=abb34(6)
      acd34(7)=spvak5k1(iv1)
      acd34(8)=abb34(5)
      acd34(9)=spvak5k4(iv1)
      acd34(10)=dotproduct(qshift,spvak6k1)
      acd34(11)=abb34(14)
      acd34(12)=spvak6k1(iv1)
      acd34(13)=dotproduct(qshift,spvak5k4)
      acd34(14)=abb34(13)
      acd34(15)=spvak6k2(iv1)
      acd34(16)=abb34(8)
      acd34(17)=-acd34(9)*acd34(10)
      acd34(18)=-acd34(12)*acd34(13)
      acd34(17)=acd34(17)+acd34(18)
      acd34(17)=acd34(11)*acd34(17)
      acd34(18)=acd34(15)*acd34(16)
      acd34(19)=acd34(7)*acd34(8)
      acd34(20)=acd34(5)*acd34(6)
      acd34(21)=acd34(3)*acd34(4)
      acd34(22)=acd34(1)*acd34(2)
      acd34(23)=acd34(12)*acd34(14)
      brack=acd34(17)+acd34(18)+acd34(19)+acd34(20)-2.0_ki*acd34(21)+acd34(22)+&
      &acd34(23)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(9) :: acd34
      complex(ki) :: brack
      acd34(1)=d(iv1,iv2)
      acd34(2)=abb34(12)
      acd34(3)=spvak5k4(iv1)
      acd34(4)=spvak6k1(iv2)
      acd34(5)=abb34(14)
      acd34(6)=spvak5k4(iv2)
      acd34(7)=spvak6k1(iv1)
      acd34(8)=acd34(4)*acd34(3)
      acd34(9)=acd34(7)*acd34(6)
      acd34(8)=acd34(9)+acd34(8)
      acd34(8)=acd34(5)*acd34(8)
      acd34(9)=acd34(2)*acd34(1)
      brack=acd34(8)+2.0_ki*acd34(9)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd34h0
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
      qshift = k2
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
!---#[ subroutine reconstruct_d34:
   subroutine     reconstruct_d34(coeffs)
      use p12_cbbar_hepneg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_34:
      coeffs%coeffs_34%c0 = derivative(czip)
      coeffs%coeffs_34%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_34%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_34%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_34%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_34%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_34%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_34%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_34%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_34%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_34%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_34%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_34%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_34%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_34%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_34:
   end subroutine reconstruct_d34
!---#] subroutine reconstruct_d34:
end module     p12_cbbar_hepneg_d34h0l1d
