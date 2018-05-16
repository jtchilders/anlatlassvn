module     p12_cbbar_hepneg_d94h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d94h0l1d.f90
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
   public :: derivative , reconstruct_d94
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(24) :: acd94
      complex(ki) :: brack
      acd94(1)=dotproduct(k1,qshift)
      acd94(2)=dotproduct(e6,qshift)
      acd94(3)=abb94(17)
      acd94(4)=abb94(7)
      acd94(5)=dotproduct(k6,qshift)
      acd94(6)=abb94(8)
      acd94(7)=dotproduct(qshift,spvak6k1)
      acd94(8)=abb94(9)
      acd94(9)=abb94(11)
      acd94(10)=dotproduct(qshift,qshift)
      acd94(11)=abb94(10)
      acd94(12)=abb94(6)
      acd94(13)=dotproduct(qshift,spvak1k6)
      acd94(14)=abb94(12)
      acd94(15)=dotproduct(qshift,spvae6k1)
      acd94(16)=abb94(15)
      acd94(17)=abb94(5)
      acd94(18)=acd94(7)*acd94(8)
      acd94(19)=acd94(1)*acd94(3)
      acd94(18)=acd94(19)-acd94(9)+acd94(18)
      acd94(18)=acd94(2)*acd94(18)
      acd94(19)=-acd94(15)*acd94(16)
      acd94(20)=-acd94(13)*acd94(14)
      acd94(21)=acd94(10)*acd94(11)
      acd94(22)=-acd94(5)*acd94(6)
      acd94(23)=-acd94(7)*acd94(12)
      acd94(24)=-acd94(1)*acd94(4)
      brack=acd94(17)+acd94(18)+acd94(19)+acd94(20)+acd94(21)+acd94(22)+acd94(2&
      &3)+acd94(24)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(26) :: acd94
      complex(ki) :: brack
      acd94(1)=k1(iv1)
      acd94(2)=dotproduct(e6,qshift)
      acd94(3)=abb94(17)
      acd94(4)=abb94(7)
      acd94(5)=k6(iv1)
      acd94(6)=abb94(8)
      acd94(7)=e6(iv1)
      acd94(8)=dotproduct(k1,qshift)
      acd94(9)=dotproduct(qshift,spvak6k1)
      acd94(10)=abb94(9)
      acd94(11)=abb94(11)
      acd94(12)=qshift(iv1)
      acd94(13)=abb94(10)
      acd94(14)=spvak6k1(iv1)
      acd94(15)=abb94(6)
      acd94(16)=spvak1k6(iv1)
      acd94(17)=abb94(12)
      acd94(18)=spvae6k1(iv1)
      acd94(19)=abb94(15)
      acd94(20)=acd94(10)*acd94(9)
      acd94(21)=acd94(3)*acd94(8)
      acd94(20)=acd94(21)-acd94(11)+acd94(20)
      acd94(20)=acd94(7)*acd94(20)
      acd94(21)=acd94(2)*acd94(10)
      acd94(21)=acd94(21)-acd94(15)
      acd94(21)=acd94(14)*acd94(21)
      acd94(22)=-acd94(18)*acd94(19)
      acd94(23)=-acd94(16)*acd94(17)
      acd94(24)=acd94(12)*acd94(13)
      acd94(25)=-acd94(5)*acd94(6)
      acd94(26)=acd94(2)*acd94(3)
      acd94(26)=-acd94(4)+acd94(26)
      acd94(26)=acd94(1)*acd94(26)
      brack=acd94(20)+acd94(21)+acd94(22)+acd94(23)+2.0_ki*acd94(24)+acd94(25)+&
      &acd94(26)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(13) :: acd94
      complex(ki) :: brack
      acd94(1)=d(iv1,iv2)
      acd94(2)=abb94(10)
      acd94(3)=k1(iv1)
      acd94(4)=e6(iv2)
      acd94(5)=abb94(17)
      acd94(6)=k1(iv2)
      acd94(7)=e6(iv1)
      acd94(8)=spvak6k1(iv2)
      acd94(9)=abb94(9)
      acd94(10)=spvak6k1(iv1)
      acd94(11)=acd94(9)*acd94(8)
      acd94(12)=acd94(5)*acd94(6)
      acd94(11)=acd94(12)+acd94(11)
      acd94(11)=acd94(7)*acd94(11)
      acd94(12)=acd94(9)*acd94(10)
      acd94(13)=acd94(5)*acd94(3)
      acd94(12)=acd94(12)+acd94(13)
      acd94(12)=acd94(4)*acd94(12)
      acd94(13)=acd94(1)*acd94(2)
      brack=acd94(11)+acd94(12)+2.0_ki*acd94(13)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd94h0
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
      qshift = -k3-k6-k5-k4
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
!---#[ subroutine reconstruct_d94:
   subroutine     reconstruct_d94(coeffs)
      use p12_cbbar_hepneg_groups, only: tensrec_info_group3
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group3), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_94:
      coeffs%coeffs_94%c0 = derivative(czip)
      coeffs%coeffs_94%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_94%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_94%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_94%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_94%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_94%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_94%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_94%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_94%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_94%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_94%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_94%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_94%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_94%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_94:
   end subroutine reconstruct_d94
!---#] subroutine reconstruct_d94:
end module     p12_cbbar_hepneg_d94h0l1d
