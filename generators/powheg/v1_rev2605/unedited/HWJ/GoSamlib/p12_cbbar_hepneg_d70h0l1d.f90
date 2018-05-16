module     p12_cbbar_hepneg_d70h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d70h0l1d.f90
   ! generator: buildfortran_d.py
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d70
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd70h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd70
      complex(ki) :: brack
      acd70(1)=abb70(8)
      brack=acd70(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd70h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd70
      complex(ki) :: brack
      acd70(1)=k6(iv1)
      acd70(2)=abb70(20)
      acd70(3)=e6(iv1)
      acd70(4)=abb70(14)
      acd70(5)=spvak2k6(iv1)
      acd70(6)=abb70(21)
      acd70(7)=spvak5k6(iv1)
      acd70(8)=abb70(22)
      acd70(9)=spvak6k1(iv1)
      acd70(10)=abb70(16)
      acd70(11)=spvak6k4(iv1)
      acd70(12)=abb70(18)
      acd70(13)=spvae6k1(iv1)
      acd70(14)=abb70(15)
      acd70(15)=spvak2e6(iv1)
      acd70(16)=abb70(13)
      acd70(17)=spvae6k2(iv1)
      acd70(18)=abb70(10)
      acd70(19)=spval3e6(iv1)
      acd70(20)=abb70(32)
      acd70(21)=spvae6l3(iv1)
      acd70(22)=abb70(25)
      acd70(23)=spvae6k4(iv1)
      acd70(24)=abb70(26)
      acd70(25)=spvak5e6(iv1)
      acd70(26)=abb70(28)
      acd70(27)=spvak6e6(iv1)
      acd70(28)=abb70(31)
      acd70(29)=spvae6k6(iv1)
      acd70(30)=abb70(29)
      acd70(31)=-acd70(2)*acd70(1)
      acd70(32)=-acd70(4)*acd70(3)
      acd70(33)=-acd70(6)*acd70(5)
      acd70(34)=-acd70(8)*acd70(7)
      acd70(35)=-acd70(10)*acd70(9)
      acd70(36)=-acd70(12)*acd70(11)
      acd70(37)=-acd70(14)*acd70(13)
      acd70(38)=-acd70(16)*acd70(15)
      acd70(39)=-acd70(18)*acd70(17)
      acd70(40)=-acd70(20)*acd70(19)
      acd70(41)=-acd70(22)*acd70(21)
      acd70(42)=-acd70(24)*acd70(23)
      acd70(43)=-acd70(26)*acd70(25)
      acd70(44)=-acd70(28)*acd70(27)
      acd70(45)=-acd70(30)*acd70(29)
      brack=acd70(31)+acd70(32)+acd70(33)+acd70(34)+acd70(35)+acd70(36)+acd70(3&
      &7)+acd70(38)+acd70(39)+acd70(40)+acd70(41)+acd70(42)+acd70(43)+acd70(44)&
      &+acd70(45)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd70h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(17) :: acd70
      complex(ki) :: brack
      acd70(1)=d(iv1,iv2)
      acd70(2)=abb70(17)
      acd70(3)=e6(iv1)
      acd70(4)=spvak2k1(iv2)
      acd70(5)=abb70(12)
      acd70(6)=spvak2k4(iv2)
      acd70(7)=abb70(9)
      acd70(8)=spvak5k1(iv2)
      acd70(9)=abb70(11)
      acd70(10)=e6(iv2)
      acd70(11)=spvak2k1(iv1)
      acd70(12)=spvak2k4(iv1)
      acd70(13)=spvak5k1(iv1)
      acd70(14)=acd70(9)*acd70(13)
      acd70(15)=acd70(7)*acd70(12)
      acd70(16)=acd70(5)*acd70(11)
      acd70(14)=acd70(16)+acd70(14)+acd70(15)
      acd70(14)=acd70(10)*acd70(14)
      acd70(15)=acd70(9)*acd70(8)
      acd70(16)=acd70(7)*acd70(6)
      acd70(17)=acd70(5)*acd70(4)
      acd70(15)=acd70(17)+acd70(15)+acd70(16)
      acd70(15)=acd70(3)*acd70(15)
      acd70(16)=acd70(1)*acd70(2)
      brack=acd70(14)+acd70(15)+2.0_ki*acd70(16)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd70h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd70
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd70h0
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
!---#[ subroutine reconstruct_d70:
   subroutine     reconstruct_d70(coeffs)
      use p12_cbbar_hepneg_groups, only: tensrec_info_group0
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group0), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_70:
      coeffs%coeffs_70%c0 = derivative(czip)
      coeffs%coeffs_70%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_70%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_70%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_70%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_70%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_70%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_70%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_70%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_70%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_70%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_70%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_70%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_70%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_70%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_70%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_70%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_70%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_70%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_70%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_70%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_70%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_70%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_70%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_70%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_70%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_70%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_70%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_70%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_70%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_70%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_70%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_70%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_70%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_70%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_70:
   end subroutine reconstruct_d70
!---#] subroutine reconstruct_d70:
end module     p12_cbbar_hepneg_d70h0l1d
