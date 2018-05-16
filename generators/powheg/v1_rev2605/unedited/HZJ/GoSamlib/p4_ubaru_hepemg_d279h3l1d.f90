module     p4_ubaru_hepemg_d279h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity3d279h3l1d.f90
   ! generator: buildfortran_d.py
   use p4_ubaru_hepemg_config, only: ki
   use p4_ubaru_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d279
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd279h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd279
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd279h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(36) :: acd279
      complex(ki) :: brack
      acd279(1)=k1(iv1)
      acd279(2)=abb279(42)
      acd279(3)=k2(iv1)
      acd279(4)=abb279(28)
      acd279(5)=k6(iv1)
      acd279(6)=abb279(31)
      acd279(7)=e6(iv1)
      acd279(8)=abb279(22)
      acd279(9)=spvak2k1(iv1)
      acd279(10)=abb279(12)
      acd279(11)=spvak2k5(iv1)
      acd279(12)=abb279(26)
      acd279(13)=spvak4k1(iv1)
      acd279(14)=abb279(46)
      acd279(15)=spvak4k5(iv1)
      acd279(16)=abb279(54)
      acd279(17)=spvak4k6(iv1)
      acd279(18)=abb279(52)
      acd279(19)=spvak6k5(iv1)
      acd279(20)=abb279(43)
      acd279(21)=spvae6k1(iv1)
      acd279(22)=abb279(48)
      acd279(23)=spvak2e6(iv1)
      acd279(24)=abb279(25)
      acd279(25)=-acd279(2)*acd279(1)
      acd279(26)=-acd279(4)*acd279(3)
      acd279(27)=-acd279(6)*acd279(5)
      acd279(28)=-acd279(8)*acd279(7)
      acd279(29)=-acd279(10)*acd279(9)
      acd279(30)=-acd279(12)*acd279(11)
      acd279(31)=-acd279(14)*acd279(13)
      acd279(32)=-acd279(16)*acd279(15)
      acd279(33)=-acd279(18)*acd279(17)
      acd279(34)=-acd279(20)*acd279(19)
      acd279(35)=-acd279(22)*acd279(21)
      acd279(36)=-acd279(24)*acd279(23)
      brack=acd279(25)+acd279(26)+acd279(27)+acd279(28)+acd279(29)+acd279(30)+a&
      &cd279(31)+acd279(32)+acd279(33)+acd279(34)+acd279(35)+acd279(36)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd279h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(75) :: acd279
      complex(ki) :: brack
      acd279(1)=k1(iv1)
      acd279(2)=e6(iv2)
      acd279(3)=abb279(40)
      acd279(4)=spvae6k1(iv2)
      acd279(5)=abb279(30)
      acd279(6)=spvak2e6(iv2)
      acd279(7)=abb279(34)
      acd279(8)=spvak4e6(iv2)
      acd279(9)=abb279(45)
      acd279(10)=spvae6k5(iv2)
      acd279(11)=abb279(18)
      acd279(12)=k1(iv2)
      acd279(13)=e6(iv1)
      acd279(14)=spvae6k1(iv1)
      acd279(15)=spvak2e6(iv1)
      acd279(16)=spvak4e6(iv1)
      acd279(17)=spvae6k5(iv1)
      acd279(18)=k2(iv1)
      acd279(19)=k2(iv2)
      acd279(20)=k6(iv1)
      acd279(21)=spvak2k1(iv2)
      acd279(22)=abb279(32)
      acd279(23)=spvak2k5(iv2)
      acd279(24)=abb279(55)
      acd279(25)=spvak4k1(iv2)
      acd279(26)=abb279(57)
      acd279(27)=k6(iv2)
      acd279(28)=spvak2k1(iv1)
      acd279(29)=spvak2k5(iv1)
      acd279(30)=spvak4k1(iv1)
      acd279(31)=abb279(27)
      acd279(32)=abb279(21)
      acd279(33)=abb279(23)
      acd279(34)=spvak2k6(iv2)
      acd279(35)=abb279(60)
      acd279(36)=spvak6k1(iv2)
      acd279(37)=abb279(50)
      acd279(38)=spvak2k6(iv1)
      acd279(39)=spvak6k1(iv1)
      acd279(40)=abb279(15)
      acd279(41)=spvak4k5(iv2)
      acd279(42)=abb279(37)
      acd279(43)=spvak4k5(iv1)
      acd279(44)=abb279(39)
      acd279(45)=abb279(53)
      acd279(46)=abb279(35)
      acd279(47)=abb279(19)
      acd279(48)=spvak4k6(iv2)
      acd279(49)=abb279(29)
      acd279(50)=spvak6k5(iv2)
      acd279(51)=abb279(17)
      acd279(52)=spvak1e6(iv2)
      acd279(53)=abb279(16)
      acd279(54)=spvae6k2(iv2)
      acd279(55)=abb279(20)
      acd279(56)=spvak4k6(iv1)
      acd279(57)=spvak6k5(iv1)
      acd279(58)=spvak1e6(iv1)
      acd279(59)=spvae6k2(iv1)
      acd279(60)=abb279(51)
      acd279(61)=abb279(56)
      acd279(62)=acd279(55)*acd279(54)
      acd279(63)=acd279(53)*acd279(52)
      acd279(64)=acd279(51)*acd279(50)
      acd279(65)=acd279(49)*acd279(48)
      acd279(66)=acd279(10)*acd279(47)
      acd279(67)=acd279(8)*acd279(46)
      acd279(68)=acd279(27)*acd279(22)
      acd279(69)=acd279(6)*acd279(44)
      acd279(70)=acd279(4)*acd279(40)
      acd279(71)=acd279(2)*acd279(31)
      acd279(62)=acd279(71)+acd279(70)+acd279(69)+acd279(68)+acd279(67)+acd279(&
      &66)+acd279(65)+acd279(64)+acd279(62)+acd279(63)
      acd279(62)=acd279(28)*acd279(62)
      acd279(63)=acd279(55)*acd279(59)
      acd279(64)=acd279(53)*acd279(58)
      acd279(65)=acd279(51)*acd279(57)
      acd279(66)=acd279(49)*acd279(56)
      acd279(67)=acd279(17)*acd279(47)
      acd279(68)=acd279(16)*acd279(46)
      acd279(69)=acd279(20)*acd279(22)
      acd279(70)=acd279(15)*acd279(44)
      acd279(71)=acd279(14)*acd279(40)
      acd279(72)=acd279(13)*acd279(31)
      acd279(63)=acd279(72)+acd279(71)+acd279(70)+acd279(69)+acd279(68)+acd279(&
      &67)+acd279(66)+acd279(65)+acd279(63)+acd279(64)
      acd279(63)=acd279(21)*acd279(63)
      acd279(64)=acd279(12)+acd279(19)
      acd279(65)=-acd279(3)*acd279(64)
      acd279(66)=acd279(36)*acd279(37)
      acd279(67)=acd279(34)*acd279(35)
      acd279(68)=acd279(25)*acd279(33)
      acd279(69)=acd279(23)*acd279(32)
      acd279(65)=acd279(69)+acd279(68)+acd279(66)+acd279(67)+acd279(65)
      acd279(65)=acd279(13)*acd279(65)
      acd279(66)=acd279(1)+acd279(18)
      acd279(67)=-acd279(3)*acd279(66)
      acd279(68)=acd279(37)*acd279(39)
      acd279(69)=acd279(35)*acd279(38)
      acd279(70)=acd279(30)*acd279(33)
      acd279(71)=acd279(29)*acd279(32)
      acd279(67)=acd279(71)+acd279(70)+acd279(68)+acd279(69)+acd279(67)
      acd279(67)=acd279(2)*acd279(67)
      acd279(68)=acd279(36)*acd279(61)
      acd279(69)=acd279(34)*acd279(60)
      acd279(70)=acd279(6)*acd279(45)
      acd279(71)=acd279(4)*acd279(42)
      acd279(68)=acd279(71)+acd279(70)+acd279(68)+acd279(69)
      acd279(68)=acd279(43)*acd279(68)
      acd279(69)=acd279(39)*acd279(61)
      acd279(70)=acd279(38)*acd279(60)
      acd279(71)=acd279(15)*acd279(45)
      acd279(72)=acd279(14)*acd279(42)
      acd279(69)=acd279(72)+acd279(71)+acd279(69)+acd279(70)
      acd279(69)=acd279(41)*acd279(69)
      acd279(70)=acd279(11)*acd279(17)
      acd279(71)=acd279(9)*acd279(16)
      acd279(70)=acd279(70)-acd279(71)
      acd279(70)=acd279(70)*acd279(64)
      acd279(71)=acd279(11)*acd279(10)
      acd279(72)=acd279(9)*acd279(8)
      acd279(71)=acd279(71)-acd279(72)
      acd279(71)=acd279(71)*acd279(66)
      acd279(72)=acd279(26)*acd279(30)
      acd279(73)=acd279(24)*acd279(29)
      acd279(72)=acd279(72)+acd279(73)
      acd279(72)=acd279(27)*acd279(72)
      acd279(73)=acd279(25)*acd279(26)
      acd279(74)=acd279(23)*acd279(24)
      acd279(73)=acd279(73)+acd279(74)
      acd279(73)=acd279(20)*acd279(73)
      acd279(64)=acd279(27)+acd279(64)
      acd279(74)=-acd279(15)*acd279(64)
      acd279(66)=acd279(20)+acd279(66)
      acd279(75)=-acd279(6)*acd279(66)
      acd279(74)=acd279(75)+acd279(74)
      acd279(74)=acd279(7)*acd279(74)
      acd279(64)=acd279(14)*acd279(64)
      acd279(66)=acd279(4)*acd279(66)
      acd279(64)=acd279(66)+acd279(64)
      acd279(64)=acd279(5)*acd279(64)
      brack=acd279(62)+acd279(63)+acd279(64)+acd279(65)+acd279(67)+acd279(68)+a&
      &cd279(69)+acd279(70)+acd279(71)+acd279(72)+acd279(73)+acd279(74)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd279h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd279
      complex(ki) :: brack
      acd279(1)=d(iv1,iv2)
      acd279(2)=e6(iv3)
      acd279(3)=abb279(36)
      acd279(4)=spvak2k1(iv3)
      acd279(5)=abb279(33)
      acd279(6)=spvak2k5(iv3)
      acd279(7)=abb279(14)
      acd279(8)=spvak4k1(iv3)
      acd279(9)=abb279(58)
      acd279(10)=spvae6k1(iv3)
      acd279(11)=abb279(30)
      acd279(12)=spvak2e6(iv3)
      acd279(13)=abb279(34)
      acd279(14)=spvak4e6(iv3)
      acd279(15)=abb279(24)
      acd279(16)=spvae6k5(iv3)
      acd279(17)=abb279(44)
      acd279(18)=d(iv1,iv3)
      acd279(19)=e6(iv2)
      acd279(20)=spvak2k1(iv2)
      acd279(21)=spvak2k5(iv2)
      acd279(22)=spvak4k1(iv2)
      acd279(23)=spvae6k1(iv2)
      acd279(24)=spvak2e6(iv2)
      acd279(25)=spvak4e6(iv2)
      acd279(26)=spvae6k5(iv2)
      acd279(27)=d(iv2,iv3)
      acd279(28)=e6(iv1)
      acd279(29)=spvak2k1(iv1)
      acd279(30)=spvak2k5(iv1)
      acd279(31)=spvak4k1(iv1)
      acd279(32)=spvae6k1(iv1)
      acd279(33)=spvak2e6(iv1)
      acd279(34)=spvak4e6(iv1)
      acd279(35)=spvae6k5(iv1)
      acd279(36)=-acd279(17)*acd279(35)
      acd279(37)=-acd279(15)*acd279(34)
      acd279(38)=acd279(13)*acd279(33)
      acd279(39)=-acd279(11)*acd279(32)
      acd279(40)=-acd279(9)*acd279(31)
      acd279(41)=-acd279(7)*acd279(30)
      acd279(42)=-acd279(5)*acd279(29)
      acd279(43)=-acd279(3)*acd279(28)
      acd279(36)=acd279(43)+acd279(42)+acd279(41)+acd279(40)+acd279(39)+acd279(&
      &38)+acd279(36)+acd279(37)
      acd279(36)=acd279(27)*acd279(36)
      acd279(37)=-acd279(17)*acd279(26)
      acd279(38)=-acd279(15)*acd279(25)
      acd279(39)=acd279(13)*acd279(24)
      acd279(40)=-acd279(11)*acd279(23)
      acd279(41)=-acd279(9)*acd279(22)
      acd279(42)=-acd279(7)*acd279(21)
      acd279(43)=-acd279(5)*acd279(20)
      acd279(44)=-acd279(3)*acd279(19)
      acd279(37)=acd279(44)+acd279(43)+acd279(42)+acd279(41)+acd279(40)+acd279(&
      &39)+acd279(37)+acd279(38)
      acd279(37)=acd279(18)*acd279(37)
      acd279(38)=-acd279(17)*acd279(16)
      acd279(39)=-acd279(15)*acd279(14)
      acd279(40)=acd279(13)*acd279(12)
      acd279(41)=-acd279(11)*acd279(10)
      acd279(42)=-acd279(9)*acd279(8)
      acd279(43)=-acd279(7)*acd279(6)
      acd279(44)=-acd279(5)*acd279(4)
      acd279(45)=-acd279(3)*acd279(2)
      acd279(38)=acd279(45)+acd279(44)+acd279(43)+acd279(42)+acd279(41)+acd279(&
      &40)+acd279(38)+acd279(39)
      acd279(38)=acd279(1)*acd279(38)
      acd279(36)=acd279(38)+acd279(36)+acd279(37)
      brack=2.0_ki*acd279(36)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd279h3
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
!---#[ subroutine reconstruct_d279:
   subroutine     reconstruct_d279(coeffs)
      use p4_ubaru_hepemg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_279:
      coeffs%coeffs_279%c0 = derivative(czip)
      coeffs%coeffs_279%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_279%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_279%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_279%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_279%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_279%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_279%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_279%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_279%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_279%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_279%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_279%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_279%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_279%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_279%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_279%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_279%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_279%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_279%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_279%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_279%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_279%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_279%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_279%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_279%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_279%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_279%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_279%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_279%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_279%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_279%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_279%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_279%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_279%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_279:
   end subroutine reconstruct_d279
!---#] subroutine reconstruct_d279:
end module     p4_ubaru_hepemg_d279h3l1d
