module     p4_ubaru_hepemg_d533h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity3d533h3l1d.f90
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
   public :: derivative , reconstruct_d533
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(48) :: acd533
      complex(ki) :: brack
      acd533(1)=dotproduct(k1,qshift)
      acd533(2)=dotproduct(e6,qshift)
      acd533(3)=abb533(6)
      acd533(4)=dotproduct(qshift,spvak2k5)
      acd533(5)=abb533(22)
      acd533(6)=dotproduct(qshift,spvak4k5)
      acd533(7)=abb533(25)
      acd533(8)=dotproduct(qshift,spvae6k1)
      acd533(9)=abb533(21)
      acd533(10)=dotproduct(qshift,spvak2e6)
      acd533(11)=abb533(19)
      acd533(12)=abb533(4)
      acd533(13)=dotproduct(k6,qshift)
      acd533(14)=abb533(16)
      acd533(15)=dotproduct(qshift,qshift)
      acd533(16)=abb533(18)
      acd533(17)=abb533(13)
      acd533(18)=dotproduct(qshift,spvak2k1)
      acd533(19)=abb533(11)
      acd533(20)=abb533(8)
      acd533(21)=dotproduct(qshift,spvak2k6)
      acd533(22)=abb533(30)
      acd533(23)=dotproduct(qshift,spvak6k1)
      acd533(24)=abb533(12)
      acd533(25)=abb533(10)
      acd533(26)=abb533(28)
      acd533(27)=abb533(14)
      acd533(28)=abb533(20)
      acd533(29)=abb533(15)
      acd533(30)=abb533(7)
      acd533(31)=abb533(27)
      acd533(32)=abb533(26)
      acd533(33)=abb533(17)
      acd533(34)=abb533(31)
      acd533(35)=abb533(23)
      acd533(36)=abb533(5)
      acd533(37)=abb533(29)
      acd533(38)=abb533(9)
      acd533(39)=acd533(3)*acd533(1)
      acd533(40)=-acd533(16)*acd533(15)
      acd533(41)=acd533(17)*acd533(4)
      acd533(42)=acd533(18)*acd533(6)
      acd533(43)=-acd533(19)*acd533(42)
      acd533(44)=acd533(22)*acd533(21)
      acd533(45)=acd533(24)*acd533(23)
      acd533(39)=-acd533(25)+acd533(45)+acd533(44)+acd533(43)+acd533(41)+acd533&
      &(40)+acd533(39)
      acd533(39)=acd533(2)*acd533(39)
      acd533(40)=acd533(7)*acd533(1)
      acd533(41)=acd533(31)*acd533(8)
      acd533(43)=acd533(32)*acd533(10)
      acd533(44)=-acd533(34)*acd533(21)
      acd533(40)=-acd533(35)+acd533(44)+acd533(43)+acd533(41)+acd533(40)
      acd533(40)=acd533(6)*acd533(40)
      acd533(41)=-acd533(5)*acd533(4)
      acd533(43)=-acd533(9)*acd533(8)
      acd533(44)=acd533(11)*acd533(10)
      acd533(41)=acd533(44)+acd533(41)+acd533(43)
      acd533(43)=acd533(13)-acd533(1)
      acd533(41)=acd533(43)*acd533(41)
      acd533(43)=-acd533(26)*acd533(4)
      acd533(44)=-acd533(27)*acd533(8)
      acd533(45)=-acd533(28)*acd533(10)
      acd533(43)=acd533(29)+acd533(45)+acd533(44)+acd533(43)
      acd533(43)=acd533(15)*acd533(43)
      acd533(44)=acd533(20)*acd533(2)
      acd533(44)=-acd533(36)+acd533(44)
      acd533(44)=acd533(18)*acd533(44)
      acd533(45)=-acd533(12)*acd533(1)
      acd533(46)=-acd533(14)*acd533(13)
      acd533(47)=-acd533(30)*acd533(4)
      acd533(42)=acd533(33)*acd533(42)
      acd533(48)=-acd533(37)*acd533(21)
      brack=acd533(38)+acd533(39)+acd533(40)+acd533(41)+acd533(42)+acd533(43)+a&
      &cd533(44)+acd533(45)+acd533(46)+acd533(47)+acd533(48)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(61) :: acd533
      complex(ki) :: brack
      acd533(1)=k1(iv1)
      acd533(2)=dotproduct(e6,qshift)
      acd533(3)=abb533(6)
      acd533(4)=dotproduct(qshift,spvak2k5)
      acd533(5)=abb533(22)
      acd533(6)=dotproduct(qshift,spvak4k5)
      acd533(7)=abb533(25)
      acd533(8)=dotproduct(qshift,spvae6k1)
      acd533(9)=abb533(21)
      acd533(10)=dotproduct(qshift,spvak2e6)
      acd533(11)=abb533(19)
      acd533(12)=abb533(4)
      acd533(13)=k6(iv1)
      acd533(14)=abb533(16)
      acd533(15)=e6(iv1)
      acd533(16)=dotproduct(k1,qshift)
      acd533(17)=dotproduct(qshift,qshift)
      acd533(18)=abb533(18)
      acd533(19)=abb533(13)
      acd533(20)=dotproduct(qshift,spvak2k1)
      acd533(21)=abb533(11)
      acd533(22)=abb533(8)
      acd533(23)=dotproduct(qshift,spvak2k6)
      acd533(24)=abb533(30)
      acd533(25)=dotproduct(qshift,spvak6k1)
      acd533(26)=abb533(12)
      acd533(27)=abb533(10)
      acd533(28)=qshift(iv1)
      acd533(29)=abb533(28)
      acd533(30)=abb533(14)
      acd533(31)=abb533(20)
      acd533(32)=abb533(15)
      acd533(33)=spvak2k5(iv1)
      acd533(34)=dotproduct(k6,qshift)
      acd533(35)=abb533(7)
      acd533(36)=spvak4k5(iv1)
      acd533(37)=abb533(27)
      acd533(38)=abb533(26)
      acd533(39)=abb533(17)
      acd533(40)=abb533(31)
      acd533(41)=abb533(23)
      acd533(42)=spvae6k1(iv1)
      acd533(43)=spvak2e6(iv1)
      acd533(44)=spvak2k1(iv1)
      acd533(45)=abb533(5)
      acd533(46)=spvak2k6(iv1)
      acd533(47)=abb533(29)
      acd533(48)=spvak6k1(iv1)
      acd533(49)=acd533(6)*acd533(21)
      acd533(49)=acd533(49)-acd533(22)
      acd533(50)=acd533(44)*acd533(49)
      acd533(51)=-acd533(26)*acd533(48)
      acd533(52)=-acd533(46)*acd533(24)
      acd533(53)=-acd533(33)*acd533(19)
      acd533(54)=2.0_ki*acd533(28)
      acd533(55)=acd533(18)*acd533(54)
      acd533(56)=-acd533(1)*acd533(3)
      acd533(57)=acd533(36)*acd533(20)*acd533(21)
      acd533(50)=acd533(57)+acd533(56)+acd533(55)+acd533(53)+acd533(51)+acd533(&
      &52)+acd533(50)
      acd533(50)=acd533(2)*acd533(50)
      acd533(49)=acd533(20)*acd533(49)
      acd533(51)=-acd533(26)*acd533(25)
      acd533(52)=-acd533(23)*acd533(24)
      acd533(53)=acd533(17)*acd533(18)
      acd533(55)=-acd533(4)*acd533(19)
      acd533(56)=-acd533(16)*acd533(3)
      acd533(49)=acd533(56)+acd533(55)+acd533(53)+acd533(52)+acd533(27)+acd533(&
      &51)+acd533(49)
      acd533(49)=acd533(15)*acd533(49)
      acd533(51)=acd533(23)*acd533(40)
      acd533(52)=-acd533(20)*acd533(39)
      acd533(53)=-acd533(10)*acd533(38)
      acd533(55)=-acd533(8)*acd533(37)
      acd533(56)=-acd533(16)*acd533(7)
      acd533(51)=acd533(56)+acd533(55)+acd533(53)+acd533(52)+acd533(41)+acd533(&
      &51)
      acd533(51)=acd533(36)*acd533(51)
      acd533(52)=acd533(46)*acd533(40)
      acd533(53)=-acd533(44)*acd533(39)
      acd533(55)=-acd533(43)*acd533(38)
      acd533(56)=-acd533(42)*acd533(37)
      acd533(57)=-acd533(1)*acd533(7)
      acd533(52)=acd533(57)+acd533(56)+acd533(55)+acd533(52)+acd533(53)
      acd533(52)=acd533(6)*acd533(52)
      acd533(53)=acd533(10)*acd533(11)
      acd533(55)=acd533(8)*acd533(9)
      acd533(56)=acd533(4)*acd533(5)
      acd533(53)=-acd533(56)+acd533(53)-acd533(55)
      acd533(55)=acd533(14)-acd533(53)
      acd533(55)=acd533(13)*acd533(55)
      acd533(56)=acd533(10)*acd533(31)
      acd533(57)=acd533(8)*acd533(30)
      acd533(58)=acd533(4)*acd533(29)
      acd533(56)=acd533(58)+acd533(57)-acd533(32)+acd533(56)
      acd533(54)=acd533(56)*acd533(54)
      acd533(53)=acd533(12)+acd533(53)
      acd533(53)=acd533(1)*acd533(53)
      acd533(56)=acd533(11)*acd533(43)
      acd533(57)=acd533(9)*acd533(42)
      acd533(56)=acd533(56)-acd533(57)
      acd533(57)=-acd533(34)*acd533(56)
      acd533(58)=acd533(43)*acd533(31)
      acd533(59)=acd533(42)*acd533(30)
      acd533(58)=acd533(58)+acd533(59)
      acd533(58)=acd533(17)*acd533(58)
      acd533(59)=acd533(17)*acd533(29)
      acd533(60)=acd533(5)*acd533(34)
      acd533(59)=acd533(60)+acd533(35)+acd533(59)
      acd533(59)=acd533(33)*acd533(59)
      acd533(60)=-acd533(33)*acd533(5)
      acd533(56)=acd533(60)+acd533(56)
      acd533(56)=acd533(16)*acd533(56)
      acd533(60)=acd533(46)*acd533(47)
      acd533(61)=acd533(44)*acd533(45)
      brack=acd533(49)+acd533(50)+acd533(51)+acd533(52)+acd533(53)+acd533(54)+a&
      &cd533(55)+acd533(56)+acd533(57)+acd533(58)+acd533(59)+acd533(60)+acd533(&
      &61)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(60) :: acd533
      complex(ki) :: brack
      acd533(1)=d(iv1,iv2)
      acd533(2)=dotproduct(e6,qshift)
      acd533(3)=abb533(18)
      acd533(4)=dotproduct(qshift,spvak2k5)
      acd533(5)=abb533(28)
      acd533(6)=dotproduct(qshift,spvae6k1)
      acd533(7)=abb533(14)
      acd533(8)=dotproduct(qshift,spvak2e6)
      acd533(9)=abb533(20)
      acd533(10)=abb533(15)
      acd533(11)=k1(iv1)
      acd533(12)=e6(iv2)
      acd533(13)=abb533(6)
      acd533(14)=spvak2k5(iv2)
      acd533(15)=abb533(22)
      acd533(16)=spvae6k1(iv2)
      acd533(17)=abb533(21)
      acd533(18)=spvak2e6(iv2)
      acd533(19)=abb533(19)
      acd533(20)=spvak4k5(iv2)
      acd533(21)=abb533(25)
      acd533(22)=k1(iv2)
      acd533(23)=e6(iv1)
      acd533(24)=spvak2k5(iv1)
      acd533(25)=spvae6k1(iv1)
      acd533(26)=spvak2e6(iv1)
      acd533(27)=spvak4k5(iv1)
      acd533(28)=k6(iv1)
      acd533(29)=k6(iv2)
      acd533(30)=qshift(iv2)
      acd533(31)=abb533(13)
      acd533(32)=dotproduct(qshift,spvak2k1)
      acd533(33)=abb533(11)
      acd533(34)=spvak2k1(iv2)
      acd533(35)=dotproduct(qshift,spvak4k5)
      acd533(36)=abb533(8)
      acd533(37)=spvak2k6(iv2)
      acd533(38)=abb533(30)
      acd533(39)=spvak6k1(iv2)
      acd533(40)=abb533(12)
      acd533(41)=qshift(iv1)
      acd533(42)=spvak2k1(iv1)
      acd533(43)=spvak2k6(iv1)
      acd533(44)=spvak6k1(iv1)
      acd533(45)=abb533(27)
      acd533(46)=abb533(26)
      acd533(47)=abb533(17)
      acd533(48)=abb533(31)
      acd533(49)=acd533(40)*acd533(39)
      acd533(50)=acd533(37)*acd533(38)
      acd533(51)=acd533(34)*acd533(36)
      acd533(52)=2.0_ki*acd533(3)
      acd533(52)=-acd533(30)*acd533(52)
      acd533(53)=acd533(14)*acd533(31)
      acd533(54)=acd533(22)*acd533(13)
      acd533(55)=acd533(33)*acd533(34)
      acd533(56)=-acd533(35)*acd533(55)
      acd533(57)=acd533(33)*acd533(32)
      acd533(58)=-acd533(20)*acd533(57)
      acd533(49)=acd533(58)+acd533(56)+acd533(54)+acd533(53)+acd533(52)+acd533(&
      &51)+acd533(49)+acd533(50)
      acd533(49)=acd533(23)*acd533(49)
      acd533(50)=acd533(40)*acd533(44)
      acd533(51)=acd533(38)*acd533(43)
      acd533(52)=acd533(42)*acd533(36)
      acd533(53)=2.0_ki*acd533(41)
      acd533(54)=-acd533(3)*acd533(53)
      acd533(56)=acd533(24)*acd533(31)
      acd533(58)=acd533(11)*acd533(13)
      acd533(59)=acd533(33)*acd533(42)
      acd533(60)=-acd533(35)*acd533(59)
      acd533(57)=-acd533(27)*acd533(57)
      acd533(50)=acd533(57)+acd533(60)+acd533(58)+acd533(56)+acd533(54)+acd533(&
      &52)+acd533(50)+acd533(51)
      acd533(50)=acd533(12)*acd533(50)
      acd533(51)=-acd533(37)*acd533(48)
      acd533(52)=acd533(34)*acd533(47)
      acd533(54)=acd533(18)*acd533(46)
      acd533(56)=acd533(16)*acd533(45)
      acd533(57)=acd533(22)*acd533(21)
      acd533(55)=-acd533(2)*acd533(55)
      acd533(51)=acd533(55)+acd533(57)+acd533(56)+acd533(54)+acd533(51)+acd533(&
      &52)
      acd533(51)=acd533(27)*acd533(51)
      acd533(52)=-acd533(43)*acd533(48)
      acd533(54)=acd533(42)*acd533(47)
      acd533(55)=acd533(26)*acd533(46)
      acd533(56)=acd533(25)*acd533(45)
      acd533(57)=acd533(11)*acd533(21)
      acd533(58)=-acd533(2)*acd533(59)
      acd533(52)=acd533(58)+acd533(57)+acd533(56)+acd533(55)+acd533(52)+acd533(&
      &54)
      acd533(52)=acd533(20)*acd533(52)
      acd533(54)=-acd533(9)*acd533(8)
      acd533(55)=-acd533(7)*acd533(6)
      acd533(56)=-acd533(5)*acd533(4)
      acd533(57)=-acd533(2)*acd533(3)
      acd533(54)=acd533(57)+acd533(56)+acd533(55)+acd533(10)+acd533(54)
      acd533(54)=acd533(1)*acd533(54)
      acd533(55)=-acd533(26)*acd533(9)
      acd533(56)=-acd533(25)*acd533(7)
      acd533(57)=-acd533(24)*acd533(5)
      acd533(55)=acd533(57)+acd533(55)+acd533(56)
      acd533(55)=acd533(30)*acd533(55)
      acd533(54)=acd533(55)+acd533(54)
      acd533(55)=acd533(19)*acd533(26)
      acd533(56)=acd533(17)*acd533(25)
      acd533(57)=acd533(15)*acd533(24)
      acd533(55)=-acd533(57)+acd533(55)-acd533(56)
      acd533(56)=-acd533(22)+acd533(29)
      acd533(55)=acd533(55)*acd533(56)
      acd533(56)=-acd533(18)*acd533(19)
      acd533(57)=acd533(16)*acd533(17)
      acd533(58)=acd533(14)*acd533(15)
      acd533(56)=acd533(58)+acd533(56)+acd533(57)
      acd533(56)=acd533(11)*acd533(56)
      acd533(57)=-acd533(9)*acd533(53)
      acd533(58)=acd533(19)*acd533(28)
      acd533(57)=acd533(57)+acd533(58)
      acd533(57)=acd533(18)*acd533(57)
      acd533(58)=-acd533(7)*acd533(53)
      acd533(59)=-acd533(17)*acd533(28)
      acd533(58)=acd533(58)+acd533(59)
      acd533(58)=acd533(16)*acd533(58)
      acd533(53)=-acd533(5)*acd533(53)
      acd533(59)=-acd533(15)*acd533(28)
      acd533(53)=acd533(53)+acd533(59)
      acd533(53)=acd533(14)*acd533(53)
      brack=acd533(49)+acd533(50)+acd533(51)+acd533(52)+acd533(53)+2.0_ki*acd53&
      &3(54)+acd533(55)+acd533(56)+acd533(57)+acd533(58)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(32) :: acd533
      complex(ki) :: brack
      acd533(1)=d(iv1,iv2)
      acd533(2)=e6(iv3)
      acd533(3)=abb533(18)
      acd533(4)=spvak2k5(iv3)
      acd533(5)=abb533(28)
      acd533(6)=spvae6k1(iv3)
      acd533(7)=abb533(14)
      acd533(8)=spvak2e6(iv3)
      acd533(9)=abb533(20)
      acd533(10)=d(iv1,iv3)
      acd533(11)=e6(iv2)
      acd533(12)=spvak2k5(iv2)
      acd533(13)=spvae6k1(iv2)
      acd533(14)=spvak2e6(iv2)
      acd533(15)=d(iv2,iv3)
      acd533(16)=e6(iv1)
      acd533(17)=spvak2k5(iv1)
      acd533(18)=spvae6k1(iv1)
      acd533(19)=spvak2e6(iv1)
      acd533(20)=spvak2k1(iv2)
      acd533(21)=spvak4k5(iv3)
      acd533(22)=abb533(11)
      acd533(23)=spvak2k1(iv3)
      acd533(24)=spvak4k5(iv2)
      acd533(25)=spvak2k1(iv1)
      acd533(26)=spvak4k5(iv1)
      acd533(27)=acd533(9)*acd533(19)
      acd533(28)=acd533(7)*acd533(18)
      acd533(29)=acd533(5)*acd533(17)
      acd533(30)=acd533(3)*acd533(16)
      acd533(27)=acd533(30)+acd533(29)+acd533(27)+acd533(28)
      acd533(27)=acd533(15)*acd533(27)
      acd533(28)=acd533(9)*acd533(14)
      acd533(29)=acd533(7)*acd533(13)
      acd533(30)=acd533(5)*acd533(12)
      acd533(31)=acd533(3)*acd533(11)
      acd533(28)=acd533(31)+acd533(30)+acd533(28)+acd533(29)
      acd533(28)=acd533(10)*acd533(28)
      acd533(29)=acd533(9)*acd533(8)
      acd533(30)=acd533(7)*acd533(6)
      acd533(31)=acd533(5)*acd533(4)
      acd533(32)=acd533(2)*acd533(3)
      acd533(29)=acd533(32)+acd533(31)+acd533(29)+acd533(30)
      acd533(29)=acd533(1)*acd533(29)
      acd533(27)=acd533(29)+acd533(27)+acd533(28)
      acd533(28)=acd533(23)*acd533(24)
      acd533(29)=acd533(20)*acd533(21)
      acd533(28)=acd533(28)+acd533(29)
      acd533(28)=acd533(16)*acd533(28)
      acd533(29)=acd533(23)*acd533(26)
      acd533(30)=acd533(21)*acd533(25)
      acd533(29)=acd533(29)+acd533(30)
      acd533(29)=acd533(11)*acd533(29)
      acd533(30)=acd533(24)*acd533(25)
      acd533(31)=acd533(20)*acd533(26)
      acd533(30)=acd533(30)+acd533(31)
      acd533(30)=acd533(2)*acd533(30)
      acd533(28)=acd533(30)+acd533(28)+acd533(29)
      acd533(28)=acd533(22)*acd533(28)
      brack=2.0_ki*acd533(27)+acd533(28)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd533h3
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
!---#[ subroutine reconstruct_d533:
   subroutine     reconstruct_d533(coeffs)
      use p4_ubaru_hepemg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_533:
      coeffs%coeffs_533%c0 = derivative(czip)
      coeffs%coeffs_533%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_533%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_533%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_533%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_533%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_533%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_533%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_533%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_533%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_533%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_533%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_533%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_533%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_533%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_533%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_533%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_533%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_533%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_533%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_533%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_533%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_533%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_533%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_533%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_533%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_533%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_533%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_533%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_533%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_533%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_533%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_533%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_533%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_533%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_533:
   end subroutine reconstruct_d533
!---#] subroutine reconstruct_d533:
end module     p4_ubaru_hepemg_d533h3l1d
