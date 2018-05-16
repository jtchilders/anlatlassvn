module     p12_sbars_hepemg_d277h1l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity1d277h1l1d.f90
   ! generator: buildfortran_d.py
   use p12_sbars_hepemg_config, only: ki
   use p12_sbars_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d277
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd277h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(82) :: acd277
      complex(ki) :: brack
      acd277(1)=dotproduct(k1,qshift)
      acd277(2)=dotproduct(e6,qshift)
      acd277(3)=abb277(50)
      acd277(4)=dotproduct(qshift,spvae6k1)
      acd277(5)=abb277(29)
      acd277(6)=dotproduct(qshift,spvak2e6)
      acd277(7)=abb277(34)
      acd277(8)=dotproduct(qshift,spvae6k4)
      acd277(9)=abb277(56)
      acd277(10)=dotproduct(qshift,spvak5e6)
      acd277(11)=abb277(52)
      acd277(12)=abb277(37)
      acd277(13)=dotproduct(k2,qshift)
      acd277(14)=abb277(49)
      acd277(15)=dotproduct(k6,qshift)
      acd277(16)=dotproduct(qshift,spvak2k1)
      acd277(17)=abb277(32)
      acd277(18)=dotproduct(qshift,spvak2k4)
      acd277(19)=abb277(68)
      acd277(20)=dotproduct(qshift,spvak5k1)
      acd277(21)=abb277(66)
      acd277(22)=abb277(30)
      acd277(23)=dotproduct(qshift,qshift)
      acd277(24)=abb277(41)
      acd277(25)=abb277(31)
      acd277(26)=abb277(67)
      acd277(27)=abb277(65)
      acd277(28)=dotproduct(qshift,spvak2k6)
      acd277(29)=abb277(48)
      acd277(30)=dotproduct(qshift,spvak6k1)
      acd277(31)=abb277(63)
      acd277(32)=abb277(20)
      acd277(33)=abb277(21)
      acd277(34)=abb277(27)
      acd277(35)=abb277(33)
      acd277(36)=abb277(17)
      acd277(37)=abb277(28)
      acd277(38)=abb277(40)
      acd277(39)=dotproduct(qshift,spvak5k4)
      acd277(40)=abb277(59)
      acd277(41)=abb277(62)
      acd277(42)=abb277(38)
      acd277(43)=abb277(18)
      acd277(44)=abb277(43)
      acd277(45)=abb277(24)
      acd277(46)=abb277(54)
      acd277(47)=abb277(35)
      acd277(48)=abb277(42)
      acd277(49)=dotproduct(qshift,spvak5k6)
      acd277(50)=abb277(51)
      acd277(51)=dotproduct(qshift,spvak6k4)
      acd277(52)=abb277(47)
      acd277(53)=dotproduct(qshift,spvak1e6)
      acd277(54)=abb277(19)
      acd277(55)=dotproduct(qshift,spvae6k2)
      acd277(56)=abb277(23)
      acd277(57)=abb277(26)
      acd277(58)=abb277(25)
      acd277(59)=abb277(22)
      acd277(60)=abb277(36)
      acd277(61)=abb277(39)
      acd277(62)=abb277(46)
      acd277(63)=abb277(44)
      acd277(64)=abb277(45)
      acd277(65)=abb277(16)
      acd277(66)=acd277(17)*acd277(15)
      acd277(67)=acd277(25)*acd277(2)
      acd277(68)=-acd277(35)*acd277(23)
      acd277(69)=acd277(38)*acd277(4)
      acd277(70)=acd277(42)*acd277(6)
      acd277(71)=acd277(45)*acd277(8)
      acd277(72)=acd277(47)*acd277(10)
      acd277(73)=acd277(50)*acd277(49)
      acd277(74)=acd277(52)*acd277(51)
      acd277(75)=acd277(54)*acd277(53)
      acd277(76)=acd277(56)*acd277(55)
      acd277(66)=-acd277(57)+acd277(76)+acd277(75)+acd277(74)+acd277(73)+acd277&
      &(72)+acd277(71)+acd277(70)+acd277(69)+acd277(68)+acd277(67)+acd277(66)
      acd277(66)=acd277(16)*acd277(66)
      acd277(67)=-acd277(24)*acd277(23)
      acd277(68)=acd277(26)*acd277(18)
      acd277(69)=acd277(27)*acd277(20)
      acd277(70)=acd277(29)*acd277(28)
      acd277(71)=acd277(31)*acd277(30)
      acd277(67)=-acd277(32)+acd277(71)+acd277(70)+acd277(69)+acd277(68)+acd277&
      &(67)
      acd277(67)=acd277(2)*acd277(67)
      acd277(68)=-acd277(33)*acd277(8)
      acd277(69)=-acd277(34)*acd277(10)
      acd277(70)=-acd277(36)*acd277(18)
      acd277(71)=-acd277(37)*acd277(20)
      acd277(68)=acd277(71)+acd277(70)+acd277(69)+acd277(68)
      acd277(68)=acd277(23)*acd277(68)
      acd277(69)=acd277(40)*acd277(4)
      acd277(70)=acd277(43)*acd277(6)
      acd277(71)=acd277(60)*acd277(28)
      acd277(72)=acd277(61)*acd277(30)
      acd277(69)=-acd277(64)+acd277(72)+acd277(71)+acd277(70)+acd277(69)
      acd277(69)=acd277(39)*acd277(69)
      acd277(70)=-acd277(3)*acd277(2)
      acd277(71)=acd277(9)*acd277(8)
      acd277(72)=-acd277(11)*acd277(10)
      acd277(70)=acd277(70)+acd277(72)+acd277(71)
      acd277(71)=acd277(13)+acd277(1)
      acd277(70)=acd277(71)*acd277(70)
      acd277(72)=acd277(5)*acd277(4)
      acd277(73)=-acd277(7)*acd277(6)
      acd277(72)=acd277(72)+acd277(73)
      acd277(71)=acd277(71)+acd277(15)-acd277(23)
      acd277(71)=acd277(71)*acd277(72)
      acd277(72)=acd277(19)*acd277(18)
      acd277(73)=acd277(21)*acd277(20)
      acd277(72)=-acd277(22)+acd277(73)+acd277(72)
      acd277(72)=acd277(15)*acd277(72)
      acd277(73)=-acd277(12)*acd277(1)
      acd277(74)=-acd277(14)*acd277(13)
      acd277(75)=-acd277(41)*acd277(4)
      acd277(76)=-acd277(44)*acd277(6)
      acd277(77)=-acd277(46)*acd277(8)
      acd277(78)=-acd277(48)*acd277(10)
      acd277(79)=-acd277(58)*acd277(18)
      acd277(80)=-acd277(59)*acd277(20)
      acd277(81)=-acd277(62)*acd277(49)
      acd277(82)=-acd277(63)*acd277(51)
      brack=acd277(65)+acd277(66)+acd277(67)+acd277(68)+acd277(69)+acd277(70)+a&
      &cd277(71)+acd277(72)+acd277(73)+acd277(74)+acd277(75)+acd277(76)+acd277(&
      &77)+acd277(78)+acd277(79)+acd277(80)+acd277(81)+acd277(82)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd277h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(106) :: acd277
      complex(ki) :: brack
      acd277(1)=k1(iv1)
      acd277(2)=dotproduct(e6,qshift)
      acd277(3)=abb277(50)
      acd277(4)=dotproduct(qshift,spvae6k1)
      acd277(5)=abb277(29)
      acd277(6)=dotproduct(qshift,spvak2e6)
      acd277(7)=abb277(34)
      acd277(8)=dotproduct(qshift,spvae6k4)
      acd277(9)=abb277(56)
      acd277(10)=dotproduct(qshift,spvak5e6)
      acd277(11)=abb277(52)
      acd277(12)=abb277(37)
      acd277(13)=k2(iv1)
      acd277(14)=abb277(49)
      acd277(15)=k6(iv1)
      acd277(16)=dotproduct(qshift,spvak2k1)
      acd277(17)=abb277(32)
      acd277(18)=dotproduct(qshift,spvak2k4)
      acd277(19)=abb277(68)
      acd277(20)=dotproduct(qshift,spvak5k1)
      acd277(21)=abb277(66)
      acd277(22)=abb277(30)
      acd277(23)=e6(iv1)
      acd277(24)=dotproduct(k1,qshift)
      acd277(25)=dotproduct(k2,qshift)
      acd277(26)=dotproduct(qshift,qshift)
      acd277(27)=abb277(41)
      acd277(28)=abb277(31)
      acd277(29)=abb277(67)
      acd277(30)=abb277(65)
      acd277(31)=dotproduct(qshift,spvak2k6)
      acd277(32)=abb277(48)
      acd277(33)=dotproduct(qshift,spvak6k1)
      acd277(34)=abb277(63)
      acd277(35)=abb277(20)
      acd277(36)=qshift(iv1)
      acd277(37)=abb277(21)
      acd277(38)=abb277(27)
      acd277(39)=abb277(33)
      acd277(40)=abb277(17)
      acd277(41)=abb277(28)
      acd277(42)=spvae6k1(iv1)
      acd277(43)=dotproduct(k6,qshift)
      acd277(44)=abb277(40)
      acd277(45)=dotproduct(qshift,spvak5k4)
      acd277(46)=abb277(59)
      acd277(47)=abb277(62)
      acd277(48)=spvak2e6(iv1)
      acd277(49)=abb277(38)
      acd277(50)=abb277(18)
      acd277(51)=abb277(43)
      acd277(52)=spvae6k4(iv1)
      acd277(53)=abb277(24)
      acd277(54)=abb277(54)
      acd277(55)=spvak5e6(iv1)
      acd277(56)=abb277(35)
      acd277(57)=abb277(42)
      acd277(58)=spvak2k1(iv1)
      acd277(59)=dotproduct(qshift,spvak5k6)
      acd277(60)=abb277(51)
      acd277(61)=dotproduct(qshift,spvak6k4)
      acd277(62)=abb277(47)
      acd277(63)=dotproduct(qshift,spvak1e6)
      acd277(64)=abb277(19)
      acd277(65)=dotproduct(qshift,spvae6k2)
      acd277(66)=abb277(23)
      acd277(67)=abb277(26)
      acd277(68)=spvak2k4(iv1)
      acd277(69)=abb277(25)
      acd277(70)=spvak5k1(iv1)
      acd277(71)=abb277(22)
      acd277(72)=spvak2k6(iv1)
      acd277(73)=abb277(36)
      acd277(74)=spvak6k1(iv1)
      acd277(75)=abb277(39)
      acd277(76)=spvak5k6(iv1)
      acd277(77)=abb277(46)
      acd277(78)=spvak6k4(iv1)
      acd277(79)=abb277(44)
      acd277(80)=spvak1e6(iv1)
      acd277(81)=spvae6k2(iv1)
      acd277(82)=spvak5k4(iv1)
      acd277(83)=abb277(45)
      acd277(84)=acd277(66)*acd277(81)
      acd277(85)=acd277(64)*acd277(80)
      acd277(86)=acd277(62)*acd277(78)
      acd277(87)=acd277(60)*acd277(76)
      acd277(88)=acd277(55)*acd277(56)
      acd277(89)=acd277(52)*acd277(53)
      acd277(90)=acd277(15)*acd277(17)
      acd277(91)=acd277(48)*acd277(49)
      acd277(92)=acd277(42)*acd277(44)
      acd277(93)=2.0_ki*acd277(36)
      acd277(94)=-acd277(39)*acd277(93)
      acd277(95)=acd277(23)*acd277(28)
      acd277(84)=acd277(95)+acd277(94)+acd277(92)+acd277(91)+acd277(90)+acd277(&
      &89)+acd277(88)+acd277(87)+acd277(86)+acd277(84)+acd277(85)
      acd277(84)=acd277(16)*acd277(84)
      acd277(85)=acd277(66)*acd277(65)
      acd277(86)=acd277(64)*acd277(63)
      acd277(87)=acd277(62)*acd277(61)
      acd277(88)=acd277(60)*acd277(59)
      acd277(89)=acd277(10)*acd277(56)
      acd277(90)=acd277(8)*acd277(53)
      acd277(91)=acd277(43)*acd277(17)
      acd277(92)=acd277(6)*acd277(49)
      acd277(94)=acd277(4)*acd277(44)
      acd277(95)=-acd277(26)*acd277(39)
      acd277(96)=acd277(2)*acd277(28)
      acd277(85)=acd277(96)+acd277(95)+acd277(94)+acd277(92)+acd277(91)+acd277(&
      &90)+acd277(89)+acd277(88)+acd277(87)+acd277(86)-acd277(67)+acd277(85)
      acd277(85)=acd277(58)*acd277(85)
      acd277(86)=acd277(1)+acd277(13)
      acd277(87)=-acd277(3)*acd277(86)
      acd277(88)=acd277(34)*acd277(74)
      acd277(89)=acd277(32)*acd277(72)
      acd277(90)=acd277(70)*acd277(30)
      acd277(91)=acd277(68)*acd277(29)
      acd277(92)=-acd277(27)*acd277(93)
      acd277(87)=acd277(92)+acd277(91)+acd277(90)+acd277(88)+acd277(89)+acd277(&
      &87)
      acd277(87)=acd277(2)*acd277(87)
      acd277(88)=acd277(24)+acd277(25)
      acd277(89)=-acd277(3)*acd277(88)
      acd277(90)=acd277(33)*acd277(34)
      acd277(91)=acd277(31)*acd277(32)
      acd277(92)=acd277(20)*acd277(30)
      acd277(94)=acd277(18)*acd277(29)
      acd277(95)=-acd277(26)*acd277(27)
      acd277(89)=acd277(95)+acd277(94)+acd277(92)+acd277(91)-acd277(35)+acd277(&
      &90)+acd277(89)
      acd277(89)=acd277(23)*acd277(89)
      acd277(90)=acd277(33)*acd277(75)
      acd277(91)=acd277(31)*acd277(73)
      acd277(92)=acd277(6)*acd277(50)
      acd277(94)=acd277(4)*acd277(46)
      acd277(90)=acd277(94)+acd277(92)+acd277(91)-acd277(83)+acd277(90)
      acd277(90)=acd277(82)*acd277(90)
      acd277(91)=-acd277(20)*acd277(41)
      acd277(92)=-acd277(18)*acd277(40)
      acd277(94)=-acd277(10)*acd277(38)
      acd277(95)=-acd277(8)*acd277(37)
      acd277(91)=acd277(95)+acd277(94)+acd277(91)+acd277(92)
      acd277(91)=acd277(91)*acd277(93)
      acd277(92)=-acd277(70)*acd277(41)
      acd277(94)=-acd277(68)*acd277(40)
      acd277(95)=-acd277(55)*acd277(38)
      acd277(96)=-acd277(52)*acd277(37)
      acd277(92)=acd277(96)+acd277(95)+acd277(92)+acd277(94)
      acd277(92)=acd277(26)*acd277(92)
      acd277(94)=acd277(55)*acd277(11)
      acd277(95)=acd277(52)*acd277(9)
      acd277(94)=acd277(94)-acd277(95)
      acd277(94)=-acd277(94)*acd277(88)
      acd277(95)=acd277(74)*acd277(75)
      acd277(96)=acd277(72)*acd277(73)
      acd277(95)=acd277(95)+acd277(96)
      acd277(95)=acd277(45)*acd277(95)
      acd277(96)=acd277(70)*acd277(21)
      acd277(97)=acd277(68)*acd277(19)
      acd277(96)=acd277(96)+acd277(97)
      acd277(96)=acd277(43)*acd277(96)
      acd277(97)=acd277(20)*acd277(21)
      acd277(98)=acd277(18)*acd277(19)
      acd277(97)=acd277(98)-acd277(22)+acd277(97)
      acd277(97)=acd277(15)*acd277(97)
      acd277(98)=acd277(10)*acd277(11)
      acd277(99)=acd277(8)*acd277(9)
      acd277(98)=acd277(98)-acd277(99)
      acd277(99)=-acd277(14)-acd277(98)
      acd277(99)=acd277(13)*acd277(99)
      acd277(98)=-acd277(12)-acd277(98)
      acd277(98)=acd277(1)*acd277(98)
      acd277(86)=-acd277(93)+acd277(15)+acd277(86)
      acd277(93)=-acd277(6)*acd277(86)
      acd277(88)=-acd277(26)+acd277(88)+acd277(43)
      acd277(100)=-acd277(48)*acd277(88)
      acd277(93)=acd277(93)+acd277(100)
      acd277(93)=acd277(7)*acd277(93)
      acd277(86)=acd277(4)*acd277(86)
      acd277(88)=acd277(42)*acd277(88)
      acd277(86)=acd277(86)+acd277(88)
      acd277(86)=acd277(5)*acd277(86)
      acd277(88)=-acd277(78)*acd277(79)
      acd277(100)=-acd277(76)*acd277(77)
      acd277(101)=-acd277(70)*acd277(71)
      acd277(102)=-acd277(68)*acd277(69)
      acd277(103)=-acd277(55)*acd277(57)
      acd277(104)=-acd277(52)*acd277(54)
      acd277(105)=acd277(45)*acd277(50)
      acd277(105)=-acd277(51)+acd277(105)
      acd277(105)=acd277(48)*acd277(105)
      acd277(106)=acd277(45)*acd277(46)
      acd277(106)=-acd277(47)+acd277(106)
      acd277(106)=acd277(42)*acd277(106)
      brack=acd277(84)+acd277(85)+acd277(86)+acd277(87)+acd277(88)+acd277(89)+a&
      &cd277(90)+acd277(91)+acd277(92)+acd277(93)+acd277(94)+acd277(95)+acd277(&
      &96)+acd277(97)+acd277(98)+acd277(99)+acd277(100)+acd277(101)+acd277(102)&
      &+acd277(103)+acd277(104)+acd277(105)+acd277(106)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd277h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(95) :: acd277
      complex(ki) :: brack
      acd277(1)=d(iv1,iv2)
      acd277(2)=dotproduct(e6,qshift)
      acd277(3)=abb277(41)
      acd277(4)=dotproduct(qshift,spvak2k1)
      acd277(5)=abb277(33)
      acd277(6)=dotproduct(qshift,spvak2k4)
      acd277(7)=abb277(17)
      acd277(8)=dotproduct(qshift,spvak5k1)
      acd277(9)=abb277(28)
      acd277(10)=dotproduct(qshift,spvae6k1)
      acd277(11)=abb277(29)
      acd277(12)=dotproduct(qshift,spvak2e6)
      acd277(13)=abb277(34)
      acd277(14)=dotproduct(qshift,spvae6k4)
      acd277(15)=abb277(21)
      acd277(16)=dotproduct(qshift,spvak5e6)
      acd277(17)=abb277(27)
      acd277(18)=k1(iv1)
      acd277(19)=e6(iv2)
      acd277(20)=abb277(50)
      acd277(21)=spvae6k1(iv2)
      acd277(22)=spvak2e6(iv2)
      acd277(23)=spvae6k4(iv2)
      acd277(24)=abb277(56)
      acd277(25)=spvak5e6(iv2)
      acd277(26)=abb277(52)
      acd277(27)=k1(iv2)
      acd277(28)=e6(iv1)
      acd277(29)=spvae6k1(iv1)
      acd277(30)=spvak2e6(iv1)
      acd277(31)=spvae6k4(iv1)
      acd277(32)=spvak5e6(iv1)
      acd277(33)=k2(iv1)
      acd277(34)=k2(iv2)
      acd277(35)=k6(iv1)
      acd277(36)=spvak2k1(iv2)
      acd277(37)=abb277(32)
      acd277(38)=spvak2k4(iv2)
      acd277(39)=abb277(68)
      acd277(40)=spvak5k1(iv2)
      acd277(41)=abb277(66)
      acd277(42)=k6(iv2)
      acd277(43)=spvak2k1(iv1)
      acd277(44)=spvak2k4(iv1)
      acd277(45)=spvak5k1(iv1)
      acd277(46)=qshift(iv2)
      acd277(47)=abb277(31)
      acd277(48)=abb277(67)
      acd277(49)=abb277(65)
      acd277(50)=spvak2k6(iv2)
      acd277(51)=abb277(48)
      acd277(52)=spvak6k1(iv2)
      acd277(53)=abb277(63)
      acd277(54)=qshift(iv1)
      acd277(55)=spvak2k6(iv1)
      acd277(56)=spvak6k1(iv1)
      acd277(57)=abb277(40)
      acd277(58)=abb277(38)
      acd277(59)=abb277(24)
      acd277(60)=abb277(35)
      acd277(61)=spvak5k6(iv2)
      acd277(62)=abb277(51)
      acd277(63)=spvak6k4(iv2)
      acd277(64)=abb277(47)
      acd277(65)=spvak1e6(iv2)
      acd277(66)=abb277(19)
      acd277(67)=spvae6k2(iv2)
      acd277(68)=abb277(23)
      acd277(69)=spvak5k6(iv1)
      acd277(70)=spvak6k4(iv1)
      acd277(71)=spvak1e6(iv1)
      acd277(72)=spvae6k2(iv1)
      acd277(73)=spvak5k4(iv2)
      acd277(74)=abb277(59)
      acd277(75)=spvak5k4(iv1)
      acd277(76)=abb277(18)
      acd277(77)=abb277(36)
      acd277(78)=abb277(39)
      acd277(79)=acd277(68)*acd277(67)
      acd277(80)=acd277(66)*acd277(65)
      acd277(81)=acd277(64)*acd277(63)
      acd277(82)=acd277(62)*acd277(61)
      acd277(83)=acd277(25)*acd277(60)
      acd277(84)=acd277(23)*acd277(59)
      acd277(85)=acd277(42)*acd277(37)
      acd277(86)=acd277(22)*acd277(58)
      acd277(87)=acd277(21)*acd277(57)
      acd277(88)=2.0_ki*acd277(46)
      acd277(89)=-acd277(5)*acd277(88)
      acd277(90)=acd277(19)*acd277(47)
      acd277(79)=acd277(90)+acd277(89)+acd277(87)+acd277(86)+acd277(85)+acd277(&
      &84)+acd277(83)+acd277(82)+acd277(81)+acd277(79)+acd277(80)
      acd277(79)=acd277(43)*acd277(79)
      acd277(80)=acd277(68)*acd277(72)
      acd277(81)=acd277(66)*acd277(71)
      acd277(82)=acd277(64)*acd277(70)
      acd277(83)=acd277(62)*acd277(69)
      acd277(84)=acd277(32)*acd277(60)
      acd277(85)=acd277(31)*acd277(59)
      acd277(86)=acd277(35)*acd277(37)
      acd277(87)=acd277(30)*acd277(58)
      acd277(89)=acd277(29)*acd277(57)
      acd277(90)=2.0_ki*acd277(54)
      acd277(91)=-acd277(5)*acd277(90)
      acd277(92)=acd277(28)*acd277(47)
      acd277(80)=acd277(92)+acd277(91)+acd277(89)+acd277(87)+acd277(86)+acd277(&
      &85)+acd277(84)+acd277(83)+acd277(82)+acd277(80)+acd277(81)
      acd277(80)=acd277(36)*acd277(80)
      acd277(81)=acd277(27)+acd277(34)
      acd277(82)=-acd277(20)*acd277(81)
      acd277(83)=acd277(52)*acd277(53)
      acd277(84)=acd277(50)*acd277(51)
      acd277(85)=acd277(40)*acd277(49)
      acd277(86)=acd277(38)*acd277(48)
      acd277(87)=-acd277(3)*acd277(88)
      acd277(82)=acd277(87)+acd277(86)+acd277(85)+acd277(83)+acd277(84)+acd277(&
      &82)
      acd277(82)=acd277(28)*acd277(82)
      acd277(83)=acd277(18)+acd277(33)
      acd277(84)=-acd277(20)*acd277(83)
      acd277(85)=acd277(53)*acd277(56)
      acd277(86)=acd277(51)*acd277(55)
      acd277(87)=acd277(45)*acd277(49)
      acd277(89)=acd277(44)*acd277(48)
      acd277(91)=-acd277(3)*acd277(90)
      acd277(84)=acd277(91)+acd277(89)+acd277(87)+acd277(85)+acd277(86)+acd277(&
      &84)
      acd277(84)=acd277(19)*acd277(84)
      acd277(85)=-acd277(17)*acd277(16)
      acd277(86)=-acd277(15)*acd277(14)
      acd277(87)=-acd277(9)*acd277(8)
      acd277(89)=-acd277(7)*acd277(6)
      acd277(91)=-acd277(5)*acd277(4)
      acd277(92)=-acd277(3)*acd277(2)
      acd277(85)=acd277(92)+acd277(91)+acd277(89)+acd277(87)+acd277(85)+acd277(&
      &86)
      acd277(86)=2.0_ki*acd277(1)
      acd277(85)=acd277(85)*acd277(86)
      acd277(87)=acd277(52)*acd277(78)
      acd277(89)=acd277(50)*acd277(77)
      acd277(91)=acd277(22)*acd277(76)
      acd277(92)=acd277(21)*acd277(74)
      acd277(87)=acd277(92)+acd277(91)+acd277(87)+acd277(89)
      acd277(87)=acd277(75)*acd277(87)
      acd277(89)=acd277(56)*acd277(78)
      acd277(91)=acd277(55)*acd277(77)
      acd277(92)=acd277(30)*acd277(76)
      acd277(93)=acd277(29)*acd277(74)
      acd277(89)=acd277(93)+acd277(92)+acd277(89)+acd277(91)
      acd277(89)=acd277(73)*acd277(89)
      acd277(91)=-acd277(9)*acd277(40)
      acd277(92)=-acd277(7)*acd277(38)
      acd277(93)=-acd277(25)*acd277(17)
      acd277(94)=-acd277(23)*acd277(15)
      acd277(91)=acd277(94)+acd277(93)+acd277(91)+acd277(92)
      acd277(91)=acd277(91)*acd277(90)
      acd277(92)=-acd277(9)*acd277(45)
      acd277(93)=-acd277(7)*acd277(44)
      acd277(94)=-acd277(32)*acd277(17)
      acd277(95)=-acd277(31)*acd277(15)
      acd277(92)=acd277(95)+acd277(94)+acd277(92)+acd277(93)
      acd277(92)=acd277(92)*acd277(88)
      acd277(88)=-acd277(88)+acd277(81)+acd277(42)
      acd277(93)=-acd277(30)*acd277(88)
      acd277(90)=-acd277(90)+acd277(83)+acd277(35)
      acd277(94)=-acd277(22)*acd277(90)
      acd277(95)=acd277(12)*acd277(86)
      acd277(93)=acd277(95)+acd277(94)+acd277(93)
      acd277(93)=acd277(13)*acd277(93)
      acd277(88)=acd277(29)*acd277(88)
      acd277(90)=acd277(21)*acd277(90)
      acd277(86)=-acd277(10)*acd277(86)
      acd277(86)=acd277(86)+acd277(90)+acd277(88)
      acd277(86)=acd277(11)*acd277(86)
      acd277(88)=acd277(26)*acd277(32)
      acd277(90)=acd277(24)*acd277(31)
      acd277(88)=acd277(88)-acd277(90)
      acd277(81)=-acd277(88)*acd277(81)
      acd277(88)=acd277(25)*acd277(26)
      acd277(90)=acd277(23)*acd277(24)
      acd277(88)=acd277(88)-acd277(90)
      acd277(83)=-acd277(88)*acd277(83)
      acd277(88)=acd277(45)*acd277(41)
      acd277(90)=acd277(44)*acd277(39)
      acd277(88)=acd277(88)+acd277(90)
      acd277(88)=acd277(42)*acd277(88)
      acd277(90)=acd277(40)*acd277(41)
      acd277(94)=acd277(38)*acd277(39)
      acd277(90)=acd277(90)+acd277(94)
      acd277(90)=acd277(35)*acd277(90)
      brack=acd277(79)+acd277(80)+acd277(81)+acd277(82)+acd277(83)+acd277(84)+a&
      &cd277(85)+acd277(86)+acd277(87)+acd277(88)+acd277(89)+acd277(90)+acd277(&
      &91)+acd277(92)+acd277(93)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd277h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd277
      complex(ki) :: brack
      acd277(1)=d(iv1,iv2)
      acd277(2)=e6(iv3)
      acd277(3)=abb277(41)
      acd277(4)=spvak2k1(iv3)
      acd277(5)=abb277(33)
      acd277(6)=spvak2k4(iv3)
      acd277(7)=abb277(17)
      acd277(8)=spvak5k1(iv3)
      acd277(9)=abb277(28)
      acd277(10)=spvae6k1(iv3)
      acd277(11)=abb277(29)
      acd277(12)=spvak2e6(iv3)
      acd277(13)=abb277(34)
      acd277(14)=spvae6k4(iv3)
      acd277(15)=abb277(21)
      acd277(16)=spvak5e6(iv3)
      acd277(17)=abb277(27)
      acd277(18)=d(iv1,iv3)
      acd277(19)=e6(iv2)
      acd277(20)=spvak2k1(iv2)
      acd277(21)=spvak2k4(iv2)
      acd277(22)=spvak5k1(iv2)
      acd277(23)=spvae6k1(iv2)
      acd277(24)=spvak2e6(iv2)
      acd277(25)=spvae6k4(iv2)
      acd277(26)=spvak5e6(iv2)
      acd277(27)=d(iv2,iv3)
      acd277(28)=e6(iv1)
      acd277(29)=spvak2k1(iv1)
      acd277(30)=spvak2k4(iv1)
      acd277(31)=spvak5k1(iv1)
      acd277(32)=spvae6k1(iv1)
      acd277(33)=spvak2e6(iv1)
      acd277(34)=spvae6k4(iv1)
      acd277(35)=spvak5e6(iv1)
      acd277(36)=-acd277(17)*acd277(35)
      acd277(37)=-acd277(15)*acd277(34)
      acd277(38)=acd277(13)*acd277(33)
      acd277(39)=-acd277(11)*acd277(32)
      acd277(40)=-acd277(9)*acd277(31)
      acd277(41)=-acd277(7)*acd277(30)
      acd277(42)=-acd277(5)*acd277(29)
      acd277(43)=-acd277(3)*acd277(28)
      acd277(36)=acd277(43)+acd277(42)+acd277(41)+acd277(40)+acd277(39)+acd277(&
      &38)+acd277(36)+acd277(37)
      acd277(36)=acd277(27)*acd277(36)
      acd277(37)=-acd277(17)*acd277(26)
      acd277(38)=-acd277(15)*acd277(25)
      acd277(39)=acd277(13)*acd277(24)
      acd277(40)=-acd277(11)*acd277(23)
      acd277(41)=-acd277(9)*acd277(22)
      acd277(42)=-acd277(7)*acd277(21)
      acd277(43)=-acd277(5)*acd277(20)
      acd277(44)=-acd277(3)*acd277(19)
      acd277(37)=acd277(44)+acd277(43)+acd277(42)+acd277(41)+acd277(40)+acd277(&
      &39)+acd277(37)+acd277(38)
      acd277(37)=acd277(18)*acd277(37)
      acd277(38)=-acd277(17)*acd277(16)
      acd277(39)=-acd277(15)*acd277(14)
      acd277(40)=acd277(13)*acd277(12)
      acd277(41)=-acd277(11)*acd277(10)
      acd277(42)=-acd277(9)*acd277(8)
      acd277(43)=-acd277(7)*acd277(6)
      acd277(44)=-acd277(5)*acd277(4)
      acd277(45)=-acd277(3)*acd277(2)
      acd277(38)=acd277(45)+acd277(44)+acd277(43)+acd277(42)+acd277(41)+acd277(&
      &40)+acd277(38)+acd277(39)
      acd277(38)=acd277(1)*acd277(38)
      acd277(36)=acd277(38)+acd277(36)+acd277(37)
      brack=2.0_ki*acd277(36)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd277h1
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
      qshift = k6+k5+k4
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
!---#[ subroutine reconstruct_d277:
   subroutine     reconstruct_d277(coeffs)
      use p12_sbars_hepemg_groups, only: tensrec_info_group5
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group5), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_277:
      coeffs%coeffs_277%c0 = derivative(czip)
      coeffs%coeffs_277%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_277%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_277%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_277%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_277%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_277%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_277%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_277%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_277%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_277%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_277%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_277%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_277%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_277%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_277%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_277%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_277%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_277%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_277%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_277%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_277%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_277%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_277%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_277%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_277%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_277%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_277%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_277%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_277%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_277%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_277%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_277%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_277%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_277%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_277:
   end subroutine reconstruct_d277
!---#] subroutine reconstruct_d277:
end module     p12_sbars_hepemg_d277h1l1d
