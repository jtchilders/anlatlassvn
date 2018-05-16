module     p12_sbars_hepemg_d45h1l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity1d45h1l1d.f90
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
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d45
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(90) :: acd45
      complex(ki) :: brack
      acd45(1)=dotproduct(k1,qshift)
      acd45(2)=dotproduct(e6,qshift)
      acd45(3)=abb45(35)
      acd45(4)=dotproduct(qshift,spvak2k1)
      acd45(5)=abb45(14)
      acd45(6)=dotproduct(qshift,spvak2k4)
      acd45(7)=abb45(77)
      acd45(8)=dotproduct(qshift,spvak5k1)
      acd45(9)=abb45(41)
      acd45(10)=dotproduct(qshift,spvae6k1)
      acd45(11)=abb45(67)
      acd45(12)=dotproduct(qshift,spvak2e6)
      acd45(13)=abb45(65)
      acd45(14)=dotproduct(qshift,spvae6k4)
      acd45(15)=abb45(50)
      acd45(16)=dotproduct(qshift,spvak5e6)
      acd45(17)=abb45(43)
      acd45(18)=abb45(30)
      acd45(19)=dotproduct(k2,qshift)
      acd45(20)=abb45(34)
      acd45(21)=dotproduct(l3,qshift)
      acd45(22)=abb45(51)
      acd45(23)=dotproduct(k6,qshift)
      acd45(24)=abb45(55)
      acd45(25)=abb45(37)
      acd45(26)=abb45(48)
      acd45(27)=dotproduct(qshift,qshift)
      acd45(28)=abb45(39)
      acd45(29)=abb45(9)
      acd45(30)=abb45(44)
      acd45(31)=abb45(28)
      acd45(32)=dotproduct(qshift,spvak2k6)
      acd45(33)=abb45(8)
      acd45(34)=dotproduct(qshift,spvak6k1)
      acd45(35)=abb45(70)
      acd45(36)=abb45(6)
      acd45(37)=abb45(11)
      acd45(38)=abb45(76)
      acd45(39)=abb45(13)
      acd45(40)=abb45(66)
      acd45(41)=abb45(64)
      acd45(42)=abb45(63)
      acd45(43)=abb45(59)
      acd45(44)=abb45(32)
      acd45(45)=abb45(16)
      acd45(46)=abb45(15)
      acd45(47)=dotproduct(qshift,spvak5k6)
      acd45(48)=abb45(20)
      acd45(49)=dotproduct(qshift,spvak6k4)
      acd45(50)=abb45(17)
      acd45(51)=abb45(7)
      acd45(52)=abb45(74)
      acd45(53)=abb45(24)
      acd45(54)=dotproduct(qshift,spvak5k4)
      acd45(55)=abb45(75)
      acd45(56)=abb45(12)
      acd45(57)=abb45(73)
      acd45(58)=abb45(52)
      acd45(59)=abb45(40)
      acd45(60)=abb45(25)
      acd45(61)=abb45(18)
      acd45(62)=abb45(26)
      acd45(63)=dotproduct(qshift,spvak1k4)
      acd45(64)=abb45(21)
      acd45(65)=abb45(72)
      acd45(66)=abb45(69)
      acd45(67)=abb45(19)
      acd45(68)=dotproduct(qshift,spvak5k2)
      acd45(69)=abb45(23)
      acd45(70)=abb45(10)
      acd45(71)=-acd45(28)*acd45(2)
      acd45(72)=-acd45(37)*acd45(4)
      acd45(73)=-acd45(38)*acd45(6)
      acd45(74)=-acd45(39)*acd45(8)
      acd45(75)=-acd45(40)*acd45(10)
      acd45(76)=acd45(41)*acd45(12)
      acd45(77)=-acd45(42)*acd45(14)
      acd45(78)=-acd45(43)*acd45(16)
      acd45(71)=acd45(44)+acd45(78)+acd45(77)+acd45(76)+acd45(75)+acd45(74)+acd&
      &45(73)+acd45(72)+acd45(71)
      acd45(71)=acd45(27)*acd45(71)
      acd45(72)=acd45(1)-acd45(19)
      acd45(72)=acd45(3)*acd45(72)
      acd45(73)=acd45(29)*acd45(4)
      acd45(74)=acd45(30)*acd45(6)
      acd45(75)=acd45(31)*acd45(8)
      acd45(76)=acd45(33)*acd45(32)
      acd45(77)=acd45(35)*acd45(34)
      acd45(72)=acd45(72)-acd45(36)+acd45(77)+acd45(76)+acd45(75)+acd45(74)+acd&
      &45(73)
      acd45(72)=acd45(2)*acd45(72)
      acd45(73)=-acd45(7)*acd45(6)
      acd45(74)=-acd45(9)*acd45(8)
      acd45(75)=acd45(15)*acd45(14)
      acd45(76)=acd45(17)*acd45(16)
      acd45(73)=acd45(76)+acd45(75)+acd45(74)+acd45(73)
      acd45(74)=acd45(19)+acd45(1)
      acd45(73)=acd45(74)*acd45(73)
      acd45(75)=acd45(45)*acd45(14)
      acd45(76)=acd45(46)*acd45(16)
      acd45(77)=acd45(48)*acd45(47)
      acd45(78)=acd45(50)*acd45(49)
      acd45(75)=-acd45(51)+acd45(78)+acd45(77)+acd45(76)+acd45(75)
      acd45(75)=acd45(4)*acd45(75)
      acd45(76)=acd45(55)*acd45(10)
      acd45(77)=acd45(57)*acd45(12)
      acd45(78)=acd45(61)*acd45(32)
      acd45(79)=acd45(62)*acd45(34)
      acd45(76)=-acd45(67)+acd45(79)+acd45(78)+acd45(77)+acd45(76)
      acd45(76)=acd45(54)*acd45(76)
      acd45(77)=acd45(5)*acd45(4)
      acd45(78)=acd45(11)*acd45(10)
      acd45(79)=-acd45(13)*acd45(12)
      acd45(77)=acd45(79)+acd45(77)+acd45(78)
      acd45(74)=acd45(74)-acd45(23)
      acd45(74)=acd45(74)*acd45(77)
      acd45(77)=acd45(24)*acd45(6)
      acd45(78)=acd45(25)*acd45(8)
      acd45(77)=-acd45(26)+acd45(78)+acd45(77)
      acd45(77)=acd45(23)*acd45(77)
      acd45(78)=-acd45(18)*acd45(1)
      acd45(79)=-acd45(20)*acd45(19)
      acd45(80)=-acd45(22)*acd45(21)
      acd45(81)=-acd45(52)*acd45(6)
      acd45(82)=-acd45(53)*acd45(8)
      acd45(83)=-acd45(56)*acd45(10)
      acd45(84)=-acd45(58)*acd45(12)
      acd45(85)=-acd45(59)*acd45(14)
      acd45(86)=-acd45(60)*acd45(16)
      acd45(87)=-acd45(64)*acd45(63)
      acd45(88)=-acd45(65)*acd45(47)
      acd45(89)=-acd45(66)*acd45(49)
      acd45(90)=-acd45(69)*acd45(68)
      brack=acd45(70)+acd45(71)+acd45(72)+acd45(73)+acd45(74)+acd45(75)+acd45(7&
      &6)+acd45(77)+acd45(78)+acd45(79)+acd45(80)+acd45(81)+acd45(82)+acd45(83)&
      &+acd45(84)+acd45(85)+acd45(86)+acd45(87)+acd45(88)+acd45(89)+acd45(90)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(111) :: acd45
      complex(ki) :: brack
      acd45(1)=k1(iv1)
      acd45(2)=dotproduct(e6,qshift)
      acd45(3)=abb45(35)
      acd45(4)=dotproduct(qshift,spvak2k1)
      acd45(5)=abb45(14)
      acd45(6)=dotproduct(qshift,spvak2k4)
      acd45(7)=abb45(77)
      acd45(8)=dotproduct(qshift,spvak5k1)
      acd45(9)=abb45(41)
      acd45(10)=dotproduct(qshift,spvae6k1)
      acd45(11)=abb45(67)
      acd45(12)=dotproduct(qshift,spvak2e6)
      acd45(13)=abb45(65)
      acd45(14)=dotproduct(qshift,spvae6k4)
      acd45(15)=abb45(50)
      acd45(16)=dotproduct(qshift,spvak5e6)
      acd45(17)=abb45(43)
      acd45(18)=abb45(30)
      acd45(19)=k2(iv1)
      acd45(20)=abb45(34)
      acd45(21)=l3(iv1)
      acd45(22)=abb45(51)
      acd45(23)=k6(iv1)
      acd45(24)=abb45(55)
      acd45(25)=abb45(37)
      acd45(26)=abb45(48)
      acd45(27)=e6(iv1)
      acd45(28)=dotproduct(k1,qshift)
      acd45(29)=dotproduct(k2,qshift)
      acd45(30)=dotproduct(qshift,qshift)
      acd45(31)=abb45(39)
      acd45(32)=abb45(9)
      acd45(33)=abb45(44)
      acd45(34)=abb45(28)
      acd45(35)=dotproduct(qshift,spvak2k6)
      acd45(36)=abb45(8)
      acd45(37)=dotproduct(qshift,spvak6k1)
      acd45(38)=abb45(70)
      acd45(39)=abb45(6)
      acd45(40)=qshift(iv1)
      acd45(41)=abb45(11)
      acd45(42)=abb45(76)
      acd45(43)=abb45(13)
      acd45(44)=abb45(66)
      acd45(45)=abb45(64)
      acd45(46)=abb45(63)
      acd45(47)=abb45(59)
      acd45(48)=abb45(32)
      acd45(49)=spvak2k1(iv1)
      acd45(50)=dotproduct(k6,qshift)
      acd45(51)=abb45(16)
      acd45(52)=abb45(15)
      acd45(53)=dotproduct(qshift,spvak5k6)
      acd45(54)=abb45(20)
      acd45(55)=dotproduct(qshift,spvak6k4)
      acd45(56)=abb45(17)
      acd45(57)=abb45(7)
      acd45(58)=spvak2k4(iv1)
      acd45(59)=abb45(74)
      acd45(60)=spvak5k1(iv1)
      acd45(61)=abb45(24)
      acd45(62)=spvae6k1(iv1)
      acd45(63)=dotproduct(qshift,spvak5k4)
      acd45(64)=abb45(75)
      acd45(65)=abb45(12)
      acd45(66)=spvak2e6(iv1)
      acd45(67)=abb45(73)
      acd45(68)=abb45(52)
      acd45(69)=spvae6k4(iv1)
      acd45(70)=abb45(40)
      acd45(71)=spvak5e6(iv1)
      acd45(72)=abb45(25)
      acd45(73)=spvak2k6(iv1)
      acd45(74)=abb45(18)
      acd45(75)=spvak6k1(iv1)
      acd45(76)=abb45(26)
      acd45(77)=spvak1k4(iv1)
      acd45(78)=abb45(21)
      acd45(79)=spvak5k6(iv1)
      acd45(80)=abb45(72)
      acd45(81)=spvak6k4(iv1)
      acd45(82)=abb45(69)
      acd45(83)=spvak5k4(iv1)
      acd45(84)=abb45(19)
      acd45(85)=spvak5k2(iv1)
      acd45(86)=abb45(23)
      acd45(87)=acd45(16)*acd45(47)
      acd45(88)=acd45(14)*acd45(46)
      acd45(89)=-acd45(12)*acd45(45)
      acd45(90)=acd45(10)*acd45(44)
      acd45(91)=acd45(8)*acd45(43)
      acd45(92)=acd45(6)*acd45(42)
      acd45(93)=acd45(2)*acd45(31)
      acd45(87)=acd45(93)+acd45(92)+acd45(91)+acd45(90)+acd45(89)+acd45(88)-acd&
      &45(48)+acd45(87)
      acd45(88)=2.0_ki*acd45(40)
      acd45(87)=acd45(87)*acd45(88)
      acd45(89)=-acd45(19)+acd45(23)
      acd45(89)=acd45(5)*acd45(89)
      acd45(90)=-acd45(56)*acd45(81)
      acd45(91)=-acd45(54)*acd45(79)
      acd45(92)=-acd45(71)*acd45(52)
      acd45(93)=-acd45(69)*acd45(51)
      acd45(88)=acd45(41)*acd45(88)
      acd45(94)=-acd45(27)*acd45(32)
      acd45(88)=acd45(94)+acd45(88)+acd45(93)+acd45(92)+acd45(90)+acd45(91)+acd&
      &45(89)
      acd45(88)=acd45(4)*acd45(88)
      acd45(89)=acd45(28)+acd45(29)
      acd45(90)=acd45(50)-acd45(89)
      acd45(90)=acd45(5)*acd45(90)
      acd45(91)=-acd45(56)*acd45(55)
      acd45(92)=-acd45(54)*acd45(53)
      acd45(93)=-acd45(16)*acd45(52)
      acd45(94)=-acd45(14)*acd45(51)
      acd45(95)=acd45(30)*acd45(41)
      acd45(96)=-acd45(2)*acd45(32)
      acd45(90)=acd45(96)+acd45(95)+acd45(94)+acd45(93)+acd45(92)+acd45(57)+acd&
      &45(91)+acd45(90)
      acd45(90)=acd45(49)*acd45(90)
      acd45(91)=acd45(71)*acd45(17)
      acd45(92)=acd45(69)*acd45(15)
      acd45(93)=acd45(60)*acd45(9)
      acd45(94)=acd45(58)*acd45(7)
      acd45(95)=acd45(13)*acd45(66)
      acd45(96)=acd45(11)*acd45(62)
      acd45(91)=-acd45(91)-acd45(92)+acd45(93)+acd45(94)+acd45(95)-acd45(96)
      acd45(89)=acd45(91)*acd45(89)
      acd45(91)=acd45(71)*acd45(47)
      acd45(92)=acd45(69)*acd45(46)
      acd45(93)=-acd45(66)*acd45(45)
      acd45(94)=acd45(62)*acd45(44)
      acd45(95)=acd45(60)*acd45(43)
      acd45(96)=acd45(58)*acd45(42)
      acd45(91)=acd45(96)+acd45(95)+acd45(94)+acd45(93)+acd45(91)+acd45(92)
      acd45(91)=acd45(30)*acd45(91)
      acd45(92)=-acd45(28)+acd45(29)
      acd45(92)=acd45(3)*acd45(92)
      acd45(93)=-acd45(37)*acd45(38)
      acd45(94)=-acd45(35)*acd45(36)
      acd45(95)=-acd45(8)*acd45(34)
      acd45(96)=-acd45(6)*acd45(33)
      acd45(97)=acd45(30)*acd45(31)
      acd45(92)=acd45(97)+acd45(96)+acd45(95)+acd45(94)+acd45(39)+acd45(93)+acd&
      &45(92)
      acd45(92)=acd45(27)*acd45(92)
      acd45(93)=acd45(16)*acd45(17)
      acd45(94)=acd45(14)*acd45(15)
      acd45(95)=acd45(8)*acd45(9)
      acd45(96)=acd45(6)*acd45(7)
      acd45(97)=acd45(13)*acd45(12)
      acd45(98)=acd45(11)*acd45(10)
      acd45(93)=-acd45(93)-acd45(94)+acd45(95)+acd45(96)+acd45(97)-acd45(98)
      acd45(94)=acd45(2)*acd45(3)
      acd45(95)=acd45(94)+acd45(20)+acd45(93)
      acd45(95)=acd45(19)*acd45(95)
      acd45(96)=-acd45(4)*acd45(5)
      acd45(93)=acd45(96)-acd45(94)+acd45(18)+acd45(93)
      acd45(93)=acd45(1)*acd45(93)
      acd45(94)=-acd45(37)*acd45(76)
      acd45(96)=-acd45(35)*acd45(74)
      acd45(97)=-acd45(12)*acd45(67)
      acd45(98)=-acd45(10)*acd45(64)
      acd45(94)=acd45(98)+acd45(97)+acd45(96)+acd45(84)+acd45(94)
      acd45(94)=acd45(83)*acd45(94)
      acd45(96)=-acd45(38)*acd45(75)
      acd45(97)=-acd45(36)*acd45(73)
      acd45(98)=-acd45(60)*acd45(34)
      acd45(99)=-acd45(58)*acd45(33)
      acd45(96)=acd45(99)+acd45(98)+acd45(96)+acd45(97)
      acd45(96)=acd45(2)*acd45(96)
      acd45(97)=-acd45(75)*acd45(76)
      acd45(98)=-acd45(73)*acd45(74)
      acd45(97)=acd45(97)+acd45(98)
      acd45(97)=acd45(63)*acd45(97)
      acd45(98)=-acd45(8)*acd45(25)
      acd45(99)=-acd45(6)*acd45(24)
      acd45(98)=acd45(99)+acd45(26)+acd45(98)
      acd45(98)=acd45(23)*acd45(98)
      acd45(99)=-acd45(66)*acd45(50)
      acd45(100)=-acd45(23)*acd45(12)
      acd45(99)=acd45(99)+acd45(100)
      acd45(99)=acd45(13)*acd45(99)
      acd45(100)=acd45(62)*acd45(50)
      acd45(101)=acd45(23)*acd45(10)
      acd45(100)=acd45(100)+acd45(101)
      acd45(100)=acd45(11)*acd45(100)
      acd45(101)=acd45(85)*acd45(86)
      acd45(102)=acd45(77)*acd45(78)
      acd45(103)=acd45(21)*acd45(22)
      acd45(104)=acd45(81)*acd45(82)
      acd45(105)=acd45(79)*acd45(80)
      acd45(106)=acd45(71)*acd45(72)
      acd45(107)=acd45(69)*acd45(70)
      acd45(108)=-acd45(63)*acd45(67)
      acd45(108)=acd45(68)+acd45(108)
      acd45(108)=acd45(66)*acd45(108)
      acd45(109)=-acd45(63)*acd45(64)
      acd45(109)=acd45(65)+acd45(109)
      acd45(109)=acd45(62)*acd45(109)
      acd45(110)=-acd45(50)*acd45(25)
      acd45(110)=acd45(61)+acd45(110)
      acd45(110)=acd45(60)*acd45(110)
      acd45(111)=-acd45(50)*acd45(24)
      acd45(111)=acd45(59)+acd45(111)
      acd45(111)=acd45(58)*acd45(111)
      brack=acd45(87)+acd45(88)+acd45(89)+acd45(90)+acd45(91)+acd45(92)+acd45(9&
      &3)+acd45(94)+acd45(95)+acd45(96)+acd45(97)+acd45(98)+acd45(99)+acd45(100&
      &)+acd45(101)+acd45(102)+acd45(103)+acd45(104)+acd45(105)+acd45(106)+acd4&
      &5(107)+acd45(108)+acd45(109)+acd45(110)+acd45(111)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(93) :: acd45
      complex(ki) :: brack
      acd45(1)=d(iv1,iv2)
      acd45(2)=dotproduct(e6,qshift)
      acd45(3)=abb45(39)
      acd45(4)=dotproduct(qshift,spvak2k1)
      acd45(5)=abb45(11)
      acd45(6)=dotproduct(qshift,spvak2k4)
      acd45(7)=abb45(76)
      acd45(8)=dotproduct(qshift,spvak5k1)
      acd45(9)=abb45(13)
      acd45(10)=dotproduct(qshift,spvae6k1)
      acd45(11)=abb45(66)
      acd45(12)=dotproduct(qshift,spvak2e6)
      acd45(13)=abb45(64)
      acd45(14)=dotproduct(qshift,spvae6k4)
      acd45(15)=abb45(63)
      acd45(16)=dotproduct(qshift,spvak5e6)
      acd45(17)=abb45(59)
      acd45(18)=abb45(32)
      acd45(19)=k1(iv1)
      acd45(20)=e6(iv2)
      acd45(21)=abb45(35)
      acd45(22)=spvak2k1(iv2)
      acd45(23)=abb45(14)
      acd45(24)=spvak2k4(iv2)
      acd45(25)=abb45(77)
      acd45(26)=spvak5k1(iv2)
      acd45(27)=abb45(41)
      acd45(28)=spvae6k1(iv2)
      acd45(29)=abb45(67)
      acd45(30)=spvak2e6(iv2)
      acd45(31)=abb45(65)
      acd45(32)=spvae6k4(iv2)
      acd45(33)=abb45(50)
      acd45(34)=spvak5e6(iv2)
      acd45(35)=abb45(43)
      acd45(36)=k1(iv2)
      acd45(37)=e6(iv1)
      acd45(38)=spvak2k1(iv1)
      acd45(39)=spvak2k4(iv1)
      acd45(40)=spvak5k1(iv1)
      acd45(41)=spvae6k1(iv1)
      acd45(42)=spvak2e6(iv1)
      acd45(43)=spvae6k4(iv1)
      acd45(44)=spvak5e6(iv1)
      acd45(45)=k2(iv1)
      acd45(46)=k2(iv2)
      acd45(47)=k6(iv1)
      acd45(48)=abb45(55)
      acd45(49)=abb45(37)
      acd45(50)=k6(iv2)
      acd45(51)=qshift(iv2)
      acd45(52)=abb45(9)
      acd45(53)=abb45(44)
      acd45(54)=abb45(28)
      acd45(55)=spvak2k6(iv2)
      acd45(56)=abb45(8)
      acd45(57)=spvak6k1(iv2)
      acd45(58)=abb45(70)
      acd45(59)=qshift(iv1)
      acd45(60)=spvak2k6(iv1)
      acd45(61)=spvak6k1(iv1)
      acd45(62)=abb45(16)
      acd45(63)=abb45(15)
      acd45(64)=spvak5k6(iv2)
      acd45(65)=abb45(20)
      acd45(66)=spvak6k4(iv2)
      acd45(67)=abb45(17)
      acd45(68)=spvak5k6(iv1)
      acd45(69)=spvak6k4(iv1)
      acd45(70)=spvak5k4(iv2)
      acd45(71)=abb45(75)
      acd45(72)=spvak5k4(iv1)
      acd45(73)=abb45(73)
      acd45(74)=abb45(18)
      acd45(75)=abb45(26)
      acd45(76)=-acd45(17)*acd45(16)
      acd45(77)=-acd45(15)*acd45(14)
      acd45(78)=acd45(13)*acd45(12)
      acd45(79)=-acd45(11)*acd45(10)
      acd45(80)=-acd45(9)*acd45(8)
      acd45(81)=-acd45(7)*acd45(6)
      acd45(82)=-acd45(5)*acd45(4)
      acd45(83)=-acd45(3)*acd45(2)
      acd45(76)=acd45(83)+acd45(82)+acd45(81)+acd45(80)+acd45(79)+acd45(78)+acd&
      &45(77)+acd45(18)+acd45(76)
      acd45(76)=acd45(1)*acd45(76)
      acd45(77)=acd45(36)+acd45(46)-acd45(50)
      acd45(77)=acd45(23)*acd45(77)
      acd45(78)=acd45(67)*acd45(66)
      acd45(79)=acd45(65)*acd45(64)
      acd45(80)=acd45(34)*acd45(63)
      acd45(81)=acd45(32)*acd45(62)
      acd45(82)=2.0_ki*acd45(51)
      acd45(83)=-acd45(5)*acd45(82)
      acd45(84)=acd45(20)*acd45(52)
      acd45(77)=acd45(84)+acd45(83)+acd45(81)+acd45(80)+acd45(78)+acd45(79)+acd&
      &45(77)
      acd45(77)=acd45(38)*acd45(77)
      acd45(78)=acd45(19)+acd45(45)-acd45(47)
      acd45(78)=acd45(23)*acd45(78)
      acd45(79)=acd45(67)*acd45(69)
      acd45(80)=acd45(65)*acd45(68)
      acd45(81)=acd45(44)*acd45(63)
      acd45(83)=acd45(43)*acd45(62)
      acd45(84)=2.0_ki*acd45(59)
      acd45(85)=-acd45(5)*acd45(84)
      acd45(86)=acd45(37)*acd45(52)
      acd45(78)=acd45(86)+acd45(85)+acd45(83)+acd45(81)+acd45(79)+acd45(80)+acd&
      &45(78)
      acd45(78)=acd45(22)*acd45(78)
      acd45(79)=-acd45(34)*acd45(17)
      acd45(80)=-acd45(32)*acd45(15)
      acd45(81)=acd45(30)*acd45(13)
      acd45(83)=-acd45(28)*acd45(11)
      acd45(85)=-acd45(26)*acd45(9)
      acd45(86)=-acd45(24)*acd45(7)
      acd45(79)=acd45(86)+acd45(85)+acd45(83)+acd45(81)+acd45(79)+acd45(80)
      acd45(79)=acd45(79)*acd45(84)
      acd45(80)=-acd45(44)*acd45(17)
      acd45(81)=-acd45(43)*acd45(15)
      acd45(83)=acd45(42)*acd45(13)
      acd45(85)=-acd45(41)*acd45(11)
      acd45(86)=-acd45(40)*acd45(9)
      acd45(87)=-acd45(39)*acd45(7)
      acd45(80)=acd45(87)+acd45(86)+acd45(85)+acd45(83)+acd45(80)+acd45(81)
      acd45(80)=acd45(80)*acd45(82)
      acd45(81)=acd45(35)*acd45(44)
      acd45(83)=acd45(33)*acd45(43)
      acd45(85)=acd45(40)*acd45(27)
      acd45(86)=acd45(39)*acd45(25)
      acd45(87)=acd45(31)*acd45(42)
      acd45(88)=acd45(29)*acd45(41)
      acd45(81)=-acd45(85)-acd45(86)-acd45(87)+acd45(88)+acd45(81)+acd45(83)
      acd45(83)=acd45(46)*acd45(81)
      acd45(85)=acd45(34)*acd45(35)
      acd45(86)=acd45(32)*acd45(33)
      acd45(87)=acd45(26)*acd45(27)
      acd45(88)=acd45(24)*acd45(25)
      acd45(89)=acd45(31)*acd45(30)
      acd45(90)=acd45(29)*acd45(28)
      acd45(85)=-acd45(87)-acd45(88)-acd45(89)+acd45(90)+acd45(85)+acd45(86)
      acd45(86)=acd45(45)*acd45(85)
      acd45(87)=acd45(57)*acd45(58)
      acd45(88)=acd45(55)*acd45(56)
      acd45(89)=acd45(26)*acd45(54)
      acd45(90)=acd45(24)*acd45(53)
      acd45(82)=-acd45(3)*acd45(82)
      acd45(91)=-acd45(46)*acd45(21)
      acd45(82)=acd45(91)+acd45(82)+acd45(90)+acd45(89)+acd45(87)+acd45(88)
      acd45(82)=acd45(37)*acd45(82)
      acd45(87)=acd45(37)*acd45(21)
      acd45(81)=acd45(87)+acd45(81)
      acd45(81)=acd45(36)*acd45(81)
      acd45(87)=acd45(58)*acd45(61)
      acd45(88)=acd45(56)*acd45(60)
      acd45(89)=acd45(40)*acd45(54)
      acd45(90)=acd45(39)*acd45(53)
      acd45(84)=-acd45(3)*acd45(84)
      acd45(91)=-acd45(45)*acd45(21)
      acd45(84)=acd45(91)+acd45(84)+acd45(90)+acd45(89)+acd45(87)+acd45(88)
      acd45(84)=acd45(20)*acd45(84)
      acd45(87)=acd45(20)*acd45(21)
      acd45(85)=acd45(87)+acd45(85)
      acd45(85)=acd45(19)*acd45(85)
      acd45(87)=acd45(57)*acd45(75)
      acd45(88)=acd45(55)*acd45(74)
      acd45(89)=acd45(30)*acd45(73)
      acd45(90)=acd45(28)*acd45(71)
      acd45(87)=acd45(90)+acd45(89)+acd45(87)+acd45(88)
      acd45(87)=acd45(72)*acd45(87)
      acd45(88)=acd45(61)*acd45(75)
      acd45(89)=acd45(60)*acd45(74)
      acd45(90)=acd45(42)*acd45(73)
      acd45(91)=acd45(41)*acd45(71)
      acd45(88)=acd45(91)+acd45(90)+acd45(88)+acd45(89)
      acd45(88)=acd45(70)*acd45(88)
      acd45(89)=acd45(40)*acd45(49)
      acd45(90)=acd45(39)*acd45(48)
      acd45(89)=acd45(90)+acd45(89)
      acd45(89)=acd45(50)*acd45(89)
      acd45(90)=acd45(26)*acd45(49)
      acd45(91)=acd45(24)*acd45(48)
      acd45(90)=acd45(91)+acd45(90)
      acd45(90)=acd45(47)*acd45(90)
      acd45(91)=acd45(42)*acd45(50)
      acd45(92)=acd45(30)*acd45(47)
      acd45(91)=acd45(91)+acd45(92)
      acd45(91)=acd45(31)*acd45(91)
      acd45(92)=-acd45(41)*acd45(50)
      acd45(93)=-acd45(28)*acd45(47)
      acd45(92)=acd45(92)+acd45(93)
      acd45(92)=acd45(29)*acd45(92)
      brack=2.0_ki*acd45(76)+acd45(77)+acd45(78)+acd45(79)+acd45(80)+acd45(81)+&
      &acd45(82)+acd45(83)+acd45(84)+acd45(85)+acd45(86)+acd45(87)+acd45(88)+ac&
      &d45(89)+acd45(90)+acd45(91)+acd45(92)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd45
      complex(ki) :: brack
      acd45(1)=d(iv1,iv2)
      acd45(2)=e6(iv3)
      acd45(3)=abb45(39)
      acd45(4)=spvak2k1(iv3)
      acd45(5)=abb45(11)
      acd45(6)=spvak2k4(iv3)
      acd45(7)=abb45(76)
      acd45(8)=spvak5k1(iv3)
      acd45(9)=abb45(13)
      acd45(10)=spvae6k1(iv3)
      acd45(11)=abb45(66)
      acd45(12)=spvak2e6(iv3)
      acd45(13)=abb45(64)
      acd45(14)=spvae6k4(iv3)
      acd45(15)=abb45(63)
      acd45(16)=spvak5e6(iv3)
      acd45(17)=abb45(59)
      acd45(18)=d(iv1,iv3)
      acd45(19)=e6(iv2)
      acd45(20)=spvak2k1(iv2)
      acd45(21)=spvak2k4(iv2)
      acd45(22)=spvak5k1(iv2)
      acd45(23)=spvae6k1(iv2)
      acd45(24)=spvak2e6(iv2)
      acd45(25)=spvae6k4(iv2)
      acd45(26)=spvak5e6(iv2)
      acd45(27)=d(iv2,iv3)
      acd45(28)=e6(iv1)
      acd45(29)=spvak2k1(iv1)
      acd45(30)=spvak2k4(iv1)
      acd45(31)=spvak5k1(iv1)
      acd45(32)=spvae6k1(iv1)
      acd45(33)=spvak2e6(iv1)
      acd45(34)=spvae6k4(iv1)
      acd45(35)=spvak5e6(iv1)
      acd45(36)=acd45(17)*acd45(35)
      acd45(37)=acd45(15)*acd45(34)
      acd45(38)=-acd45(13)*acd45(33)
      acd45(39)=acd45(11)*acd45(32)
      acd45(40)=acd45(9)*acd45(31)
      acd45(41)=acd45(7)*acd45(30)
      acd45(42)=acd45(5)*acd45(29)
      acd45(43)=acd45(3)*acd45(28)
      acd45(36)=acd45(43)+acd45(42)+acd45(41)+acd45(40)+acd45(39)+acd45(38)+acd&
      &45(36)+acd45(37)
      acd45(36)=acd45(27)*acd45(36)
      acd45(37)=acd45(17)*acd45(26)
      acd45(38)=acd45(15)*acd45(25)
      acd45(39)=-acd45(13)*acd45(24)
      acd45(40)=acd45(11)*acd45(23)
      acd45(41)=acd45(9)*acd45(22)
      acd45(42)=acd45(7)*acd45(21)
      acd45(43)=acd45(5)*acd45(20)
      acd45(44)=acd45(3)*acd45(19)
      acd45(37)=acd45(44)+acd45(43)+acd45(42)+acd45(41)+acd45(40)+acd45(39)+acd&
      &45(37)+acd45(38)
      acd45(37)=acd45(18)*acd45(37)
      acd45(38)=acd45(17)*acd45(16)
      acd45(39)=acd45(15)*acd45(14)
      acd45(40)=-acd45(13)*acd45(12)
      acd45(41)=acd45(11)*acd45(10)
      acd45(42)=acd45(9)*acd45(8)
      acd45(43)=acd45(7)*acd45(6)
      acd45(44)=acd45(5)*acd45(4)
      acd45(45)=acd45(3)*acd45(2)
      acd45(38)=acd45(45)+acd45(44)+acd45(43)+acd45(42)+acd45(41)+acd45(40)+acd&
      &45(38)+acd45(39)
      acd45(38)=acd45(1)*acd45(38)
      acd45(36)=acd45(38)+acd45(36)+acd45(37)
      brack=2.0_ki*acd45(36)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd45
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd45h1
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      integer, intent(in), optional :: i3
      integer, intent(in), optional :: i4
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = k3
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
      if(present(i4)) then
          iv4=i4
          deg=4
      else
          iv4=1
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
      if(deg.eq.4) then
         numerator = cond(epspow.eq.t1,brack_5,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
!---#[ subroutine reconstruct_d45:
   subroutine     reconstruct_d45(coeffs)
      use p12_sbars_hepemg_groups, only: tensrec_info_group5
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group5), intent(out) :: coeffs
      ! rank 4 case :
      !---[# reconstruct coeffs%coeffs_45:
      coeffs%coeffs_45%c0 = derivative(czip)
      coeffs%coeffs_45%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_45%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_45%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_45%c1(1,4) = derivative(czip,1,1,1,1)/ 24.0_ki
      coeffs%coeffs_45%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_45%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_45%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_45%c1(2,4) = derivative(czip,2,2,2,2)/ 24.0_ki
      coeffs%coeffs_45%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_45%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_45%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_45%c1(3,4) = derivative(czip,3,3,3,3)/ 24.0_ki
      coeffs%coeffs_45%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_45%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_45%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_45%c1(4,4) = derivative(czip,4,4,4,4)/ 24.0_ki
      coeffs%coeffs_45%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_45%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_45%c2(1,3) = -derivative(czip,1,2,2,2)/ 6.0_ki
      coeffs%coeffs_45%c2(1,4) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_45%c2(1,5) = derivative(czip,1,1,2,2)/ 4.0_ki
      coeffs%coeffs_45%c2(1,6) = -derivative(czip,1,1,1,2)/ 6.0_ki
      coeffs%coeffs_45%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_45%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_45%c2(2,3) = -derivative(czip,1,3,3,3)/ 6.0_ki
      coeffs%coeffs_45%c2(2,4) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_45%c2(2,5) = derivative(czip,1,1,3,3)/ 4.0_ki
      coeffs%coeffs_45%c2(2,6) = -derivative(czip,1,1,1,3)/ 6.0_ki
      coeffs%coeffs_45%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_45%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_45%c2(3,3) = -derivative(czip,1,4,4,4)/ 6.0_ki
      coeffs%coeffs_45%c2(3,4) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_45%c2(3,5) = derivative(czip,1,1,4,4)/ 4.0_ki
      coeffs%coeffs_45%c2(3,6) = -derivative(czip,1,1,1,4)/ 6.0_ki
      coeffs%coeffs_45%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_45%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_45%c2(4,3) = derivative(czip,2,3,3,3)/ 6.0_ki
      coeffs%coeffs_45%c2(4,4) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_45%c2(4,5) = derivative(czip,2,2,3,3)/ 4.0_ki
      coeffs%coeffs_45%c2(4,6) = derivative(czip,2,2,2,3)/ 6.0_ki
      coeffs%coeffs_45%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_45%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_45%c2(5,3) = derivative(czip,2,4,4,4)/ 6.0_ki
      coeffs%coeffs_45%c2(5,4) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_45%c2(5,5) = derivative(czip,2,2,4,4)/ 4.0_ki
      coeffs%coeffs_45%c2(5,6) = derivative(czip,2,2,2,4)/ 6.0_ki
      coeffs%coeffs_45%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_45%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_45%c2(6,3) = derivative(czip,3,4,4,4)/ 6.0_ki
      coeffs%coeffs_45%c2(6,4) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_45%c2(6,5) = derivative(czip,3,3,4,4)/ 4.0_ki
      coeffs%coeffs_45%c2(6,6) = derivative(czip,3,3,3,4)/ 6.0_ki
      coeffs%coeffs_45%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_45%c3(1,2) = -derivative(czip,1,2,3,3)/ 2.0_ki
      coeffs%coeffs_45%c3(1,3) = -derivative(czip,1,2,2,3)/ 2.0_ki
      coeffs%coeffs_45%c3(1,4) = derivative(czip,1,1,2,3)/ 2.0_ki
      coeffs%coeffs_45%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_45%c3(2,2) = -derivative(czip,1,2,4,4)/ 2.0_ki
      coeffs%coeffs_45%c3(2,3) = -derivative(czip,1,2,2,4)/ 2.0_ki
      coeffs%coeffs_45%c3(2,4) = derivative(czip,1,1,2,4)/ 2.0_ki
      coeffs%coeffs_45%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_45%c3(3,2) = -derivative(czip,1,3,4,4)/ 2.0_ki
      coeffs%coeffs_45%c3(3,3) = -derivative(czip,1,3,3,4)/ 2.0_ki
      coeffs%coeffs_45%c3(3,4) = derivative(czip,1,1,3,4)/ 2.0_ki
      coeffs%coeffs_45%c3(4,1) = -derivative(czip,2,3,4)
      coeffs%coeffs_45%c3(4,2) = derivative(czip,2,3,4,4)/ 2.0_ki
      coeffs%coeffs_45%c3(4,3) = derivative(czip,2,3,3,4)/ 2.0_ki
      coeffs%coeffs_45%c3(4,4) = derivative(czip,2,2,3,4)/ 2.0_ki
      coeffs%coeffs_45%c4(1,1) = -derivative(czip,1,2,3,4)
      !---#] reconstruct coeffs%coeffs_45:
   end subroutine reconstruct_d45
!---#] subroutine reconstruct_d45:
end module     p12_sbars_hepemg_d45h1l1d
