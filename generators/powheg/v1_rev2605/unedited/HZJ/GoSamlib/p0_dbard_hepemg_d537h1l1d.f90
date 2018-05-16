module     p0_dbard_hepemg_d537h1l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbard_hepemg/helicity1d537h1l1d.f90
   ! generator: buildfortran_d.py
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d537
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd537h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(46) :: acd537
      complex(ki) :: brack
      acd537(1)=dotproduct(k1,qshift)
      acd537(2)=dotproduct(e6,qshift)
      acd537(3)=abb537(5)
      acd537(4)=abb537(4)
      acd537(5)=dotproduct(k2,qshift)
      acd537(6)=dotproduct(qshift,spvak2k4)
      acd537(7)=abb537(24)
      acd537(8)=dotproduct(qshift,spvak5k1)
      acd537(9)=abb537(16)
      acd537(10)=dotproduct(qshift,spvak5k4)
      acd537(11)=abb537(25)
      acd537(12)=dotproduct(qshift,spvak2e6)
      acd537(13)=abb537(11)
      acd537(14)=abb537(18)
      acd537(15)=dotproduct(k6,qshift)
      acd537(16)=abb537(17)
      acd537(17)=abb537(22)
      acd537(18)=abb537(29)
      acd537(19)=dotproduct(qshift,spvak2k1)
      acd537(20)=abb537(14)
      acd537(21)=abb537(13)
      acd537(22)=dotproduct(qshift,spvak2k6)
      acd537(23)=abb537(19)
      acd537(24)=abb537(12)
      acd537(25)=dotproduct(qshift,qshift)
      acd537(26)=abb537(26)
      acd537(27)=abb537(6)
      acd537(28)=abb537(27)
      acd537(29)=abb537(21)
      acd537(30)=abb537(20)
      acd537(31)=abb537(28)
      acd537(32)=abb537(15)
      acd537(33)=dotproduct(qshift,spvak6k1)
      acd537(34)=abb537(7)
      acd537(35)=abb537(23)
      acd537(36)=abb537(10)
      acd537(37)=abb537(8)
      acd537(38)=abb537(9)
      acd537(39)=-acd537(33)*acd537(34)
      acd537(40)=acd537(22)*acd537(32)
      acd537(41)=-acd537(25)*acd537(28)
      acd537(42)=acd537(12)*acd537(31)
      acd537(43)=acd537(5)*acd537(11)
      acd537(39)=acd537(43)+acd537(42)+acd537(41)+acd537(40)-acd537(35)+acd537(&
      &39)
      acd537(39)=acd537(10)*acd537(39)
      acd537(40)=acd537(10)*acd537(20)
      acd537(40)=acd537(40)+acd537(21)
      acd537(40)=acd537(19)*acd537(40)
      acd537(41)=acd537(5)+acd537(1)
      acd537(41)=acd537(3)*acd537(41)
      acd537(42)=acd537(22)*acd537(23)
      acd537(43)=acd537(8)*acd537(18)
      acd537(44)=acd537(6)*acd537(17)
      acd537(40)=acd537(44)+acd537(43)+acd537(42)-acd537(24)+acd537(41)+acd537(&
      &40)
      acd537(40)=acd537(2)*acd537(40)
      acd537(41)=acd537(8)*acd537(9)
      acd537(42)=acd537(6)*acd537(7)
      acd537(41)=acd537(41)-acd537(42)
      acd537(42)=-acd537(16)-acd537(41)
      acd537(42)=acd537(15)*acd537(42)
      acd537(43)=-acd537(8)*acd537(27)
      acd537(44)=-acd537(6)*acd537(26)
      acd537(43)=acd537(44)+acd537(30)+acd537(43)
      acd537(43)=acd537(25)*acd537(43)
      acd537(44)=-acd537(15)*acd537(13)
      acd537(45)=-acd537(25)*acd537(29)
      acd537(44)=acd537(45)-acd537(36)+acd537(44)
      acd537(44)=acd537(12)*acd537(44)
      acd537(45)=acd537(12)*acd537(13)
      acd537(41)=acd537(45)-acd537(14)+acd537(41)
      acd537(41)=acd537(5)*acd537(41)
      acd537(45)=-acd537(1)*acd537(4)
      acd537(46)=-acd537(22)*acd537(37)
      brack=acd537(38)+acd537(39)+acd537(40)+acd537(41)+acd537(42)+acd537(43)+a&
      &cd537(44)+acd537(45)+acd537(46)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd537h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(62) :: acd537
      complex(ki) :: brack
      acd537(1)=k1(iv1)
      acd537(2)=dotproduct(e6,qshift)
      acd537(3)=abb537(5)
      acd537(4)=abb537(4)
      acd537(5)=k2(iv1)
      acd537(6)=dotproduct(qshift,spvak2k4)
      acd537(7)=abb537(24)
      acd537(8)=dotproduct(qshift,spvak5k1)
      acd537(9)=abb537(16)
      acd537(10)=dotproduct(qshift,spvak5k4)
      acd537(11)=abb537(25)
      acd537(12)=dotproduct(qshift,spvak2e6)
      acd537(13)=abb537(11)
      acd537(14)=abb537(18)
      acd537(15)=k6(iv1)
      acd537(16)=abb537(17)
      acd537(17)=e6(iv1)
      acd537(18)=dotproduct(k1,qshift)
      acd537(19)=dotproduct(k2,qshift)
      acd537(20)=abb537(22)
      acd537(21)=abb537(29)
      acd537(22)=dotproduct(qshift,spvak2k1)
      acd537(23)=abb537(14)
      acd537(24)=abb537(13)
      acd537(25)=dotproduct(qshift,spvak2k6)
      acd537(26)=abb537(19)
      acd537(27)=abb537(12)
      acd537(28)=qshift(iv1)
      acd537(29)=abb537(26)
      acd537(30)=abb537(6)
      acd537(31)=abb537(27)
      acd537(32)=abb537(21)
      acd537(33)=abb537(20)
      acd537(34)=spvak2k4(iv1)
      acd537(35)=dotproduct(k6,qshift)
      acd537(36)=dotproduct(qshift,qshift)
      acd537(37)=spvak5k1(iv1)
      acd537(38)=spvak5k4(iv1)
      acd537(39)=abb537(28)
      acd537(40)=abb537(15)
      acd537(41)=dotproduct(qshift,spvak6k1)
      acd537(42)=abb537(7)
      acd537(43)=abb537(23)
      acd537(44)=spvak2e6(iv1)
      acd537(45)=abb537(10)
      acd537(46)=spvak2k1(iv1)
      acd537(47)=spvak2k6(iv1)
      acd537(48)=abb537(8)
      acd537(49)=spvak6k1(iv1)
      acd537(50)=acd537(10)*acd537(23)
      acd537(50)=acd537(50)+acd537(24)
      acd537(51)=acd537(46)*acd537(50)
      acd537(52)=acd537(5)+acd537(1)
      acd537(52)=acd537(3)*acd537(52)
      acd537(53)=acd537(47)*acd537(26)
      acd537(54)=acd537(37)*acd537(21)
      acd537(55)=acd537(34)*acd537(20)
      acd537(56)=acd537(38)*acd537(22)*acd537(23)
      acd537(51)=acd537(56)+acd537(55)+acd537(54)+acd537(53)+acd537(52)+acd537(&
      &51)
      acd537(51)=acd537(2)*acd537(51)
      acd537(52)=-acd537(42)*acd537(41)
      acd537(53)=acd537(25)*acd537(40)
      acd537(54)=-acd537(36)*acd537(31)
      acd537(55)=acd537(12)*acd537(39)
      acd537(56)=acd537(19)*acd537(11)
      acd537(52)=acd537(56)+acd537(55)+acd537(54)+acd537(53)-acd537(43)+acd537(&
      &52)
      acd537(52)=acd537(38)*acd537(52)
      acd537(53)=-acd537(42)*acd537(49)
      acd537(54)=acd537(47)*acd537(40)
      acd537(55)=acd537(44)*acd537(39)
      acd537(56)=2.0_ki*acd537(28)
      acd537(57)=-acd537(31)*acd537(56)
      acd537(58)=acd537(5)*acd537(11)
      acd537(53)=acd537(58)+acd537(57)+acd537(55)+acd537(53)+acd537(54)
      acd537(53)=acd537(10)*acd537(53)
      acd537(50)=acd537(22)*acd537(50)
      acd537(54)=acd537(19)+acd537(18)
      acd537(54)=acd537(3)*acd537(54)
      acd537(55)=acd537(25)*acd537(26)
      acd537(57)=acd537(8)*acd537(21)
      acd537(58)=acd537(6)*acd537(20)
      acd537(50)=acd537(58)+acd537(57)-acd537(27)+acd537(55)+acd537(54)+acd537(&
      &50)
      acd537(50)=acd537(17)*acd537(50)
      acd537(54)=acd537(12)*acd537(13)
      acd537(55)=acd537(8)*acd537(9)
      acd537(57)=acd537(6)*acd537(7)
      acd537(54)=-acd537(57)+acd537(54)+acd537(55)
      acd537(55)=-acd537(16)-acd537(54)
      acd537(55)=acd537(15)*acd537(55)
      acd537(57)=-acd537(12)*acd537(32)
      acd537(58)=-acd537(8)*acd537(30)
      acd537(59)=-acd537(6)*acd537(29)
      acd537(57)=acd537(59)+acd537(58)+acd537(33)+acd537(57)
      acd537(56)=acd537(57)*acd537(56)
      acd537(54)=-acd537(14)+acd537(54)
      acd537(54)=acd537(5)*acd537(54)
      acd537(57)=acd537(9)*acd537(37)
      acd537(58)=acd537(7)*acd537(34)
      acd537(57)=acd537(57)-acd537(58)
      acd537(58)=-acd537(35)*acd537(57)
      acd537(59)=-acd537(37)*acd537(30)
      acd537(60)=-acd537(34)*acd537(29)
      acd537(59)=acd537(60)+acd537(59)
      acd537(59)=acd537(36)*acd537(59)
      acd537(60)=-acd537(36)*acd537(32)
      acd537(61)=-acd537(13)*acd537(35)
      acd537(60)=acd537(61)-acd537(45)+acd537(60)
      acd537(60)=acd537(44)*acd537(60)
      acd537(61)=acd537(44)*acd537(13)
      acd537(57)=acd537(61)+acd537(57)
      acd537(57)=acd537(19)*acd537(57)
      acd537(61)=-acd537(1)*acd537(4)
      acd537(62)=-acd537(47)*acd537(48)
      brack=acd537(50)+acd537(51)+acd537(52)+acd537(53)+acd537(54)+acd537(55)+a&
      &cd537(56)+acd537(57)+acd537(58)+acd537(59)+acd537(60)+acd537(61)+acd537(&
      &62)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd537h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(60) :: acd537
      complex(ki) :: brack
      acd537(1)=d(iv1,iv2)
      acd537(2)=dotproduct(qshift,spvak2k4)
      acd537(3)=abb537(26)
      acd537(4)=dotproduct(qshift,spvak5k1)
      acd537(5)=abb537(6)
      acd537(6)=dotproduct(qshift,spvak5k4)
      acd537(7)=abb537(27)
      acd537(8)=dotproduct(qshift,spvak2e6)
      acd537(9)=abb537(21)
      acd537(10)=abb537(20)
      acd537(11)=k1(iv1)
      acd537(12)=e6(iv2)
      acd537(13)=abb537(5)
      acd537(14)=k1(iv2)
      acd537(15)=e6(iv1)
      acd537(16)=k2(iv1)
      acd537(17)=spvak2k4(iv2)
      acd537(18)=abb537(24)
      acd537(19)=spvak5k1(iv2)
      acd537(20)=abb537(16)
      acd537(21)=spvak5k4(iv2)
      acd537(22)=abb537(25)
      acd537(23)=spvak2e6(iv2)
      acd537(24)=abb537(11)
      acd537(25)=k2(iv2)
      acd537(26)=spvak2k4(iv1)
      acd537(27)=spvak5k1(iv1)
      acd537(28)=spvak5k4(iv1)
      acd537(29)=spvak2e6(iv1)
      acd537(30)=k6(iv1)
      acd537(31)=k6(iv2)
      acd537(32)=abb537(22)
      acd537(33)=abb537(29)
      acd537(34)=dotproduct(qshift,spvak2k1)
      acd537(35)=abb537(14)
      acd537(36)=spvak2k1(iv2)
      acd537(37)=abb537(13)
      acd537(38)=spvak2k6(iv2)
      acd537(39)=abb537(19)
      acd537(40)=spvak2k1(iv1)
      acd537(41)=spvak2k6(iv1)
      acd537(42)=qshift(iv1)
      acd537(43)=qshift(iv2)
      acd537(44)=abb537(28)
      acd537(45)=dotproduct(e6,qshift)
      acd537(46)=abb537(15)
      acd537(47)=spvak6k1(iv2)
      acd537(48)=abb537(7)
      acd537(49)=spvak6k1(iv1)
      acd537(50)=-acd537(48)*acd537(47)
      acd537(51)=acd537(38)*acd537(46)
      acd537(52)=2.0_ki*acd537(7)
      acd537(52)=-acd537(43)*acd537(52)
      acd537(53)=acd537(23)*acd537(44)
      acd537(54)=acd537(25)*acd537(22)
      acd537(55)=acd537(35)*acd537(45)
      acd537(56)=acd537(36)*acd537(55)
      acd537(50)=acd537(56)+acd537(54)+acd537(53)+acd537(52)+acd537(50)+acd537(&
      &51)
      acd537(50)=acd537(28)*acd537(50)
      acd537(51)=-acd537(48)*acd537(49)
      acd537(52)=acd537(41)*acd537(46)
      acd537(53)=2.0_ki*acd537(42)
      acd537(54)=-acd537(7)*acd537(53)
      acd537(56)=acd537(29)*acd537(44)
      acd537(57)=acd537(16)*acd537(22)
      acd537(55)=acd537(40)*acd537(55)
      acd537(51)=acd537(55)+acd537(57)+acd537(56)+acd537(54)+acd537(51)+acd537(&
      &52)
      acd537(51)=acd537(21)*acd537(51)
      acd537(52)=acd537(35)*acd537(6)
      acd537(52)=acd537(52)+acd537(37)
      acd537(54)=acd537(36)*acd537(52)
      acd537(55)=acd537(25)+acd537(14)
      acd537(55)=acd537(13)*acd537(55)
      acd537(56)=acd537(38)*acd537(39)
      acd537(57)=acd537(19)*acd537(33)
      acd537(58)=acd537(17)*acd537(32)
      acd537(59)=acd537(35)*acd537(34)
      acd537(60)=acd537(21)*acd537(59)
      acd537(54)=acd537(60)+acd537(58)+acd537(57)+acd537(56)+acd537(55)+acd537(&
      &54)
      acd537(54)=acd537(15)*acd537(54)
      acd537(52)=acd537(40)*acd537(52)
      acd537(55)=acd537(16)+acd537(11)
      acd537(55)=acd537(13)*acd537(55)
      acd537(56)=acd537(39)*acd537(41)
      acd537(57)=acd537(27)*acd537(33)
      acd537(58)=acd537(26)*acd537(32)
      acd537(59)=acd537(28)*acd537(59)
      acd537(52)=acd537(59)+acd537(58)+acd537(57)+acd537(56)+acd537(55)+acd537(&
      &52)
      acd537(52)=acd537(12)*acd537(52)
      acd537(55)=-acd537(9)*acd537(8)
      acd537(56)=-acd537(6)*acd537(7)
      acd537(57)=-acd537(5)*acd537(4)
      acd537(58)=-acd537(3)*acd537(2)
      acd537(55)=acd537(58)+acd537(57)+acd537(56)+acd537(10)+acd537(55)
      acd537(55)=acd537(1)*acd537(55)
      acd537(56)=-acd537(29)*acd537(9)
      acd537(57)=-acd537(27)*acd537(5)
      acd537(58)=-acd537(26)*acd537(3)
      acd537(56)=acd537(58)+acd537(56)+acd537(57)
      acd537(56)=acd537(43)*acd537(56)
      acd537(55)=acd537(56)+acd537(55)
      acd537(56)=acd537(24)*acd537(29)
      acd537(57)=acd537(20)*acd537(27)
      acd537(58)=acd537(18)*acd537(26)
      acd537(56)=-acd537(58)+acd537(56)+acd537(57)
      acd537(57)=acd537(25)-acd537(31)
      acd537(56)=acd537(56)*acd537(57)
      acd537(57)=acd537(23)*acd537(24)
      acd537(58)=acd537(19)*acd537(20)
      acd537(59)=-acd537(17)*acd537(18)
      acd537(57)=acd537(59)+acd537(57)+acd537(58)
      acd537(57)=acd537(16)*acd537(57)
      acd537(58)=-acd537(9)*acd537(53)
      acd537(59)=-acd537(24)*acd537(30)
      acd537(58)=acd537(58)+acd537(59)
      acd537(58)=acd537(23)*acd537(58)
      acd537(59)=-acd537(5)*acd537(53)
      acd537(60)=-acd537(20)*acd537(30)
      acd537(59)=acd537(59)+acd537(60)
      acd537(59)=acd537(19)*acd537(59)
      acd537(53)=-acd537(3)*acd537(53)
      acd537(60)=acd537(18)*acd537(30)
      acd537(53)=acd537(53)+acd537(60)
      acd537(53)=acd537(17)*acd537(53)
      brack=acd537(50)+acd537(51)+acd537(52)+acd537(53)+acd537(54)+2.0_ki*acd53&
      &7(55)+acd537(56)+acd537(57)+acd537(58)+acd537(59)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd537h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(32) :: acd537
      complex(ki) :: brack
      acd537(1)=d(iv1,iv2)
      acd537(2)=spvak2k4(iv3)
      acd537(3)=abb537(26)
      acd537(4)=spvak5k1(iv3)
      acd537(5)=abb537(6)
      acd537(6)=spvak5k4(iv3)
      acd537(7)=abb537(27)
      acd537(8)=spvak2e6(iv3)
      acd537(9)=abb537(21)
      acd537(10)=d(iv1,iv3)
      acd537(11)=spvak2k4(iv2)
      acd537(12)=spvak5k1(iv2)
      acd537(13)=spvak5k4(iv2)
      acd537(14)=spvak2e6(iv2)
      acd537(15)=d(iv2,iv3)
      acd537(16)=spvak2k4(iv1)
      acd537(17)=spvak5k1(iv1)
      acd537(18)=spvak5k4(iv1)
      acd537(19)=spvak2e6(iv1)
      acd537(20)=e6(iv1)
      acd537(21)=spvak2k1(iv3)
      acd537(22)=abb537(14)
      acd537(23)=spvak2k1(iv2)
      acd537(24)=e6(iv2)
      acd537(25)=spvak2k1(iv1)
      acd537(26)=e6(iv3)
      acd537(27)=-acd537(9)*acd537(19)
      acd537(28)=-acd537(7)*acd537(18)
      acd537(29)=-acd537(5)*acd537(17)
      acd537(30)=-acd537(3)*acd537(16)
      acd537(27)=acd537(30)+acd537(29)+acd537(27)+acd537(28)
      acd537(27)=acd537(15)*acd537(27)
      acd537(28)=-acd537(9)*acd537(14)
      acd537(29)=-acd537(7)*acd537(13)
      acd537(30)=-acd537(5)*acd537(12)
      acd537(31)=-acd537(3)*acd537(11)
      acd537(28)=acd537(31)+acd537(30)+acd537(28)+acd537(29)
      acd537(28)=acd537(10)*acd537(28)
      acd537(29)=-acd537(9)*acd537(8)
      acd537(30)=-acd537(6)*acd537(7)
      acd537(31)=-acd537(5)*acd537(4)
      acd537(32)=-acd537(3)*acd537(2)
      acd537(29)=acd537(32)+acd537(31)+acd537(29)+acd537(30)
      acd537(29)=acd537(1)*acd537(29)
      acd537(27)=acd537(29)+acd537(27)+acd537(28)
      acd537(28)=acd537(23)*acd537(26)
      acd537(29)=acd537(21)*acd537(24)
      acd537(28)=acd537(28)+acd537(29)
      acd537(28)=acd537(18)*acd537(28)
      acd537(29)=acd537(25)*acd537(26)
      acd537(30)=acd537(20)*acd537(21)
      acd537(29)=acd537(29)+acd537(30)
      acd537(29)=acd537(13)*acd537(29)
      acd537(30)=acd537(24)*acd537(25)
      acd537(31)=acd537(20)*acd537(23)
      acd537(30)=acd537(30)+acd537(31)
      acd537(30)=acd537(6)*acd537(30)
      acd537(28)=acd537(30)+acd537(28)+acd537(29)
      acd537(28)=acd537(22)*acd537(28)
      brack=2.0_ki*acd537(27)+acd537(28)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p0_dbard_hepemg_globalsl1, only: epspow
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_abbrevd537h1
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
      qshift = -k2
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
!---#[ subroutine reconstruct_d537:
   subroutine     reconstruct_d537(coeffs)
      use p0_dbard_hepemg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_537:
      coeffs%coeffs_537%c0 = derivative(czip)
      coeffs%coeffs_537%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_537%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_537%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_537%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_537%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_537%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_537%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_537%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_537%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_537%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_537%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_537%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_537%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_537%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_537%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_537%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_537%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_537%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_537%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_537%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_537%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_537%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_537%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_537%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_537%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_537%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_537%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_537%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_537%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_537%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_537%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_537%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_537%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_537%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_537:
   end subroutine reconstruct_d537
!---#] subroutine reconstruct_d537:
end module     p0_dbard_hepemg_d537h1l1d
