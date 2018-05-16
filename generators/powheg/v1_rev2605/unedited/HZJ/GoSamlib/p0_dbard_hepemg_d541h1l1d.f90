module     p0_dbard_hepemg_d541h1l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbard_hepemg/helicity1d541h1l1d.f90
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
   public :: derivative , reconstruct_d541
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd541h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(46) :: acd541
      complex(ki) :: brack
      acd541(1)=dotproduct(k1,qshift)
      acd541(2)=dotproduct(e6,qshift)
      acd541(3)=abb541(23)
      acd541(4)=dotproduct(qshift,spvak2e6)
      acd541(5)=abb541(26)
      acd541(6)=abb541(12)
      acd541(7)=dotproduct(k2,qshift)
      acd541(8)=dotproduct(qshift,spvak2k4)
      acd541(9)=abb541(24)
      acd541(10)=dotproduct(qshift,spvak5k4)
      acd541(11)=abb541(28)
      acd541(12)=abb541(22)
      acd541(13)=dotproduct(k6,qshift)
      acd541(14)=abb541(20)
      acd541(15)=abb541(19)
      acd541(16)=dotproduct(qshift,qshift)
      acd541(17)=abb541(6)
      acd541(18)=abb541(4)
      acd541(19)=dotproduct(qshift,spvak2k1)
      acd541(20)=abb541(29)
      acd541(21)=abb541(15)
      acd541(22)=abb541(5)
      acd541(23)=dotproduct(qshift,spvak2k6)
      acd541(24)=abb541(10)
      acd541(25)=dotproduct(qshift,spvak5k1)
      acd541(26)=abb541(14)
      acd541(27)=abb541(13)
      acd541(28)=abb541(25)
      acd541(29)=abb541(21)
      acd541(30)=abb541(27)
      acd541(31)=abb541(7)
      acd541(32)=abb541(17)
      acd541(33)=abb541(11)
      acd541(34)=abb541(18)
      acd541(35)=abb541(8)
      acd541(36)=abb541(3)
      acd541(37)=abb541(9)
      acd541(38)=acd541(25)*acd541(26)
      acd541(39)=acd541(19)*acd541(22)
      acd541(40)=acd541(23)*acd541(24)
      acd541(41)=acd541(1)*acd541(3)
      acd541(42)=acd541(13)*acd541(14)
      acd541(43)=acd541(8)*acd541(18)
      acd541(44)=-acd541(16)*acd541(17)
      acd541(45)=acd541(19)*acd541(20)
      acd541(45)=acd541(21)+acd541(45)
      acd541(45)=acd541(10)*acd541(45)
      acd541(38)=acd541(45)+acd541(44)+acd541(43)+acd541(42)+acd541(41)+acd541(&
      &40)+acd541(39)-acd541(27)+acd541(38)
      acd541(38)=acd541(2)*acd541(38)
      acd541(39)=acd541(13)+acd541(7)
      acd541(39)=acd541(11)*acd541(39)
      acd541(40)=acd541(23)*acd541(34)
      acd541(41)=acd541(4)*acd541(32)
      acd541(42)=acd541(16)*acd541(30)
      acd541(39)=acd541(42)+acd541(41)-acd541(35)+acd541(40)+acd541(39)
      acd541(39)=acd541(10)*acd541(39)
      acd541(40)=-acd541(8)*acd541(29)
      acd541(41)=acd541(4)*acd541(28)
      acd541(40)=acd541(41)+acd541(31)+acd541(40)
      acd541(40)=acd541(16)*acd541(40)
      acd541(41)=-acd541(23)*acd541(36)
      acd541(42)=-acd541(7)*acd541(12)
      acd541(43)=acd541(1)*acd541(6)
      acd541(44)=-acd541(13)*acd541(15)
      acd541(45)=-acd541(7)+acd541(13)
      acd541(45)=acd541(8)*acd541(9)*acd541(45)
      acd541(46)=acd541(1)*acd541(5)
      acd541(46)=-acd541(33)+acd541(46)
      acd541(46)=acd541(4)*acd541(46)
      brack=acd541(37)+acd541(38)+acd541(39)+acd541(40)+acd541(41)+acd541(42)+a&
      &cd541(43)+acd541(44)+acd541(45)+acd541(46)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd541h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(59) :: acd541
      complex(ki) :: brack
      acd541(1)=k1(iv1)
      acd541(2)=dotproduct(e6,qshift)
      acd541(3)=abb541(23)
      acd541(4)=dotproduct(qshift,spvak2e6)
      acd541(5)=abb541(26)
      acd541(6)=abb541(12)
      acd541(7)=k2(iv1)
      acd541(8)=dotproduct(qshift,spvak2k4)
      acd541(9)=abb541(24)
      acd541(10)=dotproduct(qshift,spvak5k4)
      acd541(11)=abb541(28)
      acd541(12)=abb541(22)
      acd541(13)=k6(iv1)
      acd541(14)=abb541(20)
      acd541(15)=abb541(19)
      acd541(16)=e6(iv1)
      acd541(17)=dotproduct(k1,qshift)
      acd541(18)=dotproduct(k6,qshift)
      acd541(19)=dotproduct(qshift,qshift)
      acd541(20)=abb541(6)
      acd541(21)=abb541(4)
      acd541(22)=dotproduct(qshift,spvak2k1)
      acd541(23)=abb541(29)
      acd541(24)=abb541(15)
      acd541(25)=abb541(5)
      acd541(26)=dotproduct(qshift,spvak2k6)
      acd541(27)=abb541(10)
      acd541(28)=dotproduct(qshift,spvak5k1)
      acd541(29)=abb541(14)
      acd541(30)=abb541(13)
      acd541(31)=qshift(iv1)
      acd541(32)=abb541(25)
      acd541(33)=abb541(21)
      acd541(34)=abb541(27)
      acd541(35)=abb541(7)
      acd541(36)=spvak2e6(iv1)
      acd541(37)=abb541(17)
      acd541(38)=abb541(11)
      acd541(39)=spvak2k4(iv1)
      acd541(40)=dotproduct(k2,qshift)
      acd541(41)=spvak5k4(iv1)
      acd541(42)=abb541(18)
      acd541(43)=abb541(8)
      acd541(44)=spvak2k1(iv1)
      acd541(45)=spvak2k6(iv1)
      acd541(46)=abb541(3)
      acd541(47)=spvak5k1(iv1)
      acd541(48)=-acd541(29)*acd541(28)
      acd541(49)=-acd541(26)*acd541(27)
      acd541(50)=-acd541(3)*acd541(17)
      acd541(51)=-acd541(22)*acd541(25)
      acd541(52)=-acd541(18)*acd541(14)
      acd541(53)=acd541(19)*acd541(20)
      acd541(54)=-acd541(8)*acd541(21)
      acd541(55)=acd541(22)*acd541(23)
      acd541(55)=acd541(55)+acd541(24)
      acd541(56)=-acd541(10)*acd541(55)
      acd541(48)=acd541(56)+acd541(54)+acd541(53)+acd541(52)+acd541(51)+acd541(&
      &50)+acd541(49)+acd541(30)+acd541(48)
      acd541(48)=acd541(16)*acd541(48)
      acd541(49)=-acd541(10)*acd541(23)
      acd541(49)=acd541(49)-acd541(25)
      acd541(49)=acd541(44)*acd541(49)
      acd541(50)=-acd541(29)*acd541(47)
      acd541(51)=-acd541(45)*acd541(27)
      acd541(52)=-acd541(1)*acd541(3)
      acd541(53)=-acd541(39)*acd541(21)
      acd541(54)=-acd541(13)*acd541(14)
      acd541(56)=2.0_ki*acd541(31)
      acd541(57)=acd541(20)*acd541(56)
      acd541(55)=-acd541(41)*acd541(55)
      acd541(49)=acd541(55)+acd541(57)+acd541(54)+acd541(53)+acd541(52)+acd541(&
      &51)+acd541(50)+acd541(49)
      acd541(49)=acd541(2)*acd541(49)
      acd541(50)=-acd541(26)*acd541(42)
      acd541(51)=-acd541(4)*acd541(37)
      acd541(52)=-acd541(19)*acd541(34)
      acd541(53)=-acd541(40)-acd541(18)
      acd541(53)=acd541(11)*acd541(53)
      acd541(50)=acd541(53)+acd541(52)+acd541(51)+acd541(43)+acd541(50)
      acd541(50)=acd541(41)*acd541(50)
      acd541(51)=-acd541(45)*acd541(42)
      acd541(52)=-acd541(36)*acd541(37)
      acd541(53)=-acd541(7)-acd541(13)
      acd541(53)=acd541(11)*acd541(53)
      acd541(54)=-acd541(34)*acd541(56)
      acd541(51)=acd541(54)+acd541(53)+acd541(51)+acd541(52)
      acd541(51)=acd541(10)*acd541(51)
      acd541(52)=acd541(40)-acd541(18)
      acd541(52)=acd541(39)*acd541(52)
      acd541(53)=acd541(7)-acd541(13)
      acd541(53)=acd541(8)*acd541(53)
      acd541(52)=acd541(53)+acd541(52)
      acd541(52)=acd541(9)*acd541(52)
      acd541(53)=acd541(39)*acd541(33)
      acd541(54)=-acd541(36)*acd541(32)
      acd541(53)=acd541(53)+acd541(54)
      acd541(53)=acd541(19)*acd541(53)
      acd541(54)=-acd541(4)*acd541(32)
      acd541(55)=acd541(8)*acd541(33)
      acd541(54)=acd541(55)-acd541(35)+acd541(54)
      acd541(54)=acd541(54)*acd541(56)
      acd541(55)=acd541(45)*acd541(46)
      acd541(56)=acd541(7)*acd541(12)
      acd541(57)=-acd541(4)*acd541(5)
      acd541(57)=-acd541(6)+acd541(57)
      acd541(57)=acd541(1)*acd541(57)
      acd541(58)=-acd541(5)*acd541(17)
      acd541(58)=acd541(38)+acd541(58)
      acd541(58)=acd541(36)*acd541(58)
      acd541(59)=acd541(13)*acd541(15)
      brack=acd541(48)+acd541(49)+acd541(50)+acd541(51)+acd541(52)+acd541(53)+a&
      &cd541(54)+acd541(55)+acd541(56)+acd541(57)+acd541(58)+acd541(59)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd541h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(57) :: acd541
      complex(ki) :: brack
      acd541(1)=d(iv1,iv2)
      acd541(2)=dotproduct(e6,qshift)
      acd541(3)=abb541(6)
      acd541(4)=dotproduct(qshift,spvak2k4)
      acd541(5)=abb541(21)
      acd541(6)=dotproduct(qshift,spvak5k4)
      acd541(7)=abb541(27)
      acd541(8)=dotproduct(qshift,spvak2e6)
      acd541(9)=abb541(25)
      acd541(10)=abb541(7)
      acd541(11)=k1(iv1)
      acd541(12)=e6(iv2)
      acd541(13)=abb541(23)
      acd541(14)=spvak2e6(iv2)
      acd541(15)=abb541(26)
      acd541(16)=k1(iv2)
      acd541(17)=e6(iv1)
      acd541(18)=spvak2e6(iv1)
      acd541(19)=k2(iv1)
      acd541(20)=spvak2k4(iv2)
      acd541(21)=abb541(24)
      acd541(22)=spvak5k4(iv2)
      acd541(23)=abb541(28)
      acd541(24)=k2(iv2)
      acd541(25)=spvak2k4(iv1)
      acd541(26)=spvak5k4(iv1)
      acd541(27)=k6(iv1)
      acd541(28)=abb541(20)
      acd541(29)=k6(iv2)
      acd541(30)=qshift(iv2)
      acd541(31)=abb541(4)
      acd541(32)=dotproduct(qshift,spvak2k1)
      acd541(33)=abb541(29)
      acd541(34)=abb541(15)
      acd541(35)=spvak2k1(iv2)
      acd541(36)=abb541(5)
      acd541(37)=spvak2k6(iv2)
      acd541(38)=abb541(10)
      acd541(39)=spvak5k1(iv2)
      acd541(40)=abb541(14)
      acd541(41)=qshift(iv1)
      acd541(42)=spvak2k1(iv1)
      acd541(43)=spvak2k6(iv1)
      acd541(44)=spvak5k1(iv1)
      acd541(45)=abb541(17)
      acd541(46)=abb541(18)
      acd541(47)=acd541(33)*acd541(6)
      acd541(47)=acd541(47)+acd541(36)
      acd541(48)=acd541(35)*acd541(47)
      acd541(49)=acd541(40)*acd541(39)
      acd541(50)=acd541(37)*acd541(38)
      acd541(51)=acd541(13)*acd541(16)
      acd541(52)=acd541(29)*acd541(28)
      acd541(53)=2.0_ki*acd541(30)
      acd541(54)=-acd541(3)*acd541(53)
      acd541(55)=acd541(20)*acd541(31)
      acd541(56)=acd541(33)*acd541(32)
      acd541(56)=acd541(56)+acd541(34)
      acd541(57)=acd541(22)*acd541(56)
      acd541(48)=acd541(57)+acd541(55)+acd541(54)+acd541(52)+acd541(51)+acd541(&
      &49)+acd541(50)+acd541(48)
      acd541(48)=acd541(17)*acd541(48)
      acd541(47)=acd541(42)*acd541(47)
      acd541(49)=acd541(40)*acd541(44)
      acd541(50)=acd541(38)*acd541(43)
      acd541(51)=acd541(11)*acd541(13)
      acd541(52)=acd541(27)*acd541(28)
      acd541(54)=2.0_ki*acd541(41)
      acd541(55)=-acd541(3)*acd541(54)
      acd541(57)=acd541(25)*acd541(31)
      acd541(56)=acd541(26)*acd541(56)
      acd541(47)=acd541(56)+acd541(57)+acd541(55)+acd541(52)+acd541(51)+acd541(&
      &49)+acd541(50)+acd541(47)
      acd541(47)=acd541(12)*acd541(47)
      acd541(49)=acd541(37)*acd541(46)
      acd541(50)=acd541(14)*acd541(45)
      acd541(51)=acd541(7)*acd541(53)
      acd541(52)=acd541(24)+acd541(29)
      acd541(52)=acd541(23)*acd541(52)
      acd541(55)=acd541(33)*acd541(2)
      acd541(56)=acd541(35)*acd541(55)
      acd541(49)=acd541(56)+acd541(52)+acd541(51)+acd541(49)+acd541(50)
      acd541(49)=acd541(26)*acd541(49)
      acd541(50)=acd541(43)*acd541(46)
      acd541(51)=acd541(18)*acd541(45)
      acd541(52)=acd541(7)*acd541(54)
      acd541(56)=acd541(19)+acd541(27)
      acd541(56)=acd541(23)*acd541(56)
      acd541(55)=acd541(42)*acd541(55)
      acd541(50)=acd541(55)+acd541(56)+acd541(52)+acd541(50)+acd541(51)
      acd541(50)=acd541(22)*acd541(50)
      acd541(51)=acd541(9)*acd541(8)
      acd541(52)=acd541(6)*acd541(7)
      acd541(55)=-acd541(5)*acd541(4)
      acd541(56)=-acd541(2)*acd541(3)
      acd541(51)=acd541(56)+acd541(55)+acd541(52)+acd541(10)+acd541(51)
      acd541(51)=acd541(1)*acd541(51)
      acd541(52)=acd541(11)*acd541(15)
      acd541(55)=acd541(9)*acd541(54)
      acd541(52)=acd541(55)+acd541(52)
      acd541(52)=acd541(14)*acd541(52)
      acd541(55)=acd541(9)*acd541(18)
      acd541(56)=-acd541(25)*acd541(5)
      acd541(55)=acd541(56)+acd541(55)
      acd541(53)=acd541(53)*acd541(55)
      acd541(54)=-acd541(5)*acd541(54)
      acd541(55)=-acd541(19)+acd541(27)
      acd541(55)=acd541(21)*acd541(55)
      acd541(54)=acd541(54)+acd541(55)
      acd541(54)=acd541(20)*acd541(54)
      acd541(55)=acd541(18)*acd541(15)*acd541(16)
      acd541(56)=-acd541(24)+acd541(29)
      acd541(56)=acd541(21)*acd541(25)*acd541(56)
      brack=acd541(47)+acd541(48)+acd541(49)+acd541(50)+2.0_ki*acd541(51)+acd54&
      &1(52)+acd541(53)+acd541(54)+acd541(55)+acd541(56)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd541h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(29) :: acd541
      complex(ki) :: brack
      acd541(1)=d(iv1,iv2)
      acd541(2)=e6(iv3)
      acd541(3)=abb541(6)
      acd541(4)=spvak2k4(iv3)
      acd541(5)=abb541(21)
      acd541(6)=spvak5k4(iv3)
      acd541(7)=abb541(27)
      acd541(8)=spvak2e6(iv3)
      acd541(9)=abb541(25)
      acd541(10)=d(iv1,iv3)
      acd541(11)=e6(iv2)
      acd541(12)=spvak2k4(iv2)
      acd541(13)=spvak5k4(iv2)
      acd541(14)=spvak2e6(iv2)
      acd541(15)=d(iv2,iv3)
      acd541(16)=e6(iv1)
      acd541(17)=spvak2k4(iv1)
      acd541(18)=spvak5k4(iv1)
      acd541(19)=spvak2e6(iv1)
      acd541(20)=spvak2k1(iv3)
      acd541(21)=abb541(29)
      acd541(22)=spvak2k1(iv2)
      acd541(23)=spvak2k1(iv1)
      acd541(24)=-acd541(9)*acd541(19)
      acd541(25)=-acd541(7)*acd541(18)
      acd541(26)=acd541(5)*acd541(17)
      acd541(27)=acd541(3)*acd541(16)
      acd541(24)=acd541(27)+acd541(26)+acd541(24)+acd541(25)
      acd541(24)=acd541(15)*acd541(24)
      acd541(25)=-acd541(9)*acd541(14)
      acd541(26)=-acd541(7)*acd541(13)
      acd541(27)=acd541(5)*acd541(12)
      acd541(28)=acd541(3)*acd541(11)
      acd541(25)=acd541(28)+acd541(27)+acd541(25)+acd541(26)
      acd541(25)=acd541(10)*acd541(25)
      acd541(26)=-acd541(9)*acd541(8)
      acd541(27)=-acd541(6)*acd541(7)
      acd541(28)=acd541(5)*acd541(4)
      acd541(29)=acd541(2)*acd541(3)
      acd541(26)=acd541(29)+acd541(28)+acd541(26)+acd541(27)
      acd541(26)=acd541(1)*acd541(26)
      acd541(24)=acd541(26)+acd541(24)+acd541(25)
      acd541(25)=-acd541(13)*acd541(16)
      acd541(26)=-acd541(11)*acd541(18)
      acd541(25)=acd541(25)+acd541(26)
      acd541(25)=acd541(20)*acd541(25)
      acd541(26)=-acd541(16)*acd541(22)
      acd541(27)=-acd541(11)*acd541(23)
      acd541(26)=acd541(26)+acd541(27)
      acd541(26)=acd541(6)*acd541(26)
      acd541(27)=-acd541(18)*acd541(22)
      acd541(28)=-acd541(13)*acd541(23)
      acd541(27)=acd541(27)+acd541(28)
      acd541(27)=acd541(2)*acd541(27)
      acd541(25)=acd541(27)+acd541(26)+acd541(25)
      acd541(25)=acd541(21)*acd541(25)
      brack=2.0_ki*acd541(24)+acd541(25)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p0_dbard_hepemg_globalsl1, only: epspow
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_abbrevd541h1
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
      qshift = k6
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
!---#[ subroutine reconstruct_d541:
   subroutine     reconstruct_d541(coeffs)
      use p0_dbard_hepemg_groups, only: tensrec_info_group0
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group0), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_541:
      coeffs%coeffs_541%c0 = derivative(czip)
      coeffs%coeffs_541%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_541%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_541%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_541%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_541%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_541%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_541%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_541%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_541%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_541%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_541%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_541%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_541%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_541%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_541%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_541%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_541%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_541%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_541%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_541%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_541%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_541%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_541%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_541%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_541%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_541%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_541%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_541%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_541%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_541%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_541%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_541%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_541%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_541%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_541:
   end subroutine reconstruct_d541
!---#] subroutine reconstruct_d541:
end module     p0_dbard_hepemg_d541h1l1d
