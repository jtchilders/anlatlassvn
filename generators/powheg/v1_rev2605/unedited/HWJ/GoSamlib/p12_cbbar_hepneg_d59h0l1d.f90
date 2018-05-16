module     p12_cbbar_hepneg_d59h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d59h0l1d.f90
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
   public :: derivative , reconstruct_d59
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd59h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(44) :: acd59
      complex(ki) :: brack
      acd59(1)=dotproduct(k1,qshift)
      acd59(2)=dotproduct(e6,qshift)
      acd59(3)=abb59(11)
      acd59(4)=dotproduct(qshift,spvak2k4)
      acd59(5)=abb59(14)
      acd59(6)=dotproduct(qshift,spvak5k1)
      acd59(7)=abb59(27)
      acd59(8)=dotproduct(qshift,spvak5k4)
      acd59(9)=abb59(13)
      acd59(10)=dotproduct(qshift,spvae6k1)
      acd59(11)=abb59(20)
      acd59(12)=abb59(15)
      acd59(13)=dotproduct(k2,qshift)
      acd59(14)=abb59(4)
      acd59(15)=dotproduct(k6,qshift)
      acd59(16)=abb59(8)
      acd59(17)=abb59(9)
      acd59(18)=dotproduct(qshift,spvak2k1)
      acd59(19)=abb59(12)
      acd59(20)=abb59(16)
      acd59(21)=dotproduct(qshift,spvak6k1)
      acd59(22)=abb59(21)
      acd59(23)=abb59(10)
      acd59(24)=dotproduct(qshift,qshift)
      acd59(25)=abb59(28)
      acd59(26)=abb59(5)
      acd59(27)=abb59(24)
      acd59(28)=abb59(19)
      acd59(29)=abb59(17)
      acd59(30)=abb59(25)
      acd59(31)=abb59(6)
      acd59(32)=abb59(26)
      acd59(33)=dotproduct(qshift,spvak2k6)
      acd59(34)=abb59(22)
      acd59(35)=abb59(23)
      acd59(36)=abb59(18)
      acd59(37)=abb59(7)
      acd59(38)=-acd59(13)-acd59(1)
      acd59(38)=acd59(38)*acd59(3)
      acd59(39)=acd59(16)*acd59(4)
      acd59(40)=acd59(17)*acd59(6)
      acd59(41)=acd59(18)*acd59(8)
      acd59(42)=-acd59(19)*acd59(41)
      acd59(43)=acd59(20)*acd59(18)
      acd59(44)=acd59(22)*acd59(21)
      acd59(38)=acd59(43)+acd59(38)-acd59(23)+acd59(44)+acd59(42)+acd59(40)+acd&
      &59(39)
      acd59(38)=acd59(2)*acd59(38)
      acd59(39)=acd59(9)*acd59(1)
      acd59(40)=-acd59(27)*acd59(24)
      acd59(42)=acd59(30)*acd59(10)
      acd59(43)=acd59(32)*acd59(21)
      acd59(44)=acd59(34)*acd59(33)
      acd59(39)=-acd59(35)+acd59(44)+acd59(43)+acd59(42)+acd59(40)+acd59(39)
      acd59(39)=acd59(8)*acd59(39)
      acd59(40)=acd59(5)*acd59(4)
      acd59(42)=-acd59(7)*acd59(6)
      acd59(43)=acd59(11)*acd59(10)
      acd59(40)=acd59(12)+acd59(43)+acd59(42)+acd59(40)
      acd59(42)=acd59(15)-acd59(1)
      acd59(40)=acd59(42)*acd59(40)
      acd59(42)=-acd59(25)*acd59(4)
      acd59(43)=-acd59(26)*acd59(6)
      acd59(44)=-acd59(28)*acd59(10)
      acd59(42)=acd59(29)+acd59(44)+acd59(43)+acd59(42)
      acd59(42)=acd59(24)*acd59(42)
      acd59(43)=-acd59(14)*acd59(13)
      acd59(41)=acd59(31)*acd59(41)
      acd59(44)=-acd59(36)*acd59(10)
      brack=acd59(37)+acd59(38)+acd59(39)+acd59(40)+acd59(41)+acd59(42)+acd59(4&
      &3)+acd59(44)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd59h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(59) :: acd59
      complex(ki) :: brack
      acd59(1)=k1(iv1)
      acd59(2)=dotproduct(e6,qshift)
      acd59(3)=abb59(11)
      acd59(4)=dotproduct(qshift,spvak2k4)
      acd59(5)=abb59(14)
      acd59(6)=dotproduct(qshift,spvak5k1)
      acd59(7)=abb59(27)
      acd59(8)=dotproduct(qshift,spvak5k4)
      acd59(9)=abb59(13)
      acd59(10)=dotproduct(qshift,spvae6k1)
      acd59(11)=abb59(20)
      acd59(12)=abb59(15)
      acd59(13)=k2(iv1)
      acd59(14)=abb59(4)
      acd59(15)=k6(iv1)
      acd59(16)=e6(iv1)
      acd59(17)=dotproduct(k1,qshift)
      acd59(18)=dotproduct(k2,qshift)
      acd59(19)=abb59(8)
      acd59(20)=abb59(9)
      acd59(21)=dotproduct(qshift,spvak2k1)
      acd59(22)=abb59(12)
      acd59(23)=abb59(16)
      acd59(24)=dotproduct(qshift,spvak6k1)
      acd59(25)=abb59(21)
      acd59(26)=abb59(10)
      acd59(27)=qshift(iv1)
      acd59(28)=abb59(28)
      acd59(29)=abb59(5)
      acd59(30)=abb59(24)
      acd59(31)=abb59(19)
      acd59(32)=abb59(17)
      acd59(33)=spvak2k4(iv1)
      acd59(34)=dotproduct(k6,qshift)
      acd59(35)=dotproduct(qshift,qshift)
      acd59(36)=spvak5k1(iv1)
      acd59(37)=spvak5k4(iv1)
      acd59(38)=abb59(25)
      acd59(39)=abb59(6)
      acd59(40)=abb59(26)
      acd59(41)=dotproduct(qshift,spvak2k6)
      acd59(42)=abb59(22)
      acd59(43)=abb59(23)
      acd59(44)=spvae6k1(iv1)
      acd59(45)=abb59(18)
      acd59(46)=spvak2k1(iv1)
      acd59(47)=spvak6k1(iv1)
      acd59(48)=spvak2k6(iv1)
      acd59(49)=-acd59(42)*acd59(48)
      acd59(50)=-acd59(40)*acd59(47)
      acd59(51)=-acd59(46)*acd59(39)
      acd59(52)=-acd59(44)*acd59(38)
      acd59(53)=2.0_ki*acd59(27)
      acd59(54)=acd59(30)*acd59(53)
      acd59(55)=-acd59(1)*acd59(9)
      acd59(56)=acd59(21)*acd59(22)
      acd59(57)=acd59(16)*acd59(56)
      acd59(49)=acd59(57)+acd59(55)+acd59(54)+acd59(52)+acd59(51)+acd59(49)+acd&
      &59(50)
      acd59(49)=acd59(8)*acd59(49)
      acd59(50)=-acd59(42)*acd59(41)
      acd59(51)=-acd59(24)*acd59(40)
      acd59(52)=acd59(35)*acd59(30)
      acd59(54)=-acd59(21)*acd59(39)
      acd59(55)=-acd59(10)*acd59(38)
      acd59(57)=-acd59(17)*acd59(9)
      acd59(50)=acd59(57)+acd59(55)+acd59(54)+acd59(52)+acd59(51)+acd59(43)+acd&
      &59(50)
      acd59(50)=acd59(37)*acd59(50)
      acd59(51)=acd59(8)*acd59(22)
      acd59(51)=acd59(51)-acd59(23)
      acd59(51)=acd59(46)*acd59(51)
      acd59(52)=acd59(1)+acd59(13)
      acd59(52)=acd59(3)*acd59(52)
      acd59(54)=-acd59(25)*acd59(47)
      acd59(55)=-acd59(36)*acd59(20)
      acd59(57)=-acd59(33)*acd59(19)
      acd59(56)=acd59(37)*acd59(56)
      acd59(51)=acd59(56)+acd59(57)+acd59(55)+acd59(54)+acd59(52)+acd59(51)
      acd59(51)=acd59(2)*acd59(51)
      acd59(52)=acd59(17)+acd59(18)
      acd59(52)=acd59(3)*acd59(52)
      acd59(54)=-acd59(24)*acd59(25)
      acd59(55)=-acd59(21)*acd59(23)
      acd59(56)=-acd59(6)*acd59(20)
      acd59(57)=-acd59(4)*acd59(19)
      acd59(52)=acd59(57)+acd59(56)+acd59(55)+acd59(26)+acd59(54)+acd59(52)
      acd59(52)=acd59(16)*acd59(52)
      acd59(54)=acd59(10)*acd59(11)
      acd59(55)=acd59(6)*acd59(7)
      acd59(56)=acd59(4)*acd59(5)
      acd59(54)=acd59(54)-acd59(55)+acd59(56)+acd59(12)
      acd59(55)=acd59(1)-acd59(15)
      acd59(54)=acd59(54)*acd59(55)
      acd59(55)=acd59(10)*acd59(31)
      acd59(56)=acd59(6)*acd59(29)
      acd59(57)=acd59(4)*acd59(28)
      acd59(55)=acd59(57)+acd59(56)-acd59(32)+acd59(55)
      acd59(53)=acd59(55)*acd59(53)
      acd59(55)=acd59(7)*acd59(36)
      acd59(56)=acd59(5)*acd59(33)
      acd59(55)=acd59(55)-acd59(56)
      acd59(56)=acd59(34)*acd59(55)
      acd59(57)=acd59(36)*acd59(29)
      acd59(58)=acd59(33)*acd59(28)
      acd59(57)=acd59(58)+acd59(57)
      acd59(57)=acd59(35)*acd59(57)
      acd59(58)=acd59(35)*acd59(31)
      acd59(59)=-acd59(11)*acd59(34)
      acd59(58)=acd59(59)+acd59(45)+acd59(58)
      acd59(58)=acd59(44)*acd59(58)
      acd59(59)=acd59(44)*acd59(11)
      acd59(55)=acd59(59)-acd59(55)
      acd59(55)=acd59(17)*acd59(55)
      acd59(59)=acd59(13)*acd59(14)
      brack=acd59(49)+acd59(50)+acd59(51)+acd59(52)+acd59(53)+acd59(54)+acd59(5&
      &5)+acd59(56)+acd59(57)+acd59(58)+acd59(59)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd59h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(62) :: acd59
      complex(ki) :: brack
      acd59(1)=d(iv1,iv2)
      acd59(2)=dotproduct(qshift,spvak2k4)
      acd59(3)=abb59(28)
      acd59(4)=dotproduct(qshift,spvak5k1)
      acd59(5)=abb59(5)
      acd59(6)=dotproduct(qshift,spvak5k4)
      acd59(7)=abb59(24)
      acd59(8)=dotproduct(qshift,spvae6k1)
      acd59(9)=abb59(19)
      acd59(10)=abb59(17)
      acd59(11)=k1(iv1)
      acd59(12)=e6(iv2)
      acd59(13)=abb59(11)
      acd59(14)=spvak2k4(iv2)
      acd59(15)=abb59(14)
      acd59(16)=spvak5k1(iv2)
      acd59(17)=abb59(27)
      acd59(18)=spvak5k4(iv2)
      acd59(19)=abb59(13)
      acd59(20)=spvae6k1(iv2)
      acd59(21)=abb59(20)
      acd59(22)=k1(iv2)
      acd59(23)=e6(iv1)
      acd59(24)=spvak2k4(iv1)
      acd59(25)=spvak5k1(iv1)
      acd59(26)=spvak5k4(iv1)
      acd59(27)=spvae6k1(iv1)
      acd59(28)=k2(iv1)
      acd59(29)=k2(iv2)
      acd59(30)=k6(iv1)
      acd59(31)=k6(iv2)
      acd59(32)=abb59(8)
      acd59(33)=abb59(9)
      acd59(34)=dotproduct(qshift,spvak2k1)
      acd59(35)=abb59(12)
      acd59(36)=spvak2k1(iv2)
      acd59(37)=abb59(16)
      acd59(38)=spvak6k1(iv2)
      acd59(39)=abb59(21)
      acd59(40)=spvak2k1(iv1)
      acd59(41)=spvak6k1(iv1)
      acd59(42)=qshift(iv1)
      acd59(43)=qshift(iv2)
      acd59(44)=abb59(25)
      acd59(45)=dotproduct(e6,qshift)
      acd59(46)=abb59(6)
      acd59(47)=abb59(26)
      acd59(48)=spvak2k6(iv2)
      acd59(49)=abb59(22)
      acd59(50)=spvak2k6(iv1)
      acd59(51)=acd59(49)*acd59(50)
      acd59(52)=acd59(41)*acd59(47)
      acd59(53)=2.0_ki*acd59(42)
      acd59(54)=-acd59(7)*acd59(53)
      acd59(55)=acd59(40)*acd59(46)
      acd59(56)=acd59(27)*acd59(44)
      acd59(57)=acd59(11)*acd59(19)
      acd59(58)=acd59(35)*acd59(40)
      acd59(59)=-acd59(45)*acd59(58)
      acd59(60)=acd59(35)*acd59(34)
      acd59(61)=-acd59(23)*acd59(60)
      acd59(51)=acd59(61)+acd59(59)+acd59(57)+acd59(56)+acd59(55)+acd59(54)+acd&
      &59(51)+acd59(52)
      acd59(51)=acd59(18)*acd59(51)
      acd59(52)=acd59(49)*acd59(48)
      acd59(54)=acd59(38)*acd59(47)
      acd59(55)=2.0_ki*acd59(7)
      acd59(55)=-acd59(43)*acd59(55)
      acd59(56)=acd59(36)*acd59(46)
      acd59(57)=acd59(20)*acd59(44)
      acd59(59)=acd59(22)*acd59(19)
      acd59(61)=acd59(35)*acd59(36)
      acd59(62)=-acd59(45)*acd59(61)
      acd59(52)=acd59(62)+acd59(59)+acd59(57)+acd59(56)+acd59(55)+acd59(52)+acd&
      &59(54)
      acd59(52)=acd59(26)*acd59(52)
      acd59(54)=-acd59(11)-acd59(28)
      acd59(54)=acd59(13)*acd59(54)
      acd59(55)=acd59(39)*acd59(41)
      acd59(56)=acd59(40)*acd59(37)
      acd59(57)=acd59(25)*acd59(33)
      acd59(59)=acd59(24)*acd59(32)
      acd59(58)=-acd59(6)*acd59(58)
      acd59(60)=-acd59(26)*acd59(60)
      acd59(54)=acd59(60)+acd59(58)+acd59(59)+acd59(57)+acd59(55)+acd59(56)+acd&
      &59(54)
      acd59(54)=acd59(12)*acd59(54)
      acd59(55)=-acd59(22)-acd59(29)
      acd59(55)=acd59(13)*acd59(55)
      acd59(56)=acd59(38)*acd59(39)
      acd59(57)=acd59(36)*acd59(37)
      acd59(58)=acd59(16)*acd59(33)
      acd59(59)=acd59(14)*acd59(32)
      acd59(60)=-acd59(6)*acd59(61)
      acd59(55)=acd59(60)+acd59(59)+acd59(58)+acd59(56)+acd59(57)+acd59(55)
      acd59(55)=acd59(23)*acd59(55)
      acd59(56)=-acd59(9)*acd59(8)
      acd59(57)=-acd59(6)*acd59(7)
      acd59(58)=-acd59(5)*acd59(4)
      acd59(59)=-acd59(3)*acd59(2)
      acd59(56)=acd59(59)+acd59(58)+acd59(57)+acd59(10)+acd59(56)
      acd59(56)=acd59(1)*acd59(56)
      acd59(57)=-acd59(27)*acd59(9)
      acd59(58)=-acd59(25)*acd59(5)
      acd59(59)=-acd59(24)*acd59(3)
      acd59(57)=acd59(59)+acd59(57)+acd59(58)
      acd59(57)=acd59(43)*acd59(57)
      acd59(56)=acd59(57)+acd59(56)
      acd59(57)=acd59(21)*acd59(27)
      acd59(58)=acd59(17)*acd59(25)
      acd59(59)=acd59(15)*acd59(24)
      acd59(57)=acd59(59)+acd59(57)-acd59(58)
      acd59(58)=-acd59(22)+acd59(31)
      acd59(57)=acd59(57)*acd59(58)
      acd59(58)=-acd59(20)*acd59(21)
      acd59(59)=acd59(16)*acd59(17)
      acd59(60)=-acd59(14)*acd59(15)
      acd59(58)=acd59(60)+acd59(58)+acd59(59)
      acd59(58)=acd59(11)*acd59(58)
      acd59(59)=-acd59(9)*acd59(53)
      acd59(60)=acd59(21)*acd59(30)
      acd59(59)=acd59(59)+acd59(60)
      acd59(59)=acd59(20)*acd59(59)
      acd59(60)=-acd59(5)*acd59(53)
      acd59(61)=-acd59(17)*acd59(30)
      acd59(60)=acd59(60)+acd59(61)
      acd59(60)=acd59(16)*acd59(60)
      acd59(53)=-acd59(3)*acd59(53)
      acd59(61)=acd59(15)*acd59(30)
      acd59(53)=acd59(53)+acd59(61)
      acd59(53)=acd59(14)*acd59(53)
      brack=acd59(51)+acd59(52)+acd59(53)+acd59(54)+acd59(55)+2.0_ki*acd59(56)+&
      &acd59(57)+acd59(58)+acd59(59)+acd59(60)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd59h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(32) :: acd59
      complex(ki) :: brack
      acd59(1)=d(iv1,iv2)
      acd59(2)=spvak2k4(iv3)
      acd59(3)=abb59(28)
      acd59(4)=spvak5k1(iv3)
      acd59(5)=abb59(5)
      acd59(6)=spvak5k4(iv3)
      acd59(7)=abb59(24)
      acd59(8)=spvae6k1(iv3)
      acd59(9)=abb59(19)
      acd59(10)=d(iv1,iv3)
      acd59(11)=spvak2k4(iv2)
      acd59(12)=spvak5k1(iv2)
      acd59(13)=spvak5k4(iv2)
      acd59(14)=spvae6k1(iv2)
      acd59(15)=d(iv2,iv3)
      acd59(16)=spvak2k4(iv1)
      acd59(17)=spvak5k1(iv1)
      acd59(18)=spvak5k4(iv1)
      acd59(19)=spvae6k1(iv1)
      acd59(20)=e6(iv1)
      acd59(21)=spvak2k1(iv3)
      acd59(22)=abb59(12)
      acd59(23)=spvak2k1(iv2)
      acd59(24)=e6(iv2)
      acd59(25)=spvak2k1(iv1)
      acd59(26)=e6(iv3)
      acd59(27)=acd59(9)*acd59(19)
      acd59(28)=acd59(7)*acd59(18)
      acd59(29)=acd59(5)*acd59(17)
      acd59(30)=acd59(3)*acd59(16)
      acd59(27)=acd59(30)+acd59(29)+acd59(27)+acd59(28)
      acd59(27)=acd59(15)*acd59(27)
      acd59(28)=acd59(9)*acd59(14)
      acd59(29)=acd59(7)*acd59(13)
      acd59(30)=acd59(5)*acd59(12)
      acd59(31)=acd59(3)*acd59(11)
      acd59(28)=acd59(31)+acd59(30)+acd59(28)+acd59(29)
      acd59(28)=acd59(10)*acd59(28)
      acd59(29)=acd59(9)*acd59(8)
      acd59(30)=acd59(6)*acd59(7)
      acd59(31)=acd59(5)*acd59(4)
      acd59(32)=acd59(3)*acd59(2)
      acd59(29)=acd59(32)+acd59(31)+acd59(29)+acd59(30)
      acd59(29)=acd59(1)*acd59(29)
      acd59(27)=acd59(29)+acd59(27)+acd59(28)
      acd59(28)=acd59(23)*acd59(26)
      acd59(29)=acd59(21)*acd59(24)
      acd59(28)=acd59(28)+acd59(29)
      acd59(28)=acd59(18)*acd59(28)
      acd59(29)=acd59(25)*acd59(26)
      acd59(30)=acd59(20)*acd59(21)
      acd59(29)=acd59(29)+acd59(30)
      acd59(29)=acd59(13)*acd59(29)
      acd59(30)=acd59(24)*acd59(25)
      acd59(31)=acd59(20)*acd59(23)
      acd59(30)=acd59(30)+acd59(31)
      acd59(30)=acd59(6)*acd59(30)
      acd59(28)=acd59(30)+acd59(28)+acd59(29)
      acd59(28)=acd59(22)*acd59(28)
      brack=2.0_ki*acd59(27)+acd59(28)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd59h0
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
!---#[ subroutine reconstruct_d59:
   subroutine     reconstruct_d59(coeffs)
      use p12_cbbar_hepneg_groups, only: tensrec_info_group3
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group3), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_59:
      coeffs%coeffs_59%c0 = derivative(czip)
      coeffs%coeffs_59%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_59%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_59%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_59%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_59%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_59%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_59%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_59%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_59%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_59%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_59%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_59%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_59%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_59%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_59%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_59%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_59%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_59%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_59%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_59%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_59%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_59%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_59%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_59%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_59%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_59%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_59%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_59%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_59%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_59%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_59%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_59%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_59%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_59%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_59:
   end subroutine reconstruct_d59
!---#] subroutine reconstruct_d59:
end module     p12_cbbar_hepneg_d59h0l1d
