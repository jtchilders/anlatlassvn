module     p5_usbar_hepneg_d60h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p5_usbar_hepneg/helicity0d60h0l1d.f90
   ! generator: buildfortran_d.py
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d60
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(43) :: acd60
      complex(ki) :: brack
      acd60(1)=dotproduct(k2,qshift)
      acd60(2)=dotproduct(e6,qshift)
      acd60(3)=abb60(5)
      acd60(4)=dotproduct(qshift,spvak5k1)
      acd60(5)=abb60(24)
      acd60(6)=dotproduct(qshift,spvak5k4)
      acd60(7)=abb60(27)
      acd60(8)=dotproduct(qshift,spvae6k1)
      acd60(9)=abb60(20)
      acd60(10)=dotproduct(qshift,spvak2e6)
      acd60(11)=abb60(19)
      acd60(12)=abb60(4)
      acd60(13)=dotproduct(k6,qshift)
      acd60(14)=abb60(16)
      acd60(15)=dotproduct(qshift,qshift)
      acd60(16)=abb60(18)
      acd60(17)=abb60(13)
      acd60(18)=dotproduct(qshift,spvak2k1)
      acd60(19)=abb60(11)
      acd60(20)=abb60(17)
      acd60(21)=dotproduct(qshift,spvak2k6)
      acd60(22)=abb60(12)
      acd60(23)=dotproduct(qshift,spvak6k1)
      acd60(24)=abb60(23)
      acd60(25)=abb60(10)
      acd60(26)=abb60(7)
      acd60(27)=abb60(21)
      acd60(28)=abb60(14)
      acd60(29)=abb60(15)
      acd60(30)=abb60(6)
      acd60(31)=abb60(30)
      acd60(32)=abb60(28)
      acd60(33)=abb60(8)
      acd60(34)=abb60(25)
      acd60(35)=abb60(22)
      acd60(36)=abb60(9)
      acd60(37)=acd60(6)*acd60(19)
      acd60(37)=acd60(37)+acd60(20)
      acd60(37)=acd60(18)*acd60(37)
      acd60(38)=acd60(21)*acd60(22)
      acd60(39)=acd60(23)*acd60(24)
      acd60(40)=-acd60(15)*acd60(16)
      acd60(41)=acd60(4)*acd60(17)
      acd60(42)=acd60(1)*acd60(3)
      acd60(37)=acd60(42)+acd60(41)+acd60(40)+acd60(39)-acd60(25)+acd60(38)+acd&
      &60(37)
      acd60(37)=acd60(2)*acd60(37)
      acd60(38)=acd60(23)*acd60(33)
      acd60(39)=acd60(10)*acd60(32)
      acd60(40)=acd60(8)*acd60(31)
      acd60(38)=acd60(40)+acd60(39)-acd60(34)+acd60(38)
      acd60(38)=acd60(6)*acd60(38)
      acd60(39)=acd60(10)*acd60(11)
      acd60(40)=acd60(8)*acd60(9)
      acd60(39)=acd60(39)-acd60(40)
      acd60(40)=-acd60(4)*acd60(5)
      acd60(41)=acd60(6)*acd60(7)
      acd60(40)=acd60(41)+acd60(40)-acd60(12)-acd60(39)
      acd60(40)=acd60(1)*acd60(40)
      acd60(39)=-acd60(14)+acd60(39)
      acd60(39)=acd60(13)*acd60(39)
      acd60(41)=-acd60(10)*acd60(28)
      acd60(42)=-acd60(8)*acd60(27)
      acd60(41)=acd60(42)+acd60(29)+acd60(41)
      acd60(41)=acd60(15)*acd60(41)
      acd60(42)=acd60(13)*acd60(5)
      acd60(43)=-acd60(15)*acd60(26)
      acd60(42)=acd60(43)-acd60(30)+acd60(42)
      acd60(42)=acd60(4)*acd60(42)
      acd60(43)=-acd60(23)*acd60(35)
      brack=acd60(36)+acd60(37)+acd60(38)+acd60(39)+acd60(40)+acd60(41)+acd60(4&
      &2)+acd60(43)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(58) :: acd60
      complex(ki) :: brack
      acd60(1)=k2(iv1)
      acd60(2)=dotproduct(e6,qshift)
      acd60(3)=abb60(5)
      acd60(4)=dotproduct(qshift,spvak5k1)
      acd60(5)=abb60(24)
      acd60(6)=dotproduct(qshift,spvak5k4)
      acd60(7)=abb60(27)
      acd60(8)=dotproduct(qshift,spvae6k1)
      acd60(9)=abb60(20)
      acd60(10)=dotproduct(qshift,spvak2e6)
      acd60(11)=abb60(19)
      acd60(12)=abb60(4)
      acd60(13)=k6(iv1)
      acd60(14)=abb60(16)
      acd60(15)=e6(iv1)
      acd60(16)=dotproduct(k2,qshift)
      acd60(17)=dotproduct(qshift,qshift)
      acd60(18)=abb60(18)
      acd60(19)=abb60(13)
      acd60(20)=dotproduct(qshift,spvak2k1)
      acd60(21)=abb60(11)
      acd60(22)=abb60(17)
      acd60(23)=dotproduct(qshift,spvak2k6)
      acd60(24)=abb60(12)
      acd60(25)=dotproduct(qshift,spvak6k1)
      acd60(26)=abb60(23)
      acd60(27)=abb60(10)
      acd60(28)=qshift(iv1)
      acd60(29)=abb60(7)
      acd60(30)=abb60(21)
      acd60(31)=abb60(14)
      acd60(32)=abb60(15)
      acd60(33)=spvak5k1(iv1)
      acd60(34)=dotproduct(k6,qshift)
      acd60(35)=abb60(6)
      acd60(36)=spvak5k4(iv1)
      acd60(37)=abb60(30)
      acd60(38)=abb60(28)
      acd60(39)=abb60(8)
      acd60(40)=abb60(25)
      acd60(41)=spvae6k1(iv1)
      acd60(42)=spvak2e6(iv1)
      acd60(43)=spvak2k1(iv1)
      acd60(44)=spvak2k6(iv1)
      acd60(45)=spvak6k1(iv1)
      acd60(46)=abb60(22)
      acd60(47)=acd60(6)*acd60(21)
      acd60(47)=acd60(47)+acd60(22)
      acd60(48)=acd60(43)*acd60(47)
      acd60(49)=acd60(24)*acd60(44)
      acd60(50)=acd60(45)*acd60(26)
      acd60(51)=acd60(33)*acd60(19)
      acd60(52)=2.0_ki*acd60(28)
      acd60(53)=-acd60(18)*acd60(52)
      acd60(54)=acd60(36)*acd60(20)*acd60(21)
      acd60(55)=acd60(1)*acd60(3)
      acd60(48)=acd60(55)+acd60(54)+acd60(53)+acd60(51)+acd60(50)+acd60(49)+acd&
      &60(48)
      acd60(48)=acd60(2)*acd60(48)
      acd60(47)=acd60(20)*acd60(47)
      acd60(49)=acd60(25)*acd60(26)
      acd60(50)=acd60(24)*acd60(23)
      acd60(51)=-acd60(17)*acd60(18)
      acd60(53)=acd60(4)*acd60(19)
      acd60(54)=acd60(16)*acd60(3)
      acd60(47)=acd60(54)+acd60(53)+acd60(51)+acd60(50)-acd60(27)+acd60(49)+acd&
      &60(47)
      acd60(47)=acd60(15)*acd60(47)
      acd60(49)=acd60(25)*acd60(39)
      acd60(50)=acd60(10)*acd60(38)
      acd60(51)=acd60(8)*acd60(37)
      acd60(53)=acd60(16)*acd60(7)
      acd60(49)=acd60(53)+acd60(51)+acd60(50)-acd60(40)+acd60(49)
      acd60(49)=acd60(36)*acd60(49)
      acd60(50)=acd60(10)*acd60(11)
      acd60(51)=acd60(8)*acd60(9)
      acd60(53)=acd60(4)*acd60(5)
      acd60(50)=acd60(53)+acd60(50)-acd60(51)
      acd60(51)=-acd60(14)+acd60(50)
      acd60(51)=acd60(13)*acd60(51)
      acd60(53)=-acd60(10)*acd60(31)
      acd60(54)=-acd60(8)*acd60(30)
      acd60(55)=-acd60(4)*acd60(29)
      acd60(53)=acd60(55)+acd60(54)+acd60(32)+acd60(53)
      acd60(52)=acd60(53)*acd60(52)
      acd60(53)=acd60(45)*acd60(39)
      acd60(54)=acd60(42)*acd60(38)
      acd60(55)=acd60(41)*acd60(37)
      acd60(53)=acd60(55)+acd60(53)+acd60(54)
      acd60(53)=acd60(6)*acd60(53)
      acd60(54)=acd60(6)*acd60(7)
      acd60(50)=acd60(54)-acd60(12)-acd60(50)
      acd60(50)=acd60(1)*acd60(50)
      acd60(54)=acd60(11)*acd60(42)
      acd60(55)=acd60(9)*acd60(41)
      acd60(54)=acd60(54)-acd60(55)
      acd60(55)=acd60(34)*acd60(54)
      acd60(56)=-acd60(42)*acd60(31)
      acd60(57)=-acd60(41)*acd60(30)
      acd60(56)=acd60(56)+acd60(57)
      acd60(56)=acd60(17)*acd60(56)
      acd60(57)=-acd60(17)*acd60(29)
      acd60(58)=acd60(5)*acd60(34)
      acd60(57)=acd60(58)-acd60(35)+acd60(57)
      acd60(57)=acd60(33)*acd60(57)
      acd60(58)=-acd60(33)*acd60(5)
      acd60(54)=acd60(58)-acd60(54)
      acd60(54)=acd60(16)*acd60(54)
      acd60(58)=-acd60(45)*acd60(46)
      brack=acd60(47)+acd60(48)+acd60(49)+acd60(50)+acd60(51)+acd60(52)+acd60(5&
      &3)+acd60(54)+acd60(55)+acd60(56)+acd60(57)+acd60(58)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(58) :: acd60
      complex(ki) :: brack
      acd60(1)=d(iv1,iv2)
      acd60(2)=dotproduct(e6,qshift)
      acd60(3)=abb60(18)
      acd60(4)=dotproduct(qshift,spvak5k1)
      acd60(5)=abb60(7)
      acd60(6)=dotproduct(qshift,spvae6k1)
      acd60(7)=abb60(21)
      acd60(8)=dotproduct(qshift,spvak2e6)
      acd60(9)=abb60(14)
      acd60(10)=abb60(15)
      acd60(11)=k2(iv1)
      acd60(12)=e6(iv2)
      acd60(13)=abb60(5)
      acd60(14)=spvak5k1(iv2)
      acd60(15)=abb60(24)
      acd60(16)=spvae6k1(iv2)
      acd60(17)=abb60(20)
      acd60(18)=spvak2e6(iv2)
      acd60(19)=abb60(19)
      acd60(20)=spvak5k4(iv2)
      acd60(21)=abb60(27)
      acd60(22)=k2(iv2)
      acd60(23)=e6(iv1)
      acd60(24)=spvak5k1(iv1)
      acd60(25)=spvae6k1(iv1)
      acd60(26)=spvak2e6(iv1)
      acd60(27)=spvak5k4(iv1)
      acd60(28)=k6(iv1)
      acd60(29)=k6(iv2)
      acd60(30)=qshift(iv2)
      acd60(31)=abb60(13)
      acd60(32)=dotproduct(qshift,spvak2k1)
      acd60(33)=abb60(11)
      acd60(34)=spvak2k1(iv2)
      acd60(35)=dotproduct(qshift,spvak5k4)
      acd60(36)=abb60(17)
      acd60(37)=spvak2k6(iv2)
      acd60(38)=abb60(12)
      acd60(39)=spvak6k1(iv2)
      acd60(40)=abb60(23)
      acd60(41)=qshift(iv1)
      acd60(42)=spvak2k1(iv1)
      acd60(43)=spvak2k6(iv1)
      acd60(44)=spvak6k1(iv1)
      acd60(45)=abb60(30)
      acd60(46)=abb60(28)
      acd60(47)=abb60(8)
      acd60(48)=acd60(33)*acd60(35)
      acd60(48)=acd60(48)+acd60(36)
      acd60(49)=acd60(34)*acd60(48)
      acd60(50)=acd60(39)*acd60(40)
      acd60(51)=acd60(38)*acd60(37)
      acd60(52)=2.0_ki*acd60(3)
      acd60(52)=-acd60(30)*acd60(52)
      acd60(53)=acd60(14)*acd60(31)
      acd60(54)=acd60(22)*acd60(13)
      acd60(55)=acd60(33)*acd60(32)
      acd60(56)=acd60(20)*acd60(55)
      acd60(49)=acd60(56)+acd60(54)+acd60(53)+acd60(52)+acd60(50)+acd60(51)+acd&
      &60(49)
      acd60(49)=acd60(23)*acd60(49)
      acd60(48)=acd60(42)*acd60(48)
      acd60(50)=acd60(40)*acd60(44)
      acd60(51)=acd60(38)*acd60(43)
      acd60(52)=2.0_ki*acd60(41)
      acd60(53)=-acd60(3)*acd60(52)
      acd60(54)=acd60(24)*acd60(31)
      acd60(56)=acd60(11)*acd60(13)
      acd60(55)=acd60(27)*acd60(55)
      acd60(48)=acd60(55)+acd60(56)+acd60(54)+acd60(53)+acd60(50)+acd60(51)+acd&
      &60(48)
      acd60(48)=acd60(12)*acd60(48)
      acd60(50)=acd60(39)*acd60(47)
      acd60(51)=acd60(18)*acd60(46)
      acd60(53)=acd60(16)*acd60(45)
      acd60(54)=acd60(22)*acd60(21)
      acd60(55)=acd60(33)*acd60(2)
      acd60(56)=acd60(34)*acd60(55)
      acd60(50)=acd60(56)+acd60(54)+acd60(53)+acd60(50)+acd60(51)
      acd60(50)=acd60(27)*acd60(50)
      acd60(51)=acd60(44)*acd60(47)
      acd60(53)=acd60(26)*acd60(46)
      acd60(54)=acd60(25)*acd60(45)
      acd60(56)=acd60(11)*acd60(21)
      acd60(55)=acd60(42)*acd60(55)
      acd60(51)=acd60(55)+acd60(56)+acd60(54)+acd60(51)+acd60(53)
      acd60(51)=acd60(20)*acd60(51)
      acd60(53)=-acd60(9)*acd60(8)
      acd60(54)=-acd60(7)*acd60(6)
      acd60(55)=-acd60(5)*acd60(4)
      acd60(56)=-acd60(2)*acd60(3)
      acd60(53)=acd60(56)+acd60(55)+acd60(54)+acd60(10)+acd60(53)
      acd60(53)=acd60(1)*acd60(53)
      acd60(54)=-acd60(26)*acd60(9)
      acd60(55)=-acd60(25)*acd60(7)
      acd60(56)=-acd60(24)*acd60(5)
      acd60(54)=acd60(56)+acd60(54)+acd60(55)
      acd60(54)=acd60(30)*acd60(54)
      acd60(53)=acd60(54)+acd60(53)
      acd60(54)=acd60(19)*acd60(26)
      acd60(55)=acd60(17)*acd60(25)
      acd60(56)=acd60(15)*acd60(24)
      acd60(54)=acd60(56)+acd60(54)-acd60(55)
      acd60(55)=-acd60(22)+acd60(29)
      acd60(54)=acd60(54)*acd60(55)
      acd60(55)=-acd60(18)*acd60(19)
      acd60(56)=acd60(16)*acd60(17)
      acd60(57)=-acd60(14)*acd60(15)
      acd60(55)=acd60(57)+acd60(55)+acd60(56)
      acd60(55)=acd60(11)*acd60(55)
      acd60(56)=-acd60(9)*acd60(52)
      acd60(57)=acd60(19)*acd60(28)
      acd60(56)=acd60(56)+acd60(57)
      acd60(56)=acd60(18)*acd60(56)
      acd60(57)=-acd60(7)*acd60(52)
      acd60(58)=-acd60(17)*acd60(28)
      acd60(57)=acd60(57)+acd60(58)
      acd60(57)=acd60(16)*acd60(57)
      acd60(52)=-acd60(5)*acd60(52)
      acd60(58)=acd60(15)*acd60(28)
      acd60(52)=acd60(52)+acd60(58)
      acd60(52)=acd60(14)*acd60(52)
      brack=acd60(48)+acd60(49)+acd60(50)+acd60(51)+acd60(52)+2.0_ki*acd60(53)+&
      &acd60(54)+acd60(55)+acd60(56)+acd60(57)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(32) :: acd60
      complex(ki) :: brack
      acd60(1)=d(iv1,iv2)
      acd60(2)=e6(iv3)
      acd60(3)=abb60(18)
      acd60(4)=spvak5k1(iv3)
      acd60(5)=abb60(7)
      acd60(6)=spvae6k1(iv3)
      acd60(7)=abb60(21)
      acd60(8)=spvak2e6(iv3)
      acd60(9)=abb60(14)
      acd60(10)=d(iv1,iv3)
      acd60(11)=e6(iv2)
      acd60(12)=spvak5k1(iv2)
      acd60(13)=spvae6k1(iv2)
      acd60(14)=spvak2e6(iv2)
      acd60(15)=d(iv2,iv3)
      acd60(16)=e6(iv1)
      acd60(17)=spvak5k1(iv1)
      acd60(18)=spvae6k1(iv1)
      acd60(19)=spvak2e6(iv1)
      acd60(20)=spvak2k1(iv2)
      acd60(21)=spvak5k4(iv3)
      acd60(22)=abb60(11)
      acd60(23)=spvak2k1(iv3)
      acd60(24)=spvak5k4(iv2)
      acd60(25)=spvak2k1(iv1)
      acd60(26)=spvak5k4(iv1)
      acd60(27)=-acd60(9)*acd60(19)
      acd60(28)=-acd60(7)*acd60(18)
      acd60(29)=-acd60(5)*acd60(17)
      acd60(30)=-acd60(3)*acd60(16)
      acd60(27)=acd60(30)+acd60(29)+acd60(27)+acd60(28)
      acd60(27)=acd60(15)*acd60(27)
      acd60(28)=-acd60(9)*acd60(14)
      acd60(29)=-acd60(7)*acd60(13)
      acd60(30)=-acd60(5)*acd60(12)
      acd60(31)=-acd60(3)*acd60(11)
      acd60(28)=acd60(31)+acd60(30)+acd60(28)+acd60(29)
      acd60(28)=acd60(10)*acd60(28)
      acd60(29)=-acd60(9)*acd60(8)
      acd60(30)=-acd60(7)*acd60(6)
      acd60(31)=-acd60(5)*acd60(4)
      acd60(32)=-acd60(2)*acd60(3)
      acd60(29)=acd60(32)+acd60(31)+acd60(29)+acd60(30)
      acd60(29)=acd60(1)*acd60(29)
      acd60(27)=acd60(29)+acd60(27)+acd60(28)
      acd60(28)=acd60(23)*acd60(24)
      acd60(29)=acd60(20)*acd60(21)
      acd60(28)=acd60(28)+acd60(29)
      acd60(28)=acd60(16)*acd60(28)
      acd60(29)=acd60(23)*acd60(26)
      acd60(30)=acd60(21)*acd60(25)
      acd60(29)=acd60(29)+acd60(30)
      acd60(29)=acd60(11)*acd60(29)
      acd60(30)=acd60(24)*acd60(25)
      acd60(31)=acd60(20)*acd60(26)
      acd60(30)=acd60(30)+acd60(31)
      acd60(30)=acd60(2)*acd60(30)
      acd60(28)=acd60(30)+acd60(28)+acd60(29)
      acd60(28)=acd60(22)*acd60(28)
      brack=2.0_ki*acd60(27)+acd60(28)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd60h0
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
!---#[ subroutine reconstruct_d60:
   subroutine     reconstruct_d60(coeffs)
      use p5_usbar_hepneg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_60:
      coeffs%coeffs_60%c0 = derivative(czip)
      coeffs%coeffs_60%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_60%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_60%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_60%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_60%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_60%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_60%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_60%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_60%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_60%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_60%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_60%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_60%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_60%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_60%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_60%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_60%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_60%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_60%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_60%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_60%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_60%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_60%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_60%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_60%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_60%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_60%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_60%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_60%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_60%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_60%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_60%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_60%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_60%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_60:
   end subroutine reconstruct_d60
!---#] subroutine reconstruct_d60:
end module     p5_usbar_hepneg_d60h0l1d
