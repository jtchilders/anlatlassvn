module     p11_csbar_hepneg_d61h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p11_csbar_hepneg/helicity0d61h0l1d.f90
   ! generator: buildfortran_d.py
   use p11_csbar_hepneg_config, only: ki
   use p11_csbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d61
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd61h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd61
      complex(ki) :: brack
      acd61(1)=dotproduct(k1,qshift)
      acd61(2)=dotproduct(qshift,spvak5k1)
      acd61(3)=abb61(10)
      acd61(4)=dotproduct(qshift,spvak5k4)
      acd61(5)=abb61(28)
      acd61(6)=abb61(20)
      acd61(7)=dotproduct(k2,qshift)
      acd61(8)=dotproduct(e6,qshift)
      acd61(9)=abb61(15)
      acd61(10)=dotproduct(qshift,spvae6k1)
      acd61(11)=abb61(21)
      acd61(12)=abb61(16)
      acd61(13)=dotproduct(k6,qshift)
      acd61(14)=abb61(27)
      acd61(15)=abb61(13)
      acd61(16)=dotproduct(qshift,qshift)
      acd61(17)=abb61(6)
      acd61(18)=abb61(4)
      acd61(19)=dotproduct(qshift,spvak2k1)
      acd61(20)=abb61(26)
      acd61(21)=abb61(24)
      acd61(22)=abb61(5)
      acd61(23)=dotproduct(qshift,spvak2k6)
      acd61(24)=abb61(11)
      acd61(25)=dotproduct(qshift,spvak6k1)
      acd61(26)=abb61(23)
      acd61(27)=abb61(14)
      acd61(28)=abb61(18)
      acd61(29)=abb61(25)
      acd61(30)=abb61(12)
      acd61(31)=abb61(7)
      acd61(32)=abb61(17)
      acd61(33)=abb61(29)
      acd61(34)=abb61(3)
      acd61(35)=abb61(8)
      acd61(36)=abb61(19)
      acd61(37)=abb61(9)
      acd61(38)=-acd61(13)-acd61(7)
      acd61(38)=acd61(9)*acd61(38)
      acd61(39)=acd61(23)*acd61(24)
      acd61(40)=acd61(19)*acd61(22)
      acd61(41)=acd61(25)*acd61(26)
      acd61(42)=acd61(2)*acd61(18)
      acd61(43)=-acd61(16)*acd61(17)
      acd61(44)=acd61(19)*acd61(20)
      acd61(44)=acd61(21)+acd61(44)
      acd61(44)=acd61(4)*acd61(44)
      acd61(38)=acd61(44)+acd61(43)+acd61(42)+acd61(41)+acd61(40)-acd61(27)+acd&
      &61(39)+acd61(38)
      acd61(38)=acd61(8)*acd61(38)
      acd61(39)=acd61(25)*acd61(33)
      acd61(40)=acd61(1)*acd61(5)
      acd61(41)=acd61(13)*acd61(14)
      acd61(42)=acd61(16)*acd61(29)
      acd61(43)=acd61(10)*acd61(32)
      acd61(39)=acd61(43)+acd61(42)+acd61(41)+acd61(40)-acd61(34)+acd61(39)
      acd61(39)=acd61(4)*acd61(39)
      acd61(40)=-acd61(13)+acd61(7)
      acd61(40)=acd61(11)*acd61(40)
      acd61(41)=acd61(16)*acd61(30)
      acd61(40)=acd61(41)-acd61(35)+acd61(40)
      acd61(40)=acd61(10)*acd61(40)
      acd61(41)=-acd61(25)*acd61(36)
      acd61(42)=-acd61(7)*acd61(12)
      acd61(43)=acd61(2)*acd61(3)
      acd61(43)=-acd61(6)+acd61(43)
      acd61(43)=acd61(1)*acd61(43)
      acd61(44)=-acd61(13)*acd61(15)
      acd61(45)=-acd61(2)*acd61(28)
      acd61(45)=acd61(31)+acd61(45)
      acd61(45)=acd61(16)*acd61(45)
      brack=acd61(37)+acd61(38)+acd61(39)+acd61(40)+acd61(41)+acd61(42)+acd61(4&
      &3)+acd61(44)+acd61(45)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd61h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(59) :: acd61
      complex(ki) :: brack
      acd61(1)=k1(iv1)
      acd61(2)=dotproduct(qshift,spvak5k1)
      acd61(3)=abb61(10)
      acd61(4)=dotproduct(qshift,spvak5k4)
      acd61(5)=abb61(28)
      acd61(6)=abb61(20)
      acd61(7)=k2(iv1)
      acd61(8)=dotproduct(e6,qshift)
      acd61(9)=abb61(15)
      acd61(10)=dotproduct(qshift,spvae6k1)
      acd61(11)=abb61(21)
      acd61(12)=abb61(16)
      acd61(13)=k6(iv1)
      acd61(14)=abb61(27)
      acd61(15)=abb61(13)
      acd61(16)=e6(iv1)
      acd61(17)=dotproduct(k2,qshift)
      acd61(18)=dotproduct(k6,qshift)
      acd61(19)=dotproduct(qshift,qshift)
      acd61(20)=abb61(6)
      acd61(21)=abb61(4)
      acd61(22)=dotproduct(qshift,spvak2k1)
      acd61(23)=abb61(26)
      acd61(24)=abb61(24)
      acd61(25)=abb61(5)
      acd61(26)=dotproduct(qshift,spvak2k6)
      acd61(27)=abb61(11)
      acd61(28)=dotproduct(qshift,spvak6k1)
      acd61(29)=abb61(23)
      acd61(30)=abb61(14)
      acd61(31)=qshift(iv1)
      acd61(32)=abb61(18)
      acd61(33)=abb61(25)
      acd61(34)=abb61(12)
      acd61(35)=abb61(7)
      acd61(36)=spvak5k1(iv1)
      acd61(37)=dotproduct(k1,qshift)
      acd61(38)=spvak5k4(iv1)
      acd61(39)=abb61(17)
      acd61(40)=abb61(29)
      acd61(41)=abb61(3)
      acd61(42)=spvae6k1(iv1)
      acd61(43)=abb61(8)
      acd61(44)=spvak2k1(iv1)
      acd61(45)=spvak2k6(iv1)
      acd61(46)=spvak6k1(iv1)
      acd61(47)=abb61(19)
      acd61(48)=-acd61(28)*acd61(29)
      acd61(49)=-acd61(27)*acd61(26)
      acd61(50)=-acd61(22)*acd61(25)
      acd61(51)=-acd61(2)*acd61(21)
      acd61(52)=acd61(19)*acd61(20)
      acd61(53)=acd61(17)+acd61(18)
      acd61(53)=acd61(9)*acd61(53)
      acd61(54)=acd61(22)*acd61(23)
      acd61(54)=acd61(54)+acd61(24)
      acd61(55)=-acd61(4)*acd61(54)
      acd61(48)=acd61(55)+acd61(53)+acd61(52)+acd61(51)+acd61(50)+acd61(49)+acd&
      &61(30)+acd61(48)
      acd61(48)=acd61(16)*acd61(48)
      acd61(49)=-acd61(4)*acd61(23)
      acd61(49)=acd61(49)-acd61(25)
      acd61(49)=acd61(44)*acd61(49)
      acd61(50)=-acd61(27)*acd61(45)
      acd61(51)=-acd61(46)*acd61(29)
      acd61(52)=-acd61(36)*acd61(21)
      acd61(53)=acd61(7)+acd61(13)
      acd61(53)=acd61(9)*acd61(53)
      acd61(55)=2.0_ki*acd61(31)
      acd61(56)=acd61(20)*acd61(55)
      acd61(54)=-acd61(38)*acd61(54)
      acd61(49)=acd61(54)+acd61(56)+acd61(53)+acd61(52)+acd61(51)+acd61(50)+acd&
      &61(49)
      acd61(49)=acd61(8)*acd61(49)
      acd61(50)=-acd61(28)*acd61(40)
      acd61(51)=-acd61(5)*acd61(37)
      acd61(52)=-acd61(18)*acd61(14)
      acd61(53)=-acd61(19)*acd61(33)
      acd61(54)=-acd61(10)*acd61(39)
      acd61(50)=acd61(54)+acd61(53)+acd61(52)+acd61(51)+acd61(41)+acd61(50)
      acd61(50)=acd61(38)*acd61(50)
      acd61(51)=-acd61(46)*acd61(40)
      acd61(52)=-acd61(1)*acd61(5)
      acd61(53)=-acd61(13)*acd61(14)
      acd61(54)=-acd61(42)*acd61(39)
      acd61(56)=-acd61(33)*acd61(55)
      acd61(51)=acd61(56)+acd61(54)+acd61(53)+acd61(51)+acd61(52)
      acd61(51)=acd61(4)*acd61(51)
      acd61(52)=-acd61(3)*acd61(37)
      acd61(53)=acd61(19)*acd61(32)
      acd61(52)=acd61(53)+acd61(52)
      acd61(52)=acd61(36)*acd61(52)
      acd61(53)=-acd61(19)*acd61(34)
      acd61(54)=-acd61(17)+acd61(18)
      acd61(54)=acd61(11)*acd61(54)
      acd61(53)=acd61(54)+acd61(43)+acd61(53)
      acd61(53)=acd61(42)*acd61(53)
      acd61(54)=acd61(2)*acd61(32)
      acd61(56)=-acd61(10)*acd61(34)
      acd61(54)=acd61(56)-acd61(35)+acd61(54)
      acd61(54)=acd61(54)*acd61(55)
      acd61(55)=acd61(46)*acd61(47)
      acd61(56)=acd61(7)*acd61(12)
      acd61(57)=-acd61(2)*acd61(3)
      acd61(57)=acd61(6)+acd61(57)
      acd61(57)=acd61(1)*acd61(57)
      acd61(58)=acd61(13)*acd61(15)
      acd61(59)=-acd61(7)+acd61(13)
      acd61(59)=acd61(10)*acd61(11)*acd61(59)
      brack=acd61(48)+acd61(49)+acd61(50)+acd61(51)+acd61(52)+acd61(53)+acd61(5&
      &4)+acd61(55)+acd61(56)+acd61(57)+acd61(58)+acd61(59)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd61h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(57) :: acd61
      complex(ki) :: brack
      acd61(1)=d(iv1,iv2)
      acd61(2)=dotproduct(e6,qshift)
      acd61(3)=abb61(6)
      acd61(4)=dotproduct(qshift,spvak5k1)
      acd61(5)=abb61(18)
      acd61(6)=dotproduct(qshift,spvak5k4)
      acd61(7)=abb61(25)
      acd61(8)=dotproduct(qshift,spvae6k1)
      acd61(9)=abb61(12)
      acd61(10)=abb61(7)
      acd61(11)=k1(iv1)
      acd61(12)=spvak5k1(iv2)
      acd61(13)=abb61(10)
      acd61(14)=spvak5k4(iv2)
      acd61(15)=abb61(28)
      acd61(16)=k1(iv2)
      acd61(17)=spvak5k1(iv1)
      acd61(18)=spvak5k4(iv1)
      acd61(19)=k2(iv1)
      acd61(20)=e6(iv2)
      acd61(21)=abb61(15)
      acd61(22)=spvae6k1(iv2)
      acd61(23)=abb61(21)
      acd61(24)=k2(iv2)
      acd61(25)=e6(iv1)
      acd61(26)=spvae6k1(iv1)
      acd61(27)=k6(iv1)
      acd61(28)=abb61(27)
      acd61(29)=k6(iv2)
      acd61(30)=qshift(iv2)
      acd61(31)=abb61(4)
      acd61(32)=dotproduct(qshift,spvak2k1)
      acd61(33)=abb61(26)
      acd61(34)=abb61(24)
      acd61(35)=spvak2k1(iv2)
      acd61(36)=abb61(5)
      acd61(37)=spvak2k6(iv2)
      acd61(38)=abb61(11)
      acd61(39)=spvak6k1(iv2)
      acd61(40)=abb61(23)
      acd61(41)=qshift(iv1)
      acd61(42)=spvak2k1(iv1)
      acd61(43)=spvak2k6(iv1)
      acd61(44)=spvak6k1(iv1)
      acd61(45)=abb61(17)
      acd61(46)=abb61(29)
      acd61(47)=acd61(33)*acd61(6)
      acd61(47)=acd61(47)+acd61(36)
      acd61(48)=acd61(35)*acd61(47)
      acd61(49)=acd61(39)*acd61(40)
      acd61(50)=acd61(38)*acd61(37)
      acd61(51)=acd61(12)*acd61(31)
      acd61(52)=2.0_ki*acd61(30)
      acd61(53)=-acd61(3)*acd61(52)
      acd61(54)=-acd61(24)-acd61(29)
      acd61(54)=acd61(21)*acd61(54)
      acd61(55)=acd61(33)*acd61(32)
      acd61(55)=acd61(55)+acd61(34)
      acd61(56)=acd61(14)*acd61(55)
      acd61(48)=acd61(56)+acd61(54)+acd61(53)+acd61(51)+acd61(49)+acd61(50)+acd&
      &61(48)
      acd61(48)=acd61(25)*acd61(48)
      acd61(47)=acd61(42)*acd61(47)
      acd61(49)=acd61(40)*acd61(44)
      acd61(50)=acd61(38)*acd61(43)
      acd61(51)=acd61(17)*acd61(31)
      acd61(53)=2.0_ki*acd61(41)
      acd61(54)=-acd61(3)*acd61(53)
      acd61(56)=-acd61(19)-acd61(27)
      acd61(56)=acd61(21)*acd61(56)
      acd61(55)=acd61(18)*acd61(55)
      acd61(47)=acd61(55)+acd61(56)+acd61(54)+acd61(51)+acd61(49)+acd61(50)+acd&
      &61(47)
      acd61(47)=acd61(20)*acd61(47)
      acd61(49)=acd61(39)*acd61(46)
      acd61(50)=acd61(15)*acd61(16)
      acd61(51)=acd61(29)*acd61(28)
      acd61(54)=acd61(7)*acd61(52)
      acd61(55)=acd61(22)*acd61(45)
      acd61(56)=acd61(33)*acd61(2)
      acd61(57)=acd61(35)*acd61(56)
      acd61(49)=acd61(57)+acd61(55)+acd61(54)+acd61(51)+acd61(49)+acd61(50)
      acd61(49)=acd61(18)*acd61(49)
      acd61(50)=acd61(44)*acd61(46)
      acd61(51)=acd61(11)*acd61(15)
      acd61(54)=acd61(27)*acd61(28)
      acd61(55)=acd61(7)*acd61(53)
      acd61(57)=acd61(26)*acd61(45)
      acd61(56)=acd61(42)*acd61(56)
      acd61(50)=acd61(56)+acd61(57)+acd61(55)+acd61(54)+acd61(50)+acd61(51)
      acd61(50)=acd61(14)*acd61(50)
      acd61(51)=acd61(9)*acd61(8)
      acd61(54)=acd61(6)*acd61(7)
      acd61(55)=-acd61(5)*acd61(4)
      acd61(56)=-acd61(2)*acd61(3)
      acd61(51)=acd61(56)+acd61(55)+acd61(54)+acd61(10)+acd61(51)
      acd61(51)=acd61(1)*acd61(51)
      acd61(54)=acd61(11)*acd61(13)
      acd61(55)=-acd61(5)*acd61(53)
      acd61(54)=acd61(55)+acd61(54)
      acd61(54)=acd61(12)*acd61(54)
      acd61(55)=-acd61(5)*acd61(17)
      acd61(56)=acd61(26)*acd61(9)
      acd61(55)=acd61(56)+acd61(55)
      acd61(52)=acd61(52)*acd61(55)
      acd61(53)=acd61(9)*acd61(53)
      acd61(55)=acd61(19)-acd61(27)
      acd61(55)=acd61(23)*acd61(55)
      acd61(53)=acd61(53)+acd61(55)
      acd61(53)=acd61(22)*acd61(53)
      acd61(55)=acd61(17)*acd61(13)*acd61(16)
      acd61(56)=acd61(24)-acd61(29)
      acd61(56)=acd61(23)*acd61(26)*acd61(56)
      brack=acd61(47)+acd61(48)+acd61(49)+acd61(50)+2.0_ki*acd61(51)+acd61(52)+&
      &acd61(53)+acd61(54)+acd61(55)+acd61(56)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd61h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(29) :: acd61
      complex(ki) :: brack
      acd61(1)=d(iv1,iv2)
      acd61(2)=e6(iv3)
      acd61(3)=abb61(6)
      acd61(4)=spvak5k1(iv3)
      acd61(5)=abb61(18)
      acd61(6)=spvak5k4(iv3)
      acd61(7)=abb61(25)
      acd61(8)=spvae6k1(iv3)
      acd61(9)=abb61(12)
      acd61(10)=d(iv1,iv3)
      acd61(11)=e6(iv2)
      acd61(12)=spvak5k1(iv2)
      acd61(13)=spvak5k4(iv2)
      acd61(14)=spvae6k1(iv2)
      acd61(15)=d(iv2,iv3)
      acd61(16)=e6(iv1)
      acd61(17)=spvak5k1(iv1)
      acd61(18)=spvak5k4(iv1)
      acd61(19)=spvae6k1(iv1)
      acd61(20)=spvak2k1(iv3)
      acd61(21)=abb61(26)
      acd61(22)=spvak2k1(iv2)
      acd61(23)=spvak2k1(iv1)
      acd61(24)=-acd61(9)*acd61(19)
      acd61(25)=-acd61(7)*acd61(18)
      acd61(26)=acd61(5)*acd61(17)
      acd61(27)=acd61(3)*acd61(16)
      acd61(24)=acd61(27)+acd61(26)+acd61(24)+acd61(25)
      acd61(24)=acd61(15)*acd61(24)
      acd61(25)=-acd61(9)*acd61(14)
      acd61(26)=-acd61(7)*acd61(13)
      acd61(27)=acd61(5)*acd61(12)
      acd61(28)=acd61(3)*acd61(11)
      acd61(25)=acd61(28)+acd61(27)+acd61(25)+acd61(26)
      acd61(25)=acd61(10)*acd61(25)
      acd61(26)=-acd61(9)*acd61(8)
      acd61(27)=-acd61(6)*acd61(7)
      acd61(28)=acd61(5)*acd61(4)
      acd61(29)=acd61(2)*acd61(3)
      acd61(26)=acd61(29)+acd61(28)+acd61(26)+acd61(27)
      acd61(26)=acd61(1)*acd61(26)
      acd61(24)=acd61(26)+acd61(24)+acd61(25)
      acd61(25)=-acd61(13)*acd61(16)
      acd61(26)=-acd61(11)*acd61(18)
      acd61(25)=acd61(25)+acd61(26)
      acd61(25)=acd61(20)*acd61(25)
      acd61(26)=-acd61(16)*acd61(22)
      acd61(27)=-acd61(11)*acd61(23)
      acd61(26)=acd61(26)+acd61(27)
      acd61(26)=acd61(6)*acd61(26)
      acd61(27)=-acd61(18)*acd61(22)
      acd61(28)=-acd61(13)*acd61(23)
      acd61(27)=acd61(27)+acd61(28)
      acd61(27)=acd61(2)*acd61(27)
      acd61(25)=acd61(27)+acd61(26)+acd61(25)
      acd61(25)=acd61(21)*acd61(25)
      brack=2.0_ki*acd61(24)+acd61(25)
   end function brack_4
!---#] function brack_4:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3) result(numerator)
      use p11_csbar_hepneg_globalsl1, only: epspow
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_abbrevd61h0
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
!---#[ subroutine reconstruct_d61:
   subroutine     reconstruct_d61(coeffs)
      use p11_csbar_hepneg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 3 case :
      !---[# reconstruct coeffs%coeffs_61:
      coeffs%coeffs_61%c0 = derivative(czip)
      coeffs%coeffs_61%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_61%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_61%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_61%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_61%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_61%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_61%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_61%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_61%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_61%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_61%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_61%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_61%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_61%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_61%c2(1,3) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_61%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_61%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_61%c2(2,3) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_61%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_61%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_61%c2(3,3) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_61%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_61%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_61%c2(4,3) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_61%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_61%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_61%c2(5,3) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_61%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_61%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_61%c2(6,3) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_61%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_61%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_61%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_61%c3(4,1) = -derivative(czip,2,3,4)
      !---#] reconstruct coeffs%coeffs_61:
   end subroutine reconstruct_d61
!---#] subroutine reconstruct_d61:
end module     p11_csbar_hepneg_d61h0l1d
