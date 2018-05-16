module     p12_sbars_hepemg_d77h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity3d77h3l1d.f90
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
   public :: derivative , reconstruct_d77
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd77h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd77
      complex(ki) :: brack
      acd77(1)=abb77(12)
      brack=acd77(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd77h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(69) :: acd77
      complex(ki) :: brack
      acd77(1)=k1(iv1)
      acd77(2)=abb77(62)
      acd77(3)=k2(iv1)
      acd77(4)=abb77(68)
      acd77(5)=l3(iv1)
      acd77(6)=abb77(42)
      acd77(7)=k6(iv1)
      acd77(8)=abb77(8)
      acd77(9)=e6(iv1)
      acd77(10)=abb77(21)
      acd77(11)=spvak2k1(iv1)
      acd77(12)=abb77(11)
      acd77(13)=spvak2l3(iv1)
      acd77(14)=abb77(37)
      acd77(15)=spvak2k5(iv1)
      acd77(16)=abb77(6)
      acd77(17)=spval3k1(iv1)
      acd77(18)=abb77(36)
      acd77(19)=spval3k5(iv1)
      acd77(20)=abb77(31)
      acd77(21)=spvak4k1(iv1)
      acd77(22)=abb77(9)
      acd77(23)=spvak4l3(iv1)
      acd77(24)=abb77(34)
      acd77(25)=spvak4k5(iv1)
      acd77(26)=abb77(88)
      acd77(27)=spvak4k6(iv1)
      acd77(28)=abb77(58)
      acd77(29)=spvak6k5(iv1)
      acd77(30)=abb77(16)
      acd77(31)=spvae6k1(iv1)
      acd77(32)=abb77(7)
      acd77(33)=spvak2e6(iv1)
      acd77(34)=abb77(13)
      acd77(35)=spval3e6(iv1)
      acd77(36)=abb77(65)
      acd77(37)=spvae6l3(iv1)
      acd77(38)=abb77(82)
      acd77(39)=spvak4e6(iv1)
      acd77(40)=abb77(80)
      acd77(41)=spvae6k5(iv1)
      acd77(42)=abb77(32)
      acd77(43)=spvak6e6(iv1)
      acd77(44)=abb77(76)
      acd77(45)=spvae6k6(iv1)
      acd77(46)=abb77(64)
      acd77(47)=-acd77(2)*acd77(1)
      acd77(48)=-acd77(4)*acd77(3)
      acd77(49)=-acd77(6)*acd77(5)
      acd77(50)=-acd77(8)*acd77(7)
      acd77(51)=-acd77(10)*acd77(9)
      acd77(52)=-acd77(12)*acd77(11)
      acd77(53)=-acd77(14)*acd77(13)
      acd77(54)=-acd77(16)*acd77(15)
      acd77(55)=-acd77(18)*acd77(17)
      acd77(56)=-acd77(20)*acd77(19)
      acd77(57)=-acd77(22)*acd77(21)
      acd77(58)=-acd77(24)*acd77(23)
      acd77(59)=-acd77(26)*acd77(25)
      acd77(60)=-acd77(28)*acd77(27)
      acd77(61)=-acd77(30)*acd77(29)
      acd77(62)=-acd77(32)*acd77(31)
      acd77(63)=-acd77(34)*acd77(33)
      acd77(64)=-acd77(36)*acd77(35)
      acd77(65)=-acd77(38)*acd77(37)
      acd77(66)=-acd77(40)*acd77(39)
      acd77(67)=-acd77(42)*acd77(41)
      acd77(68)=-acd77(44)*acd77(43)
      acd77(69)=-acd77(46)*acd77(45)
      brack=acd77(47)+acd77(48)+acd77(49)+acd77(50)+acd77(51)+acd77(52)+acd77(5&
      &3)+acd77(54)+acd77(55)+acd77(56)+acd77(57)+acd77(58)+acd77(59)+acd77(60)&
      &+acd77(61)+acd77(62)+acd77(63)+acd77(64)+acd77(65)+acd77(66)+acd77(67)+a&
      &cd77(68)+acd77(69)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd77h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(118) :: acd77
      complex(ki) :: brack
      acd77(1)=d(iv1,iv2)
      acd77(2)=abb77(57)
      acd77(3)=k1(iv1)
      acd77(4)=e6(iv2)
      acd77(5)=abb77(45)
      acd77(6)=spvae6k1(iv2)
      acd77(7)=abb77(61)
      acd77(8)=spvak2e6(iv2)
      acd77(9)=abb77(69)
      acd77(10)=spvak4e6(iv2)
      acd77(11)=abb77(48)
      acd77(12)=spvae6k5(iv2)
      acd77(13)=abb77(60)
      acd77(14)=k1(iv2)
      acd77(15)=e6(iv1)
      acd77(16)=spvae6k1(iv1)
      acd77(17)=spvak2e6(iv1)
      acd77(18)=spvak4e6(iv1)
      acd77(19)=spvae6k5(iv1)
      acd77(20)=k2(iv1)
      acd77(21)=abb77(70)
      acd77(22)=abb77(71)
      acd77(23)=abb77(85)
      acd77(24)=abb77(81)
      acd77(25)=abb77(77)
      acd77(26)=spvak2k5(iv2)
      acd77(27)=abb77(30)
      acd77(28)=spvak4k1(iv2)
      acd77(29)=abb77(78)
      acd77(30)=spvak4k5(iv2)
      acd77(31)=abb77(89)
      acd77(32)=k2(iv2)
      acd77(33)=spvak2k5(iv1)
      acd77(34)=spvak4k1(iv1)
      acd77(35)=spvak4k5(iv1)
      acd77(36)=l3(iv1)
      acd77(37)=abb77(28)
      acd77(38)=abb77(67)
      acd77(39)=l3(iv2)
      acd77(40)=k6(iv1)
      acd77(41)=spvak2k1(iv2)
      acd77(42)=abb77(15)
      acd77(43)=k6(iv2)
      acd77(44)=spvak2k1(iv1)
      acd77(45)=abb77(27)
      acd77(46)=abb77(14)
      acd77(47)=abb77(17)
      acd77(48)=abb77(90)
      acd77(49)=abb77(35)
      acd77(50)=abb77(41)
      acd77(51)=abb77(33)
      acd77(52)=abb77(25)
      acd77(53)=abb77(22)
      acd77(54)=abb77(19)
      acd77(55)=spvak2l3(iv2)
      acd77(56)=abb77(44)
      acd77(57)=spvak2k6(iv2)
      acd77(58)=abb77(39)
      acd77(59)=spval3k1(iv2)
      acd77(60)=abb77(38)
      acd77(61)=spvak6k1(iv2)
      acd77(62)=abb77(91)
      acd77(63)=spvak2l3(iv1)
      acd77(64)=spvak2k6(iv1)
      acd77(65)=spval3k1(iv1)
      acd77(66)=spvak6k1(iv1)
      acd77(67)=spval3k5(iv2)
      acd77(68)=abb77(23)
      acd77(69)=spvak4k2(iv2)
      acd77(70)=abb77(46)
      acd77(71)=spvak4l3(iv2)
      acd77(72)=abb77(18)
      acd77(73)=spvak4k6(iv2)
      acd77(74)=abb77(43)
      acd77(75)=spvak6k5(iv2)
      acd77(76)=abb77(40)
      acd77(77)=spvak1e6(iv2)
      acd77(78)=abb77(10)
      acd77(79)=spvae6k2(iv2)
      acd77(80)=abb77(29)
      acd77(81)=spval3e6(iv2)
      acd77(82)=abb77(47)
      acd77(83)=spvae6l3(iv2)
      acd77(84)=abb77(26)
      acd77(85)=spvak6e6(iv2)
      acd77(86)=abb77(24)
      acd77(87)=spvae6k6(iv2)
      acd77(88)=abb77(20)
      acd77(89)=spval3k5(iv1)
      acd77(90)=spvak4k2(iv1)
      acd77(91)=spvak4l3(iv1)
      acd77(92)=spvak4k6(iv1)
      acd77(93)=spvak6k5(iv1)
      acd77(94)=spvak1e6(iv1)
      acd77(95)=spvae6k2(iv1)
      acd77(96)=spval3e6(iv1)
      acd77(97)=spvae6l3(iv1)
      acd77(98)=spvak6e6(iv1)
      acd77(99)=spvae6k6(iv1)
      acd77(100)=acd77(88)*acd77(87)
      acd77(101)=acd77(86)*acd77(85)
      acd77(102)=acd77(84)*acd77(83)
      acd77(103)=-acd77(82)*acd77(81)
      acd77(104)=acd77(80)*acd77(79)
      acd77(105)=acd77(78)*acd77(77)
      acd77(106)=acd77(76)*acd77(75)
      acd77(107)=acd77(74)*acd77(73)
      acd77(108)=acd77(72)*acd77(71)
      acd77(109)=acd77(70)*acd77(69)
      acd77(110)=acd77(68)*acd77(67)
      acd77(111)=acd77(43)*acd77(42)
      acd77(112)=acd77(12)*acd77(53)
      acd77(113)=acd77(10)*acd77(52)
      acd77(114)=acd77(26)*acd77(54)
      acd77(115)=acd77(8)*acd77(51)
      acd77(116)=acd77(6)*acd77(49)
      acd77(117)=acd77(4)*acd77(47)
      acd77(100)=acd77(117)+acd77(116)+acd77(115)+acd77(114)+acd77(113)+acd77(1&
      &12)+acd77(111)+acd77(110)+acd77(109)+acd77(108)+acd77(107)+acd77(106)+ac&
      &d77(105)+acd77(104)+acd77(103)+acd77(102)+acd77(100)+acd77(101)
      acd77(100)=acd77(44)*acd77(100)
      acd77(101)=acd77(88)*acd77(99)
      acd77(102)=acd77(86)*acd77(98)
      acd77(103)=acd77(84)*acd77(97)
      acd77(104)=-acd77(82)*acd77(96)
      acd77(105)=acd77(80)*acd77(95)
      acd77(106)=acd77(78)*acd77(94)
      acd77(107)=acd77(76)*acd77(93)
      acd77(108)=acd77(74)*acd77(92)
      acd77(109)=acd77(72)*acd77(91)
      acd77(110)=acd77(70)*acd77(90)
      acd77(111)=acd77(68)*acd77(89)
      acd77(112)=acd77(40)*acd77(42)
      acd77(113)=acd77(19)*acd77(53)
      acd77(114)=acd77(18)*acd77(52)
      acd77(115)=acd77(33)*acd77(54)
      acd77(116)=acd77(17)*acd77(51)
      acd77(117)=acd77(16)*acd77(49)
      acd77(118)=acd77(15)*acd77(47)
      acd77(101)=acd77(118)+acd77(117)+acd77(116)+acd77(115)+acd77(114)+acd77(1&
      &13)+acd77(112)+acd77(111)+acd77(110)+acd77(109)+acd77(108)+acd77(107)+ac&
      &d77(106)+acd77(105)+acd77(104)+acd77(103)+acd77(101)+acd77(102)
      acd77(101)=acd77(41)*acd77(101)
      acd77(102)=acd77(34)*acd77(29)
      acd77(103)=acd77(19)*acd77(25)
      acd77(104)=acd77(18)*acd77(24)
      acd77(105)=acd77(33)*acd77(27)
      acd77(106)=acd77(17)*acd77(23)
      acd77(107)=-acd77(16)*acd77(22)
      acd77(108)=acd77(15)*acd77(21)
      acd77(109)=acd77(35)*acd77(31)
      acd77(102)=acd77(109)+acd77(108)+acd77(107)+acd77(106)+acd77(105)+acd77(1&
      &04)+acd77(102)+acd77(103)
      acd77(102)=acd77(32)*acd77(102)
      acd77(103)=acd77(28)*acd77(29)
      acd77(104)=acd77(12)*acd77(25)
      acd77(105)=acd77(10)*acd77(24)
      acd77(106)=acd77(26)*acd77(27)
      acd77(107)=acd77(8)*acd77(23)
      acd77(108)=-acd77(6)*acd77(22)
      acd77(109)=acd77(4)*acd77(21)
      acd77(110)=acd77(30)*acd77(31)
      acd77(103)=acd77(110)+acd77(109)+acd77(108)+acd77(107)+acd77(106)+acd77(1&
      &05)+acd77(103)+acd77(104)
      acd77(103)=acd77(20)*acd77(103)
      acd77(104)=acd77(62)*acd77(61)
      acd77(105)=acd77(60)*acd77(59)
      acd77(106)=acd77(58)*acd77(57)
      acd77(107)=acd77(56)*acd77(55)
      acd77(108)=acd77(8)*acd77(50)
      acd77(109)=acd77(6)*acd77(48)
      acd77(104)=acd77(109)+acd77(108)+acd77(107)+acd77(106)+acd77(104)+acd77(1&
      &05)
      acd77(104)=acd77(35)*acd77(104)
      acd77(105)=acd77(62)*acd77(66)
      acd77(106)=acd77(60)*acd77(65)
      acd77(107)=acd77(58)*acd77(64)
      acd77(108)=acd77(56)*acd77(63)
      acd77(109)=acd77(17)*acd77(50)
      acd77(110)=acd77(16)*acd77(48)
      acd77(105)=acd77(110)+acd77(109)+acd77(108)+acd77(107)+acd77(105)+acd77(1&
      &06)
      acd77(105)=acd77(30)*acd77(105)
      acd77(106)=acd77(13)*acd77(19)
      acd77(107)=acd77(11)*acd77(18)
      acd77(108)=acd77(17)*acd77(9)
      acd77(109)=acd77(16)*acd77(7)
      acd77(106)=-acd77(106)+acd77(107)+acd77(108)-acd77(109)
      acd77(107)=acd77(15)*acd77(5)
      acd77(107)=acd77(107)-acd77(106)
      acd77(107)=acd77(14)*acd77(107)
      acd77(108)=acd77(12)*acd77(13)
      acd77(109)=acd77(10)*acd77(11)
      acd77(110)=acd77(8)*acd77(9)
      acd77(111)=acd77(6)*acd77(7)
      acd77(108)=-acd77(108)+acd77(109)+acd77(110)-acd77(111)
      acd77(109)=acd77(4)*acd77(5)
      acd77(109)=acd77(109)-acd77(108)
      acd77(109)=acd77(3)*acd77(109)
      acd77(110)=acd77(34)*acd77(38)
      acd77(111)=acd77(33)*acd77(37)
      acd77(110)=acd77(110)-acd77(111)
      acd77(106)=acd77(110)+acd77(106)
      acd77(106)=acd77(39)*acd77(106)
      acd77(111)=acd77(28)*acd77(38)
      acd77(112)=acd77(26)*acd77(37)
      acd77(111)=acd77(111)-acd77(112)
      acd77(108)=acd77(111)+acd77(108)
      acd77(108)=acd77(36)*acd77(108)
      acd77(110)=acd77(43)*acd77(110)
      acd77(111)=acd77(40)*acd77(111)
      acd77(112)=acd77(28)*acd77(46)
      acd77(113)=acd77(26)*acd77(45)
      acd77(112)=acd77(112)+acd77(113)
      acd77(112)=acd77(15)*acd77(112)
      acd77(113)=acd77(34)*acd77(46)
      acd77(114)=acd77(33)*acd77(45)
      acd77(113)=acd77(113)+acd77(114)
      acd77(113)=acd77(4)*acd77(113)
      acd77(114)=acd77(1)*acd77(2)
      brack=acd77(100)+acd77(101)+acd77(102)+acd77(103)+acd77(104)+acd77(105)+a&
      &cd77(106)+acd77(107)+acd77(108)+acd77(109)+acd77(110)+acd77(111)+acd77(1&
      &12)+acd77(113)+2.0_ki*acd77(114)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd77h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd77
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd77h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd77
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd77h3
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
!---#[ subroutine reconstruct_d77:
   subroutine     reconstruct_d77(coeffs)
      use p12_sbars_hepemg_groups, only: tensrec_info_group4
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group4), intent(out) :: coeffs
      ! rank 4 case :
      !---[# reconstruct coeffs%coeffs_77:
      coeffs%coeffs_77%c0 = derivative(czip)
      coeffs%coeffs_77%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_77%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_77%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_77%c1(1,4) = derivative(czip,1,1,1,1)/ 24.0_ki
      coeffs%coeffs_77%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_77%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_77%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_77%c1(2,4) = derivative(czip,2,2,2,2)/ 24.0_ki
      coeffs%coeffs_77%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_77%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_77%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_77%c1(3,4) = derivative(czip,3,3,3,3)/ 24.0_ki
      coeffs%coeffs_77%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_77%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_77%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_77%c1(4,4) = derivative(czip,4,4,4,4)/ 24.0_ki
      coeffs%coeffs_77%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_77%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_77%c2(1,3) = -derivative(czip,1,2,2,2)/ 6.0_ki
      coeffs%coeffs_77%c2(1,4) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_77%c2(1,5) = derivative(czip,1,1,2,2)/ 4.0_ki
      coeffs%coeffs_77%c2(1,6) = -derivative(czip,1,1,1,2)/ 6.0_ki
      coeffs%coeffs_77%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_77%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_77%c2(2,3) = -derivative(czip,1,3,3,3)/ 6.0_ki
      coeffs%coeffs_77%c2(2,4) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_77%c2(2,5) = derivative(czip,1,1,3,3)/ 4.0_ki
      coeffs%coeffs_77%c2(2,6) = -derivative(czip,1,1,1,3)/ 6.0_ki
      coeffs%coeffs_77%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_77%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_77%c2(3,3) = -derivative(czip,1,4,4,4)/ 6.0_ki
      coeffs%coeffs_77%c2(3,4) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_77%c2(3,5) = derivative(czip,1,1,4,4)/ 4.0_ki
      coeffs%coeffs_77%c2(3,6) = -derivative(czip,1,1,1,4)/ 6.0_ki
      coeffs%coeffs_77%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_77%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_77%c2(4,3) = derivative(czip,2,3,3,3)/ 6.0_ki
      coeffs%coeffs_77%c2(4,4) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_77%c2(4,5) = derivative(czip,2,2,3,3)/ 4.0_ki
      coeffs%coeffs_77%c2(4,6) = derivative(czip,2,2,2,3)/ 6.0_ki
      coeffs%coeffs_77%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_77%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_77%c2(5,3) = derivative(czip,2,4,4,4)/ 6.0_ki
      coeffs%coeffs_77%c2(5,4) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_77%c2(5,5) = derivative(czip,2,2,4,4)/ 4.0_ki
      coeffs%coeffs_77%c2(5,6) = derivative(czip,2,2,2,4)/ 6.0_ki
      coeffs%coeffs_77%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_77%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_77%c2(6,3) = derivative(czip,3,4,4,4)/ 6.0_ki
      coeffs%coeffs_77%c2(6,4) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_77%c2(6,5) = derivative(czip,3,3,4,4)/ 4.0_ki
      coeffs%coeffs_77%c2(6,6) = derivative(czip,3,3,3,4)/ 6.0_ki
      coeffs%coeffs_77%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_77%c3(1,2) = -derivative(czip,1,2,3,3)/ 2.0_ki
      coeffs%coeffs_77%c3(1,3) = -derivative(czip,1,2,2,3)/ 2.0_ki
      coeffs%coeffs_77%c3(1,4) = derivative(czip,1,1,2,3)/ 2.0_ki
      coeffs%coeffs_77%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_77%c3(2,2) = -derivative(czip,1,2,4,4)/ 2.0_ki
      coeffs%coeffs_77%c3(2,3) = -derivative(czip,1,2,2,4)/ 2.0_ki
      coeffs%coeffs_77%c3(2,4) = derivative(czip,1,1,2,4)/ 2.0_ki
      coeffs%coeffs_77%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_77%c3(3,2) = -derivative(czip,1,3,4,4)/ 2.0_ki
      coeffs%coeffs_77%c3(3,3) = -derivative(czip,1,3,3,4)/ 2.0_ki
      coeffs%coeffs_77%c3(3,4) = derivative(czip,1,1,3,4)/ 2.0_ki
      coeffs%coeffs_77%c3(4,1) = -derivative(czip,2,3,4)
      coeffs%coeffs_77%c3(4,2) = derivative(czip,2,3,4,4)/ 2.0_ki
      coeffs%coeffs_77%c3(4,3) = derivative(czip,2,3,3,4)/ 2.0_ki
      coeffs%coeffs_77%c3(4,4) = derivative(czip,2,2,3,4)/ 2.0_ki
      coeffs%coeffs_77%c4(1,1) = -derivative(czip,1,2,3,4)
      !---#] reconstruct coeffs%coeffs_77:
   end subroutine reconstruct_d77
!---#] subroutine reconstruct_d77:
end module     p12_sbars_hepemg_d77h3l1d
