module     p16_bbarb_hepemg_d101h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity0d101h0l1d.f90
   ! generator: buildfortran_d.py
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   integer, private :: iv3
   integer, private :: iv4
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d101
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd101
      complex(ki) :: brack
      acd101(1)=abb101(20)
      brack=acd101(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(51) :: acd101
      complex(ki) :: brack
      acd101(1)=k1(iv1)
      acd101(2)=abb101(19)
      acd101(3)=k2(iv1)
      acd101(4)=abb101(17)
      acd101(5)=l3(iv1)
      acd101(6)=abb101(50)
      acd101(7)=k6(iv1)
      acd101(8)=abb101(13)
      acd101(9)=e6(iv1)
      acd101(10)=abb101(23)
      acd101(11)=spvak1k2(iv1)
      acd101(12)=abb101(9)
      acd101(13)=spvak1k4(iv1)
      acd101(14)=abb101(21)
      acd101(15)=spvak1k6(iv1)
      acd101(16)=abb101(10)
      acd101(17)=spvak5k2(iv1)
      acd101(18)=abb101(24)
      acd101(19)=spvak5k4(iv1)
      acd101(20)=abb101(39)
      acd101(21)=spvak5k6(iv1)
      acd101(22)=abb101(67)
      acd101(23)=spvak6k2(iv1)
      acd101(24)=abb101(15)
      acd101(25)=spvak6k4(iv1)
      acd101(26)=abb101(63)
      acd101(27)=spvak1e6(iv1)
      acd101(28)=abb101(26)
      acd101(29)=spvae6k2(iv1)
      acd101(30)=abb101(60)
      acd101(31)=spvae6k4(iv1)
      acd101(32)=abb101(28)
      acd101(33)=spvak5e6(iv1)
      acd101(34)=abb101(12)
      acd101(35)=acd101(33)*acd101(34)
      acd101(36)=acd101(31)*acd101(32)
      acd101(37)=acd101(29)*acd101(30)
      acd101(38)=acd101(27)*acd101(28)
      acd101(39)=acd101(25)*acd101(26)
      acd101(40)=acd101(23)*acd101(24)
      acd101(41)=acd101(21)*acd101(22)
      acd101(42)=acd101(19)*acd101(20)
      acd101(43)=acd101(17)*acd101(18)
      acd101(44)=acd101(15)*acd101(16)
      acd101(45)=acd101(13)*acd101(14)
      acd101(46)=acd101(11)*acd101(12)
      acd101(47)=acd101(9)*acd101(10)
      acd101(48)=acd101(7)*acd101(8)
      acd101(49)=acd101(5)*acd101(6)
      acd101(50)=acd101(3)*acd101(4)
      acd101(51)=acd101(1)*acd101(2)
      brack=acd101(35)+acd101(36)+acd101(37)+acd101(38)+acd101(39)+acd101(40)+a&
      &cd101(41)+acd101(42)+acd101(43)+acd101(44)+acd101(45)+acd101(46)+acd101(&
      &47)+acd101(48)+acd101(49)+acd101(50)+acd101(51)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(65) :: acd101
      complex(ki) :: brack
      acd101(1)=d(iv1,iv2)
      acd101(2)=abb101(8)
      acd101(3)=k1(iv1)
      acd101(4)=e6(iv2)
      acd101(5)=abb101(53)
      acd101(6)=k1(iv2)
      acd101(7)=e6(iv1)
      acd101(8)=k2(iv1)
      acd101(9)=abb101(58)
      acd101(10)=k2(iv2)
      acd101(11)=k6(iv1)
      acd101(12)=abb101(45)
      acd101(13)=spvak1k4(iv2)
      acd101(14)=abb101(11)
      acd101(15)=spvak5k2(iv2)
      acd101(16)=abb101(6)
      acd101(17)=spvae6k4(iv2)
      acd101(18)=abb101(47)
      acd101(19)=spvak5e6(iv2)
      acd101(20)=abb101(42)
      acd101(21)=k6(iv2)
      acd101(22)=spvak1k4(iv1)
      acd101(23)=spvak5k2(iv1)
      acd101(24)=spvae6k4(iv1)
      acd101(25)=spvak5e6(iv1)
      acd101(26)=abb101(7)
      acd101(27)=abb101(36)
      acd101(28)=spvak1k2(iv2)
      acd101(29)=abb101(14)
      acd101(30)=spvak1k6(iv2)
      acd101(31)=abb101(40)
      acd101(32)=spvak6k2(iv2)
      acd101(33)=abb101(66)
      acd101(34)=spvak1k2(iv1)
      acd101(35)=spvak1k6(iv1)
      acd101(36)=spvak6k2(iv1)
      acd101(37)=abb101(30)
      acd101(38)=abb101(18)
      acd101(39)=spvak5k6(iv2)
      acd101(40)=abb101(29)
      acd101(41)=spvak6k4(iv2)
      acd101(42)=abb101(27)
      acd101(43)=spvak5k6(iv1)
      acd101(44)=spvak6k4(iv1)
      acd101(45)=spvak5k4(iv2)
      acd101(46)=abb101(32)
      acd101(47)=spvak5k4(iv1)
      acd101(48)=abb101(65)
      acd101(49)=spvak1e6(iv2)
      acd101(50)=abb101(37)
      acd101(51)=spvae6k2(iv2)
      acd101(52)=abb101(31)
      acd101(53)=spvak1e6(iv1)
      acd101(54)=spvae6k2(iv1)
      acd101(55)=acd101(32)*acd101(33)
      acd101(56)=acd101(30)*acd101(31)
      acd101(57)=acd101(15)*acd101(27)
      acd101(58)=acd101(13)*acd101(26)
      acd101(59)=acd101(9)*acd101(10)
      acd101(60)=acd101(5)*acd101(6)
      acd101(61)=acd101(28)*acd101(29)
      acd101(62)=acd101(21)*acd101(12)
      acd101(55)=acd101(62)+acd101(61)+acd101(60)+acd101(59)+acd101(58)+acd101(&
      &57)+acd101(55)+acd101(56)
      acd101(55)=acd101(7)*acd101(55)
      acd101(56)=acd101(33)*acd101(36)
      acd101(57)=acd101(31)*acd101(35)
      acd101(58)=acd101(23)*acd101(27)
      acd101(59)=acd101(22)*acd101(26)
      acd101(60)=acd101(9)*acd101(8)
      acd101(61)=acd101(5)*acd101(3)
      acd101(62)=acd101(34)*acd101(29)
      acd101(63)=acd101(11)*acd101(12)
      acd101(56)=acd101(63)+acd101(62)+acd101(61)+acd101(60)+acd101(59)+acd101(&
      &58)+acd101(56)+acd101(57)
      acd101(56)=acd101(4)*acd101(56)
      acd101(57)=-acd101(52)*acd101(51)
      acd101(58)=acd101(50)*acd101(49)
      acd101(59)=acd101(32)*acd101(48)
      acd101(60)=acd101(30)*acd101(46)
      acd101(57)=acd101(60)+acd101(59)+acd101(57)+acd101(58)
      acd101(57)=acd101(47)*acd101(57)
      acd101(58)=-acd101(52)*acd101(54)
      acd101(59)=acd101(50)*acd101(53)
      acd101(60)=acd101(36)*acd101(48)
      acd101(61)=acd101(35)*acd101(46)
      acd101(58)=acd101(61)+acd101(60)+acd101(58)+acd101(59)
      acd101(58)=acd101(45)*acd101(58)
      acd101(59)=acd101(42)*acd101(41)
      acd101(60)=acd101(40)*acd101(39)
      acd101(61)=acd101(19)*acd101(38)
      acd101(62)=acd101(17)*acd101(37)
      acd101(59)=acd101(62)+acd101(61)+acd101(59)+acd101(60)
      acd101(59)=acd101(34)*acd101(59)
      acd101(60)=acd101(42)*acd101(44)
      acd101(61)=acd101(40)*acd101(43)
      acd101(62)=acd101(25)*acd101(38)
      acd101(63)=acd101(24)*acd101(37)
      acd101(60)=acd101(63)+acd101(62)+acd101(60)+acd101(61)
      acd101(60)=acd101(28)*acd101(60)
      acd101(61)=acd101(20)*acd101(25)
      acd101(62)=acd101(18)*acd101(24)
      acd101(63)=acd101(16)*acd101(23)
      acd101(64)=acd101(14)*acd101(22)
      acd101(61)=acd101(64)+acd101(63)+acd101(61)+acd101(62)
      acd101(61)=acd101(21)*acd101(61)
      acd101(62)=acd101(19)*acd101(20)
      acd101(63)=acd101(17)*acd101(18)
      acd101(64)=acd101(15)*acd101(16)
      acd101(65)=acd101(13)*acd101(14)
      acd101(62)=acd101(65)+acd101(64)+acd101(62)+acd101(63)
      acd101(62)=acd101(11)*acd101(62)
      acd101(63)=acd101(1)*acd101(2)
      brack=acd101(55)+acd101(56)+acd101(57)+acd101(58)+acd101(59)+acd101(60)+a&
      &cd101(61)+acd101(62)+2.0_ki*acd101(63)
   end function brack_3
!---#] function brack_3:
!---#[ function brack_4:
   pure function brack_4(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(45) :: acd101
      complex(ki) :: brack
      acd101(1)=d(iv1,iv2)
      acd101(2)=e6(iv3)
      acd101(3)=abb101(25)
      acd101(4)=spvak1k2(iv3)
      acd101(5)=abb101(16)
      acd101(6)=spvak1k4(iv3)
      acd101(7)=abb101(22)
      acd101(8)=spvak5k2(iv3)
      acd101(9)=abb101(46)
      acd101(10)=spvak1e6(iv3)
      acd101(11)=abb101(62)
      acd101(12)=spvae6k2(iv3)
      acd101(13)=abb101(61)
      acd101(14)=spvae6k4(iv3)
      acd101(15)=abb101(33)
      acd101(16)=spvak5e6(iv3)
      acd101(17)=abb101(35)
      acd101(18)=d(iv1,iv3)
      acd101(19)=e6(iv2)
      acd101(20)=spvak1k2(iv2)
      acd101(21)=spvak1k4(iv2)
      acd101(22)=spvak5k2(iv2)
      acd101(23)=spvak1e6(iv2)
      acd101(24)=spvae6k2(iv2)
      acd101(25)=spvae6k4(iv2)
      acd101(26)=spvak5e6(iv2)
      acd101(27)=d(iv2,iv3)
      acd101(28)=e6(iv1)
      acd101(29)=spvak1k2(iv1)
      acd101(30)=spvak1k4(iv1)
      acd101(31)=spvak5k2(iv1)
      acd101(32)=spvak1e6(iv1)
      acd101(33)=spvae6k2(iv1)
      acd101(34)=spvae6k4(iv1)
      acd101(35)=spvak5e6(iv1)
      acd101(36)=acd101(17)*acd101(35)
      acd101(37)=acd101(15)*acd101(34)
      acd101(38)=acd101(13)*acd101(33)
      acd101(39)=acd101(11)*acd101(32)
      acd101(40)=acd101(9)*acd101(31)
      acd101(41)=acd101(7)*acd101(30)
      acd101(42)=acd101(5)*acd101(29)
      acd101(43)=acd101(3)*acd101(28)
      acd101(36)=acd101(43)+acd101(42)+acd101(41)+acd101(40)+acd101(39)+acd101(&
      &38)+acd101(36)+acd101(37)
      acd101(36)=acd101(27)*acd101(36)
      acd101(37)=acd101(17)*acd101(26)
      acd101(38)=acd101(15)*acd101(25)
      acd101(39)=acd101(13)*acd101(24)
      acd101(40)=acd101(11)*acd101(23)
      acd101(41)=acd101(9)*acd101(22)
      acd101(42)=acd101(7)*acd101(21)
      acd101(43)=acd101(5)*acd101(20)
      acd101(44)=acd101(3)*acd101(19)
      acd101(37)=acd101(44)+acd101(43)+acd101(42)+acd101(41)+acd101(40)+acd101(&
      &39)+acd101(37)+acd101(38)
      acd101(37)=acd101(18)*acd101(37)
      acd101(38)=acd101(17)*acd101(16)
      acd101(39)=acd101(15)*acd101(14)
      acd101(40)=acd101(13)*acd101(12)
      acd101(41)=acd101(11)*acd101(10)
      acd101(42)=acd101(9)*acd101(8)
      acd101(43)=acd101(7)*acd101(6)
      acd101(44)=acd101(5)*acd101(4)
      acd101(45)=acd101(3)*acd101(2)
      acd101(38)=acd101(45)+acd101(44)+acd101(43)+acd101(42)+acd101(41)+acd101(&
      &40)+acd101(38)+acd101(39)
      acd101(38)=acd101(1)*acd101(38)
      acd101(36)=acd101(38)+acd101(36)+acd101(37)
      brack=2.0_ki*acd101(36)
   end function brack_4
!---#] function brack_4:
!---#[ function brack_5:
   pure function brack_5(Q, mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd101
      complex(ki) :: brack
      brack=0.0_ki
   end function brack_5
!---#] function brack_5:
!---#[ function derivative:
   function derivative(mu2,i1,i2,i3,i4) result(numerator)
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd101h0
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
!---#[ subroutine reconstruct_d101:
   subroutine     reconstruct_d101(coeffs)
      use p16_bbarb_hepemg_groups, only: tensrec_info_group3
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group3), intent(out) :: coeffs
      ! rank 4 case :
      !---[# reconstruct coeffs%coeffs_101:
      coeffs%coeffs_101%c0 = derivative(czip)
      coeffs%coeffs_101%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_101%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_101%c1(1,3) = derivative(czip,1,1,1)/ 6.0_ki
      coeffs%coeffs_101%c1(1,4) = derivative(czip,1,1,1,1)/ 24.0_ki
      coeffs%coeffs_101%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_101%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_101%c1(2,3) = -derivative(czip,2,2,2)/ 6.0_ki
      coeffs%coeffs_101%c1(2,4) = derivative(czip,2,2,2,2)/ 24.0_ki
      coeffs%coeffs_101%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_101%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_101%c1(3,3) = -derivative(czip,3,3,3)/ 6.0_ki
      coeffs%coeffs_101%c1(3,4) = derivative(czip,3,3,3,3)/ 24.0_ki
      coeffs%coeffs_101%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_101%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_101%c1(4,3) = -derivative(czip,4,4,4)/ 6.0_ki
      coeffs%coeffs_101%c1(4,4) = derivative(czip,4,4,4,4)/ 24.0_ki
      coeffs%coeffs_101%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_101%c2(1,2) = derivative(czip,1,2,2)/ 2.0_ki
      coeffs%coeffs_101%c2(1,3) = -derivative(czip,1,2,2,2)/ 6.0_ki
      coeffs%coeffs_101%c2(1,4) = -derivative(czip,1,1,2)/ 2.0_ki
      coeffs%coeffs_101%c2(1,5) = derivative(czip,1,1,2,2)/ 4.0_ki
      coeffs%coeffs_101%c2(1,6) = -derivative(czip,1,1,1,2)/ 6.0_ki
      coeffs%coeffs_101%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_101%c2(2,2) = derivative(czip,1,3,3)/ 2.0_ki
      coeffs%coeffs_101%c2(2,3) = -derivative(czip,1,3,3,3)/ 6.0_ki
      coeffs%coeffs_101%c2(2,4) = -derivative(czip,1,1,3)/ 2.0_ki
      coeffs%coeffs_101%c2(2,5) = derivative(czip,1,1,3,3)/ 4.0_ki
      coeffs%coeffs_101%c2(2,6) = -derivative(czip,1,1,1,3)/ 6.0_ki
      coeffs%coeffs_101%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_101%c2(3,2) = derivative(czip,1,4,4)/ 2.0_ki
      coeffs%coeffs_101%c2(3,3) = -derivative(czip,1,4,4,4)/ 6.0_ki
      coeffs%coeffs_101%c2(3,4) = -derivative(czip,1,1,4)/ 2.0_ki
      coeffs%coeffs_101%c2(3,5) = derivative(czip,1,1,4,4)/ 4.0_ki
      coeffs%coeffs_101%c2(3,6) = -derivative(czip,1,1,1,4)/ 6.0_ki
      coeffs%coeffs_101%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_101%c2(4,2) = -derivative(czip,2,3,3)/ 2.0_ki
      coeffs%coeffs_101%c2(4,3) = derivative(czip,2,3,3,3)/ 6.0_ki
      coeffs%coeffs_101%c2(4,4) = -derivative(czip,2,2,3)/ 2.0_ki
      coeffs%coeffs_101%c2(4,5) = derivative(czip,2,2,3,3)/ 4.0_ki
      coeffs%coeffs_101%c2(4,6) = derivative(czip,2,2,2,3)/ 6.0_ki
      coeffs%coeffs_101%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_101%c2(5,2) = -derivative(czip,2,4,4)/ 2.0_ki
      coeffs%coeffs_101%c2(5,3) = derivative(czip,2,4,4,4)/ 6.0_ki
      coeffs%coeffs_101%c2(5,4) = -derivative(czip,2,2,4)/ 2.0_ki
      coeffs%coeffs_101%c2(5,5) = derivative(czip,2,2,4,4)/ 4.0_ki
      coeffs%coeffs_101%c2(5,6) = derivative(czip,2,2,2,4)/ 6.0_ki
      coeffs%coeffs_101%c2(6,1) = derivative(czip,3,4)
      coeffs%coeffs_101%c2(6,2) = -derivative(czip,3,4,4)/ 2.0_ki
      coeffs%coeffs_101%c2(6,3) = derivative(czip,3,4,4,4)/ 6.0_ki
      coeffs%coeffs_101%c2(6,4) = -derivative(czip,3,3,4)/ 2.0_ki
      coeffs%coeffs_101%c2(6,5) = derivative(czip,3,3,4,4)/ 4.0_ki
      coeffs%coeffs_101%c2(6,6) = derivative(czip,3,3,3,4)/ 6.0_ki
      coeffs%coeffs_101%c3(1,1) = derivative(czip,1,2,3)
      coeffs%coeffs_101%c3(1,2) = -derivative(czip,1,2,3,3)/ 2.0_ki
      coeffs%coeffs_101%c3(1,3) = -derivative(czip,1,2,2,3)/ 2.0_ki
      coeffs%coeffs_101%c3(1,4) = derivative(czip,1,1,2,3)/ 2.0_ki
      coeffs%coeffs_101%c3(2,1) = derivative(czip,1,2,4)
      coeffs%coeffs_101%c3(2,2) = -derivative(czip,1,2,4,4)/ 2.0_ki
      coeffs%coeffs_101%c3(2,3) = -derivative(czip,1,2,2,4)/ 2.0_ki
      coeffs%coeffs_101%c3(2,4) = derivative(czip,1,1,2,4)/ 2.0_ki
      coeffs%coeffs_101%c3(3,1) = derivative(czip,1,3,4)
      coeffs%coeffs_101%c3(3,2) = -derivative(czip,1,3,4,4)/ 2.0_ki
      coeffs%coeffs_101%c3(3,3) = -derivative(czip,1,3,3,4)/ 2.0_ki
      coeffs%coeffs_101%c3(3,4) = derivative(czip,1,1,3,4)/ 2.0_ki
      coeffs%coeffs_101%c3(4,1) = -derivative(czip,2,3,4)
      coeffs%coeffs_101%c3(4,2) = derivative(czip,2,3,4,4)/ 2.0_ki
      coeffs%coeffs_101%c3(4,3) = derivative(czip,2,3,3,4)/ 2.0_ki
      coeffs%coeffs_101%c3(4,4) = derivative(czip,2,2,3,4)/ 2.0_ki
      coeffs%coeffs_101%c4(1,1) = -derivative(czip,1,2,3,4)
      !---#] reconstruct coeffs%coeffs_101:
   end subroutine reconstruct_d101
!---#] subroutine reconstruct_d101:
end module     p16_bbarb_hepemg_d101h0l1d
