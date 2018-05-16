module     p0_dbard_hepemg_d77h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbard_hepemg/helicity0d77h0l1.f90
   ! generator: buildfortran.py
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p0_dbard_hepemg_model
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_color
      use p0_dbard_hepemg_abbrevd77h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc77(92)
      complex(ki) :: Qspk6
      complex(ki) :: Qspe6
      complex(ki) :: Qspvak5l3
      complex(ki) :: Qspvae6k6
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak5e6
      complex(ki) :: Qspvae6k4
      complex(ki) :: Qspval3k4
      complex(ki) :: Qspvak6e6
      complex(ki) :: Qspval3e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvak6k4
      complex(ki) :: Qspvak5k6
      complex(ki) :: Qspvae6l3
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak2k4
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak1k4
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak1l3
      complex(ki) :: Qspval3k2
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspk1
      complex(ki) :: Qspl3
      complex(ki) :: QspQ
      Qspk6 = dotproduct(Q,k6)
      Qspe6 = dotproduct(Q,e6)
      Qspvak5l3 = dotproduct(Q,spvak5l3)
      Qspvae6k6 = dotproduct(Q,spvae6k6)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak5e6 = dotproduct(Q,spvak5e6)
      Qspvae6k4 = dotproduct(Q,spvae6k4)
      Qspval3k4 = dotproduct(Q,spval3k4)
      Qspvak6e6 = dotproduct(Q,spvak6e6)
      Qspval3e6 = dotproduct(Q,spval3e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvak6k4 = dotproduct(Q,spvak6k4)
      Qspvak5k6 = dotproduct(Q,spvak5k6)
      Qspvae6l3 = dotproduct(Q,spvae6l3)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak2k4 = dotproduct(Q,spvak2k4)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak1k4 = dotproduct(Q,spvak1k4)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspk2 = dotproduct(Q,k2)
      Qspvak1l3 = dotproduct(Q,spvak1l3)
      Qspval3k2 = dotproduct(Q,spval3k2)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspk1 = dotproduct(Q,k1)
      Qspl3 = dotproduct(Q,l3)
      QspQ = dotproduct(Q,Q)
      acc77(1)=abb77(6)
      acc77(2)=abb77(7)
      acc77(3)=abb77(8)
      acc77(4)=abb77(9)
      acc77(5)=abb77(10)
      acc77(6)=abb77(11)
      acc77(7)=abb77(12)
      acc77(8)=abb77(13)
      acc77(9)=abb77(14)
      acc77(10)=abb77(15)
      acc77(11)=abb77(16)
      acc77(12)=abb77(17)
      acc77(13)=abb77(18)
      acc77(14)=abb77(19)
      acc77(15)=abb77(20)
      acc77(16)=abb77(21)
      acc77(17)=abb77(22)
      acc77(18)=abb77(23)
      acc77(19)=abb77(24)
      acc77(20)=abb77(25)
      acc77(21)=abb77(26)
      acc77(22)=abb77(27)
      acc77(23)=abb77(28)
      acc77(24)=abb77(29)
      acc77(25)=abb77(30)
      acc77(26)=abb77(31)
      acc77(27)=abb77(32)
      acc77(28)=abb77(33)
      acc77(29)=abb77(34)
      acc77(30)=abb77(35)
      acc77(31)=abb77(36)
      acc77(32)=abb77(37)
      acc77(33)=abb77(38)
      acc77(34)=abb77(39)
      acc77(35)=abb77(40)
      acc77(36)=abb77(41)
      acc77(37)=abb77(42)
      acc77(38)=abb77(43)
      acc77(39)=abb77(44)
      acc77(40)=abb77(45)
      acc77(41)=abb77(47)
      acc77(42)=abb77(48)
      acc77(43)=abb77(53)
      acc77(44)=abb77(57)
      acc77(45)=abb77(60)
      acc77(46)=abb77(61)
      acc77(47)=abb77(62)
      acc77(48)=abb77(63)
      acc77(49)=abb77(65)
      acc77(50)=abb77(67)
      acc77(51)=abb77(68)
      acc77(52)=abb77(69)
      acc77(53)=abb77(70)
      acc77(54)=abb77(71)
      acc77(55)=abb77(72)
      acc77(56)=abb77(77)
      acc77(57)=abb77(78)
      acc77(58)=abb77(79)
      acc77(59)=abb77(80)
      acc77(60)=abb77(81)
      acc77(61)=abb77(84)
      acc77(62)=abb77(85)
      acc77(63)=abb77(87)
      acc77(64)=abb77(88)
      acc77(65)=abb77(89)
      acc77(66)=abb77(90)
      acc77(67)=acc77(6)*Qspk6
      acc77(68)=acc77(10)*Qspe6
      acc77(69)=acc77(12)*Qspvak5l3
      acc77(70)=acc77(13)*Qspvae6k6
      acc77(71)=acc77(14)*Qspvak5k2
      acc77(72)=acc77(15)*Qspvak5e6
      acc77(73)=acc77(17)*Qspvae6k4
      acc77(74)=acc77(18)*Qspval3k4
      acc77(75)=acc77(19)*Qspvak6e6
      acc77(76)=acc77(20)*Qspval3e6
      acc77(77)=acc77(21)*Qspvae6k2
      acc77(78)=acc77(28)*Qspvak1e6
      acc77(79)=acc77(33)*Qspvak6k4
      acc77(80)=acc77(35)*Qspvak5k6
      acc77(81)=acc77(41)*Qspvae6l3
      acc77(82)=Qspvak2e6*acc77(24)
      acc77(83)=Qspvae6k1*acc77(5)
      acc77(84)=Qspvak2k4*acc77(39)
      acc77(67)=acc77(84)+acc77(83)+acc77(82)+acc77(81)+acc77(80)+acc77(79)+acc&
      &77(78)+acc77(77)+acc77(76)+acc77(75)+acc77(74)+acc77(73)+acc77(72)+acc77&
      &(71)+acc77(70)+acc77(69)+acc77(68)+acc77(67)+acc77(1)
      acc77(67)=Qspvak1k2*acc77(67)
      acc77(68)=acc77(32)*Qspvak5k2
      acc77(69)=acc77(53)*Qspe6
      acc77(70)=acc77(54)*Qspvak1e6
      acc77(71)=-acc77(57)*Qspvak1k4
      acc77(72)=acc77(58)*Qspvak5e6
      acc77(73)=acc77(60)*Qspvae6k4
      acc77(74)=acc77(62)*Qspvae6k2
      acc77(75)=acc77(65)*Qspvak5k4
      acc77(68)=acc77(75)+acc77(74)+acc77(73)+acc77(72)+acc77(71)+acc77(70)+acc&
      &77(69)+acc77(51)+acc77(68)
      acc77(68)=Qspk2*acc77(68)
      acc77(69)=acc77(31)*Qspvak1l3
      acc77(70)=acc77(36)*Qspvae6k2
      acc77(71)=acc77(43)*Qspval3k2
      acc77(72)=acc77(66)*Qspvak1e6
      acc77(73)=Qspvak6k2*acc77(63)
      acc77(74)=Qspvak1k6*acc77(55)
      acc77(69)=acc77(74)+acc77(73)+acc77(72)+acc77(64)+acc77(71)+acc77(70)+acc&
      &77(69)
      acc77(69)=Qspvak5k4*acc77(69)
      acc77(70)=acc77(42)*Qspvae6k4
      acc77(71)=-acc77(45)*Qspvak5e6
      acc77(72)=-acc77(46)*Qspvak1e6
      acc77(73)=acc77(52)*Qspvae6k2
      acc77(70)=acc77(73)+acc77(72)+acc77(71)+acc77(70)
      acc77(71)=Qspk1-Qspl3
      acc77(70)=acc77(71)*acc77(70)
      acc77(71)=acc77(9)*Qspvak1k4
      acc77(72)=acc77(23)*Qspvak5k2
      acc77(73)=acc77(40)*Qspk1
      acc77(71)=acc77(73)+acc77(72)+acc77(16)+acc77(71)
      acc77(71)=Qspe6*acc77(71)
      acc77(72)=acc77(25)*Qspvak5k2
      acc77(73)=-acc77(50)*Qspvak1k4
      acc77(72)=acc77(73)+acc77(72)
      acc77(73)=Qspk6+Qspl3
      acc77(72)=acc77(73)*acc77(72)
      acc77(73)=acc77(2)*Qspvak1e6
      acc77(74)=acc77(3)*Qspk6
      acc77(75)=acc77(4)*Qspvak1k4
      acc77(76)=acc77(8)*Qspvae6k2
      acc77(77)=acc77(11)*Qspvak6k4
      acc77(78)=acc77(22)*Qspvak5k2
      acc77(79)=acc77(26)*Qspvak1l3
      acc77(80)=acc77(27)*Qspvak5e6
      acc77(81)=acc77(29)*Qspval3k4
      acc77(82)=acc77(30)*Qspval3k2
      acc77(83)=acc77(34)*Qspvak5l3
      acc77(84)=acc77(37)*Qspl3
      acc77(85)=acc77(38)*Qspvak5k6
      acc77(86)=acc77(47)*Qspk1
      acc77(87)=acc77(48)*Qspvae6k6
      acc77(88)=acc77(49)*Qspvae6l3
      acc77(89)=acc77(56)*Qspvak6e6
      acc77(90)=acc77(59)*Qspvae6k4
      acc77(91)=acc77(61)*Qspval3e6
      acc77(92)=QspQ*acc77(44)
      brack=acc77(7)+acc77(67)+acc77(68)+acc77(69)+acc77(70)+acc77(71)+acc77(72&
      &)+acc77(73)+acc77(74)+acc77(75)+acc77(76)+acc77(77)+acc77(78)+acc77(79)+&
      &acc77(80)+acc77(81)+acc77(82)+acc77(83)+acc77(84)+acc77(85)+acc77(86)+ac&
      &c77(87)+acc77(88)+acc77(89)+acc77(90)+acc77(91)+acc77(92)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p0_dbard_hepemg_groups, only: &
!           & sign => diagram77_sign, shift => diagram77_shift
      use p0_dbard_hepemg_globalsl1, only: epspow
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_abbrevd77h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d77
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(-Q_ext(4),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d77 = 0.0_ki
      d77 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d77, ki), aimag(d77), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p0_dbard_hepemg_globalsl1, only: epspow
      use p0_dbard_hepemg_kinematics
      use p0_dbard_hepemg_abbrevd77h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d77
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(-Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d77 = 0.0_ki
      d77 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d77, ki), aimag(d77), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p0_dbard_hepemg_d77h0l1
