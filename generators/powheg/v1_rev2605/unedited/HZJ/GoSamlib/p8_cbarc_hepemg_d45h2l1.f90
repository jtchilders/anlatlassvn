module     p8_cbarc_hepemg_d45h2l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p8_cbarc_hepemg/helicity2d45h2l1.f90
   ! generator: buildfortran.py
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd45h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc45(70)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak1k5
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k5
      complex(ki) :: Qspvak4e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvak4k2
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspk1
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspvak6k5
      complex(ki) :: Qspvak4k6
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak4k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspl3
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak1k5 = dotproduct(Q,spvak1k5)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k5 = dotproduct(Q,spvae6k5)
      Qspvak4e6 = dotproduct(Q,spvak4e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvak4k2 = dotproduct(Q,spvak4k2)
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspk1 = dotproduct(Q,k1)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspvak6k5 = dotproduct(Q,spvak6k5)
      Qspvak4k6 = dotproduct(Q,spvak4k6)
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspk6 = dotproduct(Q,k6)
      Qspvak4k1 = dotproduct(Q,spvak4k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspl3 = dotproduct(Q,l3)
      acc45(1)=abb45(6)
      acc45(2)=abb45(7)
      acc45(3)=abb45(8)
      acc45(4)=abb45(9)
      acc45(5)=abb45(10)
      acc45(6)=abb45(11)
      acc45(7)=abb45(12)
      acc45(8)=abb45(13)
      acc45(9)=abb45(14)
      acc45(10)=abb45(15)
      acc45(11)=abb45(16)
      acc45(12)=abb45(17)
      acc45(13)=abb45(18)
      acc45(14)=abb45(19)
      acc45(15)=abb45(20)
      acc45(16)=abb45(23)
      acc45(17)=abb45(24)
      acc45(18)=abb45(25)
      acc45(19)=abb45(26)
      acc45(20)=abb45(28)
      acc45(21)=abb45(30)
      acc45(22)=abb45(32)
      acc45(23)=abb45(34)
      acc45(24)=abb45(35)
      acc45(25)=abb45(37)
      acc45(26)=abb45(38)
      acc45(27)=abb45(39)
      acc45(28)=abb45(40)
      acc45(29)=abb45(41)
      acc45(30)=abb45(43)
      acc45(31)=abb45(44)
      acc45(32)=abb45(48)
      acc45(33)=abb45(50)
      acc45(34)=abb45(51)
      acc45(35)=abb45(52)
      acc45(36)=abb45(55)
      acc45(37)=abb45(59)
      acc45(38)=abb45(63)
      acc45(39)=abb45(64)
      acc45(40)=abb45(65)
      acc45(41)=abb45(66)
      acc45(42)=abb45(67)
      acc45(43)=abb45(69)
      acc45(44)=abb45(71)
      acc45(45)=abb45(72)
      acc45(46)=abb45(74)
      acc45(47)=abb45(76)
      acc45(48)=abb45(77)
      acc45(49)=abb45(78)
      acc45(50)=abb45(79)
      acc45(51)=acc45(6)*Qspvak1k2
      acc45(52)=acc45(8)*Qspvak1k5
      acc45(53)=acc45(27)*Qspe6
      acc45(54)=acc45(37)*Qspvae6k5
      acc45(55)=acc45(38)*Qspvak4e6
      acc45(56)=acc45(39)*Qspvae6k2
      acc45(57)=acc45(41)*Qspvak1e6
      acc45(58)=acc45(48)*Qspvak4k2
      acc45(51)=acc45(58)+acc45(57)+acc45(56)+acc45(55)+acc45(54)+acc45(53)+acc&
      &45(22)+acc45(52)+acc45(51)
      acc45(51)=QspQ*acc45(51)
      acc45(52)=Qspk2-Qspk1
      acc45(52)=acc45(24)*acc45(52)
      acc45(53)=acc45(3)*Qspvak6k2
      acc45(54)=acc45(4)*Qspvak1k2
      acc45(55)=acc45(20)*Qspvak1k5
      acc45(56)=acc45(31)*Qspvak4k2
      acc45(57)=acc45(50)*Qspvak1k6
      acc45(52)=acc45(52)+acc45(57)+acc45(56)+acc45(55)+acc45(54)+acc45(53)+acc&
      &45(1)
      acc45(52)=Qspe6*acc45(52)
      acc45(53)=-acc45(29)*Qspvak1k5
      acc45(54)=acc45(30)*Qspvae6k5
      acc45(55)=acc45(33)*Qspvak4e6
      acc45(56)=-acc45(49)*Qspvak4k2
      acc45(53)=acc45(56)+acc45(55)+acc45(54)+acc45(53)
      acc45(54)=Qspk1+Qspk2
      acc45(53)=acc45(54)*acc45(53)
      acc45(55)=acc45(10)*Qspvae6k5
      acc45(56)=acc45(11)*Qspvak4e6
      acc45(57)=acc45(12)*Qspvak6k5
      acc45(58)=acc45(15)*Qspvak4k6
      acc45(55)=acc45(58)+acc45(57)+acc45(56)+acc45(55)+acc45(2)
      acc45(55)=Qspvak1k2*acc45(55)
      acc45(56)=acc45(13)*Qspvak6k2
      acc45(57)=acc45(19)*Qspvak1k6
      acc45(58)=acc45(45)*Qspvae6k2
      acc45(59)=acc45(46)*Qspvak1e6
      acc45(56)=acc45(59)+acc45(58)+acc45(57)+acc45(14)+acc45(56)
      acc45(56)=Qspvak4k5*acc45(56)
      acc45(57)=acc45(9)*Qspvak1k2
      acc45(58)=acc45(40)*Qspvae6k2
      acc45(59)=-acc45(42)*Qspvak1e6
      acc45(57)=acc45(59)+acc45(57)+acc45(58)
      acc45(54)=acc45(54)-Qspk6
      acc45(54)=acc45(54)*acc45(57)
      acc45(57)=acc45(25)*Qspvak1k5
      acc45(58)=-acc45(36)*Qspvak4k2
      acc45(57)=acc45(58)+acc45(32)+acc45(57)
      acc45(57)=Qspk6*acc45(57)
      acc45(58)=acc45(7)*Qspvak1e6
      acc45(59)=acc45(17)*Qspvak1k5
      acc45(60)=acc45(18)*Qspvae6k5
      acc45(61)=acc45(21)*Qspk1
      acc45(62)=acc45(23)*Qspk2
      acc45(63)=acc45(28)*Qspvak4e6
      acc45(64)=acc45(35)*Qspvae6k2
      acc45(65)=acc45(43)*Qspvak6k5
      acc45(66)=acc45(44)*Qspvak4k6
      acc45(67)=acc45(47)*Qspvak4k2
      acc45(68)=Qspvak4k1*acc45(26)
      acc45(69)=Qspvak2k5*acc45(16)
      acc45(70)=Qspl3*acc45(34)
      brack=acc45(5)+acc45(51)+acc45(52)+acc45(53)+acc45(54)+acc45(55)+acc45(56&
      &)+acc45(57)+acc45(58)+acc45(59)+acc45(60)+acc45(61)+acc45(62)+acc45(63)+&
      &acc45(64)+acc45(65)+acc45(66)+acc45(67)+acc45(68)+acc45(69)+acc45(70)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p8_cbarc_hepemg_groups, only: &
!           & sign => diagram45_sign, shift => diagram45_shift
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd45h2
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d45
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k3
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d45 = 0.0_ki
      d45 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d45, ki), aimag(d45), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd45h2
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d45
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k3
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d45 = 0.0_ki
      d45 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d45, ki), aimag(d45), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p8_cbarc_hepemg_d45h2l1
