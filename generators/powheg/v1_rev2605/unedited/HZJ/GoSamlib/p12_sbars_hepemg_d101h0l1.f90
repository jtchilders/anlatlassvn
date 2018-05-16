module     p12_sbars_hepemg_d101h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity0d101h0l1.f90
   ! generator: buildfortran.py
   use p12_sbars_hepemg_config, only: ki
   use p12_sbars_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd101h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc101(65)
      complex(ki) :: Qspk1
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspvak1k4
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak1k2
      complex(ki) :: QspQ
      complex(ki) :: Qspe6
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvae6k4
      complex(ki) :: Qspvak5e6
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak5k6
      complex(ki) :: Qspvak6k4
      complex(ki) :: Qspl3
      Qspk1 = dotproduct(Q,k1)
      Qspk2 = dotproduct(Q,k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspvak1k4 = dotproduct(Q,spvak1k4)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspk6 = dotproduct(Q,k6)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      QspQ = dotproduct(Q,Q)
      Qspe6 = dotproduct(Q,e6)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvae6k4 = dotproduct(Q,spvae6k4)
      Qspvak5e6 = dotproduct(Q,spvak5e6)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak5k6 = dotproduct(Q,spvak5k6)
      Qspvak6k4 = dotproduct(Q,spvak6k4)
      Qspl3 = dotproduct(Q,l3)
      acc101(1)=abb101(6)
      acc101(2)=abb101(7)
      acc101(3)=abb101(8)
      acc101(4)=abb101(9)
      acc101(5)=abb101(10)
      acc101(6)=abb101(11)
      acc101(7)=abb101(12)
      acc101(8)=abb101(13)
      acc101(9)=abb101(14)
      acc101(10)=abb101(15)
      acc101(11)=abb101(16)
      acc101(12)=abb101(17)
      acc101(13)=abb101(18)
      acc101(14)=abb101(19)
      acc101(15)=abb101(20)
      acc101(16)=abb101(21)
      acc101(17)=abb101(22)
      acc101(18)=abb101(23)
      acc101(19)=abb101(24)
      acc101(20)=abb101(25)
      acc101(21)=abb101(26)
      acc101(22)=abb101(27)
      acc101(23)=abb101(28)
      acc101(24)=abb101(29)
      acc101(25)=abb101(30)
      acc101(26)=abb101(31)
      acc101(27)=abb101(32)
      acc101(28)=abb101(33)
      acc101(29)=abb101(35)
      acc101(30)=abb101(36)
      acc101(31)=abb101(37)
      acc101(32)=abb101(39)
      acc101(33)=abb101(40)
      acc101(34)=abb101(42)
      acc101(35)=abb101(45)
      acc101(36)=abb101(46)
      acc101(37)=abb101(47)
      acc101(38)=abb101(50)
      acc101(39)=abb101(53)
      acc101(40)=abb101(58)
      acc101(41)=abb101(60)
      acc101(42)=abb101(61)
      acc101(43)=abb101(62)
      acc101(44)=abb101(63)
      acc101(45)=abb101(65)
      acc101(46)=abb101(66)
      acc101(47)=abb101(67)
      acc101(48)=Qspk1*acc101(39)
      acc101(49)=Qspk2*acc101(40)
      acc101(50)=Qspvak1k6*acc101(33)
      acc101(51)=Qspvak6k2*acc101(46)
      acc101(52)=Qspvak1k4*acc101(2)
      acc101(53)=Qspvak5k2*acc101(30)
      acc101(54)=Qspk6*acc101(35)
      acc101(55)=Qspvak1k2*acc101(9)
      acc101(56)=QspQ*acc101(20)
      acc101(48)=acc101(56)+acc101(55)+acc101(54)+acc101(53)+acc101(52)+acc101(&
      &51)+acc101(50)+acc101(49)+acc101(18)+acc101(48)
      acc101(48)=Qspe6*acc101(48)
      acc101(49)=Qspvak1e6*acc101(43)
      acc101(50)=Qspvae6k2*acc101(42)
      acc101(51)=Qspvak1k4*acc101(17)
      acc101(52)=Qspvak5k2*acc101(36)
      acc101(53)=Qspvae6k4*acc101(28)
      acc101(54)=Qspvak5e6*acc101(29)
      acc101(55)=Qspvak1k2*acc101(11)
      acc101(49)=acc101(55)+acc101(54)+acc101(53)+acc101(52)+acc101(51)+acc101(&
      &50)+acc101(3)+acc101(49)
      acc101(49)=QspQ*acc101(49)
      acc101(50)=Qspvak1k6*acc101(27)
      acc101(51)=Qspvak6k2*acc101(45)
      acc101(52)=Qspvak1e6*acc101(31)
      acc101(53)=-Qspvae6k2*acc101(26)
      acc101(50)=acc101(53)+acc101(52)+acc101(51)+acc101(32)+acc101(50)
      acc101(50)=Qspvak5k4*acc101(50)
      acc101(51)=Qspvak1k4*acc101(6)
      acc101(52)=Qspvak5k2*acc101(1)
      acc101(53)=Qspvae6k4*acc101(37)
      acc101(54)=Qspvak5e6*acc101(34)
      acc101(51)=acc101(54)+acc101(53)+acc101(52)+acc101(8)+acc101(51)
      acc101(51)=Qspk6*acc101(51)
      acc101(52)=Qspvak5k6*acc101(24)
      acc101(53)=Qspvak6k4*acc101(22)
      acc101(54)=Qspvae6k4*acc101(25)
      acc101(55)=Qspvak5e6*acc101(13)
      acc101(52)=acc101(55)+acc101(54)+acc101(53)+acc101(4)+acc101(52)
      acc101(52)=Qspvak1k2*acc101(52)
      acc101(53)=acc101(38)*Qspl3
      acc101(54)=Qspk1*acc101(14)
      acc101(55)=Qspk2*acc101(12)
      acc101(56)=Qspvak5k6*acc101(47)
      acc101(57)=Qspvak6k4*acc101(44)
      acc101(58)=Qspvak1k6*acc101(5)
      acc101(59)=Qspvak6k2*acc101(10)
      acc101(60)=Qspvak1e6*acc101(21)
      acc101(61)=Qspvae6k2*acc101(41)
      acc101(62)=Qspvak1k4*acc101(16)
      acc101(63)=Qspvak5k2*acc101(19)
      acc101(64)=Qspvae6k4*acc101(23)
      acc101(65)=Qspvak5e6*acc101(7)
      brack=acc101(15)+acc101(48)+acc101(49)+acc101(50)+acc101(51)+acc101(52)+a&
      &cc101(53)+acc101(54)+acc101(55)+acc101(56)+acc101(57)+acc101(58)+acc101(&
      &59)+acc101(60)+acc101(61)+acc101(62)+acc101(63)+acc101(64)+acc101(65)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p12_sbars_hepemg_groups, only: &
!           & sign => diagram101_sign, shift => diagram101_shift
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd101h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d101
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(+Q_ext(4),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d101 = 0.0_ki
      d101 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d101, ki), aimag(d101), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd101h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d101
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(+Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d101 = 0.0_ki
      d101 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d101, ki), aimag(d101), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p12_sbars_hepemg_d101h0l1
