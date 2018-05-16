module     p4_ubaru_hepemg_d277h3l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity3d277h3l1.f90
   ! generator: buildfortran.py
   use p4_ubaru_hepemg_config, only: ki
   use p4_ubaru_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p4_ubaru_hepemg_model
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_color
      use p4_ubaru_hepemg_abbrevd277h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc277(63)
      complex(ki) :: Qspvak4e6
      complex(ki) :: Qspe6
      complex(ki) :: Qspk6
      complex(ki) :: QspQ
      complex(ki) :: Qspvae6k5
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak6k5
      complex(ki) :: Qspvak4k6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak4k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspk1
      complex(ki) :: Qspk2
      Qspvak4e6 = dotproduct(Q,spvak4e6)
      Qspe6 = dotproduct(Q,e6)
      Qspk6 = dotproduct(Q,k6)
      QspQ = dotproduct(Q,Q)
      Qspvae6k5 = dotproduct(Q,spvae6k5)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak6k5 = dotproduct(Q,spvak6k5)
      Qspvak4k6 = dotproduct(Q,spvak4k6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak4k1 = dotproduct(Q,spvak4k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspk1 = dotproduct(Q,k1)
      Qspk2 = dotproduct(Q,k2)
      acc277(1)=abb277(16)
      acc277(2)=abb277(17)
      acc277(3)=abb277(18)
      acc277(4)=abb277(19)
      acc277(5)=abb277(20)
      acc277(6)=abb277(21)
      acc277(7)=abb277(22)
      acc277(8)=abb277(23)
      acc277(9)=abb277(24)
      acc277(10)=abb277(25)
      acc277(11)=abb277(26)
      acc277(12)=abb277(27)
      acc277(13)=abb277(28)
      acc277(14)=abb277(29)
      acc277(15)=abb277(30)
      acc277(16)=abb277(31)
      acc277(17)=abb277(32)
      acc277(18)=abb277(33)
      acc277(19)=abb277(34)
      acc277(20)=abb277(35)
      acc277(21)=abb277(36)
      acc277(22)=abb277(37)
      acc277(23)=abb277(38)
      acc277(24)=abb277(39)
      acc277(25)=abb277(40)
      acc277(26)=abb277(41)
      acc277(27)=abb277(42)
      acc277(28)=abb277(43)
      acc277(29)=abb277(44)
      acc277(30)=abb277(45)
      acc277(31)=abb277(46)
      acc277(32)=abb277(47)
      acc277(33)=abb277(48)
      acc277(34)=abb277(49)
      acc277(35)=abb277(50)
      acc277(36)=abb277(51)
      acc277(37)=abb277(53)
      acc277(38)=abb277(55)
      acc277(39)=abb277(56)
      acc277(40)=abb277(59)
      acc277(41)=abb277(62)
      acc277(42)=abb277(63)
      acc277(43)=abb277(65)
      acc277(44)=abb277(66)
      acc277(45)=abb277(67)
      acc277(46)=abb277(68)
      acc277(47)=acc277(12)*Qspvak4e6
      acc277(48)=acc277(16)*Qspe6
      acc277(49)=acc277(17)*Qspk6
      acc277(50)=acc277(18)*QspQ
      acc277(51)=acc277(20)*Qspvae6k5
      acc277(52)=acc277(23)*Qspvak2e6
      acc277(53)=acc277(25)*Qspvae6k1
      acc277(54)=acc277(32)*Qspvak6k5
      acc277(55)=acc277(36)*Qspvak4k6
      acc277(56)=Qspvae6k2*acc277(4)
      acc277(57)=Qspvak1e6*acc277(6)
      acc277(47)=acc277(57)+acc277(56)+acc277(55)+acc277(54)+acc277(53)+acc277(&
      &52)+acc277(51)+acc277(50)+acc277(49)+acc277(48)+acc277(47)+acc277(10)
      acc277(47)=Qspvak2k1*acc277(47)
      acc277(48)=acc277(26)*QspQ
      acc277(49)=acc277(33)*Qspvak2k6
      acc277(50)=acc277(42)*Qspvak6k1
      acc277(51)=acc277(43)*Qspvak4k1
      acc277(52)=acc277(45)*Qspvak2k5
      acc277(48)=acc277(52)+acc277(51)+acc277(50)+acc277(49)+acc277(48)+acc277(&
      &7)
      acc277(48)=Qspe6*acc277(48)
      acc277(49)=acc277(1)*Qspvak2k5
      acc277(50)=acc277(2)*Qspvae6k5
      acc277(51)=acc277(11)*Qspvak4e6
      acc277(52)=acc277(13)*Qspvak4k1
      acc277(49)=acc277(52)+acc277(51)+acc277(50)+acc277(49)
      acc277(49)=QspQ*acc277(49)
      acc277(50)=acc277(9)*Qspvak2e6
      acc277(51)=acc277(21)*Qspvak2k6
      acc277(52)=acc277(24)*Qspvak6k1
      acc277(53)=acc277(40)*Qspvae6k1
      acc277(50)=acc277(53)+acc277(30)+acc277(52)+acc277(51)+acc277(50)
      acc277(50)=Qspvak4k5*acc277(50)
      acc277(51)=-acc277(35)*Qspe6
      acc277(52)=acc277(37)*Qspvae6k5
      acc277(53)=-acc277(39)*Qspvak4e6
      acc277(51)=acc277(51)+acc277(53)+acc277(52)
      acc277(52)=Qspk1+Qspk2
      acc277(51)=acc277(52)*acc277(51)
      acc277(53)=acc277(14)*Qspvae6k1
      acc277(54)=-acc277(19)*Qspvak2e6
      acc277(53)=acc277(53)+acc277(54)
      acc277(52)=acc277(52)+Qspk6+QspQ
      acc277(52)=acc277(52)*acc277(53)
      acc277(53)=acc277(44)*Qspvak4k1
      acc277(54)=acc277(46)*Qspvak2k5
      acc277(53)=acc277(54)+acc277(53)+acc277(15)
      acc277(53)=Qspk6*acc277(53)
      acc277(54)=acc277(5)*Qspvak4k1
      acc277(55)=acc277(8)*Qspvak2k5
      acc277(56)=acc277(22)*Qspk1
      acc277(57)=acc277(27)*Qspvae6k5
      acc277(58)=acc277(28)*Qspvak2e6
      acc277(59)=acc277(29)*Qspvak6k5
      acc277(60)=acc277(31)*Qspvak4k6
      acc277(61)=acc277(34)*Qspk2
      acc277(62)=acc277(38)*Qspvak4e6
      acc277(63)=acc277(41)*Qspvae6k1
      brack=acc277(3)+acc277(47)+acc277(48)+acc277(49)+acc277(50)+acc277(51)+ac&
      &c277(52)+acc277(53)+acc277(54)+acc277(55)+acc277(56)+acc277(57)+acc277(5&
      &8)+acc277(59)+acc277(60)+acc277(61)+acc277(62)+acc277(63)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p4_ubaru_hepemg_groups, only: &
!           & sign => diagram277_sign, shift => diagram277_shift
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd277h3
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d277
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k6+k5+k4
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d277 = 0.0_ki
      d277 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d277, ki), aimag(d277), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd277h3
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d277
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k6+k5+k4
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d277 = 0.0_ki
      d277 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d277, ki), aimag(d277), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p4_ubaru_hepemg_d277h3l1
