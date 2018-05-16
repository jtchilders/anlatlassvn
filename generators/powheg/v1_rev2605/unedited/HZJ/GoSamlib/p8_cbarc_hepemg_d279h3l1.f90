module     p8_cbarc_hepemg_d279h3l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p8_cbarc_hepemg/helicity3d279h3l1.f90
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
      use p8_cbarc_hepemg_abbrevd279h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc279(58)
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak6k5
      complex(ki) :: Qspvae6k5
      complex(ki) :: Qspe6
      complex(ki) :: Qspvak4k6
      complex(ki) :: Qspk6
      complex(ki) :: QspQ
      complex(ki) :: Qspvak4e6
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak4k1
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspk1
      complex(ki) :: Qspk2
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak6k5 = dotproduct(Q,spvak6k5)
      Qspvae6k5 = dotproduct(Q,spvae6k5)
      Qspe6 = dotproduct(Q,e6)
      Qspvak4k6 = dotproduct(Q,spvak4k6)
      Qspk6 = dotproduct(Q,k6)
      QspQ = dotproduct(Q,Q)
      Qspvak4e6 = dotproduct(Q,spvak4e6)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak4k1 = dotproduct(Q,spvak4k1)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspk1 = dotproduct(Q,k1)
      Qspk2 = dotproduct(Q,k2)
      acc279(1)=abb279(12)
      acc279(2)=abb279(14)
      acc279(3)=abb279(15)
      acc279(4)=abb279(16)
      acc279(5)=abb279(17)
      acc279(6)=abb279(18)
      acc279(7)=abb279(19)
      acc279(8)=abb279(20)
      acc279(9)=abb279(21)
      acc279(10)=abb279(22)
      acc279(11)=abb279(23)
      acc279(12)=abb279(24)
      acc279(13)=abb279(25)
      acc279(14)=abb279(26)
      acc279(15)=abb279(27)
      acc279(16)=abb279(28)
      acc279(17)=abb279(29)
      acc279(18)=abb279(30)
      acc279(19)=abb279(31)
      acc279(20)=abb279(32)
      acc279(21)=abb279(33)
      acc279(22)=abb279(34)
      acc279(23)=abb279(35)
      acc279(24)=abb279(36)
      acc279(25)=abb279(37)
      acc279(26)=abb279(39)
      acc279(27)=abb279(40)
      acc279(28)=abb279(42)
      acc279(29)=abb279(43)
      acc279(30)=abb279(44)
      acc279(31)=abb279(45)
      acc279(32)=abb279(46)
      acc279(33)=abb279(48)
      acc279(34)=abb279(50)
      acc279(35)=abb279(51)
      acc279(36)=abb279(52)
      acc279(37)=abb279(53)
      acc279(38)=abb279(54)
      acc279(39)=abb279(55)
      acc279(40)=abb279(56)
      acc279(41)=abb279(57)
      acc279(42)=abb279(58)
      acc279(43)=abb279(60)
      acc279(44)=acc279(3)*Qspvae6k1
      acc279(45)=acc279(5)*Qspvak6k5
      acc279(46)=acc279(7)*Qspvae6k5
      acc279(47)=acc279(15)*Qspe6
      acc279(48)=acc279(17)*Qspvak4k6
      acc279(49)=acc279(20)*Qspk6
      acc279(50)=acc279(21)*QspQ
      acc279(51)=acc279(23)*Qspvak4e6
      acc279(52)=acc279(26)*Qspvak2e6
      acc279(53)=Qspvae6k2*acc279(8)
      acc279(54)=Qspvak1e6*acc279(4)
      acc279(44)=acc279(54)+acc279(53)+acc279(52)+acc279(51)+acc279(50)+acc279(&
      &49)+acc279(48)+acc279(47)+acc279(46)+acc279(45)+acc279(44)+acc279(1)
      acc279(44)=Qspvak2k1*acc279(44)
      acc279(45)=acc279(9)*Qspvak2k5
      acc279(46)=acc279(11)*Qspvak4k1
      acc279(47)=acc279(24)*QspQ
      acc279(48)=acc279(34)*Qspvak6k1
      acc279(49)=acc279(43)*Qspvak2k6
      acc279(45)=acc279(49)+acc279(48)+acc279(47)+acc279(46)+acc279(10)+acc279(&
      &45)
      acc279(45)=Qspe6*acc279(45)
      acc279(46)=acc279(2)*Qspvak2k5
      acc279(47)=acc279(12)*Qspvak4e6
      acc279(48)=acc279(30)*Qspvae6k5
      acc279(49)=acc279(42)*Qspvak4k1
      acc279(46)=acc279(49)+acc279(48)+acc279(47)+acc279(46)
      acc279(46)=QspQ*acc279(46)
      acc279(47)=acc279(25)*Qspvae6k1
      acc279(48)=acc279(35)*Qspvak2k6
      acc279(49)=acc279(37)*Qspvak2e6
      acc279(50)=acc279(40)*Qspvak6k1
      acc279(47)=acc279(50)+acc279(38)+acc279(49)+acc279(48)+acc279(47)
      acc279(47)=Qspvak4k5*acc279(47)
      acc279(48)=acc279(6)*Qspvae6k5
      acc279(49)=-acc279(27)*Qspe6
      acc279(50)=-acc279(31)*Qspvak4e6
      acc279(48)=acc279(49)+acc279(50)+acc279(48)
      acc279(49)=Qspk1+Qspk2
      acc279(48)=acc279(49)*acc279(48)
      acc279(50)=acc279(18)*Qspvae6k1
      acc279(51)=-acc279(22)*Qspvak2e6
      acc279(50)=acc279(50)+acc279(51)
      acc279(49)=acc279(49)+Qspk6+QspQ
      acc279(49)=acc279(49)*acc279(50)
      acc279(50)=acc279(39)*Qspvak2k5
      acc279(51)=acc279(41)*Qspvak4k1
      acc279(50)=acc279(51)+acc279(50)+acc279(19)
      acc279(50)=Qspk6*acc279(50)
      acc279(51)=acc279(13)*Qspvak2e6
      acc279(52)=acc279(14)*Qspvak2k5
      acc279(53)=acc279(16)*Qspk2
      acc279(54)=acc279(28)*Qspk1
      acc279(55)=acc279(29)*Qspvak6k5
      acc279(56)=acc279(32)*Qspvak4k1
      acc279(57)=acc279(33)*Qspvae6k1
      acc279(58)=acc279(36)*Qspvak4k6
      brack=acc279(44)+acc279(45)+acc279(46)+acc279(47)+acc279(48)+acc279(49)+a&
      &cc279(50)+acc279(51)+acc279(52)+acc279(53)+acc279(54)+acc279(55)+acc279(&
      &56)+acc279(57)+acc279(58)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p8_cbarc_hepemg_groups, only: &
!           & sign => diagram279_sign, shift => diagram279_shift
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd279h3
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d279
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(-Q_ext(4),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d279 = 0.0_ki
      d279 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d279, ki), aimag(d279), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd279h3
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d279
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(-Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d279 = 0.0_ki
      d279 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d279, ki), aimag(d279), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p8_cbarc_hepemg_d279h3l1
