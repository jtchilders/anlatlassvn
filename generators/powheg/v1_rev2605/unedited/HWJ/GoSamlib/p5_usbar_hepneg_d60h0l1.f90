module     p5_usbar_hepneg_d60h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p5_usbar_hepneg/helicity0d60h0l1.f90
   ! generator: buildfortran.py
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc60(32)
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak6k1
      complex(ki) :: QspQ
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspk2
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspk6
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      QspQ = dotproduct(Q,Q)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspk2 = dotproduct(Q,k2)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspk6 = dotproduct(Q,k6)
      acc60(1)=abb60(4)
      acc60(2)=abb60(5)
      acc60(3)=abb60(6)
      acc60(4)=abb60(7)
      acc60(5)=abb60(8)
      acc60(6)=abb60(9)
      acc60(7)=abb60(10)
      acc60(8)=abb60(11)
      acc60(9)=abb60(12)
      acc60(10)=abb60(13)
      acc60(11)=abb60(14)
      acc60(12)=abb60(15)
      acc60(13)=abb60(16)
      acc60(14)=abb60(17)
      acc60(15)=abb60(18)
      acc60(16)=abb60(19)
      acc60(17)=abb60(20)
      acc60(18)=abb60(21)
      acc60(19)=abb60(22)
      acc60(20)=abb60(23)
      acc60(21)=abb60(24)
      acc60(22)=abb60(25)
      acc60(23)=abb60(27)
      acc60(24)=abb60(28)
      acc60(25)=abb60(30)
      acc60(26)=-Qspvak5k4*acc60(8)
      acc60(26)=acc60(26)+acc60(14)
      acc60(26)=Qspvak2k1*acc60(26)
      acc60(27)=acc60(9)*Qspvak2k6
      acc60(28)=Qspvak6k1*acc60(20)
      acc60(29)=QspQ*acc60(15)
      acc60(30)=Qspvak5k1*acc60(10)
      acc60(31)=Qspk2*acc60(2)
      acc60(26)=acc60(31)+acc60(30)+acc60(29)+acc60(28)+acc60(27)+acc60(7)+acc6&
      &0(26)
      acc60(26)=Qspe6*acc60(26)
      acc60(27)=Qspvak6k1*acc60(5)
      acc60(28)=Qspvae6k1*acc60(25)
      acc60(29)=Qspvak2e6*acc60(24)
      acc60(30)=Qspk2*acc60(23)
      acc60(27)=acc60(30)+acc60(29)+acc60(28)+acc60(22)+acc60(27)
      acc60(27)=Qspvak5k4*acc60(27)
      acc60(28)=Qspvae6k1*acc60(17)
      acc60(29)=Qspvak2e6*acc60(16)
      acc60(28)=acc60(28)-acc60(29)
      acc60(29)=acc60(13)-acc60(28)
      acc60(29)=Qspk6*acc60(29)
      acc60(30)=Qspvae6k1*acc60(18)
      acc60(31)=Qspvak2e6*acc60(11)
      acc60(30)=acc60(31)+acc60(12)+acc60(30)
      acc60(30)=QspQ*acc60(30)
      acc60(31)=Qspk6*acc60(21)
      acc60(32)=QspQ*acc60(4)
      acc60(31)=acc60(32)+acc60(3)+acc60(31)
      acc60(31)=Qspvak5k1*acc60(31)
      acc60(32)=-Qspvak5k1*acc60(21)
      acc60(28)=acc60(32)+acc60(1)+acc60(28)
      acc60(28)=Qspk2*acc60(28)
      acc60(32)=Qspvak6k1*acc60(19)
      brack=acc60(6)+acc60(26)+acc60(27)+acc60(28)+acc60(29)+acc60(30)+acc60(31&
      &)+acc60(32)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p5_usbar_hepneg_groups, only: &
!           & sign => diagram60_sign, shift => diagram60_shift
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d60
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d60 = 0.0_ki
      d60 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d60, ki), aimag(d60), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd60h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d60
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k2
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d60 = 0.0_ki
      d60 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d60, ki), aimag(d60), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p5_usbar_hepneg_d60h0l1
