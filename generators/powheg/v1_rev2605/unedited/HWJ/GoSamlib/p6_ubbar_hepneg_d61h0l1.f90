module     p6_ubbar_hepneg_d61h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p6_ubbar_hepneg/helicity0d61h0l1.f90
   ! generator: buildfortran.py
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p6_ubbar_hepneg_model
      use p6_ubbar_hepneg_kinematics
      use p6_ubbar_hepneg_color
      use p6_ubbar_hepneg_abbrevd61h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc61(34)
      complex(ki) :: Qspk6
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspvak6k1
      complex(ki) :: QspQ
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspe6
      complex(ki) :: Qspk1
      complex(ki) :: Qspvae6k1
      Qspk6 = dotproduct(Q,k6)
      Qspk2 = dotproduct(Q,k2)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      QspQ = dotproduct(Q,Q)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspe6 = dotproduct(Q,e6)
      Qspk1 = dotproduct(Q,k1)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      acc61(1)=abb61(3)
      acc61(2)=abb61(4)
      acc61(3)=abb61(5)
      acc61(4)=abb61(6)
      acc61(5)=abb61(7)
      acc61(6)=abb61(8)
      acc61(7)=abb61(9)
      acc61(8)=abb61(10)
      acc61(9)=abb61(11)
      acc61(10)=abb61(12)
      acc61(11)=abb61(13)
      acc61(12)=abb61(14)
      acc61(13)=abb61(15)
      acc61(14)=abb61(16)
      acc61(15)=abb61(17)
      acc61(16)=abb61(18)
      acc61(17)=abb61(19)
      acc61(18)=abb61(20)
      acc61(19)=abb61(21)
      acc61(20)=abb61(23)
      acc61(21)=abb61(24)
      acc61(22)=abb61(25)
      acc61(23)=abb61(26)
      acc61(24)=abb61(27)
      acc61(25)=abb61(28)
      acc61(26)=abb61(29)
      acc61(27)=-Qspk6-Qspk2
      acc61(27)=acc61(13)*acc61(27)
      acc61(28)=acc61(9)*Qspvak2k6
      acc61(29)=Qspvak2k1*acc61(3)
      acc61(30)=Qspvak5k1*acc61(2)
      acc61(31)=Qspvak6k1*acc61(20)
      acc61(32)=QspQ*acc61(4)
      acc61(33)=-Qspvak2k1*acc61(23)
      acc61(33)=acc61(21)+acc61(33)
      acc61(33)=Qspvak5k4*acc61(33)
      acc61(27)=acc61(33)+acc61(32)+acc61(31)+acc61(30)+acc61(29)+acc61(12)+acc&
      &61(28)+acc61(27)
      acc61(27)=Qspe6*acc61(27)
      acc61(28)=Qspk1*acc61(25)
      acc61(29)=Qspvak6k1*acc61(26)
      acc61(30)=Qspk6*acc61(24)
      acc61(31)=-QspQ*acc61(22)
      acc61(32)=Qspvae6k1*acc61(15)
      acc61(28)=acc61(32)+acc61(31)+acc61(30)+acc61(29)+acc61(1)+acc61(28)
      acc61(28)=Qspvak5k4*acc61(28)
      acc61(29)=-Qspk6+Qspk2
      acc61(29)=acc61(19)*acc61(29)
      acc61(30)=-QspQ*acc61(10)
      acc61(29)=acc61(30)+acc61(6)+acc61(29)
      acc61(29)=Qspvae6k1*acc61(29)
      acc61(30)=Qspvak5k1*acc61(8)
      acc61(30)=acc61(30)+acc61(18)
      acc61(30)=Qspk1*acc61(30)
      acc61(31)=Qspk2*acc61(14)
      acc61(32)=Qspvak6k1*acc61(17)
      acc61(33)=Qspk6*acc61(11)
      acc61(34)=Qspvak5k1*acc61(16)
      acc61(34)=acc61(5)+acc61(34)
      acc61(34)=QspQ*acc61(34)
      brack=acc61(7)+acc61(27)+acc61(28)+acc61(29)+acc61(30)+acc61(31)+acc61(32&
      &)+acc61(33)+acc61(34)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p6_ubbar_hepneg_groups, only: &
!           & sign => diagram61_sign, shift => diagram61_shift
      use p6_ubbar_hepneg_globalsl1, only: epspow
      use p6_ubbar_hepneg_kinematics
      use p6_ubbar_hepneg_abbrevd61h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d61
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k6
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d61 = 0.0_ki
      d61 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d61, ki), aimag(d61), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p6_ubbar_hepneg_globalsl1, only: epspow
      use p6_ubbar_hepneg_kinematics
      use p6_ubbar_hepneg_abbrevd61h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d61
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k6
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d61 = 0.0_ki
      d61 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d61, ki), aimag(d61), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p6_ubbar_hepneg_d61h0l1
