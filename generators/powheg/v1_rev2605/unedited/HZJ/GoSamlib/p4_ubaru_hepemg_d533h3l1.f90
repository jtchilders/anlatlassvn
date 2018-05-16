module     p4_ubaru_hepemg_d533h3l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity3d533h3l1.f90
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
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc533(35)
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak2k6
      complex(ki) :: QspQ
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspk1
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspk6
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      QspQ = dotproduct(Q,Q)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspk1 = dotproduct(Q,k1)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspk6 = dotproduct(Q,k6)
      acc533(1)=abb533(4)
      acc533(2)=abb533(5)
      acc533(3)=abb533(6)
      acc533(4)=abb533(7)
      acc533(5)=abb533(8)
      acc533(6)=abb533(9)
      acc533(7)=abb533(10)
      acc533(8)=abb533(11)
      acc533(9)=abb533(12)
      acc533(10)=abb533(13)
      acc533(11)=abb533(14)
      acc533(12)=abb533(15)
      acc533(13)=abb533(16)
      acc533(14)=abb533(17)
      acc533(15)=abb533(18)
      acc533(16)=abb533(19)
      acc533(17)=abb533(20)
      acc533(18)=abb533(21)
      acc533(19)=abb533(22)
      acc533(20)=abb533(23)
      acc533(21)=abb533(25)
      acc533(22)=abb533(26)
      acc533(23)=abb533(27)
      acc533(24)=abb533(28)
      acc533(25)=abb533(29)
      acc533(26)=abb533(30)
      acc533(27)=abb533(31)
      acc533(28)=Qspvak4k5*acc533(8)
      acc533(28)=acc533(28)+acc533(5)
      acc533(28)=Qspvak2k1*acc533(28)
      acc533(29)=acc533(9)*Qspvak6k1
      acc533(30)=Qspvak2k6*acc533(26)
      acc533(31)=QspQ*acc533(15)
      acc533(32)=Qspvak2k5*acc533(10)
      acc533(33)=Qspk1*acc533(3)
      acc533(28)=acc533(33)+acc533(32)+acc533(31)+acc533(30)+acc533(29)+acc533(&
      &7)+acc533(28)
      acc533(28)=Qspe6*acc533(28)
      acc533(29)=-Qspvak2k6*acc533(27)
      acc533(30)=Qspvak2k1*acc533(14)
      acc533(31)=Qspvae6k1*acc533(23)
      acc533(32)=Qspvak2e6*acc533(22)
      acc533(33)=Qspk1*acc533(21)
      acc533(29)=acc533(33)+acc533(32)+acc533(31)+acc533(30)+acc533(20)+acc533(&
      &29)
      acc533(29)=Qspvak4k5*acc533(29)
      acc533(30)=Qspvae6k1*acc533(18)
      acc533(31)=Qspvak2e6*acc533(16)
      acc533(30)=acc533(30)-acc533(31)
      acc533(31)=acc533(13)-acc533(30)
      acc533(31)=Qspk6*acc533(31)
      acc533(32)=Qspvae6k1*acc533(11)
      acc533(33)=Qspvak2e6*acc533(17)
      acc533(32)=acc533(33)+acc533(12)+acc533(32)
      acc533(32)=QspQ*acc533(32)
      acc533(33)=-Qspk6*acc533(19)
      acc533(34)=QspQ*acc533(24)
      acc533(33)=acc533(34)+acc533(4)+acc533(33)
      acc533(33)=Qspvak2k5*acc533(33)
      acc533(34)=Qspvak2k5*acc533(19)
      acc533(30)=acc533(34)+acc533(1)+acc533(30)
      acc533(30)=Qspk1*acc533(30)
      acc533(34)=Qspvak2k6*acc533(25)
      acc533(35)=Qspvak2k1*acc533(2)
      brack=acc533(6)+acc533(28)+acc533(29)+acc533(30)+acc533(31)+acc533(32)+ac&
      &c533(33)+acc533(34)+acc533(35)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p4_ubaru_hepemg_groups, only: &
!           & sign => diagram533_sign, shift => diagram533_shift
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d533
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d533 = 0.0_ki
      d533 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d533, ki), aimag(d533), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd533h3
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d533
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d533 = 0.0_ki
      d533 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d533, ki), aimag(d533), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p4_ubaru_hepemg_d533h3l1
