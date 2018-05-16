module     p1_dbarc_hepneg_d59h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p1_dbarc_hepneg/helicity0d59h0l1.f90
   ! generator: buildfortran.py
   use p1_dbarc_hepneg_config, only: ki
   use p1_dbarc_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p1_dbarc_hepneg_model
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_color
      use p1_dbarc_hepneg_abbrevd59h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc59(35)
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspvak6k2
      complex(ki) :: QspQ
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspk1
      complex(ki) :: Qspe6
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspk6
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      QspQ = dotproduct(Q,Q)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspk1 = dotproduct(Q,k1)
      Qspe6 = dotproduct(Q,e6)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspk6 = dotproduct(Q,k6)
      acc59(1)=abb59(4)
      acc59(2)=abb59(5)
      acc59(3)=abb59(6)
      acc59(4)=abb59(7)
      acc59(5)=abb59(8)
      acc59(6)=abb59(9)
      acc59(7)=abb59(10)
      acc59(8)=abb59(11)
      acc59(9)=abb59(12)
      acc59(10)=abb59(13)
      acc59(11)=abb59(14)
      acc59(12)=abb59(15)
      acc59(13)=abb59(16)
      acc59(14)=abb59(17)
      acc59(15)=abb59(18)
      acc59(16)=abb59(19)
      acc59(17)=abb59(20)
      acc59(18)=abb59(21)
      acc59(19)=abb59(22)
      acc59(20)=abb59(23)
      acc59(21)=abb59(24)
      acc59(22)=abb59(25)
      acc59(23)=abb59(26)
      acc59(24)=abb59(28)
      acc59(25)=abb59(29)
      acc59(26)=abb59(30)
      acc59(27)=abb59(31)
      acc59(28)=Qspvak5k4*acc59(8)
      acc59(28)=acc59(28)+acc59(5)
      acc59(28)=Qspvak1k2*acc59(28)
      acc59(29)=acc59(9)*Qspvak1k6
      acc59(30)=Qspvak6k2*acc59(21)
      acc59(31)=QspQ*acc59(15)
      acc59(32)=Qspvak5k2*acc59(10)
      acc59(33)=Qspk1*acc59(3)
      acc59(28)=acc59(33)+acc59(32)+acc59(31)+acc59(30)+acc59(29)+acc59(7)+acc5&
      &9(28)
      acc59(28)=Qspe6*acc59(28)
      acc59(29)=-Qspvak6k2*acc59(27)
      acc59(30)=Qspvak1k2*acc59(14)
      acc59(31)=Qspvak1e6*acc59(26)
      acc59(32)=Qspvae6k2*acc59(25)
      acc59(33)=Qspk1*acc59(24)
      acc59(29)=acc59(33)+acc59(32)+acc59(31)+acc59(30)+acc59(23)+acc59(29)
      acc59(29)=Qspvak5k4*acc59(29)
      acc59(30)=Qspvak1e6*acc59(18)
      acc59(31)=Qspvae6k2*acc59(16)
      acc59(30)=acc59(30)-acc59(31)
      acc59(31)=acc59(13)-acc59(30)
      acc59(31)=Qspk6*acc59(31)
      acc59(32)=Qspvak1e6*acc59(11)
      acc59(33)=Qspvae6k2*acc59(17)
      acc59(32)=acc59(33)+acc59(12)+acc59(32)
      acc59(32)=QspQ*acc59(32)
      acc59(33)=-Qspk6*acc59(19)
      acc59(34)=QspQ*acc59(22)
      acc59(33)=acc59(34)+acc59(4)+acc59(33)
      acc59(33)=Qspvak5k2*acc59(33)
      acc59(34)=Qspvak5k2*acc59(19)
      acc59(30)=acc59(34)+acc59(1)+acc59(30)
      acc59(30)=Qspk1*acc59(30)
      acc59(34)=Qspvak6k2*acc59(20)
      acc59(35)=Qspvak1k2*acc59(2)
      brack=acc59(6)+acc59(28)+acc59(29)+acc59(30)+acc59(31)+acc59(32)+acc59(33&
      &)+acc59(34)+acc59(35)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p1_dbarc_hepneg_groups, only: &
!           & sign => diagram59_sign, shift => diagram59_shift
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd59h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d59
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d59 = 0.0_ki
      d59 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d59, ki), aimag(d59), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd59h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d59
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d59 = 0.0_ki
      d59 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d59, ki), aimag(d59), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p1_dbarc_hepneg_d59h0l1
