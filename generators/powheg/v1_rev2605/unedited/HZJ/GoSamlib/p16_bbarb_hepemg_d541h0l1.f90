module     p16_bbarb_hepemg_d541h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity0d541h0l1.f90
   ! generator: buildfortran.py
   use p16_bbarb_hepemg_config, only: ki
   use p16_bbarb_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p16_bbarb_hepemg_model
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_color
      use p16_bbarb_hepemg_abbrevd541h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc541(35)
      complex(ki) :: Qspvak1k4
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspk1
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak5k2
      complex(ki) :: QspQ
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspe6
      complex(ki) :: Qspk2
      complex(ki) :: Qspvae6k2
      Qspvak1k4 = dotproduct(Q,spvak1k4)
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspk1 = dotproduct(Q,k1)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspk6 = dotproduct(Q,k6)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      QspQ = dotproduct(Q,Q)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspe6 = dotproduct(Q,e6)
      Qspk2 = dotproduct(Q,k2)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      acc541(1)=abb541(3)
      acc541(2)=abb541(4)
      acc541(3)=abb541(5)
      acc541(4)=abb541(6)
      acc541(5)=abb541(7)
      acc541(6)=abb541(8)
      acc541(7)=abb541(9)
      acc541(8)=abb541(10)
      acc541(9)=abb541(11)
      acc541(10)=abb541(12)
      acc541(11)=abb541(13)
      acc541(12)=abb541(14)
      acc541(13)=abb541(15)
      acc541(14)=abb541(16)
      acc541(15)=abb541(17)
      acc541(16)=abb541(18)
      acc541(17)=abb541(19)
      acc541(18)=abb541(20)
      acc541(19)=abb541(21)
      acc541(20)=abb541(22)
      acc541(21)=abb541(23)
      acc541(22)=abb541(24)
      acc541(23)=abb541(25)
      acc541(24)=abb541(26)
      acc541(25)=abb541(27)
      acc541(26)=abb541(28)
      acc541(27)=acc541(12)*Qspvak1k4
      acc541(28)=Qspvak1k2*acc541(3)
      acc541(29)=Qspk1*acc541(21)
      acc541(30)=Qspvak6k2*acc541(8)
      acc541(31)=Qspk6*acc541(18)
      acc541(32)=Qspvak5k2*acc541(14)
      acc541(33)=QspQ*acc541(4)
      acc541(34)=-Qspvak1k2*acc541(2)
      acc541(34)=acc541(13)+acc541(34)
      acc541(34)=Qspvak5k4*acc541(34)
      acc541(27)=acc541(34)+acc541(33)+acc541(32)+acc541(31)+acc541(30)+acc541(&
      &29)+acc541(28)+acc541(27)+acc541(11)
      acc541(27)=Qspe6*acc541(27)
      acc541(28)=Qspk6+Qspk2
      acc541(28)=acc541(26)*acc541(28)
      acc541(29)=Qspvak6k2*acc541(16)
      acc541(30)=Qspvae6k2*acc541(15)
      acc541(31)=-QspQ*acc541(25)
      acc541(28)=acc541(31)+acc541(30)+acc541(29)+acc541(6)+acc541(28)
      acc541(28)=Qspvak5k4*acc541(28)
      acc541(29)=Qspvak5k2*acc541(19)
      acc541(30)=-Qspvae6k2*acc541(23)
      acc541(29)=acc541(30)+acc541(5)+acc541(29)
      acc541(29)=QspQ*acc541(29)
      acc541(30)=-Qspk1*acc541(10)
      acc541(31)=Qspk2*acc541(20)
      acc541(32)=Qspvak6k2*acc541(1)
      acc541(33)=Qspk6*acc541(17)
      acc541(34)=-Qspk2+Qspk6
      acc541(34)=Qspvak5k2*acc541(22)*acc541(34)
      acc541(35)=Qspk1*acc541(24)
      acc541(35)=acc541(9)+acc541(35)
      acc541(35)=Qspvae6k2*acc541(35)
      brack=acc541(7)+acc541(27)+acc541(28)+acc541(29)+acc541(30)+acc541(31)+ac&
      &c541(32)+acc541(33)+acc541(34)+acc541(35)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p16_bbarb_hepemg_groups, only: &
!           & sign => diagram541_sign, shift => diagram541_shift
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd541h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d541
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k6
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d541 = 0.0_ki
      d541 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d541, ki), aimag(d541), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd541h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d541
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k6
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d541 = 0.0_ki
      d541 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d541, ki), aimag(d541), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p16_bbarb_hepemg_d541h0l1
