module     p12_sbars_hepemg_d333h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity0d333h0l1.f90
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
      use p12_sbars_hepemg_abbrevd333h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc333(19)
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspk6
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak5k4
      complex(ki) :: QspQ
      complex(ki) :: Qspk1
      complex(ki) :: Qspvak1k4
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspk6 = dotproduct(Q,k6)
      Qspk2 = dotproduct(Q,k2)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      QspQ = dotproduct(Q,Q)
      Qspk1 = dotproduct(Q,k1)
      Qspvak1k4 = dotproduct(Q,spvak1k4)
      acc333(1)=abb333(5)
      acc333(2)=abb333(6)
      acc333(3)=abb333(7)
      acc333(4)=abb333(8)
      acc333(5)=abb333(10)
      acc333(6)=abb333(11)
      acc333(7)=abb333(13)
      acc333(8)=abb333(14)
      acc333(9)=abb333(15)
      acc333(10)=abb333(16)
      acc333(11)=abb333(17)
      acc333(12)=abb333(18)
      acc333(13)=acc333(6)*Qspvak1k6
      acc333(14)=Qspe6*acc333(3)
      acc333(15)=Qspvae6k2*acc333(10)
      acc333(16)=Qspk6-Qspk2
      acc333(17)=acc333(7)*acc333(16)
      acc333(13)=acc333(17)+acc333(15)+acc333(14)+acc333(8)+acc333(13)
      acc333(13)=Qspvak5k4*acc333(13)
      acc333(14)=acc333(12)*QspQ
      acc333(15)=acc333(11)*Qspk1
      acc333(17)=acc333(4)*Qspvak1k4
      acc333(18)=Qspe6*acc333(2)
      acc333(19)=Qspvae6k2*acc333(5)
      acc333(16)=-acc333(9)*acc333(16)
      brack=acc333(1)+acc333(13)+acc333(14)+acc333(15)+acc333(16)+acc333(17)+ac&
      &c333(18)+acc333(19)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p12_sbars_hepemg_groups, only: &
!           & sign => diagram333_sign, shift => diagram333_shift
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd333h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d333
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d333 = 0.0_ki
      d333 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d333, ki), aimag(d333), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd333h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d333
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d333 = 0.0_ki
      d333 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d333, ki), aimag(d333), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p12_sbars_hepemg_d333h0l1
