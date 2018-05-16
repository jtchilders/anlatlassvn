module     p16_bbarb_hepemg_d605h3l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p16_bbarb_hepemg/helicity3d605h3l1.f90
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
      use p16_bbarb_hepemg_abbrevd605h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc605(17)
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspk2
      complex(ki) :: Qspe6
      complex(ki) :: Qspvak2e6
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspk6
      complex(ki) :: QspQ
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspk2 = dotproduct(Q,k2)
      Qspe6 = dotproduct(Q,e6)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspk6 = dotproduct(Q,k6)
      QspQ = dotproduct(Q,Q)
      acc605(1)=abb605(5)
      acc605(2)=abb605(6)
      acc605(3)=abb605(7)
      acc605(4)=abb605(11)
      acc605(5)=abb605(12)
      acc605(6)=abb605(13)
      acc605(7)=abb605(14)
      acc605(8)=abb605(16)
      acc605(9)=abb605(17)
      acc605(10)=abb605(18)
      acc605(11)=acc605(2)*Qspvak2k6
      acc605(12)=acc605(8)*Qspk2
      acc605(11)=acc605(12)+acc605(5)+acc605(11)
      acc605(11)=Qspe6*acc605(11)
      acc605(12)=acc605(3)*Qspvak2k6
      acc605(13)=acc605(4)*Qspk2
      acc605(14)=Qspvak2e6*acc605(10)
      acc605(15)=Qspvak6k2*acc605(1)
      acc605(16)=Qspk6*acc605(7)
      acc605(17)=QspQ*acc605(9)
      brack=acc605(6)+acc605(11)+acc605(12)+acc605(13)+acc605(14)+acc605(15)+ac&
      &c605(16)+acc605(17)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p16_bbarb_hepemg_groups, only: &
!           & sign => diagram605_sign, shift => diagram605_shift
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd605h3
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d605
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(+Q_ext(4),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d605 = 0.0_ki
      d605 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d605, ki), aimag(d605), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p16_bbarb_hepemg_globalsl1, only: epspow
      use p16_bbarb_hepemg_kinematics
      use p16_bbarb_hepemg_abbrevd605h3
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d605
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(+Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d605 = 0.0_ki
      d605 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d605, ki), aimag(d605), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p16_bbarb_hepemg_d605h3l1
