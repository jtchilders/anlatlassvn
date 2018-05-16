module     p12_sbars_hepemg_d613h3l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity3d613h3l1.f90
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
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc613(20)
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspk1
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak2k5
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k1
      complex(ki) :: QspQ
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspk1 = dotproduct(Q,k1)
      Qspk6 = dotproduct(Q,k6)
      Qspvak2k5 = dotproduct(Q,spvak2k5)
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      QspQ = dotproduct(Q,Q)
      acc613(1)=abb613(5)
      acc613(2)=abb613(6)
      acc613(3)=abb613(7)
      acc613(4)=abb613(8)
      acc613(5)=abb613(9)
      acc613(6)=abb613(10)
      acc613(7)=abb613(11)
      acc613(8)=abb613(12)
      acc613(9)=abb613(14)
      acc613(10)=abb613(15)
      acc613(11)=abb613(16)
      acc613(12)=abb613(19)
      acc613(13)=abb613(20)
      acc613(14)=acc613(11)*Qspvak6k1
      acc613(15)=Qspk1*acc613(4)
      acc613(16)=Qspk6*acc613(13)
      acc613(17)=Qspvak2k5*acc613(7)
      acc613(18)=Qspvak4k5*acc613(10)
      acc613(14)=acc613(18)+acc613(17)+acc613(16)+acc613(15)+acc613(14)+acc613(&
      &6)
      acc613(14)=Qspe6*acc613(14)
      acc613(15)=acc613(8)*Qspvae6k1
      acc613(16)=acc613(3)*QspQ
      acc613(17)=Qspk1*acc613(1)
      acc613(18)=Qspk6*acc613(5)
      acc613(19)=Qspvak2k5*acc613(2)
      acc613(20)=Qspvak4k5*acc613(12)
      brack=acc613(9)+acc613(14)+acc613(15)+acc613(16)+acc613(17)+acc613(18)+ac&
      &c613(19)+acc613(20)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p12_sbars_hepemg_groups, only: &
!           & sign => diagram613_sign, shift => diagram613_shift
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d613
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k6-k5-k4
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d613 = 0.0_ki
      d613 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d613, ki), aimag(d613), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d613
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k3-k6-k5-k4
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d613 = 0.0_ki
      d613 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d613, ki), aimag(d613), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p12_sbars_hepemg_d613h3l1
