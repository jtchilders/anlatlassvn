module     p12_cbbar_hepneg_d94h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d94h0l1.f90
   ! generator: buildfortran.py
   use p12_cbbar_hepneg_config, only: ki
   use p12_cbbar_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p12_cbbar_hepneg_model
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_color
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc94(17)
      complex(ki) :: Qspk1
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspvak1k6
      complex(ki) :: QspQ
      complex(ki) :: Qspk6
      Qspk1 = dotproduct(Q,k1)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      QspQ = dotproduct(Q,Q)
      Qspk6 = dotproduct(Q,k6)
      acc94(1)=abb94(5)
      acc94(2)=abb94(6)
      acc94(3)=abb94(7)
      acc94(4)=abb94(8)
      acc94(5)=abb94(9)
      acc94(6)=abb94(10)
      acc94(7)=abb94(11)
      acc94(8)=abb94(12)
      acc94(9)=abb94(15)
      acc94(10)=abb94(17)
      acc94(11)=Qspk1*acc94(10)
      acc94(12)=Qspvak6k1*acc94(5)
      acc94(11)=acc94(12)+acc94(7)+acc94(11)
      acc94(11)=Qspe6*acc94(11)
      acc94(12)=acc94(9)*Qspvae6k1
      acc94(13)=acc94(8)*Qspvak1k6
      acc94(14)=acc94(6)*QspQ
      acc94(15)=acc94(4)*Qspk6
      acc94(16)=Qspk1*acc94(3)
      acc94(17)=Qspvak6k1*acc94(2)
      brack=acc94(1)+acc94(11)+acc94(12)+acc94(13)+acc94(14)+acc94(15)+acc94(16&
      &)+acc94(17)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p12_cbbar_hepneg_groups, only: &
!           & sign => diagram94_sign, shift => diagram94_shift
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d94
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k6-k5-k4
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d94 = 0.0_ki
      d94 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d94, ki), aimag(d94), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd94h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d94
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k3-k6-k5-k4
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d94 = 0.0_ki
      d94 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d94, ki), aimag(d94), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p12_cbbar_hepneg_d94h0l1
