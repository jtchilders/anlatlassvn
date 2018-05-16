module     p12_cbbar_hepneg_d30h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_cbbar_hepneg/helicity0d30h0l1.f90
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
      use p12_cbbar_hepneg_abbrevd30h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc30(20)
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k1
      complex(ki) :: Qspk6
      complex(ki) :: Qspk1
      complex(ki) :: Qspvak5k4
      complex(ki) :: QspQ
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak2k4
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k1 = dotproduct(Q,spvae6k1)
      Qspk6 = dotproduct(Q,k6)
      Qspk1 = dotproduct(Q,k1)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      QspQ = dotproduct(Q,Q)
      Qspk2 = dotproduct(Q,k2)
      Qspvak2k4 = dotproduct(Q,spvak2k4)
      acc30(1)=abb30(5)
      acc30(2)=abb30(6)
      acc30(3)=abb30(7)
      acc30(4)=abb30(8)
      acc30(5)=abb30(10)
      acc30(6)=abb30(11)
      acc30(7)=abb30(12)
      acc30(8)=abb30(13)
      acc30(9)=abb30(14)
      acc30(10)=abb30(15)
      acc30(11)=abb30(16)
      acc30(12)=abb30(20)
      acc30(13)=abb30(21)
      acc30(14)=acc30(9)*Qspvak2k6
      acc30(15)=acc30(6)*Qspvak2k1
      acc30(16)=Qspe6*acc30(7)
      acc30(17)=Qspvae6k1*acc30(5)
      acc30(18)=Qspk6-Qspk1
      acc30(19)=acc30(8)*acc30(18)
      acc30(14)=acc30(19)+acc30(17)+acc30(16)+acc30(3)+acc30(14)+acc30(15)
      acc30(14)=Qspvak5k4*acc30(14)
      acc30(15)=acc30(13)*QspQ
      acc30(16)=acc30(12)*Qspk2
      acc30(17)=acc30(10)*Qspvak2k4
      acc30(19)=Qspe6*acc30(4)
      acc30(20)=Qspvae6k1*acc30(1)
      acc30(18)=-acc30(11)*acc30(18)
      brack=acc30(2)+acc30(14)+acc30(15)+acc30(16)+acc30(17)+acc30(18)+acc30(19&
      &)+acc30(20)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p12_cbbar_hepneg_groups, only: &
!           & sign => diagram30_sign, shift => diagram30_shift
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd30h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d30
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d30 = 0.0_ki
      d30 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d30, ki), aimag(d30), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p12_cbbar_hepneg_globalsl1, only: epspow
      use p12_cbbar_hepneg_kinematics
      use p12_cbbar_hepneg_abbrevd30h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d30
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k2
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d30 = 0.0_ki
      d30 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d30, ki), aimag(d30), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p12_cbbar_hepneg_d30h0l1
