module     p4_ubaru_hepemg_d323h1l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity1d323h1l1.f90
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
      use p4_ubaru_hepemg_abbrevd323h1
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc323(19)
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspvak2k1
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak2k4
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspk6
      complex(ki) :: Qspk1
      complex(ki) :: QspQ
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspvak2k1 = dotproduct(Q,spvak2k1)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak2k4 = dotproduct(Q,spvak2k4)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspk6 = dotproduct(Q,k6)
      Qspk1 = dotproduct(Q,k1)
      QspQ = dotproduct(Q,Q)
      acc323(1)=abb323(5)
      acc323(2)=abb323(6)
      acc323(3)=abb323(7)
      acc323(4)=abb323(8)
      acc323(5)=abb323(9)
      acc323(6)=abb323(10)
      acc323(7)=abb323(14)
      acc323(8)=abb323(15)
      acc323(9)=abb323(16)
      acc323(10)=abb323(17)
      acc323(11)=abb323(18)
      acc323(12)=acc323(10)*Qspvak5k4
      acc323(12)=acc323(12)+acc323(9)
      acc323(12)=Qspvak2k6*acc323(12)
      acc323(13)=acc323(3)*Qspvak5k4
      acc323(13)=acc323(13)+acc323(5)
      acc323(13)=Qspvak2k1*acc323(13)
      acc323(14)=Qspvak6k1*acc323(8)
      acc323(15)=Qspvak2k4*acc323(11)
      acc323(16)=Qspvak1k6*acc323(1)
      acc323(17)=Qspk6*acc323(2)
      acc323(18)=Qspk1*acc323(4)
      acc323(19)=QspQ*acc323(7)
      brack=acc323(6)+acc323(12)+acc323(13)+acc323(14)+acc323(15)+acc323(16)+ac&
      &c323(17)+acc323(18)+acc323(19)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p4_ubaru_hepemg_groups, only: &
!           & sign => diagram323_sign, shift => diagram323_shift
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd323h1
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d323
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k2
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d323 = 0.0_ki
      d323 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d323, ki), aimag(d323), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd323h1
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d323
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k2
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d323 = 0.0_ki
      d323 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d323, ki), aimag(d323), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p4_ubaru_hepemg_d323h1l1
