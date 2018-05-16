module     p4_ubaru_hepemg_d413h2l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p4_ubaru_hepemg/helicity2d413h2l1.f90
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
      use p4_ubaru_hepemg_abbrevd413h2
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc413(13)
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak4k5
      complex(ki) :: Qspvak4k2
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspk6
      complex(ki) :: Qspk1
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak4k5 = dotproduct(Q,spvak4k5)
      Qspvak4k2 = dotproduct(Q,spvak4k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspk6 = dotproduct(Q,k6)
      Qspk1 = dotproduct(Q,k1)
      acc413(1)=abb413(5)
      acc413(2)=abb413(6)
      acc413(3)=abb413(7)
      acc413(4)=abb413(8)
      acc413(5)=abb413(9)
      acc413(6)=abb413(12)
      acc413(7)=abb413(15)
      acc413(8)=Qspvak6k1*acc413(5)
      acc413(9)=Qspvak4k5*acc413(4)
      acc413(10)=Qspvak4k2*acc413(2)
      acc413(11)=Qspvak1k6*acc413(1)
      acc413(12)=Qspk6*acc413(6)
      acc413(13)=Qspk1*acc413(7)
      brack=acc413(3)+acc413(8)+acc413(9)+acc413(10)+acc413(11)+acc413(12)+acc4&
      &13(13)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p4_ubaru_hepemg_groups, only: &
!           & sign => diagram413_sign, shift => diagram413_shift
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd413h2
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d413
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k5-k4
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d413 = 0.0_ki
      d413 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d413, ki), aimag(d413), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p4_ubaru_hepemg_globalsl1, only: epspow
      use p4_ubaru_hepemg_kinematics
      use p4_ubaru_hepemg_abbrevd413h2
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d413
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k3-k5-k4
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d413 = 0.0_ki
      d413 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d413, ki), aimag(d413), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p4_ubaru_hepemg_d413h2l1
