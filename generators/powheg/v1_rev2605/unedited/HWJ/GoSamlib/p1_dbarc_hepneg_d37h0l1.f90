module     p1_dbarc_hepneg_d37h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p1_dbarc_hepneg/helicity0d37h0l1.f90
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
      use p1_dbarc_hepneg_abbrevd37h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc37(13)
      complex(ki) :: Qspk1
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak6k1
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspvak1k6
      Qspk1 = dotproduct(Q,k1)
      Qspk6 = dotproduct(Q,k6)
      Qspvak6k1 = dotproduct(Q,spvak6k1)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      acc37(1)=abb37(5)
      acc37(2)=abb37(6)
      acc37(3)=abb37(7)
      acc37(4)=abb37(8)
      acc37(5)=abb37(10)
      acc37(6)=abb37(12)
      acc37(7)=abb37(15)
      acc37(8)=acc37(7)*Qspk1
      acc37(9)=acc37(6)*Qspk6
      acc37(10)=acc37(5)*Qspvak6k1
      acc37(11)=acc37(4)*Qspvak5k4
      acc37(12)=acc37(2)*Qspvak5k2
      acc37(13)=acc37(1)*Qspvak1k6
      brack=acc37(3)+acc37(8)+acc37(9)+acc37(10)+acc37(11)+acc37(12)+acc37(13)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p1_dbarc_hepneg_groups, only: &
!           & sign => diagram37_sign, shift => diagram37_shift
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd37h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d37
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = -k3-k5-k4
      Q(1)  =cmplx(real(-Q_ext(4)  -qshift(0),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3)-qshift(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d37 = 0.0_ki
      d37 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d37, ki), aimag(d37), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd37h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d37
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = -k3-k5-k4
      Q(:)  =cmplx(real(-Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d37 = 0.0_ki
      d37 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d37, ki), aimag(d37), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p1_dbarc_hepneg_d37h0l1
