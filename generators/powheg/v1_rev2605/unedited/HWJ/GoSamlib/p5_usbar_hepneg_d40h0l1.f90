module     p5_usbar_hepneg_d40h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p5_usbar_hepneg/helicity0d40h0l1.f90
   ! generator: buildfortran.py
   use p5_usbar_hepneg_config, only: ki
   use p5_usbar_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p5_usbar_hepneg_model
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_color
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc40(11)
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspvak5k4
      complex(ki) :: Qspvak5k1
      complex(ki) :: Qspk6
      complex(ki) :: Qspk2
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspvak5k4 = dotproduct(Q,spvak5k4)
      Qspvak5k1 = dotproduct(Q,spvak5k1)
      Qspk6 = dotproduct(Q,k6)
      Qspk2 = dotproduct(Q,k2)
      acc40(1)=abb40(5)
      acc40(2)=abb40(6)
      acc40(3)=abb40(7)
      acc40(4)=abb40(8)
      acc40(5)=abb40(9)
      acc40(6)=abb40(10)
      acc40(7)=Qspvak6k2*acc40(2)
      acc40(8)=Qspvak5k4*acc40(6)
      acc40(9)=Qspvak5k1*acc40(1)
      acc40(10)=Qspk6*acc40(5)
      acc40(11)=Qspk2*acc40(3)
      brack=acc40(4)+acc40(7)+acc40(8)+acc40(9)+acc40(10)+acc40(11)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p5_usbar_hepneg_groups, only: &
!           & sign => diagram40_sign, shift => diagram40_shift
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d40
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k6
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d40 = 0.0_ki
      d40 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d40, ki), aimag(d40), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p5_usbar_hepneg_globalsl1, only: epspow
      use p5_usbar_hepneg_kinematics
      use p5_usbar_hepneg_abbrevd40h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d40
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k6
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d40 = 0.0_ki
      d40 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d40, ki), aimag(d40), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p5_usbar_hepneg_d40h0l1
