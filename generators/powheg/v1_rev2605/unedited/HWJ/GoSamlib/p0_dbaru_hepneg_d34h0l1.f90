module     p0_dbaru_hepneg_d34h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbaru_hepneg/helicity0d34h0l1.f90
   ! generator: buildfortran.py
   use p0_dbaru_hepneg_config, only: ki
   use p0_dbaru_hepneg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p0_dbaru_hepneg_model
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_color
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc34(19)
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
      acc34(1)=abb34(5)
      acc34(2)=abb34(6)
      acc34(3)=abb34(7)
      acc34(4)=abb34(8)
      acc34(5)=abb34(10)
      acc34(6)=abb34(11)
      acc34(7)=abb34(13)
      acc34(8)=abb34(14)
      acc34(9)=abb34(15)
      acc34(10)=abb34(16)
      acc34(11)=abb34(17)
      acc34(12)=abb34(18)
      acc34(13)=acc34(6)*Qspvak1k6
      acc34(14)=Qspe6*acc34(3)
      acc34(15)=Qspvae6k2*acc34(10)
      acc34(16)=Qspk6-Qspk2
      acc34(17)=acc34(7)*acc34(16)
      acc34(13)=acc34(17)+acc34(15)+acc34(14)+acc34(8)+acc34(13)
      acc34(13)=Qspvak5k4*acc34(13)
      acc34(14)=acc34(12)*QspQ
      acc34(15)=acc34(11)*Qspk1
      acc34(17)=acc34(4)*Qspvak1k4
      acc34(18)=Qspe6*acc34(2)
      acc34(19)=Qspvae6k2*acc34(5)
      acc34(16)=-acc34(9)*acc34(16)
      brack=acc34(1)+acc34(13)+acc34(14)+acc34(15)+acc34(16)+acc34(17)+acc34(18&
      &)+acc34(19)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p0_dbaru_hepneg_groups, only: &
!           & sign => diagram34_sign, shift => diagram34_shift
      use p0_dbaru_hepneg_globalsl1, only: epspow
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d34
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k2
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d34 = 0.0_ki
      d34 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d34, ki), aimag(d34), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p0_dbaru_hepneg_globalsl1, only: epspow
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d34
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k2
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d34 = 0.0_ki
      d34 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d34, ki), aimag(d34), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p0_dbaru_hepneg_d34h0l1
