module     p1_dbarc_hepneg_d70h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p1_dbarc_hepneg/helicity0d70h0l1.f90
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
      use p1_dbarc_hepneg_abbrevd70h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc70(36)
      complex(ki) :: Qspvak1k2
      complex(ki) :: Qspvak1k4
      complex(ki) :: Qspvak5k2
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6l3
      complex(ki) :: Qspvak6e6
      complex(ki) :: Qspvae6k4
      complex(ki) :: Qspvae6k6
      complex(ki) :: Qspvak5e6
      complex(ki) :: Qspval3e6
      complex(ki) :: Qspvak6k4
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak5k6
      complex(ki) :: QspQ
      complex(ki) :: Qspvak1k6
      complex(ki) :: Qspvak1e6
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak2e6
      Qspvak1k2 = dotproduct(Q,spvak1k2)
      Qspvak1k4 = dotproduct(Q,spvak1k4)
      Qspvak5k2 = dotproduct(Q,spvak5k2)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6l3 = dotproduct(Q,spvae6l3)
      Qspvak6e6 = dotproduct(Q,spvak6e6)
      Qspvae6k4 = dotproduct(Q,spvae6k4)
      Qspvae6k6 = dotproduct(Q,spvae6k6)
      Qspvak5e6 = dotproduct(Q,spvak5e6)
      Qspval3e6 = dotproduct(Q,spval3e6)
      Qspvak6k4 = dotproduct(Q,spvak6k4)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspk6 = dotproduct(Q,k6)
      Qspvak5k6 = dotproduct(Q,spvak5k6)
      QspQ = dotproduct(Q,Q)
      Qspvak1k6 = dotproduct(Q,spvak1k6)
      Qspvak1e6 = dotproduct(Q,spvak1e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak2e6 = dotproduct(Q,spvak2e6)
      acc70(1)=abb70(8)
      acc70(2)=abb70(9)
      acc70(3)=abb70(10)
      acc70(4)=abb70(11)
      acc70(5)=abb70(12)
      acc70(6)=abb70(13)
      acc70(7)=abb70(14)
      acc70(8)=abb70(15)
      acc70(9)=abb70(16)
      acc70(10)=abb70(17)
      acc70(11)=abb70(18)
      acc70(12)=abb70(20)
      acc70(13)=abb70(21)
      acc70(14)=abb70(22)
      acc70(15)=abb70(25)
      acc70(16)=abb70(26)
      acc70(17)=abb70(27)
      acc70(18)=abb70(28)
      acc70(19)=abb70(29)
      acc70(20)=abb70(32)
      acc70(21)=acc70(5)*Qspvak1k2
      acc70(22)=acc70(4)*Qspvak1k4
      acc70(23)=acc70(2)*Qspvak5k2
      acc70(21)=acc70(23)+acc70(22)+acc70(7)+acc70(21)
      acc70(21)=Qspe6*acc70(21)
      acc70(22)=acc70(20)*Qspvae6l3
      acc70(23)=acc70(19)*Qspvak6e6
      acc70(24)=acc70(18)*Qspvae6k4
      acc70(25)=acc70(17)*Qspvae6k6
      acc70(26)=acc70(16)*Qspvak5e6
      acc70(27)=acc70(15)*Qspval3e6
      acc70(28)=acc70(14)*Qspvak6k4
      acc70(29)=acc70(13)*Qspvak6k2
      acc70(30)=acc70(12)*Qspk6
      acc70(31)=acc70(11)*Qspvak5k6
      acc70(32)=acc70(10)*QspQ
      acc70(33)=acc70(9)*Qspvak1k6
      acc70(34)=acc70(8)*Qspvak1e6
      acc70(35)=acc70(6)*Qspvae6k2
      acc70(36)=acc70(3)*Qspvak2e6
      brack=acc70(1)+acc70(21)+acc70(22)+acc70(23)+acc70(24)+acc70(25)+acc70(26&
      &)+acc70(27)+acc70(28)+acc70(29)+acc70(30)+acc70(31)+acc70(32)+acc70(33)+&
      &acc70(34)+acc70(35)+acc70(36)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p1_dbarc_hepneg_groups, only: &
!           & sign => diagram70_sign, shift => diagram70_shift
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd70h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d70
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(-Q_ext(4),  ki_sam),aimag(-Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(-Q_ext(1:3),ki_sam),aimag(-Q_ext(1:3)),ki)
      d70 = 0.0_ki
      d70 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d70, ki), aimag(d70), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd70h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d70
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(-Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d70 = 0.0_ki
      d70 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d70, ki), aimag(d70), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p1_dbarc_hepneg_d70h0l1
