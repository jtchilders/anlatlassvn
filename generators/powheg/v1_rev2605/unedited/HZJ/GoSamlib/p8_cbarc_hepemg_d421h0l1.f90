module     p8_cbarc_hepemg_d421h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p8_cbarc_hepemg/helicity0d421h0l1.f90
   ! generator: buildfortran.py
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_util, only: cond
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   public :: numerator_samurai
   public :: numerator_golem95
contains
!---#[ function brack_1:
   pure function brack_1(Q,mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd421h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc421(11)
      complex(ki) :: Qspvae6k2
      complex(ki) :: Qspvak2k6
      complex(ki) :: Qspe6
      complex(ki) :: Qspk6
      complex(ki) :: Qspk2
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      Qspe6 = dotproduct(Q,e6)
      Qspk6 = dotproduct(Q,k6)
      Qspk2 = dotproduct(Q,k2)
      acc421(1)=abb421(5)
      acc421(2)=abb421(6)
      acc421(3)=abb421(7)
      acc421(4)=abb421(8)
      acc421(5)=abb421(9)
      acc421(6)=abb421(10)
      acc421(7)=Qspvae6k2*acc421(3)
      acc421(8)=Qspvak2k6*acc421(1)
      acc421(9)=Qspe6*acc421(2)
      acc421(10)=Qspk6*acc421(5)
      acc421(11)=Qspk2*acc421(4)
      brack=acc421(6)+acc421(7)+acc421(8)+acc421(9)+acc421(10)+acc421(11)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p8_cbarc_hepemg_groups, only: &
!           & sign => diagram421_sign, shift => diagram421_shift
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd421h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d421
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(0:3) :: qshift
      qshift = k6
      Q(1)  =cmplx(real(+Q_ext(4)  -qshift(0),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3)-qshift(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d421 = 0.0_ki
      d421 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d421, ki), aimag(d421), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd421h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d421
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      real(ki), dimension(4) :: qshift
      qshift = k6
      Q(:)  =cmplx(real(+Q_ext(:)  -qshift(:),  ki_gol), 0.0_ki_gol, ki)
      d421 = 0.0_ki
      d421 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d421, ki), aimag(d421), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p8_cbarc_hepemg_d421h0l1
