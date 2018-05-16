module     p1_dbarc_hepneg_d92h0l1
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p1_dbarc_hepneg/helicity0d92h0l1.f90
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
      use p1_dbarc_hepneg_abbrevd92h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki) :: brack
      complex(ki) :: acc92(17)
      complex(ki) :: Qspk2
      complex(ki) :: Qspvak6k2
      complex(ki) :: Qspe6
      complex(ki) :: Qspvae6k2
      complex(ki) :: QspQ
      complex(ki) :: Qspk6
      complex(ki) :: Qspvak2k6
      Qspk2 = dotproduct(Q,k2)
      Qspvak6k2 = dotproduct(Q,spvak6k2)
      Qspe6 = dotproduct(Q,e6)
      Qspvae6k2 = dotproduct(Q,spvae6k2)
      QspQ = dotproduct(Q,Q)
      Qspk6 = dotproduct(Q,k6)
      Qspvak2k6 = dotproduct(Q,spvak2k6)
      acc92(1)=abb92(5)
      acc92(2)=abb92(6)
      acc92(3)=abb92(9)
      acc92(4)=abb92(11)
      acc92(5)=abb92(12)
      acc92(6)=abb92(13)
      acc92(7)=abb92(14)
      acc92(8)=abb92(16)
      acc92(9)=abb92(17)
      acc92(10)=abb92(18)
      acc92(11)=Qspk2*acc92(8)
      acc92(12)=Qspvak6k2*acc92(2)
      acc92(11)=acc92(12)+acc92(5)+acc92(11)
      acc92(11)=Qspe6*acc92(11)
      acc92(12)=acc92(10)*Qspvae6k2
      acc92(13)=acc92(9)*QspQ
      acc92(14)=acc92(7)*Qspk6
      acc92(15)=acc92(1)*Qspvak2k6
      acc92(16)=Qspk2*acc92(4)
      acc92(17)=Qspvak6k2*acc92(3)
      brack=acc92(6)+acc92(11)+acc92(12)+acc92(13)+acc92(14)+acc92(15)+acc92(16&
      &)+acc92(17)
   end  function brack_1
!---#] function brack_1:
!---#[ numerator interfaces:
   !------#[ function numerator_samurai:
   function numerator_samurai(ncut,Q_ext, mu2_ext) result(numerator)
      use precision, only: ki_sam => ki
!      use p1_dbarc_hepneg_groups, only: &
!           & sign => diagram92_sign, shift => diagram92_shift
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd92h0
      implicit none
      integer, intent(in) :: ncut
      complex(ki_sam), dimension(4), intent(in) :: Q_ext
      complex(ki_sam), intent(in) :: mu2_ext
      complex(ki_sam) :: numerator
      complex(ki) :: d92
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(1)  =cmplx(real(+Q_ext(4),  ki_sam),aimag(+Q_ext(4)),  ki)
      Q(2:4)=cmplx(real(+Q_ext(1:3),ki_sam),aimag(+Q_ext(1:3)),ki)
      d92 = 0.0_ki
      d92 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d92, ki), aimag(d92), ki_sam)
   end function numerator_samurai
   !------#] function numerator_samurai:
   !------#[ function numerator_golem95:
   function numerator_golem95(Q_ext, mu2_ext) result(numerator)
      use precision_golem, only: ki_gol => ki
      use p1_dbarc_hepneg_globalsl1, only: epspow
      use p1_dbarc_hepneg_kinematics
      use p1_dbarc_hepneg_abbrevd92h0
      implicit none
      real(ki_gol), dimension(0:3), intent(in) :: Q_ext
      real(ki_gol), intent(in) :: mu2_ext
      complex(ki_gol) :: numerator
      complex(ki) :: d92
      ! The Q that goes into the diagram
      complex(ki), dimension(4) :: Q
      complex(ki) :: mu2
      Q(:)  =cmplx(real(+Q_ext(:),  ki_gol), 0.0_ki_gol, ki)
      d92 = 0.0_ki
      d92 = (cond(epspow.eq.0,brack_1,Q,mu2))
      numerator = cmplx(real(d92, ki), aimag(d92), ki_gol)
   end function numerator_golem95
   !------#] function numerator_golem95:
!---#] numerator interfaces:
end module p1_dbarc_hepneg_d92h0l1
