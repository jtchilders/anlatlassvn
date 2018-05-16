module     p12_sbars_hepemg_d613h3l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p12_sbars_hepemg/helicity3d613h3l1d.f90
   ! generator: buildfortran_d.py
   use p12_sbars_hepemg_config, only: ki
   use p12_sbars_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d613
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(28) :: acd613
      complex(ki) :: brack
      acd613(1)=dotproduct(k1,qshift)
      acd613(2)=dotproduct(e6,qshift)
      acd613(3)=abb613(8)
      acd613(4)=abb613(5)
      acd613(5)=dotproduct(k6,qshift)
      acd613(6)=abb613(20)
      acd613(7)=abb613(9)
      acd613(8)=dotproduct(qshift,spvak2k5)
      acd613(9)=abb613(11)
      acd613(10)=dotproduct(qshift,spvak4k5)
      acd613(11)=abb613(15)
      acd613(12)=dotproduct(qshift,spvak6k1)
      acd613(13)=abb613(16)
      acd613(14)=abb613(10)
      acd613(15)=dotproduct(qshift,qshift)
      acd613(16)=abb613(7)
      acd613(17)=abb613(6)
      acd613(18)=abb613(19)
      acd613(19)=dotproduct(qshift,spvae6k1)
      acd613(20)=abb613(12)
      acd613(21)=abb613(14)
      acd613(22)=acd613(12)*acd613(13)
      acd613(23)=acd613(10)*acd613(11)
      acd613(24)=acd613(8)*acd613(9)
      acd613(25)=acd613(5)*acd613(6)
      acd613(26)=acd613(1)*acd613(3)
      acd613(22)=acd613(26)+acd613(25)+acd613(24)+acd613(23)-acd613(14)+acd613(&
      &22)
      acd613(22)=acd613(2)*acd613(22)
      acd613(23)=-acd613(19)*acd613(20)
      acd613(24)=acd613(15)*acd613(16)
      acd613(25)=-acd613(10)*acd613(18)
      acd613(26)=-acd613(8)*acd613(17)
      acd613(27)=-acd613(5)*acd613(7)
      acd613(28)=-acd613(1)*acd613(4)
      brack=acd613(21)+acd613(22)+acd613(23)+acd613(24)+acd613(25)+acd613(26)+a&
      &cd613(27)+acd613(28)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(34) :: acd613
      complex(ki) :: brack
      acd613(1)=k1(iv1)
      acd613(2)=dotproduct(e6,qshift)
      acd613(3)=abb613(8)
      acd613(4)=abb613(5)
      acd613(5)=k6(iv1)
      acd613(6)=abb613(20)
      acd613(7)=abb613(9)
      acd613(8)=e6(iv1)
      acd613(9)=dotproduct(k1,qshift)
      acd613(10)=dotproduct(k6,qshift)
      acd613(11)=dotproduct(qshift,spvak2k5)
      acd613(12)=abb613(11)
      acd613(13)=dotproduct(qshift,spvak4k5)
      acd613(14)=abb613(15)
      acd613(15)=dotproduct(qshift,spvak6k1)
      acd613(16)=abb613(16)
      acd613(17)=abb613(10)
      acd613(18)=qshift(iv1)
      acd613(19)=abb613(7)
      acd613(20)=spvak2k5(iv1)
      acd613(21)=abb613(6)
      acd613(22)=spvak4k5(iv1)
      acd613(23)=abb613(19)
      acd613(24)=spvak6k1(iv1)
      acd613(25)=spvae6k1(iv1)
      acd613(26)=abb613(12)
      acd613(27)=acd613(16)*acd613(24)
      acd613(28)=acd613(14)*acd613(22)
      acd613(29)=acd613(12)*acd613(20)
      acd613(30)=acd613(5)*acd613(6)
      acd613(31)=acd613(1)*acd613(3)
      acd613(27)=acd613(31)+acd613(30)+acd613(29)+acd613(27)+acd613(28)
      acd613(27)=acd613(2)*acd613(27)
      acd613(28)=acd613(16)*acd613(15)
      acd613(29)=acd613(14)*acd613(13)
      acd613(30)=acd613(12)*acd613(11)
      acd613(31)=acd613(6)*acd613(10)
      acd613(32)=acd613(3)*acd613(9)
      acd613(28)=acd613(32)+acd613(31)+acd613(30)+acd613(29)-acd613(17)+acd613(&
      &28)
      acd613(28)=acd613(8)*acd613(28)
      acd613(29)=-acd613(25)*acd613(26)
      acd613(30)=acd613(18)*acd613(19)
      acd613(31)=-acd613(22)*acd613(23)
      acd613(32)=-acd613(20)*acd613(21)
      acd613(33)=-acd613(5)*acd613(7)
      acd613(34)=-acd613(1)*acd613(4)
      brack=acd613(27)+acd613(28)+acd613(29)+2.0_ki*acd613(30)+acd613(31)+acd61&
      &3(32)+acd613(33)+acd613(34)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p12_sbars_hepemg_model
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_color
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(25) :: acd613
      complex(ki) :: brack
      acd613(1)=d(iv1,iv2)
      acd613(2)=abb613(7)
      acd613(3)=k1(iv1)
      acd613(4)=e6(iv2)
      acd613(5)=abb613(8)
      acd613(6)=k1(iv2)
      acd613(7)=e6(iv1)
      acd613(8)=k6(iv1)
      acd613(9)=abb613(20)
      acd613(10)=k6(iv2)
      acd613(11)=spvak2k5(iv2)
      acd613(12)=abb613(11)
      acd613(13)=spvak4k5(iv2)
      acd613(14)=abb613(15)
      acd613(15)=spvak6k1(iv2)
      acd613(16)=abb613(16)
      acd613(17)=spvak2k5(iv1)
      acd613(18)=spvak4k5(iv1)
      acd613(19)=spvak6k1(iv1)
      acd613(20)=acd613(16)*acd613(15)
      acd613(21)=acd613(14)*acd613(13)
      acd613(22)=acd613(12)*acd613(11)
      acd613(23)=acd613(9)*acd613(10)
      acd613(24)=acd613(5)*acd613(6)
      acd613(20)=acd613(24)+acd613(23)+acd613(22)+acd613(20)+acd613(21)
      acd613(20)=acd613(7)*acd613(20)
      acd613(21)=acd613(16)*acd613(19)
      acd613(22)=acd613(14)*acd613(18)
      acd613(23)=acd613(12)*acd613(17)
      acd613(24)=acd613(9)*acd613(8)
      acd613(25)=acd613(5)*acd613(3)
      acd613(21)=acd613(25)+acd613(24)+acd613(23)+acd613(21)+acd613(22)
      acd613(21)=acd613(4)*acd613(21)
      acd613(22)=acd613(1)*acd613(2)
      brack=acd613(20)+acd613(21)+2.0_ki*acd613(22)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p12_sbars_hepemg_globalsl1, only: epspow
      use p12_sbars_hepemg_kinematics
      use p12_sbars_hepemg_abbrevd613h3
      implicit none
      complex(ki), intent(in) :: mu2
      integer, intent(in), optional :: i1
      integer, intent(in), optional :: i2
      complex(ki) :: numerator
      complex(ki) :: loc
      integer :: t1
      integer :: deg
      complex(ki), dimension(4), parameter :: Q = (/ (0.0_ki,0.0_ki),(0.0_ki,0.&
      &0_ki),(0.0_ki,0.0_ki),(0.0_ki,0.0_ki)/)
      qshift = -k3-k6-k5-k4
      numerator = 0.0_ki
      deg = 0
      if(present(i1)) then
          iv1=i1
          deg=1
      else
          iv1=1
      end if
      if(present(i2)) then
          iv2=i2
          deg=2
      else
          iv2=1
      end if
      t1 = 0
      if(deg.eq.0) then
         numerator = cond(epspow.eq.t1,brack_1,Q,mu2)
         return
      end if
      if(deg.eq.1) then
         numerator = cond(epspow.eq.t1,brack_2,Q,mu2)
         return
      end if
      if(deg.eq.2) then
         numerator = cond(epspow.eq.t1,brack_3,Q,mu2)
         return
      end if
   end function derivative
!---#] function derivative:
!---#[ subroutine reconstruct_d613:
   subroutine     reconstruct_d613(coeffs)
      use p12_sbars_hepemg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_613:
      coeffs%coeffs_613%c0 = derivative(czip)
      coeffs%coeffs_613%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_613%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_613%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_613%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_613%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_613%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_613%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_613%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_613%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_613%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_613%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_613%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_613%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_613%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_613:
   end subroutine reconstruct_d613
!---#] subroutine reconstruct_d613:
end module     p12_sbars_hepemg_d613h3l1d
