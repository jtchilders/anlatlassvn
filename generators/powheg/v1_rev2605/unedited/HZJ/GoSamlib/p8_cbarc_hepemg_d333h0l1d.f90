module     p8_cbarc_hepemg_d333h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HZJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p8_cbarc_hepemg/helicity0d333h0l1d.f90
   ! generator: buildfortran_d.py
   use p8_cbarc_hepemg_config, only: ki
   use p8_cbarc_hepemg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d333
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd333h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(28) :: acd333
      complex(ki) :: brack
      acd333(1)=dotproduct(k1,qshift)
      acd333(2)=abb333(17)
      acd333(3)=dotproduct(k2,qshift)
      acd333(4)=dotproduct(qshift,spvak5k4)
      acd333(5)=abb333(13)
      acd333(6)=abb333(15)
      acd333(7)=dotproduct(k6,qshift)
      acd333(8)=dotproduct(e6,qshift)
      acd333(9)=abb333(7)
      acd333(10)=abb333(6)
      acd333(11)=dotproduct(qshift,qshift)
      acd333(12)=abb333(18)
      acd333(13)=dotproduct(qshift,spvak1k6)
      acd333(14)=abb333(11)
      acd333(15)=dotproduct(qshift,spvae6k2)
      acd333(16)=abb333(16)
      acd333(17)=abb333(14)
      acd333(18)=dotproduct(qshift,spvak1k4)
      acd333(19)=abb333(8)
      acd333(20)=abb333(10)
      acd333(21)=abb333(5)
      acd333(22)=acd333(9)*acd333(8)
      acd333(23)=acd333(14)*acd333(13)
      acd333(24)=acd333(16)*acd333(15)
      acd333(22)=-acd333(17)+acd333(24)+acd333(23)+acd333(22)
      acd333(22)=acd333(4)*acd333(22)
      acd333(23)=acd333(7)-acd333(3)
      acd333(24)=acd333(5)*acd333(4)
      acd333(24)=acd333(24)+acd333(6)
      acd333(23)=acd333(24)*acd333(23)
      acd333(24)=-acd333(2)*acd333(1)
      acd333(25)=-acd333(10)*acd333(8)
      acd333(26)=acd333(12)*acd333(11)
      acd333(27)=-acd333(19)*acd333(18)
      acd333(28)=-acd333(20)*acd333(15)
      brack=acd333(21)+acd333(22)+acd333(23)+acd333(24)+acd333(25)+acd333(26)+a&
      &cd333(27)+acd333(28)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd333h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(34) :: acd333
      complex(ki) :: brack
      acd333(1)=k1(iv1)
      acd333(2)=abb333(17)
      acd333(3)=k2(iv1)
      acd333(4)=dotproduct(qshift,spvak5k4)
      acd333(5)=abb333(13)
      acd333(6)=abb333(15)
      acd333(7)=k6(iv1)
      acd333(8)=e6(iv1)
      acd333(9)=abb333(7)
      acd333(10)=abb333(6)
      acd333(11)=qshift(iv1)
      acd333(12)=abb333(18)
      acd333(13)=spvak5k4(iv1)
      acd333(14)=dotproduct(k2,qshift)
      acd333(15)=dotproduct(k6,qshift)
      acd333(16)=dotproduct(e6,qshift)
      acd333(17)=dotproduct(qshift,spvak1k6)
      acd333(18)=abb333(11)
      acd333(19)=dotproduct(qshift,spvae6k2)
      acd333(20)=abb333(16)
      acd333(21)=abb333(14)
      acd333(22)=spvak1k4(iv1)
      acd333(23)=abb333(8)
      acd333(24)=spvak1k6(iv1)
      acd333(25)=spvae6k2(iv1)
      acd333(26)=abb333(10)
      acd333(27)=-acd333(20)*acd333(25)
      acd333(28)=-acd333(18)*acd333(24)
      acd333(29)=-acd333(8)*acd333(9)
      acd333(30)=acd333(3)-acd333(7)
      acd333(31)=acd333(5)*acd333(30)
      acd333(27)=acd333(31)+acd333(29)+acd333(27)+acd333(28)
      acd333(27)=acd333(4)*acd333(27)
      acd333(28)=-acd333(20)*acd333(19)
      acd333(29)=-acd333(18)*acd333(17)
      acd333(31)=-acd333(9)*acd333(16)
      acd333(32)=-acd333(15)+acd333(14)
      acd333(32)=acd333(5)*acd333(32)
      acd333(28)=acd333(32)+acd333(31)+acd333(29)+acd333(21)+acd333(28)
      acd333(28)=acd333(13)*acd333(28)
      acd333(29)=acd333(6)*acd333(30)
      acd333(30)=acd333(22)*acd333(23)
      acd333(31)=acd333(11)*acd333(12)
      acd333(32)=acd333(1)*acd333(2)
      acd333(33)=acd333(25)*acd333(26)
      acd333(34)=acd333(8)*acd333(10)
      brack=acd333(27)+acd333(28)+acd333(29)+acd333(30)-2.0_ki*acd333(31)+acd33&
      &3(32)+acd333(33)+acd333(34)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p8_cbarc_hepemg_model
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_color
      use p8_cbarc_hepemg_abbrevd333h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(23) :: acd333
      complex(ki) :: brack
      acd333(1)=d(iv1,iv2)
      acd333(2)=abb333(18)
      acd333(3)=k2(iv1)
      acd333(4)=spvak5k4(iv2)
      acd333(5)=abb333(13)
      acd333(6)=k2(iv2)
      acd333(7)=spvak5k4(iv1)
      acd333(8)=k6(iv1)
      acd333(9)=k6(iv2)
      acd333(10)=e6(iv1)
      acd333(11)=abb333(7)
      acd333(12)=e6(iv2)
      acd333(13)=spvak1k6(iv2)
      acd333(14)=abb333(11)
      acd333(15)=spvae6k2(iv2)
      acd333(16)=abb333(16)
      acd333(17)=spvak1k6(iv1)
      acd333(18)=spvae6k2(iv1)
      acd333(19)=acd333(16)*acd333(15)
      acd333(20)=acd333(14)*acd333(13)
      acd333(21)=acd333(11)*acd333(12)
      acd333(22)=acd333(9)-acd333(6)
      acd333(22)=acd333(5)*acd333(22)
      acd333(19)=acd333(22)+acd333(21)+acd333(19)+acd333(20)
      acd333(19)=acd333(7)*acd333(19)
      acd333(20)=acd333(16)*acd333(18)
      acd333(21)=acd333(14)*acd333(17)
      acd333(22)=acd333(11)*acd333(10)
      acd333(23)=acd333(8)-acd333(3)
      acd333(23)=acd333(5)*acd333(23)
      acd333(20)=acd333(23)+acd333(22)+acd333(20)+acd333(21)
      acd333(20)=acd333(4)*acd333(20)
      acd333(21)=acd333(1)*acd333(2)
      brack=acd333(19)+acd333(20)+2.0_ki*acd333(21)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p8_cbarc_hepemg_globalsl1, only: epspow
      use p8_cbarc_hepemg_kinematics
      use p8_cbarc_hepemg_abbrevd333h0
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
      qshift = k2
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
!---#[ subroutine reconstruct_d333:
   subroutine     reconstruct_d333(coeffs)
      use p8_cbarc_hepemg_groups, only: tensrec_info_group1
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group1), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_333:
      coeffs%coeffs_333%c0 = derivative(czip)
      coeffs%coeffs_333%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_333%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_333%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_333%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_333%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_333%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_333%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_333%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_333%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_333%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_333%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_333%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_333%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_333%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_333:
   end subroutine reconstruct_d333
!---#] subroutine reconstruct_d333:
end module     p8_cbarc_hepemg_d333h0l1d
