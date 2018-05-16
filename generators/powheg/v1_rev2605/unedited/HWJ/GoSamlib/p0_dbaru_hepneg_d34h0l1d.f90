module     p0_dbaru_hepneg_d34h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p0_dbaru_hepneg/helicity0d34h0l1d.f90
   ! generator: buildfortran_d.py
   use p0_dbaru_hepneg_config, only: ki
   use p0_dbaru_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d34
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p0_dbaru_hepneg_model
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_color
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(28) :: acd34
      complex(ki) :: brack
      acd34(1)=dotproduct(k1,qshift)
      acd34(2)=abb34(17)
      acd34(3)=dotproduct(k2,qshift)
      acd34(4)=dotproduct(qshift,spvak5k4)
      acd34(5)=abb34(13)
      acd34(6)=abb34(15)
      acd34(7)=dotproduct(k6,qshift)
      acd34(8)=dotproduct(e6,qshift)
      acd34(9)=abb34(7)
      acd34(10)=abb34(6)
      acd34(11)=dotproduct(qshift,qshift)
      acd34(12)=abb34(18)
      acd34(13)=dotproduct(qshift,spvak1k6)
      acd34(14)=abb34(11)
      acd34(15)=dotproduct(qshift,spvae6k2)
      acd34(16)=abb34(16)
      acd34(17)=abb34(14)
      acd34(18)=dotproduct(qshift,spvak1k4)
      acd34(19)=abb34(8)
      acd34(20)=abb34(10)
      acd34(21)=abb34(5)
      acd34(22)=acd34(3)-acd34(7)
      acd34(23)=-acd34(5)*acd34(22)
      acd34(24)=acd34(13)*acd34(14)
      acd34(25)=acd34(15)*acd34(16)
      acd34(26)=acd34(8)*acd34(9)
      acd34(23)=acd34(26)+acd34(25)-acd34(17)+acd34(24)+acd34(23)
      acd34(23)=acd34(4)*acd34(23)
      acd34(22)=-acd34(6)*acd34(22)
      acd34(24)=-acd34(18)*acd34(19)
      acd34(25)=acd34(11)*acd34(12)
      acd34(26)=-acd34(1)*acd34(2)
      acd34(27)=-acd34(15)*acd34(20)
      acd34(28)=-acd34(8)*acd34(10)
      brack=acd34(21)+acd34(22)+acd34(23)+acd34(24)+acd34(25)+acd34(26)+acd34(2&
      &7)+acd34(28)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p0_dbaru_hepneg_model
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_color
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(34) :: acd34
      complex(ki) :: brack
      acd34(1)=k1(iv1)
      acd34(2)=abb34(17)
      acd34(3)=k2(iv1)
      acd34(4)=dotproduct(qshift,spvak5k4)
      acd34(5)=abb34(13)
      acd34(6)=abb34(15)
      acd34(7)=k6(iv1)
      acd34(8)=e6(iv1)
      acd34(9)=abb34(7)
      acd34(10)=abb34(6)
      acd34(11)=qshift(iv1)
      acd34(12)=abb34(18)
      acd34(13)=spvak5k4(iv1)
      acd34(14)=dotproduct(k2,qshift)
      acd34(15)=dotproduct(k6,qshift)
      acd34(16)=dotproduct(e6,qshift)
      acd34(17)=dotproduct(qshift,spvak1k6)
      acd34(18)=abb34(11)
      acd34(19)=dotproduct(qshift,spvae6k2)
      acd34(20)=abb34(16)
      acd34(21)=abb34(14)
      acd34(22)=spvak1k4(iv1)
      acd34(23)=abb34(8)
      acd34(24)=spvak1k6(iv1)
      acd34(25)=spvae6k2(iv1)
      acd34(26)=abb34(10)
      acd34(27)=-acd34(20)*acd34(25)
      acd34(28)=-acd34(18)*acd34(24)
      acd34(29)=-acd34(8)*acd34(9)
      acd34(30)=acd34(3)-acd34(7)
      acd34(31)=acd34(5)*acd34(30)
      acd34(27)=acd34(31)+acd34(29)+acd34(27)+acd34(28)
      acd34(27)=acd34(4)*acd34(27)
      acd34(28)=-acd34(20)*acd34(19)
      acd34(29)=-acd34(18)*acd34(17)
      acd34(31)=-acd34(9)*acd34(16)
      acd34(32)=-acd34(15)+acd34(14)
      acd34(32)=acd34(5)*acd34(32)
      acd34(28)=acd34(32)+acd34(31)+acd34(29)+acd34(21)+acd34(28)
      acd34(28)=acd34(13)*acd34(28)
      acd34(29)=acd34(6)*acd34(30)
      acd34(30)=acd34(22)*acd34(23)
      acd34(31)=acd34(11)*acd34(12)
      acd34(32)=acd34(1)*acd34(2)
      acd34(33)=acd34(25)*acd34(26)
      acd34(34)=acd34(8)*acd34(10)
      brack=acd34(27)+acd34(28)+acd34(29)+acd34(30)-2.0_ki*acd34(31)+acd34(32)+&
      &acd34(33)+acd34(34)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p0_dbaru_hepneg_model
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_color
      use p0_dbaru_hepneg_abbrevd34h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(23) :: acd34
      complex(ki) :: brack
      acd34(1)=d(iv1,iv2)
      acd34(2)=abb34(18)
      acd34(3)=k2(iv1)
      acd34(4)=spvak5k4(iv2)
      acd34(5)=abb34(13)
      acd34(6)=k2(iv2)
      acd34(7)=spvak5k4(iv1)
      acd34(8)=k6(iv1)
      acd34(9)=k6(iv2)
      acd34(10)=e6(iv1)
      acd34(11)=abb34(7)
      acd34(12)=e6(iv2)
      acd34(13)=spvak1k6(iv2)
      acd34(14)=abb34(11)
      acd34(15)=spvae6k2(iv2)
      acd34(16)=abb34(16)
      acd34(17)=spvak1k6(iv1)
      acd34(18)=spvae6k2(iv1)
      acd34(19)=acd34(16)*acd34(15)
      acd34(20)=acd34(14)*acd34(13)
      acd34(21)=acd34(11)*acd34(12)
      acd34(22)=acd34(9)-acd34(6)
      acd34(22)=acd34(5)*acd34(22)
      acd34(19)=acd34(22)+acd34(21)+acd34(19)+acd34(20)
      acd34(19)=acd34(7)*acd34(19)
      acd34(20)=acd34(16)*acd34(18)
      acd34(21)=acd34(14)*acd34(17)
      acd34(22)=acd34(11)*acd34(10)
      acd34(23)=acd34(8)-acd34(3)
      acd34(23)=acd34(5)*acd34(23)
      acd34(20)=acd34(23)+acd34(22)+acd34(20)+acd34(21)
      acd34(20)=acd34(4)*acd34(20)
      acd34(21)=acd34(1)*acd34(2)
      brack=acd34(19)+acd34(20)+2.0_ki*acd34(21)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p0_dbaru_hepneg_globalsl1, only: epspow
      use p0_dbaru_hepneg_kinematics
      use p0_dbaru_hepneg_abbrevd34h0
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
!---#[ subroutine reconstruct_d34:
   subroutine     reconstruct_d34(coeffs)
      use p0_dbaru_hepneg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_34:
      coeffs%coeffs_34%c0 = derivative(czip)
      coeffs%coeffs_34%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_34%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_34%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_34%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_34%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_34%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_34%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_34%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_34%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_34%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_34%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_34%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_34%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_34%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_34:
   end subroutine reconstruct_d34
!---#] subroutine reconstruct_d34:
end module     p0_dbaru_hepneg_d34h0l1d
