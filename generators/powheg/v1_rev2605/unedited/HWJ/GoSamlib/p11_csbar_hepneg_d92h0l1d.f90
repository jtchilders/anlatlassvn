module     p11_csbar_hepneg_d92h0l1d
   ! file: /home/gionata/Documenti/Lavoro/GoSamPowheg/POWHEG-BOX/HWJ_tmp/GoSam_ &
   ! &POWHEG/Virtual/p11_csbar_hepneg/helicity0d92h0l1d.f90
   ! generator: buildfortran_d.py
   use p11_csbar_hepneg_config, only: ki
   use p11_csbar_hepneg_util, only: cond, d => metric_tensor
   implicit none
   private
   complex(ki), parameter :: i_ = (0.0_ki, 1.0_ki)
   integer, private :: iv0
   integer, private :: iv1
   integer, private :: iv2
   real(ki), dimension(4), private :: qshift
   public :: derivative , reconstruct_d92
contains
!---#[ function brack_1:
   pure function brack_1(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd92h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(1) :: acd92
      complex(ki) :: brack
      acd92(1)=abb92(14)
      brack=acd92(1)
   end function brack_1
!---#] function brack_1:
!---#[ function brack_2:
   pure function brack_2(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd92h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(18) :: acd92
      complex(ki) :: brack
      acd92(1)=k2(iv1)
      acd92(2)=abb92(9)
      acd92(3)=k6(iv1)
      acd92(4)=abb92(16)
      acd92(5)=e6(iv1)
      acd92(6)=abb92(13)
      acd92(7)=spvak5k1(iv1)
      acd92(8)=abb92(7)
      acd92(9)=spvak5k4(iv1)
      acd92(10)=abb92(10)
      acd92(11)=spvak2e6(iv1)
      acd92(12)=abb92(6)
      acd92(13)=acd92(2)*acd92(1)
      acd92(14)=acd92(4)*acd92(3)
      acd92(15)=acd92(6)*acd92(5)
      acd92(16)=acd92(8)*acd92(7)
      acd92(17)=acd92(10)*acd92(9)
      acd92(18)=acd92(12)*acd92(11)
      brack=acd92(13)+acd92(14)+acd92(15)+acd92(16)+acd92(17)+acd92(18)
   end function brack_2
!---#] function brack_2:
!---#[ function brack_3:
   pure function brack_3(Q, mu2) result(brack)
      use p11_csbar_hepneg_model
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_color
      use p11_csbar_hepneg_abbrevd92h0
      implicit none
      complex(ki), dimension(4), intent(in) :: Q
      complex(ki), intent(in) :: mu2
      complex(ki), dimension(25) :: acd92
      complex(ki) :: brack
      acd92(1)=d(iv1,iv2)
      acd92(2)=abb92(15)
      acd92(3)=k2(iv1)
      acd92(4)=e6(iv2)
      acd92(5)=abb92(11)
      acd92(6)=k2(iv2)
      acd92(7)=e6(iv1)
      acd92(8)=k6(iv1)
      acd92(9)=abb92(17)
      acd92(10)=k6(iv2)
      acd92(11)=spvak2k6(iv2)
      acd92(12)=abb92(8)
      acd92(13)=spvak5k1(iv2)
      acd92(14)=abb92(5)
      acd92(15)=spvak5k4(iv2)
      acd92(16)=abb92(18)
      acd92(17)=spvak2k6(iv1)
      acd92(18)=spvak5k1(iv1)
      acd92(19)=spvak5k4(iv1)
      acd92(20)=acd92(16)*acd92(15)
      acd92(21)=acd92(14)*acd92(13)
      acd92(22)=acd92(12)*acd92(11)
      acd92(23)=acd92(9)*acd92(10)
      acd92(24)=acd92(5)*acd92(6)
      acd92(20)=acd92(24)+acd92(23)+acd92(22)+acd92(20)+acd92(21)
      acd92(20)=acd92(7)*acd92(20)
      acd92(21)=acd92(16)*acd92(19)
      acd92(22)=acd92(14)*acd92(18)
      acd92(23)=acd92(12)*acd92(17)
      acd92(24)=acd92(9)*acd92(8)
      acd92(25)=acd92(5)*acd92(3)
      acd92(21)=acd92(25)+acd92(24)+acd92(23)+acd92(21)+acd92(22)
      acd92(21)=acd92(4)*acd92(21)
      acd92(22)=acd92(1)*acd92(2)
      brack=acd92(20)+acd92(21)+2.0_ki*acd92(22)
   end function brack_3
!---#] function brack_3:
!---#[ function derivative:
   function derivative(mu2,i1,i2) result(numerator)
      use p11_csbar_hepneg_globalsl1, only: epspow
      use p11_csbar_hepneg_kinematics
      use p11_csbar_hepneg_abbrevd92h0
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
      qshift = 0
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
!---#[ subroutine reconstruct_d92:
   subroutine     reconstruct_d92(coeffs)
      use p11_csbar_hepneg_groups, only: tensrec_info_group2
      implicit none
      complex(ki), parameter :: czip = (0.0_ki, 0.0_ki)
      complex(ki), parameter :: cone = (1.0_ki, 0.0_ki)
      type(tensrec_info_group2), intent(out) :: coeffs
      ! rank 2 case :
      !---[# reconstruct coeffs%coeffs_92:
      coeffs%coeffs_92%c0 = derivative(czip)
      coeffs%coeffs_92%c1(1,1) = derivative(czip,1)
      coeffs%coeffs_92%c1(1,2) = derivative(czip,1,1)/ 2.0_ki
      coeffs%coeffs_92%c1(2,1) = -derivative(czip,2)
      coeffs%coeffs_92%c1(2,2) = derivative(czip,2,2)/ 2.0_ki
      coeffs%coeffs_92%c1(3,1) = -derivative(czip,3)
      coeffs%coeffs_92%c1(3,2) = derivative(czip,3,3)/ 2.0_ki
      coeffs%coeffs_92%c1(4,1) = -derivative(czip,4)
      coeffs%coeffs_92%c1(4,2) = derivative(czip,4,4)/ 2.0_ki
      coeffs%coeffs_92%c2(1,1) = -derivative(czip,1,2)
      coeffs%coeffs_92%c2(2,1) = -derivative(czip,1,3)
      coeffs%coeffs_92%c2(3,1) = -derivative(czip,1,4)
      coeffs%coeffs_92%c2(4,1) = derivative(czip,2,3)
      coeffs%coeffs_92%c2(5,1) = derivative(czip,2,4)
      coeffs%coeffs_92%c2(6,1) = derivative(czip,3,4)
      !---#] reconstruct coeffs%coeffs_92:
   end subroutine reconstruct_d92
!---#] subroutine reconstruct_d92:
end module     p11_csbar_hepneg_d92h0l1d
