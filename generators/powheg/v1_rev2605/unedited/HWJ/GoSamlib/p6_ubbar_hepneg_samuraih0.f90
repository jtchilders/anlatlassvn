module     p6_ubbar_hepneg_samuraih0
   ! This file has been generated for samurai version 2.1.1
   ! Please, not the interface changes:
   ! 2.0 -> 2.1   : mu2 has changed from real to complex.
   ! 2.1 -> 2.1.1 : samurai_cm and samurai_rm have been made public
   !                we call them directly instead of the generic routine
   !                in order to avoid problems with some older versions of
   !                gfortran.
   !              + passing of invariants has been added.
   use precision, only: ki_sam => ki
   use p6_ubbar_hepneg_config, only: ki
   use p6_ubbar_hepneg_scalar_cache
   implicit none
   private
   public :: reduce_group0
   public :: reduce_group1
   public :: reduce_group2
   public :: reduce_group3
contains
!---#[ grouped numerators for samurai:
!-----#[ function numeval_group0:
function     numeval_group0(icut, Q, mu2) result(num)
   use p6_ubbar_hepneg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d70h0l1, only: numerator_d70 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:1-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case default
      R = real(k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom2 = Q2 - mT*mT
      R = real(-k3, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 70:
   if(nonzero(0)) then
      num = num + numerator_d70(icut, Q, mu2)
   end if
   !-------#] Diagram 70:
end function numeval_group0
!-----#] function numeval_group0:
!-----#[ function numeval_group1:
function     numeval_group1(icut, Q, mu2) result(num)
   use p6_ubbar_hepneg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d61h0l1, only: numerator_d61 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:1-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case default
      R = real(-k2, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      R = real(-k6, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(k3-k2+k5+k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 61:
   if(nonzero(0)) then
      num = num + numerator_d61(icut, Q, mu2)
   end if
   !-------#] Diagram 61:
end function numeval_group1
!-----#] function numeval_group1:
!-----#[ function numeval_group2:
function     numeval_group2(icut, Q, mu2) result(num)
   use p6_ubbar_hepneg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d34h0l1, only: numerator_d34 => numerator_samurai
   use p6_ubbar_hepneg_d40h0l1, only: numerator_d40 => numerator_samurai
   use p6_ubbar_hepneg_d60h0l1, only: numerator_d60 => numerator_samurai
   use p6_ubbar_hepneg_d92h0l1, only: numerator_d92 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:4-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case(1)
      nonzero(0) = .false.
      nonzero(1) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(10)
      nonzero(0) = .false.
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(21)
      nonzero(0) = .false.
      nonzero(1) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(210)
      nonzero(0) = .false.
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(3)
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(30)
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = Q2
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(31)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(310)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(32)
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(320)
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = Q2
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(321)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(3210)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case default
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 34:
   if(nonzero(0)) then
      num = num + numerator_d34(icut, Q, mu2) * denom2
   end if
   !-------#] Diagram 34:
   !-------#[ Diagram 40:
   if(nonzero(1)) then
      num = num + numerator_d40(icut, Q, mu2) * denom2 * denom4
   end if
   !-------#] Diagram 40:
   !-------#[ Diagram 60:
   if(nonzero(2)) then
      num = num + numerator_d60(icut, Q, mu2)
   end if
   !-------#] Diagram 60:
   !-------#[ Diagram 92:
   if(nonzero(3)) then
      num = num + numerator_d92(icut, Q, mu2) * denom4
   end if
   !-------#] Diagram 92:
end function numeval_group2
!-----#] function numeval_group2:
!-----#[ function numeval_group3:
function     numeval_group3(icut, Q, mu2) result(num)
   use p6_ubbar_hepneg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d30h0l1, only: numerator_d30 => numerator_samurai
   use p6_ubbar_hepneg_d37h0l1, only: numerator_d37 => numerator_samurai
   use p6_ubbar_hepneg_d59h0l1, only: numerator_d59 => numerator_samurai
   use p6_ubbar_hepneg_d94h0l1, only: numerator_d94 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:4-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case(1)
      nonzero(0) = .false.
      nonzero(1) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = Q2
   case(10)
      nonzero(0) = .false.
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = Q2
   case(21)
      nonzero(0) = .false.
      nonzero(1) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = Q2
   case(210)
      nonzero(0) = .false.
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = Q2
   case(3)
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(30)
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(31)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(310)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(32)
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(320)
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(321)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(3210)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(3) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case default
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k3-k6-k5-k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = Q2
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 30:
   if(nonzero(0)) then
      num = num + numerator_d30(icut, Q, mu2) * denom2
   end if
   !-------#] Diagram 30:
   !-------#[ Diagram 37:
   if(nonzero(1)) then
      num = num + numerator_d37(icut, Q, mu2) * denom2 * denom4
   end if
   !-------#] Diagram 37:
   !-------#[ Diagram 59:
   if(nonzero(2)) then
      num = num + numerator_d59(icut, Q, mu2)
   end if
   !-------#] Diagram 59:
   !-------#[ Diagram 94:
   if(nonzero(3)) then
      num = num + numerator_d94(icut, Q, mu2) * denom4
   end if
   !-------#] Diagram 94:
end function numeval_group3
!-----#] function numeval_group3:
!---#] grouped numerators for samurai:
!---#[ reduce groups with samurai:
!-----#[ subroutine reduce_group0:
subroutine     reduce_group0(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p6_ubbar_hepneg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p6_ubbar_hepneg_kinematics
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d70h0l1, only: numerator_diagram70 => numerator_samurai
   use p6_ubbar_hepneg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 3
   !-----------#[ invariants for samurai:
   complex(ki_sam), dimension(3, 3) :: g_mat
   !-----------#] initialize invariants:
   real(ki_sam), dimension(3) :: msq
   real(ki_sam), dimension(3,4) :: Vi

   if(samurai_test.eq.1 .or. samurai_test.eq.3) then
      istopm = 1
      istop0 = 1
   else
      istopm = samurai_istop
      istop0 = max(2,samurai_istop)
   end if
   msq(1) = real(mT*mT, ki_sam)
   Vi(1,:) = real(k6((/2,3,4,1/)), ki_sam)
   msq(2) = real(mT*mT, ki_sam)
   Vi(2,:) = real(0, ki_sam)
   msq(3) = real(mT*mT, ki_sam)
   Vi(3,:) = real(-k3((/2,3,4,1/)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(1, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(-es123+mH**2-2.0_ki*mT**2+es45-es345+es12, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(2, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 3) = real(mH**2-2.0_ki*mT**2, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(3, 3) = real(-2.0_ki*mT**2, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group0:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group0"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(3, 3))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group0, tot, totr, Vi, msq, 3, &
         & effective_group_rank, istopm, scale2, ok, &
         & samurai_cache_flag_g0, samurai_cache_g0)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group0:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='70'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram70"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,3/), (/1,2,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram70, acc, accr, &
            & Vi((/1,2,3/),:), msq((/1,2,3/)), 3, &
            & 3, istopm, scale2, ok, &
            & samurai_cache_flag_d70, samurai_cache_d70)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", -real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", -real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", -real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", -real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot =  + acc
         totr =  + accr
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group0
!-----#] subroutine reduce_group0:
!-----#[ subroutine reduce_group1:
subroutine     reduce_group1(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p6_ubbar_hepneg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p6_ubbar_hepneg_kinematics
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d61h0l1, only: numerator_diagram61 => numerator_samurai
   use p6_ubbar_hepneg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 3
   !-----------#[ invariants for samurai:
   complex(ki_sam), dimension(4, 4) :: g_mat
   !-----------#] initialize invariants:
   real(ki_sam), dimension(4) :: msq
   real(ki_sam), dimension(4,4) :: Vi

   if(samurai_test.eq.1 .or. samurai_test.eq.3) then
      istopm = 1
      istop0 = 1
   else
      istopm = samurai_istop
      istop0 = max(2,samurai_istop)
   end if
   msq(1) = 0.0_ki_sam
   Vi(1,:) = real(-k2((/2,3,4,1/)), ki_sam)
   msq(2) = 0.0_ki_sam
   Vi(2,:) = real(0, ki_sam)
   msq(3) = 0.0_ki_sam
   Vi(3,:) = real(-k6((/2,3,4,1/)), ki_sam)
   msq(4) = 0.0_ki_sam
   Vi(4,:) = real(k3((/2,3,4,1/))-k2((/2,3,4,1/))+k5((/2,3,4,1/))+k4((/2,3,4,1/&
   &)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(0.0_ki, ki_sam)
   g_mat(1, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(-es12-es61+es345, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(es345, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(es61, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 4) = real(0.0_ki, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(0.0_ki, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group1:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group1"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group1, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istop0, scale2, ok, &
         & samurai_cache_flag_g1, samurai_cache_g1)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group1:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='61'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram61"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram61, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, ok, &
            & samurai_cache_flag_d61, samurai_cache_d61)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot =  + acc
         totr =  + accr
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group1
!-----#] subroutine reduce_group1:
!-----#[ subroutine reduce_group2:
subroutine     reduce_group2(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p6_ubbar_hepneg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p6_ubbar_hepneg_kinematics
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d34h0l1, only: numerator_diagram34 => numerator_samurai
   use p6_ubbar_hepneg_d40h0l1, only: numerator_diagram40 => numerator_samurai
   use p6_ubbar_hepneg_d60h0l1, only: numerator_diagram60 => numerator_samurai
   use p6_ubbar_hepneg_d92h0l1, only: numerator_diagram92 => numerator_samurai
   use p6_ubbar_hepneg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 3
   !-----------#[ invariants for samurai:
   complex(ki_sam), dimension(4, 4) :: g_mat
   !-----------#] initialize invariants:
   real(ki_sam), dimension(4) :: msq
   real(ki_sam), dimension(4,4) :: Vi

   if(samurai_test.eq.1 .or. samurai_test.eq.3) then
      istopm = 1
      istop0 = 1
   else
      istopm = samurai_istop
      istop0 = max(2,samurai_istop)
   end if
   msq(1) = 0.0_ki_sam
   Vi(1,:) = real(-k6((/2,3,4,1/)), ki_sam)
   msq(2) = 0.0_ki_sam
   Vi(2,:) = real(0, ki_sam)
   msq(3) = 0.0_ki_sam
   Vi(3,:) = real(-k2((/2,3,4,1/)), ki_sam)
   msq(4) = 0.0_ki_sam
   Vi(4,:) = real(-k3((/2,3,4,1/))-k6((/2,3,4,1/))-k5((/2,3,4,1/))-k4((/2,3,4,1&
   &/)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(0.0_ki, ki_sam)
   g_mat(1, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(-es12-es61+es345, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(es345, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(es12, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 4) = real(0.0_ki, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(0.0_ki, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group2:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group2"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group2, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istop0, scale2, ok, &
         & samurai_cache_flag_g2, samurai_cache_g2)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group2:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='34'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram34"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,3,4/), (/1,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram34, acc, accr, &
            & Vi((/1,3,4/),:), msq((/1,3,4/)), 3, &
            & 2, istop0, scale2, ok, &
            & samurai_cache_flag_d34, samurai_cache_d34)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot =  + acc
         totr =  + accr
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='40'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram40"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(2, 2))
         s_mat(:,:) = g_mat( (/1,3/), (/1,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram40, acc, accr, &
            & Vi((/1,3/),:), msq((/1,3/)), 2, &
            & 1, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d40, samurai_cache_d40)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='60'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram60"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram60, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d60, samurai_cache_d60)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='92'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram92"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,3/), (/1,2,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram92, acc, accr, &
            & Vi((/1,2,3/),:), msq((/1,2,3/)), 3, &
            & 2, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d92, samurai_cache_d92)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group2
!-----#] subroutine reduce_group2:
!-----#[ subroutine reduce_group3:
subroutine     reduce_group3(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p6_ubbar_hepneg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p6_ubbar_hepneg_kinematics
   use p6_ubbar_hepneg_model
   use p6_ubbar_hepneg_d30h0l1, only: numerator_diagram30 => numerator_samurai
   use p6_ubbar_hepneg_d37h0l1, only: numerator_diagram37 => numerator_samurai
   use p6_ubbar_hepneg_d59h0l1, only: numerator_diagram59 => numerator_samurai
   use p6_ubbar_hepneg_d94h0l1, only: numerator_diagram94 => numerator_samurai
   use p6_ubbar_hepneg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 3
   !-----------#[ invariants for samurai:
   complex(ki_sam), dimension(4, 4) :: g_mat
   !-----------#] initialize invariants:
   real(ki_sam), dimension(4) :: msq
   real(ki_sam), dimension(4,4) :: Vi

   if(samurai_test.eq.1 .or. samurai_test.eq.3) then
      istopm = 1
      istop0 = 1
   else
      istopm = samurai_istop
      istop0 = max(2,samurai_istop)
   end if
   msq(1) = 0.0_ki_sam
   Vi(1,:) = real(-k3((/2,3,4,1/))-k5((/2,3,4,1/))-k4((/2,3,4,1/)), ki_sam)
   msq(2) = 0.0_ki_sam
   Vi(2,:) = real(-k3((/2,3,4,1/))-k6((/2,3,4,1/))-k5((/2,3,4,1/))-k4((/2,3,4,1&
   &/)), ki_sam)
   msq(3) = 0.0_ki_sam
   Vi(3,:) = real(-k2((/2,3,4,1/)), ki_sam)
   msq(4) = 0.0_ki_sam
   Vi(4,:) = real(0, ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(0.0_ki, ki_sam)
   g_mat(1, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(es61, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(es345, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(0.0_ki, ki_sam)
   g_mat(2, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(es12, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(0.0_ki, ki_sam)
   g_mat(3, 4) = real(0.0_ki, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(0.0_ki, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group3:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group3"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group3, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istop0, scale2, ok, &
         & samurai_cache_flag_g3, samurai_cache_g3)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group3:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='30'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram30"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,3,4/), (/1,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram30, acc, accr, &
            & Vi((/1,3,4/),:), msq((/1,3,4/)), 3, &
            & 2, istop0, scale2, ok, &
            & samurai_cache_flag_d30, samurai_cache_d30)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot =  + acc
         totr =  + accr
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='37'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram37"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(2, 2))
         s_mat(:,:) = g_mat( (/1,3/), (/1,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram37, acc, accr, &
            & Vi((/1,3/),:), msq((/1,3/)), 2, &
            & 1, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d37, samurai_cache_d37)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='59'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram59"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram59, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d59, samurai_cache_d59)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='94'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram94"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,3/), (/1,2,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram94, acc, accr, &
            & Vi((/1,2,3/),:), msq((/1,2,3/)), 3, &
            & 2, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d94, samurai_cache_d94)
         !-----------#[ deallocate invariants:
         deallocate(s_mat)
         !-----------#] deallocate invariants:
         if(debug_nlo_diagrams) then
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-finite' re='", +real(acc(0), ki), &
               & "' im='", aimag(acc(0)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-single' re='", +real(acc(-1), ki), &
               & "' im='", aimag(acc(-1)), "'/>"
            write(logfile,'(A30,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-double' re='", +real(acc(-2), ki), &
               & "' im='", aimag(acc(-2)), "'/>"
            write(logfile,'(A32,E24.16,A6,E24.16,A3)') &
               & "<result kind='nlo-rational' re='", +real(accr, ki), &
               & "' im='", aimag(accr), "'/>"
            if(ok) then
               write(logfile,'(A30)') "<flag name='ok' status='yes'/>"
            else
               write(logfile,'(A29)') "<flag name='ok' status='no'/>"
            end if
            write(logfile,*) "</diagram>"
         end if

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group3
!-----#] subroutine reduce_group3:
!---#] reduce groups with samurai:
end module p6_ubbar_hepneg_samuraih0
