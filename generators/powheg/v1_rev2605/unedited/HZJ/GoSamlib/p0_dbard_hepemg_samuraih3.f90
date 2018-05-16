module     p0_dbard_hepemg_samuraih3
   ! This file has been generated for samurai version 2.1.1
   ! Please, not the interface changes:
   ! 2.0 -> 2.1   : mu2 has changed from real to complex.
   ! 2.1 -> 2.1.1 : samurai_cm and samurai_rm have been made public
   !                we call them directly instead of the generic routine
   !                in order to avoid problems with some older versions of
   !                gfortran.
   !              + passing of invariants has been added.
   use precision, only: ki_sam => ki
   use p0_dbard_hepemg_config, only: ki
   use p0_dbard_hepemg_scalar_cache
   implicit none
   private
   public :: reduce_group0
   public :: reduce_group1
   public :: reduce_group2
   public :: reduce_group3
   public :: reduce_group4
   public :: reduce_group5
contains
!---#[ grouped numerators for samurai:
!-----#[ function numeval_group0:
function     numeval_group0(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d541h3l1, only: numerator_d541 => numerator_samurai
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
   !-------#[ Diagram 541:
   if(nonzero(0)) then
      num = num + numerator_d541(icut, Q, mu2)
   end if
   !-------#] Diagram 541:
end function numeval_group0
!-----#] function numeval_group0:
!-----#[ function numeval_group1:
function     numeval_group1(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d279h3l1, only: numerator_d279 => numerator_samurai
   use p0_dbard_hepemg_d333h3l1, only: numerator_d333 => numerator_samurai
   use p0_dbard_hepemg_d421h3l1, only: numerator_d421 => numerator_samurai
   use p0_dbard_hepemg_d537h3l1, only: numerator_d537 => numerator_samurai
   use p0_dbard_hepemg_d605h3l1, only: numerator_d605 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:5-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case(1)
      nonzero(1) = .false.
      nonzero(2) = .false.
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
      nonzero(1) = .false.
      nonzero(2) = .false.
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
   case(2)
      nonzero(0) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      denom3 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(20)
      nonzero(0) = .false.
      denom1 = 0.0_ki
      denom2 = Q2
      denom3 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(21)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(2) = .false.
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
      nonzero(2) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(-k3-k6-k5-k4, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
   case(3)
      nonzero(2) = .false.
      nonzero(4) = .false.
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
      nonzero(2) = .false.
      nonzero(4) = .false.
      denom1 = 0.0_ki
      denom2 = Q2
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(31)
      nonzero(1) = .false.
      nonzero(2) = .false.
      nonzero(4) = .false.
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
      nonzero(1) = .false.
      nonzero(2) = .false.
      nonzero(4) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      R = real(-k2, ki_sam)
      denom3 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom4 = 0.0_ki
   case(32)
      nonzero(0) = .false.
      nonzero(2) = .false.
      nonzero(4) = .false.
      R = real(-k6, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)
      denom2 = Q2
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(320)
      nonzero(0) = .false.
      nonzero(2) = .false.
      nonzero(4) = .false.
      denom1 = 0.0_ki
      denom2 = Q2
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(321)
      nonzero(0) = .false.
      nonzero(1) = .false.
      nonzero(2) = .false.
      nonzero(4) = .false.
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
      nonzero(2) = .false.
      nonzero(4) = .false.
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
   !-------#[ Diagram 279:
   if(nonzero(0)) then
      num = num + numerator_d279(icut, Q, mu2) * denom3
   end if
   !-------#] Diagram 279:
   !-------#[ Diagram 333:
   if(nonzero(1)) then
      num = num + numerator_d333(icut, Q, mu2) * denom2
   end if
   !-------#] Diagram 333:
   !-------#[ Diagram 421:
   if(nonzero(2)) then
      num = num + numerator_d421(icut, Q, mu2) * denom2 * denom4
   end if
   !-------#] Diagram 421:
   !-------#[ Diagram 537:
   if(nonzero(3)) then
      num = num + numerator_d537(icut, Q, mu2)
   end if
   !-------#] Diagram 537:
   !-------#[ Diagram 605:
   if(nonzero(4)) then
      num = num + numerator_d605(icut, Q, mu2) * denom4
   end if
   !-------#] Diagram 605:
end function numeval_group1
!-----#] function numeval_group1:
!-----#[ function numeval_group2:
function     numeval_group2(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d323h3l1, only: numerator_d323 => numerator_samurai
   use p0_dbard_hepemg_d413h3l1, only: numerator_d413 => numerator_samurai
   use p0_dbard_hepemg_d533h3l1, only: numerator_d533 => numerator_samurai
   use p0_dbard_hepemg_d613h3l1, only: numerator_d613 => numerator_samurai
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
   !-------#[ Diagram 323:
   if(nonzero(0)) then
      num = num + numerator_d323(icut, Q, mu2) * denom2
   end if
   !-------#] Diagram 323:
   !-------#[ Diagram 413:
   if(nonzero(1)) then
      num = num + numerator_d413(icut, Q, mu2) * denom2 * denom4
   end if
   !-------#] Diagram 413:
   !-------#[ Diagram 533:
   if(nonzero(2)) then
      num = num + numerator_d533(icut, Q, mu2)
   end if
   !-------#] Diagram 533:
   !-------#[ Diagram 613:
   if(nonzero(3)) then
      num = num + numerator_d613(icut, Q, mu2) * denom4
   end if
   !-------#] Diagram 613:
end function numeval_group2
!-----#] function numeval_group2:
!-----#[ function numeval_group3:
function     numeval_group3(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d101h3l1, only: numerator_d101 => numerator_samurai
   use p0_dbard_hepemg_d553h3l1, only: numerator_d553 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:2-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case(0)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(-k3, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = Q2 - mT*mT
      R = real(k6, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(10)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = Q2 - mT*mT
      R = real(k6, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(20)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(-k3, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      R = real(k6, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(210)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(k6, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(30)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(-k3, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = Q2 - mT*mT
      denom4 = 0.0_ki
   case(310)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = Q2 - mT*mT
      denom4 = 0.0_ki
   case(320)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(-k3, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(3210)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case default
      R = real(-k3-k5-k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      R = real(-k3, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = Q2 - mT*mT
      R = real(k6, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 101:
   if(nonzero(0)) then
      num = num + numerator_d101(icut, Q, mu2)
   end if
   !-------#] Diagram 101:
   !-------#[ Diagram 553:
   if(nonzero(1)) then
      num = num + numerator_d553(icut, Q, mu2) * denom1
   end if
   !-------#] Diagram 553:
end function numeval_group3
!-----#] function numeval_group3:
!-----#[ function numeval_group4:
function     numeval_group4(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d77h3l1, only: numerator_d77 => numerator_samurai
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
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      R = real(k6, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = Q2 - mT*mT
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 77:
   if(nonzero(0)) then
      num = num + numerator_d77(icut, Q, mu2)
   end if
   !-------#] Diagram 77:
end function numeval_group4
!-----#] function numeval_group4:
!-----#[ function numeval_group5:
function     numeval_group5(icut, Q, mu2) result(num)
   use p0_dbard_hepemg_kinematics, only: k1, k2, k3, k4, k5, k6
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d45h3l1, only: numerator_d45 => numerator_samurai
   use p0_dbard_hepemg_d277h3l1, only: numerator_d277 => numerator_samurai
   implicit none
   integer, intent(in) :: icut
   complex(ki_sam), dimension(4), intent(in) :: Q
   complex(ki_sam), intent(in) :: mu2
   complex(ki_sam) :: num

   logical, dimension(0:2-1) :: nonzero
   real(ki_sam), dimension(0:3) :: R
   complex(ki_sam) :: Q2
   complex(ki_sam) ::denom1,denom2,denom3,denom4

   nonzero(:) = .true.
   Q2 = Q(4)*Q(4) - Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) - mu2

   select case(icut)
   case(2)
      nonzero(1) = .false.
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      R = real(k5+k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(20)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(k5+k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(21)
      nonzero(1) = .false.
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(210)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   case(32)
      nonzero(1) = .false.
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      R = real(k5+k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(320)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      R = real(k5+k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(321)
      nonzero(1) = .false.
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case(3210)
      nonzero(1) = .false.
      denom1 = 0.0_ki
      denom2 = 0.0_ki
      denom3 = 0.0_ki
      denom4 = 0.0_ki
   case default
      R = real(k6+k5+k4, ki_sam)
      denom1 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      R = real(k5+k4, ki_sam)
      denom2 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
      denom3 = Q2 - mT*mT
      R = real(-k3, ki_sam)
      denom4 = Q2 + (Q(4) + Q(4) + R(0))*R(0) &
                 &    - (Q(1) + Q(1) + R(1))*R(1) &
                 &    - (Q(2) + Q(2) + R(2))*R(2) &
                 &    - (Q(3) + Q(3) + R(3))*R(3)&
                 &    - mT*mT
   end select

   num = (0.0_ki_sam, 0.0_ki_sam)
   !-------#[ Diagram 45:
   if(nonzero(0)) then
      num = num + numerator_d45(icut, Q, mu2)
   end if
   !-------#] Diagram 45:
   !-------#[ Diagram 277:
   if(nonzero(1)) then
      num = num + numerator_d277(icut, Q, mu2) * denom3
   end if
   !-------#] Diagram 277:
end function numeval_group5
!-----#] function numeval_group5:
!---#] grouped numerators for samurai:
!---#[ reduce groups with samurai:
!-----#[ subroutine reduce_group0:
subroutine     reduce_group0(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d541h3l1, only: numerator_diagram541 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

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
      !------#[ reduce numerator numeval_group0:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group0"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group0, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istop0, scale2, ok, &
         & samurai_cache_flag_g0, samurai_cache_g0)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group0:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='541'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram541"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram541, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, ok, &
            & samurai_cache_flag_d541, samurai_cache_d541)
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
end subroutine reduce_group0
!-----#] subroutine reduce_group0:
!-----#[ subroutine reduce_group1:
subroutine     reduce_group1(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d279h3l1, only: numerator_diagram279 => numerator_samurai
   use p0_dbard_hepemg_d333h3l1, only: numerator_diagram333 => numerator_samurai
   use p0_dbard_hepemg_d421h3l1, only: numerator_diagram421 => numerator_samurai
   use p0_dbard_hepemg_d537h3l1, only: numerator_diagram537 => numerator_samurai
   use p0_dbard_hepemg_d605h3l1, only: numerator_diagram605 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4
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
            write(logfile,*) "<diagram index='279'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram279"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,4/), (/1,2,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram279, acc, accr, &
            & Vi((/1,2,4/),:), msq((/1,2,4/)), 3, &
            & 3, istop0, scale2, ok, &
            & samurai_cache_flag_d279, samurai_cache_d279)
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

         tot =  +  acc
         totr =  +  accr
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='333'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram333"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,3,4/), (/1,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram333, acc, accr, &
            & Vi((/1,3,4/),:), msq((/1,3,4/)), 3, &
            & 2, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d333, samurai_cache_d333)
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
            write(logfile,*) "<diagram index='421'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram421"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(2, 2))
         s_mat(:,:) = g_mat( (/1,3/), (/1,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram421, acc, accr, &
            & Vi((/1,3/),:), msq((/1,3/)), 2, &
            & 1, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d421, samurai_cache_d421)
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
            write(logfile,*) "<diagram index='537'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram537"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram537, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d537, samurai_cache_d537)
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
            write(logfile,*) "<diagram index='605'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram605"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,3/), (/1,2,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram605, acc, accr, &
            & Vi((/1,2,3/),:), msq((/1,2,3/)), 3, &
            & 2, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d605, samurai_cache_d605)
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
end subroutine reduce_group1
!-----#] subroutine reduce_group1:
!-----#[ subroutine reduce_group2:
subroutine     reduce_group2(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d323h3l1, only: numerator_diagram323 => numerator_samurai
   use p0_dbard_hepemg_d413h3l1, only: numerator_diagram413 => numerator_samurai
   use p0_dbard_hepemg_d533h3l1, only: numerator_diagram533 => numerator_samurai
   use p0_dbard_hepemg_d613h3l1, only: numerator_diagram613 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

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
            write(logfile,*) "<diagram index='323'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram323"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,3,4/), (/1,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram323, acc, accr, &
            & Vi((/1,3,4/),:), msq((/1,3,4/)), 3, &
            & 2, istop0, scale2, ok, &
            & samurai_cache_flag_d323, samurai_cache_d323)
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
            write(logfile,*) "<diagram index='413'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram413"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(2, 2))
         s_mat(:,:) = g_mat( (/1,3/), (/1,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram413, acc, accr, &
            & Vi((/1,3/),:), msq((/1,3/)), 2, &
            & 1, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d413, samurai_cache_d413)
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
            write(logfile,*) "<diagram index='533'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram533"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram533, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 3, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d533, samurai_cache_d533)
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
            write(logfile,*) "<diagram index='613'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram613"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,3/), (/1,2,3/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram613, acc, accr, &
            & Vi((/1,2,3/),:), msq((/1,2,3/)), 3, &
            & 2, istop0, scale2, acc_ok, &
            & samurai_cache_flag_d613, samurai_cache_d613)
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
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d101h3l1, only: numerator_diagram101 => numerator_samurai
   use p0_dbard_hepemg_d553h3l1, only: numerator_diagram553 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4
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
   msq(1) = real(mT*mT, ki_sam)
   Vi(1,:) = real(-k3((/2,3,4,1/))-k5((/2,3,4,1/))-k4((/2,3,4,1/)), ki_sam)
   msq(2) = real(mT*mT, ki_sam)
   Vi(2,:) = real(-k3((/2,3,4,1/)), ki_sam)
   msq(3) = real(mT*mT, ki_sam)
   Vi(3,:) = real(0, ki_sam)
   msq(4) = real(mT*mT, ki_sam)
   Vi(4,:) = real(k6((/2,3,4,1/)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(1, 2) = real(es45-2.0_ki*mT**2, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(es345-2.0_ki*mT**2, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(-2.0_ki*mT**2+es12, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 3) = real(mH**2-2.0_ki*mT**2, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(-es123+mH**2-2.0_ki*mT**2+es45-es345+es12, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(3, 4) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(-2.0_ki*mT**2, ki_sam)
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
         & effective_group_rank, istopm, scale2, ok, &
         & samurai_cache_flag_g3, samurai_cache_g3)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group3:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='101'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram101"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram101, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 4, istopm, scale2, ok, &
            & samurai_cache_flag_d101, samurai_cache_d101)
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='553'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram553"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/2,3,4/), (/2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram553, acc, accr, &
            & Vi((/2,3,4/),:), msq((/2,3,4/)), 3, &
            & 3, istopm, scale2, acc_ok, &
            & samurai_cache_flag_d553, samurai_cache_d553)
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

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group3
!-----#] subroutine reduce_group3:
!-----#[ subroutine reduce_group4:
subroutine     reduce_group4(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d77h3l1, only: numerator_diagram77 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4
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
   msq(1) = real(mT*mT, ki_sam)
   Vi(1,:) = real(k6((/2,3,4,1/))+k5((/2,3,4,1/))+k4((/2,3,4,1/)), ki_sam)
   msq(2) = real(mT*mT, ki_sam)
   Vi(2,:) = real(k6((/2,3,4,1/)), ki_sam)
   msq(3) = real(mT*mT, ki_sam)
   Vi(3,:) = real(0, ki_sam)
   msq(4) = real(mT*mT, ki_sam)
   Vi(4,:) = real(-k3((/2,3,4,1/)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(1, 2) = real(es45-2.0_ki*mT**2, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(-2.0_ki*mT**2+es123, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(-2.0_ki*mT**2+es12, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 3) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(-es123+mH**2-2.0_ki*mT**2+es45-es345+es12, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(3, 4) = real(mH**2-2.0_ki*mT**2, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(-2.0_ki*mT**2, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group4:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group4"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group4, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istopm, scale2, ok, &
         & samurai_cache_flag_g4, samurai_cache_g4)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group4:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='77'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram77"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram77, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 4, istopm, scale2, ok, &
            & samurai_cache_flag_d77, samurai_cache_d77)
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
end subroutine reduce_group4
!-----#] subroutine reduce_group4:
!-----#[ subroutine reduce_group5:
subroutine     reduce_group5(scale2,tot,totr,ok)
   use msamurai, only: samurai, samurai_rm, samurai_cm
   use options, only: samurai_out => iout
   use madds, only: s_mat
   use p0_dbard_hepemg_config, only: samurai_group_numerators, &
      & samurai_verbosity, samurai_istop, samurai_test, &
      & debug_nlo_diagrams, logfile
   use p0_dbard_hepemg_kinematics
   use p0_dbard_hepemg_model
   use p0_dbard_hepemg_d45h3l1, only: numerator_diagram45 => numerator_samurai
   use p0_dbard_hepemg_d277h3l1, only: numerator_diagram277 => numerator_samurai
   use p0_dbard_hepemg_globalsl1, only: epspow

   implicit none
   real(ki_sam), intent(in) :: scale2
   complex(ki_sam), dimension(-2:0), intent(out) :: tot
   complex(ki_sam), intent(out) :: totr
   logical, intent(out) :: ok

   complex(ki_sam), dimension(-2:0) :: acc
   complex(ki_sam) :: accr
   logical :: acc_ok

   integer :: istopm, istop0

   integer, parameter :: effective_group_rank = 4
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
   msq(1) = real(mT*mT, ki_sam)
   Vi(1,:) = real(k6((/2,3,4,1/))+k5((/2,3,4,1/))+k4((/2,3,4,1/)), ki_sam)
   msq(2) = real(mT*mT, ki_sam)
   Vi(2,:) = real(k5((/2,3,4,1/))+k4((/2,3,4,1/)), ki_sam)
   msq(3) = real(mT*mT, ki_sam)
   Vi(3,:) = real(0, ki_sam)
   msq(4) = real(mT*mT, ki_sam)
   Vi(4,:) = real(-k3((/2,3,4,1/)), ki_sam)
   !-----------#[ initialize invariants:
   g_mat(1, 1) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(1, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 1) = g_mat(1, 2)
   g_mat(1, 3) = real(-2.0_ki*mT**2+es123, ki_sam)
   g_mat(3, 1) = g_mat(1, 3)
   g_mat(1, 4) = real(-2.0_ki*mT**2+es12, ki_sam)
   g_mat(4, 1) = g_mat(1, 4)
   g_mat(2, 2) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(2, 3) = real(es45-2.0_ki*mT**2, ki_sam)
   g_mat(3, 2) = g_mat(2, 3)
   g_mat(2, 4) = real(es345-2.0_ki*mT**2, ki_sam)
   g_mat(4, 2) = g_mat(2, 4)
   g_mat(3, 3) = real(-2.0_ki*mT**2, ki_sam)
   g_mat(3, 4) = real(mH**2-2.0_ki*mT**2, ki_sam)
   g_mat(4, 3) = g_mat(3, 4)
   g_mat(4, 4) = real(-2.0_ki*mT**2, ki_sam)
   !-----------#] initialize invariants:

   if(samurai_group_numerators) then
      !------#[ reduce numerator numeval_group5:
      if(samurai_verbosity > 0) then
         write(samurai_out,*) "[golem-2.0] numeval_group5"
         write(samurai_out,*) "[golem-2.0] epspow=", epspow
      end if
      !-----------#[ initialize invariants:
      allocate(s_mat(4, 4))
      s_mat(:,:) = g_mat(:,:)
      !-----------#] initialize invariants:
      call samurai_rm(numeval_group5, tot, totr, Vi, msq, 4, &
         & effective_group_rank, istopm, scale2, ok, &
         & samurai_cache_flag_g5, samurai_cache_g5)
      !-----------#[ deallocate invariants:
      deallocate(s_mat)
      !-----------#] deallocate invariants:

      !------#] reduce numerator numeval_group5:
   else
      !------#[ sum over reduction of single diagrams:
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='45'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram45"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(4, 4))
         s_mat(:,:) = g_mat( (/1,2,3,4/), (/1,2,3,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram45, acc, accr, &
            & Vi((/1,2,3,4/),:), msq((/1,2,3,4/)), 4, &
            & 4, istopm, scale2, ok, &
            & samurai_cache_flag_d45, samurai_cache_d45)
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
         if(debug_nlo_diagrams) then
            write(logfile,*) "<diagram index='277'>"
         end if
         if(samurai_verbosity > 0) then
            write(samurai_out,*) "[golem-2.0] numerator_diagram277"
            write(samurai_out,*) "[golem-2.0] epspow=", epspow
         end if
         !-----------#[ initialize invariants:
         allocate(s_mat(3, 3))
         s_mat(:,:) = g_mat( (/1,2,4/), (/1,2,4/) )
         !-----------#] initialize invariants:
         call samurai_rm(numerator_diagram277, acc, accr, &
            & Vi((/1,2,4/),:), msq((/1,2,4/)), 3, &
            & 3, istopm, scale2, acc_ok, &
            & samurai_cache_flag_d277, samurai_cache_d277)
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

         tot = tot  + acc
         totr = totr  + accr
         ok = ok .and. acc_ok
      !------#] sum over reduction of single diagrams:
   end if
end subroutine reduce_group5
!-----#] subroutine reduce_group5:
!---#] reduce groups with samurai:
end module p0_dbard_hepemg_samuraih3
