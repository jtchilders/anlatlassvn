module     p4_ubaru_hepemg_golem95h2
   use precision_golem, only: ki_gol => ki
   use p4_ubaru_hepemg_config, only: ki
   implicit none
   private
   interface reconstruct_group
      module procedure reconstruct_group0
      module procedure reconstruct_group1
      module procedure reconstruct_group2
      module procedure reconstruct_group3
      module procedure reconstruct_group4
      module procedure reconstruct_group5
   end interface

   public :: reconstruct_group
contains
!---#[ subroutine reconstruct_group0:
subroutine     reconstruct_group0(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group0
   use p4_ubaru_hepemg_d541h2l1, only: numerator_d541 => numerator_golem95
   use p4_ubaru_hepemg_d541h2l1d, only: reconstruct_d541
   implicit none
   type(tensrec_info_group0), intent(out) :: coeffs
   !------#[ Diagram 541:
      if (tens_rec_by_derivatives) then
         call reconstruct_d541(coeffs)
      else
         call reconstruct3(numerator_d541, coeffs%coeffs_541)
      end if
   !------#] Diagram 541:
end subroutine reconstruct_group0
!---#] subroutine reconstruct_group0:
!---#[ subroutine reconstruct_group1:
subroutine     reconstruct_group1(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group1
   use p4_ubaru_hepemg_d279h2l1, only: numerator_d279 => numerator_golem95
   use p4_ubaru_hepemg_d279h2l1d, only: reconstruct_d279
   use p4_ubaru_hepemg_d333h2l1, only: numerator_d333 => numerator_golem95
   use p4_ubaru_hepemg_d333h2l1d, only: reconstruct_d333
   use p4_ubaru_hepemg_d421h2l1, only: numerator_d421 => numerator_golem95
   use p4_ubaru_hepemg_d421h2l1d, only: reconstruct_d421
   use p4_ubaru_hepemg_d537h2l1, only: numerator_d537 => numerator_golem95
   use p4_ubaru_hepemg_d537h2l1d, only: reconstruct_d537
   use p4_ubaru_hepemg_d605h2l1, only: numerator_d605 => numerator_golem95
   use p4_ubaru_hepemg_d605h2l1d, only: reconstruct_d605
   implicit none
   type(tensrec_info_group1), intent(out) :: coeffs
   !------#[ Diagram 279:
      if (tens_rec_by_derivatives) then
         call reconstruct_d279(coeffs)
      else
         call reconstruct3(numerator_d279, coeffs%coeffs_279)
      end if
   !------#] Diagram 279:
   !------#[ Diagram 333:
      if (tens_rec_by_derivatives) then
         call reconstruct_d333(coeffs)
      else
         call reconstruct2(numerator_d333, coeffs%coeffs_333)
      end if
   !------#] Diagram 333:
   !------#[ Diagram 421:
      if (tens_rec_by_derivatives) then
         call reconstruct_d421(coeffs)
      else
         call reconstruct1(numerator_d421, coeffs%coeffs_421)
      end if
   !------#] Diagram 421:
   !------#[ Diagram 537:
      if (tens_rec_by_derivatives) then
         call reconstruct_d537(coeffs)
      else
         call reconstruct3(numerator_d537, coeffs%coeffs_537)
      end if
   !------#] Diagram 537:
   !------#[ Diagram 605:
      if (tens_rec_by_derivatives) then
         call reconstruct_d605(coeffs)
      else
         call reconstruct2(numerator_d605, coeffs%coeffs_605)
      end if
   !------#] Diagram 605:
end subroutine reconstruct_group1
!---#] subroutine reconstruct_group1:
!---#[ subroutine reconstruct_group2:
subroutine     reconstruct_group2(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group2
   use p4_ubaru_hepemg_d323h2l1, only: numerator_d323 => numerator_golem95
   use p4_ubaru_hepemg_d323h2l1d, only: reconstruct_d323
   use p4_ubaru_hepemg_d413h2l1, only: numerator_d413 => numerator_golem95
   use p4_ubaru_hepemg_d413h2l1d, only: reconstruct_d413
   use p4_ubaru_hepemg_d533h2l1, only: numerator_d533 => numerator_golem95
   use p4_ubaru_hepemg_d533h2l1d, only: reconstruct_d533
   use p4_ubaru_hepemg_d613h2l1, only: numerator_d613 => numerator_golem95
   use p4_ubaru_hepemg_d613h2l1d, only: reconstruct_d613
   implicit none
   type(tensrec_info_group2), intent(out) :: coeffs
   !------#[ Diagram 323:
      if (tens_rec_by_derivatives) then
         call reconstruct_d323(coeffs)
      else
         call reconstruct2(numerator_d323, coeffs%coeffs_323)
      end if
   !------#] Diagram 323:
   !------#[ Diagram 413:
      if (tens_rec_by_derivatives) then
         call reconstruct_d413(coeffs)
      else
         call reconstruct1(numerator_d413, coeffs%coeffs_413)
      end if
   !------#] Diagram 413:
   !------#[ Diagram 533:
      if (tens_rec_by_derivatives) then
         call reconstruct_d533(coeffs)
      else
         call reconstruct3(numerator_d533, coeffs%coeffs_533)
      end if
   !------#] Diagram 533:
   !------#[ Diagram 613:
      if (tens_rec_by_derivatives) then
         call reconstruct_d613(coeffs)
      else
         call reconstruct2(numerator_d613, coeffs%coeffs_613)
      end if
   !------#] Diagram 613:
end subroutine reconstruct_group2
!---#] subroutine reconstruct_group2:
!---#[ subroutine reconstruct_group3:
subroutine     reconstruct_group3(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group3
   use p4_ubaru_hepemg_d101h2l1, only: numerator_d101 => numerator_golem95
   use p4_ubaru_hepemg_d101h2l1d, only: reconstruct_d101
   use p4_ubaru_hepemg_d553h2l1, only: numerator_d553 => numerator_golem95
   use p4_ubaru_hepemg_d553h2l1d, only: reconstruct_d553
   implicit none
   type(tensrec_info_group3), intent(out) :: coeffs
   !------#[ Diagram 101:
      if (tens_rec_by_derivatives) then
         call reconstruct_d101(coeffs)
      else
         call reconstruct4(numerator_d101, coeffs%coeffs_101)
      end if
   !------#] Diagram 101:
   !------#[ Diagram 553:
      if (tens_rec_by_derivatives) then
         call reconstruct_d553(coeffs)
      else
         call reconstruct3(numerator_d553, coeffs%coeffs_553)
      end if
   !------#] Diagram 553:
end subroutine reconstruct_group3
!---#] subroutine reconstruct_group3:
!---#[ subroutine reconstruct_group4:
subroutine     reconstruct_group4(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group4
   use p4_ubaru_hepemg_d77h2l1, only: numerator_d77 => numerator_golem95
   use p4_ubaru_hepemg_d77h2l1d, only: reconstruct_d77
   implicit none
   type(tensrec_info_group4), intent(out) :: coeffs
   !------#[ Diagram 77:
      if (tens_rec_by_derivatives) then
         call reconstruct_d77(coeffs)
      else
         call reconstruct4(numerator_d77, coeffs%coeffs_77)
      end if
   !------#] Diagram 77:
end subroutine reconstruct_group4
!---#] subroutine reconstruct_group4:
!---#[ subroutine reconstruct_group5:
subroutine     reconstruct_group5(coeffs)
   use tens_rec
   use p4_ubaru_hepemg_config
   use p4_ubaru_hepemg_groups, only: tensrec_info_group5
   use p4_ubaru_hepemg_d45h2l1, only: numerator_d45 => numerator_golem95
   use p4_ubaru_hepemg_d45h2l1d, only: reconstruct_d45
   use p4_ubaru_hepemg_d277h2l1, only: numerator_d277 => numerator_golem95
   use p4_ubaru_hepemg_d277h2l1d, only: reconstruct_d277
   implicit none
   type(tensrec_info_group5), intent(out) :: coeffs
   !------#[ Diagram 45:
      if (tens_rec_by_derivatives) then
         call reconstruct_d45(coeffs)
      else
         call reconstruct4(numerator_d45, coeffs%coeffs_45)
      end if
   !------#] Diagram 45:
   !------#[ Diagram 277:
      if (tens_rec_by_derivatives) then
         call reconstruct_d277(coeffs)
      else
         call reconstruct3(numerator_d277, coeffs%coeffs_277)
      end if
   !------#] Diagram 277:
end subroutine reconstruct_group5
!---#] subroutine reconstruct_group5:
end module p4_ubaru_hepemg_golem95h2
