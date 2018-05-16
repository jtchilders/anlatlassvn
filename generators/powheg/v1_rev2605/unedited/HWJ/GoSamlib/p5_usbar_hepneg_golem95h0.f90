module     p5_usbar_hepneg_golem95h0
   use precision_golem, only: ki_gol => ki
   use p5_usbar_hepneg_config, only: ki
   implicit none
   private
   interface reconstruct_group
      module procedure reconstruct_group0
      module procedure reconstruct_group1
      module procedure reconstruct_group2
      module procedure reconstruct_group3
   end interface

   public :: reconstruct_group
contains
!---#[ subroutine reconstruct_group0:
subroutine     reconstruct_group0(coeffs)
   use tens_rec
   use p5_usbar_hepneg_config
   use p5_usbar_hepneg_groups, only: tensrec_info_group0
   use p5_usbar_hepneg_d70h0l1, only: numerator_d70 => numerator_golem95
   use p5_usbar_hepneg_d70h0l1d, only: reconstruct_d70
   implicit none
   type(tensrec_info_group0), intent(out) :: coeffs
   !------#[ Diagram 70:
      if (tens_rec_by_derivatives) then
         call reconstruct_d70(coeffs)
      else
         call reconstruct3(numerator_d70, coeffs%coeffs_70)
      end if
   !------#] Diagram 70:
end subroutine reconstruct_group0
!---#] subroutine reconstruct_group0:
!---#[ subroutine reconstruct_group1:
subroutine     reconstruct_group1(coeffs)
   use tens_rec
   use p5_usbar_hepneg_config
   use p5_usbar_hepneg_groups, only: tensrec_info_group1
   use p5_usbar_hepneg_d61h0l1, only: numerator_d61 => numerator_golem95
   use p5_usbar_hepneg_d61h0l1d, only: reconstruct_d61
   implicit none
   type(tensrec_info_group1), intent(out) :: coeffs
   !------#[ Diagram 61:
      if (tens_rec_by_derivatives) then
         call reconstruct_d61(coeffs)
      else
         call reconstruct3(numerator_d61, coeffs%coeffs_61)
      end if
   !------#] Diagram 61:
end subroutine reconstruct_group1
!---#] subroutine reconstruct_group1:
!---#[ subroutine reconstruct_group2:
subroutine     reconstruct_group2(coeffs)
   use tens_rec
   use p5_usbar_hepneg_config
   use p5_usbar_hepneg_groups, only: tensrec_info_group2
   use p5_usbar_hepneg_d34h0l1, only: numerator_d34 => numerator_golem95
   use p5_usbar_hepneg_d34h0l1d, only: reconstruct_d34
   use p5_usbar_hepneg_d40h0l1, only: numerator_d40 => numerator_golem95
   use p5_usbar_hepneg_d40h0l1d, only: reconstruct_d40
   use p5_usbar_hepneg_d60h0l1, only: numerator_d60 => numerator_golem95
   use p5_usbar_hepneg_d60h0l1d, only: reconstruct_d60
   use p5_usbar_hepneg_d92h0l1, only: numerator_d92 => numerator_golem95
   use p5_usbar_hepneg_d92h0l1d, only: reconstruct_d92
   implicit none
   type(tensrec_info_group2), intent(out) :: coeffs
   !------#[ Diagram 34:
      if (tens_rec_by_derivatives) then
         call reconstruct_d34(coeffs)
      else
         call reconstruct2(numerator_d34, coeffs%coeffs_34)
      end if
   !------#] Diagram 34:
   !------#[ Diagram 40:
      if (tens_rec_by_derivatives) then
         call reconstruct_d40(coeffs)
      else
         call reconstruct1(numerator_d40, coeffs%coeffs_40)
      end if
   !------#] Diagram 40:
   !------#[ Diagram 60:
      if (tens_rec_by_derivatives) then
         call reconstruct_d60(coeffs)
      else
         call reconstruct3(numerator_d60, coeffs%coeffs_60)
      end if
   !------#] Diagram 60:
   !------#[ Diagram 92:
      if (tens_rec_by_derivatives) then
         call reconstruct_d92(coeffs)
      else
         call reconstruct2(numerator_d92, coeffs%coeffs_92)
      end if
   !------#] Diagram 92:
end subroutine reconstruct_group2
!---#] subroutine reconstruct_group2:
!---#[ subroutine reconstruct_group3:
subroutine     reconstruct_group3(coeffs)
   use tens_rec
   use p5_usbar_hepneg_config
   use p5_usbar_hepneg_groups, only: tensrec_info_group3
   use p5_usbar_hepneg_d30h0l1, only: numerator_d30 => numerator_golem95
   use p5_usbar_hepneg_d30h0l1d, only: reconstruct_d30
   use p5_usbar_hepneg_d37h0l1, only: numerator_d37 => numerator_golem95
   use p5_usbar_hepneg_d37h0l1d, only: reconstruct_d37
   use p5_usbar_hepneg_d59h0l1, only: numerator_d59 => numerator_golem95
   use p5_usbar_hepneg_d59h0l1d, only: reconstruct_d59
   use p5_usbar_hepneg_d94h0l1, only: numerator_d94 => numerator_golem95
   use p5_usbar_hepneg_d94h0l1d, only: reconstruct_d94
   implicit none
   type(tensrec_info_group3), intent(out) :: coeffs
   !------#[ Diagram 30:
      if (tens_rec_by_derivatives) then
         call reconstruct_d30(coeffs)
      else
         call reconstruct2(numerator_d30, coeffs%coeffs_30)
      end if
   !------#] Diagram 30:
   !------#[ Diagram 37:
      if (tens_rec_by_derivatives) then
         call reconstruct_d37(coeffs)
      else
         call reconstruct1(numerator_d37, coeffs%coeffs_37)
      end if
   !------#] Diagram 37:
   !------#[ Diagram 59:
      if (tens_rec_by_derivatives) then
         call reconstruct_d59(coeffs)
      else
         call reconstruct3(numerator_d59, coeffs%coeffs_59)
      end if
   !------#] Diagram 59:
   !------#[ Diagram 94:
      if (tens_rec_by_derivatives) then
         call reconstruct_d94(coeffs)
      else
         call reconstruct2(numerator_d94, coeffs%coeffs_94)
      end if
   !------#] Diagram 94:
end subroutine reconstruct_group3
!---#] subroutine reconstruct_group3:
end module p5_usbar_hepneg_golem95h0
