module p17_bbarg_hepembbar_matrix
   use p16_bbarb_hepemg_config, only: ki
   implicit none
   private

   public :: initgolem, exitgolem, samplitude
   public :: samplitudel0, samplitudel1, ir_subtraction
   ! TODO:
   ! public :: color_correlated_lo2, spin_correlated_lo2


contains

   pure function prefactor()
      use p16_bbarb_hepemg_config, only: &
      & include_color_avg_factor, include_helicity_avg_factor, &
      & include_symmetry_factor
      use p16_bbarb_hepemg_model, only: NC
      use p16_bbarb_hepemg_kinematics, only: &
      & in_helicities, symmetry_factor
      use p16_bbarb_hepemg_color, only: incolors
      implicit none
      real(ki) :: prefactor, NA

      prefactor = 1.0_ki
      NA = (NC + 1.0_ki) * (NC - 1.0_ki)

      if (include_color_avg_factor) then
         prefactor = prefactor * real(incolors, ki)/NC/NA
      end if
      if (include_helicity_avg_factor) then
         prefactor = prefactor * real(in_helicities, ki)/2.0_ki/2.0_ki
      end if
      if (include_symmetry_factor) then
         prefactor = prefactor * real(symmetry_factor, ki)/1.0_ki
      end if
   end function prefactor

   subroutine     initgolem(is_first)
      use p16_bbarb_hepemg_matrix, only: orig_initgolem => initgolem
      implicit none
      logical, optional, intent(in) :: is_first

      if (present(is_first)) then
         call orig_initgolem(is_first)
      else
         call orig_initgolem()
      endif
   end subroutine initgolem

   subroutine     exitgolem(is_last)
      use p16_bbarb_hepemg_matrix, only: orig_exitgolem => exitgolem
      implicit none
      logical, optional, intent(in) :: is_last

      if (present(is_last)) then
         call orig_exitgolem(is_last)
      else
         call orig_exitgolem()
      endif
   end subroutine exitgolem

   subroutine     samplitude(vecs, scale2, amp, ok, h)
      use p16_bbarb_hepemg_matrix, only: orig_sub => samplitude
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      real(ki), dimension(4), intent(out) :: amp
      logical, intent(out), optional :: ok
      integer, intent(in), optional :: h

      real(ki), dimension(6, 4) :: new_vecs

      call twist_momenta(vecs, new_vecs)

      if (present(ok)) then
         if (present(h)) then
            call orig_sub(new_vecs, scale2, amp, ok, h)
         else
            call orig_sub(new_vecs, scale2, amp, ok)
         end if
      else
         call orig_sub(new_vecs, scale2, amp, ok)
      end if

      amp = amp * prefactor()
   end subroutine samplitude

   function     samplitudel0(vecs, h) result(amp)
      use p16_bbarb_hepemg_matrix, only: orig_func => samplitudel0
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      integer, optional, intent(in) :: h
      real(ki) :: amp
      real(ki), dimension(6, 4) :: new_vecs

      call twist_momenta(vecs, new_vecs)

      if (present(h)) then
         amp = orig_func(new_vecs, h)
      else
         amp = orig_func(new_vecs)
      end if
      amp = amp * prefactor()
   end function samplitudel0

   function     samplitudel1(vecs,scale2,ok,rat2,h) result(amp)
      use p16_bbarb_hepemg_matrix, only: orig_func => samplitudel1
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      logical, intent(out) :: ok
      real(ki), intent(in) :: scale2
      real(ki), intent(out) :: rat2
      integer, optional, intent(in) :: h
      real(ki), dimension(-2:0) :: amp
      real(ki), dimension(6, 4) :: new_vecs

      call twist_momenta(vecs, new_vecs)

      if (present(h)) then
         amp = orig_func(new_vecs, scale2, ok, rat2, h)
      else
         amp = orig_func(new_vecs, scale2, ok, rat2)
      end if
      amp = amp * prefactor()
   end function samplitudel1

   subroutine     ir_subtraction(vecs,scale2,amp)
      use p16_bbarb_hepemg_matrix, &
      & only: orig_sub => ir_subtraction
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), intent(in) :: scale2
      real(ki), dimension(2), intent(out) :: amp

      real(ki), dimension(6, 4) :: new_vecs

      call twist_momenta(vecs, new_vecs)

      call orig_sub(new_vecs, scale2, amp)
      amp = amp * prefactor()
   end subroutine ir_subtraction

   pure subroutine twist_momenta(vecs, new_vecs)
      implicit none
      real(ki), dimension(6, 4), intent(in) :: vecs
      real(ki), dimension(6, 4), intent(out) :: new_vecs
      new_vecs(1,:) = + vecs(1,:)
      new_vecs(6,:) = - vecs(2,:)
      new_vecs(3,:) = + vecs(3,:)
      new_vecs(4,:) = + vecs(4,:)
      new_vecs(5,:) = + vecs(5,:)
      new_vecs(2,:) = - vecs(6,:)
   end  subroutine twist_momenta

end module p17_bbarg_hepembbar_matrix
