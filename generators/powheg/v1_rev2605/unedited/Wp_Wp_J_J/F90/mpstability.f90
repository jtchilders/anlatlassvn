!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECstability.f90.

!!! module which contains routine and global variables needed to find
!!! instabilities.
!
module mpstability 
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types; use consts_mp 
  implicit none 
  private 

  ! logical variable set to false at the beginning of any event, and set to 
  ! true if somewhere we divided by a small number 
  ! at the very end we check if division_by_small is true we redo the event 
  ! in mp
  
  public :: check_if_small,check_res_fails 


  type(mp_real), public :: abs_smallest_den 
  logical, public    :: division_by_small, poles_dont_match  
  logical, public    :: division_by_zero
  logical, public    :: rescheck_fails(4) = .false. 
  logical, public    :: rescheck_fails_loose(4) = .false. 
  logical, public    :: rescheck_fails_veryloose(4) = .false. 
  integer, public    :: n_rescheck_fails(4)
  integer, public    :: n_rescheck_good(4)

  logical, public    :: tricheck_sltr_fails,tricheck_lltr_fails  
  logical, public    :: bubcheck1_sltr_fails, bubcheck2_lltr_fails  
  logical, public    :: deemed_unstable 
  logical, public    :: deemed_unstable1 
  logical, public    :: deemed_unstable2 
  logical, public    :: deemed_unstable3 
  real(dp), public   :: abs_rel_err 

  logical, private :: verbose_stability = .false. 

contains

    subroutine check_if_small(den,call_from)
      type(mp_complex), intent(in) :: den 
      character(len=*), intent(in), optional :: call_from 
      type(mp_real) :: aden 

      !if (verbose_stability) write(*,*) 'in check_if_small' 
      aden = abs(den) 
      if (aden < sq2tol) then 
         division_by_small = .true. 
      endif
      if (aden < abs_smallest_den) then 
         abs_smallest_den = aden
         if (verbose_stability) write(*,*) 'smallest den is:', abs_smallest_den
         if (present(call_from) .and. verbose_stability) &
              &write(*,*) 'call from ', call_from 
      endif
      if (aden == zero) division_by_zero = .true. 
      

    end subroutine check_if_small

    subroutine check_res_fails(res_v1,res_v2,i)
      type(mp_complex), intent(in) :: res_v1,res_v2 
      integer :: i 
      logical :: fails 

      fails = .false. 
      if (abs(res_v2) > sq2tol) then 
          if (abs((res_v1-res_v2)/res_v2) > 100*sq2tol) fails = .true. 
       else
          if (abs((res_v1-res_v2)) > 100*sq2tol) fails = .true. 
       endif

      if (fails) then 
         if (i ==1) then 
            !write(*,*) 'tri:large ltr ZERO?'
            rescheck_fails(1) = .true. 
         elseif (i ==2) then 
            !write(*,*) 'tri:small ltr ZERO?'
            rescheck_fails(2) = .true. 
         elseif (i ==3) then 
            !write(*,*) 'bub:large ltr ZERO?'
            rescheck_fails(3) = .true. 
         elseif (i ==4) then 
            !write(*,*) 'bub:small ltr ZERO?'
            rescheck_fails(4) = .true. 
         else
            write(*,*) 'Check number',i 
            stop 'check_res_fails: i out of range'
         end if
         !write(*,*) mptodp(abs((res_v1-res_v2)/res_v1)),mptodp(abs(res_v1-res_v2)),mptodp(abs(res_v1)),mptodp(abs(res_v2))
      end if

      fails = .false. 
      if (abs(res_v2) > sq2tol) then 
          if (abs((res_v1-res_v2)/res_v2) > 10000*sq2tol) fails = .true. 
       else
          if (abs((res_v1-res_v2)) > 10000*sq2tol) fails = .true. 
       endif

      if (fails) then 
         if (i ==1) then 
            !write(*,*) 'tri:large ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_loose(1) = .true. 
         elseif (i ==2) then 
            !write(*,*) 'tri:small ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_loose(2) = .true. 
         elseif (i ==3) then 
            !write(*,*) 'bub:large ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_loose(3) = .true. 
         elseif (i ==4) then 
            !write(*,*) 'bub:small ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_loose(4) = .true. 
         else
            write(*,*) 'Check number',i 
            stop 'check_res_fails: i out of range'
         end if
      end if

      fails = .false. 
      if (abs(res_v2) > sq2tol) then 
          if (abs((res_v1-res_v2)/res_v2) > mpreal('0.01')) fails = .true. 
       else
          if (abs((res_v1-res_v2)) > mpreal('0.01')) fails = .true. 
       endif

      if (fails) then 
         if (i ==1) then 
            !write(*,*) 'tri:large ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_veryloose(1) = .true. 
         elseif (i ==2) then 
            !write(*,*) 'tri:small ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_veryloose(2) = .true. 
         elseif (i ==3) then 
            !write(*,*) 'bub:large ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_veryloose(3) = .true. 
         elseif (i ==4) then 
            !write(*,*) 'bub:small ltr ZERO?',mptodp(res_v1-res_v2)
            rescheck_fails_veryloose(4) = .true. 
         else
            write(*,*) 'Check number',i 
            stop 'check_res_fails: i out of range'
         end if
      end if

    end subroutine check_res_fails



end module mpstability

