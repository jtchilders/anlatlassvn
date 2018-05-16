! ===================================================================
!! Module containing the routines to change the multiple precision 
!! (SetPrecision) and to get the current precision (GetPrecision). 
!! This module is adapted specifically to Baileys MP package.  
module mpprecision 
  use types; use consts_dp 
  use mpmodule 
  implicit none 
  private 

  public :: GetPrecision, SetPrecision, MaxPrecision


contains 

  subroutine GetPrecision(precision) 
    integer, intent(out) :: precision 

    precision = int((mpnw-1)*24*log(two)/Log(10._dp))

  end subroutine GetPrecision

  subroutine SetPrecision(precision) 
    integer, intent(in) :: precision 

    ! The correct relation uses int (see. pag. 7 of mpf90.ps ) 
    ! but maybe safer to use ceiling
    ! mpnw = int(precision*Log(10._dp)/(24._dp*log(two)))+1
    mpnw = ceiling(precision*Log(10._dp)/(24._dp*log(two)))+1

    if (mpnw >= mpwds+1) stop 'SetPrecision: Desired precision &
         &in machine words exceeds the maximal allowed precision'

  end subroutine SetPrecision

  function MaxPrecision() 
    integer :: MaxPrecision 
    integer :: mpnw_lcl

    mpnw_lcl = mpwds-1 
    MaxPrecision = int((mpnw_lcl-1)*24*log(two)/Log(10._dp))    

  end function MaxPrecision


end module mpprecision
