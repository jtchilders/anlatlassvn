!----------------------------------------------------------------------
! Defines kind parameters
module types
  implicit none
  integer, parameter  :: qp = selected_real_kind(30)
  integer, parameter  :: dp = selected_real_kind(15)
  integer, parameter  :: sp = selected_real_kind(6)

  integer, parameter  :: dp15 = selected_real_kind(15)

end module types



