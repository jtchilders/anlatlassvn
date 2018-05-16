!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECinitialization.f90.

module qpinitialization
  use types 
  use define_ampl; use qpglobal; use qpmemory; use qpspinors 
  use sub_defs_io 
  implicit none 
  private 

  public :: qpinitialize


contains
  
  subroutine qpinitialize(np) 
    integer, intent(in) :: np 

    write(*,*) '# INITILIZATION OF PARAMETERS & MEMORY'  
    write(*,*) '# ',trim(command_line())
    ! -- initialize here all commond-line parameters 
    call initialize_params(np)
    ! -- initialization 
    call init_spinors 
    call fix_bvec 
    
    ! -- memory 
    call allocate_mem
    MCphi = 0._qp 

  end subroutine qpinitialize
  
end module qpinitialization
