!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECinitialization.f90.

module dpinitialization
  use types 
  use define_ampl; use dpglobal; use dpmemory; use dpspinors 
  use sub_defs_io 
  implicit none 
  private 

  public :: dpinitialize


contains
  
  subroutine dpinitialize(np) 
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
    MCphi = 0._dp 

  end subroutine dpinitialize
  
end module dpinitialization
