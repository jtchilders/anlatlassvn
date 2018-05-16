!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECinitialization.f90.

module mpinitialization
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types 
  use define_ampl; use mpglobal; use mpmemory; use mpspinors 
  use sub_defs_io 
  implicit none 
  private 

  public :: mpinitialize


contains
  
  subroutine mpinitialize(np) 
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
    MCphi = mpreal('0.') 

  end subroutine mpinitialize
  
end module mpinitialization
