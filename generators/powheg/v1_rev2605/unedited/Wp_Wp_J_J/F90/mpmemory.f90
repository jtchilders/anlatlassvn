!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECmemory.f90.

module mpmemory 
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types; use consts_mp 
  use define_ampl 
  use sub_defs_io 
  implicit none
  private 
  
  integer, allocatable, public       :: vgluon_int(:) 

  type(mp_complex), allocatable, public :: ExtCurrent(:,:) 
  logical, allocatable, public :: log_ExtCurrent(:) 

  type(mp_complex), allocatable, public :: IntCurrent(:,:,:) 
  logical, allocatable, public :: log_IntCurrent(:,:) 

  type(mp_complex), allocatable, public :: IntQCurrent(:,:,:) 
  logical, allocatable, public :: log_IntQCurrent(:,:) 



  public :: allocate_mem, riinitialize_mem, riinitialize_Extmem, memory_check, store_result 

  logical :: firsttime = .true.   

contains 

  subroutine allocate_mem

    if (.not. allocated(IntCurrent)) then 
       allocate(IntCurrent(16,8,0:2**npoint-1))
       allocate(log_IntCurrent(8,0:2**npoint-1))
       IntCurrent = czero 
       log_IntCurrent = .false. 
       
       
       allocate(IntQCurrent(16,8,0:2**npoint-1))
       allocate(log_IntQCurrent(8,0:2**npoint-1))
       IntQCurrent = czero 
       log_IntQCurrent = .false. 
       
       
       allocate(ExtCurrent(4,2**npoint-1))
       allocate(log_ExtCurrent(2**npoint-1))
       ExtCurrent = czero 
       log_ExtCurrent = .false. 
    endif

  end subroutine allocate_mem

  subroutine riinitialize_mem
    
    log_IntCurrent = .false. 
    log_IntQCurrent = .false. 
    
  end subroutine riinitialize_mem

  subroutine riinitialize_Extmem
    
    log_ExtCurrent = .false. 

  end subroutine riinitialize_Extmem


  subroutine memory_check(pol_int,res,done,giarray,qiarray,wid1,wid2)
    integer, intent(in)           :: giarray(:),pol_int
    integer, intent(in), optional :: qiarray(:),wid1,wid2
    type(mp_complex), intent(out)    :: res(:) 
    logical, intent(out)          :: done 
    integer  :: tot
    logical :: ext, int_quark 

    if (firsttime) then 
       cashing = log_val_opt('-cashing',.false.) 
       if (cashing) then 
          write(*,*) 'Cashing of previouly computed scaler integrals'
       else
          write(*,*) 'No cashing of previouly computed scaler integrals'
       endif
       firsttime = .false. 
    endif
    if (.not. cashing) then 
       done = .false.
       return
    endif


    ! -- GZ need to change cashing for ferm-loops Z because of
    ! possible different ordering of gluons wrt fermion pair
    if ((ferm_loops_Z) .or. ferm_loops_Z_sbs) then 
       done = .false.
       return 
    endif

    tot = sum(giarray)
    if (present(qiarray)) tot = tot+sum(qiarray)

    if (tot < 0) then 
       done = .false. 
       return ! means is was a call without memory implemented 
    endif
    if (present(wid1)) tot = tot+wid1
    if (present(wid2)) tot = tot+wid2


    ext = .true. 
    if (size(giarray) >0) then 
       if (giarray(size(giarray)) == 0) then 
          ext = .false. 
          int_quark = .false. 
       endif
    end if
    if (present(qiarray)) then 
       if (qiarray(size(qiarray)) == 0) then 
          ext = .false. 
          int_quark = .true. 
       endif
    endif

    if (tot > size(log_ExtCurrent)) then 
       write(*,*) 'giarray ',giarray 
       if (present(qiarray)) write(*,*) 'qiarray',qiarray 
       if (present(wid1)) write(*,*) 'wid1',wid1
       if (present(wid2)) write(*,*) 'wid2',wid2
       write(*,*) 'tot', tot, size(log_ExtCurrent)
       stop 'tot too large' 
    endif

    if (ext) then 
       ! --  checking ext current 
       if (log_ExtCurrent(tot)) then 
          res(:4) = ExtCurrent(:4,tot)
          res(5:) = czero 
          done = .true. 
       else
          done = .false.
          res(:) = czero 
       end if
       
    else

       if (int_quark) then 
          ! -- checking internal current 
          if (log_IntQCurrent(pol_int,tot)) then 
             res = IntQCurrent(:size(res),pol_int,tot)
             done = .true. 
          else
             done = .false.
             res(:) = czero 
          endif
       else
          ! -- checking internal current 
          if (log_IntCurrent(pol_int,tot)) then 
             res = IntCurrent(:size(res),pol_int,tot)
             done = .true. 
          else
             done = .false.
             res(:) = czero 
          endif
       end if
    endif
  end subroutine memory_check


  subroutine store_result(pol_int,res,giarray,qiarray,wid1,wid2)
    integer, intent(in)           :: giarray(:),pol_int
    integer, intent(in), optional :: qiarray(:),wid1,wid2
    type(mp_complex), intent(in)     :: res(:) 
    integer  :: tot
    logical :: ext, int_quark  

    tot = sum(giarray)
    if (present(qiarray)) tot = tot+sum(qiarray)
    if (present(wid1)) tot = tot+wid1
    if (present(wid2)) tot = tot+wid2


    ext = .true. 
    if (size(giarray) >0) then 
       if (giarray(size(giarray)) == 0) then 
          ext = .false. 
          int_quark = .false. 
       endif
    end if

    if (present(qiarray)) then 
       if (qiarray(size(qiarray)) == 0) then 
          ext = .false. 
          int_quark = .true. 
       endif
    endif

    if (ext) then 
       ExtCurrent(:4,tot) = res(:4) 
       log_ExtCurrent(tot) = .true. 
    else
       if (int_quark) then 
          IntQCurrent(:size(res),pol_int,tot) = res(:) 
          log_IntQCurrent(pol_int,tot) = .true. 
       else
          IntCurrent(:size(res),pol_int,tot) = res(:) 
          log_IntCurrent(pol_int,tot) = .true. 
       endif
    endif

  end subroutine store_result




end module mpmemory
