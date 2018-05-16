!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECauxiliary_functions.f90.

module mpauxiliary_functions 
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types; use consts_mp 
  use warnings_and_errors 
  implicit none 
  private 
  
  interface dot 
     module procedure dotrr,dotrc,dotcr,dotcc
  end interface

  interface quad 
     module procedure quadc,quadr
  end interface

  interface pmass 
     module procedure pmassc,pmassr
  end interface

  public :: dot,quad,pmass

contains

  function dotrr(p1,p2) result(dot) 
    type(mp_real), intent(in) :: p1(:),p2(:)
    type(mp_real) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotrr

  function dotrc(p1,p2) result(dot) 
    type(mp_real), intent(in) :: p1(:)
    type(mp_complex), intent(in) :: p2(:)
    type(mp_complex) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:)) 

  end function dotrc

  function dotcr(p1,p2) result(dot) 
    type(mp_real), intent(in) :: p2(:)
    type(mp_complex), intent(in) :: p1(:)
    type(mp_complex) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcr

  function dotcc(p1,p2) result(dot) 
    type(mp_complex), intent(in) :: p1(:),p2(:)
    type(mp_complex) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcc

  function quadc(p) result(quad) 
    type(mp_complex), intent(in) :: p(:)
    type(mp_complex) :: quad

    quad = p(1)**2-sum(p(2:)**2)

   end function quadc

  function quadr(p) result(quad) 
    type(mp_real), intent(in) :: p(:)
    type(mp_real) :: quad

    quad = p(1)**2-sum(p(2:)**2)

  end function quadr

  function pmassc(p) result(res) 
    type(mp_complex), intent(in) :: p(:)
    type(mp_complex) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassc

  function pmassr(p) result(res) 
    type(mp_real), intent(in) :: p(:)
    type(mp_real) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassr

 
end module mpauxiliary_functions



