!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECauxiliary_functions.f90.

module dpauxiliary_functions 
  use types; use consts_dp 
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
    real(dp), intent(in) :: p1(:),p2(:)
    real(dp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotrr

  function dotrc(p1,p2) result(dot) 
    real(dp), intent(in) :: p1(:)
    complex(dp), intent(in) :: p2(:)
    complex(dp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:)) 

  end function dotrc

  function dotcr(p1,p2) result(dot) 
    real(dp), intent(in) :: p2(:)
    complex(dp), intent(in) :: p1(:)
    complex(dp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcr

  function dotcc(p1,p2) result(dot) 
    complex(dp), intent(in) :: p1(:),p2(:)
    complex(dp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcc

  function quadc(p) result(quad) 
    complex(dp), intent(in) :: p(:)
    complex(dp) :: quad

    quad = p(1)**2-sum(p(2:)**2)

   end function quadc

  function quadr(p) result(quad) 
    real(dp), intent(in) :: p(:)
    real(dp) :: quad

    quad = p(1)**2-sum(p(2:)**2)

  end function quadr

  function pmassc(p) result(res) 
    complex(dp), intent(in) :: p(:)
    complex(dp) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassc

  function pmassr(p) result(res) 
    real(dp), intent(in) :: p(:)
    real(dp) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassr

 
end module dpauxiliary_functions



