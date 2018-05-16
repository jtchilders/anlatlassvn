!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECauxiliary_functions.f90.

module qpauxiliary_functions 
  use types; use consts_qp 
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
    real(qp), intent(in) :: p1(:),p2(:)
    real(qp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotrr

  function dotrc(p1,p2) result(dot) 
    real(qp), intent(in) :: p1(:)
    complex(qp), intent(in) :: p2(:)
    complex(qp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:)) 

  end function dotrc

  function dotcr(p1,p2) result(dot) 
    real(qp), intent(in) :: p2(:)
    complex(qp), intent(in) :: p1(:)
    complex(qp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcr

  function dotcc(p1,p2) result(dot) 
    complex(qp), intent(in) :: p1(:),p2(:)
    complex(qp) :: dot 

    dot = p1(1)*p2(1)
    dot = dot - sum(p1(2:)*p2(2:))

  end function dotcc

  function quadc(p) result(quad) 
    complex(qp), intent(in) :: p(:)
    complex(qp) :: quad

    quad = p(1)**2-sum(p(2:)**2)

   end function quadc

  function quadr(p) result(quad) 
    real(qp), intent(in) :: p(:)
    real(qp) :: quad

    quad = p(1)**2-sum(p(2:)**2)

  end function quadr

  function pmassc(p) result(res) 
    complex(qp), intent(in) :: p(:)
    complex(qp) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassc

  function pmassr(p) result(res) 
    real(qp), intent(in) :: p(:)
    real(qp) :: res

    res = sqrt(p(1)**2-sum(p(2:)**2))

  end function pmassr

 
end module qpauxiliary_functions



