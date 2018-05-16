!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECspinors.f90.

! using WEYL REPRESENTATION
module qpspinors
  use types; use consts_qp 
  implicit none 
  private 

  ! --- common data to be used to construct spinors

  real(qp), public :: hc(16)
  ! --- dirac spinors
  complex(qp), public :: u(8,16),uz(8,16),ux(8,16)
  !-----Weyl spinors
  complex(qp), public :: wz(8,16),wx(8,16),wy(8,16)
  complex(qp), public :: bwz(8,16), bwx(8,16), bwy(8,16)


  public :: init_spinors 

  contains 


    subroutine init_spinors

      u = czero;  uz = czero;  ux = czero 

      u(1,1)  = cone; u(2,2)  = cone; u(3,5) = cone 
      u(4,6)  = cone; u(5,9)  = cone; u(6,10)= cone
      u(7,13) = cone; u(8,14) = cone

      uz(1,1)  = cone; uz(1,3)  =  cone; 
      uz(2,2)  = cone; uz(2,4)  = -cone; 
      uz(3,5)  = cone; uz(3,7)  =  cone; 
      uz(4,6)  = cone; uz(4,8)  = -cone;
      uz(5,9)  = cone; uz(5,11) =  cone; 
      uz(6,10) = cone; uz(6,12) = -cone;
      uz(7,13) = cone; uz(7,15) =  cone; 
      uz(8,14) = cone; uz(8,16) = -cone


      ux(1,1)  =  csqrt2; ux(1,2)  =  csqrt2;
      ux(1,3)  =  csqrt2; ux(1,4)  =  csqrt2; 
      ux(2,1)  = -csqrt2; ux(2,2)  =  csqrt2;
      ux(2,3)  =  csqrt2; ux(2,4)  = -csqrt2; 
      ux(3,5)  =  csqrt2; ux(3,6)  =  csqrt2;
      ux(3,7)  =  csqrt2; ux(3,8)  =  csqrt2; 
      ux(4,5)  = -csqrt2; ux(4,6)  =  csqrt2;
      ux(4,7)  =  csqrt2; ux(4,8)  = -csqrt2; 
      ux(5,9)  =  csqrt2; ux(5,10) =  csqrt2; 
      ux(5,11) =  csqrt2; ux(5,12) =  csqrt2; 
      ux(6,9)  = -csqrt2; ux(6,10) =  csqrt2; 
      ux(6,11) =  csqrt2; ux(6,12) = -csqrt2; 
      ux(7,13) =  csqrt2; ux(7,14) =  csqrt2; 
      ux(7,15) =  csqrt2; ux(7,16) =  csqrt2; 
      ux(8,13) = -csqrt2; ux(8,14) =  csqrt2; 
      ux(8,15) =  csqrt2; ux(8,16) = -csqrt2  


      hc(1) =   one; hc(2)  =  one; 
      hc(3)  = -one; hc(4)  = -one; hc(5) = one;  hc(6)  = one; 
      hc(7)  = -one; hc(8)  = -one; hc(9) = one;  hc(10) = one; 
      hc(11) = -one; hc(12) = -one; hc(13) = one; hc(14) = one; 
      hc(15) = -one; hc(16) = -one


      wz = czero; wx = czero; wy = czero 
      bwz = czero; bwx = czero; bwy = czero  


      wz(1,1)  =  cone;   
      wz(2,4)  = -cone; 
      wz(3,5)  =  cone;  
      wz(4,8)  = -cone; 
      wz(5,9)  =  cone;  
      wz(6,12) = -cone;
      wz(7,13) =  cone; 
      wz(8,16) = -cone;

      wx(1,1)  = cone; wx(1,2)  =  cone;  
      wx(2,3)  = cone; wx(2,4)  = -cone;
      wx(3,5)  = cone; wx(3,6)  =  cone; 
      wx(4,7)  = cone; wx(4,8)  = -cone;
      wx(5,9)  = cone; wx(5,10) =  cone; 
      wx(6,11) = cone; wx(6,12) = -cone;
      wx(7,13) = cone; wx(7,14) =  cone;
      wx(8,15) = cone; wx(8,16) = -cone 

      wy(1,1)  = cone; wy(1,2)  =  ci;  
      wy(2,3)  = cone; wy(2,4)  = -ci;
      wy(3,5)  = cone; wy(3,6)  =  ci; 
      wy(4,7)  = cone; wy(4,8)  = -ci;
      wy(5,9)  = cone; wy(5,10) =  ci; 
      wy(6,11) = cone; wy(6,12) = -ci;
      wy(7,13) = cone; wy(7,14) =  ci;
      wy(8,15) = cone; wy(8,16) = -ci; 


      bwz(1,3)  =  cone; 
      bwz(2,2)  = -cone;
      bwz(3,7)  =  cone; 
      bwz(4,6)  = -cone;
      bwz(5,11) =  cone;
      bwz(6,10) = -cone;
      bwz(7,15) =  cone;
      bwz(8,14) = -cone 

      bwx(1,3)  = cone; bwx(1,4)  =  cone;  
      bwx(2,1)  = cone; bwx(2,2)  = -cone; 
      bwx(3,7)  = cone; bwx(3,8)  =  cone; 
      bwx(4,5)  = cone; bwx(4,6)  = -cone; 
      bwx(5,11) = cone; bwx(5,12) =  cone; 
      bwx(6,9)  = cone; bwx(6,10) = -cone; 
      bwx(7,15) = cone; bwx(7,16) =  cone; 
      bwx(8,13) = cone; bwx(8,14) = -cone

      bwy(1,3)  = cone; bwy(1,4)  = -ci;  
      bwy(2,1)  = cone; bwy(2,2)  =  ci; 
      bwy(3,7)  = cone; bwy(3,8)  = -ci; 
      bwy(4,5)  = cone; bwy(4,6)  =  ci; 
      bwy(5,11) = cone; bwy(5,12) = -ci; 
      bwy(6,9)  = cone; bwy(6,10) =  ci; 
      bwy(7,15) = cone; bwy(7,16) = -ci; 
      bwy(8,13) = cone; bwy(8,14) =  ci

    end subroutine init_spinors

end module qpspinors 
