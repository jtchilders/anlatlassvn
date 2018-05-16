!======================================================================
! A file containing definitions of constants in a range of precisions
!

!----------------------------------------------------------------------
module consts_qp
  use types
  implicit none
  private

  ! --- numerical constants 
  real(qp), public, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_qp
  real(qp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_qp
  real(qp), public, parameter :: zeta2 =&
       & 1.644934066848226436472415166646025189219_qp
  real(qp), public, parameter :: zeta3 =&
       & 1.202056903159594285399738161511449990765_qp
  real(qp), parameter, public :: pisq =&
       & 9.869604401089358618834490999876151135314_qp
  real(qp), parameter, public :: eulergamma =&
       & 0.577215664901532860606512090082402431042_qp
  real(qp), public, parameter :: half  = 0.5_qp, two = 2.0_qp
  real(qp), public, parameter :: zero  = 0.0_qp, one = 1.0_qp
  real(qp), public, parameter :: mone = -1._qp 
  real(qp), public, parameter :: three = 3._qp 
  real(qp), public, parameter :: four  = 4._qp 
  real(qp), public, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_qp 
  real(qp), public, parameter :: msqrt2 = &
       &-1.4142135623730950488016887242096980785696718753769_qp 
  real(qp), public, parameter :: sqrt3 = &
       &1.7320508075688772935274463415058723669428_qp

  complex(qp), parameter, public :: czero = (0.0_qp,0.0_qp)
  complex(qp), parameter, public :: cone = (1.0_qp,0.0_qp)
  complex(qp), parameter, public :: ctwo = (2.0_qp,0.0_qp)
  complex(qp), parameter, public :: chalf = (0.5_qp,0.0_qp)
  complex(qp), parameter, public :: ci=(0.0_qp,1.0_qp)
  complex(qp), parameter, public :: csqrt2 = &
       &(1.4142135623730950488016887242096980785696718753769_qp,0.0_qp)


  ! --- QCD constants 
  real(qp), parameter, public :: Nc = 3
  real(qp), parameter, public :: Cf = four/three
  real(qp), public :: TR = 5._qp/two 
  real(qp), parameter, public :: nf = 5._qp 
  real(qp), parameter, public :: V = Nc**2-1 

  ! --- numerical tolerance 
  complex(qp), public         ::  iunit 
  real(qp), parameter, public :: tol = 1e-28_qp 
  real(qp), parameter, public :: sq2tol = 1e-14_qp 
  real(qp), parameter, public :: propcut = 1e-14_qp  

  ! -- masses of various particles
  real(qp), parameter, public :: mt=0.0_qp
  real(qp), parameter, public :: mb=0.0_qp
  real(qp), parameter, public :: mw=0.8_qp
  real(qp), public            :: mz ! dynamical mass 


end module consts_qp
