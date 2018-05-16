!======================================================================
! A file containing definitions of constants in a range of precisions
!

!----------------------------------------------------------------------
module consts_dp
  use types
  implicit none
  private

  ! --- numerical constants 
  real(dp), public, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), public, parameter :: zeta2 =&
       & 1.644934066848226436472415166646025189219_dp
  real(dp), public, parameter :: zeta3 =&
       & 1.202056903159594285399738161511449990765_dp
  real(dp), parameter, public :: pisq =&
       & 9.869604401089358618834490999876151135314_dp
  real(dp), parameter, public :: pisqo6 = pisq/6.0_dp
  real(dp), parameter, public :: eulergamma =&
       & 0.577215664901532860606512090082402431042_dp
  real(dp), public, parameter :: half  = 0.5_dp, two = 2.0_dp
  real(dp), public, parameter :: zero  = 0.0_dp, one = 1.0_dp
  real(dp), public, parameter :: mone = -1._dp 
  real(dp), public, parameter :: three = 3._dp 
  real(dp), public, parameter :: four  = 4._dp 
  real(dp), public, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_dp 
  real(dp), public, parameter :: msqrt2 = &
       &-1.4142135623730950488016887242096980785696718753769_dp 
  real(dp), public, parameter :: sqrt3 = &
       &1.7320508075688772935274463415058723669428_dp
  real(dp), public, parameter :: twothirds = &
	&0.666666666666666666666666666666666666666_dp 
  real(dp), public, parameter :: onethird = &
	&0.333333333333333333333333333333333333333_dp 

  complex(dp), parameter, public :: czero = (0.0_dp,0.0_dp)
  complex(dp), parameter, public :: cone = (1.0_dp,0.0_dp)
  complex(dp), parameter, public :: ctwo = (2.0_dp,0.0_dp)
  complex(dp), parameter, public :: chalf = (0.5_dp,0.0_dp)
  complex(dp), parameter, public :: ci=(0.0_dp,1.0_dp)
  complex(dp), parameter, public :: csqrt2 = &
       &(1.4142135623730950488016887242096980785696718753769_dp,0.0_dp)
  complex(dp), parameter, public :: impi = ci*pi



  ! --- QCD constants 
  real(dp), parameter, public :: Nc = three
  real(dp), parameter, public :: Cf = four/three
  real(dp), public :: TR = 5._dp/two 
  real(dp), parameter, public :: nf = 5._dp 
  real(dp), parameter, public :: V = Nc**2-1 

  ! --- numerical tolerance 
  complex(dp), public         ::  iunit 
  real(dp), parameter, public :: tol = 1e-14_dp 
  real(dp), parameter, public :: sq2tol = 1e-7_dp 
  real(dp), parameter, public :: propcut = 1e-10_dp  

  ! -- masses of various particles
  real(dp), parameter, public :: mt=0.0_dp
  real(dp), parameter, public :: mb=0.0_dp
  real(dp), parameter, public :: mw=0.8_dp
  real(dp), public            :: mz ! dynamical mass 
  real(dp), public            :: mur, as_mur 

end module consts_dp
