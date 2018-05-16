!======================================================================
! A file containing MCFM constants
!

!----------------------------------------------------------------------
module consts_MCFM
  use types
  use consts_dp 
  implicit none
  private

  ! --- numerical constants 
  real(dp), public :: MCFMwmass
  real(dp), public :: MCFMwwidth
  real(dp), public :: MCFMzmass
  real(dp), public :: MCFMzwidth
  real(dp), public :: MCFMgw
  real(dp), public :: MCFMgsq
  integer, public  :: MCFMmaxd
  integer, public  :: MCFMndmax
  complex(dp), public :: coupl_gamz(2)  

  character(len=3), public  :: MCFMchn 

  public :: set_coupl_gamz

  contains 

    ! -- sets globallay the coupling to the quark line
    ! 
    subroutine set_coupl_gamz(pZ,itype) 
      use dpaux_functions 
      use define_ampl, only : include_gamZ  
      integer, parameter :: nf = 5 
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      real(dp), intent(in) :: pZ(4) 
      integer, intent(in) :: itype 
      ! -----------------------------------
      complex(dp) :: propZ, pZc(4) 
      integer :: minus,mplus,mp(2)
      data minus,mplus/1,2/
      data mp/-1d0,+1d0/


      pZc = dcmplx(pZ)
      propZ = sc(pZc,pZc)/(sc(pZc,pZc)-MCFMzmass**2 + &
           &(0._dp,1._dp)*MCFMzmass*MCFMzwidth)
      !-- couplings according to hep-ph/9803250 Eqs. 3.4 and 3.6
      if (include_gamZ) then 
         ! Q(1) = -1/3 (down) Q(2)=2/3 (up) 
         ! mp(1) = -1  (down) mp(2)=1 (up) (=2*T3) 
         ! itype = 1 -> down; 2 -> up  ===> 3-itype swaps 1 <-> 2 
         coupl_gamZ(minus)=(mp(3-itype)*2d0*Q(3-itype)*xw+&
              &propZ*(1d0-mp(3-itype)*2d0*Q(3-itype)*xw))
         coupl_gamZ(mplus) = mp(3-itype)*2d0*Q(3-itype)*xw*(1d0-propZ)
      else
         coupl_gamZ = 0d0 ! switch gam/Z off 
      endif

    end subroutine set_coupl_gamz
    


end module consts_MCFM
