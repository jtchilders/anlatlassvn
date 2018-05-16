!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECpol_int.f90.

! using WEYL REPRESENTATION
module dppol_int
  use types; use consts_dp; 
  use dpspinors
  use dpaux_functions
  implicit none 
  private

  public :: give_usp, give_barusp

contains


  ! --- spin eigenstates for on-shell momenta
  ! --- Weyl representation
  subroutine give_usp(Nf,Dv,Ds,p,m,r)
    integer, intent(in)      :: Nf,Dv,Ds
    real(dp), intent(in)     :: m
    complex(dp), intent(in)  :: p(:)
    complex(dp), intent(out) :: r(:,:)
    ! -----------------------------------
    complex(dp) :: u1(Ds),f(Ds),fc 
    real(dp)    :: rmax 
    integer :: i, case

    !if (size(p) /= Dv) stop 'give_usp: size(p) /= Dv' 
    !if (size(r,dim=1) /= Ds) stop 'give_usp: size(r,dim=1) /= Ds' 
    !if (size(r,dim=2) /= 8) stop 'give_usp: size(r,dim=2) /= 8' 

    r=czero

    if (abs(m) > 0.0_dp ) &
         &stop 'give_usp: Error: massive Weyl spinors not implemented'

    case = 4
    rmax = abs(p(1)-p(4))
    if (abs(p(1)-p(2)) > rmax) then 
       case = 2
       rmax = abs(p(1)-p(2))
    endif
    if (abs(p(1)-p(3)) > rmax) then 
       case = 3 
       rmax = abs(p(1)-p(3))
    endif
    if (rmax < sq2tol) then 
       write(*,*) 'p', p 
       write(*,*) 'give_usp: massless fermion fails?' 
    endif

    if (case == 2) then 
       fc = sqrt(two*(p(1)-p(2)))
    elseif (case == 3) then 
       fc = sqrt(two*(p(1)-p(3)))
    elseif (case == 4) then 
       fc = sqrt((p(1)-p(4)))
    else
       write(*,*) 'case', case 
       write(*,*) 'give_usp: undefined case' 
    endif


    ! --- massless case
    do i=1,Nf

       if (case == 4) then 
          u1=wz(i,1:Ds)
       elseif (case == 2) then 
          u1=wx(i,1:Ds)
       elseif (case == 3) then 
          u1=wy(i,1:Ds)
       else
          fc = czero 
          u1 = czero 
       endif

       f = spi2(p,u1)
       r(:,i)=f/fc

    enddo


  end subroutine give_usp

  subroutine give_barusp(Nf,Dv,Ds,p,m,r)
    integer, intent(in)      :: Nf,Dv,Ds
    real(dp), intent(in)     :: m
    complex(dp), intent(in)  :: p(:)
    complex(dp), intent(out) :: r(:,:)
    ! --------------------------------------
    complex(dp) :: u1(Ds),f(Ds)
    complex(dp) ::  fc
    real(dp)    :: rmax 
    integer :: i,j, case 

    !if (size(p) /= Dv) stop 'give_barusp: size(p) /= Dv' 
    !if (size(r,dim=1) /= Ds) stop 'give_barusp: size(r,dim=1) /= Ds' 
    !if (size(r,dim=2) /= 8) stop 'give_barusp: size(r,dim=2) /= 8' 
    r = czero

    if (abs(m) > 0.0_dp) stop 'give_barusp: Error: &
         &massive Weyl fermions not implemented'

    case = 4
    rmax = abs(p(1)-p(4))
    if (abs(p(1)-p(2)) > rmax) then 
       case = 2
       rmax = abs(p(1)-p(2))
    endif
    if (abs(p(1)-p(3)) > rmax) then 
       case = 3 
       rmax = abs(p(1)-p(3))
    endif
    if (rmax < sq2tol) then 
       write(*,*) 'p', p 
       write(*,*) 'give_barusp: massless fermion fails?' 
    endif

    if (case == 2) then 
       fc = sqrt(two*(p(1)-p(2)))
    elseif (case == 3) then 
       fc = sqrt(two*(p(1)-p(3)))
    elseif (case == 4) then 
       fc = sqrt((p(1)-p(4)))
    else
       write(*,*) 'case', case 
       write(*,*) 'give_barusp: undefined case' 
    endif


    ! --- massless case
    do i=1,Nf

       if (case == 4) then 
          u1=bwz(i,1:Ds)
       elseif (case == 2) then 
          u1=bwx(i,1:Ds)
       elseif (case == 3) then 
          u1=bwy(i,1:Ds)
       else
          fc = czero 
          u1 = czero 
       endif

       f = spb2(u1,p)
       r(:,i)=f/fc

    enddo


  end subroutine give_barusp


end module dppol_int
