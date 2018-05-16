module mpconverter
  use types; use mpmodule; use mpfunmod
  private


  interface mptodp
     module procedure mptodp_fun_0d,   mptodp_fun_1d
     module procedure mptodp_fun_0d_c, mptodp_fun_1d_c
     module procedure mptodp_fun_2d,   mptodp_fun_3d
     module procedure mptodp_fun_2d_c, mptodp_fun_3d_c
  end interface
  interface dptomp
     module procedure dptomp_fun_0d,   dptomp_fun_1d
     module procedure dptomp_fun_2d,   dptomp_fun_3d
  end interface
  
  public :: mptodp, dptomp

contains
  
  function mptodp_fun_0d(mpnumber) result(mptodp_fun)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
         type (mp_real) (q)
    real(dp) :: mptodp_fun
    type(mp_real), intent(in) :: mpnumber
    
    call mpmdc (mpnumber%mpr, da, ia)
    mptodp_fun = da * 2.d0 ** ia
  end function mptodp_fun_0d


  function dptomp_fun_0d(dpnumber) result(dptomp_fun)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
         type (mp_real) (q)
    type(mp_real) :: dptomp_fun
    real(dp), intent(in) :: dpnumber
    
    dptomp_fun = mpreal(real(dpnumber,kind=dp15))
  end function dptomp_fun_0d


  function dptomp_fun_1d(dpnumber)  result(dptomp_fun)
    real(dp)      :: dpnumber(:)
    type(mp_real) :: dptomp_fun(size(dpnumber))
    integer :: I

    do i = 1, size(dpnumber)
       dptomp_fun(i) = dptomp_fun_0d(dpnumber(i))
    end do
  end function dptomp_fun_1d
  

  function mptodp_fun_1d(mpnumber)  result(mptodp_fun)
    type(mp_real), intent(in) :: mpnumber(:)
    real(dp)                  :: mptodp_fun(size(mpnumber))
    integer :: I

    do i = 1, size(mpnumber)
       mptodp_fun(i) = mptodp_fun_0d(mpnumber(i))
    end do
  end function mptodp_fun_1d

  function dptomp_fun_2d(dpnumber)  result(dptomp_fun)
    real(dp)      :: dpnumber(:,:)
    type(mp_real) :: dptomp_fun(size(dpnumber,dim=1),size(dpnumber,dim=2))
    integer :: i1,i2

    do i1 = 1, size(dpnumber,dim=1)
       do i2 = 1, size(dpnumber,dim=2)
          dptomp_fun(i1,i2) = dptomp_fun_0d(dpnumber(i1,i2))
       end do
    end do
  end function dptomp_fun_2d
  

  function mptodp_fun_2d(mpnumber)  result(mptodp_fun)
    type(mp_real), intent(in) :: mpnumber(:,:)
    real(dp)                  :: mptodp_fun(size(mpnumber,dim=1),size(mpnumber,dim=2))
    integer :: i1,i2

    do i1 = 1, size(mpnumber,dim=1)
       do i2 = 1, size(mpnumber,dim=2)
          mptodp_fun(i1,i2) = mptodp_fun_0d(mpnumber(i1,i2))
       end do
    end do

  end function mptodp_fun_2d

  function dptomp_fun_3d(dpnumber)  result(dptomp_fun)
    real(dp)      :: dpnumber(:,:,:)
    type(mp_real) :: dptomp_fun(size(dpnumber,dim=1),&
         &size(dpnumber,dim=2),size(dpnumber,dim=3))
    integer :: i1,i2,i3

    do i1 = 1, size(dpnumber,dim=1)
       do i2 = 1, size(dpnumber,dim=2)
          do i3 = 1, size(dpnumber,dim=3)
             dptomp_fun(i1,i2,i3) = dptomp_fun_0d(dpnumber(i1,i2,i3))
          end do
       end do
    enddo
  end function dptomp_fun_3d
  

  function mptodp_fun_3d(mpnumber)  result(mptodp_fun)
    type(mp_real), intent(in) :: mpnumber(:,:,:)
    real(dp)                  :: mptodp_fun(size(mpnumber,dim=1),&
         &size(mpnumber,dim=2),size(mpnumber,dim=3))
    integer :: i1,i2,i3

    do i1 = 1, size(mpnumber,dim=1)
       do i2 = 1, size(mpnumber,dim=2)
          do i3 = 1, size(mpnumber,dim=3)
             mptodp_fun(i1,i2,i3) = mptodp_fun_0d(mpnumber(i1,i2,i3))
          end do
       end do
    enddo
  end function mptodp_fun_3d




  function mptodp_fun_0d_c(mpnumber) result(mptodp_fun)
    implicit complex (c), double precision (d), type (mp_integer) (j), &
         type (mp_real) (q)
    complex(dp) :: mptodp_fun
    type(mp_complex), intent(in) :: mpnumber
    
    ! intrinsic conversion
    mptodp_fun = mpnumber 
  end function mptodp_fun_0d_c



  function mptodp_fun_1d_c(mpnumber)  result(mptodp_fun)
    type(mp_complex), intent(in) :: mpnumber(:)
    complex(dp)                  :: mptodp_fun(size(mpnumber))
    integer :: I

    do i = 1, size(mpnumber)
       mptodp_fun(i) = mptodp_fun_0d_c(mpnumber(i))
    end do
  end function mptodp_fun_1d_c


  function mptodp_fun_2d_c(mpnumber)  result(mptodp_fun)
    type(mp_complex), intent(in) :: mpnumber(:,:)
    complex(dp)                  :: mptodp_fun(size(mpnumber,dim=1),size(mpnumber,dim=2))
    integer :: i1,i2

    do i1 = 1, size(mpnumber,dim=1)
       do i2 = 1, size(mpnumber,dim=2)
          mptodp_fun(i1,i2) = mptodp_fun_0d_c(mpnumber(i1,i2))
       end do
    end do

  end function mptodp_fun_2d_c


  function mptodp_fun_3d_c(mpnumber)  result(mptodp_fun)
    type(mp_complex), intent(in) :: mpnumber(:,:,:)
    complex(dp)                  :: mptodp_fun(size(mpnumber,dim=1),&
         &size(mpnumber,dim=2),size(mpnumber,dim=3))
    integer :: i1,i2,i3

    do i1 = 1, size(mpnumber,dim=1)
       do i2 = 1, size(mpnumber,dim=2)
          do i3 = 1, size(mpnumber,dim=3)
             mptodp_fun(i1,i2,i3) = mptodp_fun_0d_c(mpnumber(i1,i2,i3))
          end do
       end do
    enddo
  end function mptodp_fun_3d_c


end module mpconverter


