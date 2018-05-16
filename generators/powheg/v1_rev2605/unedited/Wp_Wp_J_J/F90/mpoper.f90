module mpoperations_hidden 
  use mpmodule
  implicit none 
  private 

  ! some assignments  -> = in mpoperations 
  public :: mpsetequal1,mpsetequal2,mpsetequal3,mpsetequal4
  public :: mpfixrealvalue1,mpfixrealvalue2,mpfixrealvalue3,mpfixrealvalue4
  public :: mpfixintegervalue1,mpfixintegervalue2,mpfixintegervalue3
  public :: mpfixintegervalue4
  public :: mpfixintegerarray1,mpfixintegerarray2,mpfixintegerarray3
  public :: mpfixintegerarray4
  public :: mpfixmpvalue1, mpfixmpvalue2, mpfixmpvalue3, mpfixmpvalue4


  ! some elementary operations: +,-,/,* in mpoperations 
  public :: mpscalproda, mpscalprodb
  public :: mpscalproda2, mpscalprodb2
  public :: mpintsuma2, mpintsumb2
  public :: mpscalsuma, mpscalsumb,mpscalsuma2, mpscalsumb2
  public :: mpscaldifa, mpscaldifb
  public :: mpscaldiv1, mpscaldiv2
  public :: mpplus1,mpplus2,mpplus3
  public :: mpplusintegerarray1a,mpplusintegerarray2a
  public :: mpplusintegerarray3a,mpplusintegerarray1b
  public :: mpplusintegerarray2b,mpplusintegerarray3b
  public :: mpminus1,mpminus2,mpminus3
  public :: mpminusintegerarray1a,mpminusintegerarray2a
  public :: mpminusintegerarray3a,mpminusintegerarray1b
  public :: mpminusintegerarray2b,mpminusintegerarray3b
  public :: mpminusu1 
  public :: mptimes1,mptimes2,mptimes3
  public :: mpdivision1,mpdivision2,mpdivision3
  public :: mpintdivision1,mpintdivision2,mpintdivision3
  public :: mpexpo
  public :: mpsquareroot1,mpsquareroot2,mpsquareroot3
  public :: mplog1,mplog2,mplog3
  public :: mpexp1,mpexp2,mpexp3
  public :: mpabs1,mpabs2,mpabs3
  public :: mpsin1,mpsin2,mpsin3
  public :: mpcos1,mpcos2,mpcos3
  public :: mpatan21,mpatan22,mpatan23
  public :: mpnint1,mpnint2,mpnint3
  public :: mpsign1, mpsign01, mpsign2, mpsign3

  !    == , < , > <=, >= 
  !  public :: mpequal,mpequal1,mpequal2,mpequal3
  ! AB
  public :: mpequal10,mpequal11,mpequal2,mpequal3
  !  public :: mpequal,mpequal2,mpequal3
  public :: mplessthen1, mplessthen2, mplessthen3
  public :: mplessequal1,mplessequal2,mplessequal3
  public :: mpgreaterthen1,mpgreaterthen2,mpgreaterthen3
  public :: mpgreaterequal1,mpgreaterequal2,mpgreaterequal3


contains 

  ! set a constant mp real value for an mp matrix
  subroutine mpfixmpvalue1(a,b) 
    type (mp_real),intent(in) ::  b
    type (mp_real),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       a(i)=b
    enddo
  end subroutine mpfixmpvalue1

  subroutine mpfixmpvalue2(a,b) 
    type (mp_real),intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       end do
    enddo
  end subroutine mpfixmpvalue2

  subroutine mpfixmpvalue3(a,b) 
    type (mp_real),intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:)
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             a(i,j,k)=b
          end do
       end do
    end do
  end subroutine mpfixmpvalue3

  subroutine mpfixmpvalue4(a,b) 
    type (mp_real),intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:,:)
    integer :: i,j,k,l 


    do l=lbound(a,dim=4),ubound(a,dim=4)
       do k=lbound(a,dim=3),ubound(a,dim=3)
          do j=lbound(a,dim=2),ubound(a,dim=2)
             do i=lbound(a,dim=1),ubound(a,dim=1)
                a(i,j,k,l)=b
             end do
          end do
       end do
    end do
  end subroutine mpfixmpvalue4

  subroutine mpfixrealvalue1(a,b) 
    real,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       a(i)=b
    enddo
  end subroutine mpfixrealvalue1

  subroutine mpfixrealvalue2(a,b) 
    real,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       end do
    end do
  end subroutine mpfixrealvalue2

  subroutine mpfixrealvalue3(a,b) 
    real,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:)
    integer :: i,j,k

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             a(i,j,k)=b
          end do
       end do
    end do
  end subroutine mpfixrealvalue3

  subroutine mpfixrealvalue4(a,b) 
    real,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:,:)
    integer :: i,j,k,l 

    do l=lbound(a,dim=4),ubound(a,dim=4)
       do k=lbound(a,dim=3),ubound(a,dim=3)
          do j=lbound(a,dim=2),ubound(a,dim=2)
             do i=lbound(a,dim=1),ubound(a,dim=1)
                a(i,j,k,l)=b
             end do
          end do
       end do
    end do
  end subroutine mpfixrealvalue4


  subroutine mpfixintegervalue1(a,b) 
    integer,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)      
       a(i)=b
    enddo
  end subroutine mpfixintegervalue1

  subroutine mpfixintegervalue2(a,b) 
    integer,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       enddo
    end do
  end subroutine mpfixintegervalue2

  subroutine mpfixintegervalue3(a,b) 
    integer,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:)
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)      
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             a(i,j,k)=b
          enddo
       end do
    end do
  end subroutine mpfixintegervalue3

  subroutine mpfixintegervalue4(a,b) 
    integer,intent(in) ::  b
    type (mp_real),intent(inout) :: a(:,:,:,:)
    integer :: i,j,k,l 

    do l=lbound(a,dim=4),ubound(a,dim=4)
       do k=lbound(a,dim=3),ubound(a,dim=3)
          do j=lbound(a,dim=2),ubound(a,dim=2)
             do i=lbound(a,dim=1),ubound(a,dim=1)
                a(i,j,k,l)=b
             enddo
          end do
       end do
    end do
  end subroutine mpfixintegervalue4

  subroutine mpfixintegerarray1(a,b) 
    integer,intent(in) ::  b(:)
    type (mp_real),intent(inout) :: a(:)
    integer :: i 

    !if (size(a).ne.size(b)) stop 'Size 1 differ in  mpfixintegerarray'
    do i=lbound(a,dim=1),ubound(a,dim=1)      
       a(i)=b(i)
    enddo
  end subroutine mpfixintegerarray1

  subroutine mpfixintegerarray2(a,b) 
    integer,intent(in) ::  b(:,:)
    type (mp_real),intent(inout) :: a(:,:)
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpfixintegerarry' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpfixintegerarry' 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b(i,j)
       enddo
    end do
  end subroutine mpfixintegerarray2

  subroutine mpfixintegerarray3(a,b) 
    integer,intent(in) ::  b(:,:,:)
    type (mp_real),intent(inout) :: a(:,:,:)
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpfixintegerarry' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpfixintegerarry' 
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpfixintegerarry' 
    do k=lbound(a,dim=3),ubound(a,dim=3)      
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             a(i,j,k)=b(i,j,k)
          enddo
       end do
    end do
  end subroutine mpfixintegerarray3

  subroutine mpfixintegerarray4(a,b) 
    integer,intent(in) ::  b(:,:,:,:)
    type (mp_real),intent(inout) :: a(:,:,:,:)
    integer :: i,j,k,l 


    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpfixintegerarry' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpfixintegerarry' 
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpfixintegerarry' 
    !if (size(a,dim=4).ne.size(b,dim=4)) stop 'Size 4 differ in mpfixintegerarry' 
    do l=lbound(a,dim=4),ubound(a,dim=4)
       do k=lbound(a,dim=3),ubound(a,dim=3)
          do j=lbound(a,dim=2),ubound(a,dim=2)
             do i=lbound(a,dim=1),ubound(a,dim=1)
                a(i,j,k,l)=b(i,j,k,l)
             enddo
          end do
       end do
    end do
  end subroutine mpfixintegerarray4

  ! set two mparrays equal 
  subroutine mpsetequal1(a, b) 
    type (mp_real),intent(in) ::  b(:)
    type (mp_real), intent(inout) :: a(size(b,dim=1))
    integer :: i 

    do i=lbound(b,dim=1),ubound(b,dim=1) 
       a(i)=b(i)
    enddo
  end subroutine mpsetequal1

  subroutine mpsetequal2(a, b) 
    type (mp_real),intent(in) ::  b(:,:)
    type (mp_real), intent(inout) :: a(size(b,dim=1),size(b,dim=2))
    integer :: i,j 

    do j=lbound(b,dim=2),ubound(b,dim=2) 
       do i=lbound(b,dim=1),ubound(b,dim=1) 
          a(i,j)=b(i,j)
       end do
    end do
  end subroutine mpsetequal2

  subroutine mpsetequal3(a, b) 
    type (mp_real),intent(in) ::  b(:,:,:)
    type (mp_real), intent(inout) :: a(size(b,dim=1),size(b,dim=2),size(b,dim=3))

    integer :: i,j,k 

    do k=lbound(b,dim=3),ubound(b,dim=3) 
       do j=lbound(b,dim=2),ubound(b,dim=2) 
          do i=lbound(b,dim=1),ubound(b,dim=1) 
             a(i,j,k)=b(i,j,k)
          end do
       end do
    end do
  end subroutine mpsetequal3

  subroutine mpsetequal4(a, b) 
    type (mp_real),intent(in) ::  b(:,:,:,:)
    type (mp_real), intent(inout) :: a(size(b,dim=1),size(b,dim=2),size(b,dim=3),size(b,dim=4))
    integer :: i,j,k,l 

    do l=lbound(b,dim=4),ubound(b,dim=4) 
       do k=lbound(b,dim=3),ubound(b,dim=3) 
          do j=lbound(b,dim=2),ubound(b,dim=2) 
             do i=lbound(b,dim=1),ubound(b,dim=1) 
                a(i,j,k,l)=b(i,j,k,l)
             end do
          end do
       end do
    end do
  end subroutine mpsetequal4

  ! scalar multiplication 
  function mpscalproda(a,b) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real)            ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    end do
  end function mpscalproda


  function mpscalprodb(b,a) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    enddo
  end function mpscalprodb

  function mpscalproda2(a,b) result(c)
    type (mp_real),intent(in) :: a(:,:), b
    type (mp_real) ::  c(lbound(a,dim=1):ubound(a,dim=1),&
         &lbound(a,dim=2):ubound(a,dim=2))
    integer :: i,j

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)*b
       end do
    end do
  end function mpscalproda2

  function mpscalprodb2(b,a) result(c)
    type (mp_real),intent(in) :: a(:,:), b
    type (mp_real) ::  c(lbound(a,dim=1):ubound(a,dim=1),&
         &lbound(a,dim=2):ubound(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)*b
       enddo
    end do
  end function mpscalprodb2

  function mpscaldiv1(a,b) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(ubound(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)/b
    enddo
  end function mpscaldiv1

  function mpscaldiv2(a,b) result(c)
    type (mp_real),intent(in) :: a(:,:), b
    type (mp_real) ::  c(ubound(a,dim=1),ubound(a,dim=2))
    integer :: i, j

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)/b
       end do
    enddo
  end function mpscaldiv2

  ! scalar sum
  function mpscalsuma(a,b) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)+b
    enddo
  end function mpscalsuma


  function mpscalsumb(b,a) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=b+a(i)
    enddo
  end function mpscalsumb

  function mpscalsuma2(a,b) result(c)
    type (mp_real),intent(in) :: a(:,:), b
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       enddo
    end do
  end function mpscalsuma2


  function mpscalsumb2(b,a) result(c)
    type (mp_real),intent(in) :: a(:,:), b
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       end do
    enddo
  end function mpscalsumb2

  function mpintsuma2(a,b) result(c)
    type (mp_real),intent(in) :: a(:,:)
    integer, intent(in) :: b
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       enddo
    end do
  end function mpintsuma2


  function mpintsumb2(b,a) result(c)
    type (mp_real),intent(in) :: a(:,:)
    integer, intent(in) :: b
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       end do
    enddo
  end function mpintsumb2


  ! scalar difference
  function mpscaldifa(a,b) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)-b
    enddo
  end function mpscaldifa


  function mpscaldifb(b,a) result(c)
    type (mp_real),intent(in) :: a(:), b
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=b-a(i)
    enddo
  end function mpscaldifb


  ! some binary operators with pure vectors now
  function mpplus1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:),b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1) /= size(b,dim=1) ) stop "Sizes differ in mpplus"

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)+b(i) 
    enddo
  end function mpplus1

  function mpplus2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:),b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1) /= size(b,dim=1) ) stop "Sizes differ in mpplus"
    !if (size(a,dim=2) /= size(b,dim=2) ) stop "Sizes differ in mpplus"

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)+b(i,j)
       enddo
    enddo
  end function mpplus2

  function mpplus3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1) /= size(b,dim=1) ) stop "Size 1 differ in mpplus"
    !if (size(a,dim=2) /= size(b,dim=2) ) stop "Size 2 differ in mpplus"
    !if (size(a,dim=3) /= size(b,dim=3) ) stop "Size 3 differ in mpplus"

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=a(i,j,k)+b(i,j,k)
          enddo
       enddo
    end do
  end function mpplus3


  function mpminusu1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=-a(i) 
    enddo
  end function mpminusu1

  function mpminus1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:),b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminus1

  function mpminus2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:),b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpminus' 


    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)-b(i,j) 
       end do
    end do
  end function mpminus2

  function mpminus3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpminus'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpminus' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)-b(i,j,k) 
          end do
       end do
    end do
  end function mpminus3

  function mpminusintegerarray1a(a,b) result(c)
    type (mp_real),intent(in) ::  b(:)
    integer,intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminusintegerarray1a

  function mpminusintegerarray1b(a,b) result(c)
    type (mp_real),intent(in) ::  a(:)
    integer,intent(in) ::  b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminusintegerarray1b

  function mpminusintegerarray2a(a,b) result(c)
    type (mp_real),intent(in) ::  b(:,:)
    integer,intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 1 differ in mpminus' 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)-b(i,j) 
       enddo
    end do
  end function mpminusintegerarray2a

  function mpminusintegerarray2b(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    integer,intent(in) ::  b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 1 differ in mpminus' 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)-b(i,j) 
       enddo
    end do
  end function mpminusintegerarray2b


  function mpminusintegerarray3a(a,b) result(c)
    integer,intent(in) ::  a(:,:,:)
    type (mp_real),intent(in) ::  b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpminus'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpminus' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)-b(i,j,k) 
          end do
       end do
    end do
  end function mpminusintegerarray3a

  function mpminusintegerarray3b(a,b) result(c)
    integer,intent(in) ::  b(:,:,:)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpminus'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpminus' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)-b(i,j,k) 
          end do
       end do
    end do
  end function mpminusintegerarray3b


  function mpplusintegerarray1a(a,b) result(c)
    type (mp_real),intent(in) ::  b(:)
    integer,intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpplusintegerarray1a

  function mpplusintegerarray1b(a,b) result(c)
    type (mp_real),intent(in) ::  a(:)
    integer,intent(in) ::  b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpplusintegerarray1b

  function mpplusintegerarray2a(a,b) result(c)
    type (mp_real),intent(in) ::  b(:,:)
    integer,intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 1 differ in mpplus' 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)-b(i,j) 
       enddo
    end do
  end function mpplusintegerarray2a

  function mpplusintegerarray2b(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    integer,intent(in) ::  b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 1 differ in mpplus' 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)-b(i,j) 
       enddo
    end do
  end function mpplusintegerarray2b


  function mpplusintegerarray3a(a,b) result(c)
    integer,intent(in) ::  a(:,:,:)
    type (mp_real),intent(in) ::  b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpplus'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpplus' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)-b(i,j,k) 
          end do
       end do
    end do
  end function mpplusintegerarray3a

  function mpplusintegerarray3b(a,b) result(c)
    integer,intent(in) ::  b(:,:,:)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpplus'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpplus' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)-b(i,j,k) 
          end do
       end do
    end do
  end function mpplusintegerarray3b



  function mptimes1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:),b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes1' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)*b(i) 
    enddo
  end function mptimes1

  function mptimes2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:),b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes2' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mptimes2'

    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)*b(i,j) 
       end do
    end do
  end function mptimes2

  function mptimes3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes3' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mptimes3'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mptimes3' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)*b(i,j,k) 
          end do
       end do
    end do
  end function mptimes3

  function mpdivision1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:),b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdivisions' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)/b(i) 
    enddo
  end function mpdivision1

  function mpdivision2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:),b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdivision' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpdivision'

    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)/b(i,j) 
       end do
    end do
  end function mpdivision2

  function mpdivision3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdivision' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mpdivision'
    !if (size(a,dim=3).ne.size(b,dim=3)) stop 'Size 3 differ in mpdivision' 


    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)/b(i,j,k) 
          end do
       end do
    end do
  end function mpdivision3

  function mpintdivision1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:)
    integer, intent(in) :: b 
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)/b
    enddo
  end function mpintdivision1

  function mpintdivision2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    integer, intent(in) :: b 
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 


    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)/b
       end do
    end do
  end function mpintdivision2

  function mpintdivision3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    integer, intent(in) :: b 
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1) ,ubound(a,dim=1) 
             c(i,j,k)=a(i,j,k)/b
          end do
       end do
    end do
  end function mpintdivision3


  function mpexpo(a,n) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer, intent(in) :: n 
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)**n
    enddo

  end function mpexpo


  function mpsquareroot1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i

    do i=lbound(a,dim=1),ubound(a,dim=1)
       if (a(i)<mpreal('0.0')) stop 'Negative argument of squareroot'
       c(i)=sqrt(a(i))
    enddo

  end function mpsquareroot1

  function mpsquareroot2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j)<mpreal('0.0')) stop 'Negative argument of squareroot'
          c(i,j)=sqrt(a(i,j))
       enddo
    end do
  end function mpsquareroot2

  function mpsquareroot3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k)<mpreal('0.0')) stop 'Negative argument of squareroot'
             c(i,j,k)=sqrt(a(i,j,k))
          enddo
       end do
    end do

  end function mpsquareroot3

  function mplog1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       if (a(i)<mpreal('0.0')) stop 'Negative argument of Log'
       c(i)=log(a(i))
    enddo

  end function mplog1

  function mplog2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j)<mpreal('0.0')) stop 'Negative argument of Log'
          c(i,j)=log(a(i,j))
       enddo
    end do
  end function mplog2

  function mplog3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k)<mpreal('0.0')) stop 'Negative argument of Log'
             c(i,j,k)=log(a(i,j,k))
          enddo
       end do
    end do

  end function mplog3

  function mpexp1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=exp(a(i))
    enddo

  end function mpexp1

  function mpexp2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=exp(a(i,j))
       enddo
    end do
  end function mpexp2

  function mpexp3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=exp(a(i,j,k))
          enddo
       end do
    end do

  end function mpexp3

  function mpabs1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=abs(a(i))
    enddo

  end function mpabs1

  function mpabs2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=abs(a(i,j))
       enddo
    end do
  end function mpabs2

  function mpabs3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=abs(a(i,j,k))
          enddo
       end do
    end do

  end function mpabs3


  function mpsin1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=sin(a(i))
    enddo

  end function mpsin1

  function mpsin2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=sin(a(i,j))
       enddo
    end do
  end function mpsin2

  function mpsin3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=sin(a(i,j,k))
          enddo
       end do
    end do

  end function mpsin3

  function mpcos1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=cos(a(i))
    enddo

  end function mpcos1

  function mpcos2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=cos(a(i,j))
       enddo
    end do
  end function mpcos2

  function mpcos3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=cos(a(i,j,k))
          enddo
       end do
    end do

  end function mpcos3

  function mpatan21(a,b) result(c)
    type (mp_real),intent(in) ::  a(:),b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=atan2(a(i),b(i))
    enddo

  end function mpatan21

  function mpatan22(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:),b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=atan2(a(i,j),b(i,j))
       enddo
    end do
  end function mpatan22

  function mpatan23(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=atan2(a(i,j,k),b(i,j,k))
          enddo
       end do
    end do

  end function mpatan23

  function mpnint1(a) result(c)
    type (mp_real),intent(in) ::  a(:)
    integer ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=nint(a(i))
    enddo

  end function mpnint1

  function mpnint2(a) result(c)
    type (mp_real),intent(in) ::  a(:,:)
    integer ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=nint(a(i,j))
       enddo
    end do
  end function mpnint2

  function mpnint3(a) result(c)
    type (mp_real),intent(in) ::  a(:,:,:)
    integer ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=nint(a(i,j,k))
          enddo
       end do
    end do

  end function mpnint3


  function mpequal11 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b(:)
    logical :: t_f(size(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i).ne.b(i)) t_f(i)= .false.
    end do

  end function mpequal11
  function mpequal10 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i).ne.b) t_f(i)= .false.
    end do

  end function mpequal10

  function mpequal2 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    t_f = .true.                      

    do j = lbound(a,dim=2),ubound(a,dim=2)
       do i = lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j).ne.b) t_f(i,j)= .false.
       end do
    end do
  end function mpequal2

  function mpequal3 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 
    t_f = .true.                      

    do k = lbound(a,dim=3),ubound(a,dim=3)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          do i = lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k).ne.b) t_f(i,j,k)= .false.
          end do
       end do
    end do
  end function mpequal3


  function mplessthen1 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b      
    integer :: i 
    logical :: t_f(lbound(a,dim=1):ubound(a,dim=1))
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i) .ge. b) t_f(i)= .false.
    end do

  end function mplessthen1

  function mplessthen2 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:)       
    type (mp_real), intent(in) :: b      
    integer :: i,j  
    logical :: t_f(size(a,dim=1),size(a,dim=2))
    t_f = .true.                      

    do j = lbound(a,dim=2),ubound(a,dim=2)
       do i = lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j) .ge. b) t_f(i,j)= .false.
       end do
    end do
  end function mplessthen2

  function mplessthen3 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:,:)       
    type (mp_real), intent(in) :: b      
    integer :: i,j,k  
    logical :: t_f(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    t_f = .true.                      

    do k = lbound(a,dim=3),ubound(a,dim=3)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          do i = lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k) .ge. b) t_f(i,j,k)= .false.
          end do
       end do
    end do
  end function mplessthen3

  function mplessequal1 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(lbound(a,dim=1):ubound(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i) .gt. b) t_f(i)= .false.
    end do

  end function mplessequal1
  function mplessequal2 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    t_f = .true.                      

    do j = lbound(a,dim=2),ubound(a,dim=2)
       do i = lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j) .gt. b) t_f(i,j)= .false.
       end do
    end do
  end function mplessequal2

  function mplessequal3 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 
    t_f = .true.                      

    do k = lbound(a,dim=3),ubound(a,dim=3)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          do i = lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k) .gt. b) t_f(i,j,k)= .false.
          end do
       end do
    end do
  end function mplessequal3



  function mpgreaterthen1 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b      
    integer :: i 
    logical :: t_f(lbound(a,dim=1):ubound(a,dim=1))
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i) .le. b) t_f(i)= .false.
    end do

  end function mpgreaterthen1

  function mpgreaterthen2 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:)       
    type (mp_real), intent(in) :: b      
    integer :: i,j  
    logical :: t_f(size(a,dim=1),size(a,dim=2))
    t_f = .true.                      

    do j = lbound(a,dim=2),ubound(a,dim=2)
       do i = lbound(a,dim=1),ubound(a,dim=1)
          if (a(i,j) .le. b) t_f(i,j)= .false.
       end do
    end do
  end function mpgreaterthen2

  function mpgreaterthen3 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:,:)       
    type (mp_real), intent(in) :: b      
    integer :: i,j,k  
    logical :: t_f(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    t_f = .true.                      

    do k = lbound(a,dim=3),ubound(a,dim=3)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          do i = lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k) .le. b) t_f(i,j,k)= .false.
          end do
       end do
    end do
  end function mpgreaterthen3

  function mpgreaterequal1 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(lbound(a,dim=1):ubound(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i) .lt. b) t_f(i)= .false.
    end do

  end function mpgreaterequal1
  function mpgreaterequal2 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          if (a(i,j) .lt. b) t_f(i,j)= .false.
       end do
    end do
  end function mpgreaterequal2

  function mpgreaterequal3 (a, b) result (t_f) 
    type (mp_real), intent(in) :: a(:,:,:)       
    type (mp_real), intent(in) :: b      
    logical :: t_f(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 
    t_f = .true.                      

    do k = lbound(a,dim=3),ubound(a,dim=3)
       do j = lbound(a,dim=2),ubound(a,dim=2)
          do i = lbound(a,dim=1),ubound(a,dim=1)
             if (a(i,j,k) .lt. b) t_f(i,j,k)= .false.
          end do
       end do
    end do
  end function mpgreaterequal3


  function mpsign1(a,b) result(c)
    type (mp_real),intent(in) ::  a(:), b(:)
    type (mp_real) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=sign(a(i),b(i))
    enddo

  end function mpsign1

  function mpsign01(a,b) result(c)
    type (mp_real),intent(in) ::  a, b(:)
    type (mp_real) ::  c(size(b,dim=1))
    integer :: i,j 

    do i=lbound(b,dim=1),ubound(b,dim=1)
       c(i)=sign(a,b(i))
    enddo
  end function mpsign01

  function mpsign2(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:), b(:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=sign(a(i,j),b(i,j))
       enddo
    end do
  end function mpsign2

  function mpsign3(a,b) result(c)
    type (mp_real),intent(in) ::  a(:,:,:), b(:,:,:)
    type (mp_real) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=sign(a(i,j,k),b(i,j,k))
          enddo
       end do
    end do
  end function mpsign3

end module mpoperations_hidden

module mpsimpleoperations
  use mpoperations_hidden  
  implicit none 

! NOTE: Here everything is public 

 interface operator (+)
   module procedure mpplus1,mpplus2,mpplus3,mpscalsuma,mpscalsumb
   module procedure mpscalsuma2,mpscalsumb2, mpintsuma2, mpintsumb2
   module procedure mpplusintegerarray1a,mpplusintegerarray2a
   module procedure mpplusintegerarray3a,mpplusintegerarray1b
   module procedure mpplusintegerarray2b,mpplusintegerarray3b
 end interface
 public :: operator (+)
 private mpplus1,mpplus2,mpplus3,mpscalsuma,mpscalsumb,mpscalsuma2,&
      &mpscalsumb2,mpintsuma2, mpintsumb2,mpplusintegerarray1a,&
      &mpplusintegerarray2a,mpplusintegerarray3a,mpplusintegerarray1b,&
      &mpplusintegerarray2b,mpplusintegerarray3b

  interface operator(-)
     module procedure mpminus1,mpminus2,mpminus3,mpscaldifa,mpscaldifb
     module procedure mpminusintegerarray1a,mpminusintegerarray2a
     module procedure mpminusintegerarray3a,mpminusintegerarray1b
     module procedure mpminusintegerarray2b,mpminusintegerarray3b
     module procedure mpminusu1
  end interface            
  public :: operator(-)
  private mpminus1,mpminus2,mpminus3,mpscaldifa,mpscaldifb,&
       &mpminusintegerarray1a,mpminusintegerarray2a,mpminusintegerarray3a,&
       &mpminusintegerarray1b,mpminusintegerarray2b,mpminusintegerarray3b,&
       &mpminusu1

  interface operator(*)
     module procedure mptimes1, mptimes2,mptimes3
     module procedure mpscalproda, mpscalprodb, mpscalproda2, mpscalprodb2
  end interface            
  public :: operator(*)
  private mptimes1, mptimes2,mptimes3, mpscalproda, mpscalprodb,&
       &mpscalproda2, mpscalprodb2

  interface operator(/)
     module procedure mpdivision1,mpdivision2,mpdivision3,mpscaldiv1, &
          & mpscaldiv2,mpintdivision1,mpintdivision2,mpintdivision3
  end interface            
  public :: operator(/)
  private mpdivision1,mpdivision2,mpdivision3,mpscaldiv1,mpscaldiv2
  private  mpintdivision1,mpintdivision2,mpintdivision3

  interface operator(==)
     module procedure mpequal10,mpequal11,mpequal2,mpequal3
  end interface
  public :: operator(==)
  private mpequal10,mpequal11,mpequal2,mpequal3

  interface operator(<)
     module procedure mplessthen1,mplessthen2,mplessthen3
  end interface
  public :: operator(<)
  private mplessthen1,mplessthen2,mplessthen3

  interface operator(<=)
     module procedure mplessequal1,mplessequal2,mplessequal3
  end interface
  public :: operator(<=)
  private mplessequal1,mplessequal2,mplessequal3

  interface operator(>)
     module procedure mpgreaterthen1,mpgreaterthen2,mpgreaterthen3
  end interface
  public :: operator(>)
  private mpgreaterthen1,mpgreaterthen2,mpgreaterthen3

  interface operator(>=)
     module procedure mpgreaterequal1,mpgreaterequal2,mpgreaterequal3
  end interface
  public :: operator(>=)
  private mpgreaterequal1,mpgreaterequal2,mpgreaterequal3

  interface operator(**)
     module procedure mpexpo
  end interface            
  public :: operator(**)
  private mpexpo

  interface sqrt
     module procedure mpsquareroot1,mpsquareroot2,mpsquareroot3
  end interface            
  public :: sqrt
  private mpsquareroot1,mpsquareroot2,mpsquareroot3

  interface log
     module procedure mplog1,mplog2,mplog3
  end interface
  public :: log
  private mplog1,mplog2,mplog3 

  interface exp
     module procedure mpexp1,mpexp2,mpexp3
  end interface
  public :: exp
  private mpexp1,mpexp2,mpexp3 

  interface abs
     module procedure mpabs1,mpabs2,mpabs3
  end interface
  public :: abs
  private mpabs1,mpabs2,mpabs3 


  interface sin
     module procedure mpsin1,mpsin2,mpsin3
  end interface
  public :: sin
  private mpsin1,mpsin2,mpsin3 

  interface cos
     module procedure mpcos1,mpcos2,mpcos3
  end interface
  public :: cos
  private mpcos1,mpcos2,mpcos3

  interface atan2
     module procedure mpatan21,mpatan22,mpatan23
  end interface
  public :: atan2
  private mpatan21,mpatan22,mpatan23

  interface sign
     module procedure mpsign1,mpsign01,mpsign2,mpsign3
  end interface
  public :: sign
  private mpsign1,mpsign01,mpsign2,mpsign3

  interface nint
     module procedure mpnint1,mpnint2,mpnint3
  end interface
  public :: nint
  private mpnint1,mpnint2,mpnint3

 
  interface assignment (=)
     module procedure mpfixmpvalue1, mpfixmpvalue2, mpfixmpvalue3, &
          &mpfixmpvalue4,mpfixintegervalue1,mpfixintegervalue2,&
          &mpfixintegervalue3,mpfixintegervalue4,mpfixrealvalue1,&
          &mpfixrealvalue2,mpfixrealvalue3,mpfixrealvalue4, &
          &mpfixintegerarray1,mpfixintegerarray2,mpfixintegerarray3,&
          &mpfixintegerarray4
  end interface            
  public :: assignment (=)
  private mpfixmpvalue1, mpfixmpvalue2,mpfixmpvalue3,mpfixmpvalue4
  private mpfixintegervalue1,mpfixintegervalue2,mpfixintegervalue3,&
       &mpfixintegervalue4
  private mpfixrealvalue1,mpfixrealvalue2,mpfixrealvalue3,mpfixrealvalue4
  private mpsetequal1,mpsetequal2,mpsetequal3,mpsetequal4
  private mpfixintegerarray1,mpfixintegerarray2,mpfixintegerarray3,&
       &mpfixintegerarray4

end module mpsimpleoperations

module mpadvancedoperations  
  use mpmodule
  implicit none 
  private 

! stupid routine to print (:), (:,:) and (:,:,:) mpmatrices  
  interface mpwritearray
     module procedure mpwritearray1,mpwritearray2,mpwritearray3
  end interface            
  public :: mpwritearray

  interface max
     module procedure mpmax1,mpmax2,mpmax3
  end interface
  public :: max

  interface maxloc
     module procedure mpmaxlocation1,mpmaxlocation2,mpmaxlocation3
  end interface
  public :: maxloc

  interface minloc
     module procedure mpminlocation1,mpminlocation2,mpminlocation3 
  end interface
  public :: minloc

  interface maxval
     module procedure mpmaxval1,mpmaxval2,mpmaxval3
  end interface
  public :: maxval

  interface minval
     module procedure mpminval1,mpminval2,mpminval3
  end interface
  public :: minval

  interface dot_product
     module procedure mpdot_product
  end interface
  public :: dot_product

 interface sum
     module procedure mpsum, mpsum_2ddim
  end interface            
  public :: sum 

  interface floor
     module procedure mpfloor
  end interface
  public :: floor

  interface ceiling
     module procedure mpceiling
  end interface
  public :: ceiling

  contains 

    function mpmax1(a,b) result(c)
      type(mp_real), intent(in) :: a(:),b(:)
      type(mp_real)             :: c(size(a))
      integer :: i

      do i=1,size(a)
         c(i) = max(a(i),b(i))
      end do

    end function mpmax1

    function mpmax2(a,b) result(c)
      type(mp_real), intent(in) :: a(:,:),b(:,:)
      type(mp_real)             :: c(size(a,1),size(a,2))
      integer :: i,j

         do j=1,size(a,2)
      do i=1,size(a,1)
            c(i,j) = max(a(i,j),b(i,j))
         end do
      end do

    end function mpmax2

    function mpmax3(a,b) result(c)
      type(mp_real), intent(in) :: a(:,:,:),b(:,:,:)
      type(mp_real)             :: c(size(a,1),size(a,2),size(a,3))
      integer :: i,j,k

            do k=1,size(a,3)
         do j=1,size(a,2)
      do i=1,size(a,1)
               c(i,j,k) = max(a(i,j,k),b(i,j,k))
            end do
         end do
      end do

    end function mpmax3

    function mpmaxlocation1(x) result(imax)
      type(mp_real), intent(in):: x(:)
      type(mp_real) :: maxvalue
      integer :: imax(1)
      integer :: i 

      imax = lbound(x,dim=1)
      maxvalue = x(lbound(x,dim=1))
      do i=lbound(x,dim=1), ubound(x,dim=1)
         if (x(i).gt.maxvalue) then  
            imax = i
            maxvalue = x(i)
         endif
      end do

    end function mpmaxlocation1

    function mpmaxlocation2(x) result(imax)
      type(mp_real), intent(in):: x(:,:)
      type(mp_real) :: maxvalue
      integer :: i,j 
      integer :: imax(2)

      imax(1)= lbound(x,dim=1)
      imax(2)= lbound(x,dim=2)
      maxvalue = x(lbound(x,dim=1),lbound(x,dim=2))
         do j=lbound(x,dim=2), ubound(x,dim=2)
      do i=lbound(x,dim=1), ubound(x,dim=1)
            if (x(i,j).lt.maxvalue) then  
               imax(1) = i
               imax(2) = j
               maxvalue = x(i,j)
            endif
         end do
      end do
    end function mpmaxlocation2

  
    function mpmaxlocation3(x) result(imax)
      type(mp_real), intent(in):: x(:,:,:)
      type(mp_real) :: maxvalue
      integer :: i,j,k 
      integer :: imax(3)

      imax(1)= lbound(x,dim=1)
      imax(2)= lbound(x,dim=2)
      imax(3)= lbound(x,dim=3)
      maxvalue = x(lbound(x,dim=1),lbound(x,dim=2),lbound(x,dim=3))
            do k=lbound(x,dim=3), ubound(x,dim=3)
         do j=lbound(x,dim=2), ubound(x,dim=2)
      do i=lbound(x,dim=1), ubound(x,dim=1)
               if (x(i,j,k).lt.maxvalue) then  
                  imax(1) = i
                  imax(2) = j
                  imax(3) = k
                  maxvalue = x(i,j,k)
               endif
            end do
         end do
      end do
    end function mpmaxlocation3

    function mpminlocation1(x) result(imin)
      type(mp_real), intent(in):: x(:)
      type(mp_real) :: minvalue
      integer :: i 
      integer :: imin(1)

      imin(1)= lbound(x,dim=1)
      minvalue = x(lbound(x,dim=1))
      do i=lbound(x,dim=1), ubound(x,dim=1)
         if (x(i).lt.minvalue) then  
            imin(1) = i
            minvalue = x(i)
         endif
      end do

    end function mpminlocation1

    function mpminlocation2(x) result(imin)
      type(mp_real), intent(in):: x(:,:)
      type(mp_real) :: minvalue
      integer :: i,j 
      integer :: imin(2)

      imin(1)= lbound(x,dim=1)
      imin(2)= lbound(x,dim=2)
      minvalue = x(lbound(x,dim=1),lbound(x,dim=2))
         do j=lbound(x,dim=2), ubound(x,dim=2)
      do i=lbound(x,dim=1), ubound(x,dim=1)
            if (x(i,j).lt.minvalue) then  
               imin(1) = i
               imin(2) = j
               minvalue = x(i,j)
            endif
         end do
      end do
    end function mpminlocation2

  
    function mpminlocation3(x) result(imin)
      type(mp_real), intent(in):: x(:,:,:)
      type(mp_real) :: minvalue
      integer :: i,j,k 
      integer :: imin(3)

      imin(1)= lbound(x,dim=1)
      imin(2)= lbound(x,dim=2)
      imin(3)= lbound(x,dim=3)
      minvalue = x(lbound(x,dim=1),lbound(x,dim=2),lbound(x,dim=3))
            do k=lbound(x,dim=3), ubound(x,dim=3)
         do j=lbound(x,dim=2), ubound(x,dim=2)
      do i=lbound(x,dim=1), ubound(x,dim=1)
               if (x(i,j,k).lt.minvalue) then  
                  imin(1) = i
                  imin(2) = j
                  imin(3) = k
                  minvalue = x(i,j,k)
               endif
            end do
         end do
      end do
    end function mpminlocation3     

    function mpmaxval1(x) result(maxvalue)
      type(mp_real), intent(in):: x(:)
      type(mp_real) :: maxvalue
      integer :: i 

      maxvalue = x(lbound(x,dim=1))
      do i=lbound(x,dim=1),ubound(x,dim=1)
         if (x(i).gt.maxvalue) then  
            maxvalue = x(i)
         endif
      end do

    end function mpmaxval1

    function mpmaxval2(x) result(maxvalue)
      type(mp_real), intent(in):: x(:,:)
      type(mp_real) :: maxvalue
      integer :: i,j 

      maxvalue = x(lbound(x,dim=1),lbound(x,dim=2))
         do j=lbound(x,dim=2),ubound(x,dim=2)
      do i=lbound(x,dim=1),ubound(x,dim=1)
            if (x(i,j).gt.maxvalue) then  
               maxvalue = x(i,j)
            endif
         end do
      end do

    end function mpmaxval2

    function mpmaxval3(x) result(maxvalue)
      type(mp_real), intent(in):: x(:,:,:)
      type(mp_real) :: maxvalue
      integer :: i,j,k 

      maxvalue = x(lbound(x,dim=1),lbound(x,dim=2),lbound(x,dim=3))
            do k=lbound(x,dim=3),ubound(x,dim=3)
         do j=lbound(x,dim=2),ubound(x,dim=2)
      do i=lbound(x,dim=1),ubound(x,dim=1)
               if (x(i,j,k).gt.maxvalue) then  
                  maxvalue = x(i,j,k)
               endif
            end do
         end do
      end do
    end function mpmaxval3

    function mpminval1(x) result(minvalue)
      type(mp_real), intent(in):: x(:)
      type(mp_real) :: minvalue
      integer :: i 

      minvalue = x(lbound(x,dim=1))
      do i=lbound(x,dim=1),ubound(x,dim=1)
         if (x(i).lt.minvalue) then  
            minvalue = x(i)
         endif
      end do

    end function mpminval1


    function mpminval2(x) result(minvalue)
      type(mp_real), intent(in):: x(:,:)
      type(mp_real) :: minvalue
      integer :: i,j 

      minvalue = x(lbound(x,dim=1),lbound(x,dim=2))
         do j=lbound(x,dim=2),ubound(x,dim=2)
      do i=lbound(x,dim=1),ubound(x,dim=1)
            if (x(i,j).gt.minvalue) then  
               minvalue = x(i,j)
            endif
         end do
      end do

    end function mpminval2

    function mpminval3(x) result(minvalue)
      type(mp_real), intent(in):: x(:,:,:)
      type(mp_real) :: minvalue
      integer :: i,j,k 

      minvalue = x(lbound(x,dim=1),lbound(x,dim=2),lbound(x,dim=3))
            do k=lbound(x,dim=3),ubound(x,dim=3)
         do j=lbound(x,dim=2),ubound(x,dim=2)
      do i=lbound(x,dim=1),ubound(x,dim=1)
               if (x(i,j,k).gt.minvalue) then  
                  minvalue = x(i,j,k)
               endif
            end do
         end do
      end do
    end function mpminval3

    function mpdot_product(a,b) result(c)
      type (mp_real),intent(in) ::  a(:),b(:)
      type (mp_real) ::  c
      integer :: i 

      !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdot_product' 

      c=mpreal('0.0')
      do i=lbound(a,dim=1),ubound(a,dim=1)
         c=c+a(i)*b(i) 
      enddo
    end function mpdot_product

    function mpsum(a) result(c)
      type (mp_real),intent(in) ::  a(:)
      type (mp_real) ::  c
      integer :: i
      
      c = mpreal('0.0')
      do i=lbound(a,dim=1),ubound(a,dim=1)
         c=c+a(i)
      enddo
    end function mpsum

    function mpsum_2ddim(a,dim) result(c)
      type (mp_real),intent(in) ::  a(:,:)
      integer,       intent(in) :: dim
      type (mp_real) ::  c(size(a,dim=3-dim))
      integer :: i,j
      
      if (dim == 1) then
         do j = lbound(a,dim=2),ubound(a,dim=2)
            c(j) = mpreal('0.0')
         end do
         do j = lbound(a,dim=2),ubound(a,dim=2)
            do i = lbound(a,dim=1),ubound(a,dim=1)
               c(j) = c(j) + a(i,j)
            end do
         end do
      else if (dim == 2) then
         do i = lbound(a,dim=1),ubound(a,dim=1)
            c(i) = mpreal('0.0')
         end do

            do j = lbound(a,dim=2),ubound(a,dim=2)
         do i = lbound(a,dim=1),ubound(a,dim=1)
               c(i) = c(i) + a(i,j)
            end do
         end do
      else
         stop 'mpsum_2ddim: dim was not 1 or 2'
      end if
    end function mpsum_2ddim

    subroutine mpwritearray1(a)
      type (mp_real) a(:)
      integer :: i 

      do i=lbound(a,dim=1),ubound(a,dim=1) 
         call mpwrite(6,a(i))  
      enddo
    end subroutine mpwritearray1


    subroutine mpwritearray2(a)
      type (mp_real) a(:,:)
      integer :: i,j 

         do j=lbound(a,dim=2),ubound(a,dim=2) 
      do i=lbound(a,dim=1),ubound(a,dim=1) 
            call mpwrite(6,a(i,j))  
         enddo
         write(*,*)
      enddo
    end subroutine mpwritearray2

    subroutine mpwritearray3(a)
      type (mp_real) a(:,:,:)
      integer :: i,j ,k

            do k=lbound(a,dim=3),ubound(a,dim=3) 
         do j=lbound(a,dim=2),ubound(a,dim=2) 
      do i=lbound(a,dim=1),ubound(a,dim=1) 
               call mpwrite(6,a(i,j,k))  
            enddo
            write(*,*)
         enddo
         write(*,*); write(*,*) 
      enddo
    end subroutine mpwritearray3

    function mpfloor(a) result(res)
      type(mp_real), intent(in) :: a
      integer :: res
      
      res = int(a)
      if (res > a) res = res-1

    end function mpfloor

    function mpceiling(a) result(res)
      type(mp_real), intent(in) :: a
      integer :: res

      res = int(a)
      if (res < a) res = res+1

    end function mpceiling

    


end module mpadvancedoperations

