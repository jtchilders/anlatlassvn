module mpoperations_hidden_c 
  use mpmodule
  use mpconverter 
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
  public :: mpscalprodarc, mpscalprodbrc, mpscalprodacr, mpscalprodbcr
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
  public :: mptimescr1,mptimescr2,mptimescr3
  public :: mptimesrc1,mptimesrc2,mptimesrc3
  public :: mpdivision1,mpdivision2,mpdivision3
  public :: mpintdivision1,mpintdivision2,mpintdivision3
  public :: mpexpo
  public :: mpsquareroot1,mpsquareroot2,mpsquareroot3
  public :: mplog1,mplog2,mplog3
  public :: mpexp1,mpexp2,mpexp3
  public :: mpabs1,mpabs2,mpabs3

  public :: mpequal10,mpequal11,mpequal2,mpequal3

contains 

  ! set a constant mp real value for an mp matrix
  subroutine mpfixmpvalue1(a,b) 
    type (mp_complex),intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       a(i)=b
    enddo
  end subroutine mpfixmpvalue1

  subroutine mpfixmpvalue2(a,b) 
    type (mp_complex),intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       end do
    enddo
  end subroutine mpfixmpvalue2

  subroutine mpfixmpvalue3(a,b) 
    type (mp_complex),intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:,:)
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
    type (mp_complex),intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:,:,:)
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
    type (mp_complex),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       a(i)=b
    enddo
  end subroutine mpfixrealvalue1

  subroutine mpfixrealvalue2(a,b) 
    real,intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       end do
    end do
  end subroutine mpfixrealvalue2

  subroutine mpfixrealvalue3(a,b) 
    real,intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:,:)
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
    type (mp_complex),intent(inout) :: a(:,:,:,:)
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
    type (mp_complex),intent(inout) :: a(:)
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)      
       a(i)=b
    enddo
  end subroutine mpfixintegervalue1

  subroutine mpfixintegervalue2(a,b) 
    integer,intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:)
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          a(i,j)=b
       enddo
    end do
  end subroutine mpfixintegervalue2

  subroutine mpfixintegervalue3(a,b) 
    integer,intent(in) ::  b
    type (mp_complex),intent(inout) :: a(:,:,:)
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
    type (mp_complex),intent(inout) :: a(:,:,:,:)
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
    type (mp_complex),intent(inout) :: a(:)
    integer :: i 

    !if (size(a).ne.size(b)) stop 'Size 1 differ in  mpfixintegerarray'
    do i=lbound(a,dim=1),ubound(a,dim=1)      
       a(i)=b(i)
    enddo
  end subroutine mpfixintegerarray1

  subroutine mpfixintegerarray2(a,b) 
    integer,intent(in) ::  b(:,:)
    type (mp_complex),intent(inout) :: a(:,:)
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
    type (mp_complex),intent(inout) :: a(:,:,:)
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
    type (mp_complex),intent(inout) :: a(:,:,:,:)
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
    type (mp_complex),intent(in) ::  b(:)
    type (mp_complex), intent(inout) :: a(size(b,dim=1))
    integer :: i 

    do i=lbound(b,dim=1),ubound(b,dim=1) 
       a(i)=b(i)
    enddo
  end subroutine mpsetequal1

  subroutine mpsetequal2(a, b) 
    type (mp_complex),intent(in) ::  b(:,:)
    type (mp_complex), intent(inout) :: a(size(b,dim=1),size(b,dim=2))
    integer :: i,j 

    do j=lbound(b,dim=2),ubound(b,dim=2) 
       do i=lbound(b,dim=1),ubound(b,dim=1) 
          a(i,j)=b(i,j)
       end do
    end do
  end subroutine mpsetequal2

  subroutine mpsetequal3(a, b) 
    type (mp_complex),intent(in) ::  b(:,:,:)
    type (mp_complex), intent(inout) :: a(size(b,dim=1),size(b,dim=2),size(b,dim=3))

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
    type (mp_complex),intent(in) ::  b(:,:,:,:)
    type (mp_complex), intent(inout) :: a(size(b,dim=1),size(b,dim=2),size(b,dim=3),size(b,dim=4))
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
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex)            ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    end do
  end function mpscalproda


  function mpscalprodb(b,a) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    enddo
  end function mpscalprodb

  function mpscalprodarc(a,b) result(c)
    type (mp_real),intent(in) :: a(:)
    type (mp_complex),intent(in) :: b
    type (mp_complex)            ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    end do
  end function mpscalprodarc


  function mpscalprodbrc(b,a) result(c)
    type (mp_real),intent(in) :: a(:)
    type (mp_complex),intent(in) :: b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    enddo
  end function mpscalprodbrc

  function mpscalprodacr(a,b) result(c)
    type (mp_complex),intent(in) :: a(:)
    type (mp_real),intent(in) :: b
    type (mp_complex)            ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    end do
  end function mpscalprodacr


  function mpscalprodbcr(b,a) result(c)
    type (mp_complex),intent(in) :: a(:)
    type (mp_real),intent(in) :: b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)*b
    enddo
  end function mpscalprodbcr




  function mpscalproda2(a,b) result(c)
    type (mp_complex),intent(in) :: a(:,:), b
    type (mp_complex) ::  c(lbound(a,dim=1):ubound(a,dim=1),&
         &lbound(a,dim=2):ubound(a,dim=2))
    integer :: i,j

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)*b
       end do
    end do
  end function mpscalproda2

  function mpscalprodb2(b,a) result(c)
    type (mp_complex),intent(in) :: a(:,:), b
    type (mp_complex) ::  c(lbound(a,dim=1):ubound(a,dim=1),&
         &lbound(a,dim=2):ubound(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)*b
       enddo
    end do
  end function mpscalprodb2

  function mpscaldiv1(a,b) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(ubound(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)/b
    enddo
  end function mpscaldiv1

  function mpscaldiv2(a,b) result(c)
    type (mp_complex),intent(in) :: a(:,:), b
    type (mp_complex) ::  c(ubound(a,dim=1),ubound(a,dim=2))
    integer :: i, j

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       do j=lbound(a,dim=2) ,ubound(a,dim=2) 
          c(i,j)=a(i,j)/b
       end do
    enddo
  end function mpscaldiv2

  ! scalar sum
  function mpscalsuma(a,b) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)+b
    enddo
  end function mpscalsuma


  function mpscalsumb(b,a) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=b+a(i)
    enddo
  end function mpscalsumb

  function mpscalsuma2(a,b) result(c)
    type (mp_complex),intent(in) :: a(:,:), b
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       enddo
    end do
  end function mpscalsuma2


  function mpscalsumb2(b,a) result(c)
    type (mp_complex),intent(in) :: a(:,:), b
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       end do
    enddo
  end function mpscalsumb2

  function mpintsuma2(a,b) result(c)
    type (mp_complex),intent(in) :: a(:,:)
    integer, intent(in) :: b
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       enddo
    end do
  end function mpintsuma2


  function mpintsumb2(b,a) result(c)
    type (mp_complex),intent(in) :: a(:,:)
    integer, intent(in) :: b
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 
    do j=lbound(a,dim=2) ,ubound(a,dim=2) 
       do i=lbound(a,dim=1) ,ubound(a,dim=1) 
          c(i,j)=a(i,j)+b
       end do
    enddo
  end function mpintsumb2


  ! scalar difference
  function mpscaldifa(a,b) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=a(i)-b
    enddo
  end function mpscaldifa


  function mpscaldifb(b,a) result(c)
    type (mp_complex),intent(in) :: a(:), b
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1) ,ubound(a,dim=1) 
       c(i)=b-a(i)
    enddo
  end function mpscaldifb


  ! some binary operators with pure vectors now
  function mpplus1(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:),b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1) /= size(b,dim=1) ) stop "Sizes differ in mpplus"

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)+b(i) 
    enddo
  end function mpplus1

  function mpplus2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:),b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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


  function mpminus1(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:),b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminus1

  function mpminusu1(a) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=-a(i)
    enddo
  end function mpminusu1

  function mpminus2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:),b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  b(:)
    integer,intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminusintegerarray1a

  function mpminusintegerarray1b(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:)
    integer,intent(in) ::  b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpminus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpminusintegerarray1b

  function mpminusintegerarray2a(a,b) result(c)
    type (mp_complex),intent(in) ::  b(:,:)
    integer,intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:)
    integer,intent(in) ::  b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  b(:)
    integer,intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpplusintegerarray1a

  function mpplusintegerarray1b(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:)
    integer,intent(in) ::  b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpplus' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)-b(i) 
    enddo
  end function mpplusintegerarray1b

  function mpplusintegerarray2a(a,b) result(c)
    type (mp_complex),intent(in) ::  b(:,:)
    integer,intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:)
    integer,intent(in) ::  b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:),b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes1' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)*b(i) 
    enddo
  end function mptimes1

  function mptimes2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:),b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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


  function mptimescr1(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_real),intent(in) ::  b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes1' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)*b(i) 
    enddo
  end function mptimescr1

  function mptimescr2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_real),intent(in) ::  b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes2' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mptimes2'

    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)*b(i,j) 
       end do
    end do
  end function mptimescr2

  function mptimescr3(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_real),intent(in) ::  b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
  end function mptimescr3


  function mptimesrc1(b,a) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_real),intent(in) ::  b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes1' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)*b(i) 
    enddo
  end function mptimesrc1

  function mptimesrc2(b,a) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_real),intent(in) ::  b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mptimes2' 
    !if (size(a,dim=2).ne.size(b,dim=2)) stop 'Size 2 differ in mptimes2'

    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)*b(i,j) 
       end do
    end do
  end function mptimesrc2

  function mptimesrc3(b,a) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_real),intent(in) ::  b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
  end function mptimesrc3


  function mpdivision1(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:),b(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdivisions' 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)/b(i) 
    enddo
  end function mpdivision1

  function mpdivision2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:),b(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
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
    type (mp_complex),intent(in) ::  a(:,:,:),b(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:)
    integer, intent(in) :: b 
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)/b
    enddo
  end function mpintdivision1

  function mpintdivision2(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    integer, intent(in) :: b 
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 


    do j=lbound(a,dim=1),ubound(a,dim=1)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=a(i,j)/b
       end do
    end do
  end function mpintdivision2

  function mpintdivision3(a,b) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    integer, intent(in) :: b 
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer, intent(in) :: n 
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=a(i)**n
    enddo

  end function mpexpo


  function mpsquareroot1(a) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=sqrt(a(i))
    enddo

  end function mpsquareroot1

  function mpsquareroot2(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=sqrt(a(i,j))
       enddo
    end do
  end function mpsquareroot2

  function mpsquareroot3(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=sqrt(a(i,j,k))
          enddo
       end do
    end do

  end function mpsquareroot3

  function mplog1(a) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=log(a(i))
    enddo

  end function mplog1

  function mplog2(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=log(a(i,j))
       enddo
    end do
  end function mplog2

  function mplog3(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=log(a(i,j,k))
          enddo
       end do
    end do

  end function mplog3

  function mpexp1(a) result(c)
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i,j 

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=exp(a(i))
    enddo

  end function mpexp1

  function mpexp2(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2) 
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=exp(a(i,j))
       enddo
    end do
  end function mpexp2

  function mpexp3(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
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
    type (mp_complex),intent(in) ::  a(:)
    type (mp_complex) ::  c(size(a,dim=1))
    integer :: i

    do i=lbound(a,dim=1),ubound(a,dim=1)
       c(i)=abs(a(i))
    enddo

  end function mpabs1

  function mpabs2(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2))
    integer :: i,j 

    do j=lbound(a,dim=2),ubound(a,dim=2)
       do i=lbound(a,dim=1),ubound(a,dim=1)
          c(i,j)=abs(a(i,j))
       enddo
    end do
  end function mpabs2

  function mpabs3(a) result(c)
    type (mp_complex),intent(in) ::  a(:,:,:)
    type (mp_complex) ::  c(size(a,dim=1),size(a,dim=2),size(a,dim=3))
    integer :: i,j,k 

    do k=lbound(a,dim=3),ubound(a,dim=3)
       do j=lbound(a,dim=2),ubound(a,dim=2)
          do i=lbound(a,dim=1),ubound(a,dim=1)
             c(i,j,k)=abs(a(i,j,k))
          enddo
       end do
    end do

  end function mpabs3

  function mpequal11 (a, b) result (t_f) 
    type (mp_complex), intent(in) :: a(:)       
    type (mp_complex), intent(in) :: b(:)
    logical :: t_f(size(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i).ne.b(i)) t_f(i)= .false.
    end do

  end function mpequal11
  function mpequal10 (a, b) result (t_f) 
    type (mp_complex), intent(in) :: a(:)       
    type (mp_complex), intent(in) :: b      
    logical :: t_f(size(a,dim=1))
    integer :: i 
    t_f = .true.                      

    do i = lbound(a,dim=1),ubound(a,dim=1)
       if (a(i).ne.b) t_f(i)= .false.
    end do

  end function mpequal10

  function mpequal2 (a, b) result (t_f) 
    type (mp_complex), intent(in) :: a(:,:)       
    type (mp_complex), intent(in) :: b      
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
    type (mp_complex), intent(in) :: a(:,:,:)       
    type (mp_complex), intent(in) :: b      
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

  !
end module mpoperations_hidden_c

module mpsimpleoperations_c
  use mpoperations_hidden_c  
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
     module procedure mptimescr1, mptimescr2,mptimescr3
     module procedure mptimesrc1, mptimesrc2,mptimesrc3
     module procedure mpscalproda, mpscalprodb, mpscalproda2, mpscalprodb2
     module procedure mpscalprodarc, mpscalprodbrc, mpscalprodacr, mpscalprodbcr
  end interface            
  public :: operator(*)
  private mptimes1, mptimes2,mptimes3, mpscalproda, mpscalprodb,&
       &mpscalproda2, mpscalprodb2,mpscalprodarc, mpscalprodbrc, &
       &mpscalprodacr, mpscalprodbcr
  private mptimescr1, mptimescr2,mptimescr3
  private mptimesrc1, mptimesrc2,mptimesrc3

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

end module mpsimpleoperations_c

module mpadvancedoperations_c  
  use mpmodule!; use consts_mp 
  implicit none 
  private 

! stupid routine to print (:), (:,:) and (:,:,:) mpmatrices  
  interface mpwritearray
     module procedure mpwritearray1,mpwritearray2,mpwritearray3
  end interface            
  public :: mpwritearray

  interface dot_product
     module procedure mpdot_product
  end interface
  public :: dot_product

 interface sum
     module procedure mpsum, mpsum_2ddim
  end interface            
  public :: sum 

  contains 


    function mpdot_product(a,b) result(c)
      type (mp_complex),intent(in) ::  a(:),b(:)
      type (mp_complex) ::  c
      integer :: i 

      !if (size(a,dim=1).ne.size(b,dim=1)) stop 'Size 1 differ in mpdot_product' 

      c=mpcmpl('0.0','0.0')
      do i=lbound(a,dim=1),ubound(a,dim=1)
         c=c+a(i)*b(i) 
      enddo
    end function mpdot_product

    function mpsum(a) result(c)
      type (mp_complex),intent(in) ::  a(:)
      type (mp_complex) ::  c
      integer :: i
      
      c = mpcmpl('0.0','0.0')
      do i=lbound(a,dim=1),ubound(a,dim=1)
         c=c+a(i)
      enddo
    end function mpsum

    function mpsum_2ddim(a,dim) result(c)
      type (mp_complex),intent(in) ::  a(:,:)
      integer,       intent(in) :: dim
      type (mp_complex) ::  c(size(a,dim=3-dim))
      integer :: i,j
      
      if (dim == 1) then
         do j = lbound(a,dim=2),ubound(a,dim=2)
            c(j) = mpcmpl('0.0','0.0') 
         end do
            do i = lbound(a,dim=1),ubound(a,dim=1)
         do j = lbound(a,dim=2),ubound(a,dim=2)
               c(j) = c(j) + a(i,j)
            end do
         end do
      else if (dim == 2) then
         do i = lbound(a,dim=1),ubound(a,dim=1)
            c(i) = mpcmpl('0.0','0.0') 
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
      type (mp_complex) a(:)
      integer :: i 

      do i=lbound(a,dim=1),ubound(a,dim=1) 
         call mpwrite(6,a(i))  
      enddo
    end subroutine mpwritearray1


    subroutine mpwritearray2(a)
      type (mp_complex) a(:,:)
      integer :: i,j 

         do j=lbound(a,dim=2),ubound(a,dim=2) 
      do i=lbound(a,dim=1),ubound(a,dim=1) 
            call mpwrite(6,a(i,j))  
         enddo
         write(*,*)
      enddo
    end subroutine mpwritearray2

    subroutine mpwritearray3(a)
      type (mp_complex) a(:,:,:)
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


end module mpadvancedoperations_c

