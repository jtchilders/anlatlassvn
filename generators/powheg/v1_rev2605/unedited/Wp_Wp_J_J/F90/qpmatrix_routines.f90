!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECmatrix_routines.f90.

module qpmatrix_routines
  use types; use consts_qp; use warnings_and_errors
  implicit none

  private

  interface mydet 
     module procedure mydetc
  end interface

  interface complement 
     module procedure complementc
  end interface


  public :: determinant_old, mydet

contains


  recursive function mydetc(A) result (det)
    use types; use consts_qp
    complex(qp),intent(in) :: A(:,:)
    complex(qp) :: det
    ! dummy variables
    integer     :: m,n,i
    complex(qp) ::  r11,r12,r13,r14,r15,r21,r22,r23,r24,r25
    complex(qp) ::  r31,r32,r33,r34,r35,r41,r42,r43,r44,r45
    complex(qp) ::  r51,r52,r53,r54,r55

    n = size(A,dim=1) 
    m = size(A,dim=2) 
    if (n/=m) call wae_error('mydet','Not squared matrix')
    det = zero 
    if (n==1) then 
       det = A(1,1)
    elseif (n == 2) then 
       det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    elseif (n == 3) then 
       det  =A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+&
            &A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)-&
            &A(1,2)*A(2,1)*A(3,3)-A(1,1)*A(2,3)*A(3,2)
    elseif (n == 4) then 
       det = &
            &      A(1,4)*A(2,3)*A(3,2)*A(4,1) - A(1,3)*A(2,4)*A(3,2)*A(4,1) - &
            &      A(1,4)*A(2,2)*A(3,3)*A(4,1) + A(1,2)*A(2,4)*A(3,3)*A(4,1) + &
            &      A(1,3)*A(2,2)*A(3,4)*A(4,1) - A(1,2)*A(2,3)*A(3,4)*A(4,1) - & 
            &      A(1,4)*A(2,3)*A(3,1)*A(4,2) + A(1,3)*A(2,4)*A(3,1)*A(4,2) + &
            &      A(1,4)*A(2,1)*A(3,3)*A(4,2) - A(1,1)*A(2,4)*A(3,3)*A(4,2) - &
            &      A(1,3)*A(2,1)*A(3,4)*A(4,2) + A(1,1)*A(2,3)*A(3,4)*A(4,2) + &
            &      A(1,4)*A(2,2)*A(3,1)*A(4,3) - A(1,2)*A(2,4)*A(3,1)*A(4,3) - &
            &      A(1,4)*A(2,1)*A(3,2)*A(4,3) + A(1,1)*A(2,4)*A(3,2)*A(4,3) + &
            &      A(1,2)*A(2,1)*A(3,4)*A(4,3) - A(1,1)*A(2,2)*A(3,4)*A(4,3) - &
            &      A(1,3)*A(2,2)*A(3,1)*A(4,4) + A(1,2)*A(2,3)*A(3,1)*A(4,4) + &
            &      A(1,3)*A(2,1)*A(3,2)*A(4,4) - A(1,1)*A(2,3)*A(3,2)*A(4,4) - &
            &      A(1,2)*A(2,1)*A(3,3)*A(4,4) + A(1,1)*A(2,2)*A(3,3)*A(4,4)      

    elseif (n == 5) then 
       r11=A(1,1);r12=A(1,2);r13=A(1,3);r14=A(1,4);r15=A(1,5)
       r21=A(2,1);r22=A(2,2);r23=A(2,3);r24=A(2,4);r25=A(2,5)
       r31=A(3,1);r32=A(3,2);r33=A(3,3);r34=A(3,4);r35=A(3,5)
       r41=A(4,1);r42=A(4,2);r43=A(4,3);r44=A(4,4);r45=A(4,5)
       r51=A(5,1);r52=A(5,2);r53=A(5,3);r54=A(5,4);r55=A(5,5)

       det = &
            &      r15*r24*r33*r42*r51 - r14*r25*r33*r42*r51 - &
            &      r15*r23*r34*r42*r51 + r13*r25*r34*r42*r51 + &
            &      r14*r23*r35*r42*r51 - r13*r24*r35*r42*r51 - &
            &      r15*r24*r32*r43*r51 + r14*r25*r32*r43*r51 + &
            &      r15*r22*r34*r43*r51 - r12*r25*r34*r43*r51 - &
            &      r14*r22*r35*r43*r51 + r12*r24*r35*r43*r51 + &
            &      r15*r23*r32*r44*r51 - r13*r25*r32*r44*r51 - &
            &      r15*r22*r33*r44*r51 + r12*r25*r33*r44*r51 + &
            &      r13*r22*r35*r44*r51 - r12*r23*r35*r44*r51 - &
            &      r14*r23*r32*r45*r51 + r13*r24*r32*r45*r51 + &
            &      r14*r22*r33*r45*r51 - r12*r24*r33*r45*r51 - &
            &      r13*r22*r34*r45*r51 + r12*r23*r34*r45*r51 - &
            &      r15*r24*r33*r41*r52 + r14*r25*r33*r41*r52 + &
            &      r15*r23*r34*r41*r52 - r13*r25*r34*r41*r52 - &
            &      r14*r23*r35*r41*r52 + r13*r24*r35*r41*r52 + &
            &      r15*r24*r31*r43*r52 - r14*r25*r31*r43*r52 - &
            &      r15*r21*r34*r43*r52 + r11*r25*r34*r43*r52 + &
            &      r14*r21*r35*r43*r52 - r11*r24*r35*r43*r52 - &
            &      r15*r23*r31*r44*r52 + r13*r25*r31*r44*r52 + &
            &      r15*r21*r33*r44*r52 - r11*r25*r33*r44*r52 - &
            &      r13*r21*r35*r44*r52 + r11*r23*r35*r44*r52 + &
            &      r14*r23*r31*r45*r52 - r13*r24*r31*r45*r52 - &
            &      r14*r21*r33*r45*r52 + r11*r24*r33*r45*r52 + &
            &      r13*r21*r34*r45*r52 - r11*r23*r34*r45*r52 + &
            &      r15*r24*r32*r41*r53 - r14*r25*r32*r41*r53 - &
            &      r15*r22*r34*r41*r53 + r12*r25*r34*r41*r53 + &
            &      r14*r22*r35*r41*r53 - r12*r24*r35*r41*r53 - &
            &      r15*r24*r31*r42*r53 + r14*r25*r31*r42*r53 + &
            &      r15*r21*r34*r42*r53 - r11*r25*r34*r42*r53 - &
            &      r14*r21*r35*r42*r53 + r11*r24*r35*r42*r53 + &
            &      r15*r22*r31*r44*r53 - r12*r25*r31*r44*r53 - &
            &      r15*r21*r32*r44*r53 + r11*r25*r32*r44*r53 + &
            &      r12*r21*r35*r44*r53 - r11*r22*r35*r44*r53 - &
            &      r14*r22*r31*r45*r53 + r12*r24*r31*r45*r53 + &
            &      r14*r21*r32*r45*r53 - r11*r24*r32*r45*r53 - &
            &      r12*r21*r34*r45*r53 + r11*r22*r34*r45*r53 - &
            &      r15*r23*r32*r41*r54 + r13*r25*r32*r41*r54 + &
            &      r15*r22*r33*r41*r54 - r12*r25*r33*r41*r54 - &
            &      r13*r22*r35*r41*r54 + r12*r23*r35*r41*r54 + &
            &      r15*r23*r31*r42*r54 - r13*r25*r31*r42*r54 - &
            &      r15*r21*r33*r42*r54 + r11*r25*r33*r42*r54 + &
            &      r13*r21*r35*r42*r54 - r11*r23*r35*r42*r54 - &
            &      r15*r22*r31*r43*r54 + r12*r25*r31*r43*r54 + &
            &      r15*r21*r32*r43*r54 - r11*r25*r32*r43*r54 - &
            &      r12*r21*r35*r43*r54 + r11*r22*r35*r43*r54 + &
            &      r13*r22*r31*r45*r54 - r12*r23*r31*r45*r54 - &
            &      r13*r21*r32*r45*r54 + r11*r23*r32*r45*r54 + &
            &      r12*r21*r33*r45*r54 - r11*r22*r33*r45*r54 + &
            &      r14*r23*r32*r41*r55 - r13*r24*r32*r41*r55 - &
            &      r14*r22*r33*r41*r55 + r12*r24*r33*r41*r55 + &
            &      r13*r22*r34*r41*r55 - r12*r23*r34*r41*r55 - &
            &      r14*r23*r31*r42*r55 + r13*r24*r31*r42*r55 + &
            &      r14*r21*r33*r42*r55 - r11*r24*r33*r42*r55 - &
            &      r13*r21*r34*r42*r55 + r11*r23*r34*r42*r55 + &
            &      r14*r22*r31*r43*r55 - r12*r24*r31*r43*r55 - &
            &      r14*r21*r32*r43*r55 + r11*r24*r32*r43*r55 + &
            &      r12*r21*r34*r43*r55 - r11*r22*r34*r43*r55 - &
            &      r13*r22*r31*r44*r55 + r12*r23*r31*r44*r55 + &
            &      r13*r21*r32*r44*r55 - r11*r23*r32*r44*r55 - &
            &      r12*r21*r33*r44*r55 + r11*r22*r33*r44*r55

    else
       do i=1,n
          if (A(1,i) /= zero) then 
             det = det + (-1)**(i+1)*A(1,i)*mydetc(Complement(A,1,i))
          end if
       end do
    end if

  end function mydetc



  function determinant_old(matin) result(det)
    use types; use consts_qp
    complex(qp), intent(in) :: matin(:,:)
    integer     :: n,i,j
    complex(qp) :: det
    complex(qp) :: temp(size(matin,dim=2))
    complex(qp) :: mat(size(matin,dim=1),size(matin,dim=2))
    integer :: k,kmax
    real(qp) :: maxel 

    mat=matin
    n=size(matin,dim=2)

    outer: do j=1,n-1
       inner: do i=n,j+1,-1
          if (mat(j,i-1) == zero) then
             temp(:)   = mat(:,i-1)
             mat(:,i-1)= mat(:,i)
             mat(:,i)  = -temp(:) ! need - sign when interchanging adjacent rows 
          else

             if (abs(mat(j,i-1)) < sq2tol) then
                ! find largest value to the right...
                maxel = abs(mat(j,i-1)) 
                kmax = i-1
                do k=i-2,j,-1
                   if (abs(mat(j,k)) > maxel) then 
                      kmax = k
                   end if
                end do
                if (kmax == i-1 .or. abs(mat(j,kmax)) < sq2tol) then 
                   det = mydet(matin)
                   return 
                end if
                temp(:) = mat(:,i-1)
                mat(:,i-1) = mat(:,kmax)
                mat(:,kmax) = -temp(:)
             end if
             mat(:,i)=mat(:,i)-mat(j,i)*mat(:,i-1)/mat(j,i-1)
          end if
       enddo inner
    enddo outer


    det=one
    do i=1,n
       det=mat(i,i)*det 
    end do

  end function determinant_old

  function complementc(mat,i,j) result(Aij)
    use types; use consts_qp
    complex(qp), intent(in) :: mat(:,:)
    integer, intent(in)     :: i,j 
    complex(qp) :: Aij(size(mat,dim=1)-1,size(mat,dim=2)-1)

    Aij=czero
    Aij(:i-1,:j-1)= mat(:i-1,:j-1)
    Aij(:i-1,j:)  = mat(:i-1,j+1:)
    Aij(i:,:j-1)  = mat(i+1:,:j-1)
    Aij(i:,j:)    = mat(i+1:,j+1:)

  end function complementc



end module qpmatrix_routines

