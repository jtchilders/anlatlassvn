!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECvvnbasis.f90.

! --- module for vermasseren-van neerven basis
module qpvvn
  use types; use consts_qp
  use qpgram; use qpaux_functions
  implicit none 
  private 

  public :: give4to4vect,give3to4vect,give2to4vect,give1to4vect
  public :: give1to4vect_light

contains


  !-------------------------------------------------------------------------
  ! The procedure that returns 4 vectors from 4 vectors; 
  ! relevant for ``pentagon'' basis
  !------------------------------------------------------------------------
  subroutine give4to4vect(k1,k2,k3,k4,v)
    complex(qp), intent(in) :: k1(:),k2(:),k3(:),k4(:)
    complex(qp), intent(out) :: v(:,:)
    ! -------------------------------------------------      
    complex(qp) :: v1(4),v2(4),v3(4),v4(4)
    complex(qp) :: d4

    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give4to4vect: wrong size of v'
    !if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /=4.or. size(k4) /=4 ) &
    !     &stop 'give1to4vect: wrong size of k1, k2, k3 or. k4'


    call give4vect(k2,k3,k4,k1,v1)
    call give4vect(k1,k3,k4,k2,v2)
    call give4vect(k1,k2,k4,k3,v3)
    call give4vect(k1,k2,k3,k4,v4)
    d4 = sc(k1,v1)

    if (abs(d4) > tol) then 
       v(1,:) = v1/d4
       v(2,:) = v2/d4
       v(3,:) = v3/d4
       v(4,:) = v4/d4
    elseif (abs(d4) < tol) then 
       write(*,*) 'give4to4vect: d4< tol ?', (d4), (tol)
    else
       write(*,*) 'give4to4vect: d4=NaN ?', (d4)
    endif


  end subroutine give4to4vect



  !-------------------------------------------------------------------------
  ! The procedure that returns 4 vectors from 3 vectors; `box'' basis
  !-------------------------------------------------------------------------
  subroutine give3to4vect(k1,k2,k3,v)
    complex(qp), intent(in) :: k1(:),k2(:),k3(:)
    complex(qp), intent(out) :: v(:,:) 
    ! -------------------------------------------------      
    complex(qp) :: v1(4),v2(4),v3(4),v4(4)
    complex(qp) :: vaux(4),vauxa(4),v3a(4)
    complex(qp) :: d3,d2,ko1,ko2,ko3
    complex(qp) :: a11,a12,a22,a13,a23,dnorm,a2v,a1v
    complex(qp) :: a3v

    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give3to4vect: wrong size of v'
    !if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /=4) &
    !     &stop 'give3to4vect: wrong size of k1, k2 or. k3'


    call give3vect(k2,k3,k1,v1)
    call give3vect(k1,k3,k2,v2)
    call give3vect(k1,k2,k3,v3)
    d3 = sc(k1,v1)

    if (abs(d3) > tol) then 
       v1 = v1/d3
       v2 = v2/d3
       v3 = v3/d3
    elseif (abs(d3) < tol) then 
       write(*,*) 'give3to4vect: d3< tol ?', (d3)
    else
       write(*,*) 'give3to4vect: d3=NaN ?', (d3)
    endif

    ! I need to construct a vector in a plane transverse to k1,k1,k3
    ! To do that, I first create the metric tensor of the plane, 
    ! transverse to (k1,k2):
    ! this is given by three coefficients:

    a12 = sc(k1,k2)
    a11 = sc(k1,k1)
    a22 = sc(k2,k2)

    a13 = sc(k1,k3)
    a23 = sc(k2,k3)

    d2=a12**2 -a11*a22

    if (abs(d2) > tol) then 
       ko1=a22/d2
       ko2=a11/d2
       ko3=-a12/d2
    elseif (abs(d2) < tol) then 
       write(*,*) 'give3to4vect: d2< tol ?', (d2)
    else
       write(*,*) 'give3to4vect: d2=NaN ?', (d2)
    endif


    ! after that, I take a k3 vector an convolute it with that metric

    vaux=k3+ko1*a13*k1+ko2*a23*k2 + ko3*(a13*k2+a23*k1)


    dnorm=sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    if (abs(dnorm) > tol) then 
       v3a=vaux/dnorm
    elseif (abs(dnorm) < tol) then 
       write(*,*) 'give3to4vect: dnorm <tol ?', (dnorm)
    else
       write(*,*) 'give3to4vect: dnorm=NaN ?', (dnorm)
    endif

!^^^IFmp
!    vaux(1)=cone*0.5_qp
!    vaux(2)=cone*2.7_qp
!    vaux(3)=cone*3.2_qp
!    vaux(4)=cone*4.3_qp
!^^^ELSE
    vaux(1)=cmplx(0.5_qp,0.0_qp)
    vaux(2)=cmplx(2.7_qp,0.0_qp)
    vaux(3)=cmplx(3.2_qp,0.0_qp)
    vaux(4)=cmplx(4.3_qp,0.0_qp)
!^^^END

    a1v=sc(k1,vaux)
    a2v=sc(k2,vaux)

    vauxa=vaux+ko1*a1v*k1+ko2*a2v*k2+ ko3*(a1v*k2+a2v*k1)


    a3v = sc(vauxa,v3a)

    vaux = vauxa - a3v*v3a


    dnorm =  sc(vaux,vaux)
    dnorm=sqrt(dnorm)

    if (abs(dnorm) > tol) then 
       v4=vaux/dnorm
    elseif (abs(dnorm) < tol) then 
       write(*,*) 'give3to4vect: dnorm <tol ?', (dnorm)
    else
       write(*,*) 'give3to4vect: dnorm=NaN ?', (dnorm)
    endif



    v(1,:) = v1
    v(2,:) = v2
    v(3,:) = v3
    v(4,:) = v4

  end subroutine give3to4vect



  !----------------------------------------------------------------------------
  !  the procedure returns 4 v-vectors from two k-vectors, triangle basis
  !--------------------------------------------------------------------------

  subroutine give2to4vect(k1,k2,v)
    complex(qp),  intent(in) :: k1(:),k2(:)
    complex(qp),  intent(out) ::  v(:,:)
    ! ---------------------------------------------
    complex(qp) :: v1(4),v2(4),v3(4),v4(4)
    complex(qp) :: vaux(4)
    complex(qp) :: d2,ko1,ko2,ko3
    complex(qp) :: a11,a12,a22,dnorm
    complex(qp) :: a1v,a2v,a3v,a33

    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give3to4vect: wrong size of v'
    !if (size(k1) /= 4 .or. size(k2) /= 4) &
    !     &stop 'give1to4vect: wrong size of k1 or k2'


    call give2vect(k2,k1,v1)
    call give2vect(k1,k2,v2)
    d2 = sc(k1,v1)

    v1=v1/d2
    v2=v2/d2


    a11 = sc(k1,k1)
    a12 = sc(k1,k2)
    a22 = sc(k2,k2)

    d2=a12**2-a11*a22
    ko1=a22/d2
    ko2=a11/d2
    ko3=-a12/d2

    vaux(1)=cmplx(1.3_qp,0.0_qp)
    vaux(2)=cmplx(1.7_qp,0.0_qp)
    vaux(3)=cmplx(2.4_qp,0.0_qp)
    vaux(4)=cmplx(3.5_qp,0.0_qp)

    a1v = sc(k1,vaux)
    a2v = sc(k2,vaux)


    vaux=vaux+ko1*a1v*k1+ko2*a2v*k2+ko3*(a1v*k2 +a2v*k1)


    dnorm = sc(vaux,vaux)
    dnorm = sqrt(dnorm)

    v3=vaux/dnorm


    ! obtaining the forth vector
!^^^IFmp
!    vaux(1)=czero+2.1_qp
!    vaux(2)=czero+1.2_qp
!    vaux(3)=czero+3.4_qp
!    vaux(4)=czero+0.5_qp
!^^^ELSE
    vaux(1)=cmplx(2.1_qp,0.0_qp)
    vaux(2)=cmplx(1.2_qp,0.0_qp)
    vaux(3)=cmplx(3.4_qp,0.0_qp)
    vaux(4)=cmplx(0.5_qp,0.0_qp)
!^^^END

    a1v = sc(k1,vaux)
    a2v = sc(k2,vaux)

    vaux=vaux+ko1*a1v*k1+ko2*a2v*k2+ko3*(a1v*k2+a2v*k1)


    a3v =  sc(vaux,v3)
    a33 = sc(v3,v3)

    v3=vaux/dnorm
    vaux = vaux - a3v/a33*v3


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    v4 = vaux/dnorm

    v(1,:) = v1
    v(2,:) = v2
    v(3,:) = v3
    v(4,:) = v4

  end subroutine give2to4vect


  !--------------------------------------------------------------------
  ! the procedure returns 4 vectors from a not-light like vector 
  ! the first vector, v1, is parallel to k1 and it is not normalized
  !--------------------------------------------------------------------
  subroutine give1to4vect(k1,v)
    complex(qp), intent(in) :: k1(:)
    complex(qp), intent(out):: v(:,:)
    ! ----------------------------------------
    complex(qp) :: v1(4),v2(4),v3(4),v4(4)
    complex(qp) :: vaux(4),vaux1(4)
    complex(qp) :: v2aux(4)
    complex(qp) :: d2,d1,ko1,ko2,ko3
    complex(qp) :: a11,a12,a22,dnorm
    complex(qp) :: a1v,a2v,a3v,a33
    integer i

    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give3to4vect: wrong size of v'
    !if (size(k1) /= 4) stop 'give1to4vect: wrong size of k1'


    d1 = sc(k1,k1)
    do i=1,4
       v1(i)=k1(i)/d1
    enddo

!^^^IFmp
!    v2aux(1)=czero+1.3_qp
!    v2aux(2)=czero+1.7_qp
!    v2aux(3)=czero+2.3_qp
!    v2aux(4)=czero+3.4_qp
!^^^ELSE
    v2aux(1)=cmplx(1.3_qp,0.0_qp)
    v2aux(2)=cmplx(1.7_qp,0.0_qp)
    v2aux(3)=cmplx(2.3_qp,0.0_qp)
    v2aux(4)=cmplx(3.4_qp,0.0_qp)
!^^^END
    a1v = sc(v2aux,v1)

    do i=1,4
       v2aux(i)=v2aux(i)-a1v*k1(i)
    enddo

    dnorm = sc(v2aux,v2aux)
    dnorm=sqrt(dnorm)
    do i=1,4
       v2(i) = v2aux(i)/dnorm
    enddo

    a12=sc(v1,v2)
    a11=sc(v1,v1)
    a22=sc(v2,v2)
    d2=a12**2 - a11*a22

    ko1= a22/d2
    ko2= a11/d2
    ko3=-a12/d2

!^^^IFmp
!    vaux(1)=czero + 1.5_qp
!    vaux(2)=czero + 13.2_qp
!    vaux(3)=czero + 14.1_qp
!    vaux(4)=czero + 9.1_qp
!^^^ELSE
    vaux(1)=cmplx(1.5_qp,0.0_qp)
    vaux(2)=cmplx(13.2_qp,0.0_qp)
    vaux(3)=cmplx(14.1_qp,0.0_qp)
    vaux(4)=cmplx(9.1_qp,0.0_qp)  
!^^^END
    a1v = sc(v1,vaux)
    a2v = sc(v2,vaux)


    vaux=vaux+ko1*a1v*v1+ko2*a2v*v2+ko3*(a1v*v2+a2v*v1)


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    v3=vaux/dnorm

!^^^IFmp
!    vaux1(1)=czero + 2.4_qp
!    vaux1(2)=czero + 1.3_qp
!    vaux1(3)=czero + 3.45_qp
!    vaux1(4)=czero + 0.56_qp
!^^^ELSE
    vaux1(1)=cmplx(2.4_qp,0.0_qp)
    vaux1(2)=cmplx(1.3_qp,0.0_qp)
    vaux1(3)=cmplx(3.45_qp,0.0_qp)
    vaux1(4)=cmplx(0.56_qp,0.0_qp)
!^^^END

    a1v = sc(v1,vaux1)
    a2v = sc(v2,vaux1)


    vaux=vaux1+ko1*a1v*v1+ko2*a2v*v2 + ko3*(a1v*v2+a2v*v1)


    a3v =sc(vaux,v3)
    a33 =sc(v3,v3)

    vaux=vaux-a3v/a33*v3


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    v4=vaux/dnorm

    v(1,:) = v1
    v(2,:) = v2
    v(3,:) = v3
    v(4,:) = v4

  end subroutine give1to4vect

  !-------------------------------------------------------------------------
  ! for a light-like vector, this procedure returns the ``dual'' light-like 
  ! vector -- k1d and two unit vectors (v3,v4) in the plane transverse 
  ! to k1 & k1d 
  ! the first vector that is returned is k1
  !----------------------------------------------------------------------------
  subroutine give1to4vect_light(k1,v)
    complex(qp), intent(in)  ::  k1(:)
    complex(qp), intent(out) ::  v(:,:)
    ! ------------------------------------------       
    complex(qp) :: v3(4),v4(4)
    complex(qp) :: k1d(4)
    complex(qp) :: vaux(4),vaux1(4)
    complex(qp) :: d2,ko3
    complex(qp) :: dnorm
    complex(qp) :: a1v,a2v,a33,av3


    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give1to4vectlight: wrong size of v'
    !if (size(k1) /= 4) stop 'give1to4vectlight: wrong size of k1'

    k1d=-k1
    k1d(1)=k1(1)

    d2 = sc(k1d,k1)

    if (abs(d2) < tol) then 
       k1d(1)=1.0_qp+ci*0.0_qp
       k1d(2)=1.0_qp/sqrt3 + ci*0.0_qp
       k1d(3)=1.0_qp/sqrt3 + ci*0.0_qp
       k1d(4)=1.0_qp/sqrt3 + ci*0.0_qp
    endif

    d2 = sc(k1d,k1)

    if (abs(d2).lt. sq2tol) then 
       k1d(1)=1.0_qp + ci*0.0_qp
       k1d(2)=1.0_qp/sqrt3 + ci*0.0_qp
       k1d(3)=sqrt2/sqrt3 + ci*0.0_qp
       k1d(4)=0.0_qp + ci*0.0_qp
    endif

    d2 = sc(k1d,k1)

    ko3=-1.0_qp/d2

    vaux(1)=(1.5_qp + ci*0.0_qp)
    vaux(2)=(13.2_qp + ci*0.0_qp)
    vaux(3)=(14.1_qp + ci*0.0_qp)
    vaux(4)=(9.1_qp + ci*0.0_qp)

    a1v = sc(k1,vaux)
    a2v =  sc(k1d,vaux)


    vaux=vaux+ko3*(a1v*k1d+a2v*k1)


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)

    v3=vaux/dnorm 


    ! obtaining the forth vector

    vaux1(1)=(2.4_qp + ci*0.0_qp)
    vaux1(2)=(1.35_qp + ci*0.0_qp)
    vaux1(3)=(3.61_qp + ci*0.0_qp)
    vaux1(4)=(0.54_qp + ci*0.0_qp)  


    a1v=sc(k1,vaux1)
    a2v=sc(k1d,vaux1)         

    vaux=vaux1+ko3*(a1v*k1d+a2v*k1)

    av3 =sc(vaux,v3)
    a33 =sc(v3,v3)
    vaux=vaux-av3/a33*v3

    dnorm=sc(vaux,vaux)
    dnorm=sqrt(dnorm)

    v4=vaux/dnorm

    v(1,:)=k1
    v(2,:)=k1d
    v(3,:)=v3
    v(4,:)=v4

  end subroutine give1to4vect_light

end module qpvvn
