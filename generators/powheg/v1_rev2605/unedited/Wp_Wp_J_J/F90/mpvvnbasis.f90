!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECvvnbasis.f90.

! --- module for vermasseren-van neerven basis
module mpvvn
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use types; use consts_mp
  use mpgram; use mpaux_functions
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
    type(mp_complex), intent(in) :: k1(:),k2(:),k3(:),k4(:)
    type(mp_complex), intent(out) :: v(:,:)
    ! -------------------------------------------------      
    type(mp_complex) :: v1(4),v2(4),v3(4),v4(4)
    type(mp_complex) :: d4

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
       write(*,*) 'give4to4vect: d4< tol ?', mptodp(d4), mptodp(tol)
    else
       write(*,*) 'give4to4vect: d4=NaN ?', mptodp(d4)
    endif

  end subroutine give4to4vect



  !-------------------------------------------------------------------------
  ! The procedure that returns 4 vectors from 3 vectors; `box'' basis
  !-------------------------------------------------------------------------
  subroutine give3to4vect(k1,k2,k3,v)
    type(mp_complex), intent(in) :: k1(:),k2(:),k3(:)
    type(mp_complex), intent(out) :: v(:,:) 
    ! -------------------------------------------------      
    type(mp_complex) :: v1(4),v2(4),v3(4),v4(4)
    type(mp_complex) :: vaux(4),vauxa(4),v3a(4)
    type(mp_complex) :: d3,d2,ko1,ko2,ko3
    type(mp_complex) :: a11,a12,a22,a13,a23,dnorm,a2v,a1v
    type(mp_complex) :: a3v

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
       write(*,*) 'give3to4vect: d3< tol ?', mptodp(d3)
    else
       write(*,*) 'give3to4vect: d3=NaN ?', mptodp(d3)
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
       write(*,*) 'give3to4vect: d2< tol ?', mptodp(d2)
    else
       write(*,*) 'give3to4vect: d2=NaN ?', mptodp(d2)
    endif


    ! after that, I take a k3 vector an convolute it with that metric

    vaux=k3+ko1*a13*k1+ko2*a23*k2 + ko3*(a13*k2+a23*k1)


    dnorm=sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    if (abs(dnorm) > tol) then 
       v3a=vaux/dnorm
    elseif (abs(dnorm) < tol) then 
       write(*,*) 'give3to4vect: dnorm <tol ?', mptodp(dnorm)
    else
       write(*,*) 'give3to4vect: dnorm=NaN ?', mptodp(dnorm)
    endif

!^^^IFmp
    vaux(1)=cone*mpreal('0.5')
    vaux(2)=cone*mpreal('2.7')
    vaux(3)=cone*mpreal('3.2')
    vaux(4)=cone*mpreal('4.3')
!^^^ELSE
!    vaux(1)=cmplx(mpreal('0.5'),mpreal('0.0'))
!    vaux(2)=cmplx(mpreal('2.7'),mpreal('0.0'))
!    vaux(3)=cmplx(mpreal('3.2'),mpreal('0.0'))
!    vaux(4)=cmplx(mpreal('4.3'),mpreal('0.0'))
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
       write(*,*) 'give3to4vect: dnorm <tol ?', mptodp(dnorm)
    else
       write(*,*) 'give3to4vect: dnorm=NaN ?', mptodp(dnorm)
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
    type(mp_complex),  intent(in) :: k1(:),k2(:)
    type(mp_complex),  intent(out) ::  v(:,:)
    ! ---------------------------------------------
    type(mp_complex) :: v1(4),v2(4),v3(4),v4(4)
    type(mp_complex) :: vaux(4)
    type(mp_complex) :: d2,ko1,ko2,ko3
    type(mp_complex) :: a11,a12,a22,dnorm
    type(mp_complex) :: a1v,a2v,a3v,a33

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

    vaux(1)=cmplx(mpreal('1.3'),mpreal('0.0'))
    vaux(2)=cmplx(mpreal('1.7'),mpreal('0.0'))
    vaux(3)=cmplx(mpreal('2.4'),mpreal('0.0'))
    vaux(4)=cmplx(mpreal('3.5'),mpreal('0.0'))

    a1v = sc(k1,vaux)
    a2v = sc(k2,vaux)


    vaux=vaux+ko1*a1v*k1+ko2*a2v*k2+ko3*(a1v*k2 +a2v*k1)


    dnorm = sc(vaux,vaux)
    dnorm = sqrt(dnorm)

    v3=vaux/dnorm


    ! obtaining the forth vector
!^^^IFmp
    vaux(1)=czero+mpreal('2.1')
    vaux(2)=czero+mpreal('1.2')
    vaux(3)=czero+mpreal('3.4')
    vaux(4)=czero+mpreal('0.5')
!^^^ELSE
!    vaux(1)=cmplx(mpreal('2.1'),mpreal('0.0'))
!    vaux(2)=cmplx(mpreal('1.2'),mpreal('0.0'))
!    vaux(3)=cmplx(mpreal('3.4'),mpreal('0.0'))
!    vaux(4)=cmplx(mpreal('0.5'),mpreal('0.0'))
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
    type(mp_complex), intent(in) :: k1(:)
    type(mp_complex), intent(out):: v(:,:)
    ! ----------------------------------------
    type(mp_complex) :: v1(4),v2(4),v3(4),v4(4)
    type(mp_complex) :: vaux(4),vaux1(4)
    type(mp_complex) :: v2aux(4)
    type(mp_complex) :: d2,d1,ko1,ko2,ko3
    type(mp_complex) :: a11,a12,a22,dnorm
    type(mp_complex) :: a1v,a2v,a3v,a33
    integer i

    if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
         &stop 'give3to4vect: wrong size of v'
    if (size(k1) /= 4) stop 'give1to4vect: wrong size of k1'


    d1 = sc(k1,k1)
    do i=1,4
       v1(i)=k1(i)/d1
    enddo

!^^^IFmp
    v2aux(1)=czero+mpreal('1.3')
    v2aux(2)=czero+mpreal('1.7')
    v2aux(3)=czero+mpreal('2.3')
    v2aux(4)=czero+mpreal('3.4')
!^^^ELSE
!    v2aux(1)=cmplx(mpreal('1.3'),mpreal('0.0'))
!    v2aux(2)=cmplx(mpreal('1.7'),mpreal('0.0'))
!    v2aux(3)=cmplx(mpreal('2.3'),mpreal('0.0'))
!    v2aux(4)=cmplx(mpreal('3.4'),mpreal('0.0'))
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
    vaux(1)=czero + mpreal('1.5')
    vaux(2)=czero + mpreal('13.2')
    vaux(3)=czero + mpreal('14.1')
    vaux(4)=czero + mpreal('9.1')
!^^^ELSE
!    vaux(1)=cmplx(mpreal('1.5'),mpreal('0.0'))
!    vaux(2)=cmplx(mpreal('13.2'),mpreal('0.0'))
!    vaux(3)=cmplx(mpreal('14.1'),mpreal('0.0'))
!    vaux(4)=cmplx(mpreal('9.1'),mpreal('0.0'))  
!^^^END
    a1v = sc(v1,vaux)
    a2v = sc(v2,vaux)


    vaux=vaux+ko1*a1v*v1+ko2*a2v*v2+ko3*(a1v*v2+a2v*v1)


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)
    v3=vaux/dnorm

!^^^IFmp
    vaux1(1)=czero + mpreal('2.4')
    vaux1(2)=czero + mpreal('1.3')
    vaux1(3)=czero + mpreal('3.45')
    vaux1(4)=czero + mpreal('0.56')
!^^^ELSE
!    vaux1(1)=cmplx(mpreal('2.4'),mpreal('0.0'))
!    vaux1(2)=cmplx(mpreal('1.3'),mpreal('0.0'))
!    vaux1(3)=cmplx(mpreal('3.45'),mpreal('0.0'))
!    vaux1(4)=cmplx(mpreal('0.56'),mpreal('0.0'))
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
    type(mp_complex), intent(in)  ::  k1(:)
    type(mp_complex), intent(out) ::  v(:,:)
    ! ------------------------------------------       
    type(mp_complex) :: v3(4),v4(4)
    type(mp_complex) :: k1d(4)
    type(mp_complex) :: vaux(4),vaux1(4)
    type(mp_complex) :: d2,ko3
    type(mp_complex) :: dnorm
    type(mp_complex) :: a1v,a2v,a33,av3


    !if (size(v,dim=1) /= 4 .or. size(v,dim=2) /= 4) &
    !     &stop 'give1to4vectlight: wrong size of v'
    !if (size(k1) /= 4) stop 'give1to4vectlight: wrong size of k1'

    k1d=-k1
    k1d(1)=k1(1)

    d2 = sc(k1d,k1)

    if (abs(d2) < tol) then 
       k1d(1)=mpreal('1.0')+ci*mpreal('0.0')
       k1d(2)=mpreal('1.0')/sqrt3 + ci*mpreal('0.0')
       k1d(3)=mpreal('1.0')/sqrt3 + ci*mpreal('0.0')
       k1d(4)=mpreal('1.0')/sqrt3 + ci*mpreal('0.0')
    endif

    d2 = sc(k1d,k1)

    if (abs(d2).lt. sq2tol) then 
       k1d(1)=mpreal('1.0') + ci*mpreal('0.0')
       k1d(2)=mpreal('1.0')/sqrt3 + ci*mpreal('0.0')
       k1d(3)=sqrt2/sqrt3 + ci*mpreal('0.0')
       k1d(4)=mpreal('0.0') + ci*mpreal('0.0')
    endif

    d2 = sc(k1d,k1)

    ko3= mpreal('-1.0')/d2

    vaux(1)=(mpreal('1.5') + ci*mpreal('0.0'))
    vaux(2)=(mpreal('13.2') + ci*mpreal('0.0'))
    vaux(3)=(mpreal('14.1') + ci*mpreal('0.0'))
    vaux(4)=(mpreal('9.1') + ci*mpreal('0.0'))

    a1v = sc(k1,vaux)
    a2v =  sc(k1d,vaux)


    vaux=vaux+ko3*(a1v*k1d+a2v*k1)


    dnorm = sc(vaux,vaux)
    dnorm=sqrt(dnorm)

    v3=vaux/dnorm 


    ! obtaining the forth vector

    vaux1(1)=(mpreal('2.4') + ci*mpreal('0.0'))
    vaux1(2)=(mpreal('1.35') + ci*mpreal('0.0'))
    vaux1(3)=(mpreal('3.61') + ci*mpreal('0.0'))
    vaux1(4)=(mpreal('0.54') + ci*mpreal('0.0'))  


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

end module mpvvn
