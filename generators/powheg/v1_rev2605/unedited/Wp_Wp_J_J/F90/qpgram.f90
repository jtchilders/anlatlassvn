!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECgram.f90.

! -- module for gram-determinant and vvn basis
module qpgram
  use types; use consts_qp
  use qpaux_functions 
  implicit none
  private

  public :: delta2,delta3,delta4
  public :: give2vect,give3vect,give4vect

contains


  subroutine delta2(k1,k2,d2)
    complex(qp), intent(in)  :: k1(:), k2(:)
    complex(qp), intent(out) :: d2             
    ! ---------------------------------------

    if (size(k1) /= 4 .or. size(k2) /= 4) stop 'delta2: wrong size of k1 or k2'


    d2 = -(-(k1(1)*k2(1)) + k1(2)*k2(2) + k1(3)*k2(3) + &
         &      k1(4)*k2(4))**2 + &
         &  (k1(1)**2 - k1(2)**2 - k1(3)**2 - k1(4)**2)*&
         &   (k2(1)**2 - k2(2)**2 - k2(3)**2 - k2(4)**2)

  end subroutine delta2



  subroutine give2vect(k2,k1,v)
    complex(qp), intent(in)  :: k1(:), k2(:)
    complex(qp), intent(out) :: v(:)
    ! --------------------------------------

    if (size(k1) /= 4 .or. size(k2) /=4) stop 'give2vect: wrong size of k1 or k2'
    if (size(v) /= 4 ) stop 'give2vect: wrong size of v'


    v = sc(k2,k2)*k1-sc(k1,k2)*k2

  end subroutine give2vect


  subroutine delta3(k1,k2,k3,d3)
    complex(qp), intent(in)  :: k1(:), k2(:), k3(:)
    complex(qp), intent(out) :: d3

    if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /= 4) &
         &stop 'delta3: wrong size of k1, k2 or k3'


    d3 =    -((k2(1)**2 - k2(2)**2 - k2(3)**2 - k2(4)**2)*&
         &     (-(k1(1)*k3(1)) + k1(2)*k3(2) + k1(3)*k3(3) +    &
         &        k1(4)*k3(4))**2) +                            &
         &  2*(k1(1)*k2(1) - k1(2)*k2(2) - k1(3)*k2(3) -        &
         &     k1(4)*k2(4))*                                    &
         &   (k1(1)*k3(1) - k1(2)*k3(2) - k1(3)*k3(3) -         &
         &     k1(4)*k3(4))*                                    &
         &   (k2(1)*k3(1) - k2(2)*k3(2) - k2(3)*k3(3) -         &
         &     k2(4)*k3(4)) -                                   &
         &  (k1(1)**2 - k1(2)**2 - k1(3)**2 - k1(4)**2)*        &
         &   (-(k2(1)*k3(1)) + k2(2)*k3(2) + k2(3)*k3(3) +      &
         &      k2(4)*k3(4))**2 -                               &
         &  (-(k1(1)*k2(1)) + k1(2)*k2(2) + k1(3)*k2(3) +       &
         &      k1(4)*k2(4))**2*                                &
         &   (k3(1)**2 - k3(2)**2 - k3(3)**2 - k3(4)**2) +      &
         &  (k1(1)**2 - k1(2)**2 - k1(3)**2 - k1(4)**2)*        &
         &   (k2(1)**2 - k2(2)**2 - k2(3)**2 - k2(4)**2)*       &
         &   (k3(1)**2 - k3(2)**2 - k3(3)**2 - k3(4)**2)

  end subroutine delta3


  subroutine give3vect(k2,k3,k1,v)
    complex(qp), intent(in)   :: k1(:), k2(:), k3(:)
    complex(qp),  intent(out) :: v(:)
    ! ---------------------------------------------------       

    if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /= 4) &
         &stop 'delta3: wrong size of k1, k2 or k3'
    if (size(v) /= 4 ) stop 'give3vect: wrong size of v'


    v =  k3*(sc(k2,k1)*sc(k2,k3) - sc(k2,k2)*sc(k3,k1)) + &
         &  k2*(sc(k2,k3)*sc(k3,k1) - sc(k2,k1)*sc(k3,k3)) + &
         &  k1*(-sc(k2,k3)**2 + sc(k2,k2)*sc(k3,k3))

  end subroutine give3vect

  subroutine give4vect(k2,k3,k4,k1,v)
    complex(qp), intent(in) :: k1(:), k2(:), k3(:), k4(:)
    complex(qp), intent(out) :: v(:)

    if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /= 4 .or. size(k4) /= 4) &
         &stop 'delta3: wrong size of k1, k2,k3 or k4'
    if (size(v) /= 4 ) stop 'give3vect: wrong size of v'


    v = -(k1*sc(k2,k4)**2*sc(k3,k3)) +                          &
         &  two*k1*sc(k2,k3)*sc(k2,k4)*sc(k3,k4) -                       &
         &  k2*sc(k2,k4)*sc(k3,k1)*sc(k3,k4) +                         &
         &  k2*sc(k2,k1)*sc(k3,k4)**2 - k1*sc(k2,k2)*sc(k3,k4)**2 +    &
         &  k2*sc(k2,k4)*sc(k3,k3)*sc(k4,k1) -                         &
         &  k2*sc(k2,k3)*sc(k3,k4)*sc(k4,k1) +                         &
         &  k4*(sc(k2,k1)*sc(k2,k4)*sc(k3,k3) -                        &
         &     sc(k2,k3)*(sc(k2,k4)*sc(k3,k1) + sc(k2,k1)*sc(k3,k4)) + &
         &     sc(k2,k3)**2*sc(k4,k1) +                                &
         &     sc(k2,k2)*(sc(k3,k1)*sc(k3,k4) - sc(k3,k3)*sc(k4,k1))) -& 
         &  k1*sc(k2,k3)**2*sc(k4,k4) +                                &
         &  k2*sc(k2,k3)*sc(k3,k1)*sc(k4,k4) -                         &
         &  k2*sc(k2,k1)*sc(k3,k3)*sc(k4,k4) +                         &
         &  k1*sc(k2,k2)*sc(k3,k3)*sc(k4,k4) +                         &
         &  k3*(sc(k2,k4)**2*sc(k3,k1) -                               &
         &     sc(k2,k4)*(sc(k2,k1)*sc(k3,k4) + sc(k2,k3)*sc(k4,k1)) + &
         &     sc(k2,k1)*sc(k2,k3)*sc(k4,k4) +                         &
         &     sc(k2,k2)*(sc(k3,k4)*sc(k4,k1) - sc(k3,k1)*sc(k4,k4)))


  end subroutine give4vect

  subroutine delta4(k1,k2,k3,k4,d4)
    complex(qp), intent(in)  :: k1(:), k2(:), k3(:), k4(:)
    complex(qp), intent(out) ::  d4

    if (size(k1) /= 4 .or. size(k2) /= 4 .or. size(k3) /= 4 .or. size(k4) /= 4) &
         &stop 'delta4: size k1,k2,k3 or k4 wrong'


    d4 =  -(k1(2)*k2(4)*k3(3)*k4(1) -              &
         &     k1(2)*k2(3)*k3(4)*k4(1) -                  &
         &     k1(1)*k2(4)*k3(3)*k4(2) +                  &
         &     k1(1)*k2(3)*k3(4)*k4(2) -                  &
         &     k1(2)*k2(4)*k3(1)*k4(3) +                  &
         &     k1(1)*k2(4)*k3(2)*k4(3) +                  &
         &     k1(2)*k2(1)*k3(4)*k4(3) -                  &
         &     k1(1)*k2(2)*k3(4)*k4(3) +                  &
         &     k1(4)*(-(k2(2)*k3(3)*k4(1)) +              &
         &        k2(1)*k3(3)*k4(2) +                     &
         &        k2(3)*(k3(2)*k4(1) - k3(1)*k4(2)) +     &
         &        k2(2)*k3(1)*k4(3) - k2(1)*k3(2)*k4(3)) +&
         &     k1(2)*k2(3)*k3(1)*k4(4) -                  &
         &     k1(1)*k2(3)*k3(2)*k4(4) -                  &
         &     k1(2)*k2(1)*k3(3)*k4(4) +                  &
         &     k1(1)*k2(2)*k3(3)*k4(4) +                  &
         &     k1(3)*(-(k2(4)*k3(2)*k4(1)) +              &
         &        k2(2)*k3(4)*k4(1) + k2(4)*k3(1)*k4(2) - &
         &        k2(1)*k3(4)*k4(2) - k2(2)*k3(1)*k4(4) + &
         &        k2(1)*k3(2)*k4(4)))**2


  end subroutine delta4

end module qpgram
