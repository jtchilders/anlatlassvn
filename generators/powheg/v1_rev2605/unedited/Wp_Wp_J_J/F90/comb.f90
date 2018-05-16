module comb
  implicit none
  private


  public :: getallcuts,getcutnumb,getallcuts_nf,getcutnumb_nf

contains


  function fac(N)
    integer :: N, fac,i

    if (N == 0) then 
       fac = 1
    else 
       fac = 1
       do i=1,N
          fac = i*fac
       enddo
    endif

  end function fac


  subroutine getcutnumb(lncuts,npoint)
    integer, intent(in)  :: npoint
    integer, intent(out) :: lncuts(5)
    
    if (Npoint >= 5) then 
       lncuts(5) = fac(Npoint)/fac(5)/fac(Npoint-5)
    else 
       lncuts(5) = 0
    endif
    
    lncuts(4) = fac(Npoint)/fac(4)/fac(Npoint-4)
    lncuts(3) = fac(Npoint)/fac(3)/fac(Npoint-3)
    lncuts(2) = fac(Npoint)/fac(2)/fac(Npoint-2)
    lncuts(1) = Npoint
    
    
  end subroutine getcutnumb

  subroutine getcutnumb_nf(lncuts,nterm,npoint)
    integer, intent(in)  :: nterm,npoint
    integer, intent(out) :: lncuts(5)
    integer ::  par,i
    
    par = Npoint - (Nterm - 1)
    
    do i=1,5
       lncuts(i) = fac(par)/fac(i)/fac(par-i)
    enddo
    
  end subroutine getcutnumb_nf
  




  subroutine getallcuts(lncuts,npoint,Ac5,Ac4,Ac3,Ac2,Ac1)
    integer, intent(in) :: lncuts(5),npoint 
    integer, intent(out) :: Ac5(lncuts(5),5) 
    integer, intent(out) :: Ac4(lncuts(4),4),Ac3(lncuts(3),3)
    integer, intent(out) :: Ac2(lncuts(2),2), Ac1(lncuts(1),1)
    integer :: count,i5,i4,i3,i2,i1


    count = 1

    do i1=1,Npoint-4
       do i2 = i1+1,Npoint-3
          do i3 = i2+1, Npoint-2
             do i4 = i3+1,Npoint-1
                do  i5 = i4+1,Npoint 
                   Ac5(count,1) = i1
                   Ac5(count,2) = i2
                   Ac5(count,3) = i3
                   Ac5(count,4) = i4
                   Ac5(count,5) = i5
                   
                   count = count + 1
                enddo
             enddo
          enddo
       enddo
    enddo
    if (count-1 > lncuts(5)) stop 'getallcuts count > lncuts(5)'


    count = 1

    do i1=1,Npoint-3
       do i2 = i1+1,Npoint-2
          do i3 = i2+1, Npoint-1
             do i4 = i3+1,Npoint
                Ac4(count,1) = i1
                Ac4(count,2) = i2
                Ac4(count,3) = i3
                Ac4(count,4) = i4

                count = count + 1
             enddo
          enddo
       enddo
    enddo
    if (count-1 > lncuts(4)) stop 'getallcuts count > lncuts(4)'

    count = 1

    do i1=1,Npoint-2
       do i2 = i1+1,Npoint-1
          do i3 = i2+1, Npoint

             Ac3(count,1) = i1
             Ac3(count,2) = i2
             Ac3(count,3) = i3
             count = count + 1
          enddo
       enddo
    enddo
    if (count-1 > lncuts(3)) stop 'getallcuts count > lncuts(3)'

    count = 1

    do i1=1,Npoint-1
       do i2 = i1+1,Npoint

          Ac2(count,1) = i1
          Ac2(count,2) = i2

          count = count + 1
       enddo
    enddo
    if (count-1 > lncuts(2)) stop 'getallcuts count > lncuts(2)'

    count =1
    do i1 = 1,Npoint 
       Ac1(count,1) = i1
       count = count + 1
    enddo
    if (count-1 > lncuts(1)) stop 'getallcuts count > lncuts(1)'

  end subroutine getallcuts

  
  subroutine getallcuts_nf(lncuts,npoint,nterm,Ac5,Ac4,Ac3,Ac2,Ac1)
    integer, intent(in) :: lncuts(5),npoint,nterm 
    integer, intent(out) :: Ac5(lncuts(5),5) 
    integer, intent(out) :: Ac4(lncuts(4),4),Ac3(lncuts(3),3)
    integer, intent(out) :: Ac2(lncuts(2),2), Ac1(lncuts(1),1)
    integer :: count,i5,i4,i3,i2,i1,par

    count = 1

    do i1=1,Npoint-4
       do i2 = i1+1,Npoint-3
          do i3 = i2+1, Npoint-2
             do i4 = i3+1,Npoint-1
                do  i5 = i4+1,Npoint 

                   if ((i1.eq.1.or.i1.gt.Nterm).and. &
                        & (i2.eq.1.or.i2.gt.Nterm).and. &
                        & (i3.eq.1.or.i3.gt.Nterm).and. &
                        & (i4.eq.1.or.i4.gt.Nterm).and. &
                        & (i5.eq.1.or.i5.gt.Nterm) ) then 

                      Ac5(count,1) = i1
                      Ac5(count,2) = i2
                      Ac5(count,3) = i3
                      Ac5(count,4) = i4
                      Ac5(count,5) = i5

                      count = count + 1

                   endif

                enddo
             enddo
          enddo
       enddo
    enddo
    if (count-1 > lncuts(5)) stop 'getallcuts_nf count > lncuts(5)'

    count = 1

    do i1=1,Npoint-3
       do i2 = i1+1,Npoint-2
          do i3 = i2+1, Npoint-1
             do i4 = i3+1,Npoint


                if ((i1.eq.1.or.i1.gt.Nterm).and. &
                   &(i2.eq.1.or.i2.gt.Nterm).and. &
                   &(i3.eq.1.or.i3.gt.Nterm).and. &
                   &(i4.eq.1.or.i4.gt.Nterm) ) then 

                   Ac4(count,1) = i1
                   Ac4(count,2) = i2
                   Ac4(count,3) = i3
                   Ac4(count,4) = i4

                   count = count + 1

                endif

             enddo
          enddo
       enddo
    enddo
    if (count-1 > lncuts(4)) stop 'getallcuts_nf count > lncuts(4)'

    count = 1

    do i1=1,Npoint-2
       do i2 = i1+1,Npoint-1
          do i3 = i2+1, Npoint


             if ((i1.eq.1.or.i1.gt.Nterm).and. &
                &(i2.eq.1.or.i2.gt.Nterm).and. &
                &(i3.eq.1.or.i3.gt.Nterm) )then 

                Ac3(count,1) = i1
                Ac3(count,2) = i2
                Ac3(count,3) = i3

                count = count + 1

             endif


          enddo
       enddo
    enddo
    if (count-1 > lncuts(3)) stop 'getallcuts_nf count > lncuts(3)'

    count = 1

    do i1=1,Npoint-1
       do i2 = i1+1,Npoint

          if ((i1.eq.1.or.i1.gt.Nterm).and. &
             &(i2.eq.1.or.i2.gt.Nterm) ) then 

             Ac2(count,1) = i1
             Ac2(count,2) = i2

             count = count + 1

          endif


       enddo
    enddo
    if (count-1 > lncuts(2)) stop 'getallcuts_nf count > lncuts(2)'

    count =1
    do i1 = Nterm,Npoint 
       Ac1(count,1) = i1
       count = count + 1
    enddo
    if (count-1 > lncuts(1)) stop 'getallcuts_nf count > lncuts(1)'

  end subroutine getallcuts_nf



end module comb




