!  procedures to match various cuts

module match1
  implicit none
  private

  public :: match, mismatch

contains



  subroutine match(m1,m2,lcut,fl,N5,Lc5,F5,N4,Lc4,F4,&
       &     N3,Lc3,F3,N2,Lc2,F2,ns,lmatch)
    integer, intent(in)   :: m1,m2,lcut(:)
    integer, intent(in)   :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in) :: F5(:,:)*3,F4(:,:)*3,F3(:,:)*3,F2(:,:)*3
    character, intent(in) :: fl(:)*3
    integer, intent(in)   :: N5,N4,N3,N2
    integer, intent(out)  :: ns, lmatch(:)
    ! ----------------------------------------------------------------      
    integer :: i, j1,j2,xtot, Nmax

    ns = 0

    if (m2==5) then 
       Nmax = N5
    elseif (m2==4) then 
       Nmax = N4
    elseif (m2==3) then 
       Nmax = N3
    elseif (m2==2) then 
       Nmax = N2
    else
       Nmax = 0 
       stop 'match: m2 out of range'
    endif 

    do i=1,Nmax

       xtot=0

       do j1=1,m1 
          do j2=1,m2

             if (m2 == 5) then 
                if ((lcut(j1) == Lc5(i,j2)).and.(fl(j1) == F5(i,j2))) then
                   xtot= xtot + 1
                endif
             elseif (m2 == 4) then  
                if ((lcut(j1) == Lc4(i,j2)).and.(fl(j1) == F4(i,j2))) then
                   xtot= xtot + 1
                endif
             elseif (m2 == 3) then 
                if ((lcut(j1) == Lc3(i,j2)).and.(fl(j1) == F3(i,j2))) then
                   xtot= xtot + 1
                endif
             elseif (m2 == 2) then 
                if ((lcut(j1) == Lc2(i,j2)).and.(fl(j1) == F2(i,j2))) then
                   xtot= xtot + 1
                endif
             else
                stop 'match: m2 out of range' 
             endif

          enddo
       enddo

       if (xtot == m1) then 
          ns = ns+1
          lmatch(ns)=i
       endif

    enddo

  end subroutine match


  subroutine mismatch(m1,m2,lcut,ia,N5,Lc5,N4,Lc4,N3,Lc3,N2,Lc2,lpos)
    integer, intent(in)  :: m1,m2, lcut(:), ia
    integer, intent(in)  :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    integer, intent(in)  :: N5,N4,N3,N2 ! not used.... 
    integer, intent(out) :: lpos(m2-m1)
    ! -------------------------------------------------------------       
    integer :: j1,j2,xtot,ns, Lt(m2)


    ns=0

    if (m2==5) then 
       Lt = Lc5(ia,:)
    elseif (m2==4) then 
       Lt = Lc4(ia,:)
    elseif (m2==3) then 
       Lt = Lc3(ia,:)
    elseif (m2==2) then 
       Lt = Lc2(ia,:)
    else
       stop 'mismatch: m2 out of range'
    endif

    do j1=1,m2
       xtot = 1
       do j2=1,m1

          if (Lt(j1)==lcut(j2)) then 
             xtot = 0*xtot
          else
             xtot = xtot
          endif


       enddo
       if (xtot==1) then 
          ns=ns+1
          lpos(ns) = j1
       endif

    enddo

  end subroutine mismatch


end module match1
