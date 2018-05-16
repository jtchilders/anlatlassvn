!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECinitWpWp.f90.

module qpinitWpWp

  use types
  use consts_qp 
  use qpaux_functions
  use comb
  use define_ampl
  use qpamplitude
  use qpglobal

  implicit none

  private

  public :: pampl_count_WpWp, pamplWpWp, oneloopdiv


contains


  subroutine pamplWpWp(pext,pl,pa,hl,lncuts,N5,N4,N3,tN2,N2,N1,&
       &Lc5,F5,Yc5,Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2,tree,mu,div2,div1,div3,&
       &outgoing,crossed_gl)
    complex(qp), intent(in)  :: pext(:,:),pl(:,:),pa(:,:)
    integer, intent(in)    :: hl(:)
    integer, intent(in)    :: lncuts(5)
    integer, intent(in)    :: N1,N2,N3,N4,N5,tN2
    logical, intent(in)    :: crossed_gl
    logical, intent(in), optional :: outgoing(:) 
    ! ----------------------------------------      
    integer, intent(out)   :: Lc5(:,:)
    character, intent(out) :: F5(:,:)*3
    integer, intent(out)   :: Yc5(:,:)
    integer, intent(out)   :: Lc4(:,:)
    character, intent(out) :: F4(:,:)*3
    integer, intent(out)   :: Yc4(:,:)
    integer, intent(out)   :: Lc3(:,:)
    character, intent(out) :: F3(:,:)*3
    integer, intent(out)   :: Yc3(:,:)
    integer, intent(out)   :: Lc2(:,:)
    character, intent(out) :: F2(:,:)*3
    integer, intent(out)   :: Yc2(:,:)
    logical:: nodum
    complex(qp), intent(out) :: tree
    real(qp), intent(in)     :: mu
    complex(qp), intent(out) :: div2, div1, div3
    ! -------------------------------------------------
    integer :: i,j,Np, jglob, ij
    integer :: j1, j2, ic, i1, i2, j3, pup, icount, pos_wwp, pos_wwm,pupp
    integer :: xxx, iq, iw
    integer :: aN2, aN3, aN4, aN5
    integer :: xN2, xN3, xN4, xN5
    character :: FFc5(Ncut,5)*3, FFc4(Ncut,4)*3
    character :: FFc3(Ncut,3)*3, FFc2(Ncut,2)*3 
    integer :: Ac5(lncuts(5),5)
    integer :: Ac4(lncuts(4),4), Ac3(lncuts(3),3) 
    integer :: Ac2(lncuts(2),2), Ac1(lncuts(1),1)

    character, dimension(:,:), allocatable :: Lf*3
    character, dimension(:,:), allocatable :: Lg*3
    integer, dimension(:,:), allocatable :: Yf
    integer, dimension(:,:), allocatable :: Yg
    integer :: Lc2a(tN2,2)
    character :: F2a(tN2,2)*3
    integer :: Yc2a(tN2,1)
    complex(qp) :: e(size(pext,dim=1),4)
    complex(qp) :: sp(4,Npoint), p(4,Npoint)
    character :: fl1(Npoint)*3
    integer :: ng,nq,nW
    complex(qp) :: vg(4), vs(4)


    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn


    ! Case where gluons cross inbetween quark lines -- swap ubar and v and 
    ! helicities. Momenta swapped in qqqqampl_v.

    if (crossed_gl) then 
       e(5,:) = -v0(pext(5,:), hl(6))
       e(6,:) = ubar0(pext(6,:),hl(5))
    else
       e(5,:) = v0(pext(5,:), hl(5))
       e(6,:) = ubar0(pext(6,:),hl(6))
    endif

    e(3,:) = pol_dk2mom(pl(1,:),pa(1,:),hl(3),outgoing(3))  !w
    e(4,:) = pol_dk2mom(pl(2,:),pa(2,:), hl(4), outgoing(4))


    ! -- labels for external particles. In this case, with Ncut always =1,
    ! it is simplest to do this by hand.


    if (case_b1) then
       Lab_ex(Ncut,:) = (/'top', 'wwp', 'bot', 'tp1', 'wwp', 'bt1'/)
       Lab_in(Ncut,:) = (/'glu', 'top', 'bot', 'glu', 'top', 'bot'/) 


    elseif (case_b2) then
       Lab_ex(Ncut,:) = (/'top', 'wwp', 'bt1','wwp', 'tp1', 'bot'/)
       Lab_in(Ncut,:) = (/'glu','top', 'bot', 'dum','dum','bot'/) 
    endif
    ! -- external momenta ordered

    do i=1,Ncut

       mom(i,1,:) = pext(1,:) 
       hel(i,1,:) = e(1,:)
       identity(i,1) = 1


       ic = 0
       iq = 0
       iw = 0

       do j=2,Npoint

          if (Lab_ex(i,j) == 'top' .or. Lab_ex(i,j) == 'bot' .or. Lab_ex(i,j) == 'str') then 
             mom(i,j,:) = pext(2,:)
             hel(i,j,:) = e(2,:)          
             identity(i,j) = 2
          elseif (Lab_ex(i,j) == 'tp1') then
             mom(i,j,:) = pext(5,:)
             hel(i,j,:) = e(5,:)
             identity(i,j) = 16
             iq = iq+1
             Lab_ex(i,j) ='top'

          elseif (Lab_ex(i,j) == 'bt1') then
             mom(i,j,:) = pext(6,:)
             hel(i,j,:) = e(6,:)
             identity(i,j) = 32
             iq = iq+1
             Lab_ex(i,j) ='bot'
          elseif (Lab_ex(i,j) == 'wwp') then 
             mom(i,j,:) = pext(3+iw,:)
             hel(i,j,:) = e(3+iw,:) 
             identity(i,j) = 2**(2+iw)
             iw = iw + 1
          else
             write(*,*) Lab_ex(i,j)
             stop 'pampl: undefined Lab_ex' 
          endif




       enddo
    enddo


    !--- calculate the tree-level amplitude -----------------------------
    !    for this, we need one particular configuration of external lines

    do i=1,Npoint
       sp(:,i) = hel(1,i,:)
       p(:,i)  = mom(1,i,:)
       fl1 = Lab_ex(1,:)
    enddo
    nq = 4
    nw = 2
    ng = Npoint - nq-nw

    call calc_ampl(4,4,ng,nq,nw,sp,p,fl1,vg,vs,identity(1,:),0)

    tree = psp1(vs,sp(:,1))         

    call oneloopdiv(mu,p,div2,div1, div3) ! divergence of the amplitude 
    div3  = czero

    !       the list of ``propagator momenta'' ; for the two cases

    momline = czero

    do i=1,Ncut
       momline(i,1,:)=czero

       do j=1,Npoint-1
          do j1=1,j
             momline(i,j+1,:)  = momline(i,j+1,:)+mom(i,j1,:)   
          enddo
       enddo

    enddo


    call  getallcuts(lncuts,npoint,Ac5,Ac4,Ac3,Ac2,Ac1)

    xN5 = 0
    xN4 = 0
    xN3 = 0
    xN2 = 0

    do i =1,lncuts(5)
       nodum = .true.
       do j =1,5

          if (Lab_in(Ncut,Ac5(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          xN5 = xN5+1
          Lc5(xN5,:) = Ac5(i,:)         
          F5(xN5,:) = Lab_in(Ncut, Ac5(i,:))
          Yc5(xN5,1) = 1
       endif
    enddo


    if (xN5 /= N5) then
       write(*,*) 'WARNING: Error occurred in pentuple cuts (pamplWW4q)'
       stop
    endif

    do i =1,lncuts(4)
       nodum = .true.
       do j =1,4

          if (Lab_in(Ncut,Ac4(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          xN4 = xN4+1
          Lc4(xN4,:) = Ac4(i,:)
          F4(xN4,:) = Lab_in(Ncut, Ac4(i,:))
          Yc4(xN4,1) = 1
       endif
    enddo

    if (xN4 /= N4) then
       write(*,*) 'WARNING: Error occurred in quadruple cuts (pamplWW4q)'
       stop
    endif

    do i =1,lncuts(3)
       nodum = .true.
       do j =1,3

          if (Lab_in(Ncut,Ac3(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          xN3 = xN3+1
          Lc3(xN3,:) = Ac3(i,:)
          F3(xN3,:) = Lab_in(Ncut, Ac3(i,:))
          Yc3(xN3,1) = 1
       endif
    enddo

    if (xN3 /= N3) then
       write(*,*) 'WARNING: Error occurred in triple cuts (pamplWW4q)'
       stop
    endif


    do i =1,lncuts(2)
       nodum = .true.
       do j =1,2

          if (Lab_in(Ncut,Ac2(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo
       if (nodum) then
          xN2 = xN2+1
          Lc2a(xN2,:) = Ac2(i,:)
          F2a(xN2,:) = Lab_in(Ncut, Ac2(i,:))
          Yc2a(xN2,1) = 1
       endif


    enddo

    ! Massless two-cuts must be removed.

    jglob = 0

    aN2 = lncuts(2)

    icount = 0


    do i=1,xN2


       xxx = 0

       i1 = Lc2a(i,1)
       i2 = Lc2a(i,2)



       if (i2.ne.i1 + 1) then 
          xxx = xxx + 1

          if (i1.eq.1.and.i2.eq.Npoint) then 
             xxx = xxx - 1
          endif

       else

          if (Lab_ex(Ncut,i1).ne.'top'.and.Lab_ex(Ncut,i1).ne.'bot'.and.&
               &Lab_ex(Ncut,i1).ne.'glu' .and. Lab_ex(Ncut,i1) .ne. 'str') then 
             xxx = xxx + 1
          endif
       endif



       if (xxx.eq.1) then 
          icount = icount + 1
          Lc2(icount,:) = Lc2a(i,:)
          F2(icount,:) = F2a(i,:)
          Yc2(icount,1) = 1

       endif


    enddo

    if (icount /= N2) then
       write(*,*) 'possible error in removing massless two cuts',icount, N2

    endif

    !----------- end disregarding massless 2-points

  end subroutine pamplWpWp


  !----------


  subroutine pampl_count_WpWp(lncuts,npoint,N5,N4,N3,xN2,N2,N1)
    ! This function is now only used to remove massless 2-cuts.
    ! The total number of 2-cuts, as well as the labels, are now INPUTS.
    ! Returned are the number of massive two cuts, and the number of 1-cuts (0).

    integer, intent(in):: lncuts(5), npoint ! input
    integer, intent(out):: N5, N4, N3,N1, xN2
    integer:: icount, i, i1,i2,j,N2
    integer::Acc5(lncuts(5),5), Acc4(lncuts(4),4)
    integer::Acc3(lncuts(3),3), Acc2(lncuts(2),2), Acc1(lncuts(1),1)
    character::Lab_ex(Ncut, Npoint)*3, Lab_in(Ncut,Npoint)*3
    logical:: nodum


    N5 = 0
    N4 = 0
    N3 = 0
    N2 = 0
    call getallcuts(lncuts, Npoint, Acc5, Acc4, Acc3, Acc2, Acc1)
    icount = 0


    if (case_b1) then
       Lab_ex(Ncut,:) = (/'top', 'wwp', 'bot', 'top', 'wwm', 'bot'/)
       Lab_in(Ncut,:) = (/'glu', 'top', 'bot', 'glu', 'top', 'bot'/) 
    elseif (case_b2) then
       Lab_ex(Ncut,:) = (/'top', 'wwp', 'bot','wwm', 'top', 'bot'/)
       Lab_in(Ncut,:) = (/'glu','top', 'bot', 'dum','dum', 'bot'/)
    endif



    do i =1,lncuts(5)
       nodum = .true.
       do j =1,5

          if (Lab_in(Ncut,Acc5(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          N5 = N5+1
       endif
    enddo


    do i =1,lncuts(4)
       nodum = .true.
       do j =1,4

          if (Lab_in(Ncut,Acc4(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          N4 = N4+1
       endif
    enddo


    do i =1,lncuts(3)
       nodum = .true.
       do j =1,3

          if (Lab_in(Ncut,Acc3(i,j)) .eq. 'dum') then
             nodum = .false.
          endif

       enddo

       if (nodum) then
          N3 = N3+1
       endif
    enddo

    xN2 = 0
    do i=1,lncuts(2)

       i1 = Acc2(i,1)
       i2 = Acc2(i,2)

       if ((Lab_in(Ncut,i1) == 'dum') .or. (Lab_in(Ncut,i2) == 'dum')) then

       else
          xN2 = xN2+1
          if (i2.ne.i1 + 1) then 

             icount = icount + 1

             if (i1.eq.1.and.i2.eq.Npoint) then 
                icount = icount - 1
             endif

          else

             if (Lab_ex(Ncut,i1).ne.'top'.and.Lab_ex(Ncut,i1).ne.'bot'&
                  & .and.Lab_ex(Ncut,i1).ne.'glu' .and. Lab_ex(Ncut, i1) .ne. 'str') then 

                icount = icount + 1

             endif

          endif
       endif

    enddo

    N2 = icount


    N1 = 0   ! for amplitudes with massless particles, no tadpoles

  end subroutine pampl_count_WpWp

  subroutine  oneloopdiv(mu,p,rdiv2,rdiv1,rdiv3)
    complex(qp), intent(in)  :: p(:,:)
    real(qp), intent(in)     :: mu
    complex(qp), intent(out) :: rdiv2, rdiv1,rdiv3
    ! ---------------------------------------------
    complex(qp) :: ap1(size(p,dim=1)), ap2(size(p,dim=1))
    complex(qp) :: s
    real(qp)    :: mu2
    integer       :: i, Nmax

    mu2 = mu**2
    rdiv3 = 0d0 

    if (case_a2) rdiv2 = -one
    if (case_b1 .or. case_b2) rdiv2 = -two
    rdiv1 = czero
    if (case_a1) then
       ap1 = p(:,4)
       ap2 = p(:,5)



       s = two*sc(ap1,ap2)
       if (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif

       ap1 = p(:,1)
       ap2 = p(:,6)
       s = two*sc(ap1,ap2)

       if (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif

       ! and finally
       rdiv1 = rdiv1  + two/three


    elseif (case_a2 .or. case_a4) then

       rdiv2 = -one
       ap1 = p(:,1)
       ap2 = p(:,6)
       s = two*sc(ap1,ap2)

       if (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif

       ! and finally
       rdiv1 = rdiv1  - three/two



    elseif (case_a3) then

       rdiv2=-one
       ap1 = p(:,2)
       ap2 = p(:,3)
       s = two*sc(ap1,ap2)

       if  (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif

       ! and finally
       rdiv1 = rdiv1-three/two

    elseif (case_b1) then
       ap1 = p(:,4)
       ap2 = p(:,3)

       s = two*sc(ap1,ap2)
       if (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif
       ap1 = p(:,1)
       ap2 = p(:,6)
       s = two*sc(ap1,ap2)
       
       if (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif
       ! and finally
       rdiv1 = rdiv1  + two/three

    elseif (case_b2) then
       rdiv2 = -one
       ap1 = p(:,1)
       ap2 = p(:,6)
       s = two*sc(ap1,ap2)
       if  (real(s) < zero) then 
          rdiv1 = rdiv1 -log(mu2/(-s))
       else
          rdiv1 = rdiv1 -log(mu2/s)-ci*pi 
       endif
       rdiv1  = rdiv1 -three/two
    endif

  end subroutine oneloopdiv


end module qpinitWpWp
