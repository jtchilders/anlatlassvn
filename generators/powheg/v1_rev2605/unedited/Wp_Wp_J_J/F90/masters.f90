module masters
  use types; use consts_dp 
  use dpaux_functions; use dpvvn
  use match1; use define_ampl
  use dpglobal 
  use dpstability 
  implicit none 
  private 

  public :: MasterSubs
  
contains


  subroutine MasterSubs(mu,Lc5,F5,Yc5,Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2,res)
    real(dp), intent(in)     :: mu
    integer, intent(in)      :: Lc5(:,:),Lc4(:,:),Lc3(:,:),Lc2(:,:)
    character, intent(in)    :: F5(:,:)*3,F4(:,:)*3
    character, intent(in)    :: F3(:,:)*3, F2(:,:)*3
    integer, intent(in)      :: Yc5(:,:),Yc4(:,:),Yc3(:,:)
    integer, intent(in)      :: Yc2(:,:)
    complex(dp), intent(out) :: res(-2:1)
    ! -------------------------------------------------------------       
    integer :: i,j,j1,xe,ia,ib,ij
    complex(dp) :: d1,d5,c1,c8
    complex(dp) :: b1,b10,b2,b5
    complex(dp) :: k1d(16)
    complex(dp) :: k1(16),k2(16),pr1,pr2
    complex(dp) :: p1(4),p2(4),p3(4),p4(4)
    complex(dp) :: p12(4),p23(4)
    complex(dp) :: cs11,cs22,cs33,cs44,cs12,cs23
    complex(dp) ::  qlI4,qlI3,qlI2
    complex(dp) :: d(5,1)
    complex(dp) :: dres_tri(-2:1),dres_bub(-2:1)
    complex(dp) :: res_tri(-2:1),res_bub(-2:1),res_box(-2:1)
    real(dp) :: s11,s22,s33,s44,s12,s23
    real(dp) ::  m12,m22,m32,m42
    real(dp) :: argk(5,6),argm(5,4)
    real(dp) :: mu2
    integer       :: N5,N4,N3,N2
    logical, save :: first_time = .true. 
    real(dp) :: est_err, abs_res

    N5 = size(Lc5,dim=1)
    N4 = size(Lc4,dim=1)
    N3 = size(Lc3,dim=1)
    N2 = size(Lc2,dim=1)

    mu2=mu**2
    
    if (first_time) then 
       call qlinit
       first_time = .false. 
    endif
    res(-2)=czero
    res(-1)=czero
    res(0)=czero
    res(1)=czero

    dres_tri(-2)=czero
    dres_tri(-1)=czero
    dres_tri(0)=czero
    dres_tri(1)=czero

    dres_bub(-2)=czero
    dres_bub(-1)=czero
    dres_bub(0)=czero
    dres_bub(1)=czero

    res_tri(-2)=czero
    res_tri(-1)=czero
    res_tri(0)=czero
    res_tri(1)=czero

    res_box(-2)=czero
    res_box(-1)=czero
    res_box(0)=czero
    res_box(1)=czero

    res_bub(-2)=czero
    res_bub(-1)=czero
    res_bub(0)=czero
    res_bub(1)=czero

    ! 5cuts are dropped because of the parameterization of the 5-cut
    ! there is no contribution to the rational part as well


    !     4 cut

    do i=1,N4

       d1=coeff4(i,1)
       d5=coeff4(i,5)
       ij = Yc4(i,1)  ! GPS


       !------ finding the momenta

       do j=1,4
          p1(j)=czero
          p2(j)=czero
          p3(j)=czero
          p4(j)=czero
       enddo

       ia = Lc4(i,1)
       ib = Lc4(i,2)-1

       do j=ia,ib 
          do j1=1,4
             p1(j1) = p1(j1) + mom(ij,j,j1)
          enddo
       enddo


       ia = Lc4(i,2)
       ib = Lc4(i,3)-1

       do j=ia,ib 
          do j1=1,4
             p2(j1) = p2(j1) + mom(ij,j,j1)
          enddo
       enddo



       ia = Lc4(i,3)
       ib = Lc4(i,4)-1



       do j=ia,ib 
          do j1=1,4
             p3(j1) = p3(j1) + mom(ij,j,j1)
          enddo
       enddo


       do j1=1,4
          p4(j1) = -p1(j1)-p2(j1)-p3(j1)
       enddo


       do j=1,4
          p12(j)=p1(j)+p2(j)
          p23(j)=p2(j)+p3(j)
       enddo

       cs11 =  sc(p1,p1)
       cs22 =  sc(p2,p2)
       cs33 =  sc(p3,p3)
       cs44 =  sc(p4,p4)

       cs12 = sc(p12,p12)
       cs23 = sc(p23,p23)

       s11=real(cs11,dp)



       s22=real(cs22,dp)
       s33=real(cs33,dp)
       s44=real(cs44,dp)

       s12=real(cs12,dp)
       s23=real(cs23,dp)

       m12=mass4(i,1)**2
       m22=mass4(i,2)**2
       m32=mass4(i,3)**2
       m42=mass4(i,4)**2


       if (abs(s11-mw**2).lt.propcut) then 
          s11=mw**2
       endif

       if (abs(s22-mw**2).lt.propcut) then 
          s22=mw**2
       endif

       if (abs(s33-mw**2).lt.propcut) then 
          s33=mw**2
       endif

       if (abs(s44-mw**2).lt.propcut) then 
          s44=mw**2
       endif


       if (abs(s11).lt.propcut) then 
          s11=zero
       endif

       if (abs(s22).lt.propcut) then 
          s22=zero
       endif

       if (abs(s33).lt.propcut) then 
          s33=zero
       endif

       if (abs(s44).lt.propcut) then 
          s44=zero
       endif


       do xe=-2,0
          res(xe) = res(xe) + &
               &d1*qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,xe)
          res_box(xe) = res_box(xe) + &
               &d1*qlI4(s11,s22,s33,s44,s12,s23,m12,m22,m32,m42,mu2,xe)
       enddo

       res(1) = res(1)-d5/6._dp
       res_box(1) = res_box(1)-d5/6._dp

       !-------------------- for the Lc4 loop
    enddo


    !------------------------3 cut



    do i=1,N3

       c1=coeff3(i,1)
       c8=coeff3(i,8)

       ij = Yc3(i,1) ! GPS

       !------finding the momenta 

       do j=1,4
          p1(j)=czero
          p2(j)=czero
          p3(j)=czero
       enddo

       ia = Lc3(i,1)
       ib = Lc3(i,2)-1

       do j=ia,ib 
          do j1=1,4
             p1(j1) = p1(j1) + mom(ij,j,j1)
          enddo
       enddo


       ia = Lc3(i,2)
       ib = Lc3(i,3)-1

       do j=ia,ib 
          do j1=1,4
             p2(j1) = p2(j1) + mom(ij,j,j1)
          enddo
       enddo


       do j1=1,4
          p3(j1) = -p1(j1)-p2(j1)
       enddo

       cs11 = sc(p1,p1)
       cs22 = sc(p2,p2)
       cs33 = sc(p3,p3)



       s11=real(cs11,dp)
       s22=real(cs22,dp)
       s33=real(cs33,dp)

       m12=mass3(i,1)**2
       m22=mass3(i,2)**2
       m32=mass3(i,3)**2


       if (abs(s11-mw**2).lt.propcut) then 
          s11=mw**2
       endif

       if (abs(s22-mw**2).lt.propcut) then 
          s22=mw**2
       endif

       if (abs(s33-mw**2).lt.propcut) then 
          s33=mw**2
       endif

       if (abs(s11).lt.propcut) then 
          s11=zero
       endif

       if (abs(s22).lt.propcut) then 
          s22=zero
       endif

       if (abs(s33).lt.propcut) then 
          s33=zero
       endif


       do xe=-2,0
          res(xe) = res(xe) + c1*qlI3(s11,s22,s33,m12,m22,m32,mu2,xe)
       enddo

       res(1) = res(1) -half*c8


       do xe=-2,0
          dres_tri(xe) = dres_tri(xe) + dcoeff3(i,1)*qlI3(s11,s22,s33,m12,m22,m32,mu2,xe)
          res_tri(xe) = res_tri(xe) +coeff3(i,1)*qlI3(s11,s22,s33,m12,m22,m32,mu2,xe)
       enddo
       dres_tri(1) = dres_tri(1) -half*dcoeff3(i,8)
       res_tri(1) = res_tri(1) -half*c8


       !cccccccccccccccccccccccc for the Lc3 loop
    enddo

    !     2 cut


    do i=1,N2

       ij = Yc2(i,1)


       !-----buble off the light cone
       if (tagdcut(i,1).eq.666) then 

          b1=coeff2(i,1)
          b10=coeff2(i,10)


          !      finding the momenta -ccccccccccccccccccccc!      

          do j=1,4
             p1(j)=czero
             p2(j)=czero
          enddo

          ia = Lc2(i,1)
          ib = Lc2(i,2)-1

          do j=ia,ib 
             do j1=1,4
                p1(j1) = p1(j1) + mom(ij,j,j1)
             enddo
          enddo

          do j1=1,4
             p2(j1) = -p1(j1)
          enddo

          cs11 =  sc(p1,p1)

          s11=real(cs11,dp)

          m12=mass2(i,1)**2
          m22=mass2(i,2)**2


          do xe=-2,0
             res(xe) = res(xe) + b1*qlI2(s11,m12,m22,mu2,xe)
          enddo

          do xe=-2,0
             dres_bub(xe) = dres_bub(xe)+dcoeff2(i,1)*qlI2(s11,m12,m22,mu2,xe)
             res_bub(xe) = res_bub(xe)+b1*qlI2(s11,m12,m22,mu2,xe)
          enddo

       elseif (tagdcut(i,1).eq.999) then 
          !-----bubble on a light-cone
          
          b1=coeff2(i,1)
          b2=coeff2(i,2)
          b5=coeff2(i,5)
          b10=coeff2(i,10)

          do j=1,4
             k1(j)=propv2(i,j)
             k2(j)=propv2(i,4+j)
             k1d(j)=refvect2(i,j)   
          enddo


          !      finding the momenta -cccccccccccccccccccccc      

          do j=1,4
             p1(j)=czero
             p2(j)=czero
          enddo

          ia = Lc2(i,1)
          ib = Lc2(i,2)-1

          do j=ia,ib 
             do j1=1,4
                p1(j1) = p1(j1) + mom(ij,j,j1)
             enddo
          enddo

          do j1=1,4
             p2(j1) = -p1(j1)
          enddo

          cs11 =  sc(p1,p1)

          s11=real(cs11,dp)

          m12=mass2(i,1)**2
          m22=mass2(i,2)**2


          pr1 =  sc(k1,k1d)
          pr2 =  sc(k2,k1d)
          pr1=-pr1
          pr2=-pr2

          do xe=-2,0
             res(xe) = res(xe) + b1*qlI2(s11,m12,m22,mu2,xe) &
                  &+ b2*vBub(pr1,pr2,m12,m22,mu2,xe) &
                  &+ b5*tBub(pr1,pr2,m12,m22,mu2,xe)
          enddo

          do xe=-2,0
             dres_bub(xe) = dres_bub(xe) + dcoeff2(i,1)*qlI2(s11,m12,m22,mu2,xe) &
                  &+ dcoeff2(i,2)*vBub(pr1,pr2,m12,m22,mu2,xe) &
                  &+ dcoeff2(i,5)*tBub(pr1,pr2,m12,m22,mu2,xe)
             res_bub(xe) = res_bub(xe) + b1*qlI2(s11,m12,m22,mu2,xe) &
                  &+ b2*vBub(pr1,pr2,m12,m22,mu2,xe) &
                  &+ b5*tBub(pr1,pr2,m12,m22,mu2,xe)
          enddo


       endif

       !--------contribution to the rational part

       res(1) = res(1) +b10*(-half)*(m12+m22-s11/three)
       dres_bub(1) = dres_bub(1)+dcoeff2(i,10)*(-half)*(m12+m22-s11/three)
       res_bub(1) = res_bub(1)+b10*(-half)*(m12+m22-s11/three)

       !cccccccccccccccccccccccc for the Lc2 loop
    enddo
    
    abs_res  = abs(res(0)+res(1))
    est_err = max(abs(dres_tri(0)+dres_tri(1)),abs(dres_bub(0)+dres_bub(1)))

    abs_rel_err = est_err/abs_res 

  end subroutine MasterSubs


  !-----auxiliary master integrals (vector & tensor bubbles on a light-cone)
  !---- currently, for equal masses only

  function vBub(pr1,pr2,m12,m22,mu2,xe)
    integer, intent(in)     :: xe
    complex(dp), intent(in) ::  pr1,pr2
    real(dp), intent(in)    ::  m12,m22,mu2
    complex(dp) :: vBub

    if (xe == -2) then 
       vBub=czero
    elseif (xe == -1) then 
       vBub =-half*pr2-half*pr1
    elseif (xe == 0) then 
       vBub =half*(pr2+pr1)*log(m12/mu2)
    else
       stop 'vBub: xe out of range' 
    endif

  end function Vbub


  function tBub(pr1,pr2,m12,m22,mu2,xe)
    integer, intent(in)     ::  xe
    complex(dp), intent(in) ::  pr1,pr2
    real(dp), intent(in)    ::  m12,m22,mu2
    complex(dp) :: tBub

    if (xe == -2) then 
       tBub=czero
    elseif (xe == -1) then 
       tBub =-one/three*(pr1-pr2)**2-pr2*(pr1-pr2)-pr2**2
    elseif (xe == 0) then 
       tBub = (pr2**2+ pr2*(pr1-pr2)+one/three*(pr1-pr2)**2)*log(m12/mu2)
    else
       stop 'tBub: xe out of range' 
    endif

  end function tBub


end module masters
