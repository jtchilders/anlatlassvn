!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECaux_functions_weyl.f90.

!     AUXILIARY FUNCTIONS
!     WEYL REPRESENTATION FOR gamma matrices

module qpaux_functions
  use types; use consts_qp 
  use define_ampl, only : use_pol_mass, gauge_check, MChel
  use qpglobal
  implicit none 
  private 

  public :: pol_qqgluons,pol_qqWgluons,pol_qqWqqgluons,pol_qqWWgluons
  public :: pol_qqZgluons
  
  interface pol_mless 
     module procedure pol_mless,pol_mless2 
  end interface

  interface pol_mass 
     module procedure pol_mass,pol_mass2 
  end interface

  public :: sc,spb2,spi2,psp1,pol_mless,pol_mass,pol_decay,ubar0,v0,pol_dk2mom
  public :: spi5,spb5

  public :: quark_flavour 
  public :: set_pol_ext 
contains
  
  
  subroutine set_pol_ext(pext)
    use qpauxiliary_functions  
    complex(qp), intent(in) :: pext(:,:)
    complex(qp) :: pWp(4),pWm(4),pZ(4),eWp(4),eWm(4)
    complex(qp) :: pp(size(pext,dim=1),size(pext,dim=2))
    complex(qp) :: pol_plus(4),pol_minus(4) 
    integer       :: ihel, i,j  
    
    MCFMhardscale =  dot(pext(1,:)+pext(2,:),pext(1,:)+pext(2,:))
    MCFMhardscale = sqrt(abs(MCFMhardscale))*2d0
    MCFMsqrthardscale = sqrt(MCFMhardscale)
    do i=1,size(pp,dim=1)
       do j=1,size(pp,dim=2)
          pp(i,j) = pext(i,j)/MCFMhardscale 
       enddo
    enddo
    do ihel=-1,1,2
       MCFMpol_ext_qq(1,:,ihel)  = v0(pp(1,:),ihel) 
       MCFMpol_ext_qq(2,:,ihel)  = ubar0(pp(2,:),ihel) 
       MCFMpol_ext_ww(1,:,ihel)  = pol_dk2mom(pp(3,:),pp(4,:),ihel,.true.) 
       MCFMpol_ext_ww(2,:,ihel)  = pol_dk2mom(pp(5,:),pp(6,:),ihel,.true.) 
       MCFMpol_ext_gg(1,:,ihel)  = pol_mless(pp(7,:),ihel,.true.) 
       MCFMpol_ext_gg(2,:,ihel)  = pol_mless(pp(8,:),ihel,.true.) 
    enddo
    pWp = pp(3,:)+pp(4,:)
    pWm = pp(5,:)+pp(6,:)
    pZ = pWp+pWm 
    eWp = MCFMpol_ext_ww(1,:,-1)
    eWm = MCFMpol_ext_ww(2,:,-1)
    ! Z->WW->4l 
    MCFMpol_ext_z(:) = ((pWp-pWm)*dot(eWp,eWm) &
         & + 2d0*dot(pWm,eWp)*eWm& 
         & - 2d0*dot(pWp,eWm)*eWp)/dot(pZ,pZ)
    
    if (MChel) then 
       do i = 1,2 
          pol_plus  = MCFMpol_ext_gg(i,:,1)
          pol_minus = MCFMpol_ext_gg(i,:,-1)
          MCFMpol_ext_gg(i,:,-1)  = (exp(-ci*twopi*MCphi(i))*pol_minus+&
               &exp(ci*twopi*MCphi(i))*pol_plus)
       enddo
    endif
    
    
  end subroutine set_pol_ext


  subroutine pol_qqWgluons(pext,pl,pa,hl,e,outgoing)
    integer, intent(in)           :: hl(:)
    complex(qp), intent(in)     :: pext(:,:),pl(:),pa(:)
    logical, intent(in), optional :: outgoing(:)
    complex(qp), intent(out)    :: e(:,:)
    integer :: i 

    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn

    if (use_pol_mass) then 
       e(3,:) = pol_mass(pext(3,:),mw,hl(3),outgoing(3))  !w
    else
       e(3,:) = pol_dk2mom(pl,pa,hl(3),outgoing(3))  !w
    endif
    
    if (MChel) then 
       do i=4,size(pext,dim=1)
          e(i,:) = (exp(-ci*twopi*MCphi(i-3))*pol_mless(pext(i,:),-1,outgoing(i))+&
               &exp(ci*twopi*MCphi(i-3))*pol_mless(pext(i,:),1,outgoing(i)))
       enddo
    else
       do i=4,size(pext,dim=1)
          e(i,:) = pol_mless(pext(i,:),hl(i),outgoing(i))     !gluons
       enddo
    endif
    if (gauge_check) e(3,:) = pext(3,:)

  end subroutine pol_qqWgluons


  subroutine pol_qqZgluons(pext,pl,pa,pl2,pa2,hl,e,outgoing)
    use qpauxiliary_functions 
    integer, intent(in)           :: hl(:)
    complex(qp), intent(in)     :: pext(:,:),pl(:),pa(:),pl2(:),pa2(:)
    logical, intent(in), optional :: outgoing(:)
    complex(qp), intent(out)    :: e(:,:)
    integer :: i 
    complex(qp) :: pWp(4),pWm(4),pZ(4),eWp(4),eWm(4)

    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn
    
    pWp = pl+pa
    pWm = pl2+pa2
    pZ = pWp+pWm 
    eWp = pol_dk2mom(pl,pa,-1,outgoing(3))
    eWm = pol_dk2mom(pl2,pa2,-1,outgoing(3))

    e(3,:) = ((pWp-pWm)*dot(eWp,eWm) + 2d0*dot(pWm,eWp)*eWm& ! Z->WW->4l 
         &- 2d0*dot(pWp,eWm)*eWp)/dot(pZ,pZ)

    if (MChel) then 
       do i=4,size(pext,dim=1)
          e(i,:) = (exp(-ci*twopi*MCphi(i-3))*pol_mless(pext(i,:),-1,outgoing(i))+&
               &exp(ci*twopi*MCphi(i-3))*pol_mless(pext(i,:),1,outgoing(i)))
       enddo
    else
       do i=4,size(pext,dim=1)
          e(i,:) = pol_mless(pext(i,:),hl(i),outgoing(i))     !gluons
       enddo
    endif

  end subroutine pol_qqZgluons

 subroutine pol_qqWWgluons(pext,pl,pa,hl,e,outgoing)
    integer, intent(in)           :: hl(:)
    complex(qp), intent(in)     :: pext(:,:),pl(:,:),pa(:,:)
    logical, intent(in), optional :: outgoing(:)
    complex(qp), intent(out)    :: e(:,:)
    integer :: i 
    complex(qp)	 :: factp, factm


    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn

    e(3,:) = pol_dk2mom(pl(1,:),pa(1,:),hl(3),outgoing(3))  !wm
    e(4,:) = pol_dk2mom(pl(2,:),pa(2,:),hl(4),outgoing(4))  !wp

    if (MChel) then 
       do i=5,size(pext,dim=1)
          e(i,:) = (exp(-ci*twopi*MCphi(i-4))*pol_mless(pext(i,:),-1,outgoing(i))+&
               &exp(ci*twopi*MCphi(i-4))*pol_mless(pext(i,:),1,outgoing(i)))
       enddo
    else
       do i=5,size(pext,dim=1)
          e(i,:) = pol_mless(pext(i,:),hl(i),outgoing(i))   !gluons
       enddo
    endif

    if (gauge_check) e(3,:) = pext(3,:)

  end subroutine pol_qqWWgluons


  subroutine pol_qqgluons(pext,hl,e,outgoing)
    integer, intent(in)           :: hl(:)
    complex(qp), intent(in)     :: pext(:,:)
    logical, intent(in), optional :: outgoing(:)
    complex(qp), intent(out)    :: e(:,:)
    integer :: i 
    logical, save :: first_time = .true. 

    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn

    if (MChel) then 
       do i=3,size(pext,dim=1)
          e(i,:) = (exp(-ci*twopi*MCphi(i-2))*pol_mless(pext(i,:),-1,outgoing(i))+&
               &exp(ci*twopi*MCphi(i-2))*pol_mless(pext(i,:),1,outgoing(i)))
       enddo
    else
       do i=3,size(pext,dim=1)
          e(i,:) = pol_mless(pext(i,:),hl(i),outgoing(i))     !gluons
       enddo
    endif
    if (gauge_check) then 
       e(3,:) = pext(3,:)
       if (first_time) then 
          write(*,*) 'doing gauge check: set e3=p3'
          first_time = .false.
       endif
    end if

  end subroutine pol_qqgluons


  subroutine pol_qqWqqgluons(pext,pl,pa,hl,e,vub,outgoing,swap)
    integer, intent(in)           :: hl(:)
    complex(qp), intent(in)     :: pext(:,:),pl(:),pa(:)
    logical, intent(in)           :: vub
    logical, intent(in), optional :: outgoing(:),swap 
    complex(qp), intent(out)    :: e(:,:)
    integer :: i 
    logical, save :: first_time = .true. 


    e(1,:) = v0(pext(1,:),hl(1))           !bu
    e(2,:) = ubar0(pext(2,:),hl(2))        !dn

    if (use_pol_mass) then 
       e(3,:) = pol_mass(pext(3,:),mw,hl(3),outgoing(3))  !w
    else
       e(3,:) = pol_dk2mom(pl,pa,hl(3),outgoing(3))  !w
    endif
    if (swap) then 
       if (vub) then 
          e(4,:) = v0(pext(4,:),hl(5))           !s
          e(5,:) = ubar0(pext(5,:),hl(4))        !bs
       else
          e(4,:) = ubar0(pext(4,:),hl(5))        !bs
          e(5,:) = v0(pext(5,:),hl(4))           !s
       endif
    else
       if (vub) then 
          e(4,:) = v0(pext(4,:),hl(4))           !s
          e(5,:) = ubar0(pext(5,:),hl(5))        !bs
       else
          e(4,:) = ubar0(pext(4,:),hl(4))        !bs
          e(5,:) = v0(pext(5,:),hl(5))           !s
       endif
    endif

    if (MChel) then 
       do i=6,size(pext,dim=1)
          e(i,:) = (exp(-ci*twopi*MCphi(i-5))*pol_mless(pext(i,:),-1,outgoing(i))+&
               &exp(ci*twopi*MCphi(i-5))*pol_mless(pext(i,:),1,outgoing(i)))
       enddo
    else
       do i=6,size(pext,dim=1)
          e(i,:) = pol_mless(pext(i,:),hl(i),outgoing(i))     !gluons
       enddo
    endif
    if (gauge_check) e(3,:) = pext(3,:)

  end subroutine pol_qqWqqgluons


  
  !     ubar spinor, massless
  
  function ubar0(p,i)
    complex(qp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------      
    complex(qp) :: ubar0(4)
    complex(qp) :: fc, fc2
    real(qp)    :: p0,px,py,pz
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),qp)
    px=real(p(2),qp)
    py=real(p(3),qp)
    pz=real(p(4),qp)
!^^^END

    fc2 = p0 + pz 
    fc=sqrt(fc2)

    if (.not.( (abs(px) == 0d0) .and. (abs(py) == 0d0) .and. &
         &((p0>0.and. pz<0) .or. (p0<0d0.and. pz>0d0)))) then 
       
       if (i.eq.1) then 
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then 
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range' 
       endif
       
    else
       if (i.eq.1) then 
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then 
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range' 
       endif
    endif
    
    
  end function ubar0
  
  
  
  ! -- v0  spinor, massless
  function v0(p,i)
    complex(qp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------      
    complex(qp) :: v0(4)
    complex(qp) :: fc2, fc
    real(qp)    :: p0,px,py,pz
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),qp)
    px=real(p(2),qp)
    py=real(p(3),qp)
    pz=real(p(4),qp)
!^^^END

    fc2 = p0 + pz 
    fc=sqrt(fc2)

!    if (abs(fc2).gt. tol) then 
    if (.not.( (abs(px) == 0d0) .and. (abs(py) == 0d0) .and. &
         &((p0>0.and. pz<0) .or. (p0<0d0.and. pz>0d0)))) then 
       
       if (i.eq.1) then 
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range' 
       endif

    else
       
       if (i.eq.1) then 
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range' 
       endif
       
    endif
    
  end function v0
  

  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    complex(qp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------      
    integer :: pol
    real(qp) :: p0,px,py,pz
    real(qp) :: pv,ct,st,cphi,sphi
    complex(qp) :: pol_mless(4)
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),qp)
    px=real(p(2),qp)
    py=real(p(3),qp)
    pz=real(p(4),qp)
!^^^END   


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0_qp-ct**2))
    
    if (st < tol) then 
       cphi=1.0_qp
       sphi=0.0_qp
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif
    

    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_qp) then  
       pol=i
    else
       pol=-i
    endif
    
    ! -- take complex conjugate for outgoing 
    if (present(outgoing)) then 
       if (outgoing) pol = -pol 
    endif
    
    pol_mless(1)=czero 
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless
  

  function pol_mless2(p,i,out)
    integer, intent(in)       :: i
    complex(qp), intent(in) :: p(4) 
    character(len=*), intent(in):: out
    complex(qp)             :: pol_mless2(4)
    ! -------------------------------------

    if (out == 'out') then 
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2



!--------massive vector boson polarization routine
  
  function pol_mass(p,m,i,outgoing)
    integer, intent(in)       :: i
    complex(qp), intent(in) :: p(4) 
    real(qp),  intent(in)   :: m
    logical, intent(in),optional :: outgoing
    complex(qp)             :: pol_mass(4)
    ! -------------------------------------
    real(qp) :: p0,px,py,pz, pv
    real(qp) :: ct,st,cphi,sphi
    integer :: pol
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),qp)
    px=real(p(2),qp)
    py=real(p(3),qp)
    pz=real(p(4),qp)
!^^^END   
   
    pv= sqrt(abs(p0**2 - m**2))
    ct= pz/pv
    st= sqrt(abs(one-ct**2))
    
    if (st < tol) then 
       cphi=one
       sphi=zero
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif


    ! i=0 is longitudinal polarization          

    ! -- distinguish between positive and negative energies
    if ( p0 > zero) then  
       pol=i
    else
       pol=-i
    endif
    
    ! -- take complex conjugate for outgoing 
    if (present(outgoing)) then 
       if (outgoing) pol = -pol 
    endif

    if(pol == -1.or.pol == 1) then     
       pol_mass(1)=czero
       pol_mass(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
       pol_mass(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
       pol_mass(4)=-st/sqrt2
    elseif (pol == 0) then 
       pol_mass(1)= pv/m
       pol_mass(2)= p0/m/pv*px
       pol_mass(3)= p0/m/pv*py
       pol_mass(4)= p0/m/pv*pz
       else
          stop 'pol_mass: pol out of range'
       endif

  end function pol_mass

  function pol_mass2(p,m,i,out)
    integer, intent(in)       :: i
    complex(qp), intent(in) :: p(4) 
    real(qp),  intent(in)   :: m
    character(len=*), intent(in):: out
    complex(qp)             :: pol_mass2(4)
    ! -------------------------------------

    if (out == 'out') then 
       pol_mass2 = pol_mass(p,m,i,outgoing=.true.)
    else
       pol_mass2 = pol_mass(p,m,i,outgoing=.false.)
    endif
  end function pol_mass2


  !--------massive vector boson with decay polarization
  ! send by Keith with email on 12Set08
  function pol_dk2mom(plepton,antilepton,i,outgoing)
    integer, intent(in) :: i
    integer :: j
    complex(qp), intent(in) :: plepton(:),antilepton(:)
    logical, intent(in),optional :: outgoing
    complex(qp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq
    
    
    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2
    
    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)
    
    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))
    
    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo

  end function pol_dk2mom
  

  

!--------massive vector boson with decay polarization

  function pol_decay(p,m,i)
    complex(qp), intent(in) :: p(4) 
    real(qp),  intent(in)   :: m
    integer, intent(in)       :: i
    ! -------------------------------------
    complex(qp) :: pol_decay(4),zl(4),za(4),Ub(4),V(4)
    real(qp) :: p0,px,py,pz,q(4),lcm(4),l(4),a(4)
    real(qp) :: theta

!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=real(p(1),qp)
    px=real(p(2),qp)
    py=real(p(3),qp)
    pz=real(p(4),qp)
!^^^END   
    

    ! -- define momentum of vector boson in RKE notation
      q(4)=p0
      q(1)=px
      q(2)=py
      q(3)=pz
      
      ! choose  momentum of lepton in rest frame of vector boson in RKE notation
      ! with arbitrary angle theta
      theta=0.3324455_qp
    
      lcm(4)=0.5_qp*m
      lcm(1)=0.0_qp
      lcm(2)=lcm(4)*sin(theta)
      lcm(3)=lcm(4)*cos(theta)

      ! -- Calculate l form lcm (in CM frame of q)
      call boost(m,q,lcm,l)
      
      a(:)=q(:)-l(:)
      
      ! -- Complexify in Kirill's notation
      zl(1)=cone*l(4)
      zl(2)=cone*l(1)
      zl(3)=cone*l(2)
      zl(4)=cone*l(3)
      za(1)=cone*a(4)
      za(2)=cone*a(1)
      za(3)=cone*a(2)
      za(4)=cone*a(3)

      Ub(:)=ubar0(zl,i)
      V(:)=v0(za,-i)


      !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
      pol_decay(1)=Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3)
      pol_decay(2)=-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3)
      pol_decay(3)=ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
      pol_decay(4)=Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3)
     
    end function pol_decay


    ! take momenta p_in in frame in which particle one is at rest with mass 
    ! "mass" and convert to frame in which particle one has fourvector p1
    subroutine boost(mass,p1,p_in,p_out)
      real(qp), intent(in)  ::  mass,p1(4),p_in(4)
      real(qp), intent(out) ::  p_out(4)
      ! ------------------------------------       
      real(qp) ::  gam,beta(1:3),bdotp
      integer :: j,k

      gam=p1(4)/mass
      bdotp=zero 
      do j=1,3
         beta(j)=-p1(j)/p1(4)
         bdotp=bdotp+p_in(j)*beta(j)
      enddo
      p_out(4)=gam*(p_in(4)-bdotp)
      do k=1,3
         p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(4))
      enddo
    end subroutine boost
    
    
    
    function sc(p1,p2)
      complex(qp), intent(in) :: p1(:)
      complex(qp), intent(in) :: p2(:)
      complex(qp)             :: sc
      integer :: sizemin 
      
      sizemin=min(size(p1),size(p2))

      sc = p1(1)*p2(1)
      sc = sc - sum(p1(2:sizemin)*p2(2:sizemin))
      
    end function sc
      

    function spb2(sp,v) 
      complex(qp), intent(in) :: sp(:),v(:)
      complex(qp) :: spb2(size(sp))
      complex(qp) :: x0(4,4),xx(4,4),xy(4,4)
      complex(qp) :: xz(4,4),x5(4,4)
      complex(qp) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,Dv,Ds,imax

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spb2:Dv not allowed'
      endif

      imax = Ds/4

      do i=1,imax
         i1= 1+4*(i-1)
         i2=i1+3

         y1=sp(i1)
         y2=sp(i1+1)
         y3=sp(i1+2)
         y4=sp(i1+3)

         x0(1,i)=y3
         x0(2,i)=y4
         x0(3,i)=y1
         x0(4,i)=y2

         xx(1,i) = y4
         xx(2,i) = y3
         xx(3,i) = -y2
         xx(4,i) = -y1

         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3           
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=y3
         xz(2,i)=-y4
         xz(3,i)=-y1
         xz(4,i)=y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4
      enddo

      if (Dv.eq.4) then 

         do i=1,4
            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then 
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))

         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
                 &          -v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)
            
            i1 = i+4
            spb2(i1)= v(1)*x0(i,2)-v(2)*xx(i,2) &
                 &            -v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)
         enddo
      elseif (Dv.eq.8) then 
         bp=(v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))
         
         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &          -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &         +bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4
 
            spb2(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) &
            &            -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &            -bp*x5(i,1)+cm*x5(i,4)

            i2 = i1+4

            spb2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3) &
            &             -v(3)*xy(i,3)-v(4)*xz(i,3) & 
            &             +bm*x5(i,4)+cp*x5(i,1)

            i3=i2+4

            spb2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4) &
            &             -v(3)*xy(i,4)-v(4)*xz(i,4) &
            &             -bp*x5(i,3)-cp*x5(i,2)

         enddo
      else
         stop 'spb2: Dv out of bound' 
      endif

    end function spb2



    function spi2(v,sp)
      complex(qp), intent(in) :: sp(:),v(:)
      complex(qp) :: spi2(size(sp))
      complex(qp) :: x0(4,4),xx(4,4),xy(4,4)
      complex(qp) :: xz(4,4),x5(4,4)
      complex(qp) ::  y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,imax,Dv,Ds

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spi2:Dv not allowed'
      endif

      imax = Ds/4

      do i=1,imax
         i1= 1+4*(i-1)
         i2=i1+3

         y1=sp(i1)
         y2=sp(i1+1)
         y3=sp(i1+2)
         y4=sp(i1+3)

         x0(1,i)=y3
         x0(2,i)=y4
         x0(3,i)=y1
         x0(4,i)=y2


         xx(1,i) = -y4
         xx(2,i) = -y3
         xx(3,i) = y2
         xx(4,i) = y1


         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3           
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=-y3
         xz(2,i)=y4
         xz(3,i)=y1
         xz(4,i)=-y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4

      enddo

      if(Dv.eq.4) then 

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then 
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))


         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &           -bp*x5(i,2)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2) &
            &             -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &             +bm*x5(i,1)

         enddo

      elseif (Dv.eq.8) then 
         
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)&
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)&
            &           -bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)&
            &             -v(3)*xy(i,2)-v(4)*xz(i,2)&
            &             +bm*x5(i,1)-cp*x5(i,4)

            i2=i1+4

            spi2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)&
            &           -v(3)*xy(i,3)-v(4)*xz(i,3)&
            &          -bp*x5(i,4)-cm*x5(i,1)

            i3=i2+4

            spi2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)&
            &           -v(3)*xy(i,4)-v(4)*xz(i,4)&
            &           +bm*x5(i,3)+cm*x5(i,2)


         enddo

      else
         stop 'spi2: Dv out of bounds' 
      end if

    end function spi2


    function  psp1(sp1,sp2)
      complex(qp), intent(in) :: sp1(:)
      complex(qp), intent(in) :: sp2(:)
      complex(qp) :: psp1
      
      psp1 = sum(sp1(1:)*sp2(1:))
      
    end function psp1



    ! -- multiplication of spinor with gamma_5 on the left  
    function spi5(sp)
      complex(qp), intent(in) :: sp(:)
      complex(qp) :: spi5(size(sp))
      integer :: i,j,imax, Ds  

      Ds = size(sp)
      imax = Ds/4
      do i=1,imax
         spi5(4*(i-1)+1) = sp(4*(i-1)+1) 
         spi5(4*(i-1)+2) = sp(4*(i-1)+2) 
         spi5(4*(i-1)+3) = -sp(4*(i-1)+3) 
         spi5(4*(i-1)+4) = -sp(4*(i-1)+4) 
      enddo

    end function spi5


    ! -- multiplication of bspinor with gamma_5 on the right  
    function spb5(sp)
      complex(qp), intent(in) :: sp(:)
      complex(qp) :: spb5(size(sp))
      
      spb5 = spi5(sp) 

    end function spb5
    
 logical function quark_flavour(fl)
    character(len=3), intent(in) :: fl 

    if ((fl == 'top') .or. (fl == 'bot') .or. (fl == 'str') .or. &
         &(fl == 'chr') .or. (fl == 'ch1') .or. (fl == 'ch2') .or. &
         &(fl == 'dwn').or. (fl == 'to1') .or. (fl == 'to2')) then 
       quark_flavour = .true. 
    elseif (fl == 'glu' .or. fl == 'wwp' .or. fl == 'wwm' .or. fl == 'www') then 
       quark_flavour = .false. 
    else
       write(*,*) 'fl',fl 
       stop 'qpampl:: quark_flavour: unrecognized flavour' 
    endif 
    
  end function quark_flavour


  
end module qpaux_functions
