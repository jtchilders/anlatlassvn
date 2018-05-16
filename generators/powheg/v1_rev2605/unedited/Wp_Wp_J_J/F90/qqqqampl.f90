module qqqqampl
  use types; use consts_dp;
  use dpaux_functions
  use dpmemory
  use dprecurrencebitsfour
  use define_ampl 
  implicit none
  private
  public :: getamplqqqq, getamplqqqq_stu, mtotqqqq, getamplqqqq_my

contains


  !----------------------------------------------------------!
  !   0 --> qb(1) q(2) ve(3) e+(4) vm(5) mu+(6) qb(7) q(8)   !
  !----------------------------------------------------------!

  !------------------------------!
  !getamplqqqqg returns          !
  !------------------------------!
  !     q(2)|            |qb(7)  !
  !    (u1) |            | (v2)  !
  !         ^            v       !
  !         |~~~~~~~~~~~~|       !
  !         |            |       !
  !         ^            v       !
  !         |            |       !
  !    qb(1)|            |q(8)   !
  !     (v1)    Amp(1)     (u2)  !
  !------------------------------!
  

  subroutine getamplqqqq_my(q,hel,plep,Amp)
    real(dp), intent(in) :: q(4,4),plep(4,4)
    integer, intent(in) :: hel(4)
    real(dp)  :: pin(8,4)
    complex(dp), intent(out)  :: Amp
    !--------------------------------------
    complex(dp) :: p(8,4)
    complex(dp) :: u1(4), v1(4), u2(4), v2(4)
    complex(dp) :: eW1(4), eW2(4), pW1(4), pW2(4)
    complex(dp) :: tmp(4)
    complex(dp) :: edum(4,0),kdum(4,0)
    complex(dp) :: eW(4,2),pW(4,2)
    complex(dp) :: pq(4,3),sp(4,3)
    character(len=3),save :: fl0,fl1(3) 
    integer,save :: giarray(0), qiarray(3), WWid(2), pol_int  
    integer :: iswap(8),i
    logical, save :: firsttime = .true.

    Amp = czero

    if ((hel(1).ne.1).or.(hel(2).ne.-1).or.(hel(3).ne.1).or.(hel(4).ne.-1)) return 
      
    giarray(:) = 0 
    qiarray(:) = (/1,2,4/)
    WWid = (/8,16/) 
    pol_int = 0 
    WWqqqq = .true. 
    cashing = .false. 
    fl1 = (/'top','bot','top'/)
    fl0 = 'bot'
    firsttime = .false.
    case_b1 = .true. 
    case_b2 = .false. 

    ! complexify momenta and swap spacetime indices (probably...)


    p(1,:) = dcmplx(q(:,1),0d0) 
    p(2,:) = dcmplx(q(:,2),0d0) 
    p(3,:) = dcmplx(plep(:,1),0d0) 
    p(4,:) = dcmplx(plep(:,2),0d0) 
    p(5,:) = dcmplx(plep(:,3),0d0) 
    p(6,:) = dcmplx(plep(:,4),0d0) 
    p(7,:) = dcmplx(q(:,3),0d0) 
    p(8,:) = dcmplx(q(:,4),0d0)

    pq(:,1) = p(2,:) 
    pq(:,2) = p(7,:) 
    pq(:,3) = p(8,:) 
    sp(:,1) = ubar0(p(2,:),-1) ! u1
    sp(:,2) = v0(p(7,:),1)     ! v2
    sp(:,3) = ubar0(p(8,:),-1) ! u2


    ! set Ws polarisation
    eW1 = pol_dk2mom(p(3,:),p(4,:),-1,.true.)
    eW2 = pol_dk2mom(p(5,:),p(6,:),-1,.true.)
    pW1 = p(3,:)+p(4,:)
    pW2 = p(5,:)+p(6,:)

    eW(:,1) = eW1
    eW(:,2) = eW2 
    pW(:,1) = pW1
    pW(:,2) = pW2 


    tmp = fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 
    
    !-- swap Ws
    eW(:,1) = eW2
    eW(:,2) = eW1 
    pW(:,1) = pW2
    pW(:,2) = pW1 

    tmp = tmp+fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp = psp1(tmp,v0(p(1,:),1))  ! psp1(tmp,v1)

  end subroutine getamplqqqq_my





  subroutine getamplqqqq(pin,i1,i2,i3,i4,Amp)
    real(dp),intent(in)       :: pin(8,4)
    integer, intent(in)       :: i1,i2,i3,i4
    complex(dp), intent(out)  :: Amp
    !--------------------------------------
    complex(dp) :: p(8,4)
    complex(dp) :: u1(4), v1(4), u2(4), v2(4)
    complex(dp) :: eW1(4), eW2(4), pW1(4), pW2(4)
    complex(dp) :: tmp(4)
    complex(dp) :: edum(4,0),kdum(4,0)
    complex(dp) :: eW(4,2),pW(4,2)
    complex(dp) :: pq(4,3),sp(4,3)
    character(len=3), save :: fl0,fl1(3) 
    integer, save :: giarray(0), qiarray(3), WWid(2), pol_int  
    integer :: iswap(8),i
    logical, save :: firsttime = .true.

    Amp = czero

    giarray(:) = 0 
    qiarray(:) = (/1,2,4/)
    WWid = (/8,16/) 
    pol_int = 0 
    fl1 = (/'top','bot','top'/)
    fl0 = 'bot'
    firsttime = .false.
    WWqqqq = .true. 
    ! -- virtual calculation might changes, make sure they are set here 
    case_b1 = .true. 
    case_b2 = .false.
       
    ! complexify momenta and swap spacetime indices (probably...)

    iswap = (/i1,i2,3,4,5,6,i3,i4/)
    do i=1,8
       p(i,1) = dcmplx(pin(iswap(i),4),0d0)
       p(i,2) = dcmplx(pin(iswap(i),1),0d0)
       p(i,3) = dcmplx(pin(iswap(i),2),0d0)
       p(i,4) = dcmplx(pin(iswap(i),3),0d0)
    enddo

    
    ! set Ws polarisation
    eW1 = pol_dk2mom(p(3,:),p(4,:),-1,.true.)
    eW2 = pol_dk2mom(p(5,:),p(6,:),-1,.true.)
    pW1 = p(3,:)+p(4,:)
    pW2 = p(5,:)+p(6,:)

    eW(:,1) = eW1
    eW(:,2) = eW2 
    pW(:,1) = pW1
    pW(:,2) = pW2 

    
    pq(:,1) = p(2,:) 
    pq(:,2) = p(7,:) 
    pq(:,3) = p(8,:) 
    sp(:,1) = ubar0(p(2,:),-1) ! u1
    sp(:,2) = v0(p(7,:),1)     ! v2
    sp(:,3) = ubar0(p(8,:),-1) ! u2

    tmp = fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    !-- swap Ws
    eW(:,1) = eW2
    eW(:,2) = eW1 
    pW(:,1) = pW2
    pW(:,2) = pW1 

    tmp = tmp+fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp = psp1(tmp,v0(p(1,:),1))  ! psp1(tmp,v1)

  end subroutine getamplqqqq


  !----------------------------------------------------------!
  !  pp(1) pp(2) --> ve(3) e+(4) vm(5) mu+(6) pp(7) pp(8)    !
  !----------------------------------------------------------!

  !----------------------------------------------------------!
  !  returns Amp for a given channel: qqb, qbq, qqq, qbb     !
  !     Amp(1:2) = 's diagram','t diagram' or 'u diagram'    !
  !----------------------------------------------------------!
  !   below now all momenta are outgoing                     !
  !----------------------------------------------------------!
 
  !------------------------------!------------------------------!
  !getamplqqqqg returns                                         !
  !------------------------------!------------------------------!
  !qqb channel                                                  !
  !                                                             !
  !     q(2)|            |qb(7)  !     q(8)|            |qb(7)  !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |~~~~~~~~~~~~|       !         |~~~~~~~~~~~~|       !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |            |       !         |            |       !
  !    qb(1)|            |q(8)   !    qb(1)|            |q(2)   !
  !             Amp(1)           !             Amp(2)           !
  !------------------------------!------------------------------!
  !qbq channel                                                  !
  !                                                             !
  !     q(1)|            |qb(7)  !     q(8)|            |qb(7)  !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |~~~~~~~~~~~~|       !         |~~~~~~~~~~~~|       !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |            |       !         |            |       !
  !    qb(2)|            |q(8)   !    qb(2)|            |q(1)   !
  !             Amp(1)           !             Amp(2)           !
  !------------------------------!------------------------------!
  !qqq channel                                                  !
  !                                                             !
  !     q(7)|            |qb(2)  !     q(8)|            |qb(2)  !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |~~~~~~~~~~~~|       !         |~~~~~~~~~~~~|       !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |            |       !         |            |       !
  !    qb(1)|            |q(8)   !    qb(1)|            |q(7)   !
  !             Amp(1)           !             Amp(2)           !
  !------------------------------!------------------------------!
  !qbb channel                                                  !
  !                                                             !
  !     q(1)|            |qb(8)  !     q(1)|            |qb(7)  !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |~~~~~~~~~~~~|       !         |~~~~~~~~~~~~|       !
  !         |            |       !         |            |       !
  !         ^            v       !         ^            v       !
  !         |            |       !         |            |       !
  !    qb(7)|            |q(2)   !    qb(8)|            |q(2)   !
  !             Amp(1)           !             Amp(2)           !
  !------------------------------!------------------------------!

  subroutine getamplqqqq_stu(pin,chn,Amp)
    real(dp),intent(in)       :: pin(8,4)
    character, intent(in)     :: chn*3
    complex(dp), intent(out)  :: Amp(2)

    Amp = czero

    ! --- now calculate the Amps
    if (chn .eq. 'qqb') then

       call getamplqqqq(pin,1,2,7,8,Amp(1))
       call getamplqqqq(pin,1,8,7,2,Amp(2))
    elseif (chn .eq. 'qbq') then

       call getamplqqqq(pin,2,1,7,8,Amp(1))
       call getamplqqqq(pin,2,8,7,1,Amp(2))

    elseif (chn .eq. 'qqq') then

       call getamplqqqq(pin,1,7,2,8,Amp(1))
       call getamplqqqq(pin,1,8,2,7,Amp(2))

    elseif (chn .eq. 'qbb') then

       call getamplqqqq(pin,7,1,8,2,Amp(1))
       call getamplqqqq(pin,7,2,8,1,Amp(2))

    endif
    Amp(2) = -Amp(2)  ! (-1) as swap fermions

  end subroutine getamplqqqq_stu


  !--------------------------------------------------------------!
  ! This calculates msq for different channels, for different    !
  ! flavour combinations.                                        !
  ! mtot(1) = both diagrams squared                              !
  ! mtot(2) = 'first' diagram squared                            !
  ! mtot(3) = 'second' diagram squared                           !
  ! e.g. in the qqb channel                                      !
  ! mtot(1) = |s+u|**2	(eg u dbar -> d ubar)                    !
  ! mtot(2) = |s|**2	(eg u dbar -> s cbar)                    !
  ! mtot(3) = |u|**2	(eg u sbar -> d cbar)                    !
  ! e.g. in the qqq channel                                      !       
  ! mtot(1) = |u+t|**2	(eg u u -> d d)                          !
  ! mtot(2) = |u|**2	(eg u c -> d s)                          !
  ! mtot(3) = 0 :  I  chose  u c -> d s to be in u channel       !
  !--------------------------------------------------------------!
  subroutine mtotqqqq(p,mtot,mtot_bits,chn)
    real(dp), intent(in)  :: p(8,4)
    character, intent(in) :: chn*3
    real(dp), intent(out) :: mtot(3)
    real(dp), intent(out) :: mtot_bits(3)
    !--------------------------------
    complex(dp) :: Amp(2)
    real(dp), parameter    :: c1=2d0, c2=-2d0/3d0

    call getamplqqqq_stu(p,chn,Amp)
    
    mtot_bits(1) = c1*abs(Amp(1))**2 
    mtot_bits(2) = c1*abs(Amp(2))**2
    mtot_bits(3) = c2*(Amp(1)*dconjg(Amp(2))+dconjg(Amp(1))*Amp(2))

    mtot(1) =   mtot_bits(1)+mtot_bits(2)+mtot_bits(3)
    mtot(2) =   mtot_bits(1)
    
    if (chn .eq. 'qqb' .or. chn .eq. 'qbq') then
       mtot(3) = mtot_bits(2) 
    elseif (chn .eq. 'qqq' .or. chn .eq. 'qbb') then
       mtot(3) = zero ! never used in qqb_wpwp_qqb.f 
    endif

  end subroutine mtotqqqq

end module qqqqampl

