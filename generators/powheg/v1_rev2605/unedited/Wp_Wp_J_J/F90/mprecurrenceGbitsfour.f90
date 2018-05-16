!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECrecurrenceGbitsfour.f90.

!! File generated automatically by autogen.pl from 
!! general precision template mpfiles/genmp/mprecurrenceGbitsfour.f90.

module mprecurrencebitsfour
  use mpmodule; use mpconverter
  use mpsimpleoperations; use mpadvancedoperations
  use mpsimpleoperations_c; use mpadvancedoperations_c

  use types; use consts_mp
  use mprecurrencebitsone
  use mprecurrencebitstwo
  use mprecurrencebitsthree
  use mpaux_functions
  use mpmemory
  use define_ampl
  use mpauxiliary_functions
  implicit none

  public :: vvww,fWW0,bfWW0,fVWW0,bfVWW0,fWW,bfWW,gWW_fbf
  public :: fWW_bffbf
  ! -- new functions 
  public :: fWWnogluon, bfWWnogluon 
  public :: bfWW_fbff, fWW_bffbf_1
  public :: gww_sbsfbf, gWW_bff

  interface gWW_fbf
     module procedure gWW_fbf_TM, gWW_fbf_RR
  end interface

  interface fWW
     module procedure fWW_TM, fWW_RR
  end interface

  interface bfWW
     module procedure bfWW_TM, bfWW_RR
  end interface

  interface fWW_bffbf
     module procedure fWW_bffbf_TM, fWW_bffbf_RR
  end interface


  private 

contains



  function vvww(eW1,kW1,eW2,kW2)  !---- The ZWW or \gammaWW vertex
    type(mp_complex), intent(in) :: eW1(:), eW2(:)
    type(mp_complex), intent(in) :: kW1(:), kW2(:)
    type(mp_complex)	 	    :: vvww(size(eW1))
    !------------------------------------------
    type(mp_complex):: sk1e2, se1e2, sk2e1, sk1e1, sk2e2

    sk1e2 = dot(kW1,eW2)
    se1e2 = dot(eW1,eW2)
    sk2e1 = dot(kW2,eW1)
    sk1e1 = dot(kW1,eW1)
    sk2e2 = dot(kW2,eW2)

    !----W1 is the negative one
    vvww = (se1e2*kW2-se1e2*kW1&
         &+(2*sk1e2+sk2e2)*eW1-(2*sk2e1+sk1e1)*eW2)

  end function vvww


  function fWW0(sp,p,eW1,kW1,eW2,kW2) !---Two Ws attach directly onto fermions
    type(mp_complex), intent(in)      :: sp(:), p(:)
    type(mp_complex), intent(in)    	 :: eW1(:), kW1(:),eW2(:), kW2(:)
    type(mp_complex)	      		 :: fWW0(size(sp))
    !----------------------------------------
    type(mp_complex)		:: k1(size(p))
    type(mp_complex)		:: sp1(size(sp))
    type(mp_complex)		:: k1sq
    type(mp_complex)		:: tmp(size(sp))


    fWW0 = czero

    !---W2 is nearest the 'onshell' fermion

    k1 = p + kW2(:)
    k1sq = dot(k1,k1)

    sp1 = vbqW(sp,eW2(:))

    tmp = spb2(sp1,k1)

    if (abs(k1sq) > propcut) then
       tmp = ci/k1sq*tmp
    else 
       tmp = czero
    endif

    fWW0 = vbqW(tmp,eW1(:))

  end function fWW0


  function bfWW0(sp,p,eW1,kW1,eW2,kW2) !---Two Ws attach directly onto fermions
    type(mp_complex), intent(in)      :: sp(:), p(:)
    type(mp_complex), intent(in)    	 :: eW1(:), kW1(:),eW2(:), kW2(:)
    type(mp_complex)	      		 :: bfWW0(size(sp))
    !----------------------------------------
    type(mp_complex)		:: k1(size(p))
    type(mp_complex)		:: sp1(size(sp))
    type(mp_complex)		:: k1sq
    type(mp_complex)		:: tmp(size(sp))


    bfWW0 = czero

    !---W2 is nearest the 'onshell' anti-fermion

    k1 = -p - kW2(:)
    k1sq = dot(k1,k1)

    sp1 = vWq(eW2(:),sp)
    tmp = spi2(k1,sp1)

    if (abs(k1sq) > propcut) then
       tmp = ci/k1sq*tmp
    else 
       tmp = czero
    endif

    bfWW0 = vWq(eW1(:),tmp)

  end function bfWW0



  function fVWW0(sp,p,eW,kW) !---WW -> vector boson -> fermion line
    type(mp_complex), intent(in)      :: sp(:), p(:)
    type(mp_complex), intent(in)    	 :: eW(:,:), kW(:,:)
    type(mp_complex)	      		 :: fVWW0(size(sp))
    !----------------------------------------
    type(mp_complex)		:: k1(size(p))
    type(mp_complex)		:: e1(size(sp))

    k1 = kW(:,1) + kW(:,2)
    e1 = vvww(eW(:,1),kW(:,1),eW(:,2),kW(:,2))

    fVWW0= vbqW(sp,e1)

  end function fVWW0



  function bfVWW0(sp,p,eW,kW) !---WW -> vector boson -> fermion line
    type(mp_complex), intent(in)      :: sp(:), p(:)
    type(mp_complex), intent(in)    	 :: eW(:,:), kW(:,:)
    type(mp_complex)	      		 :: bfVWW0(size(sp))
    !----------------------------------------
    type(mp_complex)		:: k1(size(p))
    type(mp_complex)		:: e1(size(sp))

    k1 = -kW(:,1) - kW(:,2)
    e1 = vvww(eW(:,1),kW(:,1),eW(:,2),kW(:,2))

    bfVWW0= vWq(e1,sp)

  end function bfVWW0


  recursive function fWW_RR(e,k,sp,p,flin, flout, eW,kW,ms,giarray,&
       &qiarray,WWid, pol_int) result(res)
    type(mp_complex), intent(in)     :: e(:,:),k(:,:), kW(:,:), eW(:,:)
    type(mp_complex), intent(in)     :: p(:), sp(:)
    character, intent(in)         :: flin*3, flout*3
    integer, intent(in), optional :: giarray(:),qiarray(:),WWid(:),pol_int,ms 
    ! -----------------------------------------------------------------------
    integer ngluons, ng1,ng2,m,n,ms1,i
    type(mp_complex):: e2(size(eW,dim=1)), res(size(sp,dim=1)), tmp(size(sp,dim=1)) 
    type(mp_complex):: kZ(size(kW,dim=1)),sp2(size(sp)),nm(size(sp2)),  kprop(size(p))
    type(mp_complex)::  den, den2
    type(mp_complex):: kgluon(size(k,dim=1))
    character:: flint*3
    logical:: done

    if (flin .ne. flout) then
       res = czero
       return
    endif

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'fw: giarray missing'
    endif


    res = czero
    if (flin .eq. 'bot') then
       flint = 'top'
    elseif (flin .eq. 'top') then
       flint = 'bot'
    endif

    ngluons = size(e,dim=2)

    if (ngluons ==0) then
       res = fWWnogluon(sp,p, eW,kW, giarray,qiarray,WWid, pol_int)
    else
       ng1 = ms
       ng2 = ngluons -ms


       ! A: Attach W1 to current with n gluons and W2. 
       ! Remember W1 is always the 
       ! upper W (may be W+ or W-, decided in main program).
       sp2 = fW(e,k,sp,p,flint, flout, eW(:,2),kW(:,2),ng1,&
            &giarray,qiarray,WWid(2),pol_int)

       if (ngluons >= 1) then
          kprop = sum(k(:,1:ngluons),dim=2)
       else
          kprop = czero
       endif
       kprop(:) = kprop(:)+kW(:,2) + p(:)
       sp2 = ci*(spb2(sp2,kprop) + mt*sp2)
       den = sc(kprop,kprop) - mt**2
       if (abs(den) > propcut) then
          tmp = 1/den*vbqW(sp2,eW(:,1))
       else
          tmp = czero
       endif
       res = res + tmp

       ! C: Recursion. Attach the nth gluon to a current of 2W's and n-1 gluons
       do m = 0,ng2-1

          if (ng1+m+1<=ngluons) then
             kgluon = sum(k(:,ng1+m+1:ngluons), dim=2)
          else
             kgluon = czero
          endif
          e2 = vgluon(e(:,ng1+m+1:ngluons),k(:,ng1+m+1:ngluons),&
               &giarray(ng1+m+1:ngluons),pol_int)

          if (1<=ng1+m) then
             kprop = sum(k(:,1:ng1+m), dim=2)
          else
             kprop = czero
          endif
          kprop(:) = kprop(:) + kW(:,1) +kW(:,2) +p(:) 
          den = sc(kprop,kprop) - mt**2

          sp2 = fWW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,flin,flout,&
               &eW,kW,ng1,giarray(1:ng1+m),qiarray,WWid,pol_int)

          sp2 = ci*(spb2(sp2,kprop) + mt*sp2)

          tmp = vqg(sp2,e2)

          if (ng2-m > 1) then     
             den2 = ci*sc(kgluon,kgluon)
          else
             den2 = 1                
          endif
          if (abs(den) > propcut) then
             tmp = tmp/den
          else
             tmp = czero
          endif

          if (abs(den2) > propcut) then
             tmp = tmp/den2
          else
             tmp = czero
          endif
          res = res + tmp
       enddo

       do m = 1,ng1
          if (m>1) then
             kgluon = sum(k(:,1:m), dim=2)
             den2 = ci*sc(kgluon,kgluon)
          else
             den2 = 1
          endif
          e2 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)

          if (m+1<= ngluons) then
             kprop = sum(k(:,m+1:ngluons),dim=2)
          else
             kprop = czero
          endif
          kprop(:) = kprop(:) + kW(:,1) +kW(:,2) +p(:)  !?
          den = sc(kprop,kprop) - mt**2
          ms1 = ng1-m

          sp2 = fWW(e(:,m+1:ngluons),k(:,m+1:ngluons),sp,p,&
               &flin,flout,eW,kW,ms1,giarray(m+1:ngluons),qiarray,WWid,pol_int)
          if (ngluons == 2) then

          endif
          sp2 = ci*(spb2(sp2,kprop) + mt*sp2)
          tmp = vgq(e2,sp2) 
          if (abs(den) > propcut) then
             tmp = tmp/den
          else
             tmp = czero
          endif
          if (abs(den2)>propcut) then
             tmp = tmp/den2
          else
             tmp = czero
          endif

          res = res +tmp
       enddo

    endif

    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,WWid(1),WWid(2))
  end function fWW_RR




  recursive function fWW_TM(e,k,sp,p,eW1,kW1,eW2,kW2,ms,giarray,&
       &qiarray,iWp,iWm,pol_int) result(res)
    type(mp_complex), intent(in)	:: e(:,:), k(:,:)
    type(mp_complex), intent(in)	:: sp(:), p(:)
    type(mp_complex), intent(in)	:: eW1(:), kW1(:)
    type(mp_complex), intent(in)	:: eW2(:), kW2(:)
    integer, intent(in)		:: ms
    integer, intent(in), optional :: giarray(:), qiarray(:),iWp,iWm,pol_int  
    !----------------------------------------
    type(mp_complex)			:: res(size(sp))
    type(mp_complex)			:: tmp(size(sp))
    type(mp_complex)			:: k1(size(p))
    type(mp_complex)			:: k2(size(kW1))
    type(mp_complex)			:: sp2(size(sp))
    type(mp_complex)			:: e1(size(e,dim=1))
    type(mp_complex)			:: k1sq, k2sq
    integer			:: ngluon, ng1, ng2, m, ms1
    character(len=3)		:: x
    logical                              :: done 

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,iWp,iWm)
       if (done) return 
    else
       if (i_warn < max_warn) then 
          write(*,*) 'bWW: giarray missing', i_warn 
          i_warn = i_warn+1
       endif
    endif


    ngluon = size(e,dim=2)
    ng1 = ms               !---#gluons to the left of the f-line
    ng2 = ngluon - ms	  !---#gluons to the right of the f-line

    if (ngluon==0) then
       res = fWW0(sp,p,eW1,kW1,eW2,kW2)
       return
    else

       res = czero

       do m=0, ng2-1
          if (ng1+1+m<=ngluon) then
             k1 =sum(k(:,ng1+1+m:ngluon),dim=2)
          else
             k1 = czero
          endif

          e1 = vgluon(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon),&
               &giarray(ng1+1+m:ngluon),pol_int)
          k1sq = sc(k1,k1)

          if (1<=ng1+m) then
             k2 = sum(k(:,1:ng1+m),dim=2)
          else
             k2 = czero
          endif

          k2 = k2 + p + kW1 + kW2
          k2sq = sc(k2,k2)
          sp2 = fWW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,eW1,kW1,eW2,kW2,&
               &ng1,giarray(1:ng1+m),qiarray,iWp,iWm,pol_int)

          sp2 = spb2(sp2,k2)

          tmp = vqg(sp2,e1)

          if (m < ng2-1) then
             if (abs(k1sq) > propcut) then
                tmp = -ci/k1sq*tmp
             else
                stop 'gluon prop blowing up'
             endif
          endif

          if (abs(k2sq) > propcut) then
             tmp = ci/k2sq*tmp
          else
             stop 'quark propagator blowing up'
          endif

          res = res + tmp

       enddo



       do m=1,ng1
          k1 = sum(k(:,1:m),dim=2)
          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)

          k1sq = sc(k1,k1)

          if (m+1<=ngluon) then
             k2 = sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero
          endif

          k2 = k2 + p + kW1 + kW2
          k2sq = sc(k2,k2)
          ms1 = ng1 - m
          sp2=fWW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,eW1,kW1,eW2,kW2,&
               &ms1,giarray(m+1:ngluon),qiarray,iWp,iWm,pol_int)

          sp2 = spb2(sp2,k2)

          tmp = vgq(e1,sp2)

          if (m > 1) then
             if (abs(k1sq) > propcut) then
                tmp = -ci/k1sq*tmp
             else
                stop 'gluon prop blowing up'
             endif
          endif

          if (abs(k2sq) > propcut) then
             tmp=ci/k2sq*tmp
          else
             stop 'quark prop blowing up'
          endif

          res = res + tmp

       enddo

       sp2 = fW(e,k,sp,p,'bot','top',eW2,kW2,ms,giarray,qiarray,iWm,pol_int)
       if (1<=ngluon) then
          k1 = sum(k,dim=2) + p + kW2
       else
          k1 = p + kW2
       endif

       k1sq = dot(k1,k1)

       sp2 = spb2(sp2,k1)

       tmp = vbqW(sp2,eW1)

       if (abs(k1sq) > propcut) then
          tmp = ci/k1sq*tmp
       else 
          stop 'quark prop blowing up'
       endif

       res = res + tmp

    endif

    ! -- store current 
    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,iWp,iWm)


  end function fWW_TM


  recursive function bfWW_RR(e,k,sp,p,flin, flout, eW,kW,ms,giarray,qiarray,WWid, pol_int) &
       &result(res)
    type(mp_complex), intent(in)::e(:,:),k(:,:), kW(:,:), eW(:,:)
    type(mp_complex), intent(in):: p(:), sp(:)
    character, intent(in) :: flin*3, flout*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int,ms 
    ! -----------------------------------------------------------------------
    integer ngluons, ng1,ng2,m,n,ms1,i
    type(mp_complex):: e2(size(eW,dim=1)), res(size(sp,dim=1)), tmp(size(sp,dim=1)) 
    type(mp_complex):: kZ(size(kW,dim=1)),sp2(size(sp)),nm(size(sp2)),  kprop(size(p))
    type(mp_complex)::  den, den2
    type(mp_complex):: kgluon(size(k,dim=1))
    character:: flint*3
    logical:: done

    if (flin .ne. flout) then
       res = czero
       return
    endif

    done = .false. 
    if (present(giarray)) then 
       ! XXX 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'fw: giarray missing'
    endif


    res = czero
    if (flin .eq. 'bot') then
       flint = 'top'
    elseif (flin .eq. 'top') then
       flint = 'bot'
    endif

    ngluons = size(e,dim=2)

    if (ngluons ==0) then
       res = bfWWnogluon(sp,p, eW,kW, giarray,qiarray,WWid, pol_int)
    else
       ng1 = ms
       ng2 = ngluons -ms


       ! A: Attach W1 to current with n gluons and W2. Remember W1 is always the 
       ! upper W (may be W+ or W-, decided in main program).
       sp2 = bfW(e,k,sp,p,flint, flout, eW(:,1),kW(:,1),ng1,giarray,qiarray,WWid(1),pol_int)

       if (ngluons >= 1) then
          kprop = sum(k(:,1:ngluons),dim=2)
       else
          kprop = czero
       endif
       kprop(:) = -kprop(:)-kW(:,1) - p(:)
       sp2 = ci*(spi2(kprop,sp2) + mt*sp2)
       den = sc(kprop,kprop) - mt**2
       if (abs(den) > propcut) then
          tmp = 1/den*vWq(eW(:,2),sp2)
       else
          tmp = czero
       endif
       res = res + tmp


!!$   ! C: Recursion. Attach the nth gluon to a current of 2W's and n-1 gluons
       do m = 0,ng2-1

          if (ng1+m+1<=ngluons) then
             kgluon = sum(k(:,ng1+m+1:ngluons), dim=2)
          else
             kgluon = czero
          endif
          e2 = vgluon(e(:,ng1+m+1:ngluons),k(:,ng1+m+1:ngluons),&
               &giarray(ng1+m+1:ngluons),pol_int)

          if (1<=ng1+m) then
             kprop = sum(k(:,1:ng1+m), dim=2)
          else
             kprop = czero
          endif
          kprop(:) = -kprop(:) - kW(:,1) -kW(:,2) -p(:) 
          den = sc(kprop,kprop) - mt**2

          sp2 = bfWW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,flin,flout,eW,kW,&
               &ng1,giarray(1:ng1+m),qiarray,WWid,pol_int)

          sp2 = ci*(spi2(kprop,sp2) + mt*sp2)

          tmp = vbqg(sp2,e2)

          if (ng2-m > 1) then     
             den2 = ci*sc(kgluon,kgluon)
          else
             den2 = 1                
          endif
          if (abs(den) > propcut) then
             tmp = tmp/den
          else
             tmp = czero
          endif

          if (abs(den2) > propcut) then
             tmp = tmp/den2
          else
             tmp = czero
          endif
          res = res + tmp

       enddo

       do m = 1,ng1
          if (m>1) then
             kgluon = sum(k(:,1:m), dim=2)
             den2 = ci*sc(kgluon,kgluon)
          else
             den2 = 1
          endif
          e2 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)

          if (m+1<= ngluons) then
             kprop = sum(k(:,m+1:ngluons),dim=2)
          else
             kprop = czero
          endif
          kprop(:) = -kprop(:) - kW(:,1) -kW(:,2) -p(:)  !?
          den = sc(kprop,kprop) - mt**2
          ms1 = ng1-m

          sp2 = bfWW(e(:,m+1:ngluons),k(:,m+1:ngluons),sp,p,flin,flout,&
               &eW,kW,ms1,giarray(m+1:ngluons),qiarray,WWid,pol_int)
          if (ngluons == 2) then

          endif
          sp2 = ci*(spi2(kprop,sp2) + mt*sp2)
          tmp = vgbq(e2,sp2) 
          if (abs(den) > propcut) then
             tmp = tmp/den
          else
             tmp = czero
          endif
          if (abs(den2)>propcut) then
             tmp = tmp/den2
          else
             tmp = czero
          endif

          res = res +tmp

       enddo

    endif

    if (present(giarray)) call store_result(pol_int,res,giarray,&
         &qiarray,WWid(1),WWid(2))
  end function bfWW_RR


  recursive function bfWW_TM(e,k,sp,p,eW1,kW1,eW2,kW2,ms,giarray,&
       &qiarray,iWp,iWm,pol_int) result(res)
    type(mp_complex), intent(in)	:: e(:,:), k(:,:)
    type(mp_complex), intent(in)	:: sp(:), p(:)
    type(mp_complex), intent(in)	:: eW1(:), kW1(:)
    type(mp_complex), intent(in)	:: eW2(:), kW2(:)
    integer, intent(in)		:: ms
    integer, intent(in), optional       :: giarray(:),qiarray(:),iWp,iWm,pol_int 
    !----------------------------------------
    type(mp_complex)			:: res(size(sp))
    type(mp_complex)			:: tmp(size(sp))
    type(mp_complex)			:: k1(size(p))
    type(mp_complex)			:: k2(size(kW1))
    type(mp_complex)			:: sp2(size(sp))
    type(mp_complex)			:: e1(size(e,dim=1))
    type(mp_complex)			:: k1sq, k2sq
    integer			        :: ngluon, ng1, ng2, m, ms1
    character(len=3)		        :: x
    logical                              :: done 

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,iWp,iWm)
       if (done) return 
    else
       if (i_warn < max_warn) then 
          write(*,*) 'bfWW: giarray missing', i_warn 
          i_warn = i_warn+1
       endif
    endif

    ngluon = size(e,dim=2)
    ng1 = ms               !---#gluons to the left of the f-line
    ng2 = ngluon - ms	  !---#gluons to the right of the f-line

    if (ngluon==0) then
       res = bfWW0(sp,p,eW1,kW1,eW2,kW2)
       return

    else

       res = czero

       do m=0, ng2-1
          if (ng1+1+m<=ngluon) then
             k1 =sum(k(:,ng1+1+m:ngluon),dim=2)
          else
             k1 = czero
          endif

          e1 = vgluon(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))

          k1sq = sc(k1,k1)

          if (1<=ng1+m) then
             k2 = sum(k(:,1:ng1+m),dim=2)
          else
             k2 = czero
          endif

          k2 = -k2 - p - kW1 - kW2
          k2sq = sc(k2,k2)
          sp2 = bfWW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,eW1,kW1,eW2,kW2,ng1,&
               &giarray(1:ng1+m),qiarray,iWp,iWm,pol_int)

          sp2 = spi2(k2,sp2)

          tmp = vbqg(sp2,e1)

          if (m < ng2-1) then
             if (abs(k1sq) > propcut) then
                tmp = -ci/k1sq*tmp
             else
                stop 'gluon prop blowing up'
             endif
          endif

          if (abs(k2sq) > propcut) then
             tmp = ci/k2sq*tmp
          else
             stop 'quark propagator blowing up'
          endif

          res = res + tmp

       enddo

       do m=1,ng1
          k1 = sum(k(:,1:m),dim=2)
          e1 = vgluon(e(:,1:m),k(:,1:m))

          k1sq = sc(k1,k1)

          if (m+1<=ngluon) then
             k2 = sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero
          endif

          k2 = -k2 - p - kW1 - kW2
          k2sq = sc(k2,k2)
          ms1 = ng1 - m
          sp2=bfWW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,eW1,kW1,eW2,kW2,ms1,&
               &giarray(m+1:ngluon),qiarray,iWp,iWm,pol_int)

          sp2 = spi2(k2,sp2)

          tmp = vgbq(e1,sp2)

          if (m > 1) then
             if (abs(k1sq) > propcut) then
                tmp = -ci/k1sq*tmp
             else
                stop 'gluon prop blowing up'
             endif
          endif

          if (abs(k2sq) > propcut) then
             tmp=ci/k2sq*tmp
          else
             stop 'quark prop blowing up'
          endif

          res = res + tmp

       enddo

       sp2 = bfW(e,k,sp,p,x,x,eW2,kW2,ms,giarray,qiarray,iWm,pol_int)
       if (1<=ngluon) then
          k1 = -sum(k,dim=2) - p - kW2
       else
          k1 = -p - kW2
       endif
       k1sq = dot(k1,k1)

       sp2 = spi2(k1,sp2)

       tmp = vWq(eW1,sp2)

       if (abs(k1sq) > propcut) then
          tmp = ci/k1sq*tmp
       else 
          stop 'quark prop blowing up'
       endif

       res = res + tmp

    endif


    ! -- store current 
    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,iWp,iWm)

  end function bfWW_TM


  recursive function gWW_fbf_RR(e,k,sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
       &giarray,qiarray,WWid,pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:)
    type(mp_complex), intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    type(mp_complex), intent(in) :: eW(:,:),kW(:,:)
    integer, intent(in) ::  ng1,ng2
    character, intent(in) :: fl1*3,fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int 
    type(mp_complex)             :: res(size(e,dim=1))
    ! -------------------------------------
    type(mp_complex) :: eW1(size(ew,dim=1)),eW2(size(ew,dim=1))
    type(mp_complex) :: kW1(size(ew,dim=1)),kW2(size(ew,dim=1))
    integer :: iWp, iWm 
    eW1 = ew(:,1)
    eW2 = ew(:,2)
    kW1 = kw(:,1)
    kW2 = kw(:,2)
    iWp = WWid(1)
    iWm = WWid(2)

    res= gWW_fbf(e,k,sp1,p1,fl1,sp2,p2,fl2,&
         &eW1,kW1,eW2,kW2,ng1,ng2,giarray,qiarray,iWp,iWm,pol_int)

  end function gWW_fbf_RR
  !
  recursive function gWW_fbf_TM(e,k,sp1,p1,fl1,sp2,p2,fl2,&
       &eW1,kW1,eW2,kW2,ng1,ng2,giarray,qiarray,iWp,iWm,pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:)
    type(mp_complex), intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    type(mp_complex), intent(in) :: eW1(:),kW1(:),eW2(:),kW2(:)
    integer, intent(in) ::  ng1,ng2
    character, intent(in) :: fl1*3,fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),iWp,iWm,pol_int 
    ! -----------------------------------------------------------------------
    integer             :: m1,m2, ms1a, ms2a,Dv!,m,m3
    integer :: ngluon,  ng3, ngL 
    type(mp_complex)             :: res(size(e,dim=1))
    type(mp_complex)             :: tmp(size(e,dim=1)),tmp2(size(e,dim=1))
    type(mp_complex)             :: k1(size(p1))
    type(mp_complex)             :: k2(size(p1)),k22(size(p1))
    type(mp_complex)             :: k3(size(p1))
    type(mp_complex)             :: k4(size(p1))
    type(mp_complex)             :: kW(size(p1))
    type(mp_complex)             :: sp3(size(sp1))
    type(mp_complex)             :: sp4(size(sp1)),sp5(size(sp1))
    type(mp_complex)             :: e1(size(e,dim=1))
    type(mp_complex)             :: e2(size(e,dim=1))
    type(mp_complex)             :: e3(size(e,dim=1))
    type(mp_complex)  :: k1sq,k2sq,k22sq,k3sq,k4sq
    logical                   :: done 

    if (fl1.eq.'str'.or.fl2.eq.'str') then 
       res = czero 
       return 
    endif

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,iWp,iWm)
       if (done) return 
    else
       if (i_warn < max_warn) then 
          write(*,*) 'gWW_fbf: giarray missing', i_warn 
          i_warn = i_warn+1
       endif
    endif

    kW = kW1 + kW2

    Dv = size(e,dim=1)

    ngluon = size(e,dim=2)
    ng3 = ngluon - ng1 - ng2

    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT D'

    if (ngluon == 0) then 

       res = czero

       sp4 = fWW(e,k,sp2,p2,eW1,kW1,eW2,kW2,0,giarray,qiarray(2:2),iWp,iWm,pol_int)   
       k2 = p2 + kW1 + kW2
       k2sq = sc(k2,k2)
       sp4 = ci/k2sq*(spb2(sp4,k2))!+mass*sp4)
       res = -cone*vbqq(Dv,sp4,sp1)

       sp4 = bfWW(e,k,sp1,p1,eW2,kW2,eW1,kW1,0,giarray,qiarray(1:1),iWm,iWp,pol_int)   
       k2 = -p1 - kW1 - kW2
       k2sq = sc(k2,k2)
       sp4 = ci/k2sq*(spi2(k2,sp4) )!+ mass*sp4)
       tmp = -cone*vbqq(Dv,sp2,sp4)

       sp4 = fW(e,k,sp2,p2,fl2,fl1,eW2,kW2,0,giarray,qiarray(2:2),iWm,pol_int)
       k2 = p2 + kW2
       k2sq = sc(k2,k2)
       sp4 = ci/k2sq*(spb2(sp4,k2))
       sp5 = bfW(e,k,sp1,p1,fl1,fl2,eW1,kW1,0,giarray,qiarray(1:1),iWp,pol_int)
       k22 = -p1 - kW1
       k22sq = sc(k22,k22)
       sp5 = ci/k22sq*(spi2(k22,sp5))
       tmp2 = -cone*vbqq(Dv,sp4,sp5)

       res = res + tmp + tmp2

    else

       res = czero

       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq = sc(k1,k1)

          ms1a = ng1-m1
          e2=gWW_fbf(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
               &sp1,p1,fl1,sp2,p2,fl2,eW1,kW1,eW2,kW2,ms1a,ng2,&
               &giarray(m1+1:ngluon),qiarray(1:2),iWp,iWm,pol_int)    
          if (m1+1<=ngluon) then 
             k2 = sum(k(:,m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p1 + p2+kW
          k2sq = sc(k2,k2)

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*vggg(e1,k1,e2,k2)
          else 
             tmp = czero 
          endif

          if (m1 > 1) then  
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo  !#1


       do m1=0,ng3-1 

          e1 = gWW_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
               &sp1,p1,fl1,sp2,p2,fl2,eW1,kW1,eW2,kW2,ng1,ng2,&
               &giarray(1:ng1+ng2+m1),qiarray,iWp,iWm,pol_int)    
          if (1<=ng1+ng2+m1) then 
             k1=sum(k(:,1:ng1+ng2+m1),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2+kW
          k1sq=sc(k1,k1)

          e2=vgluon(e(:,ng1+ng2+m1+1:ngluon),&
               &k(:,ng1+ng2+m1+1:ngluon),giarray(ng1+ng2+m1+1:ngluon),pol_int)    
          if (ng1+ng2+m1+1<=ngluon) then 
             k2=sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2sq = sc(k2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = -ci/k1sq*vggg(e1,k1,e2,k2)
          else
             tmp = czero 
          endif

          if (m1 < ng3-1) then 
             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo   !#2


       do m1 =1, ng1-1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = m1+1,ng1   

             e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),giarray(m1+1:m2),pol_int)    
             if (m1+1<=m2) then 
                k2=sum(k(:,m1+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)


             ms1a=ng1-m2
             e3=gWW_fbf(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW1,kW1,eW2,kW2,ms1a,ng2,&
                  &giarray(m2+1:ngluon),qiarray,iWp,iWm,pol_int)    
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3 = k3 + p1 + p2+kW
             k3sq=sc(k3,k3)

             if (abs(k3sq) > propcut) then 
                tmp = -ci/k3sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m1>1) then 
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2 > m1+1) then 
                if (abs(k2sq) > propcut) then 
                   tmp = -ci/k2sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo   !#3


       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = 0,ng3-1    

             ms1a=ng1-m1
             e2=gWW_fbf(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW1,kW1,eW2,kW2,ms1a,ng2,&
                  &giarray(m1+1:ng1+ng2+m2),qiarray,iWp,iWm,pol_int)    
             if (m1+1<=ng1+ng2+m2) then 
                k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p1 + p2 + kW
             k2sq=sc(k2,k2)


             e3=vgluon(e(:,ng1+ng2+m2+1:ngluon)&
                  &,k(:,ng1+ng2+m2+1:ngluon),&
                  &giarray(ng1+ng2+m2+1:ngluon),pol_int)    
             if (ng1+ng2+m2+1<=ngluon) then 
                k3=sum(k(:,ng1+ng2+m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)

             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m1>1) then 
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2+1+ng1+ng2 < ngluon) then 
                if (abs(k3sq) > propcut) then 
                   tmp = -ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo  !#4



       do m1=0,ng3-2

          ngL=ng1+ng2+m1

          e1=gWW_fbf(e(:,1:ngL),k(:,1:ngL),&
               &sp1,p1,fl1,sp2,p2,fl2,eW1,kW1,eW2,kW2,ng1,ng2,&
               &giarray(1:ngL),qiarray,iWp,iWm,pol_int)    
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2+kW
          k1sq=sc(k1,k1)

          do m2=ngL+1,ngluon-1

             e2=vgluon(e(:,ngL+1:m2),&
                  &k(:,ngL+1:m2),giarray(ngL+1:m2),pol_int)    
             if (ngL+1<=m2) then 
                k2=sum(k(:,ngL+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)

             e3=vgluon(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &giarray(m2+1:ngluon),pol_int)    
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)


             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m2 > ngL+1) then 
                if (abs(k2sq)> propcut) then 
                   tmp=-ci/k2sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2 < ng3-1) then 
                if (abs(k3sq) > propcut) then  
                   tmp=-ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo  !#5


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=bf(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl1,ng1,&
               &giarray(1:ms1a),qiarray(1:1),pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = -k3 - p1
          k3sq = sc(k3,k3)!-mass2

          if (ng1 > 0.or.m1 > 0) sp3 = spi2(k3,sp3)!+mass*sp3

          ms2a=ng2-m1
          sp4=fWW(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,eW1,kW1,eW2,kW2,ms2a, &
               &giarray(ms1a+1:ngluon),qiarray(2:2),iWp,iWm,pol_int)      
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 =  k4 + p2+kW
          k4sq = sc(k4,k4)!-mass2
          sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)

          if (ng1 > 0.or.m1 > 0) then 
             if (abs(k3sq) > propcut) then 
                tmp= ci/k3sq*tmp
             else 
                tmp = czero 
             endif
          endif


          if (abs(k4sq) > propcut) then 
             tmp= ci/k4sq*tmp 
          else 
             tmp = czero 
          endif


          res = res + tmp

       enddo  !#6


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=bfWW(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,eW2,kW2,eW1,kW1,ng1,&
               &giarray(1:ms1a),qiarray(1:1),iWm,iWp,pol_int) ! XXXX 
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = -k3 - p1 - kW
          k3sq = sc(k3,k3)!-mass2
          sp3 = spi2(k3,sp3)!+mass*sp3

          ms2a=ng2-m1
          sp4=f(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl2,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 =  k4 + p2
          k4sq = sc(k4,k4)!-mass2

          if (ng3 > 0.or.ng2-m1>0) sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)

          if (abs(k3sq) > propcut) then 
             tmp= ci/k3sq*tmp
          else 
             tmp = czero 
          endif


          if (ng3 > 0.or.ng2-m1>0) then 
             if (abs(k4sq) > propcut) then 
                tmp= ci/k4sq*tmp 
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo  !#7


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=bfW(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl2,eW1,kW1,ng1,&
               &giarray(1:ms1a),qiarray(1:1),iWp,pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = -k3 - p1 - kW1
          k3sq = sc(k3,k3)!-mass2
          sp3 = spi2(k3,sp3)!+mass*sp3

          ms2a=ng2-m1
          sp4=fW(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl1,eW2,kW2,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),iWm,pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 =  k4 + p2+kW2
          k4sq = sc(k4,k4)!-mass2
          sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)

          if (abs(k3sq) > propcut) then 
             tmp= ci/k3sq*tmp
          else 
             tmp = czero 
          endif

          if (abs(k4sq) > propcut) then 
             tmp= ci/k4sq*tmp 
          else 
             tmp = czero 
          endif


          res = res + tmp

       enddo  !#8

    endif   !end condition for ngluon 



    ! -- store current 
    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,iWp,iWm)



  end function gWW_fbf_TM


  !----------------------------------------------------------

  ! This is the current of fWW (no gluons) with the W's coupling directly to the fermion line.
  function fWWnogluon(sp,p, eW,kW, giarray,qiarray,WWid, pol_int) result(res)
    type(mp_complex), intent(in):: kW(:,:), eW(:,:)
    type(mp_complex), intent(in):: p(:), sp(:)
    integer, intent(in), optional:: giarray(:),qiarray(:),WWid(:),pol_int 
    ! -----------------------------------------------------------------------
    type(mp_complex):: res(size(sp,dim=1)), tmp(size(sp,dim=1)),e2(size(eW,dim=1))
    type(mp_complex)::  sp2(size(sp)), kprop(size(p))
    logical done

    !! The order of the W's will depend on the flavour of the
    !! quark. This is sorted out in a parent program, so that we can
    !! assume here that the first W entered is the top one, and the
    !! second is the bottom one, irrespective of their charges.


    done = .false.
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'fw: giarray missing'
    endif

    res = czero
    sp2 = vbqW(sp,eW(:,2))
    ! Momentum of quark propagator
    kprop(:) = p(:) + kW(:,2)
    sp2 = spb2(sp2,kprop) !+ mt*sp2
    tmp = ci*vbqW(sp2,eW(:,1))/(sc(kprop,kprop))! - mt**2)
    res = res +tmp


  end function fWWnogluon

  ! This is the current of fWW (no gluons) with the W's coupling directly to the fermion line.
  function bfWWnogluon(sp,p, eW,kW, giarray,qiarray,WWid, pol_int) result(res)
    type(mp_complex), intent(in):: kW(:,:), eW(:,:)
    type(mp_complex), intent(in):: p(:), sp(:)
    integer, intent(in), optional:: giarray(:),qiarray(:),WWid(:),pol_int 
    ! -----------------------------------------------------------------------
    type(mp_complex):: res(size(sp,dim=1)), tmp(size(sp,dim=1)),e2(size(eW,dim=1))
    type(mp_complex)::  sp2(size(sp)), kprop(size(p))
    logical done

    !! The order of the W's will depend on the flavour of the
    !! quark. This is sorted out in a parent program, so that we can
    !! assume here that the first W entered is the top one, and the
    !! second is the bottom one, irrespective of their charges.


    done = .false.
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'fw: giarray missing'
    endif

    res = czero

    sp2 = vWq(eW(:,1),sp)
    ! Momentum of quark propagator
    kprop(:) = -p(:) - kW(:,1)
    sp2 = spi2(kprop,sp2)+ mt*sp2
    tmp = ci*vWq(eW(:,2),sp2)
    tmp = tmp/(sc(kprop,kprop) - mt**2)
    res = res +tmp

  end function bfWWnogluon

  recursive function fWW_bffbf_RR(e,k,sp,p,fll,fl0,&
       &eWin,kWin,ng1,ng2,ng3,sw,giarray,qiarray,WWid,pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:)
    type(mp_complex), intent(in) :: sp(:,:),p(:,:)
    type(mp_complex), intent(in) :: eWin(:,:),kWin(:,:)
    integer, intent(in) ::  ng1,ng2,ng3,sw
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int   ! --------------------------------------------------------------------
    type(mp_complex)::sp4(size(sp,dim=1)), res(size(sp,dim=1)),tmp(size(sp,dim=1))
    type(mp_complex):: e2(size(e,dim=1)),sp1(size(sp,dim=1)),sp2(size(sp,dim=1))
    type(mp_complex)::kW(size(kWin,dim=1), size(kWin,dim=2)), eW(size(eWin,dim=1), size(eWin,dim=2))
    character::flaux*3,fl1*3,fl2*3,fl3*3, fllnew(3)*3
    integer::ngluons,ng4,m,ngL
    integer, parameter::Ndumm = 0
    type(mp_complex)             :: kdumm(size(k,dim=1),Ndumm),k5sq
    type(mp_complex)             :: edumm(size(e,dim=1),Ndumm)
    type(mp_complex)::k2sq,k4sq,k1sq,k4(size(k,dim=1)), k2(size(k,dim=1)), k5(size(k,dim=1))
    type(mp_complex)::k1(size(k,dim=1)), spnew(size(sp,dim=1),size(sp,dim=2)),pnew(size(p,dim=1),size(p,dim=2))


    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)
    kW(:,1) = kWin(:,1)
    kW(:,2) = kWin(:,2)

    eW(:,1) = eWin(:,1)
    eW(:,2) = eWin(:,2)


    ngluons = size(e,dim=2)


    ng4 = ngluons - ng1 - ng2-ng3
    res = czero
    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C: fWW_bffbf', ngluons, ng1, ng2, ng3, sw
    if (ngluons == 0) then

       if (sw == 1) then

          !--- Both W's above gluon

          sp4  = fWW(edumm, kdumm, sp(:,1),p(:,1), fl1, fl0,eW, kW, 0, giarray,&
               &qiarray(1:1), WWid, pol_int)
          k4 = p(:,1)+kW(:,1)+kW(:,2)
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4, k4)!+sp4*mass

          e2 = g_fbf(edumm, kdumm, sp(:,2), p(:,2), fl2, sp(:,3),p(:,3),fl3,0,0,&
               &giarray,qiarray(2:3),pol_int) 
          k2 = p(:,2)+p(:,3)
          k2sq = sc(k2,k2)

          tmp = vqg(sp4, e2)


          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          res = res+tmp   !  sw = 1, #1

          ! ---- 1 W below gluon, 1 above, and both below

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'
          sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,2), kW(:,2),&
               &ng1,ng2,ng3,1,giarray,qiarray,WWid(2),pol_int)    

          k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,2)
          sp4 = spb2(sp4,k4) !+ mass*sp4

          k4sq = sc(k4,k4)

          tmp = vbqW(sp4,eW(:,1))

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp
          else 
             tmp = czero
          endif

          res = res + tmp   ! sw=1,#2


       endif

       if (sw == 2) then

          ! -- for W above gluon line
          if (fl0 == 'top') flaux = 'bot'
          if (fl0 == 'bot') flaux = 'top'
          if (case_b2 .eqv. .false.) then


             sp4 = fW(edumm,kdumm, sp(:,1),p(:,1), fl1,fl0, eW(:,2),&
                  &kW(:,2),0,giarray,qiarray(1:1),WWid(2), pol_int)

             k4 = p(:,1)+kW(:,2)          
             k4sq = sc(k4,k4)

             sp4 = spb2(sp4, k4)!+mass*sp4

             e2 = gW_fbf(edumm,kdumm, sp(:,2), p(:,2), fl2, sp(:,3), p(:,3), fl3,&
                  &eW(:,1),kW(:,1),0,0, giarray,qiarray(2:3), WWid(1),pol_int)

             k2 = p(:,2)+p(:,3)+kW(:,1)
             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)

             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else 
                tmp = czero 
             endif

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp 
             else 
                tmp = czero 
             endif

             res = res + tmp  ! sw =2, #1



             ! ---- for W below gluon line          

             sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,1), kW(:,1),&
                  &ng1,ng2,ng3,3,&
                  &giarray,qiarray,WWid(1),pol_int)    

             k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,1)

             sp4 = spb2(sp4,k4) !+ mass*sp4
             k4sq = sc(k4,k4)

             tmp = vbqW(sp4,eW(:,2))


             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp   ! sw = 2, #2

          elseif (case_b2) then

             sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,1), kW(:,1),&
                  &ng1,ng2,ng3,sw,&
                  &giarray,qiarray,WWid(1),pol_int)    

             k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,1)

             sp4 = spb2(sp4,k4) !+ mass*sp4
             k4sq = sc(k4,k4)

             tmp = vbqW(sp4,eW(:,2))

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp   ! sw = 2, case_b2, #1

          endif

       endif

       if (sw ==3) then
          e2 = gWW_fbf(edumm, kdumm, sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,& 
               &eW,kW,0,0,giarray,qiarray(2:3),WWid,pol_int)

          k2 = p(:,2)+p(:,3)+kW(:,1)+kW(:,2)
          k2sq = sc(k2,k2)

          tmp = vqg(sp(:,1),e2)


          if (abs(k2sq) > propcut.and.fl0==fl1) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          res = res + tmp  ! #1
       endif

       ! This case is used for A3.
       if (sw==4) then 

          eW(:,1) = eWin(:,2)
          eW(:,2) = eWin(:,1)
          kW(:,1) = kWin(:,2)
          kW(:,2) = kWin(:,1)

          spnew(:,1) = sp(:,3)
          spnew(:,2) = sp(:,2)
          spnew(:,3) = sp(:,1)
          pnew(:,1) = p(:,3)
          pnew(:,2) = p(:,2)
          pnew(:,3)= p(:,1)
          fllnew(1) = fll(3)
          fllnew(2) = fll(2)
          fllnew(3) = fll(1)


          res = fWW_bffbf(e,k,spnew,pnew,fllnew,fl0,eW,kW,ng1,ng2,ng3,1,giarray,qiarray,WWid, pol_int)

       endif

       ! Also for A3...

       if (sw == 5) then

          eW(:,1) = eWin(:,2)
          eW(:,2) = eWin(:,1)
          kW(:,1) = kWin(:,2)
          kW(:,2) = kWin(:,1)


          spnew(:,1) = sp(:,3)
          spnew(:,2) = sp(:,2)
          spnew(:,3) = sp(:,1)
          pnew(:,1) = p(:,3)
          pnew(:,2) = p(:,2)
          pnew(:,3)= p(:,1)
          fllnew(1) = fll(3)
          fllnew(2) = fll(2)
          fllnew(3) = fll(1)
          res = fWW_bffbf(e,k,spnew,pnew,fllnew,fl0,eW,kW,ng1,ng2,ng3,3,giarray,qiarray,WWid, pol_int)

       endif

    else                     ! ngluons > 0


       if (sw == 2) then

          if (fl0 == 'bot') flaux = 'top'
          if (fl0 == 'top') flaux = 'bot'

          do m = 0,ng2

             sp4 = fW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),fl1,fl0,&
                  &eW(:,2),kW(:,2),ng1,giarray(1:ng1+m),qiarray(1:1),WWid(2),&
                  & pol_int)

             k4 = sum(k(:,1:ng1+m),dim=2)+ p(:,1) + kW(:,2)

             k4sq = sc(k4,k4) !- mass2

             sp4 = spb2(sp4,k4) !+ mass*sp4

             e2 = gW_fbf(e(:,ng1+m+1:ngluons),k(:,ng1+m+1:ngluons),sp(:,2),&
                  p(:,2),fll(2),sp(:,3),p(:,3),fll(3), eW(:,1),kW(:,1), &
                  & ng2-m, ng3, giarray(ng1+m+1:ngluons),qiarray(2:3),WWid(1),&
                  & pol_int)

             k2 = sum(k(:,ng1+m+1:ngluons),dim=2)+kW(:,1)+p(:,2)+p(:,3)

             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)

             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif

             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif

             res = res + tmp

          enddo


          do m = 1,ng1

             sp4 = fWW_bffbf(e(:,m+1:ngluons),k(:,m+1:ngluons),sp,p,fll,fl0, &
                  & eW,kW,ng1-m,ng2,ng3,sw,giarray(m+1:ngluons),qiarray,&
                  & WWid,pol_int)

             k4 = sum(k(:,m+1:ngluons),dim=2) + p(:,1)+p(:,2)+p(:,3) + &
                  & kW(:,1) + kW(:,2)

             k4sq = sc(k4,k4) !- mass2

             sp4 = spb2(sp4,k4) !+ mass*sp4

             e2 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)


             k2 = sum(k(:,1:m),dim=2)
             k2sq = sc(k2,k2)

             tmp = vgq(e2,sp4)


             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif


             if (m>1) then
                if (abs(k2sq) >  propcut) then
                   tmp = -ci*tmp/k2sq
                else
                   tmp = czero
                endif
             endif

             res = res + tmp

          enddo


          do m = 0,ng4-1

             ngL = ng1 + ng2 + ng3 + m

             sp4 = fWW_bffbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0, &
                  & eW,kW,ng1,ng2,ng3,sw,giarray(1:ngL),qiarray,&
                  & WWid,pol_int)


             k4 = sum(k(:,1:ngL),dim=2) + sum(p(:,1:3),dim=2) + &
                  & sum(kW(:,1:2),dim=2)

             sp4 = spb2(sp4,k4) !+ mass*sp4

             k4sq = sc(k4,k4) !- mass2

             e2 = vgluon(e(:,ngL+1:ngluons),k(:,ngL+1:ngluons),&
                  & giarray(ngL+1:ngluons),pol_int)

             k2 = sum(k(:,ngL+1:ngluons),dim=2)
             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)


             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif


             if (m < ng4-1) then
                if (abs(k2sq) >  propcut) then
                   tmp = -ci*tmp/k2sq
                else
                   tmp = czero
                endif
             endif

             res = res + tmp


          enddo

          sp4 = fW_bffbf(e,k,sp,p,fll,flaux,eW(:,1),kW(:,1),ng1,ng2,ng3,3,&
               & giarray, qiarray,WWid(1),pol_int)

          k4 = sum(k(:,1:ngluons),dim=2) + kW(:,1) + sum(p(:,1:3),dim=2)

          k4sq = sc(k4,k4) !- mass2

          sp4 = spb2(sp4,k4) !+ mass*sp4


          tmp = vbqW(sp4,eW(:,2))

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          res = res + tmp

       endif
    endif


  end function fWW_bffbf_RR

  recursive function fWW_bffbf_TM(e,k,sp,p,fll,fl0,&
       &eWp,kWp,eWm,kWm,ng1,ng2,ng3,sw,giarray,qiarray,iWp,iWm,pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:)
    type(mp_complex), intent(in) :: sp(:,:),p(:,:)
    type(mp_complex), intent(in) :: eWp(:),kWp(:),eWm(:),kWm(:)
    integer, intent(in) ::  ng1,ng2,ng3,sw
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),iWp,iWm,pol_int 
    ! -----------------------------------------------------------------------
    character :: flaux*3
    character :: fl1*3,fl2*3,fl3*3
    integer :: ngluon, ng4, ngL,m!,m1
    integer, parameter :: Ndumm=0
    type(mp_complex)             :: res(size(sp,dim=1))
    type(mp_complex)             :: tmp(size(sp,dim=1))
    type(mp_complex)             :: k1(size(k,dim=1))
    type(mp_complex)             :: k2(size(k,dim=1))
    type(mp_complex)             :: k4(size(k,dim=1))
    type(mp_complex)             :: sp1(size(sp,dim=1))
    type(mp_complex)             :: sp2(size(sp,dim=1))
    type(mp_complex)             :: sp3(size(sp,dim=1))
    type(mp_complex)             :: sp4(size(sp,dim=1))
    type(mp_complex)             :: sp5(size(sp,dim=1))
    type(mp_complex)             :: e1(size(e,dim=1))
    type(mp_complex)             :: e2(size(e,dim=1))
    type(mp_complex)             :: kdumm(size(k,dim=1),Ndumm)
    type(mp_complex)             :: edumm(size(e,dim=1),Ndumm)
    type(mp_complex)             :: k1sq,k2sq,k4sq
    logical                   :: done 


    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,iWp,iWm)
       if (done) return 
    else
       if (i_warn < max_warn) then 
          write(*,*) 'fWW_bffbf_2W: giarray missing', i_warn 
          i_warn = i_warn+1
       endif
    endif

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)

    if ((sw.ne.3).and.(sw.ne.1).and.(sw.ne.2)) stop 'sw in fW_bffbf_2W not implemented'

    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C'

    if (ngluon == 0) then 

       res = czero

       if (sw.eq.3) then 

          e2 = gWW_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
               &sp(:,3),p(:,3),fl3,eWp,kWp,eWm,kWm,0,0,&
               &giarray,qiarray(2:3),iWp,iWm,pol_int)  
          k1 = p(:,2)+p(:,3)+kWp+kWm
          k1sq=sc(k1,k1)

          if (abs(k1sq) > propcut.and.fl0==fl1) then 
             tmp = -ci/k1sq*vqg(sp(:,1),e2)
          else 
             tmp = czero 
	     write(*,*)'czero in fWW_fbffbf tmp'
          endif

          res = res + tmp  
       endif

       if ((sw.eq.1)) then 

          e2 = g_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
               &sp(:,3),p(:,3),fl3,0,0,&
               &giarray,qiarray(2:3),pol_int)    

          k2 = p(:,2)+p(:,3)
          k2sq = sc(k2,k2)

          sp4 = fWW(edumm,kdumm,sp(:,1),p(:,1),eWp,kWp,eWm,kWm,0,&
               &giarray,qiarray(1:1),iWp,iWm,pol_int)     

          k4  = p(:,1) + kWp + kWm
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4,k4) !+ mass*sp4

          tmp = vqg(sp4,e2)

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp 
          else 
             tmp = czero 
          endif

          res = res + tmp  ! #3

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'

          sp5 = fW_bffbf_2W(edumm,kdumm,sp,p,fll,flaux,&
               &eWm,kWm,0,0,0,1,&
               &giarray,qiarray,iWm,pol_int)  

          k4 = p(:,1)+p(:,2)+p(:,3) + kWm
          sp5 = spb2(sp5,k4) !+ mass*sp5
          k4sq = sc(k4,k4)

          tmp = vbqW(sp5,eWp)

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp
          else 
             tmp = czero
          endif

          res = res + tmp  

       endif

       if ((sw.eq.2)) then 

          if (case_b2 .eqv. .false.) then
             e2 = gW_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
                  &sp(:,3),p(:,3),fl3,eWp,kWp,0,0,&
                  &giarray,qiarray(2:3),iWp,pol_int) 

             k2 = p(:,2)+p(:,3) + kWp
             k2sq = sc(k2,k2)

             sp4 = fW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl0,eWm,kWm,0,&
                  &giarray,qiarray(1:1),iWm,pol_int)    
             k4  = p(:,1) + kWm
             k4sq = sc(k4,k4)
             sp4 = spb2(sp4,k4) !+ mass*sp4

             tmp = vqg(sp4,e2)

             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else 
                tmp = czero 
             endif

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp 
             else 
                tmp = czero 
             endif

             res = res + tmp

             if (fl0.eq.'top') flaux = 'bot'
             if (fl0.eq.'bot') flaux = 'top'

             sp5 = fW_bffbf_2W(edumm,kdumm,sp,p,fll,flaux,&
                  &eWp,kWp,0,0,0,3,&
                  &giarray,qiarray,iWp,pol_int)

             k4 = p(:,1)+p(:,2)+p(:,3) + kWp
             sp5 = spb2(sp5,k4) !+ mass*sp5
             k4sq = sc(k4,k4)

             tmp = vbqW(sp5,eWm)

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp  

          elseif (case_b2) then

             if (fl0.eq.'top') flaux = 'bot'
             if (fl0.eq.'bot') flaux = 'top'

             sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eWp, kWp,&
                  &ng1,ng2,ng3,sw,&
                  &giarray,qiarray,iWp,pol_int)    
             k4 = p(:,1)+p(:,2)+p(:,3)+kWp

             sp4 = spb2(sp4,k4) !+ mass*sp4
             k4sq = sc(k4,k4)

             tmp = vbqW(sp4,eWm)

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp   ! sw = 2, case_b2, #1

          endif


       endif


    else                              ! ngluons > 0

       res = czero
       if (sw == 1) then

          do m = 0,ng2

             sp4 = fWW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),eWp,kWp,eWm,kWm,&
                  &ng1,giarray(1:ng1+m),qiarray(1:1),iWp,iWm,pol_int)
             k4 = sum(k(:,1:ng1+m),dim=2) + p(:,1)+kWp+kWm
             k4sq = sc(k4,k4)
             sp4 = spb2(sp4,k4)

             e2 = g_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,&
                  &ng2-m,ng3,giarray(ng1+1+m:ngluon),qiarray(2:3),pol_int)
             k2 = sum(k(:,ng1+1+m:ngluon),dim=2) + p(:,2)+ p(:,3)
             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)

             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif

             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif

             res = res + tmp

          enddo


          if (fl0 == 'bot') flaux = 'top'
          if (fl0 == 'top') flaux = 'bot'


          sp4 = fW_bffbf(e,k,sp,p,fll,flaux,eWm,kWm,ng1,ng2,ng3,1,giarray,qiarray,&
               &iWm,pol_int)
          k4 = sum(k(:,1:ngluon),dim=2) + sum(p(:,1:3),dim=2) + kWm
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4,k4)

          tmp = vbqW(sp4,eWp)

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          res = res+tmp

       elseif (sw == 2) then

          if (fl0 == 'bot') flaux = 'top'
          if (fl0 == 'top') flaux = 'bot'

          do m = 0,ng2

             sp4 = fW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),fl1,fl0,&
                  &eWm,kWm,ng1,giarray(1:ng1+m),qiarray(1:1),iWm,&
                  & pol_int)

             k4 = sum(k(:,1:ng1+m),dim=2)+ p(:,1) + kWm

             k4sq = sc(k4,k4) !- mass2

             sp4 = spb2(sp4,k4) !+ mass*sp4

             e2 = gW_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),sp(:,2),&
                  p(:,2),fll(2),sp(:,3),p(:,3),fll(3), eWp,kWp, &
                  & ng2-m, ng3, giarray(ng1+m+1:ngluon),qiarray(2:3),iWp,&
                  & pol_int)

             k2 = sum(k(:,ng1+m+1:ngluon),dim=2)+kWp+p(:,2)+p(:,3)

             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)

             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif

             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif

             res = res + tmp

          enddo

          sp4 = fW_bffbf(e,k,sp,p,fll,flaux,eWp,kWp,ng1,ng2,ng3,3,&
               & giarray, qiarray,iWp,pol_int)

          k4 = sum(k(:,1:ngluon),dim=2) + kWp + sum(p(:,1:3),dim=2)

          k4sq = sc(k4,k4) !- mass2

          sp4 = spb2(sp4,k4) !+ mass*sp4


          tmp = vbqW(sp4,eWm)

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          res = res + tmp

       elseif (sw ==3) then

          do m = 0,ng2 

             sp4 = f(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),fl1,fl0,ng1,&
                  &giarray(1:ng1+m),qiarray(1:1),pol_int)
             k4 = p(:,1)+sum(k(:,1:ng1+m),dim=2)

             k4sq = sc(k4,k4)

             if (m>0 .or. ng1>0) then
                sp4 = spb2(sp4,k4)
             endif

             e2 = gWW_fbf(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon),sp(:,2),p(:,2),&
                  & fl2,sp(:,3),p(:,3),fl3,eWp,kWp,eWm,kWm,ng2-m,ng3,giarray(ng1+m+1:+ngluon),&
                  & qiarray(2:3),iWp,iWm,pol_int)
             k2 = sum(k(:,ng1+m+1:ngluon),dim=2) + p(:,2)+p(:,3)+kWp+kWm
             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)


             if (ng1+m>=1) then
                if (abs(k4sq) > propcut) then
                   tmp = ci*tmp/k4sq
                else
                   tmp = czero
                endif
             endif

             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif

             res = res+tmp

          enddo


       endif


       ! For any value of sw

       do m = 1,ng1

          sp4 = fWW_bffbf(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,fll,fl0, &
               & eWp,kWp,eWm,kWm,ng1-m,ng2,ng3,sw,giarray(m+1:ngluon),qiarray,&
               & iWp,iWm,pol_int)

          k4 = sum(k(:,m+1:ngluon),dim=2) + p(:,1)+p(:,2)+p(:,3) + &
               & kWp + kWm

          k4sq = sc(k4,k4) !- mass2

          sp4 = spb2(sp4,k4) !+ mass*sp4

          e2 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)


          k2 = sum(k(:,1:m),dim=2)
          k2sq = sc(k2,k2)

          tmp = vgq(e2,sp4)


          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif


          if (m>1) then
             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif
          endif

          res = res + tmp

       enddo

       do m = 0,ng4-1

          ngL = ng1 + ng2 + ng3 + m

          sp4 = fWW_bffbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0, &
               & eWp,kWp,eWm,kWm,ng1,ng2,ng3,sw,giarray(1:ngL),qiarray,&
               & iWp,iWm,pol_int)


          k4 = sum(k(:,1:ngL),dim=2) + sum(p(:,1:3),dim=2) + &
               & kWp+kWm

          sp4 = spb2(sp4,k4) !+ mass*sp4

          k4sq = sc(k4,k4) !- mass2

          e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
               & giarray(ngL+1:ngluon),pol_int)

          k2 = sum(k(:,ngL+1:ngluon),dim=2)
          k2sq = sc(k2,k2)

          tmp = vqg(sp4,e2)


          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif


          if (m < ng4-1) then
             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif
          endif

          res = res + tmp

       enddo

    endif

    ! -- store current 
    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,iWp,iWm)


  end function fWW_bffbf_TM




  recursive function gWW_bff(e,k,sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
       &giarray,qiarray,WWid,pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:)
    type(mp_complex), intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    type(mp_complex), intent(in) :: eW(:,:),kW(:,:)
    integer, intent(in) ::  ng1,ng2
    character, intent(in) :: fl1*3,fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int 
    ! -----------------------------------------------------------------------
    integer             :: m1,m2, ms1a, ms2a,Dv!,m,m3
    integer :: ngluon,  ng3, ngL 
    type(mp_complex)             :: res(size(e,dim=1))
    type(mp_complex)             :: tmp(size(e,dim=1))
    type(mp_complex)             :: k1(size(p1))
    type(mp_complex)             :: k2(size(p1))
    type(mp_complex)             :: k3(size(p1))
    type(mp_complex)             :: k4(size(p1))
    type(mp_complex)             :: sp3(size(sp1))
    type(mp_complex)             :: sp4(size(sp1))
    type(mp_complex)             :: e1(size(e,dim=1))
    type(mp_complex)             :: e2(size(e,dim=1))
    type(mp_complex)             :: e3(size(e,dim=1))
    type(mp_complex)  :: k1sq,k2sq,k3sq,k4sq
    logical                   :: done 
    character:: flint*3



    if (fl1.eq.'str'.or.fl2.eq.'str') then 
       res = czero 
       return 
    endif

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'giarray missing'
    endif


    if (fl1 == 'top') then
       flint = 'bot'
    elseif (fl1 == 'bot') then
       flint = 'top'
    endif

    Dv = size(e,dim=1)

    ngluon = size(e,dim=2)
    ng3 = ngluon - ng1 - ng2


    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT D'

    if (ngluon == 0) then 

       res = czero

       sp4 = fWW(e,k,sp1,p1,fl1,fl2,eW,kW,0,giarray,qiarray(2:2),WWid,pol_int)

       k2 = p1 + kW(:,1)+kW(:,2)
       k2sq = sc(k2,k2) !-mass2
       sp4 = ci/k2sq*(spb2(sp4,k2))!+mass*sp4)
       res = -cone*vbqq(Dv,sp4,sp2)

       sp4 = bfWW(e,k,sp2,p2,fl2,fl1,eW,kW,0,giarray,qiarray(1:1),WWid,pol_int) 
       k2 = -p2 - kW(:,1)-kW(:,2)
       k2sq = sc(k2,k2) !- mass2
       sp4 = ci/k2sq*(spi2(k2,sp4))! + mass*sp4)
       tmp = -cone*vbqq(Dv,sp1,sp4)

       res = res + tmp

       sp4 = fW(e,k,sp1,p1,flint,fl1, eW(:,2), kW(:,2), 0, giarray, qiarray(2:2), WWid(2), pol_int)
       k2 = p1 + kW(:,2)
       sp4 = ci*(spb2(sp4,k2))/(sc(k2,k2))!-mass2)

       sp3 = bfW(e,k,sp2,p2,flint,fl2,eW(:,1), kW(:,1), 0, giarray, qiarray(1:1), WWid(1), pol_int)
       k2 = -p2-kW(:,1)
       sp3 = ci*(spi2(k2,sp3))/(sc(k2,k2))!-mass2)

       tmp = -cone*vbqq(Dv, sp4,sp3)

       res = res+tmp
    else
       write(*,*) 'function gWW_bff not written for ngluons /=0'
    endif

  end function gWW_bff

  !----------------

  !--------------------------------

  !! Used for the case where there are strange quarks and dummy lines.
  recursive function fWW_bffbf_1(e,k,sp,p, fll, fl0, eW, kW, ng1, ng2, ng3,sw,giarray,qiarray,WWid, pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:),eW(:,:), kW(:,:)
    type(mp_complex), intent(in) :: sp(:,:), p(:,:)
    integer, intent(in) ::  ng1,ng2, ng3,sw
    character, intent(in) :: fll(:)*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int
    character:: fl1*3,fl2*3,fl3*3, fl0*3, flaux*3
    integer::ngluons,ng4
    type(mp_complex)::res(size(sp,dim=1)), e2(size(e,dim=1)), k2(size(k,dim=1)),k2sq
    type(mp_complex)::tmp(size(sp,dim=1)), e1(size(e,dim=1))
    type(mp_complex)::sp1(size(sp,dim=1)), k1(size(k,dim=1)),k1sq
    type(mp_complex)::sp4(size(sp,dim=1)), k4(size(k,dim=1)),k4sq
    integer, parameter::Ndumm = 0
    type(mp_complex)             :: kdumm(size(k,dim=1),Ndumm)
    type(mp_complex)             :: edumm(size(e,dim=1),Ndumm)


    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)

    ngluons = size(e,dim=2)
    ng4 = ngluons - ng1 - ng2-ng3
    res = czero
    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C: fWW_bffbf_1', ngluons, ng1, ng2, ng3

    if (ngluons == 0) then
       if (fl0 == 'top')  flaux = 'bot'
       if (fl0 == 'bot')  flaux = 'top'

       if ((sw == 1) .or. (sw==3)) then

          sp1 = fW_bffbf_2(e,k,sp,p,fll,flaux,eW(:,2),kW(:,2),ng1,ng2,ng3,sw,&
               &giarray,qiarray,WWid(2),pol_int)

          k1 = p(:,1)+p(:,2)+p(:,3)+kW(:,2)
          k1sq = sc(k1,k1)

          sp1 =spb2(sp1,k1)!+mass*sp1

          tmp = vbqW(sp1,eW(:,1))

          if (abs(k1sq) > propcut) then
             tmp = ci*tmp/k1sq
          else
             tmp = czero
          endif

          res = res+tmp

       endif


       if (sw == 2) then

          e1 = gWW_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,eW,kW,ng1,ng2,&
               &giarray,qiarray(1:2),WWid,pol_int)

          k1 = p(:,1)+p(:,2)+kW(:,1)+kW(:,2)
          k1sq = sc(k1,k1)


          tmp = vqg(sp(:,3),e1)

          if (abs(k1sq) >  propcut) then
             tmp = -ci*tmp/k1sq
          else
             tmp = czero
          endif

          res = res+tmp
       endif

       if (sw == 4) then

          !--- Both W's above gluon

          sp4  = fWW(edumm, kdumm, sp(:,3),p(:,3), fl3, fl0,eW, kW, 0, giarray,&
               &qiarray(3:3), WWid, pol_int)
          k4 = p(:,3)+kW(:,1)+kW(:,2)
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4, k4)!+sp4*mass

          e2 = g_bff(edumm, kdumm, sp(:,1), p(:,1), fl1, sp(:,2),p(:,2),fl2,0,0,&
               &giarray,qiarray(1:2),pol_int) 
          k2 = p(:,1)+p(:,2)
          k2sq = sc(k2,k2)

          tmp = vgq(e2, sp4)


          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          res = res+tmp   !  sw = 4, #1

          ! ---- 1 W below gluon, 1 above, and both below

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'

          sp4 = fW_bffbf_1(edumm,kdumm,sp,p,fll,flaux,eW(:,2), kW(:,2),&
               &ng1,ng2,ng3,giarray,qiarray,WWid(2),pol_int)    


          k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,2)
          sp4 = spb2(sp4,k4) !+ mass*sp4

          k4sq = sc(k4,k4)

          tmp = vbqW(sp4,eW(:,1))

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp
          else 
             tmp = czero
          endif

          res = res + tmp   ! sw=4,#2

       endif
    else
       write(*,*) 'function fWW_bffbf_1 not written for ngluons /=0!'
    endif


  end function fWW_bffbf_1

  !--------------------


!!$!! Only done for ngluons = 0 and sw = 2 or 4.
  recursive function bfWW_fbff(e,k,sp,p,fll,fl0,eW,kW,ng1,ng2,ng3,sw,giarray,&
       &qiarray,WWid, pol_int) result(res)
    type(mp_complex), intent(in) :: e(:,:), k(:,:),eW(:,:), kW(:,:)
    type(mp_complex), intent(in) :: sp(:,:), p(:,:)
    integer, intent(in) ::  ng1,ng2, ng3,sw
    character, intent(in) :: fll(:)*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int

    character:: fl1*3,fl2*3,fl3*3, fl0*3, flaux*3
    integer::ngluons,ng4
    type(mp_complex)::res(size(sp,dim=1)), e2(size(e,dim=1)), k2(size(k,dim=1)),k2sq
    type(mp_complex)::tmp(size(sp,dim=1))
    type(mp_complex)::e3(size(e,dim=1)), k3(size(k,dim=1)),k3sq, k4(size(k,dim=1)),k4sq
    type(mp_complex)::sp1(size(sp,dim=1)),sp2(size(sp,dim=1)), sp3(size(sp,dim=1))
    type(mp_complex)::k1(size(k,dim=1)), k1sq

    integer, parameter::Ndumm = 0
    type(mp_complex)             :: kdumm(size(k,dim=1),Ndumm)
    type(mp_complex)             :: edumm(size(e,dim=1),Ndumm)
    type(mp_complex):: eww(size(ew,dim=1),size(ew,dim=2))
    type(mp_complex)::kww(size(kw,dim=1),size(kw,dim=2))

    !    mass = mt

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)

    ngluons = size(e,dim=2)
    ng4 = ngluons - ng1 - ng2-ng3
    res = czero
    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C: bfWW_fbff', ngluons, ng1, ng2, ng3, sw
    if (ngluons == 0) then
       res = czero

       ! Both gluons attached to the first quark pair.
       if (sw == 2) then

          e2 = gWW_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,&
               &eW,kW,0,0,giarray,qiarray(1:2),WWid,pol_int)
          k2 = p(:,1)+p(:,2)+sum(kW,dim=2)
          k2sq = sc(k2,k2)

          if (abs(k2sq) > propcut.and.fl0==fl1) then 
             tmp = -ci/k2sq*vgbq(e2,sp(:,3))
          else 
             tmp = czero 
          endif

          res = res + tmp  ! #1
       endif

       ! Both quarks attached to the second quark pair.
       if (sw == 4) then

          kww(:,1)=kw(:,2)
          kww(:,2)=kw(:,1)
          eww(:,1)=ew(:,2)
          eww(:,2) = ew(:,1)

          sp3 = bfWW(edumm, kdumm, sp(:,3),p(:,3), fll(3), fl0, eW, kW, 0, &
               &giarray, qiarray(3:3), WWid, pol_int)
          k3 = -p(:,3)-kW(:,1)-kW(:,2)
          k3sq = sc(k3,k3)

          sp3 = spi2(k3,sp3)!+mass*sp3

          e3 = g_fbf(edumm, kdumm, sp(:,1),p(:,1), fll(1),sp(:,2),p(:,2), fll(2),&
               &0,0,giarray, qiarray(1:2),pol_int)
          k4 = p(:,1)+p(:,2)
          k4sq = sc(k4,k4)

          tmp = vgbq(e3,sp3)

          if (abs(k3sq) > propcut) then
             tmp = tmp*ci/k3sq
          else
             tmp = czero
          endif

          if (abs(k4sq)>propcut) then
             tmp = -ci*tmp/k4sq
          else
             tmp = czero
          endif

          res = res+tmp                ! sw = 4, #1

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'

          sp2 = bfW_fbff(edumm, kdumm, sp, p, fll, flaux, eW(:,1), kW(:,1), 0,0,&
               &0,2,giarray,qiarray,WWid(1),pol_int)
          k2 = -p(:,1)-p(:,2)-p(:,3)-kW(:,1)
          k2sq = sc(k2,k2)

          sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = vWq(eW(:,2),sp2)
          if (abs(k2sq)>propcut) then
             tmp = ci*tmp/k2sq
          else
             tmp = czero
          endif

          res=res+tmp                  ! sw=4, #2    

       endif


       ! One quark attached to each of the two quark lines.
       if (sw ==3) then

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'




          sp2 = bfW(edumm, kdumm, sp(:,3),p(:,3),fll(3),fl0,eW(:,1),kW(:,1),&
               &0,giarray,qiarray(3:3), WWid(1), pol_int)
          k2 = -p(:,3)-kW(:,1)
          k2sq = sc(k2,k2)

          sp2 = spi2(k2,sp2)!+mass*sp2

          e2 = gW_fbf(edumm,kdumm, sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2, &
               &eW(:,2),kW(:,2),0,0, giarray,qiarray(1:2),WWid(2),pol_int)     
          k3 = p(:,1)+p(:,2)+kW(:,2)
          k3sq = sc(k3,k3)

          tmp = vgbq(e2,sp2)

          if (abs(k2sq)>propcut) then
             tmp=ci*tmp/k2sq
          else
             tmp = czero
          endif

          if (abs(k3sq)>propcut) then
             tmp = -ci*tmp/k3sq
          else
             tmp = czero
          endif

          res = res +tmp                    !sw=3, #1


          sp1 = bfW_fbff(edumm, kdumm, sp, p, fll, flaux, eW(:,2),kW(:,2),0,0,&
               &0,2, giarray, qiarray,WWid(2), pol_int)
          k1 = -p(:,1)-p(:,2)-p(:,3)-kW(:,2)
          k1sq = sc(k1,k1)

          sp1 = spi2(k1,sp1)!+mass*sp1

          tmp = vWq(eW(:,1),sp1)

          if (abs(k1sq) > propcut) then
             tmp = ci*tmp/k1sq
          else
             tmp = czero
          endif

          !    res = res+tmp           !sw=3, #2


       endif

    else              ! If ngluons > 0

       write(*,*) 'function bfWW_fbff only written for ngluons = 0'
    endif


  end function bfWW_fbff

  recursive function gWW_sbsfbf(e,k,sp,p,fll,eW,kW,ng1,ng2,ng3,ng4,&
       &sw,giarray,qiarray,WWid, pol_int) result(res)

    type(mp_complex), intent(in) :: e(:,:), k(:,:),eW(:,:), kW(:,:)
    type(mp_complex), intent(in) :: sp(:,:), p(:,:)
    integer, intent(in) ::  ng1,ng2, ng3,ng4,sw
    character, intent(in) :: fll(:)*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int 
    ! -----------------------------------------------------
    type(mp_complex)::sp1(size(sp,dim=1)),sp2(size(sp,dim=1)),sp3(size(sp,dim=1))
    type(mp_complex)::sp4(size(sp,dim=1)),sp5(size(sp,dim=1)),sp6(size(sp,dim=1))
    type(mp_complex)::sp7(size(sp,dim=1)),e1(size(e,dim=1)),e2(size(e,dim=1))
    type(mp_complex)::k1(size(k,dim=1)),k2(size(k,dim=1)),k3(size(k,dim=1))
    type(mp_complex)::k4(size(k,dim=1)),k5(size(k,dim=1)),k6(size(k,dim=1))
    type(mp_complex)::k7(size(k,dim=1))
    type(mp_complex)::k1sq,k2sq,k3sq,k4sq,k5sq,k6sq
    integer::Dv,ngluon,ng5
    type(mp_complex) :: tmp(size(e,dim=1)),res(size(e,dim=1))
    character::fl1*3,fl2*3,fl3*3,fl4*3,flaux*3
    type(mp_complex)::edumm(4,0),kdumm(4,0)
    logical :: done 

    done = .false. 
    if (present(giarray)) then 
       call memory_check(pol_int,res,done,giarray,qiarray,WWid(1),WWid(2))
       if (done) return 
    else
       write(*,*) 'giarray missing'
    endif

    Dv = size(e,dim=1)
    ngluon = size(e,dim=2)

    if ((sw.ne.2).and.(sw.ne.4)) print *, 'W position incorrect'

    ng5 = ngluon - ng1 - ng2-ng3 - ng4
    if (ng5 <0 ) write(*,*) 'no of gluons incorrect'


    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)
    fl4 = fll(4)


    if (ngluon == 0) then

       if ((sw == 2) .or. (sw==4)) then
          sp1 = fWW_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),fll(2:4),&
               &fl1,eW(:,1),kW(:,1),eW(:,2),kW(:,2),0,0,0,sw-1,giarray,&
               &qiarray(2:4),WWid(1),WWid(2),pol_int)

          k1 = p(:,2)+p(:,3)+p(:,4)+sum(kW,dim=2)
          k1sq = sc(k1,k1)
          sp1 = spb2(sp1,k1)!+mass*sp1

          tmp = -(cone*vbqq(Dv,sp1,sp(:,1)))

          if (abs(k1sq) > propcut) then
             tmp = tmp*ci/k1sq
          else
             tmp = czero
          endif

          res = res + tmp

          sp2 = bfWW_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),fll(1:3),&
               &fl1,eW,kW,0,0,0,sw,giarray,qiarray(1:3),WWid,pol_int)

          k2 = -p(:,1)-p(:,2)-p(:,3) - sum(kW,dim=2)
          k2sq = sc(k2,k2)
          sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = -cone*vbqq(Dv,sp1,sp2)

          if (abs(k2sq) > propcut) then
             tmp = tmp*ci/k1sq
          else
             tmp = czero
          endif
          res = res + tmp       ! # 3

       endif

       if (sw == 2) then

          if (fl1 .eq. 'top') then
             flaux = 'bot'
          elseif (fl1 .eq. 'bot') then
             flaux = 'top'
          endif

          sp3 = fW_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
               &fll(2:4),flaux,eW(:,2),kW(:,2),0,0,0,1,giarray,&
               &qiarray(2:4),WWid(2),pol_int)

          k3 = p(:,2)+p(:,3)+p(:,4)+kW(:,2)

          sp3 = spb2(sp2,k3)!+mass*sp3
          k3sq = sc(k3,k3)

          sp4 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl2,&
               &eW(:,1),kW(:,1),0,giarray,qiarray(1:1),WWid(1),pol_int)
          k4 = -p(:,1)-kW(:,1)
          sp4 = spi2(k1,sp4) !+ mass*sp4
          k4sq = sc(k4,k4)

          tmp = -cone*vbqq(Dv,sp3,sp4)

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          if (abs(k3sq) >propcut) then
             tmp = ci*tmp/k3sq
          else
             tmp = czero
          endif

          res = res + tmp 

          !------

          sp5 = f_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),fll(2:4),&
               &fl2,0,0,0,giarray,qiarray(2:4),pol_int)

          k5 = p(:,4)+p(:,2)+p(:,3)
          k5sq = sc(k5,k5)
          sp5 = spb2(sp5,k5)!+mass*sp5

          sp6 = bfWW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl2,eW,kW,0,&
               &giarray,qiarray(1:1),WWid,pol_int)

          k6 = -p(:,1) - sum(kW,dim=2)
          k6sq = sc(k6,k6)
          sp6 = spi2(k5,sp6)!+mass*sp6


          tmp = -cone*vbqq(Dv,sp5,sp6)

          if (abs(k5sq) > propcut) then
             tmp = ci*tmp/k5sq
          else
             tmp = czero
          endif

          if (abs(k6sq) >propcut) then
             tmp = ci*tmp/k6sq
          else
             tmp = czero
          endif

          res = res+tmp

          !------------

          e1 = gWW_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),&
               &p(:,2),fl2,eW,kW,0,0,&
               &giarray,qiarray(1:2),WWid,pol_int)
          k1 = p(:,1) + p(:,2)+sum(kW,dim=2)
          k1sq = sc(k1,k1)

          e2 = g_fbf(edumm,kdumm,sp(:,3),p(:,3),fl3,sp(:,4),p(:,4),fl4,0,0,&
               &giarray,qiarray(3:4),pol_int)
          k2 = p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          tmp = vggg(e1,k1,e2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = -ci*tmp/k1sq
          else
             tmp = czero
          endif

          if (abs(k2sq) > propcut) then
             tmp = -ci*tmp/k2sq
          else
             tmp = czero
          endif


          res = res+tmp
       endif

       !------------

       if (sw == 4) then

          if (fl4 .eq. 'top') then
             flaux = 'bot'
          elseif (fl4 .eq. 'bot') then
             flaux = 'top'
          endif

          sp3 = bfW_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),fll(1:3),&
               &flaux,eW(:,1),kW(:,1),0,0,0,4,giarray,qiarray(1:3),&
               &WWid(1),pol_int)

          k3 = -p(:,1)-p(:,2)-p(:,3)-kW(:,1)
          k3sq = sc(k3,k3)
          sp3 = spi2(k3,sp3)!+mass*sp3


          sp4 = fW(edumm,kdumm,sp(:,4),p(:,4),fl4,fl3,eW(:,2),&
               &kW(:,2),0,giarray,qiarray(4:4),WWid(2),pol_int)

          k4 = p(:,4)+kW(:,2)
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)


          if (abs(k3sq) > propcut) then
             tmp = ci*tmp/k3sq
          else
             tmp = czero
          endif

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif
          res = res+tmp

          !---------

          sp5 = bf_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),fll(1:3),&
               &fl3,0,0,0,giarray,qiarray(1:3),pol_int)

          k5 = -p(:,1)-p(:,2)-p(:,3)
          k5sq = sc(k5,k5)
          sp5 = spb2(sp5,k5)!+mass*sp5

          sp6 = fWW(edumm,kdumm, sp(:,4),p(:,4),fl4,fl3,eW,kW,0,&
               &giarray,qiarray(4:4),WWid,pol_int)

          k6 = p(:,4)+sum(kW,dim=2)
          k6sq = sc(k6,k6)
          sp6 = spb2(sp6,k6)!+mass*sp6

          tmp = -cone*vbqq(Dv,sp6,sp5)
          if (abs(k5sq) > propcut) then
             tmp = ci*tmp/k5sq
          else
             tmp = czero
          endif

          if (abs(k6sq) > propcut) then
             tmp = ci*tmp/k6sq
          else
             tmp = czero
          endif

          res = res+tmp


          !------

          e1 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,0,0,&
               &giarray,qiarray(1:2),pol_int)
          k1 = p(:,1) + p(:,2)+sum(kW,dim=2)
          k1sq = sc(k1,k1)

          e2 = gWW_fbf(edumm,kdumm,sp(:,3),p(:,3),fl3,sp(:,4),&
               &p(:,4),fl4,eW,kW,0,0,&
               &giarray,qiarray(3:4),WWid,pol_int)
          k2 = p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          tmp = vggg(e1,k1,e2,k2)

          if (abs(k1sq) > propcut) then
             tmp = -ci*tmp/k1sq
          else
             tmp = czero
          endif

          if (abs(k2sq) > propcut) then
             tmp = -ci*tmp/k2sq
          else
             tmp = czero
          endif

          res = res+tmp


       endif

    elseif (ngluon > 0) then
       write(*,*) 'No recursion relations for ngluons >0'
    endif


  end function gww_sbsfbf



end module mprecurrencebitsfour


