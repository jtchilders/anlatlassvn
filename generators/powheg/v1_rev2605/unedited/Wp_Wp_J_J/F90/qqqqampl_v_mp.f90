
module qqqqampl_v_mp


  !------------------------------------------------------------------------!
  ! input momenta p has form						  !
  ! 0 ->  proton1(p1) + proton2(p2) + Wmp+ Wm + parton1(p7) + parton2(p8)  ! 
  !    	 	       		     |	  |    		     		  !
  !			   v(p3) l~(p4)	 l(p5) v~(p6)       		  !
  !------------------------------------------------------------------------!

  use mpmodule; 
  use mpsimpleoperations_c; use mpadvancedoperations_c
  use consts_mp 
  use mpconverter 


  use types
  use consts_MCFM
  use mprecurrence
  use mpaux_functions
  use mpmemory
  use define_ampl
  use mpauxiliary_functions
  use mpinitWpWp
  use comb
  use mpglobal
  use mpamplitude
  use match1
  use mpopp
  use masters
  use mpspinors
  use sub_defs_io 
  use mptodp_global 
  use mpprecision 

  implicit none 
  private 

  public :: redo_ampl_in_mp  

contains


  subroutine redo_ampl_in_mp(pext_dp,muu,i1,i2,i3,i4,css,tree,res)
    complex(dp), intent(in)   :: pext_dp(:,:) 
    real(dp), intent(in)      :: muu  
    integer, intent(in)       :: i1,i2,i3,i4,css   
    complex(dp), intent(out)  :: tree, res(-2:1)  
    !-------------------    
    type(mp_complex) :: pext_mp(size(pext_dp,dim=1),size(pext_dp,dim=2))
    type(mp_real) :: muu_mp 
    logical, save :: firsttime = .true. 
    integer :: init_prec 

    if (firsttime) then 
       ! -- MP initialization 
       init_prec = int_val_opt('-init_prec',20)
       call mpinit
       call SetPrecision(init_prec) 
       call fixmpconsts(init_prec) 
       call init_spinors
       call fix_bvec 
       firsttime = .false. 
    end if

    call generate_mless_mpevent(pext_dp,pext_mp) 
    muu_mp = muu 
    tree = 0d0  
    res = 0d0  

    call getamplqqqq_v_mp(pext_mp, muu_mp, i1,i2,i3,i4, css, tree, res)


  end subroutine redo_ampl_in_mp


  ! given a dpevent generate an mpevent 
  ! make sure that the Npart mass conditions and the 4 energy-momentum
  ! conservation conditions are fulfilled in multiple precision (not
  ! only up to 15 digits)
  subroutine generate_mless_mpevent(pext_dp,pext_mp) 
    use consts_mp 
    complex(dp), intent(in)       :: pext_dp(:,:)
    type(mp_complex), intent(out)      :: pext_mp(:,:) 
    ! --------------------------------------------
    type(mp_complex) :: pout(4),ptmp(4) 
    type(mp_complex) :: pin(4),po(4),pinp(4),pop(4)
    type(mp_complex) :: Ei,Eip,pzi,pzip,Eo,Eop,pzo,pzop,ptop2
    integer     :: Npart, i, ninc, iinc(2),abeam(2),j,iout,ipt 


    Npart = size(pext_dp,dim=1)
    if (Npart < 4) stop 'generate_mpevent: Npart<4'


    ! -- get the incoming partons and the direction of the beam (x,y,z ?) 
    ninc = 0 
    do i=1,Npart
       do j=2,4
          if (abs(abs(pext_dp(i,j))-abs(pext_dp(i,1))) < sq2tol) then 
             ninc = ninc+1 
             if (ninc >2) stop 'generate_massless_mpevent: more than 2 incoming?'
             iinc(ninc) = i 
             abeam(ninc) = j 
          end if
       end do
    end do
    if (abeam(1) /= abeam(2)) stop 'generate_massless_mpevent: different abeams?' 

    pext_mp = czero 
    do i=1,Npart
       if (i /= iinc(1) .and. i /= iinc(2)) then 
          pext_mp(i,2) = pext_dp(i,2) 
          pext_mp(i,3) = pext_dp(i,3) 
          pext_mp(i,4) = pext_dp(i,4) 
          pext_mp(i,1) = sqrt(pext_mp(i,2)**2+pext_mp(i,3)**2+pext_mp(i,4)**2)
       endif
    enddo
    do i=1,4 
       pout(i) = sum(pext_mp(1:Npart,i))
    end do

    ! after that have 2 * 4 components yet to be fixed and 6 more
    ! conditions to be fullfilled

    ! set iinc(1) along beam axis (fix 2 components) 
    pext_mp(iinc(1),2:4) = czero 
    do j=2,4
       if (j /= abeam(1)) then 
          pext_mp(iinc(2),j) = -pout(j) ! iinc(2) has tiny pt components now... 
       end if
    end do

    pext_mp(iinc(1),1) = (-pout(1)**2+pout(4)**2)
    do j=2,4 
       if (j /= abeam(1)) then 
          pext_mp(iinc(1),1) = pext_mp(iinc(1),1)+pext_mp(iinc(2),j)**2
       endif
    enddo
    pext_mp(iinc(1),1) = pext_mp(iinc(1),1)/(two*(pout(1)-pout(abeam(1))))

    pext_mp(iinc(1),abeam(1)) = pext_mp(iinc(1),1)


    pext_mp(iinc(2),1) = (-pout(1)**2+two*pout(1)*pout(abeam(1))-pout(abeam(1))**2)
    do j=2,4 
       if (j /= abeam(1)) then 
          pext_mp(iinc(2),1) = pext_mp(iinc(2),1)-pext_mp(iinc(2),j)**2
       endif
    enddo
    pext_mp(iinc(2),1) = pext_mp(iinc(2),1)/(two*(pout(1)-pout(abeam(1))))

    pext_mp(iinc(2),abeam(2)) = (pout(1)**2-two*pout(1)*pout(abeam(1))+pout(abeam(1))**2)

    do j=2,4 
       if (j /= abeam(1)) then 
          pext_mp(iinc(2),abeam(2)) = pext_mp(iinc(2),abeam(2))-pext_mp(iinc(2),j)**2
       endif
    enddo
    pext_mp(iinc(2),abeam(2)) = pext_mp(iinc(2),abeam(2))/(two*(pout(1)-pout(abeam(1))))

    ipt = iinc(2) 
    if (real((pext_mp(iinc(1),abeam(1)))*pext_dp(iinc(1),abeam(1))) < zero) then 
       ! - need to swap incoming ones 
       ptmp =  pext_mp(iinc(1),:)
       pext_mp(iinc(1),:) = pext_mp(iinc(2),:)
       pext_mp(iinc(2),:) = ptmp 
       ipt = iinc(1) 
    endif

    ! -- adjust kinematics -- let one outgoing parton absorb the pt components
    pin = pext_mp(ipt,:) 
    do i=1,Npart 
       if (i /= iinc(1) .and. i /= iinc(2)) then 
          po = pext_mp(i,:) 
          iout = i 
          exit 
       end if
    enddo

    pinp = czero 
    pop = czero 
    Ei = pin(1)
    Eo = po(1)
    pzi = pin(abeam(1))
    pzo = po(abeam(1))
    ptop2 = czero 
    do i=2,4
       if (i /= abeam(1)) then 
          ptop2 = ptop2 + (po(i)+pin(i))**2
       endif
    enddo


    Eip = (Ei**2 + 2*Ei*Eo + Eo**2 - ptop2 - pzi**2 - 2*pzi*pzo - pzo**2)/(2*(Ei + Eo + pzi + pzo))
    pzip = (-Ei**2 - 2*Ei*Eo - Eo**2 + ptop2 + pzi**2 + 2*pzi*pzo + pzo**2)/(2*(Ei + Eo + pzi + pzo))
    if (real(Eip*Ei) > zero .and. real(pzip*pzi)> zero) then 
       Eop = (Ei**2 + 2*Ei*Eo + Eo**2 + ptop2 + 2*Ei*pzi + 2*Eo*pzi + pzi**2 + 2*Ei*pzo + 2*Eo*pzo + &
            &2*pzi*pzo + pzo**2)/(2*(Ei + Eo + pzi + pzo))
       pzop = (Ei**2 + 2*Ei*Eo + Eo**2 - ptop2 + 2*Ei*pzi + 2*Eo*pzi + pzi**2 + 2*Ei*pzo + 2*Eo*pzo + &
            &2*pzi*pzo + pzo**2)/(2*(Ei + Eo + pzi + pzo))
    else
       Eip =  (-Ei**2 - 2*Ei*Eo - Eo**2 + ptop2 + pzi**2 + 2*pzi*pzo + pzo**2)/(2*(-Ei - Eo + pzi + pzo))
       pzip = (-Ei**2 - 2*Ei*Eo - Eo**2 + ptop2 + pzi**2 + 2*pzi*pzo + pzo**2)/(2*(-Ei - Eo + pzi + pzo)) 
       Eop = (Ei**2 + 2*Ei*Eo + Eo**2 + ptop2 - 2*Ei*pzi - 2*Eo*pzi + pzi**2 - 2*Ei*pzo - 2*Eo*pzo &
            &+ 2*pzi*pzo + pzo**2)/(2*(Ei + Eo - pzi - pzo))
       pzop = (-Ei**2 - 2*Ei*Eo - Eo**2 + ptop2 + 2*Ei*pzi + 2*Eo*pzi - pzi**2 + 2*Ei*pzo + 2*Eo*pzo &
            &- 2*pzi*pzo - pzo**2)/(2*(Ei + Eo - pzi - pzo))

       if (real(Eip*Ei) < zero .or. real(pzip*pzi) < zero) write(*,*) 'could not find right solution?'
    end if
    pinp(1) = Eip
    pinp(abeam(1)) = pzip
    pop(1) = Eop
    pop(abeam(1)) = pzop 
    do i=2,4
       if (i /= abeam(1)) then 
          pop(i) = po(i)+pin(i)
       endif
    enddo

    pext_mp(ipt,:) = pinp 
    pext_mp(iout,:) = pop 

  end subroutine generate_mless_mpevent

  !! Subroutine to compute primitive amplitudes.
  !! Inputs: momenta p(8,4), quark momentum ordering, scale muu, 
  !! css to decide which primitive to use, polesonly.
  !! Helicities fixed: do not need to be inputted.
  !! Outputs: tree level amplitude, double and single poles, and rational
  !! and cut-constructible part at 1L.

  !! Momenta conventions:
  !! 0 -> qb(i1) + q(i2) + W^+          + W^+  + qb(i3) + q(i4)
  !!                     |              |
  !!                     V              V
  !!        ( -> ve(3) + e^+(4))   (vm(5) + mu^+(6)

  subroutine getamplqqqq_v_mp(pin, muu, i1,i2,i3,i4, css, tree, res)
    type(mp_complex), intent(in)      :: pin(:,:)
    type(mp_real), intent(in)             :: muu
    integer, intent(in)          :: i1,i2,i3,i4
    integer, intent(in)          :: css
    complex(dp), intent(out)         :: tree, res(-2:1)
    !----------------------------------------------------
    type(mp_complex)                  :: div1, div2, div3
    complex(dp)                  :: resWW(-2:1)
    type(mp_complex)                  :: treeWW
    type(mp_complex)                  :: pp(6,4), pl(2,4), pa(2,4)
    type(mp_complex)                  :: pps(6,4), pls(2,4), pas(2,4)
    type(mp_real)                     :: mu, s
    type(mp_complex)                  :: hardscale
    logical, save                :: first_time = .true.
    integer                      :: Npart, lncuts(5), N5, N4, N3, xN2, N2, N1
    integer                      :: hl(size(pin,dim=1)),i,j,Nloop
    logical                      :: pole_check = .false.
    logical                      :: crossed_gl
    logical                      :: outgoing(size(pin,dim=1))


    res = 0d0
    tree = 0d0
    div1 = czero
    div2 = czero
    outgoing = .true.
    crossed_gl = .false.

    ! Momenta of quarks and W's.
    pp(1,:) = pin(i1,:)              
    pp(2,:) = pin(i2,:)          
    pp(3,:) = pin(3,:)+pin(4,:)  
    pp(4,:) = pin(5,:)+pin(6,:)  
    pp(5,:) = pin(i3,:)          
    pp(6,:) = pin(i4,:)          


    ! Lepton momenta.
    pl(1,:) = pin(3,:)      
    pa(1,:) = pin(4,:)      
    pl(2,:) = pin(5,:)      
    pa(2,:) = pin(6,:) 

    mu = muu 

    s = sc(pin(1,:)+pin(2,:),pin(1,:)+pin(2,:))
    hardscale = sqrt(abs(s))*2d0

    pp = pp/hardscale
    pl = pl/hardscale
    pa = pa/hardscale
    mu = muu/hardscale  


    Npart = Npoint + 2
    Nterm = 4
    Ncut = 1

    if (first_time) then
       call init_spinors
       call fix_bvec
       call allocate_mem
       call allocate_mominfo(Ncut, Npoint)
       call allocate_mominfo2(Ncut, Npoint)
       WWqqqq = .true.
       qbq_WW_and_gluons = .true.
       cashing = .false. 
    endif

    if (css == 1) then
       case_b1 = .true.
       case_b2 = .false.
       Nloop = 1
       crossed_gl = .false.
    elseif (css ==2 ) then
       case_b1 = .false.
       case_b2 = .true.
       Nloop = 2
       crossed_gl = .false.
    elseif (css == 3) then
       case_b1 = .true.
       case_b2 = .false.
       crossed_gl = .true.
       Nloop = 1
    else

       stop 'getamplqqqq_v : Not a valid input for css'
    endif


    call getcutnumb(lncuts, Npoint)

    call pampl_count_WpWp(lncuts, Npoint, N5, N4, N3, xN2, N2, N1)


    if (first_time) then         
       call allocate_coeffs(N5max,N4max,N3max,N2max,N1max)        
       call allocate_arrays(N5max,N4max,N3max,N2max)   
       first_time = .false.      
    endif

    ! This loop covers the swapping of the W's.
    do i = 1,2
       hl = (/1,-1,-1,-1,1,-1/)   


       ! The primitive B2 has the W only to the left of the gluon line. To get a result
       ! for the W on the right, the quark momenta and helicities are swapped in this 
       ! loop.
       ! Note that there is no swapping for case B1, as the recursion relations in this 
       ! case take care of all possible positions of the W's.


       do j = 1,Nloop

          call pamplWpWp(pp,pl,pa,hl,lncuts,N5,N4,N3,xN2,N2,N1,Lc5,F5,Yc5,&
               Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2,treeWW,mu,div2,div1,div3,outgoing,crossed_gl)


          call pentcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:))    

          call quadcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
               & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), &
               & F2(:N2,:),Yc2(:N2,:))

          call tripcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
               & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), & 
               & F2(:N2,:),Yc2(:N2,:))                                    

          call doubcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
               & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), &
               & F2(:N2,:),Yc2(:N2,:))                                   

          ! now switch to DP ned to allocate DO global and make transfer MP info in DP 
          ! then mastersubs uses DP info 
          call mptodp_glob(tagdcut,mom,hel,momline,Lab_ex, &
               &Lab_in,coeff5,propv5,coeff4,refvect4,propv4,&
               &coeff3,refvect3,propv3,coeff2,refvect2,propv2,coeff1,&
               &dcoeff5,dcoeff4,dcoeff3,dcoeff2,dcoeff1,&
               &mass5,mass4,mass3,mass2,mass1,.true.)

          call MasterSubs(mptodp(mu),Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:), &
               & F4(:N4,:),Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:), &
               & Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:),resWW)

          resWW = resWW*(0d0,1d0)

          res(:) = res(:) + resWW(:)
          tree = tree + treeWW


          ! helicity change:
          hl(1) = -hl(1)
          hl(2) = -hl(2)
          hl(5) = -hl(5)
          hl(6) = -hl(6)

          ! quark momentum swap...
          pp(1,:) = pin(i2,:)/hardscale
          pp(2,:) = pin(i1,:)/hardscale
          pp(5,:) = pin(i4,:)/hardscale
          pp(6,:) = pin(i3,:)/hardscale

       enddo

       ! Set quark momenta back to normal...
       pp(1,:) = pin(i1,:)/hardscale
       pp(2,:) = pin(i2,:)/hardscale
       pp(5,:) = pin(i3,:)/hardscale
       pp(6,:) = pin(i4,:)/hardscale

       ! Swap W momenta.
       pp(4,:) = (pin(3,:)+pin(4,:))/hardscale      
       pp(3,:) = (pin(5,:)+pin(6,:))/hardscale     
       pl(1,:) = pin(5,:)/hardscale      
       pa(1,:) = pin(6,:)/hardscale      
       pl(2,:) = pin(3,:)/hardscale
       pa(2,:) = pin(4,:)/hardscale

    enddo

    res =  res/mptodp(hardscale**4)
    tree = tree/mptodp(hardscale**4)

    if (pole_check) then
       write(*,*) 'Tree amplitude ', tree
       write(*,*) ''
       write(*,*) 'Res ', res
       write(*,*) ''
       write(*,*) 'div2 ',  div2
       write(*,*) 'res(-2)/tree ', res(-2)/tree
       write(*,*) 'Relative difference of values: ', (res(-2)/tree - div2)/div2
       write(*,*) ''
       write(*,*) 'div1 ',  div1
       write(*,*) 'res(-1)/tree ', res(-1)/tree
       write(*,*) 'Relative difference of values: ', (res(-1)/tree - div1)/div1
       write(*,*) ''
    endif


  end subroutine getamplqqqq_v_mp




end module qqqqampl_v_mp
