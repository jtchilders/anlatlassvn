module qqqqampl_v

  use types
  use consts_dp
  use consts_MCFM
  use dprecurrence
  use dpmemory
  use dpglobal
  use comb
  use dpspinors
  use define_ampl
  use dpaux_functions
  use dpinitWpWp
  use dpopp
  use masters
  use sub_defs_io
  use dpstability 
  use qqqqampl_v_qp
!  use qqqqampl_v_mp
  use counter 

  implicit none

  private

  public  :: getamplqqqq_v, getamplqqqq_v_stu, mtotqqqq_v
  
  logical :: verbose = .false.
  logical :: verbose1 = .false.

contains
  
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
  
  subroutine getamplqqqq_v(pin, muu, i1,i2,i3,i4, css, tree, treesplit, res, polesonly)
    complex(dp), intent(inout)      :: pin(:,:)
    real(dp), intent(in)             :: muu
    integer, intent(in)          :: i1,i2,i3,i4
    integer, intent(in)          :: css
    logical, intent(in)          :: polesonly
    complex(dp), intent(out)         :: tree, res(-2:1), treesplit(2)
!----------------------------------------------------
    complex(dp)                  :: div1, div2, div3
    complex(dp)                  :: resWW(-2:1), treeWW
    complex(dp)                  :: pp(6,4), pl(2,4), pa(2,4)
    complex(dp)                  :: pps(6,4), pls(2,4), pas(2,4)
    real(dp)                     :: mu, s, hardscale
!    logical                      :: case_b1, case_b2
    logical, save                :: first_time = .true.
    integer                      :: Npart, lncuts(5), N5, N4, N3, xN2, N2, N1
    integer                      :: hl(size(pin,dim=1)),i,j,Nloop,ii
    logical                      :: pole_check = .false.
    logical                      :: crossed_gl
    logical                      :: outgoing(size(pin,dim=1))
!-----------------------------------------------------
    real(dp)    :: lnprec(-2:-1) ,lnprec_qp(-2:-1), abs_rel_err_dp 
    complex(dp) :: res_qp(-2:1), tree_qp 
    complex(dp) :: res_mp(-2:1), tree_mp 
    real(dp)                :: ln10 = 2.3025850929940456840_dp 


    if (verbose)    write(*,*) 'entering getamplqqqq_v', i1,i2,i3,i4,css
 

    res = czero
    tree = czero
    div1 = czero
    div2 = czero
    treesplit = czero
    outgoing = .true.
    crossed_gl = .false.

    ! Momenta of quarks and W's.
    pp(1,:) = pin(i1,:)              
    pp(2,:) = pin(i2,:)          
    pp(3,:) = pin(3,:)+pin(4,:)  
    pp(4,:) = pin(5,:)+pin(6,:)  
    pp(5,:) = pin(i3,:)          
    pp(6,:) = pin(i4,:)          
    do  i =1,8
       if (verbose)  write(*,*) 'pin', i,pin(i,:),sc(pin(i,:),pin(i,:))
    enddo
    if (verbose) write(*,*) 'psum', sum(pin,dim=1) 

   

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
   
    ! This is not used in the gg virt case -- done at higher level???

    ! A bit confused about first_time -- ask Tom.

    if (first_time) then
       call init_spinors
       call fix_bvec
       call allocate_mem
       call allocate_mominfo(Ncut, Npoint)
       call allocate_mominfo2(Ncut, Npoint)
       recomputed_in_qp = 0 
       total_ampl = 0 
    endif
       WWqqqq = .true.
       qbq_WW_and_gluons = .true.
       cashing = .false. 


    if (css == 1) then
       case_b1 = .true.
       case_b2 = .false.
       Nloop = 1
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

    if (verbose) write(*,*) 'Primitive ', css
    
    if (verbose) write(*,*) 'Calling getcutnumb',Ncut 
    
    call getcutnumb(lncuts, Npoint)
    
    if (verbose) write(*,*) 'Calling pampl_count_WW4q', Ncut
    
    call pampl_count_WpWp(lncuts, Npoint, N5, N4, N3, xN2, N2, N1)
    
    if (verbose) write(*,*) 'Number of cuts', N5, N4, N3, N2, N1
    
    if ((N5>N5max .or.  N4>N4max .or.N3>N3max .or.  N2>N2max .or. &
         & xN2>tN2max)) then         
       write(*,*) 'n2,n3,n4,n5',n2,n3,n4,n5,xn2        
       write(*,*) 'n2,n3,n4,n5',n2max,n3max,n4max,n5max,tn2max        
       stop 'getamplqqqq_v: something went wrong in setting of N5,N4,N3,N2'     
    endif

    if (first_time) then         
         call allocate_coeffs(N5max,N4max,N3max,N2max,N1max)        
         call allocate_arrays(N5max,N4max,N3max,N2max)   
         first_time = .false.      
    endif

    ! This loop covers the swapping of the W's.
    do i = 1,2
       if (verbose) write(*,*) 'iloop', i 
       hl = (/1,-1,-1,-1,1,-1/)   
       

      ! The primitive B2 has the W only to the left of the gluon line. To get a result
      ! for the W on the right, the quark momenta and helicities are swapped in this 
      ! loop.
      ! Note that there is no swapping for case B1, as the recursion relations in this 
      ! case take care of all possible positions of the W's.
       

       do j = 1,Nloop
          if (verbose) write(*,*) 'calling pamplWpWp', Ncut 
         call pamplWpWp(pp,pl,pa,hl,lncuts,N5,N4,N3,xN2,N2,N1,Lc5,F5,Yc5,&
              Lc4,F4,Yc4,Lc3,F3,Yc3,Lc2,F2,Yc2,treeWW,mu,div2,div1,div3,outgoing,crossed_gl)
          
         if (polesonly) then
             
            tree = tree + treeWW

            res(-2) = tree*div2
            res(-1) = tree*div1
            res(0) =  tree*div3
            res(1) = czero

         else
          
            if (verbose) write(*,*) 'doing pentcut'                  
          
            call pentcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:))    
          
            if (verbose ) write(*,*) 'doing quadcut'                  
          
            call quadcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
                 & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), &
                 & F2(:N2,:),Yc2(:N2,:))
          
            if (verbose) write(*,*) 'doing tripcut'  
             
            call tripcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
                 & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), & 
                 & F2(:N2,:),Yc2(:N2,:))                                    

            if (verbose) write(*,*) 'doing doubcut'                  
            
            call doubcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:), &
                 & Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:), &
                 & F2(:N2,:),Yc2(:N2,:))                                   
             
            if (verbose) write(*,*) 'calling Masters'                  

            call MasterSubs(mu,Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:), &
                 & F4(:N4,:),Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:), &
                 & Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:),resWW)

            abs_rel_err_dp = abs_rel_err
            
            resWW = ci*resWW
          
          if (verbose)  write(*,*) 'Tree init ', treeWW
          if (verbose)  write(*,*) 'res init', resWW

            res(:) = res(:) + resWW(:)
 ! Added RR 24-03-2011 -- for fermion loops
             treesplit(i) = treesplit(i) + treeWW
             tree = tree + treeWW
             

         endif
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
       
      if (verbose) write(*,*) 'swapped Ws',i
             
   enddo
   if (verbose) write(*,*) 'tree',tree
   if (verbose) write(*,*) 'res',res

   res =  res/hardscale**4
   tree = tree/hardscale**4
   treesplit = treesplit/hardscale**4
   total_ampl = total_ampl+1 
 
   if (verbose) write(*,*) 'tree',tree
   if (verbose) write(*,*) 'res',res
   call compute_lnprec(res,div2,div1,tree,lnprec)

   ! -- if presision check fail then recompute the amplitude in higher precision 
   if (lnprec(-2) > -3d0 .or. lnprec(-1) > -3d0 .or. abs_rel_err_dp > 2d-1&
        &.or. &
        &res(0)+1d0 == res(0) .or. res(1)+1d0 == res(1) .or. & ! true is result is Infinite 
        &res(0) /= res(0) .or. res(1) /= res(1) .or. &         ! true is result is NaN 
        &abs(res(0)+res(1)) > 1000d0*abs(tree)) then           

      if (recomputed_in_qp == 0) &
           &write(*,*) 'Switch to Multiple Precision if lnprec > -3, abs_rel_err > 2d-1, |V| > 1000|B|, found NaN, or infinite'  
      
       call redo_ampl_in_qp(pin,muu,i1,i2,i3,i4,css,tree_qp,res_qp)

      if ((recomputed_in_qp / 1000)*1000 == recomputed_in_qp) &
           &write(*,*) 'Amplitudes done in HP',recomputed_in_qp
      recomputed_in_qp = recomputed_in_qp+1 
      res= res_qp
   endif

      
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


   if (verbose)   write(*,*) 'Done getamplqqq_v', res, tree  



  end subroutine getamplqqqq_v


    !-------------------------

    !---------------------------
    ! For each channel (e.g. qqb) we now return EIGHT amplitudes: four primitives 
    ! for each of s- and t-channels.
    ! Primitives labeled by 1-4.
    ! Primitive 1 : case B1
    ! Primitive 2 : case B1 with gluons crossed in the middle
    ! Primitive 3 :  case B2 with gluon loop on 1st quark line
    ! Primitive 4 : case B2 with gluon loop on 2nd quark line  

    subroutine getamplqqqq_v_stu(pin,muu, chn, polesonly, tree_s,res_s, tree_t, res_t)

      complex(dp), intent(inout)        :: pin(:,:)
      real(dp), intent(in)           :: muu
      character, intent(in)          :: chn*3
      logical, intent(in)            :: polesonly
      complex(dp), intent(out)       :: tree_s(4), tree_t(4), res_s(5,-2:1), res_t(5,-2:1)
      !---------------------------------------------------------------------------
      real(dp)                       :: rat_s, rat_t
      logical                        :: csp
      complex(dp)                    :: Ls1, Ls2,Lt1,Lt2, treesplit_s(2), treesplit_t(2), treesplit(2)
      complex(dp)                    ::pW1(4),pW2(4)
      !---------------------------------------------------------------------------

      if (verbose) write(*,*) 'entering getamplqqqq_v_stu'

      tree_s = czero
      tree_t = czero
      res_s = czero
      res_t = czero
      pW1(:) = pin(3,:) + pin(4,:)
      pW2(:) = pin(5,:) + pin(6,:)
      ! check single pole
      csp = log_val_opt('-justepinv1', .false.)

    ! --- now calculate the amps
      ! Note: treesplit is dummy variable, not used. Variables actually used: treesplit_s & treesplit_t
    if (chn .eq. 'qqb') then
      
       ! s-channel
       if (verbose1) write(*,*) 's-channel B1'
       call getamplqqqq_v(pin,muu,1,2,7,8,1,tree_s(1),treesplit_s, res_s(1,-2:1),polesonly) ! Prim #1 - B1

       if (verbose1) write(*,*) 's-channel crossed gluons'
       call getamplqqqq_v(pin,muu,1,2,8,7,3,tree_s(2),treesplit, res_s(2,-2:1),polesonly) ! Prim #2 - B1 crossed gl
       
       if (verbose1) write(*,*) 's-channel B2-1'
       call getamplqqqq_v(pin,muu,1,2,7,8,2,tree_s(3),treesplit, res_s(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 's-channel B2-2' 
       call getamplqqqq_v(pin,muu,7,8,1,2,2,tree_s(4),treesplit, res_s(4,-2:1),polesonly) ! Prim #4 - B2           

       !t-channel
       if (verbose1) write(*,*) 't-channel B1'
       call getamplqqqq_v(pin,muu,1,8,7,2,1,tree_t(1),treesplit_t, res_t(1,-2:1),polesonly) ! Prim #1 - B1


       if (verbose1) write(*,*) 't-channel crossed gluons'                              
       call getamplqqqq_v(pin,muu,1,8,2,7,3,tree_t(2),treesplit, res_t(2,-2:1),polesonly) ! Prim #2 - B1 crossed gl
       
       if (verbose1)  write(*,*) 't-channel B2-1'

       call getamplqqqq_v(pin,muu,1,8,7,2,2,tree_t(3),treesplit,res_t(3,-2:1),polesonly) ! Prim #3 - B2
       
       if (verbose1) write(*,*) 't-channel B2-2'                                                        

       call getamplqqqq_v(pin,muu,7,2,1,8,2,tree_t(4),treesplit,res_t(4,-2:1),polesonly) ! Prim #4 - B2

       if (csp) call single_pole_check(pin,muu,1,2,7,8, tree_s(1), tree_t(1))


       if (ferm_loops_WW_nf) then
        
          ! Logs take 3 momenta as arguments
          call log_mom(pin(1,:),pin(2,:),pW1,muu,Ls1)
          call log_mom(pin(1,:),pin(2,:),pW2,muu,Ls2)

          call log_mom(pin(1,:),pin(8,:),pW1,muu,Lt1)
          call log_mom(pin(1,:),pin(8,:),pW2,muu,Lt2)
       endif

       
    elseif (chn .eq. 'qbq') then

       ! s-channel
       if (verbose1) write(*,*) 's-channel B1'
       call getamplqqqq_v(pin,muu,2,1,7,8,1,tree_s(1),treesplit_s,res_s(1,-2:1),polesonly)  ! Prim #1 - B1
       
       if (verbose1) write(*,*) 's-channel crossed gluons'
       call getamplqqqq_v(pin,muu,2,1,8,7,3,tree_s(2),treesplit,res_s(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 's-channel B2-1'
       call getamplqqqq_v(pin,muu,2,1,7,8,2,tree_s(3),treesplit,res_s(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 's-channel B2-2'
       call getamplqqqq_v(pin,muu,7,8,2,1,2,tree_s(4),treesplit,res_s(4,-2:1),polesonly)  ! Prim #4 - B2


      !t-channel     
       if (verbose1) write(*,*) 't-channel B1'
       call getamplqqqq_v(pin,muu,2,8,7,1,1,tree_t(1),treesplit_t,res_t(1,-2:1),polesonly)  ! Prim #1 - B1
       
       if (verbose1) write(*,*) 't-channel crossed gluons'
       call getamplqqqq_v(pin,muu,2,8,1,7,3,tree_t(2),treesplit,res_t(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 't-channel B2-1'
       call getamplqqqq_v(pin,muu,2,8,7,1,2,tree_t(3),treesplit,res_t(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 't-channel B2-2'
       call getamplqqqq_v(pin,muu,7,1,2,8,2,tree_t(4),treesplit,res_t(4,-2:1),polesonly)  ! Prim #4 - B2

!       if (csp) call single_pole_check(pin,muu,2,1,7,8, tree_s(1), tree_t(1))

       if (ferm_loops_WW_nf) then
          call log_mom(pin(2,:),pin(1,:),pW1,muu,Ls1)
          call log_mom(pin(2,:),pin(1,:),pW2,muu,Ls2)
          call log_mom(pin(2,:),pin(8,:),pW1,muu,Lt1)
          call log_mom(pin(2,:),pin(8,:),pW2,muu,Lt2)
       endif
     
    elseif (chn .eq. 'qqq') then

       ! s-channel
       if (verbose1) write(*,*) 's-channel B1'
       call getamplqqqq_v(pin,muu,1,7,2,8,1,tree_s(1),treesplit_s,res_s(1,-2:1),polesonly)  ! Prim #1 - B1
       
       if (verbose1) write(*,*) 's-channel crossed gluons'
       call getamplqqqq_v(pin,muu,1,7,8,2,3,tree_s(2),treesplit,res_s(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 's-channel B2-1'
       call getamplqqqq_v(pin,muu,1,7,2,8,2,tree_s(3),treesplit,res_s(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 's-channel B2-2'
       call getamplqqqq_v(pin,muu,2,8,1,7,2,tree_s(4),treesplit,res_s(4,-2:1),polesonly)  ! Prim #4 - B2
       


       !t-channel
       if (verbose1) write(*,*) 't-channel B1'
       call getamplqqqq_v(pin,muu,1,8,2,7,1,tree_t(1),treesplit_t,res_t(1,-2:1),polesonly)  ! Prim #1 - B1
       
       if (verbose1) write(*,*) 't-channel crossed gluons'
       call getamplqqqq_v(pin,muu,1,8,7,2,3,tree_t(2),treesplit,res_t(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 't-channel B2-1'
       call getamplqqqq_v(pin,muu,1,8,2,7,2,tree_t(3),treesplit,res_t(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 't-channel B2-2'
       call getamplqqqq_v(pin,muu,2,7,1,8,2,tree_t(4),treesplit,res_t(4,-2:1),polesonly)  ! Prim #4 - B2

       if (csp) call single_pole_check(pin,muu,1,7,2,8,tree_s(1), tree_t(1))

       if (ferm_loops_WW_nf) then
          call log_mom(pin(1,:),pin(7,:),pW1,muu,Ls1)
          call log_mom(pin(1,:),pin(7,:),pW2,muu,Ls2)
          call log_mom(pin(1,:),pin(8,:),pW1,muu,Lt1)
          call log_mom(pin(1,:),pin(8,:),pW2,muu,Lt2)
       endif
       
    elseif (chn .eq. 'qbb') then

       ! s-channel
       if (verbose1) write(*,*) 's-channel B1'
       call getamplqqqq_v(pin,muu,7,1,8,2,1,tree_s(1),treesplit_s,res_s(1,-2:1),polesonly)  ! Prim #1 - B1

       if (verbose1) write(*,*) 's-channel crossed gluons'
       call getamplqqqq_v(pin,muu,7,1,2,8,3,tree_s(2),treesplit,res_s(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 's-channel B2'
       call getamplqqqq_v(pin,muu,7,1,8,2,2,tree_s(3),treesplit,res_s(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) ',s-channel B2-2'
       call getamplqqqq_v(pin,muu,8,2,7,1,2,tree_s(4),treesplit,res_s(4,-2:1),polesonly)  ! Prim #4 - B2


       !t-channel
       if (verbose1) write(*,*) 't-channel B1'
       call getamplqqqq_v(pin,muu,7,2,8,1,1,tree_t(1),treesplit_t,res_t(1,-2:1),polesonly)  ! Prim #1 - B1

       if (verbose1) write(*,*) 't-channel crossed gluons'
       call getamplqqqq_v(pin,muu,7,2,1,8,3,tree_t(2),treesplit,res_t(2,-2:1),polesonly)  ! Prim #2 - B1 crossed g

       if (verbose1) write(*,*) 't-channel B2'
       call getamplqqqq_v(pin,muu,7,2,8,1,2,tree_t(3),treesplit,res_t(3,-2:1),polesonly)  ! Prim #3 - B2

       if (verbose1) write(*,*) 't-channel B2-2'
       call getamplqqqq_v(pin,muu,8,1,7,2,2,tree_t(4),treesplit,res_t(4,-2:1),polesonly)  ! Prim #3 - B2

       if (csp) call single_pole_check(pin,muu,7,1,8,2, tree_s(1), tree_t(1))

       if (ferm_loops_WW_nf) then
          call log_mom(pin(7,:),pin(1,:),pW1,muu,Ls1)
          call log_mom(pin(7,:),pin(1,:),pW2,muu,Ls2)

          call log_mom(pin(7,:),pin(2,:),pW1,muu,Lt1)
          call log_mom(pin(7,:),pin(2,:),pW2,muu,Lt2)

       endif

    endif

    ! Include the fermion loop.
    ! This has the same colour factor as for prim #1 - TaTbTbTa
    ! Prefactor - see hep-ph/9305239, and notes.
    if (ferm_loops_WW_nf) then        

       res_s(5,-1) = res_s(5,-1) + Nf*(-two/three)*tree_s(1)  
       res_t(5,-1) = res_t(5,-1) + Nf*(-two/three)*tree_t(1)  
                                                                 
       ! Use Q^2 of gluon - this will depend on which W is radiated by quark
       ! line, so split amplitude into two parts.

       res_s(5,0) = res_s(5,0) + Nf*(-two/three)*(Ls1*treesplit_s(1) + Ls2*treesplit_s(2)) 
       res_t(5,0) = res_t(5,0) + Nf*(-two/three)*(Lt1*treesplit_t(1) + Lt2*treesplit_t(2)) 

       res_s(5,1) = res_s(5,1) + Nf*(-10.0d0/9.0d0)*tree_s(1) 
       res_t(5,1) = res_t(5,1) + Nf*(-10.0d0/9.0d0)*tree_t(1)                                                                  

    endif
       
    ! Tree amplitudes should be the same regardless of primitive:
    
    rat_s = abs((tree_s(1) - tree_s(2))/tree_s(1))
    rat_t = abs((tree_t(1) - tree_t(2))/tree_t(1))

!    if (rat_s > 1e-12 .or. rat_t > 1e-12) then
!       write(*,*) 'Tree level amplidue in s-channel : ', tree_s
!       write(*,*) 'Tree level amplidue in t-channel : ', tree_t

!       stop 'Error: the tree-level amplitudes are different for different primitives'
!    endif
    
    
    tree_t = -tree_t  ! (-1) as swap fermions
    res_t = -res_t

        
  end subroutine getamplqqqq_v_stu

!----------------------------

  subroutine single_pole_check(p,mu,i1,i2,i3,i4,tree_s, tree_t)

    complex(dp), intent(in)     :: p(:,:)
    complex(dp), intent(in)     :: tree_s, tree_t
    real(dp), intent(in)        :: mu
    integer, intent(in)         :: i1,i2,i3,i4
!-----------
  
    complex(dp)                 :: ap1(4),ap2(4),ap3(4),ap4(4)
    complex(dp)                 :: Ls(4),Lt(4),s
    real(dp)                    :: c(8,2), mu2, hardscale
    complex(dp)                    :: cres_s, cres_t, cres_st
    complex(dp)                    :: born_s, born_t, born_st
    integer                     :: i

    mu2 = mu**2
    
    ap1 = p(i1,:)
    ap2 = p(i2,:)
    ap3 = p(i3,:)
    ap4 = p(i4,:)

    s = sc(p(1,:)+p(2,:),p(1,:)+p(2,:))
    hardscale = sqrt(abs(s))*2d0

    ap1 = ap1/hardscale
    ap2 = ap2/hardscale
    ap3 = ap3/hardscale
    ap4 = ap4/hardscale

    mu2 = mu2/hardscale**2

    Ls = czero
    Lt = czero
    
    ! s-channel B1
    s = two*sc(ap1,ap4)
    if (real(s) < zero) then 
       Ls(1) = Ls(1) -log(mu2/(-s))
    else
       Ls(1) = Ls(1) -log(mu2/s)-ci*pi 
    endif
    
    s = two*sc(ap2,ap3)
    if (real(s) < zero) then 
       Ls(1) = Ls(1) -log(mu2/(-s))
    else
       Ls(1) = Ls(1) -log(mu2/s)-ci*pi 
    endif


    Ls(1) = Ls(1) + two/three

    !s-channel B1 crossed gluons

    s = two*sc(ap1,ap3)
    if (real(s) < zero) then 
       Ls(2) = Ls(2) -log(mu2/(-s))
    else
       Ls(2) = Ls(2) -log(mu2/s)-ci*pi 
    endif
    
    s = two*sc(ap2,ap4)
    if (real(s) < zero) then 
       Ls(2) = Ls(2) -log(mu2/(-s))
    else
       Ls(2) = Ls(2) -log(mu2/s)-ci*pi 
    endif

    Ls(2) = Ls(2) + two/three

    ! s-channel B2-1
    s = two*sc(ap1,ap2)
    if (real(s) < zero) then 
       Ls(3) = Ls(3) -log(mu2/(-s))
    else
       Ls(3) = Ls(3) -log(mu2/s)-ci*pi 
    endif

    Ls(3) = Ls(3) - three/two

    ! s-channel B2-2
    s = two*sc(ap3,ap4)
    if (real(s) < zero) then 
       Ls(4) = Ls(4) - log(mu2/(-s))
    else
       Ls(4) = Ls(4) -log(mu2/s)-ci*pi 
    endif

    Ls(4) = Ls(4) - three/two


     ! t-channel B1
    s = two*sc(ap3,ap4)

    if (real(s) < zero) then 
       Lt(1) = Lt(1) -log(mu2/(-s))
    else
       Lt(1) = Lt(1) -log(mu2/s)-ci*pi 
    endif
   

    s = two*sc(ap1,ap2)

    if (real(s) < zero) then 
       Lt(1) = Lt(1) -log(mu2/(-s))
    else
       Lt(1) = Lt(1) -log(mu2/s)-ci*pi 
    endif

    Lt(1) = Lt(1) + two/three 


    !t-channel B1 crossed gluons

    s = two*sc(ap1,ap3)
    if (real(s) < zero) then 
       Lt(2) = Lt(2) -log(mu2/(-s))
    else
       Lt(2) = Lt(2) -log(mu2/s)-ci*pi 
    endif
    
    s = two*sc(ap2,ap4)
    if (real(s) < zero) then 
       Lt(2) = Lt(2) -log(mu2/(-s))
    else
       Lt(2) = Lt(2) -log(mu2/s)-ci*pi 
    endif

    Lt(2) = Lt(2) + two/three 

    ! t-channel B2-1
    s = two*sc(ap1,ap4)
    if (real(s) < zero) then 
       Lt(3) = Lt(3) -log(mu2/(-s))
    else
       Lt(3) = Lt(3) -log(mu2/s)-ci*pi 
    endif

    Lt(3) = Lt(3) - three/two

    ! t-channel B2-2
    s = two*sc(ap2,ap3)
    if (real(s) < zero) then 
       Lt(4) = Lt(4) -log(mu2/(-s))
    else
       Lt(4) = Lt(4) -log(mu2/s)-ci*pi 
    endif

    Lt(4) = Lt(4) - three/two

    call cmatrix4q(c)

    cres_s = czero    
    cres_t = czero
    cres_st = czero

    do i =1,4
       cres_t = cres_t + c(i+4,2)*Lt(i)*abs(tree_t)**2
    enddo

    do i = 1,4
       cres_s = cres_s + c(i,1)*Ls(i)*abs(tree_s)**2
    enddo

    do i = 1,4
       cres_st = cres_st - (c(i+4,1)*Lt(i))*tree_s*conjg(tree_t) - c(i,2)*Ls(i)*tree_t*conjg(tree_s)
    enddo

    cres_st = cres_st + cres_s + cres_t 
    ! Contribution from fermion loops
    if (ferm_loops_WW_nf) then
       cres_t = cres_t - 70.0d0/27.0d0*abs(tree_t)**2
       cres_s = cres_s - 70.0d0/27.0d0*abs(tree_s)**2
       cres_st = cres_st + 10.0d0/81.0d0*tree_s*conjg(tree_t) + 10.0d0/81.0d0*tree_t*conjg(tree_s) &
            & - 70.0d0/27.0d0*abs(tree_s)**2  - 70.0d0/27.0d0*abs(tree_t)**2 
    endif


    born_s = Cf*Nc/two*abs(tree_s)**2
    born_t = Cf*Nc/two*abs(tree_t)**2
    born_st = born_s + born_t  + Cf*tree_s*conjg(tree_t)
    
    write(*,*) 'Single pole virtual and Born mixing:'
    
    write(*,*) 's-channel only: ', two*real(cres_s/born_s)

    write(*,*) 't-channel only: ', two*real(cres_t/born_t)
    
    write(*,*) 's- and t-channels mixed : ', two*real(cres_st)/real(born_st)

    write(*,*) 's- and t-channels mixed plus s-channel : ', two*real(cres_st+cres_s)/real(born_st + born_s)


  end subroutine single_pole_check
    


!-------------------------------------
  ! Colour matrix

      
  subroutine cmatrix4q(c)
    real(dp), intent(out)    :: c(8,2)
    ! -------
    

    ! s-virt * s-tree
    c(1,1) = Cf**2*Nc/two - Cf/four                  
    c(2,1) = -one*(-Cf/two)                          
    c(3,1) = -Cf/four                                
    c(4,1) = -Cf/four                                
  

    ! t-virt * s-tree
    c(5,1) = Cf/four/Nc
    c(6,1) = -one*(Cf*Nc/four + Cf/four/Nc)
    c(7,1) = Cf/four/Nc
    c(8,1) = Cf/four/Nc
    
    ! s-virt * t-tree
    c(1,2) = c(5,1)
    c(2,2) = c(6,1)
    c(3,2) = c(7,1)
    c(4,2) = c(8,1)
    
    ! t-virt * t-tree
    c(5,2) = c(1,1)
    c(6,2) = c(2,1)
    c(7,2) = c(3,1)
    c(8,2) = c(4,1)

  end subroutine cmatrix4q
    

   !---------------------------------------


    subroutine mtotqqqq_v(pin, muu, chn, polesonly, mtot)
      real(dp), intent(in)          :: pin(:,:)
      real(dp), intent(in)          :: muu
      character, intent(in)         :: chn*3
      logical, intent(in)           :: polesonly
      real(dp), intent(out)         :: mtot(3,-2:1)
!----------------
      complex(dp)                   :: pc(size(pin,dim=1), size(pin,dim=2))
      complex(dp)                   :: tree_s(4), tree_t(4), res_s(5,-2:1), res_t(5,-2:1)
      complex(dp)                   :: tree_st(2), res_st(8,-2:1)
      real(dp)                      :: c(8,2)
      integer                       :: i,j
       

      mtot = zero

      if (verbose) write(*,*) 'entering mtotqqqq_v'

      ! -- change from real to complex momenta, and put energy first.

      do i=1,8
         pc(i,1) = cmplx(pin(i,4),0d0,dp)
         pc(i,2) = cmplx(pin(i,1),0d0,dp)
         pc(i,3) = cmplx(pin(i,2),0d0,dp)
         pc(i,4) = cmplx(pin(i,3),0d0,dp)
      enddo


      do i = 1,8
         if (verbose) write(*,*) 'pc',i, pc(i,:)
      enddo
  
      call getamplqqqq_v_stu(pc,muu, chn, polesonly, tree_s, res_s, tree_t, res_t)

      ! Note that all elements of tree_s are equal; same for tree_t.
      tree_st(1) = tree_s(1)
      tree_st(2) = tree_t(1)
      do i =1,4
         res_st(i,-2:1) = res_s(i,-2:1)
         res_st(i+4,-2:1) = res_t(i,-2:1)
      enddo

      call cmatrix4q(c)
      
      ! This is the product  of the s- and t-channels.
      do j = 1,2
         do i = 1,8
            mtot(1,-2:1) = mtot(1,-2:1) + two*real(c(i,j)*dconjg(tree_st(j))*res_st(i,-2:1))
         enddo
      enddo

    
      ! This is the tree level multiplied by the virtual for s-channel only.
      do i = 1,4
         do j = 1,1
            mtot(2,-2:1) = mtot(2,-2:1) + two*real(c(i,j)*dconjg(tree_st(j))*res_st(i,-2:1))
         enddo
      enddo
            
      
      ! This is the tree level multiplied by the virtual for t-channel only.
      ! Only happens for qqb or qbq case, never qq or qbqb cases.
      if (chn .eq. 'qqb' .or. chn .eq. 'qbq') then
        do i = 5,8
         do j = 2,2
            mtot(3,-2:1) = mtot(3,-2:1) + two*real(c(i,j)*dconjg(tree_st(j))*res_st(i,-2:1))
         enddo
      enddo
      endif



      if (ferm_loops_WW_nf) then
         
         mtot(1,-2:1) = mtot(1,-2:1) + two*real(res_s(5,-2:1)*dconjg(tree_s(1))+&
              &res_t(5,-2:1)*dconjg(tree_t(1))-one/three*res_s(5,-2:1)*dconjg(tree_t(1))&
              &-one/three*res_t(5,-2:1)*dconjg(tree_s(1)))  !Colour fac isFc*N/4=1 and -CF/4
         mtot(2,-2:1) = mtot(2,-2:1) + two*real(res_s(5,-2:1)*dconjg(tree_s(1))) !Colour fac is Cf*N/4=1
         mtot(3,-2:1) = mtot(3,-2:1) + two*real(res_t(5,-2:1)*dconjg(tree_t(1))) !colour fac is CF*N/4=1         
      endif

    end subroutine mtotqqqq_v
    
!-----------------------------------------


    subroutine log_mom(p1,p2,p3,muu,L)
      complex(dp), intent(in)    :: p1(:), p2(:),p3(:)
      real(dp), intent(in)           :: muu
      complex(dp), intent(out)   :: L
!------------
      complex(dp)                :: s, psum(4)
      

      psum(:) = p1(:)+p2(:) + p3(:)
      s = sc(psum,psum)
!      s = two*sc(p1,p2)

      if (real(s) < zero) then 
          L = log(muu**2/(-s))
       else
          L = log(muu**2/s)+ci*pi 
       endif

     end subroutine log_mom
      
  subroutine compute_lnprec(res,rdiv2,rdiv1,tree,lnprec) 
    complex(dp), intent(in) :: res(-2:1),rdiv2,rdiv1,tree
    real(dp), intent(out)   :: lnprec(-2:-1) 
    real(dp)                :: ln10 = 2.3025850929940456840_dp 
    real(dp), parameter     :: min_prec = 1e-16_dp
    

    if (abs(rdiv2) > zero ) then 
       lnprec(-2) = log(min_prec+abs((res(-2)/tree-rdiv2)/rdiv2))/ln10
    else
       if (abs(res(-2)) > zero) then 
          lnprec(-2) = log(abs(res(-2))) 
       else
         lnprec(-2) = -20._dp
      endif
    end if

    if (abs(rdiv1) > zero ) then 
       lnprec(-1) = log(min_prec+abs((res(-1)/tree-rdiv1)/rdiv1))/ln10 
    else
       if (abs(res(-1)) > zero) then 
          lnprec(-1) = log(abs(res(-1)))
       else
         lnprec(-1) = -20._dp 
      endif
    end if
  end subroutine compute_lnprec



  end module qqqqampl_v
      
    
