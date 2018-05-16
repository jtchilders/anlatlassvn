module qpevent 
  use define_ampl, only : swap 
  use types;  use consts_qp 
  use warnings_and_errors 
  implicit none
  private 

  public :: redo_event_in_qp  
  logical :: verbose = .false. 

contains

  ! given a dpevent generate an qpevent 
  ! make sure that the Npart mass conditions and the 4 energy-momentum
  ! conservation conditions are fulfilled in multiple precision (not
  ! only up to 15 digits..)
  subroutine generate_massless_qpevent(pext_dp,pext_qp,rn) 
    complex(dp), intent(in)       :: pext_dp(:,:)
    complex(qp), intent(out) :: pext_qp(:,:) 
    logical, intent(in)           :: rn 
    ! --------------------------------------------
    complex(qp) :: pout(4),ptmp(4) 
    integer  :: Npart, i, ninc, iinc(2),abeam(2),j,iout,ipt 
    complex(qp) :: pin(4),po(4),pinp(4),pop(4)
    complex(qp) :: Ei,Eip,pzi,pzip,Eo,Eop,pzo,pzop,ptop2
    ! needed for fixed point 
    complex(qp) :: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
    real(qp) :: mu,en,theta,phi,al,ga,cbeta,sbeta,eg3,eg4,eg5,eg6
    

    Npart = size(pext_dp,dim=1)
    if (Npart < 4) call wae_error('generate_qpevent','Npart<4')

    if (rn) then 

       ! -- get the incoming partons and the direction of the beam (x,y,z ?) 
       ninc = 0 
       do i=1,Npart
          do j=2,4
             if (abs(abs(pext_dp(i,j))-abs(pext_dp(i,1))) < sq2tol) then 
                ninc = ninc+1 
                if (ninc >2) stop 'generate_massless_qpevent: more than 2 incoming?'
                iinc(ninc) = i 
                abeam(ninc) = j 
             end if
          end do
       end do
       if (abeam(1) /= abeam(2)) stop 'generate_massless_qpevent: different abeams?' 

       pext_qp = czero 
       do i=1,Npart
          if (i /= iinc(1) .and. i /= iinc(2)) then 
             pext_qp(i,2) = pext_dp(i,2) 
             pext_qp(i,3) = pext_dp(i,3) 
             pext_qp(i,4) = pext_dp(i,4) 
             pext_qp(i,1) = sqrt(pext_qp(i,2)**2+pext_qp(i,3)**2+pext_qp(i,4)**2)
          endif
       enddo
       do i=1,4 
          pout(i) = sum(pext_qp(1:Npart,i))
       end do

       ! after that have 2 * 4 components yet to be fixed and 6 more
       ! conditions to be fullfilled

       ! set iinc(1) along z axis (fix 2 components) 
       pext_qp(iinc(1),2:4) = czero 
       do j=2,4
          if (j /= abeam(1)) then 
             pext_qp(iinc(2),j) = -pout(j) ! iinc(2) has tiny pt components now... 
          end if
       end do

       pext_qp(iinc(1),1) = (-pout(1)**2+pout(4)**2)
       do j=2,4 
          if (j /= abeam(1)) then 
             pext_qp(iinc(1),1) = pext_qp(iinc(1),1)+pext_qp(iinc(2),j)**2
          endif
       enddo
       pext_qp(iinc(1),1) = pext_qp(iinc(1),1)/(two*(pout(1)-pout(abeam(1))))

       pext_qp(iinc(1),abeam(1)) = pext_qp(iinc(1),1)


       pext_qp(iinc(2),1) = (-pout(1)**2+two*pout(1)*pout(abeam(1))-pout(abeam(1))**2)
       do j=2,4 
          if (j /= abeam(1)) then 
             pext_qp(iinc(2),1) = pext_qp(iinc(2),1)-pext_qp(iinc(2),j)**2
          endif
       enddo
       pext_qp(iinc(2),1) = pext_qp(iinc(2),1)/(two*(pout(1)-pout(abeam(1))))

       pext_qp(iinc(2),abeam(2)) = (pout(1)**2-two*pout(1)*pout(abeam(1))+pout(abeam(1))**2)

       do j=2,4 
          if (j /= abeam(1)) then 
             pext_qp(iinc(2),abeam(2)) = pext_qp(iinc(2),abeam(2))-pext_qp(iinc(2),j)**2
          endif
       enddo
       pext_qp(iinc(2),abeam(2)) = pext_qp(iinc(2),abeam(2))/(two*(pout(1)-pout(abeam(1))))

       ipt = iinc(2) 
       if (real((pext_qp(iinc(1),abeam(1)))*pext_dp(iinc(1),abeam(1))) < zero) then 
          ! - need to swap incoming ones 
          ptmp =  pext_qp(iinc(1),:)
          pext_qp(iinc(1),:) = pext_qp(iinc(2),:)
          pext_qp(iinc(2),:) = ptmp 
          ipt = iinc(1) 
       endif
       
       ! -- adjust kinematics -- let one outgoing parton absort the pt components
       pin = pext_qp(ipt,:) 
       do i=1,Npart 
          if (i /= iinc(1) .and. i /= iinc(2)) then 
             po = pext_qp(i,:) 
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
       if (real(Eip*Ei,qp) > zero .and. real(pzip*pzi,qp)> zero) then 
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
          
          if (real(Eip*Ei,qp) < zero .or. real(pzip*pzi,qp) < zero) write(*,*) 'could not find right solution?'
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
       
       pext_qp(ipt,:) = pinp 
       pext_qp(iout,:) = pop 

    else
       include 'BDKpoint7qp.'
    endif


  end subroutine generate_massless_qpevent

  subroutine generate_W_qpevent(pext_mless,pext,pl,pa) 
    complex(qp), intent(in)  :: pext_mless(:,:)
    complex(qp), intent(out) :: pext(:,:),pl(:),pa(:)
    integer :: i, Npart 

    Npart = size(pext_mless,dim=1)
    if (size(pext,dim=1)+1 /= size(pext_mless,dim=1)) &
         &stop 'generate_W_qpevent: size(pext,dim=1)+1 /= size(pext,dim=1)'

    if (size(pext,dim=2) /= size(pext_mless,dim=2)) &
         &stop 'generate_W_qpevent: size(pext,dim=2) /= size(pext,dim=2)'


    ! change notation: first index is now the particle label 
    ! and end index the component, energy compoenent is now the 1st one            
    do i=1,2
       pext(i,:) = pext_mless(i,:)
    enddo

    ! construct W momentum as sum of two leptons         
    pa(:) = pext_mless(Npart-1,:)
    pl(:) = pext_mless(Npart,:)

    pext(3,:) = pl(:)+pa(:)
    mz = sqrt(pext(3,1)**2-pext(3,2)**2-pext(3,3)**2-pext(3,4)**2)

    do i=4,Npart-1
       pext(i,:) = pext_mless(i-1,:)
    enddo

  end subroutine generate_W_qpevent

  subroutine redo_event_in_qp(pmless_dp,mu_dp,hl,outgoing,rn,res,massless_in)
    use define_ampl; use sub_defs_io  
    use qpglobal; use qptodp_global
    use qpmemory 
    use comb; use qpinitWpWp
    use qpopp; use masters; use qpspinors 
!    use qpopp; use qpmasters; use qpspinors 
    use qpauxiliary_functions ! tmp 
    complex(dp), intent(in)   :: pmless_dp(:,:)
    real(dp), intent(in)      :: mu_dp
    integer, intent(in)       :: hl(:)
    logical, intent(in)       :: outgoing(:),rn
    complex(dp), intent(out)  :: res(-2:1)
    logical, intent(in), optional      :: massless_in
    ! -----------------------------------------------       
    complex(qp) :: pmless_qp(size(pmless_dp,dim=1),&
         &size(pmless_dp,dim=2))
    complex(qp) :: pext(size(pmless_dp,dim=1),&
         &size(pmless_dp,dim=2))
!    complex(qp) :: pext(size(pmless_dp,dim=1)-1,&
!         &size(pmless_dp,dim=2))
    complex(qp)   :: pl(4),pa(4)
    logical, save :: first_time = .true. 
    !    integer :: init_prec 
    integer,save :: lncuts(5)
    integer :: N1,N2,N3,N4,N5,tN2,tN3,tN4,tN5,xe
    complex(qp) :: tree,rdiv1, rdiv2
    real(qp) :: mu
    integer :: i,j,np 
    logical :: massless 

    if (present(massless_in)) then 
       massless = massless_in
    else
       massless = .false. 
    endif
    if (first_time) then 
       ! -- MP initialization  
       !init_prec =  int_val_opt('-init_prec',32)
       !call mpinit 
       !call SetPrecision(init_prec)
       !call fixmpconsts(init_prec)
       call init_spinors 
       call fix_bvec 
    endif

    ! generate qpevent from dp one 
    !do i=1,7
    !   write(*,*) 'DP:pi',i,pmless_dp(i,:)
    !enddo
    call generate_massless_qpevent(pmless_dp,pmless_qp,rn)
    if (verbose) write(*,*) 'generated massless qp event' 
    do i=1,size(pmless_qp,dim=1)
       if (any(abs(pmless_qp(i,:)-pmless_dp(i,:)) > sqrt(sq2tol))) then 
          write(*,*) 'momenta do not match'
          write(*,*) 'DP:pi',i,(pmless_dp(i,:))
          write(*,*) 'QP:pi',i,(pmless_qp(i,:)),dot(pmless_qp(i,:),pmless_qp(i,:))
       endif
!       if ((abs(dot(pmless_qp(i,:),pmless_qp(i,:))) > 1000*tol)) then 
!          write(*,*) 'not onshell?',i,(pmless_qp(i,:)),dot(pmless_qp(i,:),pmless_qp(i,:)) 
!       endif
    enddo

    if (any(abs(sum(pmless_qp(:,:),dim=1)) > 1000*tol)) then 
       write(*,*) 'no en conservation? ',sum(pmless_qp(:,:),dim=1)
    endif

    if (massless) then 
       pext = pmless_qp 
       np = size(pext,dim=1)
    else
       call generate_W_qpevent(pmless_qp,pext(:size(pext,dim=1)-1,:),pl,pa)         
       np = size(pext,dim=1)-1
    endif
    mu = mu_dp
    if (verbose) write(*,*) 'generated W qp event' 


    if (first_time) then 
       ! -- memory 
       call allocate_mem

       ! -- momentum flow info 
       ! args are ncut,npoint ncut = ncut_max = npoint -2 
       !call allocate_mominfo(ncut,npoint) 
       call allocate_mominfo(npoint-2,npoint) 
    endif
    log_ExtCurrent = .false. 
    ! get the number of cuts depends on npoint only 
    if (ferm_loops .or. ferm_loops_Z) then 
       call getcutnumb_nf(lncuts,nterm,npoint)  
    else
       call getcutnumb(lncuts,npoint)  
    endif
    
    call pampl_count_gen(ampl_type,lncuts,Npoint,&
         &tN5,N5,tN4,N4,tN3,N3,tN2,N2,N1)  

    if ((N5>N5max .or.  N4>N4max .or.  &
         &N3>N3max .or.  N2>N2max .or. tN2>tN2max)) then 
       write(*,*) 'ampl_type',ampl_type
       write(*,*) 'lncuts',lncuts
       write(*,*) 'Cur: n2,n3,n4,n5',n2,n3,n4,n5
       write(*,*) 'Max: n2,n3,n4,n5',n2max,n3max,n4max,n5max
       write(*,*) 'npoint,nterm,ncut',npoint,nterm,ncut 
       write(*,*) 'redo_event_in_mp: something went wrong in setting of N5,N4,N3,N2'
    endif

    if (first_time) then 
       ! Lci will contain info about cut labeling (integer) 
       ! Fi  will contain info about flavour (char) 
       ! Yci will contain info about "reference for external 
       ! momenta ordering (integer) 
       call allocate_arrays(N5max,N4max,N3max,N2max)
       call allocate_coeffs(N5max,N4max,N3max,N2max,N1max)
       first_time = .false. 
    endif

    ! returns Lci, Fi, Yci as well as tree level and poles of oneloop 
    call pampl_gen(ampl_type,pext(:np,:),pl,pa,hl,lncuts,tN5,N5,tN4,N4,tN3,N3,tN2,N2,N1,&
         &Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:),Yc4(:N4,:),&
         &Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:),&
         &tree,mu,rdiv2,rdiv1,outgoing)

    if (writeout)  write(*,*) 'QP: TREE AMPLITUDE', (tree)

    call pentcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:))

    call quadcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:),Yc4(:N4,:),&
         &Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:))

    call tripcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:),Yc4(:N4,:),&
         &Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:))

    call doubcut(Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),Lc4(:N4,:),F4(:N4,:),Yc4(:N4,:),&
         &Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:))


    ! call singlecut ! there is no single cut for W + jets

    ! now switch to DP ned to allocate DO global and make transfer QP info in DP 
    ! then mastersubs uses DP info 
    call qptodp_glob(tagdcut,mom,hel,momline,Lab_ex, &
         &Lab_in,coeff5,propv5,coeff4,refvect4,propv4,&
         &coeff3,refvect3,propv3,coeff2,refvect2,propv2,coeff1,&
         &dcoeff5,dcoeff4,dcoeff3,dcoeff2,dcoeff1,&
         &mass5,mass4,mass3,mass2,mass1,.true.)!,en,&
    !             &ptcut,etacut, eg,noalloc=.true.)


!    call MasterSubs(mu,Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),&
    call MasterSubs(mu_dp,Lc5(:N5,:),F5(:N5,:),Yc5(:N5,:),&
         &Lc4(:N4,:),F4(:N4,:),Yc4(:N4,:),Lc3(:N3,:),F3(:N3,:),Yc3(:N3,:),&
         &Lc2(:N2,:),F2(:N2,:),Yc2(:N2,:),res)

    if ((ampl_type == iextra_ferm_pair_nf) .or. (ampl_type == iferm_loops) &
         &.or. (ampl_type == iferm_loops_Z) .or. (ampl_type == iferm_loops_Z_sbs))&
         & res = -res
    !if (ampl_type == iextra_ferm_pair_nf .or. ampl_type == iferm_loops) rdiv2 = -rdiv2
    !if (ampl_type == iextra_ferm_pair_nf .or. ampl_type == iferm_loops) rdiv1 = -rdiv1
    
    ! -- include ci form loop integral here 
    res = (0.0_dp,1.0_dp)*res 

    if (writeout) then 
       print *, 'QP: FINAL RESULT'
       if (abs(tree).gt.propcut) then 
          do xe=-2,1 
             print *, (res(xe)/(tree)) ,(abs(res(xe)/(tree)))
          enddo
       else 
          print *, 'qp: tree amplitude vanishes'
          do xe=-2,1
             print *, (res(xe))
          enddo
       endif
       write(*,*) 'qp: cut+rat', ((res(0)+res(1))/(tree))

       print *, 'QP: DIVERGENT PARTS'
       print *, 'QP: xe=-2         ', rdiv2
       print *, 'QP: xe = -1       ', rdiv1
    endif

  end subroutine redo_event_in_qp


end module qpevent
