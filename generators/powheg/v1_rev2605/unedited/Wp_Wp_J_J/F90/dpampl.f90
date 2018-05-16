!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECampl.f90.

!! File generated automatically by autogen.pl from 
!! general precision template dpfiles/gendp/dpampl.f90.

module dpamplitude
  use types; use consts_dp ; use define_ampl
  use dppol_int; use dpaux_functions; use dpvvn
  use dprecurrence;  use dpglobal 
  use dpmemory; use dpcut_utils
  implicit none
  private 

  public :: givepol, ampl, calc_ampl

  logical :: verbose = .false. 
  logical :: use_mynivecs = .true. 

contains


  !-------- polarization vectors for fermions and bosons
  subroutine givepol(lab,q1,Dv,Ds,tag_pol,Nj,BPOL,POL)
    character, intent(in)      :: lab*3 
    integer, intent(in)        :: tag_pol,Ds,Dv
    complex(dp), intent(in)  :: q1(5) 
    integer, intent(out)       :: Nj
    complex(dp), intent(out) :: BPOL(8,16),POL(8,16)
    ! ---------------------------------------------------
    complex(dp) :: kuks, r12, r14
    complex(dp) :: p(Dv), q14(4),tmp
    complex(dp) :: v(4,4), v1(4),v2(4)
    complex(dp) :: ulist(Ds,8),barulist(Ds,8)
    real(dp)    :: mass 
    integer       :: Nhf,Nhg,i,j
    ! ------------------------------------
    complex(dp) :: ni(4,3), pred(4,1) 

    !----       initializing polarization lists
    POL=czero
    BPOL=czero

    if (verbose) write(*,*) 'entering givepol',(q1)

    ! tag_pol always zero for use (was used for top) 
    if (tag_pol == 0) then 

       if (Dv == 4) then 
          Nhf=2
          Nhg=2
       elseif (Dv == 6) then 
          Nhf=4
          Nhg=4
       elseif (Dv == 8) then 
          Nhf=8
          Nhg=2
       else
          stop 'givepol: Dv not allowed'
       endif

    elseif (tag_pol == 1) then 

       if (Dv == 4) then 
          Nhf=2
          Nhg=4
       elseif (Dv == 6) then 
          Nhf=4
          Nhg=6
       elseif (Dv == 8) then 
          Nhf=8
          Nhg=2
       else
          stop 'givepol: Dv not allowed'
       endif

    else
       stop 'givepol: tag_pol not recognized' 
    endif


    if (quark_flavour(lab)) then 

       if(Dv == 4) then    
          do i=1,4
             p(i)=-q1(i)
          enddo
       else
          do i=1,5
             p(i)=-q1(i)
          enddo
          do i=6,Dv
             p(i)=czero
          enddo
       endif

       if (lab == 'top') then 
          mass = mt
       elseif (lab == 'bot') then 
          mass = mb
       else
          mass = zero 
       endif
       call give_usp(Nhf,Dv,Ds,p,mass,ulist)
       call give_barusp(Nhf,Dv,Ds,p,mass,barulist)

       do i=1,Nhf
          do j=1,Ds
             BPOL(i,j)=barulist(j,i)
             POL(i,j)=ulist(j,i) 
          enddo
       enddo

       Nj=Nhf


    elseif (lab == 'glu') then 

       if (Dv == 4) then 

          call give1to4vect_light(q1(1:4),v)
          if (verbose) write(*,*) 'in give_pol:Ds=4', (v) 
          do i=1,4
             POL(1,i)=v(3,i)
             POL(2,i)=v(4,i)
             BPOL(1,i)=v(3,i) 
             BPOL(2,i)=v(4,i)
          enddo


          do i=1,4
             v1(i)=v(1,i)
             v2(i)=v(2,i)
          enddo

          r12 =  sc(v1,v2)
          kuks=1/sqrt(two*r12)

          do i=1,4
             POL(3,i)=kuks*(v1(i)+v2(i))
             BPOL(3,i)=POL(3,i)

             POL(4,i)=ci*kuks*(v1(i)-v2(i))
             BPOL(4,i)=POL(4,i)
          enddo

          if (verbose) write(*,*) 'in give_pol:Ds=4:pol1', (POL(1,:)) 

       elseif (Dv == 6) then  

          do i=1,4
             q14(i)=q1(i)
          enddo

          r14 = sc(q14,q14)
          r14=sqrt(r14)

          if (use_mynivecs) then 
             pred(:,1) = q14
             call compute_ni(pred,bvec,ni) 
             tmp = sc(q14,q14)
             v(1,:)= q14/tmp
             v(2,:)=ni(:,1)
             v(3,:)=ni(:,2)
             v(4,:)=ni(:,3)
          else
             call give1to4vect(q14,v)
          endif

          do i=1,4          
             POL(1,i)=v(2,i)
             BPOL(1,i)=v(2,i)
             POL(2,i)=v(3,i)
             BPOL(2,i)=v(3,i)
             POL(3,i)=v(4,i)
             BPOL(3,i)=v(4,i)
             POL(4,i)=czero
             BPOL(4,i)=czero
             POL(5,i)=czero
             BPOL(5,i)=czero
             POL(6,i)=q14(i)/r14
             BPOL(6,i)=POL(6,i)
          enddo

          POL(4,5)=czero
          BPOL(4,5)=czero
          POL(4,6) =ci
          BPOL(4,6)=ci
          POL(5,5) =ci
          BPOL(5,5)=ci

          POL(6,5)=czero
          BPOL(6,5)=czero
          POL(6,6)=czero
          BPOL(6,6)=czero            
          if (verbose) write(*,*) 'in give_pol:Ds=6:pol1', (POL(5,:)) 

       elseif (Dv == 8) then 

          do i=1,6
             POL(1,i)=czero
             BPOL(1,i)=czero
             POL(2,i)=czero
             BPOL(2,i)=czero
          enddo

          POL(1,7)=ci
          BPOL(1,7)=ci
          POL(2,7)=czero
          BPOL(2,7)=czero
          POL(2,8)=ci
          BPOL(2,8)=ci

       else
          stop 'givepol: gluon case: Vd not allowed' 
       endif ! end if Dv 

       Nj=Nhg

    else
       stop 'givepol: lab undefined' 
    endif ! end if lab 
    if (verbose) write(*,*) 'in givepol', &
         &lab, Dv, Ds, tag_pol, Nj, (POL(1,:))

  end subroutine givepol

  !--------- procedure to compute scattering amplitudes
  ! Dv/Ds dimensions 
  ! Nj1, Nj2 number of polarizations
  ! q1/q2 two cut momenta, POLI, POLF their polarization              
  ! lab1/lab2 their flavour       
  ! ia and ib are positions of the cut lines 
  subroutine ampl(Dv,Ds,Nj1,Nj2,POLI,q1,lab1,tg,ia,ib,BPOLF,q2,lab2,ij,mur)
    integer, intent(in)        ::  Dv, Ds,Nj1,Nj2,tg,ia,ib
    character, intent(in)      ::  lab1*3,lab2*3
    complex(dp), intent(in)  ::  POLI(8,16),BPOLF(8,16)
    complex(dp), intent(in)  ::  q1(5), q2(5)
    integer, intent(in)        ::  ij
    complex(dp), intent(out) ::  mur(10,10)
    ! -----------------------------------------------------------       
    complex(dp) :: res
    complex(dp) :: po(ib-ia+3,Ds),pe(ib-ia+3,Dv)
    complex(dp) :: p(Dv,ib-ia+3),sp(Ds,ib-ia+3)
    complex(dp) :: vg(Dv), vs(Ds),pol_glu(Dv)
    integer       :: ng,nq,nw
    integer       :: i,np,i1,i2,j,j2,j1
    character     :: Lab_ex1(ib-ia+3)*3 
    character     :: fl(ib-ia+3)*3
    integer       :: iarray(ib-ia+3) 

    i1=tg
    np = ib-ia+3

    iarray(1) = 0 ! dummy 
    do i=1,ib-ia+1
       i2=i1 + (i-1)

       if (i2 > Npoint) i2=i2 - Npoint

       do j=1,4
          po(i+1,j)=hel(ij,i2,j)
          pe(i+1,j)=mom(ij,i2,j)
       enddo
       do j=5,Ds
          po(i+1,j)=czero
       enddo
       do j=5,Dv
          pe(i+1,j)=czero
       enddo
       Lab_ex1(i+1)=Lab_ex(ij,i2)
       iarray(i+1) = identity(ij,i2)
    enddo
    iarray(np) = 0 ! internal line 

    ! -- now definining particle content of the amplitude           
    !    flavors/spins/momenta

    ng=0
    nq=0
    nw=0


    fl(1) =lab1
    fl(np)=lab2

    do i=2,np-1
       fl(i)=Lab_ex1(i)
    enddo

    do i=1,np
       if (fl(i) == 'glu') ng = ng+1
       if (quark_flavour(fl(i))) nq = nq+1
       if ((fl(i) == 'www') .or. (fl(i) == 'wwm') .or. (fl(i) == 'wwp')) nw = nw+1
    enddo

    if(Dv.eq.4) then 
       do j=1,4 
          p(j,1)=q1(j)
          p(j,np)=-q2(j)          
       enddo
    else
       do j=1,5
          p(j,1)=q1(j)
          p(j,np)=-q2(j)
       enddo
       do j=6,Dv
          p(j,1)=czero
          p(j,np)=czero
       enddo
    endif

    do j=1,Dv
       do i=2,np-1
          p(j,i)=pe(i,j)
       enddo
    enddo


    do j=1,Ds
       do i=2,np-1
          sp(j,i)=po(i,j) 
       enddo
    enddo


    ! -- memory: ri-initialize here everything        
    call riinitialize_mem

    do j2=1,Nj2

       do j=1,Ds
          sp(j,np)=BPOLF(j2,j)
       enddo

       sp(:,1) = czero

       call calc_ampl(Dv,Ds,ng,nq,nw,sp,p,fl,vg,vs,iarray,j2)

       do j1=1,Nj1

          do j=1,Ds
             sp(j,1)=POLI(j1,j)
          enddo

          if (fl(1).eq.'glu') then 
             pol_glu= sp(1:Dv,1)
             res = sc(pol_glu,vg)

            elseif(quark_flavour(fl(1))) then  
             res = psp1(vs,sp(:,1))                  
          else
             stop 'ampl: fl(1) undefined'
          endif
          if (verbose) write(*,*) 'in ampl:res', (res) 

          mur(j1,j2)=res

       enddo

    enddo
  end subroutine ampl

  ! Dv and Ds are the space and spin dimensions  
  ! ng, nq, nw denote the number of gluons, quarks and w bosons              
  ! sp and p are the spin/polarization and momentum information 
  ! fl is the flavour of everything  
  ! vg and vs are the results
  subroutine calc_ampl(Dv,Ds,ng,nq,nw,sp,p,flampl,vg,vs,iarray_in,pol_int_in)
    complex(dp), intent(in)    :: sp(:,:), p(:,:)
    character, intent(in)        :: flampl(:)*3
    integer, intent(in)          :: Dv,Ds,ng,nq,nw
    integer, intent(in),optional :: iarray_in(:)
    integer, intent(in),optional :: pol_int_in  
    complex(dp), intent(out)   :: vg(Dv),vs(Ds)
    ! -----------------------------------------------
    complex(dp) :: pol_glu(Dv,ng), p_glu(Dv,ng)
    complex(dp) :: pol_q(Ds,nq), p_q(Dv,nq)
    complex(dp) :: p_W(Dv), pol_W(Dv)
    complex(dp) :: p_Wp(Dv), pol_Wp(Dv)
    complex(dp) :: p_Wm(Dv), pol_Wm(Dv)
    complex(dp) :: p_WW(Dv,nw), pol_WW(Dv,nw), holder(Dv,nw)
    character     :: fl1(nq)*3      ! fermion flavors
    integer       :: pos1(1),pos2(2),pos3(3),pos4(4),pos6(6), pos5(5),posW
    integer       :: ng1,ng2, ng3,ng4,ng5, iW, sw, detstr, iWW(2)  
    integer       :: np             ! number of particles
    integer       :: i, ig,iq,iWp,iWm
    integer       :: giarray(ng),qiarray(nq),wiarray(2), oneWiarray
    integer :: calc
    integer :: iarray(size(p,dim=2)), pol_int
    character :: fl(size(flampl))*3
   logical::isstr
    

    fl = flampl 
    do i=1,size(flampl)
       if ((flampl(i) == 'ch1') .or. (flampl(i) =='ch2')) fl(i) = 'chr' 
       if ((flampl(i) == 'to1') .or. (flampl(i) =='to2')) fl(i) = 'top' 
    enddo


    vg = czero 
    vs = czero 
    if (present(pol_int_in)) then 
       if (size(iarray) /= size(iarray_in)) stop 'calc_ampl: array wrong size'
       pol_int = pol_int_in 
       iarray = iarray_in
    else
       pol_int = -1
       iarray = -1
    endif

    if (verbose) write(*,*) 'in calc_ampl: iarray,fl,nq,ng,nw', iarray,fl,nq,ng,nw

    calc = 0

    np=size(p,dim=2)


    if (nw == 0) then ! first case, no vector boson

       if (fl(1) == 'glu'.and.nq == 0) then    

          pol_glu(:,:) = sp(1:Dv,:)
          p_glu = p
          giarray = iarray 

          vg=vgluon(pol_glu(:,2:ng),p_glu(:,2:ng),giarray(2:ng),pol_int)
          calc = 1 

       elseif (fl(1) == 'glu'.and.nq == 2) then

          iq = 0
          ig = 0
          
          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos2,ig,iq,giarray,qiarray)

          ng1 = pos2(1) -2
          ng2 = pos2(2) - pos2(1) - 1


          pol_glu(:,ng) = sp(1:Dv,1)

          if (fl1(1) == fl1(2)) then 
             if (size(pol_glu,dim=2) + 2 /= np) stop 'wrong size' 
             vg=g_fbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                  &pol_q(:,1),p_q(:,1),fl1(1),&
                  &pol_q(:,2),p_q(:,2),fl1(2),ng1,ng2,&
                 &giarray(1:ng-1),qiarray(:nq),pol_int)
          else 
             vg = czero
          endif
          calc = 1 


       elseif (((fl(1) == 'top').or.(fl(1) == 'bot') .or.(fl(1) == 'str').or. &
            &(fl(1) == 'chr') .or. (fl(1) == 'dwn') ) .and. nq == 2) then 


          iq = 0 
          ig = 0
          
          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,fl1,pos1,ig,iq,giarray,qiarray)


          ng1 = pos1(1) -2

          if (fl(1) == fl1(1)) then 

             if (size(pol_glu,dim=2)+2 /= np) &
                  &stop 'calc_amp: size(pol_glu,dim=2)+2 /= np)'
             vs=f(pol_glu,p_glu,pol_q(:,1),p_q(:,1),fl1(1),fl(1),ng1,&
             &giarray(1:ng),qiarray(1:1),pol_int)
          else 
             vs = czero 
          endif
          calc = 1 


       elseif ((((ferm_loops_Z.eqv..true.)) &
            &.and. (fl(1) == 'chr').and.nq == 4)) then 
          
             iq = 1
             ig = 0

          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)
          pos3(1:3) = pos4(2:4) 
          
          ng1 = pos3(1) -2
          ng2 = pos3(2) - pos3(1) -1
          ng3 = pos3(3) - pos3(2) - 1
          vs=f_fbfbf_3(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
          
             calc = 1

       elseif ( ((extra_ferm_pair2.eqv..false.).and.(WWqqqq .eqv. .false.) .and.&
            &(extra_ferm_pair_nf.eqv..false.) .and. &
            &(ferm_loops_Z .eqv. .false.) ) .and.  &
            &(((fl(1) == 'top').or.(fl(1) == 'bot').or.&
               &(fl(1) == 'str').or.(fl(1) == 'chr')).and.nq == 4)) then 

             iq = 1
             ig = 0
             
          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)
          pos3(1:3) = pos4(2:4) 
          
              ng1 = pos3(1) -2
              ng2 = pos3(2) - pos3(1) -1
              ng3 = pos3(3) - pos3(2) - 1
              
              if (qbq_WW_and_gluons) then 
                 vs=f_bffbf_2(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                      &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
              else
                 vs=f_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                      &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
              endif
             calc = 1

       elseif ((extra_ferm_pair2.eqv..true.) .and. (WWqqqq .eqv. .false.) .and.&
            &((fl(1) == 'top'.or.fl(1) == 'bot'.or.&
               &fl(1) == 'str').and.nq == 4)) then 

             iq = 1
             ig = 0


          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)
          pos3(1:3) = pos4(2:4) 
          
          ng1 = pos3(1) -2
          ng2 = pos3(2) - pos3(1) -1
          ng3 = pos3(3) - pos3(2) - 1

             vs=f_bffbf_2(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray,qiarray(2:nq),pol_int)

             calc = 1


       elseif (((extra_ferm_pair_nf.eqv..true.) .and. (WWqqqq .eqv. .false.)) .and. &
            & ((fl(1) == 'top'.or.fl(1) == 'bot'.or.&
               &fl(1) == 'str').and.nq == 4)) then 

             iq = 1
             ig = 0

          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)
          pos3(1:3) = pos4(2:4) 
          
          ng1 = pos3(1) -2
          ng2 = pos3(2) - pos3(1) -1
          ng3 = pos3(3) - pos3(2) - 1
          vs=f_bffbf_3(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
          calc = 1


          !--------------------------------------------------

          elseif (nq == 4 .and. fl(1) /= 'glu' .and. WWqqqq .eqv. .true.) then
             iq = 1
             ig = 0
             
             call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                  &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)
             pos3(1:3) = pos4(2:4)

             if (((case_a2 .eqv. .true.) .or. (case_a4 .eqv. .true.)).and. &
                  &((fl(1) == 'top').or.(fl(1) == 'bot').or.&
                  & (fl(1) == 'chr')) .and. (fl1(2) /= 'str' .and. &
                  & fl1(3) /= 'str')) then 

                pos3(1:3) = pos4(2:4) 
                ng1 = pos3(1) -2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                
                vs=f_bffbf_2(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),&
                     &pol_int)

                calc = 1
              
           
             elseif((case_a1 .eqv. .true.) .or. (case_b1 .eqv. .true.) .or. (case_b2 .eqv. .true.) .or. (case_a3 .eqv. .true.)) then
                
                pos3(1:3) = pos4(2:4) 
                ng1 = pos3(1) -2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                
                vs=f_bffbf_2(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
                
                calc = 1
                
             elseif (((case_a4 .eqv. .true.) .or. (case_a2 .eqv. .true.)) .and. (fl1(2) == 'str' .and. fl1(3) == 'str')) then
                

                pos3(1:3) = pos4(2:4) 
                ng1 = pos3(1) -2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                
                vs=f_bffbf_3(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                   &fl1(2:nq),fl(1),ng1,ng2,ng3,giarray(1:ng),qiarray(2:4),pol_int)
               
                calc = 1

             endif

           
       elseif (fl(1) == 'glu'.and.nq == 4) then

          iq = 0
          ig = 0

          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos4,ig,iq,giarray,qiarray)

          ng1 = pos4(1) -2
          ng2 = pos4(2) - pos4(1) - 1
          ng3 = pos4(3) - pos4(2) - 1
          ng4 = pos4(4) - pos4(3) - 1

          pol_glu(:,ng) = sp(1:Dv,1)

          vg=g_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
               &pol_q,p_q,fl1,ng1,ng2,ng3,ng4,giarray(1:ng-1),qiarray,pol_int)

          calc = 1


          !--------------------------------------------------------------------
       elseif (((fl(1) == 'bot') .or. (fl(1) == 'top') .or. &
            &(fl(1) == 'chr')) .and.nq == 6) then 

          iq = 1
          ig = 0

          call split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos6,ig,iq,giarray,qiarray)

          ng1 = pos6(2) -2
          ng2 = pos6(3) - pos6(2) -1
          ng3 = pos6(4) - pos6(3) - 1
          ng4 = pos6(5) - pos6(4) - 1
          ng5 = pos6(6) - pos6(5) - 1

          
          vs=f_bffbffbf(pol_glu,p_glu,pol_q(:,2:nq),&
               &p_q(:,2:nq),fl1(2:nq),fl(1),ng1,ng2,ng3,ng4,ng5,&
               &giarray,qiarray(2:nq),pol_int)

          calc = 1

       else
          print *, 'UNCALCULATED';
          print *, 'nw=',nw,'nq=',nq,'ng=',ng
          print *, fl
          write(*,*)  (fl(1) == 'top'.or.fl(1) == 'bot'.or.&
               &fl(1) == 'str'.or.fl(1) == 'chr')
          write(*,*) nq == 4

          stop 'calc_ampl: nw =0: uncalculated'


       endif

       !-----------------------------  nw = 1 case

    elseif (nw == 1) then 

       if (((fl(1) == 'top') .or. (fl(1) == 'bot') .or. (fl(1) == 'chr')) .and.nq == 2) then 

          iq = 0
          ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos1,ig,iq,iW,giarray,qiarray,&
       &wiarray(1:1),posW=posW)

          if (ferm_loops_Z .or. ferm_loops_Z_sbs) then  ! tree level 
             if (pos1(1) > posW) then 
                ng1 = pos1(1) - 3   
             else
                ng1 = pos1(1) - 2   
             endif
          else
             ng1 = pos1(1) - 3   
          endif

          vs=fW(pol_glu,p_glu,pol_q(:,1),p_q(:,1),fl1(1),fl(1),&
               &pol_W,p_W,ng1,giarray,qiarray(1:1),wiarray(1),pol_int)
          calc = 1

       elseif ((nq == 4) .and. (fl(1) /= 'glu') .and. (WWqqqq .eqv. .true.)) then
          iq = 1
          ig = 0
          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                  &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
                  &giarray,qiarray,wiarray(1:1))
          isstr = .false.
          do i=1,size(fl)
             if (fl(i) == 'str') isstr = .true.
          enddo    

          if (((isstr .eqv. .false.).and. (fl(1) == 'top' .or. fl(1) == 'bot')) &
               &.or. (((case_a1 .eqv. .true.) .or. (case_b1 .eqv. .true.)) .and. (fl(1) /= 'glu'))) then 
             pos3(1:3) = pos4(2:4) 
             
             if (iW.lt.pos3(1)) then 
                ng1 = pos3(1) -3
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 1
             endif

             if ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -2
                ng3 = pos3(3) - pos3(2) - 1
                sw = 2
             endif

             if ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 2
                sw = 3
             endif

             if (iW.gt.pos3(3)) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 4
             endif
       
             vs=fW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),4,pol_int)


             calc = 1

          
          elseif ((isstr .eqv. .true.) .and. (fl(1) /= 'glu') &
               & .and. (fl1(3) /= 'str' .and. fl1(2) /= 'str')) then
      
             if (iW.lt.pos3(1)) then 
                ng1 = pos3(1) -3
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
             endif

             if ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -2
                ng3 = pos3(3) - pos3(2) - 1
             endif

             if ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 2
             endif

             if (iW.gt.pos3(3)) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
             endif
       
             vs=fW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),4,pol_int)


             calc = 1

          elseif ((isstr .eqv. .true.) .and. (fl(1) /= 'glu') &
               & .and. ((case_a4 .eqv. .true.) .or. (case_a2 .eqv. .true.)) .and. &
               & (fl1(3) == 'str') .and. (fl1(2) == 'str')) then
             
             pos3(1:3) = pos4(2:4) 
             
             if (iW.lt.pos3(1)) then 
                ng1 = pos3(1) -3
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 1
             endif

             if ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -2
                ng3 = pos3(3) - pos3(2) - 1
                sw = 2
             endif

             if ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 2
                sw = 3
             endif

             if (iW.gt.pos3(3)) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 4
             endif

             vs=fW_bffbf_2(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                     &giarray,qiarray(2:nq),4,pol_int)
             calc = 1

          elseif (case_a3) then
              pos3(1:3) = pos4(2:4) 
             
             if (iW.lt.pos3(1)) then 
                ng1 = pos3(1) -3
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 1
             endif

             if ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -2
                ng3 = pos3(3) - pos3(2) - 1
                sw = 2
             endif

             if ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 2
                sw = 3
             endif

             if (iW.gt.pos3(3)) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 4
             endif
             ! Swap quarks
           
             vs=fW_bffbf_1(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,&
                     &giarray,qiarray(2:nq),4,pol_int)
             calc = 1

             
          endif
        
       elseif ((extra_ferm_pair1.eqv. .false.) .and. (WWqqqq .eqv. .false.) &
            & .and. ((fl(1) == 'top'.or.fl(1) == 'bot'.or. fl(1) == 'str') &
               &.and.nq == 4) .and. .not.(ferm_loops .and. fl(1) == 'str') ) then 

             iq = 1
             ig = 0

             call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                  &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
                  &giarray,qiarray,wiarray(1:1))
             
             pos3(1:3) = pos4(2:4) 


             if (iW.lt.pos3(1)) then 
                ng1 = pos3(1) -3
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 1
             endif

             if ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -2
                ng3 = pos3(3) - pos3(2) - 1
                sw = 2
             endif


             if ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 2
                sw = 3
             endif

             if (iW.gt.pos3(3)) then 
                ng1 = pos3(1) - 2
                ng2 = pos3(2) - pos3(1) -1
                ng3 = pos3(3) - pos3(2) - 1
                sw = 4
             endif
             
             ! -- take into account that is flavour is different that Z must be 
             !    emitted from both fermions in the loop, so set sw = 4
             if (ferm_loops_Z_sbs .and. &
                  &((flampl(1)  == 'ch1' .and. flampl(nq+ng+nw) == 'ch2') .or. &
                  & (flampl(1)  == 'ch2' .and. flampl(nq+ng+nw) == 'ch1'))) then 
                sw = 4
             endif

          if (qbq_WW_and_gluons .and. (.not.WWqqqq)) then 
             vs = fW_bffbf_2W(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),wiarray(1),pol_int)
          else            
             vs=fW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),wiarray(1),pol_int)
          endif

             calc = 1

       elseif ((extra_ferm_pair1.eqv. .true.) .and. (WWqqqq .eqv. .false.) .and.& 
          & ((fl(1) == 'top'.or.fl(1) == 'bot') .and. nq == 4)) then 

             detstr = 0

             do i=1,np
                if (fl(i)=='str') detstr = 1
             enddo

             if (detstr == 0 ) then 

                iq = 1
                ig = 0

             call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                  &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))
             pos3(1:3) = pos4(2:4) 

                if (iW.lt.pos3(1)) then 
                   ng1 = pos3(1) -3
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 1
                elseif ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -2
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 2
                elseif ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 2
                   sw = 3
                elseif (iW.gt.pos3(3)) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 4
                endif

                vs=fW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
                     &giarray,qiarray(2:nq),wiarray(1),pol_int)

                calc = 1

             else  ! this means that there is an sbs pair


                iq = 1
                ig = 0

             call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                  &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
                  &giarray,qiarray,wiarray(1:1))
             pos3(1:3) = pos4(2:4) 


                if (iW.lt.pos3(1)) then 
                   ng1 = pos3(1) -3
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 1
                elseif ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -2
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 2
                elseif ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 2
                   sw = 3
                elseif (iW.gt.pos3(3)) then 
                   ng1 = pos3(1) - 2
                   ng2 = pos3(2) - pos3(1) -1
                   ng3 = pos3(3) - pos3(2) - 1
                   sw = 4
                endif

                !-------no sw here 

                vs=fW_bffbf_1(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                     &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,&
                     &giarray,qiarray(2:nq),wiarray(1),pol_int)


                calc = 1

             endif

       !----------------------------------------------
       elseif (fl(1) == 'glu'.and.nq == 2) then

          iq = 0
          ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos2,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1),nostr=.true.)

          ng1 = pos2(1) -2
          ng2 = pos2(2) - pos2(1) - 2

          pol_glu(:,ng) = sp(1:Dv,1)
          
          vg=gW_fbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
               &pol_q(:,1),p_q(:,1),fl1(1),&
               &pol_q(:,2),p_q(:,2),fl1(2),pol_W,p_W,ng1,ng2,&
               &giarray(1:ng-1),qiarray(1:2),wiarray(1),pol_int)
          calc = 1


       !----------------------------------------------- 
       elseif ((extra_ferm_pair1 .eqv. .false.) .and. (WWqqqq .eqv. .false.)&
            &.and.(extra_ferm_pair2 .eqv. .false.) .and. &
            &(fl(1) == 'glu'.and.nq == 4)) then

             iq = 0
             ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))

             if (iW.lt.pos4(1)) then 
                ng1 = pos4(1) -3
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 1
             elseif ((iW.gt.pos4(1)).and.(iW.lt.pos4(2))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -2
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 2
             elseif ((iW.gt.pos4(2)).and.(iW.lt.pos4(3))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 2
                ng4 = pos4(4) - pos4(3) - 1
                sw = 3
             elseif ((iW.gt.pos4(3)).and.(iW.lt.pos4(4))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 2
                sw = 4
             elseif (iW.gt.pos4(4)) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 5
             endif

             vg=gW_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                  &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,sW,&
                  &giarray(1:ng-1),qiarray,wiarray(1),pol_int)

             calc = 1

          elseif ((WWqqqq .eqv. .true.) .and.(fl(1) == 'glu').and.(nq == 4)) then
           
           iq = 0
           ig = 0

           call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
                &giarray,qiarray,wiarray(1:1))

           if (iW.lt.pos4(1)) then 
              ng1 = pos4(1) -3
              ng2 = pos4(2) - pos4(1) -1
              ng3 = pos4(3) - pos4(2) - 1
              ng4 = pos4(4) - pos4(3) - 1
              sw = 1
           elseif ((iW.gt.pos4(1)).and.(iW.lt.pos4(2))) then 
              ng1 = pos4(1) - 2
              ng2 = pos4(2) - pos4(1) -2
              ng3 = pos4(3) - pos4(2) - 1
              ng4 = pos4(4) - pos4(3) - 1
              sw = 2
           elseif ((iW.gt.pos4(2)).and.(iW.lt.pos4(3))) then 
              ng1 = pos4(1) - 2
              ng2 = pos4(2) - pos4(1) -1
              ng3 = pos4(3) - pos4(2) - 2
              ng4 = pos4(4) - pos4(3) - 1
              sw = 3
           elseif ((iW.gt.pos4(3)).and.(iW.lt.pos4(4))) then 
              ng1 = pos4(1) - 2
              ng2 = pos4(2) - pos4(1) -1
              ng3 = pos4(3) - pos4(2) - 1
              ng4 = pos4(4) - pos4(3) - 2
              sw = 4
           elseif (iW.gt.pos4(4)) then 
              ng1 = pos4(1) - 2
              ng2 = pos4(2) - pos4(1) -1
              ng3 = pos4(3) - pos4(2) - 1
              ng4 = pos4(4) - pos4(3) - 1
              sw = 5
           endif

           if (fl1(2) /= 'str' .or. fl1(3) /= 'str') then 

              vg=gW_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                   &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,sw,&
                   &giarray(1:ng-1),qiarray,4,pol_int)

              calc = 1
           elseif (fl1(2) == 'str' .and. fl1(3) == 'str') then
                vg=gW_sbsfbf_3(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                   &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,sw,&
                   &giarray(1:ng-1),qiarray,4,pol_int)

             
              calc = 1
           endif

        elseif ((extra_ferm_pair1.eqv. .true.) .and. (WWqqqq .eqv. .false.) .and. &
            &(fl(1) == 'glu'.and.nq == 4)) then

             iq = 0
             ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))


             if (iW.lt.pos4(1)) then 
                ng1 = pos4(1) -3
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 1
             elseif ((iW.gt.pos4(1)).and.(iW.lt.pos4(2))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -2
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 2
             elseif ((iW.gt.pos4(2)).and.(iW.lt.pos4(3))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 2
                ng4 = pos4(4) - pos4(3) - 1
                sw = 3
             elseif ((iW.gt.pos4(3)).and.(iW.lt.pos4(4))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 2
                sw = 4
             elseif (iW.gt.pos4(4)) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 5
             endif


             if ( (fl(1).eq.'str'.and.fl(2).eq.'str')&
                  &.or.(fl(3).eq.'str'.and.fl(4).eq.'str') ) then 
                vg=gW_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                     &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,sW,&
                     &giarray(1:ng-1),qiarray,wiarray(1),pol_int)
             else
                vg=gW_sbsfbf_1(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                     &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,&
                     &giarray(1:ng-1),qiarray,wiarray(1),pol_int)
             endif

             calc = 1

       elseif ((extra_ferm_pair2.eqv. .true.) .and. (WWqqqq .eqv. .false.) .and. &
            &(fl(1) == 'glu'.and.nq == 4)) then

             iq = 0
             ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))


             if (iW.lt.pos4(1)) then 
                ng1 = pos4(1) -3
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 1
             elseif ((iW.gt.pos4(1)).and.(iW.lt.pos4(2))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -2
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 2
             elseif ((iW.gt.pos4(2)).and.(iW.lt.pos4(3))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 2
                ng4 = pos4(4) - pos4(3) - 1
                sw = 3
             elseif ((iW.gt.pos4(3)).and.(iW.lt.pos4(4))) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 2
                sw = 4
             elseif (iW.gt.pos4(4)) then 
                ng1 = pos4(1) - 2
                ng2 = pos4(2) - pos4(1) -1
                ng3 = pos4(3) - pos4(2) - 1
                ng4 = pos4(4) - pos4(3) - 1
                sw = 5
             endif

             if ( (fl1(1).eq.'str'.and.fl1(2).eq.'str')&
                  &.or.(fl1(3).eq.'str'.and.fl1(4).eq.'str') ) then 
                vg=gW_sbsfbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                     &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,sW,&
                     &giarray(1:ng-1),qiarray,wiarray(1),pol_int)
             else
                if ( (fl1(2).eq.'str'.and.fl1(3).eq.'str') ) then 
                   vg=gW_sbsfbf_1(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                        &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,&
                        &giarray(1:ng-1),qiarray,wiarray(1),pol_int)
                else 
                   if ( (fl1(1).eq.'str'.and.fl1(4).eq.'str') ) then 
                      vg=gW_sbsfbf_2(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
                           &pol_q,p_q,fl1,pol_W,p_W,ng1,ng2,ng3,ng4,&
                           &giarray(1:ng-1),qiarray,wiarray(1),pol_int)
                   endif
                endif
             endif
             calc = 1

       !-------- new staff; only relevant for fermion loops to qbqWggg
       elseif (((fl(1) == 'str').and.nq == 4.and.ferm_loops).and. &
            (WWqqqq .eqv. .false.)) then 

          iq = 1
          ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))
          pos3(1:3) = pos4(2:4) 

          if (iW.lt.pos3(1)) then 
             ng1 = pos3(1) -3
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 1
          elseif ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -2
             ng3 = pos3(3) - pos3(2) - 1
             sw = 2
          elseif ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 2
             sw = 3
          elseif (iW.gt.pos3(3)) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 4
          endif

          !-------note that f's and bf's are in the strange order          
          vs=fsW_fbfbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
               &giarray,qiarray(2:nq),wiarray(1),pol_int)

          calc = 1

       elseif ((fl(1) == 'chr') .and. nq == 4 .and. &
            &(ferm_loops_Z)) then 


          iq = 1
          ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))
          pos3(1:3) = pos4(2:4) 

          if (iW.lt.pos3(1)) then 
             ng1 = pos3(1) -3
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 1
          elseif ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -2
             ng3 = pos3(3) - pos3(2) - 1
             sw = 2
          elseif ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 2
             sw = 3
          elseif (iW.gt.pos3(3)) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 4
          endif

          vs=fW_fbfbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
               &giarray,qiarray(2:nq),wiarray(1),pol_int)

          calc = 1

       elseif ((fl(1) == 'chr') .and. nq == 4 .and. &
            &(ferm_loops_Z_sbs)) then 


          iq = 1
          ig = 0

          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos4,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))
          pos3(1:3) = pos4(2:4) 

          if (iW.lt.pos3(1)) then 
             ng1 = pos3(1) -3
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 1
          elseif ((iW.gt.pos3(1)).and.(iW.lt.pos3(2))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -2
             ng3 = pos3(3) - pos3(2) - 1
             sw = 2
          elseif ((iW.gt.pos3(2)).and.(iW.lt.pos3(3))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 2
             sw = 3
          elseif (iW.gt.pos3(3)) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 4
          endif

          vs=fW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),pol_W,p_W,ng1,ng2,ng3,sw,&
               &giarray,qiarray(2:nq),wiarray(1),pol_int)

          calc = 1



!-------------------------------------
      
       elseif ((fl(1) == 'top' .or. fl(1) == 'bot') .and. nq == 6) then
          
          iq = 1
          ig = 0
          
          call split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos6,ig,iq,iW,&
               &giarray,qiarray,wiarray(1:1))

          pos5(1:5) = pos6(2:6) 
     
          if (iW.lt.pos5(1)) then 
             ng1 = pos5(1) -3
             ng2 = pos5(2) - pos5(1) -1
             ng3 = pos5(3) - pos5(2) - 1
             ng4 = pos5(4)-pos5(3)-1
             ng5 = pos5(5)-pos5(4)-1
             sw = 1
      
          elseif ((iW.gt.pos5(2)).and.(iW.lt.pos5(3))) then 
             
             ng1 = pos5(1) -2
             ng2 = pos5(2) - pos5(1) -1
             ng3 = pos5(3) - pos5(2) - 2
             ng4 = pos5(4)-pos5(3)-1
             ng5 = pos5(5)-pos5(4)-1
             sw = 2
             
          elseif ((iW > pos5(1)) .and. (iW < pos5(2))) then
             sw = 4
             ng1 = pos5(1) -2
             ng2 = pos5(2)-pos5(1)-2
             ng3 = pos5(3)-pos5(2)-1
             ng4 = pos5(4)-pos5(3)-1
             ng5 = pos5(5)-pos5(4)-1
             if (case_b2 .eqv. .false.) stop 'not case b2!'

          elseif ((iW.gt.pos5(4)).and.(iW.lt.pos5(5))) then 
             
             ng1 = pos5(1) -2
             ng2 = pos5(2) - pos5(1) -1
             ng3 = pos5(3) - pos5(2) - 1
             ng4 = pos5(4)-pos5(3)-1
             ng5 = pos5(5)-pos5(4)-2
             sw = 3
          endif
      
          vs=fW_bffbffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
               &fl1(2:nq),fl(1),ng1,ng2,ng3,ng4, ng5,sw,pol_W,p_W,&
               &giarray,qiarray(2:nq),iW,pol_int)

          
          calc = 1
          
       
          
!-----------------------------        


          
       else
          print *, 'UNCALCULATED',calc
          print *, 'nw=',nw,'nq=',nq,'ng=',ng
          print *, fl
          stop 'calc_ampl: Uncalculated: nw=1'
       endif
       !-------2W's
    elseif (nw ==2) then

       if (((fl(1)== 'top').or.(fl(1)=='bot')) .and. (nq==2)) then

          iq = 0
	  ig = 0

	  call split_WW_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_Wp,p_Wp,pol_Wm,p_Wm,fl1,pos1,ig,iq,iWp,iWm,&
               &giarray,qiarray,wiarray)

	  ng1 = pos1(1) - 4

          vs=fWW(pol_glu,p_glu,pol_q(:,1),p_q(:,1),&
               &pol_Wp,p_Wp,pol_Wm,p_Wm,ng1, &
               &giarray,qiarray(1:1),wiarray(1),wiarray(2),pol_int)


          calc = 1


       elseif (fl(1) == 'glu'.and.nq == 2) then

          iq = 0
          ig = 0

          call split_WW_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_Wp,p_Wp,pol_Wm,p_Wm,fl1,pos2,ig,iq,iWp,iWm,&
               &giarray,qiarray,wiarray)

          ng1 = pos2(1) -2
          ng2 = pos2(2) - pos2(1) - 3

          pol_glu(:,ng) = sp(1:Dv,1)

          vg=gWW_fbf(pol_glu(:,1:ng-1),p_glu(:,1:ng-1),&
               &pol_q(:,1),p_q(:,1),fl1(1),&
               &pol_q(:,2),p_q(:,2),fl1(2),pol_Wp,p_Wp,pol_Wm,p_Wm,ng1,ng2,&
               &giarray(1:ng-1),qiarray(1:2),wiarray(1),wiarray(2),pol_int)


          calc = 1


      elseif ((nq == 4) .and. (WWqqqq .eqv. .true.)) then
              iq = 1
              ig = 0
              call split_WW_quarks_and_gluons(Dv,sp,p,fl,iarray,&
                   &pol_glu,p_glu,pol_q,p_q,pol_WW(:,1),p_WW(:,1),&
                   &pol_WW(:,2),p_WW(:,2),fl1,pos4,ig,iq,iWW(1),iWW(2),&
                   &giarray,qiarray,wiarray)

              isstr = .false.
              do i=1,size(fl)
                 if (fl(i) == 'str') isstr = .true.
              enddo
              if (((isstr .eqv. .false.).and.&
                   & (fl(1) == 'top'.or.fl(1) == 'bot')) &
                   &.or. (((case_b1 .eqv. .true.) .or. &
                   &(case_a1 .eqv. .true.)) .and. (fl(1) /= 'glu'))) then         
                 pos3(1:3) = pos4(2:4)
                 
             if ((iWW(2) < pos3(1)) .and. (iWW(1) < pos3(1))) then
                sw = 1
                ng1 = pos3(1)-4
                ng2 = pos3(2) - pos3(1) - 1
                ng3 = pos3(3)-pos3(2) -1
             elseif ((iWW(2) < pos3(1)) .and. (iWW(1) > pos3(2)) &
                  &.and. (iWW(1) < pos3(3))) then
                sw = 2
                ng1 = pos3(1)-3
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-2
             
             elseif ((iWW(1) < pos3(1)) .and. (iWW(2) > pos3(2)) &
                  &.and. (iWW(2) <  pos3(3))) then
                sw = 2
                ng1 = pos3(1)-3
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-2
                holder(:,1) = pol_WW(:,1)
                pol_WW(:,1) = pol_WW(:,2)
                pol_WW(:,2) = holder(:,1)
                holder(:,1) = p_WW(:,1)
                p_WW(:,1) = p_WW(:,2)
                p_WW(:,2) = holder(:,1)

             elseif ((iWW(1) < pos3(1)) .and. (iWW(2) < pos3(2))) then
                sw = 2
                ng1 = pos3(1) - 3
                ng2 = pos3(2)-pos3(1)-2
                ng3 = pos3(3)-pos3(2)-1
                 holder(:,1) = pol_WW(:,1)
                pol_WW(:,1) = pol_WW(:,2)
                pol_WW(:,2) = holder(:,1)
                holder(:,1) = p_WW(:,1)
                p_WW(:,1) = p_WW(:,2)
                p_WW(:,2) = holder(:,1)

             elseif ((iWW(2) > pos3(2)) .and. (iWW(2) < pos3(3)) &
                  &.and. (iWW(1) >  pos3(2)) .and. (iWW(2) < pos3(3))) then
                sw = 3
                ng1 = pos3(1)-2
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-3
           
             else
                write(*,*) 'iWW', iWW 
                write(*,*) 'pos3', pos3
                write(*,*) 'unknown sw', sw, fl
                
             endif

             vs=fWW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_WW,p_WW,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),iWW,pol_int)

     
             calc = 1
          
          elseif ((isstr .eqv. .true.) .and. (fl(1) /= 'glu') &
               & .and. ((fl1(4) == 'str') .and. (fl1(3) == 'str'))) then

             pos3(1:3) = pos4(2:4) 
             if ((iWW(2) < pos3(1)) .and. (iWW(1) < pos3(1))) then
                sw = 1
                ng1 = pos3(1)-4
                ng2 = pos3(2) - pos3(1) - 1
                ng3 = pos3(3)-pos3(2) -1
             elseif ((iWW(2) < pos3(1)) .and. (iWW(1) > pos3(2)) &
                  &.and. (iWW(1) < pos3(3))) then
                sw = 2
                ng1 = pos3(1)-3
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-2
                write(*,*) 'should have a W swap?' 
                
             elseif ((iWW(1) < pos3(1)) .and. (iWW(2) > pos3(2)) &
                  &.and. (iWW(2) < pos3(3))) then
                sw = 2
                ng1 = pos3(1)-3
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-2
             elseif ((iWW(2) > pos3(2)) .and. (iWW(2) < pos3(3)) &
                  &.and. (iWW(1) >  pos3(2)) .and. (iWW(2) < pos3(3))) then
                sw = 3
                ng1 = pos3(1)-2
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-3
        
             elseif ((iWW(2) < iWW(1)) .and. (iWW(2) > pos3(3))) then
                
                sw = 4
                ng1 = pos3(1) - 2
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-1
                write(*,*) 'sw = 4'
             elseif ((iWW(2) > pos3(1)) .and. (iWW(2) < pos3(2))) then
                sw = 5
                ng1 = pos3(1) - 2
                ng2 = pos3(2)-pos3(1)-3
                ng3 = pos3(3)-pos3(2)-1
             
             else
                write(*,*) 'unknown sw ', fl
                write(*,*) iWW, pos3
                
             endif
             vs = fWW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_WW,p_WW,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),iWW,pol_int)
             
             calc = 1
             
          elseif (isstr .eqv. .true. .and. fl1(1) /= 'glu' .and. &
               & ((fl1(2) == 'str' .and. fl1(3) == 'str') .or. &
               & (fl1(4) == 'str' .and. fl(1) == 'str' ))) then
             
             pos3(1:3) = pos4(2:4)

             if (iWW(1) < pos3(1) .and. iWW(2) < pos3(1)) then
                sw = 1
                ng1 = pos3(1)-4
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-1
               
             elseif  (iWW(1) > pos3(1) .and. iWW(2) < pos3(2)) then
                sw = 2
                ng1 = pos3(1)-2
                ng2 = pos3(2)-pos3(1)-3
                ng3 = pos3(3)-pos3(2)-1
                holder(:,1) = pol_WW(:,1)
                pol_WW(:,1) = pol_WW(:,2)
                pol_WW(:,2) = holder(:,1)
                holder(:,1) = p_WW(:,1)
                p_WW(:,1) = p_WW(:,2)
                p_WW(:,2) = holder(:,1)
             elseif (iWW(1) < pos3(1) .and. iWW(2) > pos3(2) &
                  &.and. iWW(2) < pos3(3)) then
                sw = 3
                ng1 = pos3(1)-3
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-2
                
                
             elseif (iWW(1) > pos3(3)) then
                holder(:,1) = pol_WW(:,1)
                pol_WW(:,1) = pol_WW(:,2)
                pol_WW(:,2) = holder(:,1)
                holder(:,1) = p_WW(:,1)
                p_WW(:,1) = p_WW(:,2)
                p_WW(:,2) = holder(:,1)
                sw = 4
                ng1 = pos3(1)-2
                ng2 = pos3(2)-pos3(1)-1
                ng3 = pos3(3)-pos3(2)-1
             else
                write(*,*) 'unknown sw ', fl
                write(*,*) iWW
                write(*,*) pos3
             endif

             vs = fWW_bffbf_1(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_WW,p_WW,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),iWW,pol_int)

             calc = 1
          endif


! -------------------------------------------- RRR 
       elseif ((extra_ferm_pair1.eqv. .false.) .and. &
            & ((fl(1) == 'top'.or.fl(1) == 'bot'.or. fl(1) == 'str') &
            &.and.nq == 4) .and. .not.(ferm_loops .and. fl(1) == 'str') ) then 

          iq = 1
          ig = 0

          call split_WW_quarks_and_gluons(Dv,sp,p,fl,iarray,&
               &pol_glu,p_glu,pol_q,p_q,pol_Wp,p_Wp,pol_Wm,p_Wm,&
               & fl1,pos4,ig,iq,iWp,iWm,giarray,qiarray,wiarray)

          pos3(1:3) = pos4(2:4) 

          sw = 99

          if ((iWm.lt.pos3(1)).and.(iWp.gt.pos3(2))) then 
             ng1 = pos3(1) -3
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 2
             sw = 2
          endif

          if ((iWp.lt.pos3(1)).and.(iWm.lt.pos3(1))) then 
             ng1 = pos3(1) - 4
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 1
             sw = 1
          endif


          if ((iWm.gt.pos3(2)).and.(iWp.gt.pos3(2))) then 
             ng1 = pos3(1) - 2
             ng2 = pos3(2) - pos3(1) -1
             ng3 = pos3(3) - pos3(2) - 3
             sw = 3
          endif

          if (sw == 99) stop 'sw in nW=2 fl(1)=top nq=4 not working'

          if (qbq_WW_and_gluons) then 
             vs = fWW_bffbf(pol_glu,p_glu,pol_q(:,2:nq),p_q(:,2:nq),&
                  &fl1(2:nq),fl(1),pol_Wp,p_Wp,pol_Wm,p_Wm,ng1,ng2,ng3,sw,&
                  &giarray,qiarray(2:nq),wiarray(1),wiarray(2),pol_int) 
          else
             stop 'calling fWW_bffbf without qbq_WW_and_gluons'
             write(*,*)'NOT QQBQ_WW IN dpAMPL'

          endif


          calc = 1

       endif

    else
       stop 'calc_ampl: nw > 2 ?' 
    endif  ! endif for nw = 2 case   
    
    ! -- calc should be redundant, but keep it for now.... 
    if (calc == 0) then 
       print *, 'UNCALCULATED';
       print *, 'nw=',nw,'nq=',nq,'ng=',ng
       print *, fl(1)
       stop 
    endif

    
  end subroutine calc_ampl


  subroutine split_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,fl1,pos,ig,iq,giarray,qiarray)
    integer, intent(in)       :: Dv
    complex(dp), intent(in) :: sp(:,:), p(:,:)
    character, intent(in)     :: fl(:)*3
    integer, intent(in)       :: iarray(:)
    complex(dp), intent(out):: pol_glu(:,:),p_glu(:,:)
    complex(dp), intent(out):: pol_q(:,:),p_q(:,:)
    character, intent(out)    :: fl1(:)*3
    integer, intent(out)      :: pos(:)
    integer, intent(inout)    :: ig,iq 
    integer, intent(out)      :: giarray(:),qiarray(:)
    ! ------------------------------------------
    integer :: i


    do i=2,size(p,dim=2) 
       if (fl(i) == 'glu') then
          ig = ig + 1
          pol_glu(:,ig) = sp(1:Dv,i)
          p_glu(:,ig) = p(:,i)
          giarray(ig) = iarray(i)
       elseif (fl(i) == 'top'.or.fl(i) == 'bot'&
            &.or.fl(i)=='str'.or.fl(i) == 'chr' .or. fl(i) == 'dwn') then 
          iq = iq + 1
          pol_q(:,iq) = sp(:,i)
          p_q(:,iq) = p(:,i)
          pos(iq) = i
          fl1(iq) = fl(i)
          qiarray(iq) = iarray(i)
       else
          write(*,*) 'flavour is', fl(i)
          stop 'split_quarks_and_gluons: undefined flavour' 
       endif
    end do
    
  end subroutine split_quarks_and_gluons

  subroutine split_W_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,pol_W,p_W,fl1,pos,ig,iq,iW,giarray,qiarray,Wiarray,&
       &nostr,posW)
    integer, intent(in)       :: Dv
    complex(dp), intent(in) :: sp(:,:), p(:,:)
    character, intent(in)     :: fl(:)*3
    integer, intent(in)       :: iarray(:)
    logical, optional, intent(in) :: nostr 
    integer, optional, intent(out) :: posW
    complex(dp), intent(out) :: pol_glu(:,:),p_glu(:,:)
    complex(dp), intent(out) :: pol_q(:,:),p_q(:,:)
    complex(dp), intent(out) :: pol_W(:),p_W(:)
    character, intent(out)     :: fl1(:)*3
    integer, intent(out)       :: pos(:)
    integer, intent(inout)     :: ig,iq 
    integer, intent(out)       :: iW,giarray(:),qiarray(:),Wiarray(:)
    ! ------------------------------------------
    integer :: i

    do i=2,size(p,dim=2) 
       if (fl(i) == 'glu') then
          ig = ig + 1
          pol_glu(:,ig) = sp(1:Dv,i)
          p_glu(:,ig) = p(:,i)
          giarray(ig) = iarray(i)
       elseif (fl(i) == 'top'.or.fl(i) == 'bot'&
            &.or.fl(i)=='str' .or. fl(i) == 'chr' .or. fl(i) == 'dwn') then 
          iq = iq + 1
          pol_q(:,iq) = sp(:,i)
          p_q(:,iq) = p(:,i)
          pos(iq) = i
          fl1(iq) = fl(i)
          qiarray(iq) = iarray(i)
          if (present(nostr)) then 
             if (fl(i) =='str') stop &
                  &'split_quarks_and_gluons: strange quark not allowed here' 
          endif
       elseif (fl(i) == 'www'.or.fl(i)=='wwp'.or.fl(i)=='wwm') then 
          pol_W(:) = sp(1:Dv,i)
          p_W(:) = p(:,i)
          iW = i 
          wiarray(1) = iarray(i)
          if (present(posW)) posW = i 
       else
          write(*,*) 'flavour is', fl(i)
          stop 'split_W_quarks_and_gluons: undefined flavour' 
       endif
    end do
    
  end subroutine split_W_quarks_and_gluons

 subroutine split_WW_quarks_and_gluons(Dv,sp,p,fl,iarray,&
       &pol_glu,p_glu,pol_q,p_q,pol_Wp,p_Wp,pol_Wm,p_Wm,fl1,pos,&
       ig,iq,iWp,iWm,giarray,qiarray,wiarray,nostr)
    integer, intent(in)       :: Dv
    complex(dp), intent(in) :: sp(:,:), p(:,:)
    character, intent(in)     :: fl(:)*3
    integer, intent(in)       :: iarray(:)
    logical, optional, intent(in) :: nostr 
    complex(dp), intent(out) :: pol_glu(:,:),p_glu(:,:)
    complex(dp), intent(out) :: pol_q(:,:),p_q(:,:)
    complex(dp), intent(out) :: pol_Wp(:),p_Wp(:)
    complex(dp), intent(out) :: pol_Wm(:),p_Wm(:)
    character, intent(out)     :: fl1(:)*3
    integer, intent(out)       :: pos(:)
    integer, intent(inout)     :: ig,iq 
    integer, intent(out)       :: iWp,iWm,giarray(:),qiarray(:),wiarray(:)
    ! ------------------------------------------
    integer :: i, iw 

    iWp= 0
    iWm = 0

    iw = 0 
    do i=2,size(p,dim=2) 
       if (fl(i) == 'glu') then
          ig = ig + 1
          pol_glu(:,ig) = sp(1:Dv,i)
          p_glu(:,ig) = p(:,i)
          giarray(ig) = iarray(i)
       elseif (fl(i) == 'top'.or.fl(i) == 'bot'&
            &.or.fl(i)=='str' .or. fl(i) == 'chr') then 
          iq = iq + 1
          pol_q(:,iq) = sp(:,i)
          p_q(:,iq) = p(:,i)
          pos(iq) = i
          fl1(iq) = fl(i)
          qiarray(iq) = iarray(i)
          if (present(nostr)) then 
             if (fl(i) =='str') stop &
                  &'split_quarks_and_gluons: strange quark not allowed here' 
          endif
       elseif (fl(i) == 'wwp') then 
          iw = iw+1 
          if (iwp /= 0) then 
             iWm = i ! fill iWW(1:2) both are really Wp 
             pol_Wm(:) = sp(1:Dv,i)
             p_Wm(:) = p(:,i)
          else
             iWp = i 
             pol_Wp(:) = sp(1:Dv,i)
             p_Wp(:) = p(:,i)
          endif
          wiarray(iw) = iarray(i)
       elseif (fl(i) == 'wwm') then 
          iw = iw+1 
          if (iwm /= 0) then ! -- to be fixed for W-W- 
             pol_Wp(:) = sp(1:Dv,i)
             p_Wp(:) = p(:,i)
             iWp = i 
          else
             pol_Wm(:) = sp(1:Dv,i)
             p_Wm(:) = p(:,i)
             iWm = i 
          endif
          wiarray(iw) = iarray(i)
       else
          stop 'split_quarks_and_gluons: undefined flavour' 
       endif
    end do
    if (iw /= 2 .or. sum(wiarray) /= 12) then 
       write(*,*) 'iarray', iarray 
       write(*,*) 'iw', iw
       write(*,*) 'wiarray', wiarray 
       stop 'split_WW_quarks_and_gluons: sth failed' 
    endif
  end subroutine split_WW_quarks_and_gluons



end module dpamplitude

