!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECcut_utils.f90.

module dpcut_utils
  use types; use consts_dp 
  use warnings_and_errors 
  use dpauxiliary_functions;
  use dpmatrix_routines; use assertions 
  implicit none 
  private 

  public :: deltaM,deltavM
  public :: compute_ni, compute_pol, compute_vi, compute_bigV   
  public :: loopmom ,loopmoms 

  public :: dot, quad  

  integer, parameter :: dimmax = 4 

contains 
  
  ! deltaM(k1....km) = delta_{kv1_1..kv1_m}^{kv2_1...kv2_m}
  function deltaM(kv1,kv2,new,ind,newout) 
    complex(dp), intent(in) :: kv1(:,:),kv2(:,:)
    logical, intent(in)     :: new  
    integer, intent(in)     :: ind 
    logical, intent(in), optional     :: newout 
    !-----------------------------------------------------------------
    complex(dp), save :: delta_mat(dimmax,dimmax) != czero 
    complex(dp) :: deltaM
    integer :: i,j,M,D
    logical :: newout_lcl 

    newout_lcl = default_or_opt(.true.,newout) 

    M = size(kv1,dim=2)
    if (M > size(delta_mat,dim=1)) &
         &call wae_error('delta_mat','size delta_mat too small') 

    if (new .and. newout_lcl) then 
       do j=1,M 
          do i=1,M 
             delta_mat(i,j) = dot(kv1(:,i),kv2(:,j))
          end do
       end do
    elseif (.not. new) then 
       ! recompute only first part 
       do i=1,M 
          delta_mat(i,ind) = dot(kv1(:,i),kv2(:,ind))
       end do
    elseif (.not. newout_lcl) then 
       do i=1,M 
          delta_mat(i,ind-1) = dot(kv1(:,i),kv2(:,ind-1))
          delta_mat(i,ind)   = dot(kv1(:,i),kv2(:,ind))
       end do
    else
       call wae_error('deltaM','should never get here')
    end if

    deltaM = determinant_old(delta_mat(:M,:M))

  end function deltaM


  ! deltaM(k1....km) = delta_{kv1_1..kv1_m}^{kv2_1...kv2_m}
  function deltaMprime(kv1,kv2,new,ind,read_prev,offset) 
    complex(dp), intent(in) :: kv1(:,:),kv2(:,:)
    logical, intent(in)     :: new  
    integer, intent(in)     :: ind 
    logical, intent(in)     :: read_prev 
    integer, intent(in)     :: offset 
    !-----------------------------------------------------------------
    complex(dp), save :: delta_mat(dimmax,dimmax) != czero 
    complex(dp), save :: delta_mat_store(dimmax,dimmax) != czero 
    complex(dp) :: deltaMprime
    integer :: i,j,M,D

    M = size(kv1,dim=2)
    if (M > size(delta_mat,dim=1)) &
         &call wae_error('delta_mat','size delta_mat too small') 

    if (new .and. .not. read_prev) then 
       do j=1,M 
          do i=1,M 
             delta_mat(i,j) = dot(kv1(:,i),kv2(:,j))
             delta_mat_store(i,j) = delta_mat(i,j) 
          end do
       end do
    elseif (.not. new) then 
       ! recompute only first part 
       do i=1,M 
          delta_mat(i,ind) = dot(kv1(:,i),kv2(:,ind))
       end do
    elseif (read_prev) then 

       do i=1,M 
          do j=1,M
             if (j==1) then 
                delta_mat(i,j) = delta_mat_store(i+offset,j)
             else
                ! skip j=2 (replaced) 
                delta_mat(i,j) = delta_mat_store(i+offset,j+offset)
             end if
          end do
       end do
    else
       call wae_error('deltaMprime','should never get here')

    end if

    deltaMprime = determinant_old(delta_mat(:M,:M))

  end function deltaMprime


  ! deltaM needed to compute vectors ni k = (b_m...b1 p1..pn}
  ! deltaM(k1....km) = delta_{kv1_1..kv1_m}^{mu,kv2_1...kv2_m}
  function deltavM(kv1,kv2,read_prev,offset) 
    complex(dp), intent(in) :: kv1(:,:),kv2(:,:)
    logical, intent(in), optional :: read_prev 
    integer, intent(in), optional :: offset
    complex(dp) :: muv(size(kv1,dim=1)),kv2p(size(kv1,dim=1),size(kv2,dim=2)+1)
    complex(dp) :: deltavM(size(kv1,dim=1)) 
    integer :: i,offset_lcl
    logical :: new, read_prev_lcl  

    read_prev_lcl = default_or_opt(.false.,read_prev) 
    offset_lcl = default_or_opt(0,offset) 
    kv2p(:,2:) = kv2(:,1:)
    new = .true. 
    do i=1,size(kv1,dim=1)  
       muv = czero 
       if (i == 1) then 
          muv(i) = cone
       else
          muv(i) = -cone
       endif

       kv2p(:,1) = muv
       deltavM(i) = deltaMprime(kv1,kv2p,new,1,read_prev_lcl,offset_lcl)
       new = .false. 
    end do

  end function deltavM


!                     |p1^1 p1p2 ... p1pm |   
!                     |p1p2 p2^2 ... p2pm | 
!     Gram(p1....pm) =| .    .        .   | = delta_{p1..pm}^{p1...pm}
!                     | .    .        .   |
!                     |p1pm .........pm^2 | 
  function Gram(p)
    complex(dp), intent(in) :: p(:,:)
    complex(dp) :: Gram !, gram2 
    complex(dp) :: delta_mat(size(p,dim=2),size(p,dim=2))
    integer :: i,j,M

    M = size(p,dim=2)
    do j=1,M 
       do i=j,M 
          delta_mat(i,j) = dot(p(:,i),p(:,j))
          if (i>j) delta_mat(j,i) = delta_mat(i,j)
       end do
    end do
    gram =  determinant_old(delta_mat)

  end function Gram

  ! Given a particle momemta (in D dimensions) and it's gauge vector, computes
  ! D-2 polarization vectors 
  subroutine compute_pol(p,gv,b,pol)
    complex(dp), intent(in)       :: p(:),gv(:),b(:,:)
    complex(dp), intent(inout) :: pol(:,:)
    ! -----------------------------------------------------------
    complex(dp) :: pa(size(p,dim=1),2) 

    pa(:,1) = p(:)
    pa(:,2) = gv(:)
    call compute_nibis(pa,b,pol,computing_pol = .true.)
    pol = pol*ci

  end subroutine compute_pol

  ! Given M momenta p(1:D,1:M) in D dimensions and 
  ! (D-M) D-dimensional arbitrary vectors b(1:D,1:D-M) compute   
  ! (D-M) vectors n_i such that ni*nj=delta_ij (i,j=1...D-M) 
  ! and ni*pj = 0 (i=1...D-M, j=1,..M)
  ! NB: results independent od b_{D-M} -> only (D-M-1) arbitrary vectors 
  subroutine compute_ni(p,b,nv,computing_pol)
    complex(dp), intent(in) :: p(:,:),b(:,:)
    complex(dp), intent(inout) :: nv(:,:)
    logical, intent(in), optional :: computing_pol
    ! -------------------------------------------------------------
    complex(dp) :: kv(size(p,dim=1),size(p,dim=2)+size(b,dim=2)),tmp
    integer     :: i,j,k,D,M

    logical :: computing_pol_lcl 

    D = size(p,dim=1)
    M = size(p,dim=2)

    computing_pol_lcl = default_or_opt(.false., computing_pol)

    do i=1,D-M
       kv(:,i) = b(:,D-M-i+1)
    end do
    do i=D-M+1,D
       kv(:,i) = p(:,i-(D-M))
    end do

    nv = czero 
    do i=1,D-M
       !                    bi....b1  , p1..... .pM
       nv(:,i) = deltavM(kv(:,D-M+1-i:),kv(:,D-M+2-i:))
       tmp = sqrt(quad(nv(:,i)))
       nv(:,i) = nv(:,i)/tmp !sqrt(quad(nv(:,i)))

    end do


    

  end subroutine compute_ni

  ! Given M momenta p(1:D,1:M) in D dimensions and 
  ! (D-M) D-dimensional arbitrary vectors b(1:D,1:D-M) compute   
  ! (D-M) vectors n_i such that ni*nj=delta_ij (i,j=1...D-M) 
  ! and ni*pj = 0 (i=1...D-M, j=1,..M)
  ! NB: results independent od b_{D-M} -> only (D-M-1) arbitrary vectors 
  subroutine compute_nibis(p,b,nv,computing_pol)
    complex(dp), intent(in) :: p(:,:),b(:,:)
    complex(dp), intent(inout) :: nv(:,:)
    logical, intent(in), optional :: computing_pol
    ! -------------------------------------------------------------
    complex(dp) :: kv(size(p,dim=1),size(p,dim=2)+size(b,dim=2)),tmp
    integer     :: i,j,k,D,M

    logical :: computing_pol_lcl 
    logical :: read_prev 
    integer :: offset 

    D = size(p,dim=1)
    M = size(p,dim=2) ! M = 2

    computing_pol_lcl = default_or_opt(.false., computing_pol)

    do i=1,D-2
       kv(:,i) = b(:,D-M-i+1)
    end do
    do i=D-1,D
       kv(:,i) = p(:,i-(D-M))
    end do

    nv = czero 
    read_prev = .false. 
    offset = -1 

    do i=D-2,1,-1
       offset = offset+1

       !                    bi....b1  , p1..... .pM
       nv(:,i) = deltavM(kv(:,D-1-i:),kv(:,D-i:),read_prev,offset)

       tmp = sqrt(quad(nv(:,i)))
       nv(:,i) = nv(:,i)/tmp
       read_prev = .true. 
    end do
    


  end subroutine compute_nibis

  ! --- eq. (11) of EGK paper 
  subroutine compute_vi(ki,vi) 
    complex(dp), intent(in)  :: ki(:,:)
    complex(dp), intent(out) :: vi(:,:)
    complex(dp) :: kiprime(size(ki,dim=1),size(ki,dim=2))
    complex(dp) :: muv(size(ki,dim=1)),tmp
    integer :: i1,i2,i,j 
    logical :: new, newout


    newout = .true. 
    do i2=1,size(ki,dim=2)  
       kiprime=ki 
       new = .true. 
       do i1=1,size(ki,dim=1)  
          muv = czero 
          if (i1 == 1) then 
             muv(i1) = one ! time component 
          else
             muv(i1) = -one ! space components 
          endif

          kiprime(:,i2) = muv(:)
          vi(i1,i2) = deltaM(ki,kiprime,new,i2,newout)
          new = .false. 
          newout = .false. 
       end do
    enddo
    
       tmp = dot(ki(:,1),vi(:,1))
       do i=1,size(ki,dim=2)
       vi(:,i)=vi(:,i)/tmp 
    end do


  end subroutine compute_vi
     
  ! -- Essentially eq. 19 of EGK 
  subroutine compute_bigV(qi,vi,V)
    complex(dp), intent(in)    :: qi(:,:),vi(:,:)
    complex(dp), intent(inout) :: V(:)
    complex(dp) :: di(size(qi,dim=2))!,qisum(size(qi,dim=1))
    integer :: i,im1 

    do i=1,size(qi,dim=2)
       ! -- NO INTERNAL MASSES 
       di(i) = quad(qi(:,i)-qi(:,size(qi,dim=2)))
    end do

    V(:) = czero 
    do i=1,size(qi,dim=2)-1
       if (i>1) then 
          im1 = i-1
       else
          im1= size(qi,dim=2)
       end if
       V(:) = V(:) +(di(i)-di(im1))*vi(:,i)
    end do
    V= -half*V

  end subroutine compute_bigV



  ! -- Essentially eq. 19 of EGK 
  function loopmom(V,alphai,ni) 
    complex(dp), intent(in) :: V(:),alphai(:),ni(:,:)
    complex(dp)             :: loopmom(size(V,dim=1)),tmp
    integer :: i

    loopmom = czero 

    if (sum(alphai*alphai) ==czero) return 
    
    do i=1,size(ni,dim=2)
       loopmom(:) = loopmom(:) + alphai(i)*ni(:,i)
    end do
    
    tmp = sum(alphai*alphai)
    loopmom = sqrt(-dot(V,V)/tmp)*loopmom 
    loopmom = loopmom + V 

  end function loopmom

  ! special case when V3^2 = 0 (2 onshell legs in triangle)        
  function loopmoms(V,alphai,ni) 
    complex(dp), intent(in)    :: V(:),alphai(:),ni(:,:)
    complex(dp)                :: loopmoms(size(V,dim=1))
    integer :: i
    
    loopmoms = V 
    do i=1,size(ni,dim=2)
       loopmoms(:) = loopmoms(:) + alphai(i)*ni(:,i)
    end do

  end function loopmoms


end module dpcut_utils
