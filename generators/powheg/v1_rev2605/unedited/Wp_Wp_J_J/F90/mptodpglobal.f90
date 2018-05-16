module mptodp_global 
  use types; use consts_dp 
  use mpmodule 
  use define_ampl 
  use dpglobal 
  use mpconverter
  implicit none 
  private

  public :: mptodp_glob 
contains

  ! do not access mpglobal, pass all mpinfo as argument and assign 
  ! values to DP global 
  subroutine mptodp_glob(mptagdcut,mpmom,mphel,mpmomline,mpLab_ex, &
       &mpLab_in,mpcoeff5,mppropv5,mpcoeff4,mprefvect4,mppropv4,&
       &mpcoeff3,mprefvect3,mppropv3,mpcoeff2,mprefvect2,mppropv2,&
       &mpcoeff1,&
       &mpdcoeff5,mpdcoeff4,mpdcoeff3,mpdcoeff2,mpdcoeff1,&
       &mpmass5,mpmass4,mpmass3,mpmass2,mpmass1,noalloc)

    integer, intent(in)             :: mptagdcut(:,:)
    type(mp_complex), intent(in)    :: mpmom(:,:,:)
    type(mp_complex), intent(in)    :: mphel(:,:,:)
    type(mp_complex), intent(in)    :: mpmomline(:,:,:)
    character(len=3), intent(in)    :: mpLab_ex(:,:), mpLab_in(:,:)
    type(mp_complex), intent(in)    :: mpcoeff5(:,:),mppropv5(:,:)
    type(mp_complex), intent(in)    :: mpcoeff4(:,:)
    type(mp_complex), intent(in)    :: mprefvect4(:,:),mppropv4(:,:)
    type(mp_complex), intent(in)    :: mpcoeff3(:,:),mprefvect3(:,:)
    type(mp_complex), intent(in)    :: mppropv3(:,:)
    type(mp_complex), intent(in)    :: mpcoeff2(:,:),mprefvect2(:,:)
    type(mp_complex), intent(in)    :: mppropv2(:,:)
    type(mp_complex), intent(in)    :: mpcoeff1(:,:)	
    type(mp_complex), intent(in)    :: mpdcoeff1(:,:)	
    type(mp_complex), intent(in)    :: mpdcoeff2(:,:)	
    type(mp_complex), intent(in)    :: mpdcoeff3(:,:)	
    type(mp_complex), intent(in)    :: mpdcoeff4(:,:)	
    type(mp_complex), intent(in)    :: mpdcoeff5(:,:)	
    type(mp_real), intent(in)       :: mpmass5(:,:),mpmass4(:,:),mpmass3(:,:)
    type(mp_real), intent(in)       :: mpmass2(:,:),mpmass1(:,:)
    logical, intent(in), optional   :: noalloc
    logical, save                   :: first_time = .true. 
    logical                         :: do_allocation

    
    if (present(noalloc)) then 
       if (noalloc) then 
          do_allocation = .false.
       else
          do_allocation = .true.
       endif
    else
       do_allocation = .true.
    endif
    if (first_time .and. do_allocation) then 
       call allocate_mominfo(size(mpmom,dim=1),size(mpmom,dim=2))
       call allocate_coeffs(size(mpmass5,dim=1) ,size(mpmass4,dim=1) ,&
            &size(mpmass3,dim=1) ,size(mpmass2,dim=1) ,size(mpmass1,dim=1))
       first_time = .false.
    endif
    tagdcut(:,:)  = mptagdcut(:,:) 
    mom(:,:,:)    = mptodp(mpmom(:,:,:))
    hel(:,:,:)    = mptodp(mphel(:,:,:))   
    momline(:,:,:)= mptodp(mpmomline(:,:,:))
    Lab_ex(:,:)   = mpLab_ex(:,:)
    Lab_in(:,:)   = mpLab_in(:,:)
    coeff5(:,:)   = mptodp(mpcoeff5(:,:))  
    propv5(:,:)   = mptodp(mppropv5(:,:))  
    coeff4(:,:)   = mptodp(mpcoeff4(:,:))  
    refvect4(:,:) = mptodp(mprefvect4(:,:))
    propv4(:,:)   = mptodp(mppropv4(:,:))  
    coeff3(:,:)   = mptodp(mpcoeff3(:,:))  
    refvect3(:,:) = mptodp(mprefvect3(:,:))
    propv3(:,:)   = mptodp(mppropv3(:,:))  
    coeff2(:,:)   = mptodp(mpcoeff2(:,:))  
    refvect2(:,:) = mptodp(mprefvect2(:,:))
    propv2(:,:)   = mptodp(mppropv2(:,:))  
    coeff1(:,:)	  = mptodp(mpcoeff1(:,:))	 

    dcoeff1(:,:) = mptodp(mpdcoeff1(:,:))	 
    dcoeff2(:,:) = mptodp(mpdcoeff2(:,:))	 
    dcoeff3(:,:) = mptodp(mpdcoeff3(:,:))	 
    dcoeff4(:,:) = mptodp(mpdcoeff4(:,:))	 
    dcoeff5(:,:) = mptodp(mpdcoeff5(:,:))	 

    mass5(:,:)    = mptodp(mpmass5(:,:))   
    mass4(:,:)    = mptodp(mpmass4(:,:))   
    mass3(:,:)    = mptodp(mpmass3(:,:))   
    mass2(:,:)    = mptodp(mpmass2(:,:))   
    mass1(:,:)    = mptodp(mpmass1(:,:))   

  end subroutine mptodp_glob
  


end module mptodp_global
