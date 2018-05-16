module qptodp_global 
  use types; use consts_dp 
  use define_ampl 
  use dpglobal 
  implicit none 
  private

  public :: qptodp_glob 



contains

  ! do not access qpglobal, pass all qpinfo as argument and assign 
  ! values to DP global 
  subroutine qptodp_glob(qptagdcut,qpmom,qphel,qpmomline,qpLab_ex, &
       &qpLab_in,qpcoeff5,qppropv5,qpcoeff4,qprefvect4,qppropv4,&
       &qpcoeff3,qprefvect3,qppropv3,qpcoeff2,qprefvect2,qppropv2,&
       &qpcoeff1,&
       &qpdcoeff5,qpdcoeff4,qpdcoeff3,qpdcoeff2,qpdcoeff1,&
       &qpmass5,qpmass4,qpmass3,qpmass2,qpmass1,noalloc)

    integer, intent(in)             :: qptagdcut(:,:)
    complex(qp), intent(in)    :: qpmom(:,:,:)
    complex(qp), intent(in)    :: qphel(:,:,:)
    complex(qp), intent(in)    :: qpmomline(:,:,:)
    character(len=3), intent(in)    :: qpLab_ex(:,:), qpLab_in(:,:)
    complex(qp), intent(in)    :: qpcoeff5(:,:),qppropv5(:,:)
    complex(qp), intent(in)    :: qpcoeff4(:,:)
    complex(qp), intent(in)    :: qprefvect4(:,:),qppropv4(:,:)
    complex(qp), intent(in)    :: qpcoeff3(:,:),qprefvect3(:,:)
    complex(qp), intent(in)    :: qppropv3(:,:)
    complex(qp), intent(in)    :: qpcoeff2(:,:),qprefvect2(:,:)
    complex(qp), intent(in)    :: qppropv2(:,:)
    complex(qp), intent(in)    :: qpcoeff1(:,:)	
    complex(qp), intent(in)    :: qpdcoeff1(:,:)	
    complex(qp), intent(in)    :: qpdcoeff2(:,:)	
    complex(qp), intent(in)    :: qpdcoeff3(:,:)	
    complex(qp), intent(in)    :: qpdcoeff4(:,:)	
    complex(qp), intent(in)    :: qpdcoeff5(:,:)	
    real(qp), intent(in)       :: qpmass5(:,:),qpmass4(:,:),qpmass3(:,:)
    real(qp), intent(in)       :: qpmass2(:,:),qpmass1(:,:)
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
       call allocate_mominfo(size(qpmom,dim=1),size(qpmom,dim=2))
       call allocate_coeffs(size(qpmass5,dim=1) ,size(qpmass4,dim=1) ,&
            &size(qpmass3,dim=1) ,size(qpmass2,dim=1) ,size(qpmass1,dim=1))
       first_time = .false.
    endif

    tagdcut(:,:)  = qptagdcut(:,:) 
    mom(:,:,:)    = (qpmom(:,:,:))
    hel(:,:,:)    = (qphel(:,:,:))   
    momline(:,:,:)= (qpmomline(:,:,:))
    Lab_ex(:,:)   = qpLab_ex(:,:)
    Lab_in(:,:)   = qpLab_in(:,:)
    coeff5(:,:)   = (qpcoeff5(:,:))  
    propv5(:,:)   = (qppropv5(:,:))  
    coeff4(:,:)   = (qpcoeff4(:,:))  
    refvect4(:,:) = (qprefvect4(:,:))
    propv4(:,:)   = (qppropv4(:,:))  
    coeff3(:,:)   = (qpcoeff3(:,:))  
    refvect3(:,:) = (qprefvect3(:,:))
    propv3(:,:)   = (qppropv3(:,:))  
    coeff2(:,:)   = (qpcoeff2(:,:))  
    refvect2(:,:) = (qprefvect2(:,:))
    propv2(:,:)   = (qppropv2(:,:))  
    coeff1(:,:)	  = (qpcoeff1(:,:))	 

    dcoeff1(:,:) = (qpdcoeff1(:,:))	 
    dcoeff2(:,:) = (qpdcoeff2(:,:))	 
    dcoeff3(:,:) = (qpdcoeff3(:,:))	 
    dcoeff4(:,:) = (qpdcoeff4(:,:))	 
    dcoeff5(:,:) = (qpdcoeff5(:,:))	 

    mass5(:,:)    = (qpmass5(:,:))   
    mass4(:,:)    = (qpmass4(:,:))   
    mass3(:,:)    = (qpmass3(:,:))   
    mass2(:,:)    = (qpmass2(:,:))   
    mass1(:,:)    = (qpmass1(:,:))   

  end subroutine qptodp_glob
  


end module qptodp_global
