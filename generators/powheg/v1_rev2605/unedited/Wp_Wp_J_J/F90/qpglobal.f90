!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECglobal.f90.

! globla variables with the info about cuts, flavor, residues, essentilly everything...
module qpglobal 
  use types; use consts_qp 
  implicit none
  private 

  !----- global declarations 
  integer, allocatable         ,public   :: tagdcut(:,:)
  complex(qp), allocatable   ,public   :: mom(:,:,:)
  complex(qp), allocatable   ,public   :: mom2(:,:,:,:)
  complex(qp), allocatable   ,public   :: hel(:,:,:)
  complex(qp), allocatable   ,public   :: hel2(:,:,:,:)
  complex(qp), allocatable   ,public   :: momline(:,:,:)
  complex(qp), allocatable   ,public   :: momline2(:,:,:,:)
  integer, allocatable   ,public         :: identity(:,:)
  integer, allocatable   ,public         :: identity2(:,:,:)
  character(len=3), allocatable,public   :: Lab_ex(:,:), Lab_in(:,:)
  character(len=3), allocatable,public   :: Lab_ex2(:,:,:), Lab_in2(:,:,:)
  complex(qp), allocatable   ,public   :: coeff5(:,:),propv5(:,:)
  complex(qp), allocatable   ,public   :: coeff4(:,:)
  complex(qp), allocatable   ,public   :: refvect4(:,:),propv4(:,:)
  complex(qp), allocatable   ,public   :: coeff3(:,:),refvect3(:,:)
  complex(qp), allocatable   ,public   :: propv3(:,:)
  complex(qp), allocatable   ,public   :: coeff2(:,:),refvect2(:,:)
  complex(qp), allocatable   ,public   :: propv2(:,:)
  complex(qp), allocatable   ,public   :: coeff1(:,:)	

  complex(qp), allocatable   ,public   :: dcoeff1(:,:),dcoeff2(:,:)
  complex(qp), allocatable   ,public   :: dcoeff3(:,:),dcoeff4(:,:)
  complex(qp), allocatable   ,public   :: dcoeff5(:,:)

  real(qp), allocatable      ,public   :: mass5(:,:),mass4(:,:),mass3(:,:)
  real(qp), allocatable      ,public   :: mass2(:,:),mass1(:,:)
  real(qp)                   ,public   :: en,ptcut,etacut, eg

  integer, allocatable  , public :: Lc5(:,:), Lc4(:,:), Lc3(:,:), Lc2(:,:)
  character, allocatable, public :: F5(:,:)*3, F4(:,:)*3, F3(:,:)*3,F2(:,:)*3
  integer, allocatable  , public :: Yc5(:,:), Yc4(:,:), Yc3(:,:), Yc2(:,:)

  integer, parameter, public :: N5max = 100 !21 
  integer, parameter, public :: N4max = 100 !45
  integer, parameter, public :: N3max = 100 !50 
  integer, parameter, public :: N2max = 100 !22 
  integer, parameter, public :: tN2max = 100 !30 
  integer, parameter, public :: N1max = 1 

  ! auxiliary vectors 
  complex(qp), public                  :: bvec(4,3)  

  ! angle for MC over helicities of gluons, max of 4 gluons allowed 
  real(qp), public :: MCphi(4)

  ! -- define global external polarizations, specifi to WW2j 
  ! -- 6 external polarization, two possible helicities, compute this only once
  complex(qp), public :: MCFMpol_ext_qq(2,4,-1:1)  
  complex(qp), public :: MCFMpol_ext_gg(2,4,-1:1)  
  complex(qp), public :: MCFMpol_ext_ww(2,4,-1:1) 
  complex(qp), public :: MCFMpol_ext_z(4) 
  real(qp), public    :: MCFMhardscale, MCFMsqrthardscale 
  logical, public       :: swap_ww, swap_gg, newkinematics



  !------ end global declarations

  public :: allocate_mominfo, allocate_arrays, allocate_coeffs
  public :: deallocate_global ,deallocate_arrays 
  public :: fix_bvec,  allocate_mominfo2 
  
  contains 
    
    subroutine allocate_mominfo(ncut,npoint)
      integer, intent(in) :: ncut, npoint 
      if (.not. allocated(mom)) allocate(mom(ncut,npoint,4))
      if (.not. allocated(hel)) allocate(hel(ncut,npoint,4))
      if (.not. allocated(momline)) allocate(momline(ncut,npoint,4))
      if (.not. allocated(identity)) allocate(identity(ncut,npoint))
      if (.not. allocated(Lab_ex)) allocate(Lab_ex(ncut,npoint))
      if (.not. allocated(Lab_in)) allocate(Lab_in(ncut,npoint))

    end subroutine allocate_mominfo

   subroutine allocate_mominfo2(ncut,npoint)
      integer, intent(in) :: ncut, npoint 
      if (.not. allocated(mom2)) allocate(mom2(ncut,ncut,npoint,4))
      if (.not. allocated(hel2)) allocate(hel2(ncut,ncut,npoint,4))
      if (.not. allocated(identity2)) allocate(identity2(ncut,ncut,npoint))

      if (.not. allocated(Lab_ex2)) allocate(Lab_ex2(ncut,ncut,npoint))
      if (.not. allocated(Lab_in2)) allocate(Lab_in2(ncut,ncut,npoint))

    end subroutine allocate_mominfo2

    subroutine allocate_arrays(N5,N4,N3,N2)
      integer, intent(in) :: N5,N4,N3,N2
    
      if (.not. allocated(Lc5)) then 
         allocate(Lc5(N5,5),Lc4(N4,4),Lc3(N3,3),Lc2(N2,2))
      else
         if (size(Lc5,dim=1) < N5) stop 'allocate_arrays: size(Lc5,dim=1) < N5'
         if (size(Lc4,dim=1) < N4) stop 'allocate_arrays: size(Lc4,dim=1) < N4'
         if (size(Lc3,dim=1) < N3) stop 'allocate_arrays: size(Lc3,dim=1) < N3'
         if (size(Lc2,dim=1) < N2) stop 'allocate_arrays: size(Lc2,dim=1) < N2'
      endif

      if (.not. allocated(F5))  allocate(F5(N5,5),F4(N4,4),F3(N3,3),F2(N2,2))
      if (.not. allocated(Yc5)) allocate(Yc5(N5,1),Yc4(N4,1),Yc3(N3,1),Yc2(N2,1))
      

    end subroutine allocate_arrays

    subroutine allocate_coeffs(N5,N4,N3,N2,N1)
      integer, intent(in) :: N5,N4,N3,N2,N1
          
      if (.not. allocated(tagdcut)) then 
         allocate(tagdcut(N2,1))
         allocate(coeff5(N5,1),propv5(N5,20))
         allocate(coeff4(N4,5))
         allocate(refvect4(N4,4),propv4(N4,16))
         allocate(coeff3(N3,10),refvect3(N3,8))
         allocate(propv3(N3,12))
         allocate(coeff2(N2,10),refvect2(N2,12))
         allocate(propv2(N2,8))
         allocate(coeff1(N1,1))	
         
         allocate(dcoeff5(N5,1))
         allocate(dcoeff4(N4,5))
         allocate(dcoeff3(N3,10))
         allocate(dcoeff2(N2,10))
         allocate(dcoeff1(N1,1))	
         
         allocate(mass5(N5,5),mass4(N4,4),mass3(N3,3))
         allocate(mass2(N2,2),mass1(N1,1))
      endif

    end subroutine allocate_coeffs


    subroutine deallocate_global 

      deallocate(tagdcut)
      deallocate(mom,hel,momline,identity)
      
      deallocate(Lab_ex,Lab_in)

 
      deallocate(mom2,hel2,identity2)
      deallocate(lab_ex2,Lab_in2)
     
      deallocate(coeff5,propv5)
      deallocate(coeff4,refvect4,propv4)
      deallocate(coeff3,refvect3,propv3)
      deallocate(coeff2,refvect2,propv2)
      deallocate(coeff1)	

      deallocate(dcoeff5,dcoeff4,dcoeff3,dcoeff2,dcoeff1)	
      deallocate(mass5,mass4,mass3,mass2,mass1)

    end subroutine deallocate_global

    subroutine deallocate_arrays()
    
      deallocate(Lc5,Lc4,Lc3,Lc2)
      deallocate(F5,F4,F3,F2)
      deallocate(Yc5,Yc4,Yc3,Yc2)

    end subroutine deallocate_arrays

    subroutine fix_bvec 
      
      bvec(1,1) = 7._qp 
      bvec(2,1) = -5._qp 
      bvec(3,1) = three
      bvec(4,1) = 17._qp

      bvec(1,2) = one 
      bvec(2,2) = -two 
      bvec(3,2) = four
      bvec(4,2) = -four
      
      bvec(1,3) = one 
      bvec(2,3) = zero 
      bvec(3,3) = 19._qp
      bvec(4,3) = one

    end subroutine fix_bvec

    
end module qpglobal
