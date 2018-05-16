module define_ampl
  use types; use sub_defs_io 
  implicit none
  private

  integer, public :: Npoint      ! N-point scattering amplitude
  integer, public :: Nterm       ! where fermion line terminates
  integer, public :: Ncut        ! # of ``highest level masters''
  integer, public :: Ncut2 
  integer, public :: NcutZ
  integer, public :: NtermZ
  integer, public :: NpointZ      

  ! -- things realted to amplitudes with additional quarks 

  logical, public :: qbq_and_gluons         
  logical, public :: ferm_loops           
  logical, public :: extra_ferm_pair      
  logical, public :: extra_ferm_pair1     
  logical, public :: extra_ferm_pair2     
  logical, public :: extra_ferm_pair_nf   
  logical, public :: gluonsonly           
  logical, public :: gluons_ferm_loops           
  logical, public :: qbq_and_gluonsonly   
  integer, public :: ampl_type 
  logical, public :: ferm_loops_Z   
  logical, public :: ferm_loops_Z_sbs   
  logical, public :: qbq_WW_and_gluons         
  logical, public:: case_a2 = .false.
  logical, public:: case_a3 = .false.
  logical, public:: case_a4 = .false.
  logical, public:: case_a1 = .false.
  logical, public:: case_b1 = .false.
  logical, public:: case_b2 = .false.
  logical, public:: case_one
  logical, public:: WWqqqq
  logical, public :: ferm_loops_WW_nf 

  logical, public :: include_gamZ ! In WW decide whether to include gam/Z 

  ! variable deciding whether to swap the fermions in the 2nd fermion pair 
  ! Q and QB -- for extra_ferm_pair swap can be tru or false
  ! for extra_ferm_pair_1,extra_ferm_pair_2, extra_ferm_pair_nf swap 
  ! should be true 
  ! for qbq_and_gluons and ferm_loops swap .eq. true does not make sense (and 
  ! helicities are not properly swapped in the polarization routine) 
  logical, public :: swap 

  !-- the next line is for the amplitudes with additional QUARK pair 
  !   indicating where the this new (bar) quark enters the diagram (NQ) and 
  !   leaves (NBQ) the diagram 

  integer, public :: NQU 
  integer, public :: NQUB 

  integer, public :: Nglue 

  ! other things specifying what to do with the numbers... 
  logical,public :: doplot,writeout, rn_point   
  integer,public :: iterm_min, iterm_max, save_after, nev   

  ! decide whether to use pol_mass or pol_dk2mom for the W polarization (obsolete) 
  logical, public :: use_pol_mass 
  ! perform gauge check 
  logical, public :: gauge_check 


  ! -- amplitudes with W 
  integer, parameter, public :: iqbq_and_gluons     = 1
  integer, parameter, public :: iferm_loops         = 2
  integer, parameter, public :: iextra_ferm_pair    = 3
  integer, parameter, public :: iextra_ferm_pair1   = 4 
  integer, parameter, public :: iextra_ferm_pair2   = 5
  integer, parameter, public :: iextra_ferm_pair_nf = 6 
  ! -- amplitudes without W 
  integer, parameter, public :: igluonsonly         = 7 
  integer, parameter, public :: iqbq_and_gluonsonly = 8 
  integer, parameter, public :: igluons_ferm_loops  = 10 
  ! -- amplitudes with Z 
  integer, parameter, public :: iferm_loops_Z       = 9
  integer, parameter, public :: iferm_loops_Z_sbs   = 11
  integer, parameter, public :: iqbq_WW_and_gluons  = 12


  integer, parameter, public :: v_vector   = 1
  integer, parameter, public :: v_axial  = 2
  integer, public            :: v_coupling 

  integer, parameter, public :: max_warn = 10 
  integer, public            :: i_warn = 0  
  logical, public            :: MChel, MCchn, leadingcol, subleadingcol   
  logical, public            :: cashing 

  public :: initialize_params

contains 


  ! initialize here all commond-line parameters 
  subroutine initialize_params(np) 
    integer, intent(in) :: np 
    !-------input parameters 
    rn_point     = log_val_opt('-rn',.false.) 
    npoint       = int_val_opt('-npoint',np)
    iterm_min    = int_val_opt('-iterm_min',3)
    iterm_max    = int_val_opt('-iterm_max',npoint)
    iterm_min    = int_val_opt('-nterm',iterm_min)
    iterm_max    = int_val_opt('-nterm',iterm_max)
    nev          = int_val_opt('-nev',1)
    doplot       = log_val_opt('-doplot',.true.) 
    save_after   = int_val_opt('-saveafter',nev)
    writeout     = log_val_opt('-writeout',.false.)
    use_pol_mass = log_val_opt('-use_pol_mass',.false.)
    gauge_check  = log_val_opt('-gauge_check',.false.)

    qbq_and_gluons     = log_val_opt('-qbq_and_gluons',.false.)
    ferm_loops         = log_val_opt('-ferm_loops',.false.)
    ferm_loops_Z       = log_val_opt('-ferm_loops_Z',.false.)
    ferm_loops_Z_sbs   = log_val_opt('-ferm_loops_Z_sbs',.false.)
    extra_ferm_pair    = log_val_opt('-extra_ferm_pair',.false.)
    extra_ferm_pair1   = log_val_opt('-extra_ferm_pair1',.false.)
    extra_ferm_pair2   = log_val_opt('-extra_ferm_pair2',.false.)
    extra_ferm_pair_nf = log_val_opt('-extra_ferm_pair_nf',.false.)
    gluonsonly         = log_val_opt('-gluonsonly',.false.)
    gluons_ferm_loops  = log_val_opt('-gluons_ferm_loops',.false.)
    qbq_and_gluonsonly = log_val_opt('-qbq_and_gluonsonly',.false.)
    qbq_WW_and_gluons  = log_val_opt('-qbq_WW_and_gluons',.false.)
    WWqqqq             = log_val_opt('-WWqqqq',.false.)
    case_a1            = log_val_opt('-case_a1',.false.)
    case_a2            = log_val_opt('-case_a2',.false.)
    case_a3            = log_val_opt('-case_a3',.false.)
    case_a4            = log_val_opt('-case_a4',.false.)
    case_b1            = log_val_opt('-case_b1',.false.)
    case_b2            = log_val_opt('-case_b2',.false.)
    ferm_loops_WW_nf   = log_val_opt('-ferm_loops_WW_nf',.true.)
    cashing            = log_val_opt('-cashing',.false.)

    swap = log_val_opt('-swap',.false.)
    if ((extra_ferm_pair1 .or. extra_ferm_pair2 .or. extra_ferm_pair_nf) .and. &
         &.not.swap) write(*,*) 'WARNING NO QQB SWAP PERFORMED!'

    if ((qbq_and_gluons .or. ferm_loops .or. ferm_loops_Z) .and. &
         &swap) write(*,*) 'WARNING SWAP PERFORMED!'

    NQU   = int_val_opt('-NQU',4)
    NQUB  = int_val_opt('-NQUB',5)
    Nglue = int_val_opt('-Nglue',0)

    ! when nothing is specified do q bq W + n gluons as default 
    if (.not.ferm_loops .and. .not. extra_ferm_pair .and. .not. extra_ferm_pair1&
         &.and. .not. extra_ferm_pair2 .and. .not. extra_ferm_pair_nf &
         &.and. .not. gluonsonly .and. .not. qbq_and_gluonsonly .and. &
         &.not. ferm_loops_Z .and. .not. gluons_ferm_loops .and. .not. ferm_loops_Z_sbs) &
         &qbq_and_gluons = .true. 


    if (qbq_and_gluons    ) ampl_type = iqbq_and_gluons    
    if (ferm_loops        ) ampl_type = iferm_loops        
    if (ferm_loops_Z      ) ampl_type = iferm_loops_Z        
    if (ferm_loops_Z_sbs  ) ampl_type = iferm_loops_Z_sbs        
    if (extra_ferm_pair   ) ampl_type = iextra_ferm_pair   
    if (extra_ferm_pair1  ) ampl_type = iextra_ferm_pair1  
    if (extra_ferm_pair2  ) ampl_type = iextra_ferm_pair2  
    if (extra_ferm_pair_nf) ampl_type = iextra_ferm_pair_nf
    if (gluonsonly        ) ampl_type = igluonsonly
    if (gluons_ferm_loops ) ampl_type = igluons_ferm_loops
    if (qbq_and_gluonsonly) ampl_type = iqbq_and_gluonsonly    
    if (qbq_WW_and_gluons ) ampl_type = iqbq_WW_and_gluons


    if (qbq_and_gluons    ) write(*,*) ' ---> Doing qbq_and_gluons    '
    if (ferm_loops        ) write(*,*) ' ---> Doing ferm_loops        '
    if (ferm_loops_Z      ) write(*,*) ' ---> Doing ferm_loops_Z      ' 
    if (ferm_loops_Z_sbs  ) write(*,*) ' ---> Doing ferm_loops_Z_sbs  ' 
    if (extra_ferm_pair   ) write(*,*) ' ---> Doing extra_ferm_pair   '
    if (extra_ferm_pair1  ) write(*,*) ' ---> Doing extra_ferm_pair1  '
    if (extra_ferm_pair2  ) write(*,*) ' ---> Doing extra_ferm_pair2  '
    if (extra_ferm_pair_nf) write(*,*) ' ---> Doing extra_ferm_pair_nf'
    if (gluonsonly        ) write(*,*) ' ---> Doing gluonsonly        '
    if (gluons_ferm_loops ) write(*,*) ' ---> Doing gluons_ferm_loops '
    if (qbq_and_gluonsonly) write(*,*) ' ---> Doing qbq_and_gluonsonly'   
    if (qbq_WW_and_gluons ) write(*,*) ' ---> Doing qbq_WW_and_gluons '   

    ! -- decide which type of W/Z amplitude 
    v_coupling= int_val_opt('-v_coupling',v_vector)
    if (v_coupling == v_vector) then 
       write(*,*) 'Doing V vector coupling'
    elseif (v_coupling == v_axial) then 
       write(*,*) 'Doing V axial coupling'
    else
       stop 'initialize_params: undefined vector coupling' 
    endif

    include_gamZ = log_val_opt('-include_gamZ', .true.) 
    if ((case_a1) .or. (case_b1)) case_one = .true.

    MChel = log_val_opt('-MChel', .false.)     
    MCchn = log_val_opt('-MCchn', .false.) 
    leadingcol    = log_val_opt('-leadingcol', .false.)     
    subleadingcol = log_val_opt('-subleadingcol', .false.) 

  end subroutine initialize_params


end module define_ampl



