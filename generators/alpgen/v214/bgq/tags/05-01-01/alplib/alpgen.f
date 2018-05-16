c     Version to allow reading in a standard event, and perform the kt
C     clustering and alphas rescaling
C-----------------------------------------------------
      program alpgen
C-------------------------------------------------------------------------
C 
C driver for multi-parton matrix element generator based on ALPHA
C
C--------------------------------------------------------------------------
C
#ifdef USE_MPI
c     initialize mpi
      call initialize_mpi
c
#endif
      call alsprc
c
c     Routine location: XXX.f
c     Purpose: 
c        -  assign the hard process code for selected process XXX
c
C setup running defaults:
c
      call alsdef
c
c     Routine location: alpgen.f
c     Purpose: 
c        -  initialise the default generation parameters (e.g. beam
c           type, energy, PDF set) 
c        -  initialise default mass and couplings for particles
c
c
C setup  event generation options, bookkeeping, etc
      call alhset
c
c     Routine location: XXX.f
c     Purpose: 
c        - setup event generation options, bookkeeping, etc, for
c              the specific hard process XXX
c        - write run information on stat, unwpar files

c initialise alpha parameters
      call alinit
c     Routine location: alpgen.f
c     Purpose: 
c        - evaluate parameter-dependent quantities (e.g. Higgs width)
c	 - fill the apar array of ALPHA (which contains all parameters
c	   required by ALPHA
c
c
C setup internal bookkeeping histograms
      call alsbkk
c
c     Routine location: alpgen.f
c     Purpose: initialise histograms common to all processes, and needed
C     for internal purposes
c
c
C setup user histograms
#ifndef USE_MPI
      call alshis
#endif
c
c     Routine location: XXXusr.f
c     Purpose: initialise histograms
c
c
c
C setup integration grids, including optimization if required
c
      call alsgrd
c
c     Routine location: XXX.f
c     Purpose: setup integration grid variables
      call aligrd
c
c     Routine location: alpgen.f
c     Purpose: initialise grid with warm-up iterations, if required
c
c
C generate events
c
      Call alegen
c
c     Routine location: alpgen.f
c     Purpose: generates events, calling in a standad format the
c          the process-depepdent phase-space and flavour-selection
c          routines, contained in XXX.f

C finalise histograms
c
#ifndef USE_MPI
      call alfhis
#endif
c
c     Routine location: XXXusr.f
c     Purpose: finalize analysis and histograms
      call alfbkk
c
c     Routine location: alpgen.f
c     Purpose: finalize internal histograms
c
c
#ifdef USE_MPI
c     FINALIZE MPI
      call finalize_mpi
c
c
#endif
      end

c-------------------------------------------------------------------
      subroutine alsbkk
c     setup weight bookeeping histograms
c-------------------------------------------------------------------
      implicit none
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
      integer i
      wmin=log10(1d-20)
      wmax=log10(1d20)
      wbin=(wmax-wmin)/1000d0
      do i=1,1000
        logwgt(i)=wmin+wbin*(i-0.5)
        wgtdis(i)=0d0
      enddo
c under- and overflow bin
      wgtdis(0)=0d0
      wgtdis(1001)=0d0
      call mbook(190,'rewgt factors',0.02,0.,2.)
      end

c-------------------------------------------------------------------
      subroutine alhbkk(wgt)
c     online bookeeping of weight distributions
c-------------------------------------------------------------------
      implicit none
      double precision wgt
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
      integer i
      double precision lwgt
      if(wgt.le.0.d0) return
      lwgt=log10(wgt)
      i=int((lwgt-wmin)/wbin)+1
      if(i.lt.0) then
c underflow bin
        wgtdis(0)=wgtdis(0)+wgt
      elseif(i.le.1000) then
        wgtdis(i)=wgtdis(i)+wgt
      else
c overflow bin
        wgtdis(1001)=wgtdis(1001)+wgt
      endif
      end

c-------------------------------------------------------------------
      subroutine alfbkk
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
c locals
      double precision sigwgt(0:1001),newwmax(5)
      integer nev,iproc
      real Sq,Savgwgt
      integer i,j
      real  xnorm
      character*80 tmpstr
      integer aluisl,nfile1,nfile2

c find maxwgt thresholds for 5-4-3-2-1% of total weight
      if(imode.eq.1) then
        sigwgt(1001)=wgtdis(1001)
        do i=1000,0,-1
          sigwgt(i)=sigwgt(i+1)+wgtdis(i)
        enddo
        do i=1,1000
          do j=1,5
            if(sigwgt(i)/sigwgt(0).gt.0.01*float(j)) then
              newwmax(j)=logwgt(i)
            endif
          enddo
        enddo
        write(niopar,'(5(e12.6,1x))')(10d0**(newwmax(i)+0.5*wbin),i=1,5)
        close(niopar)
      endif

c     
      if(imode.eq.2) then
        close(niounw)
      endif
c      open(unit=99,file=topfile,err=999,status='old')
c      call aluend(99)
c compatbility with latest gfortran
#ifndef USE_MPI
      open(unit=99,file=topfile,position='append',err=999,status='old')
#endif
      if(imode.le.1) then
         xnorm=sngl(avgwgt/totwgt)
      elseif(imode.eq.2) then
         xnorm=1e0/real(unwev)
      else
         write(6,*) 'imode type not allowed, stop'
         stop
      endif
c
      do i=160,200
         if(i.ne.192) call mopera(i,'F',i,i,xnorm,1.)
         call mfinal(i)
      enddo 
c
c
#ifndef USE_MPI
      call mtop(190,99,'rewgt ','dN/d rewgt ','LIN')
      close(99)
#endif

 999  return
      end

c-------------------------------------------------------------------
      subroutine alsdef
c     setup default run parameters
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
c common declarations
c   for MPI but always included as these get initialized for the MT
      integer mpirank,mpiworldsize
      character*10 time
      character*8 date,rankstring
      common/mpi/mpirank,mpiworldsize,rankstring


c local variables
      double precision dummy,ranset,ranset2
      integer i,retval,chdir
      character*100 input_filename,worldsize
      integer iargs,input_file_status
c
c beginning of MPI block, get MPI rank and set common variables
#ifdef USE_MPI
      character*256 output_path
      include 'mpif.h'
#else
      mpirank = 0
      write( rankstring,*) ''
#endif
c end of MPI block

c print rank and time
      call date_and_time(date,time)
      write(worldsize,*) mpiworldsize
      write(6,*) ' Rank: ',rankstring,' of ',worldsize,
     $           ' Date/Time: ',date,' ',time

c list processes
      do i=1,nprocs
        chprc(i)='unavailable'
      enddo
      chprc(1)='wqq'
      chprc(2)='zqq'
      chprc(3)='wjet'
      chprc(4)='zjet'
      chprc(5)='vbjet'
      chprc(6)='2Q'
      chprc(7)='4Q'
      chprc(8)='QQh'
      chprc(9)='Njet'
      chprc(10)='wcjet'
      chprc(11)='phjet'
      chprc(12)='hjet'
      chprc(13)='top'
      chprc(14)='wphjet'
      chprc(15)='wphqq'
      chprc(16)='2Qph'
      nactprc=16
c fixed parameters:
      resc=1
c     number of calls for the spin/colour event-by-event average
      navg=1
c     beam 1   (ih=1: proton;  ih=-1: pbar)
      ih1= 1
c     lepton masses 
      mlep(1)=0.511d-3
      mlep(2)=0.10566d0
      mlep(3)=1.777d0
c     ckm mixings
      scab2=0.05
      ccab2=1-scab2
c     EW parameters:
      mw=80.419d0
      mz=91.188d0
      stw=sqrt(0.231d0)
      aem=1.d0/128.89
      gfermi=1.16639d-5
C--   To ensure the gauge invariance of the ALPHA calculations,
c     we shall use the LO relations between mW, mZ, weak couplings,
c     GFermi and alpha(em).  Which parameters can be input, and which
c     are calculated, is governed by the option IEWOPT, which can be
c     reassigned in the ALSUSR routine. The implementation of the option
c     is contained in the routine ALINIT
c
      iewopt=3
c eliminate propagating photons in (w)phjet proceses:
      if(ihrd.eq.2.or.ihrd.eq.4.or.ihrd.eq.11.
     +   or.ihrd.eq.14.or.ihrd.eq.15) resc=1d3
c scale setting parameters:
      iqopt=1
      qfac=1d0
      ktfac=1d0
      ickkw=0
c clustering parameters
c     1: pclu(1:4)=p1(1:4)+p2(1:4)  2: pclu(1:3)=p1(1:3)+p2(1:3), mclu=0
      mrgopt=1
      cluopt=1
c
c     assign Default PARameter values for usr-accessible parameters
      call aldpar(1)
c
C     MANDATORY INPUTS:
#ifdef USE_MPI
c input from file
c     rank 0 checks the command line arguments
c     the output path is appended to the weighted event data file,
c        the weighted event parameters file, and the 
c        unweighted event output files.
      if(mpirank == 0) then
         iargs = COMMAND_ARGUMENT_COUNT()
c        if only one argument passed it should be the input config file
         if(iargs.eq.1) then
           call getarg(1,input_filename)
c          set the output path to be local
           output_path = '.'
c        if two arguments passed then there is also an output path
         elseif(iargs.eq.2) then
           call getarg(1,input_filename)
           call getarg(2,output_path)
           write(6,*) 'Using output path:',output_path
c        print usage and exit if any other number of arguments passed
         else
           write(6,*) 'No input file provided'
           write(6,*) 'usage: ./executable input-card',
     $                ' [optional-output-path]'
           call exit(-1)
         endif

         write(6,*) 'Reading input card: ',input_filename
         open(9,file=input_filename,action='read',
     $          iostat=input_file_status)
         if(input_file_status /= 0) then
            write(6,*) ' Error opening input file: ',input_filename
            call exit(-1)
         endif
         read(9,*,iostat=input_file_status) imode
         if(input_file_status /= 0) then
            write(6,*) ' Error reading from file: ',input_filename
            write(6,*) '                  status: ',input_file_status
            call exit(-1)
         endif
      endif
      call MPI_BCAST(imode,1,MPI_INT,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting imode '
         call exit(-1)
      endif

      call MPI_BCAST(output_path,256,MPI_CHAR,0,
     $                    MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting output_path '
         call exit(-1)
      endif

#else
      write(6,*) 'Input RUN generation mode:'
      write(6,*) '0: generate weighted events, no evt dumps to file'
      write(6,*)
     $   '1: generate wgtd events, write to file for later unweighting'
      write(6,*)
     $   '2: read events from file for unweighting or processing'
      write(6,*)
     $ 'or documentation modes:'
      write(6,*)
     $ '3: print parameter options and defaults, then stop'
      write(6,*)
     $ '4: write to par.list parameter options and defaults, then stop'
      write(6,*)
     $ '5: write to prc.list complete list of processes, parameter ',
     $ '   options and defaults, scale choices, PDF, etc., then stop'
      read(5,*) imode
#endif
      
      write(6,*) 'imode = ',imode
#ifdef USE_MPI
      if(imode.eq.0.and.mpiworldsize.gt.1) then
         write(6,*) 'ERROR Running more than 1 MPI Rank in imode = 0'
         write(6,*) 'imode = 0 cannot be parallelized and should be'
         write(6,*) '  run with only one rank '
         write(6,*) ' rank: ',mpirank,' world size: ',mpiworldsize
         call exit(-1)
      endif
#endif
c
c documentation modes:
      if(imode.ge.3.and.imode.le.5) then
        call alppar(imode)
        stop
      endif

#ifdef USE_MPI
      if(mpirank == 0) then
         read(9,*,iostat=input_file_status) fname
         if(input_file_status /= 0) then
            write(6,*) 'error reading file base'
            call exit(-1)
         endif
      endif
      call MPI_BCAST(fname,140,MPI_CHAR,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting fname '
         call exit(-1)
      endif
#else
      write(6,*) 'Input string labeling output and input files'
      write(6,*) '(e.g. w2j to output files w2j.stat, etc.)'
      read(5,*) fname
#endif
      write(6,*) 'read fname: ',fname
c for MPI add rank number to fname
#ifdef USE_MPI
c only use the rank string for modes 1 & 2 in MPI mode
c otherwise use just the basename
      if (imode.ne.0) then
c        add rank number to the output file name base for uniqueness
         call alustc(fname,'.'//rankstring,fname_mpirank)
c        add the output path to the output file name
         call alustc(output_path,'/'//fname_mpirank,fname_mpirank)
      else
         fname_mpirank = fname
      endif
#else
      fname_mpirank = fname
#endif
      write(6,*) 'basename: ',fname
      write(6,*) 'basename for MPI: ',fname_mpirank
c     define input/output files
      call alstio
c     
      if(imode.eq.2) then
c read in parameters from fname.par
        call alrpar
      elseif(imode.lt.2) then
c read in parameters from input
        call alrusr
      else
        write(6,*) 'With imode=',imode,' we should not be here, stop'
        stop
      endif
c     save parameter values to relative variables in common blocks
      call alspar
c     deposit random number generator seeds
C     ADD MPI Rank to ISEED(2), has no affect if MPI not used
      dummy= ranset(iseed)
      if(imode.eq.2) dummy= ranset2(iseed2)
      end
      

c-------------------------------------------------------------------
      subroutine alrusr
c     input datacards superseding default parameters
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      double precision chvalue,tmp
      character*8 chparam
      integer iflag,retval
      integer aluisl,itmp,ntmp

#ifdef USE_MPI
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      include 'mpif.h'

      if(mpirank == 0) then
         read(9,*) igrid
      endif
      call MPI_BCAST(igrid,1,MPI_INT,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting igrid '
         call exit(-1)
      endif
#else
c     
c     select options for grid selection/generation
      write(6,*) ' '
      write(6,*) 'To generate new grid input 0, '
      write(6,*) 'To use grid generated by the warmup iterations of',
     +     ' the previous run, input 1,'
      write(6,*) 'To use grid generated by the event generation of',
     +     ' the previous run, input 2:'
      read(5,*) igrid
#endif

#ifdef USE_MPI
      if(mpirank == 0) then
         read(9,*) nopt,niter
      endif
      call MPI_BCAST(nopt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting nopt '
         call exit(-1)
      endif
      call MPI_BCAST(niter,1,MPI_INT,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting niter '
         call exit(-1)
      endif
#else
      write(6,*) 'Input N(events)/iteration and N(iterations) for the 
     + warmup iterations:'
      read(5,*) nopt,niter
#endif

#ifdef USE_MPI
      if(mpirank == 0) then
         read(9,*) maxev
      endif
      call MPI_BCAST(maxev,1,MPI_REAL8,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) ' Error broadcasting maxev '
         call exit(-1)
      endif
#else
      write(6,*) 'Input number evts to generate:'
      read(5,*) maxev
#endif
c     input scale options
      write(6,*) ' grid source = ',igrid
      write(6,*) ' N(events)/iterations and N(iterations) = ',nopt,
     $ ',',niter
      write(6,*) ' number of events to generate = ',maxev
      call alhsca(ihrd,6)
c     
c     
c     
#ifdef USE_MPI
 1    if(mpirank == 0) write(6,*) ' '
 2    if(mpirank == 0) then
         read(9,*,end=3,err=10) chparam,chvalue
      endif
      call MPI_BCAST(chparam,8,MPI_CHAR,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) mpirank,' Error broadcasting chparam '
         call exit(-1)
      endif
      call MPI_BCAST(chvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               retval)
      if(retval /= 0) then
         write(6,*) mpirank,' Error broadcasting chvalue '
         call exit(-1)
      endif
#else
 1    write(6,*) ' '
      write(6,*) 'Input parameters to replace defaults:  type and value'
      write(6,*)
     $    '(input ''print 1'' to display the list of parameter types and
     $ their current values)'
      write(6,*)
     $     '(input ''print 2'' to write the list to file par.list)'
      write(6,*)
     $     '(input ''ctrl-D'' to terminate the input sequence)'
 2    read(5,*,end=3,err=10) chparam,chvalue
#endif
      
 10   if(chparam(1:3).eq.'eoi') goto 3
      if(chparam(1:1).eq.'*') goto 2
      itmp=aluisl(chparam)
      call alfpar(chparam(1:itmp),chvalue,iflag)
      if(iflag.le.1) then
        goto 2
      elseif(iflag.eq.2) then
        goto 1
      elseif(iflag.eq.3) then
#ifndef USE_MPI
        write(6,*) 'param ',chparam(1:itmp),
     +   ' not available/changeable for this process or imode=',imode
#endif
        goto 1
      endif
#ifdef USE_MPI
 3    if(mpirank.eq.0) then
         call MPI_BCAST('eoi',8,MPI_CHAR,0,MPI_COMM_WORLD,retval)
          if(retval /= 0) then
            write(6,*) mpirank,' Error broadcasting eoi '
            call exit(-1)
         endif
         call MPI_BCAST(1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $                  retval)
         if(retval /= 0) then
            write(6,*) mpirank,' Error broadcasting 1 '
            call exit(-1)
         endif
      endif
      write(niopar,*) 'eoi',1
#else
 3    write(niopar,*) 'eoi',1
#endif
      return
      end

c-------------------------------------------------------------------
      subroutine alhsca(ihrd,iunit)
c     defines the hard process scales
c-------------------------------------------------------------------
      implicit none
      integer ihrd,iunit
c     qsq scale; 
      write(iunit,*)
     $     'Options for Factorization/renormalization scale Q:'
      write(iunit,*) 'iqopt=0 => Q=qfac'
      if(ihrd.eq.1) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.2) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m0^2+ sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*m0'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m0^2 + pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
        write(iunit,*) '- m0=mll' 
      elseif(ihrd.eq.3) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
      elseif(ihrd.eq.4) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m0^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*m0'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m0^2 + pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where m0=mll' 
      elseif(ihrd.eq.5) then
        write(iunit,*) 
     $ 'iqopt=1 => Q=qfac*sqrt{sum(mV)^2+sum(pt_photons^2+pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sum(mV)'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{shat}'
      elseif(ihrd.eq.6) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.7) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.8) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{mh^2+sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.9) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*s)'
      elseif(ihrd.eq.10) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
      elseif(ihrd.eq.11) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(pt_ph^2+pt_jets^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{pt_jets^2}'
      elseif(ihrd.eq.12) then 
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{nh*mh^2+pt_jets^2}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{shat}'
      elseif(ihrd.eq.13) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt(mw^2+sum(m_tr^2))'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{shat}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
        write(iunit,*) 
     $ '- the mw^2 term is present only for processes with a W'
      elseif(ihrd.eq.14) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons and jets'
      elseif(ihrd.eq.15) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons, heavy quarks and light jets'
      elseif(ihrd.eq.16) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      endif
c
c      write(iunit,*) ' '
c      write(iunit,*) 'Default is iqopt=1, qfac=1'
c      write(iunit,*) ' '
      if(ihrd.le.6.or.ihrd.eq.9.or.ihrd.eq.10.or.ihrd.eq.11.or.ihrd.eq
     $     .12.or.ihrd.eq.14.or.ihrd.eq.15.or.ihrd.eq.16) then
        write(iunit,*) 'To select CKKW scale input ''ickkw 1'' '
        write(iunit,*) '(mandatory for later use of jet matching)'
        write(iunit,*) 'In imode=2 select clustering option ''cluopt'':'
        write(iunit,*) 'cluopt=1: kperp propto pt(cluster) (default)'
        write(iunit,*) 'cluopt=2: kperp propto mt(cluster)'
        write(iunit,*) 'kperp is then rescaled by ktfac'
      endif
      end

c-------------------------------------------------------------------
      subroutine alrpar
c     read in parameters different from defaults for imode=2
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
#ifdef USE_MPI
      integer mpirank,mpiworldsize,retval
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      include 'mpif.h'
      integer filestat
#endif
      double precision chvalue
      character chparam*8
      integer aluisl,iflag,itmp
C first read in parameters from the imode=1 run
      itmp=0
 1    read(niopar,*,end=2,err=10) chparam,chvalue
 10   if(chparam(1:3).eq.'eoi') goto 2
      if(chparam(1:3).eq.'***') goto 1
#ifndef USE_MPI
      if(itmp.eq.0) then
        itmp=1
        write(6,*) ' '
        write(6,*) 'read in generation parameters:'
      endif
#endif
      call alfpar(chparam,chvalue,iflag)
      if(iflag.le.1) then
        goto 1
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam,' not recognised, stop'
        stop
      else
        write(6,*) 'unrecognized status after parameter input, stop'
        stop
      endif
 2    continue
c Restore defaults of imode=2 params that depend on imode=1 inputs:
c No option for Z decays if Z->nunu
      if((ihrd.eq.2.or.ihrd.eq.4).and.(ilep.eq.1)) then
        paruse(152,ihrd)=0
      endif
c Then read in decay parameters specific to the imode=2 run:
#ifdef USE_MPI
 15   continue
 20   if(mpirank.eq.0) then
         read(9,*,end=200,err=100) chparam,chvalue
      endif
      call MPI_BCAST(chparam,8,MPI_CHAR,0,MPI_COMM_WORLD,retval)
      if(retval /= 0) then
         write(6,*) '> Error broadcasting chparam '
         call exit(-1)
      endif
      call MPI_BCAST(chvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               retval)
      if(retval /= 0) then
         write(6,*) '> Error broadcasting chvalue '
         call exit(-1)
      endif
#else
 15   write(6,*) ' '
      write(6,*)'Input decay params to replace defaults: type and value'
      write(6,*)
     $   '(input ''print 1'' to display the list of parameter types and
     $ their current values)'
      write(6,*)
     $     '(input ''ctrl-D'' to terminate the input sequence)'
 20   read(5,*,end=200,err=100) chparam,chvalue
#endif
 100  if(chparam(1:3).eq.'eoi') goto 200
      if(chparam(1:1).eq.'*') goto 20
      itmp=aluisl(chparam)
      call alfdpa(chparam(1:itmp),chvalue,iflag)
      if(iflag.le.1) then
        goto 20
      elseif(iflag.eq.2) then
        goto 15
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam(1:itmp),' not accepted in imode=2'
        goto 15
      endif
#ifdef USE_MPI
 200  if(mpirank.eq.0) then
         call MPI_BCAST('eoi',8,MPI_CHAR,0,MPI_COMM_WORLD,retval)
          if(retval /= 0) then
            write(6,*) mpirank,' Error broadcasting eoi '
            call exit(-1)
         endif
         call MPI_BCAST(1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $                  retval)
         if(retval /= 0) then
            write(6,*) mpirank,' Error broadcasting 1 '
            call exit(-1)
         endif
      endif
      return
#else
 200  return
#endif
      end

      
c-------------------------------------------------------------------
      subroutine Alfpar(chparam,chvalue,iflag)
c     deposit parameter values into common blocks and store in fname.par
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
#ifdef USE_MPI
      integer mpirank
      character*8 rankstring
      common/mpi/mpirank,rankstring
#endif
      double precision chvalue
      character chparam*8,chtmp*8
      integer i,iflag,ipar,itmp,aluisl,istall
      data istall/0/
      iflag=0
      itmp=aluisl(chparam)
      if(chparam(1:5).eq.'print') then
        call alppar(int(chvalue))
        iflag=2
        return
      elseif(chparam(1:3).eq.'***') then
        iflag=1
        return
      else
        do i=1,npar
          chtmp=chpar(i)
          if(chparam(1:itmp).eq.chtmp(1:itmp)) then
            ipar=i
            if(paruse(ipar,ihrd).eq.0) then
              iflag=3
              return
            endif
            goto 100
          endif
        enddo
      endif
      iflag=3
      return
 100  parval(ipar)=chvalue
#ifdef USE_MPI
c     increment iseed(2) and iseed2(2) by the mpirank
c     for imode 0 and 1 only
c     don't increment for imode 2 because then it is incremented
c     twice
      if((ipar.eq.91.or.ipar.eq.191).and.imode.lt.2) then
         parval(ipar) = parval(ipar) + mpirank
         chvalue = parval(ipar)
      endif
#endif
      if(ipar.eq.4.and.chvalue.eq.-1) then
        if(istall.gt.1000) then
          write(6,*) 'PDF code unspecified, stop'
          stop
        endif
        call prntsf(6)
        write(6,*) ' enter ''ndns'' followed by value'
        iflag=1
        istall=istall+1
        return
      endif
      if(paruse(ipar,ihrd).eq.1) then
        if(imode.eq.1) write(niopar,*) chpar(ipar),chvalue
#ifndef USE_MPI
        if(imode.eq.2) write(6,*) chpar(ipar),'=',chvalue
#endif
      endif
      end

c-------------------------------------------------------------------
      subroutine Alfdpa(chparam,chvalue,iflag)
c     deposit parameter values for imode=2 into common blocks and store
C     in fname.par 
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      double precision chvalue
      character chparam*8,chtmp*8
      integer i,iflag,ipar,itmp,aluisl,istall
      data istall/0/
      iflag=0
      itmp=aluisl(chparam)
      if(chparam(1:5).eq.'print') then
        call alppar(6)
        iflag=2
        return
      elseif(chparam(1:3).eq.'***') then
        iflag=1
        return
      else
        do i=1,npar
          chtmp=chpar(i)
          if(chparam(1:itmp).eq.chtmp(1:itmp)) then
            ipar=i
            if(paruse(ipar,ihrd).eq.0) then
              iflag=3
              return
            endif
            goto 100
          endif
        enddo
      endif
      iflag=3
      return
 100  if(ipar.le.150) then
        iflag=3
        return
      endif 
      parval(ipar)=chvalue
      if(paruse(ipar,ihrd).eq.1) write(6,*) chpar(ipar),'=',chvalue
      end
c-------------------------------------------------------------------
      subroutine Aldpar(n)
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      double precision chvalue
      character chparam*8
      integer n,iunit,i,j,aluisl,itmp
      do i=1,npar
        chpar(i)='***'
        chpdes(i)='parameter not assigned'
        do itmp=1,nprocs
          paruse(i,itmp)=0
        enddo
      enddo
c     beam parameters
c
      chpar(2)='ih2'
      chpdes(2)='Select pp (1) or ppbar (-1) collisions'
      partyp(2)=1
      parval(2)=-1
      do i=1,nactprc
        paruse(2,i)=1
      enddo
c
      chpar(3)='ebeam'
      chpdes(3)='beam energy in CM frame (e.g. 7000 for LHC)'
      partyp(3)=0
      parval(3)=980
      do i=1,nactprc
        paruse(3,i)=1
      enddo
c
      chpar(4)='ndns'
      chpdes(4)='parton density set'
c currently available:
c ndns= 1      2      3       4      5      6       7      8
c pdf = cteq4m cteq4l cteq4hj cteq5m cteq5l cteq5hj cteq6m cteq6l
c ndns= 101    102        103        104        105       
c pdf = mrst99 mrst2002-1 mrst2002-2 mrst2002-3 mrst2002-4
      partyp(4)=1
      parval(4)=5
      do i=1,nactprc
        paruse(4,i)=1
      enddo
c
      chpar(5)='iqopt'
      chpdes(5)='scale option (process dependent)'
      partyp(5)=1
      parval(5)=1
      do i=1,nactprc
        paruse(5,i)=1
      enddo
c
      chpar(6)='qfac'
      chpdes(6)='Q scale rescaling factor'
      partyp(6)=0
      parval(6)=1d0
      do i=1,nactprc
        paruse(6,i)=1
      enddo
c
      chpar(7)='ickkw'
      chpdes(7)
     $     ='CKKW scale option: set to 1 to enable jet-parton matching'
      partyp(7)=1
      parval(7)=0
      do i=1,6
        paruse(7,i)=1
      enddo
      do i=9,12
        paruse(7,i)=1
      enddo
      do i=14,16
        paruse(7,i)=1
      enddo
c
      chpar(8)='ktfac'
      chpdes(8) ='scale factor for ckkw alphas scale'
      partyp(8)=0
      parval(8)=1d0
      do i=1,6
        paruse(8,i)=1
      enddo
      do i=9,12
        paruse(8,i)=1
      enddo
      do i=14,16
        paruse(8,i)=1
      enddo
c
      chpar(10)='njets'
      chpdes(10)='number of light jets'
      partyp(10)=1
      if(ihrd.eq.3.or.
     +   ihrd.eq.4.or.
     +   ihrd.eq.11.or.
     +   ihrd.eq.13) then
        parval(10)=1
      elseif(ihrd.eq.9) then
        parval(10)=2
      else
        parval(10)=0
      endif
      do i=1,nactprc
        paruse(10,i)=1
      enddo
c
      chpar(11)='ihvy'
      chpdes(11)='heavy flavour type for procs like WQQ, ZQQ, 2Q, etc'/
     $     /'(4=c, 5=b, 6=t)'
      partyp(11)=1
      parval(11)=5
      paruse(11,1)=1
      paruse(11,2)=1
      paruse(11,6)=1
      paruse(11,7)=1
      paruse(11,8)=1
      paruse(11,15)=1
      paruse(11,16)=1
c
      chpar(12)='ihvy2'
      chpdes(12)='2nd heavy flavour type for procs like 4Q'
      partyp(12)=1
      parval(12)=5
      paruse(12,7)=1
c
      chpar(13)='nw'
      chpdes(13)='number of W bosons'
      partyp(13)=1
      parval(13)=0
      if(ihrd.eq.5) parval(13)=2      
      paruse(13,5)=1
c
      chpar(14)='nz'
      chpdes(14)='number of Z bosons'
      partyp(14)=1
      parval(14)=0      
      paruse(14,5)=1
c
      chpar(15)='nh'
      chpdes(15)='number of H bosons'
      partyp(15)=1
      parval(15)=0 
      if(ihrd.eq.12) parval(15)=1
      paruse(15,5)=1
c      paruse(15,12)=1
c
      chpar(16)='nph'
      chpdes(16)='number of photons'
      partyp(16)=1
      parval(16)=0
      if(ihrd.eq.11.or.ihrd.eq.14.or.ihrd.eq.15.or.ihrd.eq.16) parval(16
     $     )=1
      paruse(16,5)=1
      paruse(16,11)=1
      paruse(16,14)=1
      paruse(16,15)=1
      paruse(16,16)=1
c
      chpar(17)='zfstate'
      chpdes(17)='Final states for Z bosons in imode=1'
      partyp(17)=1
      parval(17)=0
      paruse(17,5)=1
c
c     masses
      chpar(20)='mc'
      chpdes(20)='charm mass'
c     the charm quark is considered massless unless explicitly requested
C     (e.g. in processes like W c cbar or c cbar)
c     it is 0 in W c processes
      partyp(20)=0
      parval(20)=0
      paruse(20,1)=1
      paruse(20,2)=1
      paruse(20,6)=1
      paruse(20,7)=1
      paruse(20,15)=1
      paruse(20,16)=1
      if(paruse(20,ihrd).eq.1) parval(20)=1.5
c
      chpar(21)='mb'
      chpdes(21)='bottom mass'
      partyp(21)=0
      parval(21)=4.7d0   
      paruse(21,1)=1
      paruse(21,2)=1
      paruse(21,6)=1
      paruse(21,7)=1
      paruse(21,8)=1
      paruse(21,13)=1
      paruse(21,15)=1
      paruse(21,16)=1
c
      chpar(22)='mt'
      chpdes(22)='top mass'
      partyp(22)=0
      parval(22)=174.3d0
      paruse(22,1)=1
      paruse(22,2)=1
      paruse(22,6)=1
      paruse(22,7)=1
      paruse(22,8)=1
      paruse(22,13)=1
      paruse(22,15)=1
      paruse(22,16)=1
c
      chpar(23)='mh'
      chpdes(23)='higgs mass'
      partyp(23)=0
      parval(23)=120d0
      paruse(23,5)=1
      paruse(23,8)=1
      paruse(23,12)=1
c
c     pt cuts
      chpar(30)='ptjmin'
      chpdes(30)='minimum pt for light jets'
      partyp(30)=0
      parval(30)=20d0
      do i=1,nactprc
        paruse(30,i)=1
      enddo
c
      chpar(31)='ptbmin'
      chpdes(31)='ptmin for bottom quarks (in procs with explicit b)'
      partyp(31)=0
      parval(31)=20d0
      paruse(31,1)=1
      paruse(31,2)=1
      paruse(31,6)=1
      paruse(31,7)=1
      paruse(31,8)=1
      paruse(31,13)=1
      paruse(31,15)=1
      paruse(31,16)=1
c
      chpar(32)='ptcmin'
      chpdes(32)='ptmin for charm quarks (in procs with explicit c)'
      partyp(32)=0
      parval(32)=20d0
      paruse(32,1)=1
      paruse(32,2)=1
      paruse(32,6)=1
      paruse(32,7)=1
      paruse(32,15)=1
      paruse(32,16)=1
      paruse(32,10)=1
c
      chpar(33)='ptlmin'
      chpdes(33)='minimum pt for charged leptons'
      partyp(33)=0
      parval(33)=0d0
      do i=1,4
        paruse(33,i)=1
      enddo
      paruse(33,10)=1
      paruse(33,14)=1
      paruse(33,15)=1
c
      chpar(34)='metmin'
      chpdes(34)='minimum missing et'
      partyp(34)=0
      parval(34)=0d0
      do i=1,4
        paruse(34,i)=1
      enddo
      paruse(34,10)=1
      paruse(34,14)=1
      paruse(34,15)=1
c
      chpar(35)='ptphmin'
      chpdes(35)='minimum pt for photons'
      partyp(35)=0
      parval(35)=20d0
      paruse(35,5)=1
      paruse(35,11)=1
      paruse(35,14)=1
      paruse(35,15)=1
      paruse(35,16)=1
c
      chpar(36)='ptcen'
      chpdes(36)='min pt for central jet in VBF 3-jet final states'/
     $ /' (used if irapgap = 1 and njets= 3)'
      partyp(36)=0
      parval(36)=parval(30)
      paruse(36,4)=1
      paruse(36,5)=1
      paruse(36,12)=1
c
c     hadronic recoil cut
      chpar(37)='pthrmin'
      chpdes(37)='minimum pt for hadronic recoil (-1: off)'
      partyp(37)=0
      parval(37)=-1.
      paruse(37,1)=1
      paruse(37,2)=1
      paruse(37,3)=1
      paruse(37,4)=1
c      paruse(37,11)=1
c
      chpar(38)='pthrmax'
      chpdes(38)='maximum pt for hadronic recoil (-1: off)'
      partyp(38)=0
      parval(38)=-1.
      paruse(38,1)=1
      paruse(38,2)=1
      paruse(38,3)=1
      paruse(38,4)=1
c      paruse(38,11)=1
c
c     ptphmax cut
      chpar(39)='ptphmax'
      chpdes(39)='maximun pt for photons'
      partyp(39)=0
      parval(39)=-1.
      paruse(39,5)=1
      paruse(39,11)=1
      paruse(39,14)=1
      paruse(39,15)=1
      paruse(39,16)=1
c
c     eta cuts
      chpar(40)='etajmax'
      chpdes(40)='max|eta| for light jets'
      partyp(40)=0
      parval(40)=2.5
      do i=1,nactprc
        paruse(40,i)=1
      enddo
c
      chpar(41)='etabmax'
      chpdes(41)='max|eta| for b quarks (in procs with explicit b)'
      partyp(41)=0
      parval(41)=2.5
      paruse(41,1)=1
      paruse(41,2)=1
      paruse(41,6)=1
      paruse(41,7)=1
      paruse(41,8)=1
      paruse(41,13)=1
      paruse(41,15)=1
      paruse(41,16)=1
c
      chpar(42)='etacmax'
      chpdes(42)='max|eta| for c quarks (in procs with explicit c)'
      partyp(42)=0
      parval(42)=2.5
      paruse(42,1)=1
      paruse(42,2)=1
      paruse(42,6)=1
      paruse(42,7)=1
      paruse(42,15)=1
      paruse(42,10)=1
      paruse(42,16)=1
c
      chpar(43)='etalmax'
      chpdes(43)='max abs(eta) for charged leptons'
      partyp(43)=0
      parval(43)=10d0
      do i=1,4
        paruse(43,i)=1
      enddo
      paruse(43,10)=1
      paruse(43,14)=1
      paruse(43,15)=1
c
      chpar(44)='etaphmax'
      chpdes(44)='max abs(eta) for photons'
      partyp(44)=0
      parval(44)=2.5
      paruse(44,5)=1
      paruse(44,11)=1
      paruse(44,14)=1
      paruse(44,15)=1
      paruse(44,16)=1
c
      chpar(45)='irapgap'
      chpdes(45)='enable central rap-gap in VBF >=2 jet events'
      partyp(45)=1
      parval(45)=0
      paruse(45,4)=1
      paruse(45,5)=1
      paruse(45,12)=1
c
      chpar(46)='etagap'
      chpdes(46)='min rap for 2 "fwd" jets in VBF >=2 jet events'/
     $ /' (used if irapgap = 1)'
      partyp(46)=0
      parval(46)=2.5
      paruse(46,4)=1
      paruse(46,5)=1
      paruse(46,12)=1
c
c     Some more ATLAS-custom parameters...
      chpar(47)='ptj1min'
      chpdes(47)='minimum pt for leading jet (-1: off)'
      partyp(47)=0
      parval(47)=-1
      paruse(47,1)=1
      paruse(47,2)=1
      paruse(47,3)=1
      paruse(47,4)=1
      paruse(47,6)=1
      paruse(47,7)=1
      paruse(47,9)=1
      paruse(47,10)=1
      paruse(47,13)=1
c
      chpar(48)='ptj1max'
      chpdes(48)='maximum pt for leading jet (-1: off'
      partyp(48)=0
      parval(48)=-1
      paruse(48,1)=1
      paruse(48,2)=1
      paruse(48,3)=1
      paruse(48,4)=1
      paruse(48,6)=1
      paruse(48,7)=1
      paruse(48,9)=1
      paruse(48,10)=1
      paruse(48,13)=1
c
c     isolation cuts
      chpar(50)='drjmin'
      chpdes(50)='min deltaR(j-j), deltaR(Q-j) [j=light jet, Q=c/b]'
      partyp(50)=0
      parval(50)=0.7
      do i=1,nactprc
        paruse(50,i)=1
      enddo
c
      chpar(51)='drbmin'
      chpdes(51)='min deltaR(b-b) (procs with explicit b)'
      partyp(51)=0
      parval(51)=0.7
      paruse(51,1)=1
      paruse(51,2)=1
      paruse(51,6)=1
      paruse(51,7)=1
      paruse(51,8)=1
      paruse(51,13)=0
      paruse(51,15)=1
      paruse(51,16)=1
c
      chpar(52)='drcmin'
      chpdes(52)='min deltaR(c-c) (procs with explicit charm)'
      partyp(52)=0
      parval(52)=0.7
      paruse(52,1)=1
      paruse(52,2)=1
      paruse(52,6)=1
      paruse(52,7)=1
      paruse(52,8)=1
      paruse(52,15)=1
      paruse(52,16)=1
c      paruse(52,10)=1
c
      chpar(55)='drlmin'
      chpdes(55)='min deltaR between charged lepton and light jets'
      partyp(55)=0
      parval(55)=0d0
      do i=1,4
        paruse(55,i)=1
      enddo
      paruse(55,10)=1
      paruse(55,14)=1
      paruse(55,15)=1
c
      chpar(56)='drphjmin'
      chpdes(56)='min deltaR between photon and light jets'
      partyp(56)=0
      parval(56)=0.7
      paruse(56,5)=1
      paruse(56,11)=1
      paruse(56,14)=1
      paruse(56,15)=1
      paruse(56,16)=1
c
      chpar(57)='drphlmin'
      chpdes(57)='min deltaR between photon and charged lepton'
      partyp(57)=0
      parval(57)=0.4
      paruse(57,14)=1
      paruse(57,15)=1
c
      chpar(58)='drphmin'
      chpdes(58)='min deltaR between photons'
      partyp(58)=0
      parval(58)=0.7
      paruse(58,5)=1
      paruse(58,11)=1
      paruse(58,14)=1
      paruse(58,15)=1
      paruse(58,16)=1
c
c     dilepton cuts
      chpar(60)='ilep'
      chpdes(60)='Z*/gamma fin state: 0=lept (1 family) 1=nu (3 fam)'
      partyp(60)=1
      parval(60)=0
      paruse(60,2)=1
      paruse(60,4)=1
c
      chpar(61)='mllmin'
      chpdes(61)='min dilepton inv mass'
      partyp(61)=0
      parval(61)=40d0 
      paruse(61,2)=1
      paruse(61,4)=1
c
      chpar(62)='mllmax'
      chpdes(62)='max dilepton inv mass'
      partyp(62)=0
      parval(62)=200d0
      paruse(62,2)=1
      paruse(62,4)=1
c
c     seeds
      chpar(90)='iseed1'
      chpdes(90)='first random number seed (5-digit integer)'
      partyp(90)=1
      parval(90)=12345
      do i=1,nactprc
        paruse(90,i)=1
      enddo

c
      chpar(91)='iseed2'
      chpdes(91)='second random number seed (5-digit integer)'
      partyp(91)=1
      parval(91)=67890
      do i=1,nactprc
        paruse(91,i)=1
      enddo
c
c
c     anomalous couplings
*     changing cosvma between -1 and 1 gives a coupling of the form 
*     cosvma*(V-A) + sinvma*(V+A)   (sinvma= sqrt(1-cosvma^2))
      chpar(100)='cosvma'
      chpdes(100)='top-W coupling: cosvma*(V-A)+sinvma*(V+A)'
      partyp(100)=0
      parval(100)=1d0
      paruse(100,6)=1
      paruse(100,13)=1
      paruse(100,16)=1
c
      chpar(101)='itdec'
      chpdes(101)='forces top decays, with spin-correlations'
      partyp(101)=1
      if(ihrd.eq.6.or.ihrd.eq.8.or.ihrd.eq.13.or.ihrd.eq.16) then
        parval(101)=1
      else 
        parval(101)=0
      endif 
      paruse(101,6)=1
      paruse(101,8)=1
      paruse(101,13)=1
      paruse(101,16)=1
c
      chpar(102)='itopprc'
      chpdes(102)='Selection of single-top process'
      partyp(102)=1
      parval(102)=1
      paruse(102,13)=1
c
      chpar(110)='xlclu'
      chpdes(110)='lambda value for ckkw alpha (match shower alpha)'
      partyp(110)=0
c     default to -1. If it does not get replaced with a positive value
c     in the input file, it will be set equal to xlam in alinit
      parval(110)=-1d0
      do i=1,6
        paruse(110,i)=1
      enddo
      do i=9,12
        paruse(110,i)=1
      enddo
      do i=14,16
        paruse(110,i)=1
      enddo
c
      chpar(111)='lpclu'
      chpdes(111)='loop order for ckkw alpha (match shower alpha)'
      partyp(111)=1
c     default to -1. If it does not get replaced with a positive value
c     in the input file, it will be set equal to nloop in alinit
      parval(111)=-1
      do i=1,6
        paruse(111,i)=1
      enddo
      do i=9,12
        paruse(111,i)=1
      enddo
      do i=14,16
        paruse(111,i)=1
      enddo

c parameters for run with imode=2, npar>150
      chpar(151)='iwdecmode'
      chpdes(151)='W decay modes, in imode=2'
      partyp(151)=1
      parval(151)=1
      if(ihrd.eq.5) parval(151)=11
      paruse(151,1)=1
      paruse(151,3)=1
      paruse(151,5)=1
      paruse(151,10)=1
      paruse(151,13)=1
      paruse(151,14)=1
      paruse(151,15)=1
c
      chpar(152)='izdecmode'
      chpdes(152)='Z->l+l- (ilep=0) decay modes, in imode=2'
      partyp(152)=1
      parval(152)=1
      paruse(152,2)=1
      paruse(152,4)=1
c
      chpar(153)='itdecmode'
      chpdes(153)='top (or t-tbar) decay modes, in imode=2'
      partyp(153)=1
      parval(153)=1
      paruse(153,6)=1
      paruse(153,8)=1
      paruse(153,13)=1
      paruse(153,16)=1
c
      chpar(160)='cluopt'
      chpdes(160)
     $     ='kt scale option. 1:kt propto pt, 2:kt propto mt'
      partyp(160)=1
      parval(160)=1
      do i=1,6
        paruse(160,i)=1
      enddo
      do i=9,12
        paruse(160,i)=1
      enddo
      paruse(160,14)=1
      paruse(160,15)=1
      paruse(160,16)=1
c     seeds for unweighting
      chpar(190)='iseed3'
      chpdes(190)
     $   ='first random number seed for unweighting (5-digit integer)'
      partyp(190)=1
      parval(190)=12345
      do i=1,nactprc
        paruse(190,i)=1
      enddo

c
      chpar(191)='iseed4'
      chpdes(191)
     $   ='second random number seed for unweighting (5-digit integer)'
      partyp(191)=1
      parval(191)=67890
      do i=1,nactprc
        paruse(191,i)=1
      enddo
c
c hidden parameters
      chpar(195)='impunw'
      chpdes(195)
     $   ='impunw=2: maxwgt<0, impunw=1: ask user, impunw=0: default'
      partyp(195)=1
      parval(195)=0
      do i=1,nactprc
        paruse(195,i)=2
      enddo
c
      do i=1,npar
        parlen(i)=aluisl(chpar(i))
      enddo
      return
c
      entry Alppar(n)
c print definitions and current values of parameters
      if(n.eq.2.or.n.eq.4.or.n.eq.5) then
        call alugun(iunit)
        if(n.eq.5) then
          open(iunit,file='prc.list',status='unknown')
        else
          open(iunit,file='par.list',status='unknown')
        endif
      else
        iunit=6
      endif
      if(n.lt.5) then
        write(iunit,*) '------'
        write(iunit,*) 'hard process code (not to be changed):'
        write(iunit,*) 'ihrd=',ihrd
        do i=1,npar
          if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
            itmp=aluisl(chpdes(i))
            write(iunit,*) '------'
            write(iunit,*) chpdes(i)(1:itmp),':'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
            else
              write(iunit,*) chpar(i)(1:itmp),'=',int(parval(i))
            endif
            if(i.eq.4) call prntsf(iunit)
            if(i.eq.5) call alhsca(ihrd,iunit)
          endif
        enddo
      elseif(n.eq.5) then
c print list of all processes
        write(iunit,*) '======'
        write(iunit,*) 'List of processes'
        do i=1,nactprc
          if(chprc(i)(1:6).ne.'unavai') then
            write(iunit,'(I2,1x,a)') i,chprc(i)
          endif
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'List of parameters'
        do 490 i=1,npar
          do j=1,nactprc
            if(paruse(i,j).eq.2) goto 480
          enddo
          if(chpar(i).ne.'***') then
            write(iunit,*) '------'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) i,' ',chpar(i),parval(i)
            else
              write(iunit,*) i,' ',chpar(i),int(parval(i))
            endif
            itmp=aluisl(chpdes(i))
            write(iunit,*) chpdes(i)(1:itmp),':'
            write(iunit,'(20(I2,1x))') (paruse(i,j),j=1,nactprc)
          endif
 480      continue
 490    enddo
      write(iunit,*) '======'
        write(iunit,*) 'Scale setting choices for each process:'
        do 500 j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 495
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alhsca(j,iunit)
 495      continue
 500    enddo
        write(iunit,*) '======'
        write(iunit,*) 'imode=1 process-specific documentation'
        do 600 j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 590
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alpdoc(j,iunit,1)
 590      continue
 600    enddo
        write(iunit,*) '======'
        write(iunit,*) 'imode=2 process-specific documentation'
        do 700 j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 690
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alpdoc(j,iunit,2)
 690      continue
 700    enddo
        write(iunit,*) '======'
        write(iunit,*) 'Structure function menu'
        call prntsf(iunit)
      elseif(n.eq.6) then
        do i=151,npar
          if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
            itmp=aluisl(chpdes(i))
            write(iunit,*) '------'
            write(iunit,*) chpdes(i)(1:itmp),':'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
            else
              write(iunit,*) chpar(i)(1:itmp),'=',int(parval(i))
            endif
          endif
        enddo
      endif
      if(imode.lt.5) then
        write(6,*) ' '
        call alpdoc(ihrd,iunit,imode)
      endif
      if(iunit.ge.10) close(iunit)
      end

c-------------------------------------------------------------------
      subroutine Alpdoc(ihrd,iunit,imode)
c     print information regarding hard process ihrd
c-------------------------------------------------------------------
      implicit none
      integer ihrd,iunit,imode
      if(imode.le.1.or.imode.gt.2) then
C Z+jets and VBJET PROCESSES:
        if(ihrd.eq.4.or.ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     & 'GENERATION OF RAPIDITY-GAP CONFIGURATIONS:'
     &,'If njets>=2 the option exists to generate events with fwd jets.'
     &,'Select irapgap=1 to force two jets to satisfy the constraints:'
     &,'|eta(j1)|>etagap   |eta(j2)|>etagap  eta(j1)*eta(j2)<0'
     &,'In the case of njets=3, irapgap=1 assumes the third jet to be'
     &,'central, with abs(eta)<etagap and with pt>ptcen'
          write(iunit,*) ' '
        endif 
        if(ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     &         'SELECTION OF DECAY MODES:'
     &,'the choice of decay modes for the Z bosons should be performed'
     &,'already in imode=1 running (the EW couplings of the Z depend on'
     &,'flavour). Input the integer string "zfstate" describing  the '
     &,'decay modes of the individual Zs:'
     &,'1: Z->nu nubar (summed over all flavours)' 
     &,'2: Z->l+l- (summed over all charged leptons)'
     &,'3: Z->q qbar (summed over all quark flavours < top)'
     &,'4: Z->b bbar'
     &,'5: Z->all f fbar modes'
     &,'Example: input "zfstate 24" for ZZ -> l+l- b bar'
     &,'Example: input "zfstate 25" for ZZ -> l+l- + (Z->anything)'
     &,'Example: input "zfstate 234" for ZZZ -> l+l- q qbar b bar'
     &,'NB: The decay modes of the W boson are entered in imode=2'
        endif
        if(ihrd.eq.13) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     &         '4 single-top processes can be selected:'
          write(iunit,*)'itopprc=1: t+q (njets=0)'
          write(iunit,*)'itopprc=2: t+b (njets=0)'
          write(iunit,*)'itopprc=3: t+W(W->f fbar")+jets (njets=0,1)'
          write(iunit,*)"itopprc=4: t+b+W(W->f fbar')+jets (njets=0,1)"
        endif
      elseif(imode.eq.2.or.imode.gt.2) then
        write(iunit,*) ' '
        if(ihrd.eq.1.or.ihrd.eq.3.or.ihrd.eq.10) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
          write(iunit,*) '5: q q''bar'
          write(iunit,*) '6: fully inclusive'
        elseif(ihrd.eq.2.or.ihrd.eq.4) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2 and ilep=0:'
          write(iunit,*)
     $         'select Z->l+l- decay modes, set izdecmode to:'
          write(iunit,*) '1: e e'
          write(iunit,*) '2: mu mu'
          write(iunit,*) '3: tau tau'
          write(iunit,*) '4: ell+ ell-'
        elseif(ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*)
     $         'select decay mode for each W, set iwdecmode as:' 
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu '
          write(iunit,*) '5: q qbar'' '
          write(iunit,*) '6: fully inclusive'
          write(iunit,*)
     $         'E.g.: input iwdecmode 24 for WW-> mu l nu_mu nu_l' 
          write(iunit,*
     $         )'E.g.: input iwdecmode 45 for WW-> l nu q qbar'''  
          write(iunit,*)
     $     'E.g.: input iwdecmode 126 for WWW-> e nue mu nu(mu) ff''bar'
        elseif(ihrd.eq.6) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select top decay modes, set itdecmode to:'
          write(iunit,*) '1: e nu b bbar + 2 jets'
          write(iunit,*) '2: mu nu b bbar + 2 jets'
          write(iunit,*) '3: tau nu b bbar + 2 jets'
          write(iunit,*) '4: e/mu/tau nu b bbar + 2 jets'
          write(iunit,*) '5: l nu l'' nu  b bbar (l,l''=e/mu/tau)'
          write(iunit,*) '6: b bbar + 4 jets'
          write(iunit,*) '7: fully inclusive'
        elseif(ihrd.eq.7) then
        elseif(ihrd.eq.8) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select top decay modes, set itdecmode to:'
          write(iunit,*) '1: e nu b bbar + 2 jets'
          write(iunit,*) '2: mu nu b bbar + 2 jets'
          write(iunit,*) '3: tau nu b bbar + 2 jets'
          write(iunit,*) '4: e/mu/tau nu b bbar + 2 jets'
          write(iunit,*) '5: l nu l'' nu  b bbar (l,l''=e/mu/tau)'
          write(iunit,*) '6: b bbar + 4 jets'
          write(iunit,*) '7: fully inclusive'
        elseif(ihrd.eq.9) then
        elseif(ihrd.eq.11) then
        elseif(ihrd.eq.12) then
        elseif(ihrd.eq.13) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select top decay modes, set itdecmode to:'
          write(iunit,*) '1: e nu b'
          write(iunit,*) '2: mu nu b'
          write(iunit,*) '3: tau nu b'
          write(iunit,*) '4: e/mu/tau nu b'
          write(iunit,*) '5: b + 2 jets'
          write(iunit,*) '6: fully inclusive'
c          read(iunit,*) itdecmode
c          if(itopprc.ge.3) then
          write(iunit,*) 'if itopprc.ge.3, then'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
          write(iunit,*) '5: 2 jets'
          write(iunit,*) '6: fully inclusive'
c            read(iunit,*) iwdecmode
c          endif
        elseif(ihrd.eq.14.or.ihrd.eq.15) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
        endif
      endif
      end

c-------------------------------------------------------------------
      subroutine Alspar
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
c
      ih2=parval(2)
      ebeam=parval(3)
      ndns=parval(4)
      iqopt=parval(5)
      qfac=parval(6)
      ickkw=parval(7)
      ktfac=parval(8)
      njets=parval(10)
      ihvy=parval(11)
      ihvy2=parval(12)
      nw=parval(13)
      nz=parval(14)
      nh=parval(15)
      nph=parval(16)
      zfstate=parval(17)
      mc=parval(20)
      mb=parval(21)
      mt=parval(22)
      mh=parval(23)
      ptjmin=parval(30)
      ptbmin=parval(31)
      ptcmin=parval(32)
      ptlmin=parval(33)
      metmin=parval(34)
      ptphmin=parval(35)
      ptcen=parval(36)
      pthrmin=parval(37)
      pthrmax=parval(38)
      ptphmax=parval(39)
      etajmax=parval(40)
      etabmax=parval(41)
      etacmax=parval(42)
      etalmax=parval(43)
      etaphmax=parval(44)
      irapgap=parval(45)
      etagap=parval(46)
      ptj1min=parval(47)
      ptj1max=parval(48)
      drjmin=parval(50)
      drbmin=parval(51)
      drcmin=parval(52)
      drlmin=parval(55)
      drphjmin=parval(56)
      drphlmin=parval(57)
      drphmin=parval(58)
      ilep=parval(60)
      mllmin=parval(61)
      mllmax=parval(62)
      iseed(1)=parval(90)
      iseed(2)=parval(91)
      cosvma=parval(100)
      itdec=parval(101)
      itopprc=parval(102)
c
c      if(parval(110).gt.0) xlclu=parval(110)
c      if(parval(111).gt.0) lpclu=parval(111)
      xlclu=parval(110)
      lpclu=parval(111)
c parameters for imode=2, n>150
      iWdecmode=parval(151)
      iZdecmode=parval(152)
      itdecmode=parval(153)
      cluopt=parval(160)
      iseed2(1)=parval(190)
      iseed2(2)=parval(191)
      impunw=parval(195)
      end

c-------------------------------------------------------------------
      subroutine alinit
c     process the input parameters and fills the couplings constants
c     needed by ALPHA
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
#ifdef USE_MPI
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      include 'mpif.h'
#endif
c local variables
      double precision aluhiw,alfas,alfas_clu
      integer i,iret,iunit,aluisl,itmp,retval
      character*50 ewopt(0:3)
c
      if(iewopt.eq.0) then
c     input GF, sin-thetaW, alphaem
c     calculate the rest
         gw=sqrt(4.d0*pi*aem)/stw
         mw=gw/sqrt(4d0*sqrt(2d0)*gfermi)
         ctw=sqrt(1d0-stw**2)
         mz=mw/ctw
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.1) then
c     input mW, GF and sin_thetaW
c     calculate the rest
         ctw=sqrt(1-stw**2)
         mZ=mW/ctw
         gw=sqrt ( 4d0*sqrt(2d0)*gfermi* mw**2 )
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.2) then
c     input mZ, alpha(em) and sin_thetaW
c     calculate the rest
         ctw=sqrt(1-stw**2)
         mw=mz*ctw 
         gw=sqrt(4.d0*pi*aem)/stw
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.3) then
c     input mZ, mW and GF
c     calculate the rest
         gw=sqrt ( 4d0*sqrt(2d0)*gfermi* mw**2 )
         ctw=mw/mz
         stw=sqrt(1-ctw**2)
         gbar=gw/ctw 
         g1=gbar*stw
         aem=gw**2/(4.d0*pi)*stw**2
      else
        write(6,*)
     $       'option for selection of EW parameters not specified'
        stop
      endif
      wwid = 3d0/16d0/pi * gw**2 * mw
      zwid = 
c     leptons
     +     3*(1-2*stw**2+4*stw**4) +
c     u,d,c,s
     +     2*3*(1-2*stw**2+20d0/9d0*stw**4) +
c     b
     +     3*sqrt(1-(2*mb/mz)**2)*(1d0-(2d0*mb/mz)**2 + (-1d0+4d0/3d0
     $     *stw**2)**2*(1d0+2d0*(mb/mz)**2))/4d0
      zwid = 1/(48*pi) * mz * gbar**2 * zwid
      hwid = aluhiw(mh,mw,mz,mt)
c initialise QCD and PDF parameters
 1    call pdfconv(ndns,nmnr,pdftyp)
      call pdfpar(nmnr,ih1,xlam,sche,nloop,iret)
      if(iret.eq.1) then
         write(6,*) 'pdf set',ndns,' not available, re-input ndns:'
         call prntsf(6)
#ifdef USE_MPI
c         read(9,*) ndns
c         just exit because the input file has incorrent PDF defined
          call exit(-1)
#else
         read(5,*) ndns
#endif
         goto 1
      endif
      as=alfas(mz**2,xlam,nloop,-1)
c glu-glu-higgs coupling in mt-> infinity limit
      ggh= as/3.d0/246.d0/pi
c      ggh= as/3.d0/246.d0/pi*(1.d0+11.d0/4.d0*as/pi)
c lambda and loop order in alphas for reweighting after clustering
c     Set to xlam and nloop if default was not changed in input
      if(xlclu.lt.0) xlclu=xlam
      if(lpclu.lt.0) lpclu=nloop
c fill in  alpha common with mass and coupling parameters
      do i=1,100
         apar(i)=0d0
      enddo
c
      apar(1)=mz
      apar(2)=zwid
      apar(3)=mw
      apar(4)=wwid
      apar(5)=mh
      apar(6)=hwid
      apar(14)=mc
      apar(15)=mb
      apar(16)=mt
      apar(51)=gw
      apar(52)=g1
      apar(53)=ctw
      apar(54)=stw
      apar(55)=gbar
      apar(61)=ggh
c
      do i=1,24
        amass(i)=0
      enddo
      amass(4)=mc
      amass(5)=mb
      amass(6)=mt
      amass(11)=mlep(1)
      amass(13)=mlep(2)
      amass(15)=mlep(3)
      amass(23)=mz
      amass(24)=mw
c
      roots=2d0*ebeam
      s=roots*roots
      ptjmax=0.8*ebeam
      ptbmax=ptjmax
      ptcmax=ptjmax
c
      ewopt(0)='input GF, sin-thetaW, alphaem, calculate the rest'
      ewopt(1)='input mW, GF, sin-thetaW, calculate the rest'
      ewopt(2)='input mZ sin-thetaW, alphaem, calculate the rest'
      ewopt(3)='input mW, mZ, GF calculate the rest'
      do i=1,2
        if(i.eq.1) iunit=6
        if(imode.eq.2) goto 100
        if(i.eq.2) iunit=niosta
        write(iunit,*) ' '
        write(iunit,*) ' RUNNING PARAMETERS'
        write(iunit,*) '-------------------'
        write(iunit,*) ' '
        write(iunit,*) '        Electroweak parameters:'
        write(iunit,*) 'iewopt=',iewopt,':'
        write(iunit,*) ewopt(iewopt)
        write(iunit,*) 'M(W)=',mw,' Gamma(W)=',wwid
        write(iunit,*) 'M(Z)=',mz,' Gamma(Z)=',zwid
        write(iunit,*) 'M(H)=',mh,' Gamma(H)=',hwid
        write(iunit,*) ' gW=',gw,'; sin^2(thetaW)=',stw**2,
     $       '; 1/aem(mZ)=',1/aem
        write(iunit,*) ' '
        write(iunit,*) '        Quark masses:'
        write(iunit,*) 'm(top)=',mt,' m(b)=',mb
        write(iunit,*) ' '
        write(iunit,*) '       Beams'' parameters:'
        if(ih2.eq.1) write(iunit,*) 'beam1=proton, beam2=proton'
        if(ih2.eq.-1) write(iunit,*) 'beam1=proton, beam2=antiproton'
        write(iunit,*) 'Ebeam=',ebeam,' PDF set=',pdftyp
        write(iunit,*) 'as(MZ)[nloop=',nloop,'] = ',as
        write(iunit,*) ' '
      enddo
      return
c write out parameters for herwig
#ifdef USE_MPI
 100  if(mpirank.eq.0.or.imode.ne.2) then
         write(niosta,'(a)') '************** run parameters '
         write(niosta,9991) ihrd,' ! hard process code'
         write(niosta,9992) mc,mb,mt,mw,mz,mh,' ! mc,mb,mt,mw,mz,mh'
 9991    FORMAT(I4,A20)
 9992    FORMAT(6F8.3,A20)  
         do i=1,npar
            if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
               itmp=aluisl(chpar(i))
               write(niosta,*) i,parval(i),'  ! ',chpar(i)(1:itmp)
            endif
         enddo
         write(niosta,'(a)') '************** end parameters '
      endif
#else
 100  write(niosta,'(a)') '************** run parameters '
      write(niosta,9991) ihrd,' ! hard process code'
      write(niosta,9992) mc,mb,mt,mw,mz,mh,' ! mc,mb,mt,mw,mz,mh'
 9991 FORMAT(I4,A20)
 9992 FORMAT(6F8.3,A20)  
      do i=1,npar
        if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
          itmp=aluisl(chpar(i))
          write(niosta,*) i,parval(i),'  ! ',chpar(i)(1:itmp)
        endif
      enddo
      write(niosta,'(a)') '************** end parameters '
c
#endif

      end


c-------------------------------------------------------------------
      subroutine aligrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      integer nvar,nv
      common/psopt/nvar,nv
      real*8 totalmax,tmpfac
      common/mxfact/totalmax,tmpfac
      character*50 tmpstr
c local variables
      integer n,nit,iunit
      double precision ntmp
#ifdef USE_MPI
      include 'mpif.h'
      integer filestat
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      real*8 al1,bet1
      integer nct,nx1,nx2,nct1,init,j,k,nrd,i
      integer ncmax,maxn
      parameter (ncmax= 1000, maxn= 100)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
#endif
*
c-    If  unweighting run:
*
      if(imode.eq.2) then
         avgrew=0d0
         colgen=.true.
         evopt=.false.
         evgen=.true.
         return
      endif
*
c     input grid, if so required
      if(igrid.ne.0) then
c first of all save old grids
c grid 1:
        call alugun(iunit)
        if(igrid.eq.1) then
          call alustc(fname,'.grid1',tmpstr)
#ifdef USE_MPI
          open(unit=iunit,file=tmpstr,status='unknown',iostat=filestat)
          if(filestat /= 0) then
             write(6,*) 'Error Opening file:',tmpstr
             call exit(-1)
          endif
#else
          open(unit=iunit,file=tmpstr,status='unknown')
#endif
          call readgrid(iunit)
          close(iunit)
          call alustc(fname,'.grid1-old',tmpstr)
#ifdef USE_MPI
          open(unit=iunit,file=tmpstr,status='unknown',iostat=filestat)
          if(filestat /= 0) then
             write(6,*) 'Error Opening file:',tmpstr
             call exit(-1)
          endif
#else
          open(unit=iunit,file=tmpstr,status='unknown')
#endif
          do n= 1,nvar+1
            call grid1W(iunit,n)
          enddo
#ifdef USE_MPI
          write (iunit,20,iostat=filestat) totalmax,tmpfac
          if(filestat /= 0) then
             write(6,*) 'Error Opening file:',tmpstr
             call exit(-1)
          endif
#else
          write (iunit,20) totalmax,tmpfac
#endif
          close(iunit)
          call alustc(fname,'.grid1',tmpstr)
        elseif(igrid.eq.2) then
c grid 2
#ifndef USE_MPI
          call alustc(fname,'.grid2',tmpstr)
           open(unit=iunit,file=tmpstr,status='unknown')
          call readgrid(iunit)
          close(iunit)
          call alustc(fname,'.grid2-old',tmpstr)
          open(unit=iunit,file=tmpstr,status='unknown')

          do n= 1,nvar+1
            call grid1W(iunit,n)
          enddo
          
          write (iunit,20) totalmax,tmpfac
          close(iunit)
#endif
          call alustc(fname,'.grid2',tmpstr)
        endif
c now open the grid required for the run:
#ifdef USE_MPI

        if(mpirank.eq.0) then
          open(unit=iunit,file=tmpstr,status='unknown')
          if(filestat /= 0) then
             write(6,*) 'Error Opening file:',tmpstr
             call exit(-1)
          endif
          call readgrid(iunit)
          close (iunit)
        endif
        call MPI_BCAST(al1,ncmax*maxn,MPI_REAL8,0,MPI_COMM_WORLD,
     $                 filestat)
        if(filestat.ne.0) then
           write(6,*) 'Error broadcasting al1'
           call exit(-1)
        endif
        call MPI_BCAST(nx1,maxn,MPI_INTEGER,0,MPI_COMM_WORLD,
     $                 filestat)
        if(filestat.ne.0) then
           write(6,*) 'Error broadcasting nx1'
           call exit(-1)
        endif
        call MPI_BCAST(bet1,(1+ncmax)*maxn,MPI_REAL8,0,MPI_COMM_WORLD,
     $              filestat)
        if(filestat.ne.0) then
           write(6,*) 'Error broadcasting bet1'
           call exit(-1)
        endif
        call MPI_BCAST(totalmax,1,MPI_REAL8,0,MPI_COMM_WORLD,
     $              filestat)
        if(filestat.ne.0) then
           write(6,*) 'Error broadcasting totalmax'
           call exit(-1)
        endif
        call MPI_BCAST(tmpfac,1,MPI_REAL8,0,MPI_COMM_WORLD,
     $              filestat)
        if(filestat.ne.0) then
          write(6,*) 'Error broadcasting tmpfac'
          call exit(-1)
        endif
#else
        open(unit=iunit,file=tmpstr,status='unknown')
        call readgrid(iunit)
        close (iunit)
#endif
      endif
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
*     
c-    Optimize grid with first iterations, if required:
*
      ntmp=maxev
      if(nopt.ge.1) then
         evgen=.false.          ! do NOT store event info (e.g. histograms)
         evopt=.true.           ! perform optimization procedures
         colgen=.false.         ! do NOT generate colour flows
         maxev=nopt
         do nit=1,niter
            call alegen
         enddo
      endif
*
c-    Reset generation parameters, move to event generation
*
      maxev=ntmp
      evgen=.true.              ! store event info (e.g. histograms)
      evopt=.true.              ! perform optimization procedures
      colgen=.false.            ! do NOT generate colour flows
      end


c-------------------------------------------------------------------
      subroutine alegen
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
#ifdef USE_MPI
      include 'mpif.h'
      integer mpirank,mpiworldsize,retval
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      real timestamp
      real*8 sigerr_sum,sigerr_global
      real*8 avgwgt_sum,avgwgt_global
      integer unwev_sum
#endif

c     commons
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      real*8  psum(1000), p2sum(1000)
      common/subproc/psum,p2sum
      integer maxn,nct,nx1,nx2
      parameter (maxn= 100)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer nvar,nv
      common/psopt/nvar,nv
c     
c     local variables
      integer keeps,filestat
      real *8 djpd,djg
      real *8 totalmax,tmpfac,totalmin,ncount,ncount1,extra
      integer nmon
      data nmon/100000/
      data extra/1.5d0/
      common/mxfact/totalmax,tmpfac
      integer navgusd,maxiter,navgin
      data avgwgt/-1.d0/
      data navgin/40/,maxiter/20/
      real *8 rn1,rn2,pswgt,xlum,total,alphap(4,maxpar)
      real *8 wgcol,dummy,ranget,ranset,xrn
      integer colfl(2*maxpar)
      integer i,n,iflag,ip(maxpar),ipinv(maxpar),
     $     afl(maxpar),alphafl(maxpar),instate(2),ndummy
c     local counters
      real *8 ipscount,nev
      real *8 weff,weff1
      real *8 ipscounti(1000),ispncounti(1000),ni(1000)
      save ipscounti,ispncounti,ni
      save ipscount
      real *8 gtotal,g2total,effps,effev
      real *8 sigerr,cumavg,cumerr,cumsig
      real*8  gsum(1000),g2sum(1000),gg2sum(1000)
      data cumsig/0d0/,cumerr/0d0/,cumavg/0d0/
c     
      real *8 newwmax(5)
      real *8 rspin,evtme,evtme0,evwgt,djproc
c     local file naming for grids
      character*550 tmpstr
c
      integer iunit
c     random numbers defining the event
      integer jseed(2)
      data jseed/2*0/
c     local weight variables
      real rwgt,rx1
      double precision nevwgt
      logical unwgt
      real *8 avgwgt_in,g2total_in,maxwgt_in
      integer mxnavg
      parameter (mxnavg=100)
      real*8 refwgt,tmp,avgsig(mxnavg)
      double precision wusr
c     local spin variables
      integer hel(100)
      integer iavg,iavg_store,iavst(mxnavg)
c     clustering reweight
      double precision rewgt,shat,alfas,alfas_clu,rn,nrewgt
      integer color(2*maxpar)     !coulor string
      common/colore/color
c     local momentum variables
      integer maxdec
      parameter (maxdec=40)
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
      integer ierr

c     
c     for debugging:    
c     character*3 pname(-12:12)
c     data pname/'nb ','e+ ',5*'   ','bb ','cb ','sb ','ub ','db ','gl '
c     $     ,'d  ','u  ','s  ','c  ','b  ',5*'   ','e- ','nu '/
c     
c     reset counters
      real*8  aextra(0:8)
      integer imaxiter(0:8)  
      data imaxiter/1,4,6,10,17,28,45,49,49/
      data aextra/1.d0,1.15d0,1.3d0,1.45d0,1.6d0,
     $     1.75d0,1.9d0,2.05d0,2.2d0/
*     
      real*8 qref(4)/221.d0,0.d0,0.d0,21.d0/,qout(4),msqref/48400.d0/
      integer j1,j2
      logical flg90
      common/f90/flg90          ! to avoid spurious divergences in temporal gauge
*     
      integer jpr
      common/jp/jpr
      save qref,msqref
*     
      if(imode.le.1.and.maxev.eq.0) return
c initialise averaging over multiple color/spin configurations
      if (imode.eq.2) nprtns= 8 
      maxiter= imaxiter(nprtns)
      extra  = aextra(nprtns) 
      navgin = 1
      if (maxiter.ne.1) navgin = 2*maxiter
      if(navgin.gt.mxnavg) then
        write(6,*)
     $       'bookeeping of multiple colour/spin calls not properly',
     $       ' initialised: '
        write(6,*)
     $       'number of calls exceeds ',mxnavg 
        write(6,*) 'Redimension array ``avgsig'' in routine ``alegen'''
        stop
      endif
      gtotal=0d0
      g2total=0d0
c     reset counters for individual processes to 0
      do i=1,1000
        gsum(i)=0d0
        g2sum(i)=0d0
      enddo
      ipscount=0d0
      do i= 1,jprocmax
        ipscounti(i) = 0d0
        ispncounti(i)= 0d0
        ni(i)        = 0d0
      enddo
      maxwgt=-1d0
      unwev=0
      nevwgt=0d0
      nrewgt=0d0
c     if writing events to a file for later unweighting,
c     dump the grids first
      if(imode.eq.1.and.evgen) then
        do n= 1,nvar+1
          call grid1W(niopar,n)
        enddo
        write (niopar,20) totalmax,tmpfac
      endif
c     if reading events from a file for unweighting,
c     input the grids first
      if(imode.eq.2) then
        call readgrid(niopar)
        read(niopar,*)
        read(niopar,10) maxev,avgwgt_in,g2total_in,maxwgt_in
        read(niopar,*)
        read(niopar,'(e12.4)') tmp
c user sets the maximum weight:
        read(niopar,'(5(e12.6,1x))') newwmax
        do i=1,5
          if(newwmax(i).gt.maxwgt_in) newwmax(i)=maxwgt_in
        enddo
        if(impunw.eq.1) then
#ifdef USE_MPI
 9        if(mpirank == 0) then
            read(9,*) maxwgt
          endif
        call MPI_BCAST(maxwgt,1,MPI_REAL8,0,MPI_COMM_WORLD,
     $              retval)
        if(retval /= 0) then
          write(6,*) 'Error broadcasting maxwgt'
          call exit(-1)
        endif
#else
          write(6,*) ' '
          write(6,*)
     $         'Expected approx # of unweighted events:',int(avgwgt_in
     $         /maxwgt_in*tmp)
 8        FORMAT(1X,A,100(/,1X,A))
          write(6,8)
     $         'To improve unweigthing efficiency, you can rescale',
     $         'the true maximum weight.',
     $         'This may slightly bias a fraction of the events.',
     $  'The following factors lead to the quoted eff improvements',
     $  'affecting at most the quoted fraction of events:'
          write(6,'(a)')
     $'bias:      1%          2%           3%           4%           5%'  
          write(6,'(a,5(e12.6,1x))') 'rescale: ',(newwmax(i)/maxwgt_in,i
     $         =1,5)
          write(6,'(a,5(e12.6,1x))') 'eff impr:',(maxwgt_in
     $         /newwmax(i),i=1,5)
 9        write(6,*) ' '
          write(6,*) 'input desired increase in efficiency',
     $         ' (1 for the default max wgt)'
          read(5,*) maxwgt
#endif
c     
c     keep fraction of affected events <5% :
          if(maxwgt.gt.maxwgt_in/newwmax(5)) maxwgt=maxwgt_in/newwmax(5)
c     keep efficiency boost < 100, not to be biased by pre-unweighting:
          maxwgt=min(maxwgt,1d2)
          write(6,*) 'implemented increase in efficiency =',maxwgt
          maxwgt=maxwgt_in/maxwgt
        elseif(impunw.eq.0) then
c          write(6,'(a,5(e12.6,1x))') 'rescale: ',(newwmax(i)/maxwgt_in,i
c     $         =1,5)
c          write(6,'(a,5(e12.6,1x))') 'eff impr:',(maxwgt_in
c     $         /newwmax(i),i=1,5)
          maxwgt=maxwgt_in/newwmax(5)          
          maxwgt=min(maxwgt,1d2)
c          write(6,*) 'implemented increase in efficiency =',maxwgt
          maxwgt=maxwgt_in/maxwgt
        elseif(impunw.eq.2) then
          maxwgt=-1
        else
          write(6,*) 'unspecified choice of maximum weight setting'
          stop
        endif
      endif 
 10   format(f15.1,1x,3(e12.6,1x))
 11   format(e12.6)
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
c     generate events
      write(6,*) ' '
      if(imode.le.1) then
        write(6,*) 'starting generation of',maxev,' events'
        write(niosta,*) ' '
        write(niosta,*) 'starting generation of',maxev,' events'
      else
        write(6,*) 'starting scan/unweighting of',maxev,' events'
      endif
#ifdef USE_MPI
      call cpu_time(timestamp)
      write(6,*) 'RANK ',rankstring,' BEGIN EVTGEN AT ',timestamp
#endif
      iflag= 0
      n=0
      ncount= 0.d0
      ncount1= 0.d0
c     start loop over events, between 99 and 100
      nev=1d0
 99   continue
      n=n+1
c     reset local integer counter, to avoid overruns, Used only for
C     monitoring
      if(n.gt.1000000000) n=n-1000000000
c     generate kinematical configuration, return phase-space weight

 1    if(imode.eq.0) then
        dummy= ranget(jseed)
        continue
      elseif(imode.eq.1) then
        dummy= ranget(jseed)
      elseif(imode.eq.2) then 
        read(niowgt,15,end=1000) jseed,iavg_store,rwgt,rx1
        evwgt=dble(rwgt)
        if(.not.unwgt(evwgt,maxwgt)) goto 100
        dummy= ranset(jseed)
        unwev=unwev+1
      endif
 15   format(i12,1x,i12,1x,i4,2(1x,e12.6))

c     MC over channels
c     generate jproc in range 1-jprocmax, 
c     return inverse jacobian djproc

      call onedimbin(1,jproc,djproc,1,ndummy,dummy)
      ni(jproc)= ni(jproc)+1d0
      jpr= jproc

      total= 0.d0
      
 16   call phspace(iflag,pswgt,djpd,djg)
      if(iflag.eq.1) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        goto 16
      endif
c     evaluate parton densities
      call setpdf
c     evaluate parton luminosities, alpha_s factor,  select flavours and
c     apply possible flavour-dependent cuts
      call selflav(jproc,xlum,afl)
      if(xlum.le.0d0) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        goto 16
      endif
c     Apply user-defined selection cuts;
c     wusr allows the user to reweight the event as a function of the
c     phase-space cuts
      call usrcut(iflag,wusr)
      pswgt=pswgt*wusr
      if(iflag.eq.1) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        goto 16
      endif
c     events passed all kinematical cuts
      weff= ni(jproc)/(ni(jproc)+ipscounti(jproc))

c     check that weighted-events file was not contaminated
      if(imode.eq.2) then
        if(abs(rx1-real(x1))/(rx1+real(x1)).gt.1e-5) then
          write(6,*) 'Weighted-events file has problems:'
          write(6,*) 'x1 from phase-space=',real(x1)
          write(6,*) 'x1 from file =',rx1
          write(6,*) 'skip event',unwev
          unwev=unwev-1
          goto 100
        endif
      endif
c     setup input for alpha computation:
      call setalp(npart,afl,p,alphafl,alphap,instate,ip,ipinv)
c     generate colour      
      tmp=0d0
      navg    = navgin
      totalmin=-1.d0
      if (avgwgt.lt.0.d0.and.igrid.eq.0) then
        navg    = 1
        totalmax= 1.d25
        tmpfac  = 100.d0
      else
        if (.not.evgen) totalmax= avgwgt*tmpfac
      endif
      if(imode.eq.2) navg=navgin
      navgusd = 0
      keeps= 0
      do 111 iavg=1,navg
 2      call randa(rn1)
        call randa(rn2)
        call selcol(alphafl,instate,rn1,rn2,iflag)
        if(iflag.eq.1)   goto 2
c     generate spin
 17     call randa(rspin)
        call selspin(rspin,hel)
        call fltspn(alphafl,hel,instate,iflag)
        if(iflag.eq.1) then
          if(keeps.eq.0) ispncounti(jproc)= ispncounti(jproc)+1d0
          goto 17
        endif
c     non-zero spin configuration selected, efficiency updated
        weff1= weff*ni(jproc)/(ni(jproc)+ispncounti(jproc))     
c     
c     reevaluate alphas using ktmin as a scale
        if(ickkw.eq.1) then
c          write(6,*) ' '
c          write(6,*) (alphafl(i),i=1,npart)
c          write(6,*) (color(2*i-1),i=1,npart)
c          write(6,*) (color(2*i),i=1,npart)
          call setcol(npart,ip,ipinv,color,ifl,icu)
c          write(6,*) 'qsq=',qsq
          call cktmin
c          write(6,*) 'kres=',kres
          asmax=alfas_clu(kres,xlclu,lpclu,-1)
c     multiply alphas by a rescaling factor, to allow separation of qcd
C     and ew
c     contributions to jet production: 
c     
c     sig=1/(resc**N) * sig(as->as*resc)
c     
c     for processes where maximum power of alphas is N.
c     resc>>1 will suppress ew processes
c     resc<<1 will suppress qcd processes
          apar(56)=sqrt(4*pi*as*resc)
        endif
        
        if(imode.eq.2.and.iavg.eq.iavg_store) then 
          total=evwgt
c debug
c          write(2,*) ' '
c          write(2,*) (alphafl(i),i=1,npart)
c          write(2,*) (color(2*i-1),i=1,npart)
c          write(2,*) (color(2*i),i=1,npart)
c          write(2,*) 'ktmin=',sqrt(kres),' Q=',sqrt(qsq)
c          write(2,*) 'as(Q)/as(kres)=',alfas(qsq,xlam,nloop,-1)
c     $         /alfas(kres,xlam,nloop,-1) 
c end debug
          goto 112
        elseif(imode.le.1) then
          resonance= 'n'
          flg90=.false.
          call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     $         ,colfl)
          if (resonance.eq.'y') then
c     print*,'taglio',evtme
            evtme= 0.d0
          endif
          if(flg90.and.resonance.ne.'y') then
            do j1=1,npart
              call boost(1,msqref,qref,qout,alphap(1,j1))
              do j2=1,4
                alphap(j2,j1)=qout(j2)
              enddo
            enddo
            flg90=.false.
            call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     +           ,colfl)
          endif
          ncount= ncount+1
          ncount1= ncount1+1
c     
c     include MC sum over jprocs, gev->pb
          total = dble(jprocmax)/djproc * pbfac * evtme
c     combine with parton luminosities, phase-space weight and
c     efficiency
          total = total*xlum*pswgt*weff1
c
c     rescale alphas for the light jets
          if(ickkw.eq.1) total=total*(asmax/as)**naspow
          if ((iavg.eq.1).and.
     .         (total.le.totalmax).and. 
     .         (total.ge.totalmin)) then
            navgusd= navgusd+1
            avgsig(navgusd)=total
            iavst(navgusd)=iavg
            tmp=tmp+total
            goto 115
          else
            keeps= 1
            if ((total.ge.totalmax).or.
     .           (total.le.totalmin)) then
              navgusd= navgusd+1
              avgsig(navgusd)=total
              iavst(navgusd)=iavg
              tmp=tmp+total
              if (navgusd.ge.maxiter) goto 115 
c     else
c     cycle
            endif
          endif
        endif     
 111  continue
c     omment  navg -> navgusd
 115  total=tmp/dble(navgusd)
      if(imode.eq.1) then
        call randa(rn1)
        rn1=rn1*tmp
        tmp=0d0
        do iavg=1,navgusd
          tmp=tmp+avgsig(iavg)
          if(tmp.ge.rn1) then
            iavg_store=iavst(iavg)
            goto 112
          endif
        enddo
      endif
 112  continue
c     
c     generate color flow, if required
      if(colgen) then
        call randa(wgcol) 
        flg90=.false.
        call matrix(alphafl,instate,alphap,hel,evtme0,1
     >       ,wgcol,colfl)
        if(flg90) then
          call selcol(alphafl,instate,rn1,rn2,iflag)
          do j1=1,npart
            call boost(1,msqref,qref,qout,alphap(1,j1))
            do j2=1,4
              alphap(j2,j1)=qout(j2)
            enddo
          enddo
          flg90=.false.
          call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     +         ,colfl)
        endif
        call setcol(npart,ip,ipinv,colfl,ifl,icu)
        call colchk(npart,icu,ifl,nev,total)
      endif
      gtotal=gtotal+total                     
      g2total=g2total+total**2
      gsum(jproc) =gsum(jproc)+total
      g2sum(jproc)=g2sum(jproc)+total**2
      if(imode.le.1) maxwgt=max(maxwgt,total)
      if(evgen) then
        if(imode.eq.0) then 
          call aleana(jproc,total)
        elseif(imode.eq.1) then 

c     write events after a pre-unweighting, using as maximum weight 
c     refwgt=1d-2*maxwgt 
          refwgt=1d-2*maxwgt
          if(unwgt(total,refwgt)) then
            write(niowgt,15) jseed,iavg_store,real(max(total,refwgt))
     $           ,real(x1)
            nevwgt=nevwgt+1d0
          endif
          call aleana(jproc,total)
          call alhbkk(total)
        elseif(imode.eq.2) then ! dump event for later Herwig processing
c     First cluster and evaluate reweighting factor
          if(ickkw.eq.1) then
            rewgt=1d0
            call alpclu(rewgt)
            avgrew=avgrew+rewgt
            nrewgt=nrewgt+1d0
            call mfill(190,sngl(rewgt),1e0)
c     next calls for debugging only
c            if(rewgt.gt.1d0) then
c              write(2,*) ' '
c              write(2,*) 'nev=',nev,' rewgt:',rewgt
c              write(2,*) (ifl(i),i=1,npart)
c              write(2,*) (icu(1,i),i=1,npart)
c              write(2,*) (icu(2,i),i=1,npart)
c              write(2,*) (sqrt(clkt(2,i)),i=1,ncltot)
c            endif
c            write(2,*)
c     $         '   fl col1  col2  cmot cdau1 cdau2  p(3)   ktin  ktout' 
c            write(2,*) 'BEFORE CLUSTERING'
c            do i=1,nclext
c             write(2,"(6(I5,1x),2(f10.2,1x))") icfl(i),iccol(1,i)
c     $             ,iccol(2,i),icmot(i),icdau(1,i),icdau(2,i)
c     $             ,sqrt(clkt(1,i)),sqrt(clkt(2,i))
c            enddo
c            write(2,*) ' '
c            write(2,*) 'CLUSTERS'
c            do i=nclext+1,ncltot
c              write(2,"(6(I5,1x),2(f10.2,1x))") icfl(i),iccol(1,i)
c     $             ,iccol(2,i),icmot(i),icdau(1,i),icdau(2,i)
c     $             ,sqrt(clkt(1,i)),sqrt(clkt(2,i))
c            enddo
c            write(2,*) ' '

c     end debug
c     3: unweight w.r.t. to reweighting factor:
            if(.not.unwgt(rewgt,1d0)) then
              unwev=unwev-1
              goto 100
            endif
          endif
c     First add to the output event record decay particles, when present
          decbr=1d0
          call setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
          call aleana(jproc,avgwgt_in*decbr)
          call evdump(niounw,nwrt,jproc,iflwrt,icuwrt,pwrt,qsq
     $         ,1d0)
        endif
      endif
c     weight bookkeeping and various optimizations

      if (evopt) then
        call wgtopt(total)
      endif

c     end generation of n-th event
c     Produce monitor file, to report partial results during run
      if(imode.ne.2.and.mod(n,nmon).eq.0) then
        if (.not.evgen) then
          tmpfac= tmpfac*ncount1/(extra*nmon)
          ncount1= 0.d0
        endif
        call alustc(fname_mpirank,'.mon',tmpstr)
        call alugun(iunit)
#ifndef USE_MPI
        open(unit=iunit,file=tmpstr,status='unknown')
        write (iunit,51) nev
 51     format(3x,' processed=',f15.1,' events') 
        write (iunit,*) '       '
#endif
        totwgt= gtotal
        effps=nev/(nev+ipscount)
        effev=nev
        sigerr=sqrt(g2total-gtotal**2/(effev))/(effev)
        avgwgt=totwgt/effev
#ifndef USE_MPI
        write(iunit,*) 'average ph-space eff=',effps
        write(iunit,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr
     $       ,' maxwgt=',maxwgt
        if(maxwgt.gt.0) write(iunit,*) 'unwgt eff=',totwgt/nev
     $       /maxwgt
        write(iunit,*) '                 '
        write(iunit,*) ' sub-processes:  '
        write(iunit,*) '                 '
#endif
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
#ifndef USE_MPI
          write(iunit,*) 'jproc=',i,' total: ',gsum(i)/(effev),'+-'
     $         ,gg2sum(i) 
#endif
        enddo
c     write(iunit,*) '               '
c     write(iunit,*) 'n. of ME calls=',ncount,avgwgt,tmpfac,totalmax
c     write(iunit,*) '               '
#ifndef USE_MPI
        close (iunit)
        call monitor(n,tmpstr)
#endif
      endif
 100  continue
      nev=nev+1d0
#ifdef USE_MPI
      if(mod(nev,1d4).eq.0) then
         flush(niowgt)
      endif
#endif
      if(nev.le.maxev) goto 99
c     
c     event generation completed
c     
 1000 continue
#ifdef USE_MPI
      call cpu_time(timestamp)
      write(6,*) 'RANK ',rankstring,' END EVTGEN AT ',timestamp
#endif
c     
c     print out final statistics
      totwgt=gtotal
      if(imode.ne.2) then
        effps=maxev/(maxev+ipscount)
        effev=maxev
        sigerr=sqrt(g2total-gtotal**2/(effev))/(effev)
        avgwgt=totwgt/effev
        if (avgwgt.eq.0.d0) then
          write(6,*) 'cross section = 0, stop'
          stop
        endif
        write(6,*) 'average ph-space eff=',effps
        write(6,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr,' maxwgt=',maxwgt
        write(6,*) 'unwgt eff =',totwgt/maxev/maxwgt
        write(6,*) '                 '
        write(6,*) ' sub-processes:  '
        write(6,*) '                 '
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
          write(6,*) 'jproc=',i,' total(pb): ',gsum(i)/(effev),'+-'
     $         ,gg2sum(i) 
        enddo
c     
c     same info, written to file
        write(niosta,*) 'average ph-space eff=',effps
        write(niosta,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr,' maxwgt='
     $       ,maxwgt
        write(niosta,*) 'unwgt eff =',totwgt/maxev/maxwgt
        write(niosta,*) '                 '
        write(niosta,*) ' sub-processes:  '
        write(niosta,*) '                 '
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
          write(niosta,*) 'jproc=',i,' total(pb): ',gsum(i)/(effev)
     $         ,'+-',gg2sum(i) 
        enddo
        if(cumerr.eq.0d0) then
          cumerr=1d0/sigerr**2
        else
          cumerr=1d0/cumerr**2+1d0/sigerr**2
        endif
        cumsig=cumsig+avgwgt/sigerr**2
        cumavg=cumsig/cumerr
        cumerr=1d0/sqrt(cumerr)
        write(niosta,*) '                 '
        write(niosta,*) 'cumulated cross-section:'
        write(niosta,*) 'avgwgt(pb)=',cumavg,'+-',cumerr
        write(6,*) '                 '
        write(6,*) 'cumulated cross-section:'
        write(6,*) 'avgwgt(pb)=',cumavg,'+-',cumerr
      elseif(imode.eq.2) then
        if(ickkw.eq.1) then 
c rescale by average reweighting factor with CKKW-like scales
          avgrew=avgrew/nrewgt
        else
          avgrew=1d0
        endif
        sigerr=g2total_in*decbr*avgrew
        avgwgt=avgwgt_in*decbr*avgrew
#ifdef USE_MPI
        call MPI_Reduce(sigerr,sigerr_sum, 1,
     $       MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(ierr.ne.0) then
           write(6,*) 'Error computing aggregate sigerr:', ierr
           call MPI_Abort(MPI_COMM_WORLD,ierr)
        endif
        call MPI_Reduce(avgwgt,avgwgt_sum, 1,
     $       MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(ierr.ne.0) then
           write(6,*) 'Error computing aggregate avgwgt:', ierr
           call MPI_Abort(MPI_COMM_WORLD,ierr)
        endif
        call MPI_Reduce(unwev,unwev_sum, 1,
     $       MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(ierr.ne.0) then
           write(6,*) 'Error computing aggregate unwev:', ierr
           call MPI_Abort(MPI_COMM_WORLD,ierr)
        endif

        if(mpirank.eq.0) then
          sigerr_global=sigerr_sum/mpiworldsize
          avgwgt_global = avgwgt_sum/mpiworldsize

          write(niosta,9995) avgwgt_global,sigerr_global, 
     $     ' ! Crosssection +- error (pb)'
          write(niosta,*) unwev_sum,unwev_sum/avgwgt_global,
     $     ' ! unwtd events, lum (pb-1)'
          write(6,*)'Crosssection(pb)=',avgwgt_global,'+-',sigerr_global
          write(6,*) 'Generated ',unwev_sum,' unweighted events, lum=',
     $       unwev_sum/avgwgt_global,'pb-1' 
        endif
#else
        write(niosta,9995) avgwgt,sigerr,' ! Crosssection +- error (pb)'
        write(niosta,*) unwev,unwev/avgwgt,' ! unwtd events, lum (pb-1)'
        write(6,*)'Crosssection(pb)=',avgwgt,'+-',sigerr
        write(6,*) 'Generated ',unwev,' unweighted events, lum=',unwev
     $       /avgwgt,'pb-1' 
#endif
      endif
 9995 FORMAT(2F25.8,A30)     
c     write out updated grid, if igridw properly initialised
      if(imode.ne.2) then
        if (evgen) tmpfac= tmpfac*ncount1/(extra*maxev)
        totalmax= avgwgt*tmpfac
        if(evopt) then
          if(evgen) then
            call alustc(fname,'.grid2',tmpstr)
          else
            call alustc(fname,'.grid1',tmpstr)
          endif
          call alugun(iunit)
#ifdef USE_MPI
          if(mpirank.eq.0) then
            open(unit=iunit,file=tmpstr,status='unknown',
     $            iostat=filestat)
            if(filestat.ne.0) then
               write(6,*) 'Error reading file in alegen subroutine:',
     $                     tmpstr
            endif
            write(6,*) 'rank ',mpirank,' dumping grid file: ',tmpstr
            call dumpgrid(iunit)
            if(iunit.ne.0) close (iunit)
          else
            iunit=0
          endif
#else
          open(unit=iunit,file=tmpstr,status='unknown')
          call dumpgrid(iunit)
          close (iunit)
#endif
        endif
      endif
      if(imode.eq.1.and.evgen) then
#ifdef USE_MPI
        write(niopar,*,iostat=filestat)
     $'number wgted evts in the file,  sigma(pb), error(pb),  maxwgt'
         if(filestat /= 0) then
            write(6,*) 'Error writing to file in alegen subroutine'
         endif
         write(niopar,10,iostat=filestat) nevwgt,cumavg,cumerr,maxwgt
         if(filestat /= 0) then
            write(6,*) 'Error writing to file in alegen subroutine'
         endif
         write(niopar,*,iostat=filestat)
     $'number weights'
         if(filestat /= 0) then
            write(6,*) 'Error writing to file in alegen subroutine'
         endif
         write(niopar,'(e12.4)',iostat=filestat) real(maxev)
         if(filestat /= 0) then
            write(6,*) 'Error writing to file in alegen subroutine'
         endif
#else
        write(niopar,*)
     $'number wgted evts in the file,  sigma(pb), error(pb),  maxwgt'
         write(niopar,10) nevwgt,cumavg,cumerr,maxwgt
         write(niopar,*)
     $'number weights'
         write(niopar,'(e12.4)') real(maxev)
#endif
      endif
c next call for debugging only
c      if(imode.ne.2) call prfinal(jprocmax,effev)
c
      end

      subroutine colchk(npart,icu,ifl,nev,w)
      implicit none
c arguments
      integer npart,icu(2,10),ifl(10)
      double precision nev,w
c locals
      integer i,j,icount(0:10),itmp
c
      do i=1,10
        icount(i)=0
      enddo
      do i=1,2
        do j=1,npart
          itmp=icu(i,j)
          icount(itmp)=icount(itmp)+1
        enddo
      enddo
      itmp=0
      do i=1,10
        if(icount(i).ne.2.and.icount(i).ne.0) itmp=1
      enddo
      if(itmp.ne.0) then
        write(6,*) 'Problem with colour flow, event=',nev,' ME=',w
        write(6,*) (ifl(j),j=1,npart)
        write(6,*) (icu(1,j),j=1,npart)
        write(6,*) (icu(2,j),j=1,npart)
      endif
      end

c-------------------------------------------------------------------
      subroutine setcol(npart,ip,ipinv,colfl,ifl,icu)
c     input colour flow from alpha (colfl) , output colour flow string
c     for herwig (icu)
c-------------------------------------------------------------
      implicit none
c inputs
      integer maxpar
      parameter(maxpar=10)
      integer npart,ip(maxpar),ipinv(maxpar),colfl(2*maxpar),icu(2
     $     ,maxpar),ifl(maxpar)
c locals
      integer nfree,nfree0,iswap,iend(maxpar),ist(maxpar)
      integer i,j,k,kb,it,nev, idbg,iq(maxpar), nq,ng
c
      data nev/0/,idbg/0/
c
c      nev=nev+1
c      if(idbg.eq.1.and.nev.gt.30) stop
      nq=0
      ng=0
      do i=1,maxpar
         icu(1,i)=0
         icu(2,i)=0
         iend(i)=1
         iq(i)=0
         if(i.le.2) then
c initial state
            ist(i)=1
         else
c final state
            ist(i)=-1
         endif
      enddo
      do i=1,npart
         icu(1,i)=colfl(2*ipinv(i)-1)
         icu(2,i)=colfl(2*ipinv(i))
      enddo
      nfree=npart
      do i=1,npart
        if(ifl(i).lt.0.and.ifl(i).gt.-7) then         
c     antiquark
          iq(i)=-1
          nq=nq+1
          iend(i)=0
          nfree=nfree-1
          if(icu(2,i).eq.0) then
            icu(2,i)=icu(1,i)
            icu(1,i)=0
          endif
        elseif(ifl(i).gt.0.and.ifl(i).lt.7) then         
c     quark
          iq(i)=1
          nq=nq+1
          iend(i)=0
          nfree=nfree-1
          if(icu(1,i).eq.0) then
            icu(1,i)=icu(2,i)
            icu(2,i)=0
          endif
        elseif(ifl(i).eq.21) then
          ng=ng+1
        endif
      enddo
c
c special for gluon-only events, flip colours of initial state gluons
c      if(nq.eq.0) then
c        do 20 i=1,2
c          if(ifl(i).eq.21) then
c            k=icu(2,i)
c            icu(2,i)=icu(1,i)
c            icu(1,i)=k
c          endif
c 20     enddo
c        return
c      endif
c     

c
c flip colours for initial state gluons
      do 20 i=1,2
        if(ifl(i).eq.21) then
          k=icu(2,i)
          icu(2,i)=icu(1,i)
          icu(1,i)=k
        endif
 20   enddo
c wrap it up if no quarks
      if(nq.eq.0) return


      if(idbg.eq.1) then
         write(6,*) 'Quark status: nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
      endif
c     
c
c swap colour-anticolour for gluons whose colour line ends on antiquarks
      do 120 i=1,npart
         if(abs(iq(i)).eq.1) then
            k=icu(1,i)
            kb=icu(2,i)
c     found a quark/antiquark
            j=0
 100        j=j+1
            if(j.gt.npart) goto 110
            if(ifl(j).ne.21.or.iend(j).eq.0) goto 100
            if(k.ne.icu(1,j).and.k.ne.icu(2,j).and.kb.ne.icu(1,j).and.kb
     $           .ne.icu(2,j)) goto 100
            iend(j)=0
            nfree=nfree-1
            iswap=0
c     colour line between init-init or final-final states: swap the gluon
            if(k.eq.icu(1,j).and.ist(i)*ist(j).gt.0) iswap=1
c     colour-anticol line between init-final states: swap the gluon
            if(k.eq.icu(2,j).and.ist(i)*ist(j).lt.0) iswap=1
c     anticol line between init-init or final-final states: swap the gluon
            if(kb.eq.icu(2,j).and.ist(i)*ist(j).gt.0) iswap=1
c     colour-anticol line between init-final states: swap the gluon
            if(kb.eq.icu(1,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(iswap.eq.1) then
               k=icu(2,j)
               icu(2,j)=icu(1,j)
               icu(1,j)=k
            endif
         endif
 110     continue
 120  enddo
c
c     go through the remaining gluons
      it=0
 180  continue
      nfree0=nfree
      if(idbg.eq.1) then
         it=it+1
         write(6,*) 'Gluon status, IT=',it,' nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
      endif
      do 220 i=1,npart
         if(ifl(i).eq.21.and.iend(i).eq.0) then
            k=icu(1,i)
            kb=icu(2,i)
            j=0
 200        j=j+1
            if(j.gt.npart) goto 210
            if(ifl(j).ne.21.or.iend(j).eq.0) goto 200
            if(k.ne.icu(1,j).and.k.ne.icu(2,j).and.kb.ne.icu(1,j).and.kb
     $           .ne.icu(2,j)) goto 200
            iend(j)=0
            nfree=nfree-1
            iswap=0
            if(k.eq.icu(1,j).and.ist(i)*ist(j).gt.0) iswap=1
            if(k.eq.icu(2,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(kb.eq.icu(2,j).and.ist(i)*ist(j).gt.0) iswap=1
            if(kb.eq.icu(1,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(iswap.eq.1) then
               k=icu(2,j)
               icu(2,j)=icu(1,j)
               icu(1,j)=k
            endif
         endif
 210     continue
 220  enddo
      if(nfree.ne.nfree0) goto 180
c check for close colour clusters

c     
c
      if(idbg.eq.1) then
         write(6,*) 'Final, IT=',it,' nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
         write(6,*) ' '
      endif
      end

cc-------------------------------------------------------------------
c      subroutine prstat(n,jseed,jprocmax,jproc,evtme,xlum,pswgt
c     $     ,djproc,wgtmax)
cc     routine for the debugging of individual large-weight events
cc     Not functional to the running of the code
cc-------------------------------------------------------------------
c
c      implicit none
c      real*8  evtme,xlum,pswgt,djproc,wgtmax,total,effev
c      integer n,jprocmax,jproc,i     
c      real*8 psum(1000), p2sum(1000)
c      data psum/1000*0/,p2sum/1000*0/
c      integer jseed(2)
c      save 
c      total=evtme*xlum*pswgt
c      if(mod(n,100000).eq.0) then
c          open (unit= 31,file= 'monitor',status='unknown')
c          write (31,51) n
c   51     format('  event=',i10) 
c          close (31)
c      endif
c      if(total.gt.1d5) then
c         write(6,*) ' '
c         write(6,*) ' '
c         write(6,*) 'large weight: event=',n,' jproc=',jproc
c         write(6,*) 'total wgt=',total,' pswgt=',pswgt,' jproc jac='
c     $        ,1/djproc,'xlum=',xlum,' ME^2=',evtme*djproc
c 8       format(5(a,e12.4))
c         if(wgtmax.eq.total) write(niosta,*) 'NEW MAX WEIGHT'
c         write(niosta,*) 'large weight: event=',n,' jseeds=',jseed(1)
c     $        ,jseed(2),' jproc=',jproc
c         write(niosta,8) 'total wgt=',total,' pswgt=',pswgt,' jproc jac='
c     $        ,1/djproc,'xlum=',xlum,' ME^2=',evtme*djproc
c         write(6,*) ' '
c         call dumpwgt
c      endif
cc bookkeeping of the grandtotals
c      psum(jproc)=psum(jproc)+total
cc     bookkeeping of the phase-space weight
c      psum(jprocmax+jproc)=psum(jprocmax+jproc)+pswgt
cc bookkeeping of the jacobian of jproc selection
c      psum(2*jprocmax+jproc)=psum(2*jprocmax+jproc)+1/djproc
cc bookkeeping of the ME^2
c      psum(3*jprocmax+jproc)=psum(3*jprocmax+jproc)+evtme*djproc
c      p2sum(jproc)=p2sum(jproc)+total**2
c      p2sum(jprocmax+jproc)=p2sum(jprocmax+jproc)+pswgt**2
c      p2sum(2*jprocmax+jproc)=p2sum(2*jprocmax+jproc)+(1/djproc)**2
c      p2sum(3*jprocmax+jproc)=p2sum(3*jprocmax+jproc)+(evtme*djproc)
c     $     **2
c      return
c      entry prfinal(jprocmax,effev)
cc      write(6,*) 'bookkeeping of the grandtotals'
c      write(niosta,*) 'bookkeeping of the grandtotals'
c      do i=1,jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i,' total: ',psum(i)/(effev),'+-'
cc     $        ,p2sum(i) 
c         write(niosta,*) 'jproc=',i,' total: ',psum(i)/(effev),'+-'
c     $        ,p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the  phspace wgt'
c      write(niosta,*) 'bookkeeping of the  phspace wgt'
c      do i=jprocmax+1,2*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-jprocmax,' pswgt: ',psum(i)/(effev),'+-'
cc     $        ,p2sum(i) 
c         write(niosta,*) 'jproc=',i-jprocmax,' pswgt: ',psum(i)/(effev),'+-'
c     $        ,p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the jproc jacobian'
c      write(niosta,*) 'bookkeeping of the jproc jacobian'
c      do i=2*jprocmax+1,3*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-2*jprocmax,' prcwgt: ',psum(i)/(effev)
cc     $        ,'+-',p2sum(i) 
c         write(niosta,*) 'jproc=',i-2*jprocmax,' prcwgt: ',psum(i)/(effev)
c     $        ,'+-',p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the ME^2'
c      write(niosta,*) 'bookkeeping of the ME^2'
c      do i=3*jprocmax+1,4*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-3*jprocmax,' ME^2: ',psum(i)/(effev)
cc     $        ,'+-',p2sum(i) 
c         write(niosta,*) 'jproc=',i-3*jprocmax,' ME^2: ',psum(i)/(effev)
c     $        ,'+-',p2sum(i) 
c      enddo
c      do i=1,100 
c         psum(i)=0d0
c         p2sum(i)=0d0
c      enddo
c      end
c
      subroutine wgtopt(total)
      implicit none
c inputs
      real *8 total
c common
      integer nvar,nv
      common/psopt/nvar,nv
      integer maxn
      parameter (maxn= 100)  
      integer mask,mmask
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c locals
      real *8 dummy,wgtt
      integer n,ndummy,nn
      data ndummy/-1/
      wgtt= total

      call onedimbin(2,mmask(1),wgtt,1,ndummy,dummy)
      nn=jgrid(jproc)
      do n= (nn-1)*nv+2,nn*nv+1
        call onedimbin(2,mmask(n),wgtt,n,ndummy,dummy)
      enddo
      end

c-------------------------------------------------------------------
      subroutine evdump(Nunit,npart,iproc,ifl,icolup,p,qsq,wgt)
c     dump event details to a file, for future reading by herwig
c-------------------------------------------------------------------
      implicit none
      integer maxdec
      parameter (maxdec=40)
c inputs
      integer Nunit,npart,iproc,ifl(maxdec),icolup(2,maxdec)
      real*8 qsq,p(5,maxdec),wgt
c locals
      real Sq,Sp(3,maxdec),Sm(maxdec),Swgt
      integer nev,i,j
      data nev/0/
      save
c
      if(npart.gt.maxdec) then
         write(6,*) 'Sp(4,maxdec) in evdump underdimensioned'
         stop
      endif
      nev=nev+1
      Sq=real(sqrt(qsq))
      Swgt=real(wgt)
      do i=1,npart
         do j=1,3
            Sp(j,i)=real(p(j,i))
         enddo
         Sm(i)=real(p(5,i))
      enddo
c     Nevent, iproc, npart, wgt, Q scale 
      write(Nunit,2) nev,iproc,npart,Swgt,Sq
 2    format(i8,1x,i4,1x,i2,2(1x,e12.6))
c     flavour, colour and z-momentum of incoming partons
      write(Nunit,8) ifl(1),icolup(1,1),icolup(2,1),Sp(3,1)
      write(Nunit,8) ifl(2),icolup(1,2),icolup(2,2),Sp(3,2)
c     flavour, colour, 3-momentum and mass of outgoing partons
      do i=3,npart
      write(Nunit,9) ifl(i),icolup(1,i),icolup(2,i),Sp(1,i),Sp(2,i)
     $     ,Sp(3,i),Sm(i)
      enddo
 8    format(i8,1x,2(i4,1x),f10.3)
 9    format(i8,1x,2(i4,1x),4(1x,f10.3))
      end

c-------------------------------------------------------------------
      subroutine setalp(npart,flin,pin,flout,pout,instate,ip,ipinv)
c
c     First reorder particles according to alpha ordering
c     and then reorder the respective momenta, etc.
c     Alpha ordering:
c     i= 1    2    3     4  5    6    7    8     9 10 11  12 13 14 15 16
c        ubar dbar nubar e+ cbar sbar tbar bbar  u d  nu  e- c  s  t  b  
c        17     18   19  20  21     22
c        Higgs  Z0   W-  W+  photon gluon
c     (all particles outgoing)
c     ip(i): position of the "i-th" particle of the Alpha array,
c            in the array of momenta defined by routine alegen
c     ipinv(i): inverse of ip, i.e. where to find the "i"-th particle
c            of the alegen array in the alpha momentum array 
c
c     flin is the sequence of flavours required by alpha to fully
c     define a specific process, with the generation ordering.
c
c     flout is the sequence of flavours ordered according to the alpha
c     ordering: flout(i)=flin(ip(i))
c     
c     pout is the sequence of momenta, ordered according to the alpha
c     ordering. All energies are positive, and the pointers to the two
c     momenta corresponding to the intial state partons are instate(i),
c     i=1,2
c     
c-------------------------------------------------------------------
      implicit none
      integer i,j,k,ix,jx,kp,iff,npart,maxpar
      parameter (maxpar=10)
      integer flin(maxpar),flout(maxpar),ip(maxpar),ipinv(maxpar)
      integer instate(2)
      real *8 pin(5,maxpar),pout(4,maxpar)
c     Alpha ordering of outgoing momenta
      integer aref,aord(22)
c     Labels of the particles according to the PDG
c
      data aord/-2,-1,-12,-11,-4,-3,-6,-5,2,1,12,11,4,3,6,5,
     .          25,23,-24,24,22,0/
c
      do i=1,maxpar
       flout(i)=1001
      enddo
c
      if(npart.gt.maxpar) then
         write(6,*) 'error in setper:'
         write(6,*) 'number of particles larger than ',maxpar
         stop
      endif
c     First reorder particles according to alpha ordering
c
c     scan particle types in the Alpha order
      k=1
      do i=1,22
         aref=aord(i)
c     search for all particles of flavour aref in the alegen array, and
c     order them in the "ip" array
         do j=1,npart
            iff=flin(j)
c     turn incoming to outgoing flavours
            if(j.le.2.and.abs(iff).lt.20) iff=-iff
c
            if(iff.eq.aref) then
               ip(k)=j
               ipinv(j)=k
               flout(k)=flin(j)
               if(j.le.2) instate(j)=k
               do ix=1,4
c     alegen momentum conventions:     p(1)=Px p(2)=Py  p(3)=Pz  p(4)=E   
c     alpha momentum conventions:      p(1)=E  p(2)=Px  p(3)=Py  p(4)=Pz
                  jx=ix+1
                  if(jx.eq.5) jx=1
                  pout(jx,k)=pin(ix,j)
               enddo
               if(k.eq.npart) goto 99
               k=k+1
            endif
         enddo
      enddo
 99   continue
c reassign momentum labels to agree with alpha conventions
c                             
      do k=1,npart
        kp=ip(k)
        do ix=1,4
c     alegen momentum conventions:
c     p(1)=Px   p(2)=Py    p(3)=Pz p(4)=E   
c     alpha momentum conventions:
c     p(1)=E  p(2)=Px   p(3)=Py    p(4)=Pz
           jx=ix+1
          if(jx.eq.5) jx=1
          pout(jx,k)=pin(ix,kp)
        enddo
      enddo
c
c
      do k=1,npart
         if(abs(flout(k)).eq.24) flout(k)= -flout(k)
      enddo
      end

c************
c************
c-------------------------------------------------------------------
      subroutine setpdf
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      real *8 alfas
      real xx1,xx2,q2,tmp
c
      q2=real(qsq)                   
      xx1=real(x1)
      xx2=real(x2)
      call mlmpdf(nmnr,ih1,q2,xx1,f1,5)
      call mlmpdf(nmnr,ih2,q2,xx2,f2,5)      
c     adopt pdg numbering scheme for up and down quarks
      tmp=f1(1)
      f1(1)=f1(2)
      f1(2)=tmp
      tmp=f2(1)
      f2(1)=f2(2)
      f2(2)=tmp
      tmp=f1(-1)
      f1(-1)=f1(-2)
      f1(-2)=tmp
      tmp=f2(-1)
      f2(-1)=f2(-2)
      f2(-2)=tmp
c
      as=alfas(qsq,xlam,nloop,-1)
c     multiply alphas by a rescaling factor, to allow separation of qcd and ew
c     contributions to jet production: 
c
c                 sig=1/(resc**N) * sig(as->as*resc)
c
c     for processes where maximum power of alphas is N.
c     resc>>1 will suppress ew processes
c     resc<<1 will suppress qcd processes
      apar(56)=sqrt(4*pi*as*resc)
      end                                    

      subroutine grans(n1,nx1,x1)      
c----------------------------------------------------------------c
c                                                                c       
c     x1 is uniformely genetared in                              c
c     x1a < x1 < x1b.                                            c
c                                                                c
c     n1 is the bin, nx1 is the total number of bins             c
c                                                                c       
c----------------------------------------------------------------c

      implicit none
      real*8 x1,ran1,x1a,x1b
      integer n1,nx1   
      call rans(ran1)

      x1a= dfloat(n1-1)/dfloat(nx1)
      x1b= x1a+1.d0/dfloat(nx1)
      x1= x1b*ran1+(1.d0-ran1)*x1a

      return
      end

      subroutine photsm(lflag,cn,cxm,cxp,s,dj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Massless particle propagator:                                c
c                                                                c
c   the invariant mass squared has the distribution              c
c                                                                c
c                       1/s^cn                                   c
c                       cxm < s < cxp                            c
c                                                                c
c              INPUT                           OUTPUT            c
c                                                                c
c  lflag= 0:   cn, cxm, cxp                    s, dj             c
c                                                                c
c  lflag= 1:   cn, cxm, cxp , s                dj                c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 cn,cxm,cxp,s,dj,ran1,tjm,hj
      integer lflag
      dj= 0.d0
      if (lflag.eq.0) then
        s= tjm(0.d0,cn,cxm,cxp,1,ran1)
      else
        if (s.le.cxm.or.s.ge.cxp) return
      endif
      dj= 1.d0/(s**cn*hj(0.d0,cn,cxm,cxp))
      return
      end

      subroutine resonm(lflag,rm,ga,lim,cxm,cxp,s,dj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Resonant propagator:                                         c
c                                                                c
c   the invariant mass squared has the distribution              c
c                                                                c
c                    1/[(s-rm^2)^2+rm^2*ga^2]                    c
c                                                                c
c   if lim= 0 -inf < s < +inf                                    c
c   if lim= 1  cxm < s < cxp                                     c
c                                                                c
c              INPUT                           OUTPUT            c
c                                                                c
c  lflag= 0:   rm, ga, lim, (cxm, cxp)         s, dj             c
c                                                                c
c  lflag= 1:   rm, ga, lim, s, (cxm, cxp)      dj                c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 rm,ga,cxm,cxp,s,dj,ran1,pi,rms,ard,arp,arm,ym,ypmym
      integer lflag,lim    
      parameter(pi= 3.14159265358979324d0)
      dj= 0.d0
      rms= rm*rm
      if (lim.eq.0) then
        if (lflag.eq.0) then
          s= rms+rm*ga*tan(pi*(ran1-0.5d0))
        endif
        dj= 1.d0/(pi/rm/ga*((s-rms)*(s-rms)+rms*ga*ga))
      else
        if (lflag.eq.1) then
          if (s.lt.cxm.or.s.gt.cxp) return
        endif
        ard= cxp-cxm
        arp= (cxp-rms)/rm/ga
        arm= (cxm-rms)/rm/ga
        ym= atan(arm)
        ypmym= atan(ard/rm/ga/(1.d0+arp*arm))
        if (arp*arm.lt.-1.d0) then
          if (arp.gt.0.d0) ypmym= ypmym+pi
          if (arp.lt.0.d0) ypmym= ypmym-pi
        endif
        if (lflag.eq.0) then
          s= rms+rm*ga*tan(ran1*(ypmym)+ym)
        endif
        dj= 1.d0/(ypmym/rm/ga*((s-rms)*(s-rms)+rms*ga*ga))
      endif
      return
      end

      subroutine dec2fm(lflag,s,p,s1,s2,p1,p2,dj,ran)
c----------------------------------------------------------------c
c                                                                c
c   Isotropic 2-body decay:                                      c
c                                                                c
c                       p ----> p1 + p2                          c
c                       s1= p1^2                                 c
c                       s2= p2^2                                 c
c                                                                c
c              INPUT                          OUTPUT             c
c                                                                c
c  lflag= 0:   s, p, s1, s2                   p1, p2, dj         c
c                                                                c
c  lflag= 1:   s, s1, s2                      dj                 c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 pi,s,s1,s2,dj,sqlam,rs,p1m,ran1,ran2,ct,st,phi,cp,sp
      integer lflag,k    
      parameter (pi= 3.14159265358979324d0)
      real*8 p(0:3),p1(0:3),p2(0:3),p1h(0:3)
      real*8 ran(1:2)      
      dj= 2.d0/pi/sqlam(s,s1,s2)
      if (lflag.eq.0) then
        rs= sqrt(dabs(s))
        p1h(0)= (s+s1-s2)/rs/2.d0
        p1m= rs/dj/pi
        ran1= ran(1)
        ran2= ran(2)
        ct= 2.d0*ran1-1.d0
        st= sqrt(1.d0-ct*ct)
        phi= 2.d0*pi*ran2
        cp= cos(phi)
        sp= sin(phi)
        call qvec(p1h(0),p1m,st,ct,sp,cp,p1h)
        call boost(0,s,p,p1h,p1)
        do k= 0,3
           p2(k)= p(k)-p1(k)
        enddo
      endif
      return
      end

      function tjm(a,cn,cxm,cxp,k,ran1)
c----------------------------------------------------------------c
c                                                                c
c   tjm is produced according to the distribution                c
c                                                                c
c                    (1/tjm)^cn     where:                       c
c                                                                c
c                    tjm= a+k*ct                                 c  
c                    cxm < ct < cxp                              c
c                    k= +1 or -1                                 c  
c                                                                c  
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,ran1,ce,tjm
      integer k 
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        if (k.eq.1) then
           tjm= (ran1*(a+cxp)**ce+(1.d0-ran1)*
     *         (a+cxm)**ce)**(1.d0/ce)
        else
           tjm= (ran1*(a-cxm)**ce+(1.d0-ran1)*
     *         (a-cxp)**ce)**(1.d0/ce)
        endif
      else
        if (k.eq.1) then
          if(cxp.lt.-a) then
            tjm= -exp(ran1*log(-a-cxp)+
     *          (1.d0-ran1)*log(-a-cxm))
          else
            tjm= exp(ran1*log(a+cxp)+
     *          (1.d0-ran1)*log(a+cxm))
          endif
        else
          if(cxp.lt.a) then
            tjm= exp(ran1*log(a-cxm)+
     *          (1.d0-ran1)*log(a-cxp))
          else
            tjm= -exp(ran1*log(-a+cxm)+
     *          (1.d0-ran1)*log(-a+cxp))
          endif
        endif
      endif
      return
      end

      function rantjm(a,cn,cxm,cxp,k,tjm)
c----------------------------------------------------------------c
c                                                                c
c     given tjm produced according to the distribution           c
c                                                                c
c                    (1/tjm)^cn     where:                       c
c                                                                c
c                    tjm= a+k*ct                                 c  
c                    cxm < ct < cxp                              c
c                    k= +1 or -1                                 c  
c                                                                c  
c     rantjm is the corresp. uniform random number in [0,1]      c 
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,rantjm,ce,tjm
      integer k 
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        if (k.eq.1) then
          rantjm= (tjm**ce-(a+cxm)**ce)/((a+cxp)**ce-(a+cxm)**ce)
        else
          rantjm= (tjm**ce-(a-cxp)**ce)/((a-cxm)**ce-(a-cxp)**ce)
        endif
      else
        if (k.eq.1) then
          if(cxp.lt.-a) then
            rantjm= (log(-tjm)-log(-a-cxm))/
     *              (log(-a-cxp)-log(-a-cxm))
          else
            rantjm= (log(tjm)-log(a+cxm))/
     *              (log(a+cxp)-log(a+cxm))
          endif
        else
          if(cxp.lt.a) then
            rantjm= (log(tjm)-log(a-cxp))/
     *              (log(a-cxm)-log(a-cxp))
          else
            rantjm= (log(-tjm)-log(-a+cxp))/
     *              (log(-a+cxm)-log(-a+cxp))
          endif
        endif
      endif
      return
      end

      subroutine boost(lflag,sq,q,ph,p)
c----------------------------------------------------------------c
c                                        _                       c
c   Boost of a 4-vector ( relative speed q/q(0) ):               c
c                                                                c
c   ph is the 4-vector in the rest frame of q                    c
c   p is the corresponding 4-vector in the lab frame             c
c                                                                c  
c              INPUT                               OUTPUT        c
c                                                                c
c  lflag= 0:  sq, q, ph                            p             c
c                                                                c
c  lflag= 1:  sq, q, p                             ph            c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 sq,c1,rsq
      integer lflag,i   
      real*8 q(0:3),ph(0:3),p(0:3)
      rsq= sqrt(sq)
      if (lflag.eq.0) then
        p(0)= (q(0)*ph(0)+q(1)*ph(1)+q(2)*ph(2)+q(3)*ph(3))/rsq
        c1= (ph(0)+p(0))/(rsq+q(0))
        do i= 1,3
          p(i)= ph(i)+q(i)*c1
        enddo
      else
        ph(0)= (q(0)*p(0)-q(1)*p(1)-q(2)*p(2)-q(3)*p(3))/rsq
        c1= (p(0)+ph(0))/(rsq+q(0))
        do i= 1,3
          ph(i)= p(i)-q(i)*c1
        enddo
      endif
      return
      end

      subroutine qvec(q0,qm,st,ct,sp,cp,q)
c----------------------------------------------------------------c
c                                                                c
c   A 4-vector is built                                          c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 q0,qm,st,ct,sp,cp
      real*8 q(0:3)
      q(0)= q0
      q(1)= qm*st*sp        
      q(2)= qm*st*cp       
      q(3)= qm*ct
      return
      end 

      function sqlam(s,s1,s2)
      implicit none
      real*8 s,s1,s2,x1,x2,arg1,sqlam,arg
      x1= s1/s
      x2= s2/s
      arg= 1.d0-x1-x2
      arg1= arg*arg-4.d0*x1*x2
      if (arg1.gt.0.d0) then
        sqlam= sqrt(arg1)
      else
        sqlam= 1.d-30
      endif
      return
      end

      function hj(a,cn,cxm,cxp)
c----------------------------------------------------------------c
c                                                                c
c   Normalization factor for tj                                  c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,ce,hj
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        hj= ((a+cxp)**ce-(a+cxm)**ce)/ce
      else
        hj= log((a+cxp)/(a+cxm))
      endif
      return
      end

      subroutine onedimbin(lflag,n1,w,n,igridw,peropt)
c-
c-    Self optimizing binning
c- 
c-    n labels the variable for which one wishes to
c-    have a self optimized binning
c-
c-    if lflag= 1:
c-        n1 is generated in [1:nx1] and w= djtot1 (the inverse jacobian)
c-        is an output.
c-
c-    if lflag= 2:
c-        w is the input weight for the optimization. The optimization
c-        is performed only if n1= 1.
c-
c-    if lflag= 3:
c-        given n1 and n w= djtot1 (the inverse jacobian) is computed
c-
c-    if lflag= 0:
c-        the new set of a-priori weights is computed and
c-        written to unit= igridw if igridw.ne.0.
c-
      implicit none
      real*8 peropt,w,beta,rcho,wt1,ni
      real*8 al1,bet1,v1,aptot1,alsum1,alsum2
      integer lflag,n1,n,ncmax,maxn,modopt,igridw
      integer nx1,nx2,nct1,nct,init,i,m,j,k,lrun
      parameter (beta= 0.5d0, ncmax= 1000, maxn= 100)
c-
c-    if modopt= 0 [1] the grid optimization is done according to the
c-    variance [cross section]
c-
      parameter (modopt= 1)
c-
c-    peropt is the percentage [0:1] of allowed optimization
c-
c-    WARNING: ncmax must be .le. 6000
c-
      real*8 s1(ncmax)
      integer iopt(maxn),jad1(maxn)
      integer ntmp  
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      common/book/v1(ncmax),ni(maxn)
comment
      save 
*
c-    self-optimization of the a-priori weights   
*                
      if (lflag.eq.0) then
*
c-      only if there was optimization the new set is computed
*
        if (iopt(n).eq.1) then
           do i= nct(n)+1,nct(n)+nx1(n)
             s1(i)= v1(i)/ni(n)
           enddo
           aptot1= 0.d0                                             
           alsum1= 0.d0                                             
           alsum2= 0.d0                                             
           do m= 1,nx1(n)                                             
             if (modopt.eq.0) then
               aptot1= aptot1+al1(m,n)*s1(nct(n)+m)**beta
             else 
               if (s1(nct(n)+m).lt.1.d-10) then
                  aptot1= aptot1+1.d-10 
               else
                  aptot1= aptot1+s1(nct(n)+m)                 
               endif
             endif
           enddo                                                    
           do m= 1,nx1(n)                                             
             if (modopt.eq.0) then
               al1(m,n)= al1(m,n)*(s1(nct(n)+m)**beta)/aptot1       
             else
               if (s1(nct(n)+m).lt.1.d-10) then
                 al1(m,n)=  1.d-10/aptot1
               else
                 al1(m,n)=  s1(nct(n)+m)/aptot1
               endif 
             endif
           enddo                                                    
*
c-         rescaling  according to the chosen value of peropt
*
           do m= 1,nx1(n)                                             
             al1(m,n)= al1(m,n)*peropt
             al1(m,n)= al1(m,n)+(1.d0-peropt)/dfloat(nx1(n))
           enddo                                                    
        endif
*
        if (igridw.ne.0) then
*
c-        if igridw= 0 the a-priori weights are not written 
c-        if igridw.ne.0 they are written in unit= igridw
*
          write (igridw,20) n,(m,al1(m,n),m= 1,nx1(n))
        endif
*
c-      computes the new bet1(j) 
*
        do j= 1,nx1(n)                                             
          bet1(j,n)= 0                                             
          do k= 1,j                                              
            bet1(j,n)= bet1(j,n)+al1(k,n)                              
          enddo
        enddo
        return 
*
c-    generation of n1 and computation of w= djtot1
*
      elseif (lflag.eq.1) then
         call rans(rcho)                                            
         ntmp = nx1(n) 
         do lrun= 1,ntmp                                            
            if (rcho.le.bet1(lrun,n)) then
               jad1(n)= lrun
               goto 10                        
            endif
         enddo                                                      
 10      n1= jad1(n)
         w = ntmp*al1(n1,n)                             
         return
*
c-    bookkeeping for the self-optimization, w is an input 
*
      elseif(lflag.eq.2) then  
         if (n1.eq.1) then
            iopt(n)= 1
            if (modopt.eq.0) then
              wt1= w*w/al1(jad1(n),n)
            else
              wt1= w
            endif
            ntmp       = nct(n)+jad1(n)
            v1(ntmp)   = v1(ntmp)+wt1
            ni(n)      = ni(n)+1
         endif
*
c-    given n1 and n w= djtot1 (the inverse jacobian) is computed
*
      elseif(lflag.eq.3) then  
         w = nx1(n)*al1(n1,n) 
      endif
 20   format(/,'  a-priori weights for variable n.',i4,
     +       //,(i3,1x,d20.9))
 21   format(///,(4x,d20.9))
      return
      end
*
      subroutine grid1W(nwrt,n)
      implicit none
      real*8 al1
      integer nwrt,n,ncmax,maxn,filestat
      integer nct,nx1,nx2,nct1,m     	
      parameter (ncmax= 1000, maxn= 100)
c-
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      save 

#ifdef USE_MPI
      write (nwrt,20,iostat=filestat) n,(m,al1(m,n),m= 1,nx1(n))
      if(filestat /= 0) then
         write(6,*) 'Error writing to file in grid1W subroutine'
      endif
#else
      write (nwrt,20) n,(m,al1(m,n),m= 1,nx1(n))
#endif
 20   format(/,'  a-priori weights for variable n.',i4,
     +       //,(i3,1x,d20.9))
      end

      subroutine grid1R(nrd,n)
      implicit none 
      real*8 al1,bet1
      integer nct,nx1,nx2,nct1,init,j,k,nrd,i,n
      integer ncmax,maxn,filestat
      parameter (ncmax= 1000, maxn= 100)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      save 

c-    if nrd.ne.0
c-    the a-priori weights are read from unit= nrd
#ifdef USE_MPI
      if (nrd.ne.0) then
         read (nrd,21,iostat=filestat) (al1(i,n),i= 1,nx1(n))
         if(filestat /= 0) then
            write(6,*) 'Error reading file in grid1R subroutine'
         endif
      endif
#else
      if (nrd.ne.0) read (nrd,21) (al1(i,n),i= 1,nx1(n))
#endif
         
      do j= 1,nx1(n)                                         
        bet1(j,n)= 0                                         
        do k= 1,j                                          
          bet1(j,n)= bet1(j,n)+al1(k,n)                          
        enddo                                  
      enddo
cmlm
      if(abs(bet1(nx1(n),n)-1d0).lt.1d-6) then
         bet1(nx1(n),n)=1d0
      else
         write(6,*) 'grid entries for variable',n,' don''t sum to 1:'
         write(6,*) 'beta1=',bet1(nx1(n),n)
         write(6,*) 'error initialising grid, stop'
         stop
      endif

 21   format(///,(4x,d20.9))
      end

c-------------------------------------------------------------------
      subroutine dumpgrid(igridw)
c-------------------------------------------------------------------
      implicit none
c input variables
      integer igridw
      double precision peropt,totalmax,tmpfac
c commons
      integer nvar,nv
      common/psopt/nvar,nv
      integer maxn
      parameter (maxn= 100)  
      integer mask,mmask 
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/mxfact/totalmax,tmpfac
c locals
      integer ndummy,n
      real*8 dummy
c
      do n= 1,nvar+1
        call onedimbin(0,ndummy,dummy,n,igridw,peropt(n))
      enddo
#ifdef USE_MPI
      if(igridw.ne.0) write (igridw,20) totalmax,tmpfac
#else
      write (igridw,20) totalmax,tmpfac
#endif
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
      end

c-------------------------------------------------------------------
      subroutine readgrid(nrd)
c-------------------------------------------------------------------
      implicit none
c inputs
      integer nrd,filestat
c common
      integer nvar,nv
      common/psopt/nvar,nv
      real*8 totalmax,tmpfac
      common/mxfact/totalmax,tmpfac
c locals
      integer  n
c
c      if (nrd.ne.niopar) open (unit= nrd,status='old')
      do n= 1,nvar+1
        call grid1r(nrd,n)
      enddo
#ifdef USE_MPI
      read(nrd,20,iostat=filestat) totalmax,tmpfac
      if(filestat /= 0) then
        write(6,*) 'Error reading from file in readgrid subroutine'
      endif
#else
      read(nrd,20) totalmax,tmpfac
#endif
 20   format(///,2(d20.9))
c      if (nrd.ne.niopar) close(nrd)

      end
*
      subroutine spk(lflag,a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,
     .               etc0,etc1,eta3,wjac,ran0,lw)
c----------------------------------------------------------------c
c                                                                c
c     It returns eta3, generated within etc0 and etc1,           c
c     such that the invariant mass s23 has the distribution      c
c                                                                c
c     1/(a+s_23)^cn                                              c
c                                                                c
c     If lflag.eq.0:                                             c
c                                                                c
c     Input: a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,              c
c            etc0,etc1,ran0                                      c
c                                                                c
c     Ouput: eta3,wjac,lw                                        c
c                                                                c
c     If lflag.eq.1:                                             c
c                                                                c
c     Input: a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,eta3          c
c            etc0,etc1                                           c
c                                                                c
c     Ouput: ran0,wjac,lw                                        c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 pt(2:3),eta(2:3),etc0,etc1,xmr(2:3)
      real*8 al2,al3,xp,xm,x2p,x2m,zc0,zc1,x3p0,x3p1,x3mi,sm
      real*8 x3m,x3p,pt23s,yp,wjac,sp,spp,smm,x,al,fla,cn,wj,s
      real*8 rnn,ran0,wf,sig,sc,a
      real*8 pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,eta3
      integer lw,lflag
*
c-    Set local variables 
*      
      pt(2)  = pt2
      pt(3)  = pt3
      xmr(2) = xmr2
      xmr(3) = xmr3
      eta(2) = eta2
      if (lflag.eq.1) eta(3)= eta3
      pt23s  = pt(2)**2+pt(3)**2+2.d0*pt(2)*pt(3)*cos(phi2-phi3)
*
      al2= pt(2)**2+xmr(2)**2
      al3= pt(3)**2+xmr(3)**2
      if (eta(2).ge.0.d0) then
        xp = sqrt(pt(2)**2*(cosh(eta(2)))**2+xmr(2)**2)
     .           +pt(2)*sinh(eta(2))
        x2p= xp
        x2m= al2/x2p
      else
        xm = sqrt(pt(2)**2*(cosh(eta(2)))**2+xmr(2)**2)
     .           -pt(2)*sinh(eta(2))
        x2m= xm
        x2p= al2/x2m
      endif
      if (lflag.eq.1) then
        if (eta(3).ge.0.d0) then
          xp = sqrt(pt(3)**2*(cosh(eta(3)))**2+xmr(3)**2)
     .             +pt(3)*sinh(eta(3))
          x3p= xp
          x3m= al3/x3p
        else
          xm = sqrt(pt(3)**2*(cosh(eta(3)))**2+xmr(3)**2)
     .             -pt(3)*sinh(eta(3))
          x3m= xm
          x3p= al3/x3m
        endif
        s=-pt23s+al2+al3+x3p*x2m+x2p*x3m 
      endif
      zc0 = 2.d0*pt(3)*sinh(etc0)
      zc1 = 2.d0*pt(3)*sinh(etc1)
      if (zc0.ge.0.d0) then
        x3p0= 0.5d0*(zc0+sqrt(zc0*zc0+4.d0*al3))
      else
        x3p0=-al3/(0.5d0*(zc0-sqrt(zc0*zc0+4.d0*al3)))
      endif
      if (zc1.ge.0.d0) then
        x3p1= 0.5d0*(zc1+sqrt(zc1*zc1+4.d0*al3))
      else
        x3p1=-al3/(0.5d0*(zc1-sqrt(zc1*zc1+4.d0*al3)))
      endif
      x3mi= sqrt(al2*al3)/x2m
      sm  = -pt23s+al2+al3+x2m*x3p0+al2*al3/x2m/x3p0
      sp  = -pt23s+al2+al3+x2m*x3p1+al2*al3/x2m/x3p1
      if     (x3p0.ge.x3mi) then
        smm= sm
        spp= sp
        wf = 1.d0
        if (lflag.eq.0) then
          sig= 1.d0
          rnn= ran0
        else
          sig= 0.d0
        endif
      elseif (x3p1.le.x3mi) then
        smm= sp
        spp= sm
        wf = 1.d0
        if (lflag.eq.0) then
          sig=-1.d0
          rnn= ran0
        else
          sig= 0.d0
        endif
      elseif (x3p0.lt.x3mi.and.x3p1.gt.x3mi) then
        sc= -pt23s+al2+al3+2.d0*sqrt(al2*al3)
        wf = 2.d0
        if (lflag.eq.0) then
          if (ran0.le.0.5d0) then
            smm= sc
            spp= sm
            sig=-1.d0
            rnn= ran0*2.d0
          else
            smm= sc
            spp= sp
            sig= 1.d0
            rnn= ran0*2.d0-1.d0
          endif
        else
          if (x3p.lt.x3mi) then
            smm= sc
            spp= sm
            sig= 0.d0
          else
            smm= sc
            spp= sp
            sig= 1.d0
          endif
        endif
      else
        print*,'error'
        stop
      endif
      if (smm.le.0.d0) goto 100
      if (spp.le.0.d0) goto 100
      if (smm.ge.spp)  goto 100
      call ppeaka(lflag,a,cn,smm,spp,s,wj,rnn,lw)
      if (lw.eq.1) goto 100
      if (lflag.eq.1) then
        if (pt(3).lt.1.d-20) goto 100
        fla = dabs(x3p*x2m-x3m*x2p)
        if (fla.lt.1.d-20)   goto 100
        wjac= (x3p+x3m)/(2.d0*pt(3)*cosh(eta(3)))/
     .       fla*wf
        wjac= wjac*wj
        ran0= (rnn+sig)/wf 
        return
      endif
      al  = pt23s+s
      fla = al**2+al2**2+al3**2-2.d0*(al*al2+al*al3+al2*al3)
      fla = sqrt(dabs(fla))
      yp  = 0.5d0*((al-al2-al3)+fla)
      if (dabs(sig-1.d0).lt.1.d-5) x3p= yp/x2m
      if (dabs(sig+1.d0).lt.1.d-5) x3p= al2*al3/yp/x2m
      x3m = al3/x3p
      x   = 0.5d0/pt(3)*(x3p-al3/x3p)
      if (x.gt.0.d0) then
        eta(3)= log(sqrt(x*x+1.d0)+x)
      else
        eta(3)=-log(sqrt(x*x+1.d0)-x)
      endif
      wjac= (x3p+x3m)/(2.d0*pt(3)*cosh(eta(3)))/
     .      dabs(x3p*x2m-x3m*x2p)*wf
      wjac= wjac*wj
      eta3= eta(3)
*
      return
 100  wjac= 0.d0
      lw= 1
      end

*
      subroutine etagen(lflag,ph1,ph2,eta1,eta2m,eta2p,eta2,wjac,rnd,lw)
c---------------------------------------------------------------c
c
c     Given 2 massless 4-momenta p1 and p2, eta2 is generated 
c     in such a way that  s12 has the distribution 1/s12:
c
c             1
c        --------------
c        cosh(eta2-a)-b
c
c     When p1 is along z use a=b=0 --> eta1=0, ph1=0, ph2=pi/2.
c
c     When lflag.eq.0
c
c       input:    ph1,ph2,eta1,eta2m,eta2p,rnd
c       
c       output:   eta2,wjac
c
c     When lflag.eq.1
c
c       input:    ph1,ph2,eta1,eta2,eta2m,eta2p
c       
c       output:   wjac,rnd
c 
c---------------------------------------------------------------c

      implicit none
      real*8 ph1,ph2,eta1,eta2,eta2m,eta2p,rnd,wjac
      real*8 a,b,c,xm,xp,bb,aa,den
      integer lflag,lw
      lw= 0
*
      if (lflag.eq.1) then
        if (eta2.lt.eta2m) goto 100
        if (eta2.gt.eta2p) goto 100
      endif
*
c-    local variables:
*      
      a = eta1
      b = cos(ph1-ph2)
*
c-    Protection:
*
      if (b.ge.0.d0) then
        b= min(b,0.99d0)
      else
        b= max(b,-0.99d0)
      endif
*
      c = sqrt(dabs(1.d0-b*b))
      xm= eta2m
      xp= eta2p
*
      bb = atan((exp(xm-a)-b)/c)
      den= atan((exp(xp-a)-b)/c)-bb
      if (dabs(den).lt.1.d-12) goto 100
      aa = c/2.d0/den
      if (lflag.eq.0) then
        den= b+c*tan(0.5d0*c*rnd/aa+bb)
        if (den.le.0.d0) goto 100
        eta2= a+log(den)
      endif
      wjac= (cosh(eta2-a)-b)/aa
      if (lflag.eq.1) rnd= 2.d0*aa/c*(atan((exp(eta2-a)-b)/c)-bb)
*
      return
 100  wjac= 0.d0
      lw= 1
      return
      end

      SUBROUTINE RAMBO(LFLAG,N,ET,XM,P,DJ)
C------------------------------------------------------
C
C                       RAMBO
C
C    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
C    (MODIFIED BY R. PITTAU)
C
C                INPUT                 OUTPUT
C
C    LFLAG= 0:   N, ET, XM             P, DJ
C    LFLAG= 1:   N, ET, XM, P          DJ
C
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=100 )
C    P  = PARTICLE MOMENTA ( DIM=(4,100) )
C    DJ = 1/(WEIGHT OF THE EVENT)
C
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XM(100),P(4,100),Q(4,100),Z(100),R(4),
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(5)
      SAVE ACC,ITMAX,IBEGIN,IWARN,Z,TWOPI,PO2LOG
      DATA ACC/1.D-14/,ITMAX/10/,IBEGIN/0/,IWARN/5*0/
C
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.*DATAN(1.D0)
      PO2LOG=LOG(TWOPI/4.)
      Z(2)=PO2LOG
      DO 101 K=3,100
  101 Z(K)=Z(K-1)+PO2LOG-2.*LOG(DFLOAT(K-2))
      DO 102 K=3,100
  102 Z(K)=(Z(K)-LOG(DFLOAT(K-1)))
C
C CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      STOP
C
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+ABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP

  201 CONTINUE 
      if (lflag.eq.1) then
        w0= exp((2.*N-4.)*LOG(ET)+Z(N))
        do j= 1,N
          v(j)= sqrt(p(1,j)**2+p(2,j)**2+p(3,j)**2)
        enddo

        a1= 0.d0
        a3= 0.d0
        a2= 1.d0
        do j= 1,N
          a1= a1+v(j)/ET
          a2= a2*v(j)/p(4,j)
          a3= a3+v(j)*v(j)/p(4,j)/ET
        enddo
        wm= a1**(2*N-3)*a2/a3
        dj= 1.d0/w0/wm
        return
      endif
C
C THE PARAMETER VALUES ARE NOW ACCEPTED
C
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE

      DO 202 I=1,N
      call rans(RAN1)
      call rans(RAN2)
      call rans(RAN3)
      call rans(RAN4)
      C=2.*RAN1-1.
      S=SQRT(1.-C*C)
      F=TWOPI*RAN2
      Q(4,I)=-LOG(RAN3*RAN4)
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*COS(F)
  202 Q(1,I)=Q(4,I)*S*SIN(F)
C
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=SQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS
C
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
C
C CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*LOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
C
C RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
C
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=SQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=SQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(ABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
C
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*LOG(X)+LOG(WT2/WT3*ET)
C
C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
C
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D15.6)
      END

c RP: I have added subroutine flat
      subroutine flat(x,x0,x1,w,ran0)
      implicit none
      real*8 x,x0,x1,w,ran0
*
      x= x1*ran0+(1.d0-ran0)*x0
      w= (x1-x0)
      return
      end
*
      subroutine fflat(lflag,x,x0,x1,w,ran0,lw)
      implicit none
      real*8 x,x0,x1,w,ran0
      integer lflag,lw
      w= (x1-x0)
      if (lflag.eq.1) then
        if (x.lt.x0.or.x.gt.x1) goto 100
        ran0= (x-x0)/(x1-x0)
        return
      endif
*
      x= x1*ran0+(1.d0-ran0)*x0
      return
 100  w= 0.d0
      lw= 1
      return
      end
*

c RP: I have added function argsinh
      function argsinh(x)
      implicit none
      real*8 x,argsinh
      if (x.ge.0.d0) then
        argsinh= log(sqrt(x*x+1.d0)+x)
      else
        argsinh=-log(sqrt(x*x+1.d0)-x)
      endif
      return
      end     
*
c RP: I have added subroutine peaka
      subroutine peaka(a,cn,cxm,cxp,x,wj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Peaked distribution:                                         c
c   The variable x has the distribution                          c
c                                                                c
c                   1/(a+x)^cn                                   c
c                   cxm < x < cxp                                c
c                                                                c
c              INPUT                      OUTPUT                 c 
c                                                                c
c              a, cn, cxm, cxp            x, wj                  c
c                                                                c
c----------------------------------------------------------------c
*
      implicit none
      real*8 a,cn,cxm,cxp,x,wj,ran1,tjm,apx,hj
      apx= tjm(a,cn,cxm,cxp,1,ran1)
      x= -a+apx
      if (abs(1.d0-cn).gt.1.d-8) then
         wj= (apx)**cn*hj(a,cn,cxm,cxp)
      else
         wj= (apx*hj(a,cn,cxm,cxp))
      endif
      return
      end
* 
      subroutine triangle(x0,xm,xp,x,wj,ran0)
c----------------------------------------------------------------c
c                                                                c
c     Triangular distribution                                    c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 x0,xm,xp,x,wj,ran0,root,ran1,al1,alsum
      if (xm.lt.0.d0.and.xp.lt.0.d0) then
        root= (x0+xm)**2+ran0*(xp-xm)*(xp+xm+2.d0*x0)
        x   =-x0+sqrt(root)
        wj  = 0.5d0*(xp-xm)/(x+x0)*(xp+xm+2.d0*x0)
      elseif(xm.gt.0.d0.and.xp.gt.0.d0) then
        root= (-x0+xm)**2+ran0*(xp-xm)*(xp+xm-2.d0*x0)
        x   = x0-sqrt(root)
        wj  = 0.5d0*(xp-xm)/(x-x0)*(xp+xm-2.d0*x0)
      elseif(xm.lt.0.d0.and.xp.gt.0.d0) then
        alsum=-xm+xp
        al1  = -xm/alsum
        if (ran0.le.al1) then
          ran1= ran0/al1
          root= (x0+xm)**2-ran1*xm*(xm+2.d0*x0)
          x   =-x0+sqrt(root)
          wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
        else
          ran1= (ran0-al1)/(1.d0-al1)
          root= x0**2+ran1*xp*(xp-2.d0*x0)
          x   = x0-sqrt(root)
          wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
        endif
      else
        print*,'ERROR in subroutine triangle'
        stop 
      endif
      return
      end
*
      subroutine ppeaka(lflag,a,cn,cxm,cxp,x,wj,ran1,lw)
c----------------------------------------------------------------c
c                                                                c
c   Peaked distribution:                                         c
c   The variable x has the distribution                          c
c                                                                c
c                   1/(a+x)^cn                                   c
c                   cxm < x < cxp                                c
c                                                                c
c              INPUT                      OUTPUT                 c 
c                                                                c
c  lflag= 0:   a, cn, cxm, cxp, ran1      x, wj                  c
c                                                                c
c  lflag= 1:   a, cn, cxm, cxp, x         wj, ran1               c
c                                                                c
c----------------------------------------------------------------c
*
      implicit none
      real*8 a,cn,cxm,cxp,x,wj,ran1,tjm,apx,hj,rantjm
      integer lflag,lw
      if (lflag.eq.0) then
        apx= tjm(a,cn,cxm,cxp,1,ran1)
        x= -a+apx
      else
        apx= a+x
        if (x.lt.cxm.or.x.gt.cxp) goto 100
        ran1= rantjm(a,cn,cxm,cxp,1,apx)
      endif
      if (abs(1.d0-cn).gt.1.d-8) then
         wj= (apx)**cn*hj(a,cn,cxm,cxp)
      else
         wj= (apx*hj(a,cn,cxm,cxp))
      endif
      return
 100  wj= 0.d0
      lw= 1
      return
      end

      subroutine ttriangle(lflag,x0m,x0p,xm,xp,x,wj,ran0,lw)
c----------------------------------------------------------------c
c                                                                c
c     Triangular distribution                                    c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      integer lflag,lw
      real*8 x0,x0m,x0p,xm,xp,x,wj,ran0,root,ran1,al1,alsum,prot
      real*8 a,b
      data prot/1.d-20/
*
      if (lflag.eq.1) then
        if (x.lt.xm.or.x.gt.xp) goto 100
      endif
*
      if (xm.lt.0.d0.and.xp.lt.0.d0) then
        x0= dabs(x0m)
        if (lflag.eq.0) then
          root= (x0+xm)**2+ran0*(xp-xm)*(xp+xm+2.d0*x0)
          x   =-x0+sqrt(dabs(root))
        else
          a   =  0.5d0*(xp-xm)*(xp+xm+2.d0*x0)
          b   =  0.5d0*xm*(xm+2.d0*x0)
          ran0= (x*(0.5d0*x+x0)-b)/a
        endif
        if (dabs(x+x0).lt.prot) goto 100
        wj  = 0.5d0*(xp-xm)/(x+x0)*(xp+xm+2.d0*x0)
      elseif(xm.gt.0.d0.and.xp.gt.0.d0) then
        x0= dabs(x0p)
        if (lflag.eq.0) then
          root= (-x0+xm)**2+ran0*(xp-xm)*(xp+xm-2.d0*x0)
          x   = x0-sqrt(dabs(root))
        else
          a   =  0.5d0*(xp-xm)*(xp+xm-2.d0*x0)
          b   =  0.5d0*xm*(xm-2.d0*x0)
          ran0= (x*(0.5d0*x-x0)-b)/a
        endif
        if (dabs(x-x0).lt.prot) goto 100
        wj  = 0.5d0*(xp-xm)/(x-x0)*(xp+xm-2.d0*x0)
      elseif(xm.lt.0.d0.and.xp.gt.0.d0) then
        alsum=-xm+xp
        if (dabs(alsum).lt.prot) goto 100
        al1  = -xm/alsum
        if (lflag.eq.0) then
          if (ran0.le.al1) then
            x0= dabs(x0m)
            if (dabs(al1).lt.prot) goto 100
            ran1= ran0/al1
            root= (x0+xm)**2-ran1*xm*(xm+2.d0*x0)
            x   =-x0+sqrt(dabs(root))
            if (dabs(x+x0).lt.prot) goto 100
            wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
          else
            x0= dabs(x0p)
            if (dabs(1.d0-al1).lt.prot)    goto 100
            ran1= (ran0-al1)/(1.d0-al1)
            root= x0**2+ran1*xp*(xp-2.d0*x0)
            x   = x0-sqrt(dabs(root))
            if (dabs(x-x0).lt.prot) goto 100
            wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
          endif
        elseif (lflag.eq.1) then
          if (x.lt.0.d0) then
            x0= dabs(x0m)
            if (dabs(al1).lt.prot)  goto 100
            if (dabs(x+x0).lt.prot) goto 100
            wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
            a   =  0.5d0*(-xm)*(xm+2.d0*x0)
            b   =  0.5d0*xm*(xm+2.d0*x0)
            ran1= (x*(0.5d0*x+x0)-b)/a
            ran0= ran1*al1
          else
            x0= dabs(x0p)
            if (dabs(1.d0-al1).lt.prot)  goto 100
            if (dabs(x-x0).lt.prot)      goto 100
            wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
            a   =  0.5d0*(xp)*(xp-2.d0*x0)
            b   =  0.d0
            ran1= (x*(0.5d0*x-x0)-b)/a
            ran0= al1+ran1*(1.d0-al1)
          endif 
        else
          goto 101
        endif
      else
        goto 100
      endif
*
      return
 100  wj= 0.d0
      lw= 1
      return
 101  print*,'ERROR IN SUBROUTINE TTRIANGLE'
      stop
      end
*
      function aluhiw(mh,mw,mz,mt)
      implicit none
      double precision aluhiw,mh,mw,mz,mt,tmp
      double precision gtow,gtoz,gtot
      double precision b,x
      double precision pi,gfermi
      parameter (pi=3.14159265d0,gfermi=1.16639d-5)
c
      tmp=0d0
c     H-> WW
      if(mh.gt.2d0*mw) then
         x=mw**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtow=2d0*sqrt(2d0)*gfermi/32d0/pi * mh**3*(1d0-4d0*x+12d0*x**2)
     $        *b
         tmp=tmp+gtow
      endif
c     H-> ZZ
      if(mh.gt.2d0*mz) then
         x=mz**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtoz=sqrt(2d0)*gfermi/32d0/pi * mh**3*(1d0-4d0*x+12d0*x**2)*b
         tmp=tmp+gtoz
      endif
c     H-> t tbar
      if(mh.gt.2d0*mt) then
         x=mt**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtot= gfermi /4.d0/pi/dsqrt(2.d0)*mh*3.d0*mt**2 * b
         tmp=tmp+gtot
      endif
      aluhiw=tmp
      end


      subroutine rans(ran)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Random number generator                                    c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ran,rangen
      ran= rangen(1)
      end

      subroutine randa(rn)
      implicit none
      real *8 rn,rangen
      rn=rangen(1)
      end

*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION RANGEN(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RANGEN,RANSET,RANGET
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12345,67890/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANGEN=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
C-----------------------------------------------------------------------
      ENTRY RANSET(JSEED)
C-----------------------------------------------------------------------
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) then
         write(6,*) 'Jseeds=0, wrong settings for RANSET'
         stop
      endif
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
C-----------------------------------------------------------------------
      ENTRY RANGET(JSEED)
C-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END

*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION RANGEN2(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RANGEN2,RANSET2,RANGET2
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12347,67892/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANGEN2=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
C-----------------------------------------------------------------------
      ENTRY RANSET2(JSEED)
C-----------------------------------------------------------------------
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) then
         write(6,*) 'Jseeds=0, wrong settings for RANSET'
         stop
      endif
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
C-----------------------------------------------------------------------
      ENTRY RANGET2(JSEED)
C-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END


C-----------------------------------------------------------------------
      subroutine aluend(iunit)
c     go to end of file
C-----------------------------------------------------------------------
      ios = 0    
      dowhile(ios.eq.0)
         read(unit=iunit,fmt='(1x)',iostat=ios)
      enddo                        
      end
C-----------------------------------------------------------------------
      subroutine alugun(n)
c     look for next available unit
C-----------------------------------------------------------------------
      implicit none
      integer n,i
      logical yes
      do i=10,100
        inquire(unit=i,opened=yes)
        if(.not.yes) goto 10
      enddo
      write(6,*) 'no free units to write to available, stop'
      stop
 10   n=i
      end
C---------------- -------------------------------------------------------
      subroutine alustn(string,num)
c- writes the number num on the string string starting at the blank
c- following the last non-blank character
      character * (*) string
      character * 20 tmp
      integer aluisl
      l = len(string)
      write(tmp,'(i15)')num
      j=1
      dowhile(tmp(j:j).eq.' ')
        j=j+1
      enddo
      ipos = aluisl(string)
      ito = ipos+1+(15-j)
      if(ito.gt.l) then
         write(6,*)'error, string too short'
         write(6,*) string
         stop
      endif
      string(ipos+1:ito)=tmp(j:)
      end

      function aluisl(string)
c returns the position of the last non-blank character in string
      integer aluisl
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      aluisl = i
      end

      subroutine alustc(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      integer aluisl
      character *(*) str1,str2,str
      l1=aluisl(str1)
      l2=aluisl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(6,*) 'error: l1+l2>l in alustc'
          write(6,*) 'l1=',l1,' str1=',str1
          write(6,*) 'l2=',l2,' str2=',str2
          write(6,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end

c-------------------------------------------------------------------
      function unwgt(w,wmax)
c     routine for the unweighting of weighted events.
c     unwgt=.true.   ==> event was selected
c-------------------------------------------------------------------
      implicit none
      logical unwgt
      double precision w,wmax,rn,rangen2
      unwgt=.false.
      rn=rangen2(1)
      if(w.ge.rn*wmax) unwgt=.true.
      end



      subroutine alstio
      implicit none
      include 'alpgen.inc'
c  For MPI
#ifdef USE_MPI
      integer mpirank,filestat
      character*8 rankstring
      common/mpi/mpirank,rankstring
#endif
      character*550 tmpstr
      integer l,aluisl
      if(imode.lt.2) then
         call alustc(fname_mpirank,'.stat',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Statistics and results written to ',tmpstr(1:l)
         call alugun(niosta)
#ifndef USE_MPI
         open(unit=niosta,file=tmpstr,status='unknown')
#else
c for parallel running we don't want many files being output
c so we redirect this to the stdout instead
         niosta=6
#endif
         call alustc(fname_mpirank,'.top',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Topdrawer plots (if any) written to ',tmpstr(1:l)
         topfile=tmpstr
      endif
      if(imode.eq.1) then

         call alustc(fname_mpirank,'.par',tmpstr)
         call alugun(niopar)
#ifdef USE_MPI
         open(unit=niopar,file=tmpstr,status='unknown',iostat=filestat)
         if(filestat /= 0) then
            write(6,*) 'Error opening file: ',tmpstr
            call exit(-1)
         endif
#else
         open(unit=niopar,file=tmpstr,status='unknown')
#endif
         l=aluisl(tmpstr)
         write(6,*) 'Run parameters and diagnostics written to '
     $        ,tmpstr(1:l)
#ifdef USE_MPI
         call alustc(fname_mpirank,'.wgt',tmpstr)
#else
         call alustc(fname,'.wgt',tmpstr)
#endif
         call alugun(niowgt)
#ifdef USE_MPI
         open(unit=niowgt,file=tmpstr,status='unknown',iostat=filestat)
         if(filestat /= 0) then
            write(6,*) 'Error opening file: ',tmpstr
            call exit(-1)
         endif
#else
         open(unit=niowgt,file=tmpstr,status='unknown')
#endif
         l=aluisl(tmpstr)
         write(6,*) 'Wgted events are written to ',tmpstr(1:l)
      elseif(imode.eq.2) then
         call alustc(fname,'_unw.par',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Parameters and results written to ',tmpstr(1:l)
         call alugun(niosta)
#ifdef USE_MPI
c        only rank 0 needs to ouput a par file
         if(mpirank.eq.0) then
            open(unit=niosta,file=tmpstr,status='unknown',
     $           iostat=filestat)
            if(filestat /= 0) then
               write(6,*) 'Error opening file: ',tmpstr
               call exit(-1)
            endif
         else
            niosta = 6
         endif
#else
         open(unit=niosta,file=tmpstr,status='unknown')
#endif
         call alustc(fname_mpirank,'_unw.top',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Topdrawer plots (if any) written to ',tmpstr(1:l)
         topfile=tmpstr
         call alustc(fname_mpirank,'.par',tmpstr)
         call alugun(niopar)
#ifdef USE_MPI
         open(unit=niopar,file=tmpstr,status='unknown',iostat=filestat,
     $        action='READ')
         if(filestat /= 0) then
            write(6,*) 'Error opening file: ',tmpstr
            call exit(-1)
         endif
#else
         open(unit=niopar,file=tmpstr,status='unknown')
#endif
         call alustc(fname_mpirank,'.wgt',tmpstr)
         call alugun(niowgt)
#ifdef USE_MPI
         open(unit=niowgt,file=tmpstr,status='unknown',iostat=filestat,
     $        action='READ')
         if(filestat /= 0) then
            write(6,*) 'Error opening file: ',tmpstr
            call exit(-1)
         endif
#else
         open(unit=niowgt,file=tmpstr,status='unknown')
#endif
         l=aluisl(tmpstr)
         write(6,*) 'Wgted events are read from ',tmpstr(1:l)

         call alustc(fname_mpirank,'.unw',tmpstr)
         call alugun(niounw)
#ifdef USE_MPI
         open(unit=niounw,file=tmpstr,status='unknown',iostat=filestat)
         if(filestat /= 0) then
            write(6,*) 'Error opening file: ',tmpstr
            call exit(-1)
         endif
#else
         open(unit=niounw,file=tmpstr,status='unknown')
#endif
         l=aluisl(tmpstr)
         write(6,*) 'Unwgted events are written to ',tmpstr(1:l)
      endif
      end


      subroutine refinemom(n,p,iflag)                           
*                                                                    
c-    Refine the momenta in such a way they are                 
c-    on-shell and energy-momentum conserving up to             
c-    the computer precision.                                   
c-    n is the total number of momenta (initial+final state)    
*                                                                    
      implicit none                                             
      integer maxpar,n,i,iflag                                  
      parameter (maxpar= 10)                                    
      real*8 p(5,maxpar),ptotx,ptoty,ptotz,etot                 
*                                                               
      iflag= 0                                                  
      ptotx= 0.d0                                               
      ptoty= 0.d0                                               
      ptotz= 0.d0                                               
      do i= 3,n                                                 
         ptotx= ptotx+p(1,i)                                    
         ptoty= ptoty+p(2,i)                                    
         ptotz= ptotz+p(3,i)                                    
      enddo                                                     
*                                                                    
      etot= 0.d0                                                
      do i= 3,n                                                 
        p(1,i)= p(1,i)-ptotx/dfloat(n-2)                        
        p(2,i)= p(2,i)-ptoty/dfloat(n-2)                        
        p(4,i)= sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+p(5,i)**2)   
        etot  = etot+p(4,i)                                     
      enddo                                                     
*                                                                    
      if (dabs(ptotz).gt.etot) then                             
        iflag= 1                                                
      else                                                      
        p(4,1)= 0.5d0*(etot+ptotz)                              
        p(3,1)= p(4,1)                                          
        p(4,2)= 0.5d0*(etot-ptotz)                              
        p(3,2)=-p(4,2)                                          
      endif                                                     
      return                                                    
      end                                                       




c-------------------------------------------------------------------
      SUBROUTINE USRTIME(NOW)
c-------------------------------------------------------------------
c     this subroutine should return the current time, in the format
c     "Thu Apr  4 12:10:45 2002". 
c     The present version is based on use of the Intrinsic Time and
c     CTime functions, supported by most Unix and GNU compilers.  The
c     user may need to change them should they not be supported by her
c     compiler. The sole purpose of this routine is to provide the time
c     information when outputting the topdrawer histograms, so if you
c     don't want to bother finding a Time routine supported by your
c     compiler, you can uncomment the dummy value assigned below to NOW
c     ,and comment away the call to ctime
c     
      IMPLICIT NONE
      integer time
      character*24 CTIME,now
      now='Day Mon XX hh:mm:ss yyyy'
c      now=ctime(time())
      end


*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULB4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(1)*PS(1)+PI(2)*PS(2)
     &         +PI(3)*PS(3)+PI(4)*PS(4))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) + FN*PS(1)
        PF(2)= PI(2) + FN*PS(2)
        PF(3)= PI(3) + FN*PS(3)
        PF(4)= PF4
      END IF
      END
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULF4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN LAB) INTO PF (IN REST FRAME OF PS)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(4)*PS(4)-PI(3)*PS(3)
     &         -PI(2)*PS(2)-PI(1)*PS(1))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) - FN*PS(1)
        PF(2)= PI(2) - FN*PS(2)
        PF(3)= PI(3) - FN*PS(3)
        PF(4)= PF4
      END IF
      END
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULOB(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4,5) = (PX,PY,PZ,E,M)
C-----------------------------------------------------------------------
      DOUBLE PRECISION PS(5),PI(5),PF(5)
      CALL ALULB4(PS,PI,PF)
      PF(5)= PI(5)
      END
CDECK  ID>, ALULOF.
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULOF(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN LAB) INTO PF (IN REST FRAME OF PS)
C     N.B. P(1,2,3,4,5) = (PX,PY,PZ,E,M)
C-----------------------------------------------------------------------
      DOUBLE PRECISION PS(5),PI(5),PF(5)
      CALL ALULF4(PS,PI,PF)
      PF(5)= PI(5)
      END
CDECK  ID>, ALULOR.
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Giovanni Abbiendi & Luca Stanco
C-----------------------------------------------------------------------
      SUBROUTINE ALULOR (TRANSF,PI,PF)
C-----------------------------------------------------------------------
C     Makes the ALULOR transformation specified by TRANSF on the
C     quadrivector PI(5), giving PF(5).
C-----------------------------------------------------------------------
      DOUBLE PRECISION TRANSF(4,4),PI(5),PF(5)
      INTEGER I,J
      DO 1 I=1,5
        PF(I)=0.D0
    1 CONTINUE
      DO 3 I=1,4
       DO 2 J=1,4
         PF(I) = PF(I) + TRANSF(I,J) * PI(J)
    2  CONTINUE
    3 CONTINUE
      PF(5) = PI(5)
      RETURN
      END
CDECK  ID>, ALUMAS.
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALUMAS(P)
C-----------------------------------------------------------------------
C     PUTS INVARIANT MASS IN 5TH COMPONENT OF VECTOR
C     (NEGATIVE SIGN IF SPACELIKE)
C-----------------------------------------------------------------------
      DOUBLE PRECISION ALUSQR,P(5)
      P(5)=SQRT((P(4)+P(3))*(P(4)-P(3))-P(1)**2-P(2)**2)
      END
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION ALUPCM(EM0,EM1,EM2)
C-----------------------------------------------------------------------
C     C.M. MOMENTUM FOR DECAY MASSES EM0 -> EM1 + EM2
C     SET TO -1 BELOW THRESHOLD
C-----------------------------------------------------------------------
      DOUBLE PRECISION ALUPCM,EM0,EM1,EM2,EMS,EMD
      EMS=ABS(EM1+EM2)
      EMD=ABS(EM1-EM2)
      IF (EM0.LT.EMS.OR.EM0.LT.EMD) THEN
        ALUPCM=-1.
      ELSEIF (EM0.EQ.EMS.OR.EM0.EQ.EMD) THEN
        ALUPCM=0.
      ELSE
        ALUPCM=SQRT((EM0+EMD)*(EM0-EMD)*
     &              (EM0+EMS)*(EM0-EMS))*.5/EM0
      ENDIF
      END

      subroutine rescms(p,p1,p2,m1,m2)
      implicit none
      integer il
      double precision p(5),p1(5),p2(5),m1,m2
      double precision po1(5),po2(5),m,mo1,mo2,pcm,pcmo
      double precision alupcm
      m=p(5)
      mo1=p1(5)
      mo2=p2(5)
      call alulof(p,p1,po1)
      call alulof(p,p2,po2)
      if(max(mo1,mo2).lt.1d-6) then
        pcmo=m/2d0
      else 
        pcmo=alupcm(m,mo1,mo2)
      endif 
      pcm=alupcm(m,m1,m2)
      do il=1,3
        po1(il)=pcm/pcmo*po1(il)
        po2(il)=pcm/pcmo*po2(il)
      enddo
      po1(4)=sqrt(pcm**2+m1**2)
      po2(4)=sqrt(pcm**2+m2**2)
      po1(5)=m1
      po2(5)=m2
      call alulob(p,po1,p1)
      call alulob(p,po2,p2)
      end


c     start clustering routines
c
      subroutine alsclu
c     setup variables for clustering:
c     iopt=1: quarks and gluons only
c     iopt=2: flavour clustering, quarks/gluons and EW particles
      implicit none
      include 'alpgen.inc'
c     locals
      integer i,j,if1,ini,nfcol,nweak
      data ini/0/
c
c      if(ihrd.gt.4) return
c     count coloured particles (quarks and gluons)
      if(ini.eq.0) then
        nclext=0
        do i=1,npart
          if(abs(ifl(i)).le.6.or.ifl(i).eq.21) nclext=nclext+1
        enddo
c non strongly interacting final state particles
        nweak=nw+nh+nz+nph
c strongly interacting final state particles
        nfcol=nclext-2
        ncltot=2*nfcol+min(2,nweak)
c determine total number of clusters 
c
c w/z qq+jets
        if(ihrd.le.2) then
          ncltot=4+2*njets
        elseif(ihrd.le.4) then
          ncltot=2+2*njets
        elseif(ihrd.eq.5) then
          if(nh.eq.0) then
            ncltot=2+2*njets
          else
            if(nh.eq.nweak) then
c     qq->Hqq => purely EW
              ncltot=2*(nfcol-2)+4
            else
              ncltot=2+2*nfcol
            endif
          endif
c     2Q
        elseif(ihrd.eq.6) then
          ncltot=4+2*njets
c     4Q
        elseif(ihrd.eq.7) then
          ncltot=6+2*njets
c     QQh
        elseif(ihrd.eq.8) then
          ncltot=4+2*njets
c     Njet
        elseif(ihrd.eq.9) then
          ncltot=4+2*(njets-2)
c     W c jets
        elseif(ihrd.eq.10) then
          ncltot=3 +2*njets
c     gamma jets
        elseif(ihrd.eq.11) then
          ncltot=3+2*(njets-1)
c     H jets
        elseif(ihrd.eq.12) then
          ncltot=2+2*njets
c     W gamma jets
        elseif(ihrd.eq.14) then
          ncltot=2+2*njets
c     W gamma QQ jets
        elseif(ihrd.eq.15) then
          ncltot=4+2*njets
        endif
c
        naspow=ncltot-nclext
        ini=1
      endif
      if(nclext.eq.ncltot) return
c
      kres=qsq
c
      do i=1,nclext
        clkt(1,i)=kres
        clkt(2,i)=kres
        icst(i)=0
        icmot(i)=0
        if1=ifl(i)
        if(if1.eq.21) if1=0
        icfl(i)=if1
        icfs(i)=1
        if(i.le.2) icfs(i)=-1
        do j=1,4
          clp(j,i)=p(j,i)
        enddo
        do j=1,2
          icdau(j,i)=i
          iccol(j,i)=icu(j,i)
        enddo
      enddo
      do i=nclext+1,ncltot
        clkt(1,i)=kres
        clkt(2,i)=kres
        icst(i)=0
        icmot(i)=0
        icfl(i)=0
        icfs(i)=1
        do j=1,4
          clp(j,i)=0
        enddo
        do j=1,2
          icdau(j,i)=i
          iccol(j,i)=0
        enddo
      enddo
      end
      
      
      subroutine clumrg(p1,p2,is1,is2,pclu,mrgopt)
c     merge pseudoparticle momenta p1&p2 into new cluster momentum
      implicit none
c     inputs
      integer is1,is2,mrgopt
      double precision p1(4),p2(4),pclu(4)
c     locals
      integer k,is
      is=is1*is2
      do k=1,4
        pclu(k)=is1*p1(k)+is2*p2(k)
        if(is.lt.0) pclu(k)=-pclu(k)
      enddo
      if(mrgopt.eq.2) then
c keep massless cluster
        pclu(4)=sqrt((is*p1(1)+p2(1))**2+(is*p1(2)+p2(2))**2+(is*p1(3)
     $       +p2(3))**2)
      endif
      end


      subroutine alcclu(iflag)
c colour driven clustering
      implicit none
      include 'alpgen.inc'
c inputs
      integer iflag
c local
      integer i,j,nmax,nclmax,is,imin,jmin,istart,if1,if2,ic1(2)
     $     ,ic2(2),icol1,icol2,icmin1,icmin2
      double precision ktmin,kperp,tmp,tiny
      parameter (tiny=1d-2)
c debugging
c      integer nev
c      data nev/0/
c
c initialise
c      nev=nev+1
      iflag=0
      nclmax=ncltot
c      write(6,*) 'event=',nev
      nmax=nclext
      if(nmax.eq.nclmax) return
 1    continue
      ktmin=1d12
      do i=1,nmax-1
        if(icst(i).eq.1) goto 10
        istart=i
        if1=icfl(i)*icfs(i)
        ic1(1)=iccol(1,i)
        ic1(2)=iccol(2,i)
c inhibit clustering of initial-initial states
        if(i.eq.1) istart=2
        do j=istart+1,nmax
          if(icst(j).eq.1) goto 9
          if(icfs(i).lt.0. and .icfs(j).lt.0) goto 9
          if2=icfl(j)*icfs(j)
          ic2(1)=iccol(1,j)
          ic2(2)=iccol(2,j)
          is=icfs(i)*icfs(j)
c     check clusterability of i+j and assign colours to potential new
C     cluster
c     do not allow W-g clusters
          if( (if1.eq.0.and.if2.eq.-1000) . or.
     $        (if2.eq.0.and.if1.eq.-1000))  goto 9
          if(iccol(1,i).eq.-1) then
c     i is an EW object, cluster inherits j colours
            icol1=iccol(1,j)
            icol2=iccol(2,j)
          elseif(iccol(1,j).eq.-1) then
c     j is an EW object, cluster inherits i colours
            icol1=iccol(1,i)
            icol2=iccol(2,i)
          elseif(is.lt.0) then
c initial-final cluster first
            if((if1+if2)*if1*if2.ne.0) goto 9
            if(ic1(1).ne.0.and.ic1(1).eq.ic2(1)) then
              icol1=ic1(2)*(1+icfs(i))/2+ic2(2)*(1+icfs(j))/2
              icol2=ic1(2)*(1-icfs(i))/2+ic2(2)*(1-icfs(j))/2
c              write(6,*) 'passed 1'
            elseif(ic1(2).ne.0.and.ic1(2).eq.ic2(2)) then
              icol1=ic1(1)*(1-icfs(i))/2+ic2(1)*(1-icfs(j))/2
              icol2=ic1(1)*(1+icfs(i))/2+ic2(1)*(1+icfs(j))/2
c              write(6,*) 'passed 2'
            elseif(ic1(1).eq.ic1(2)) then
              icol1=ic2(1)*(1-icfs(j))/2+ic2(2)*(1+icfs(j))/2
              icol2=ic2(2)*(1-icfs(j))/2+ic2(1)*(1+icfs(j))/2
c              write(6,*) 'passed 3'
            elseif(ic2(1).eq.ic2(2)) then
              icol1=ic1(1)*(1-icfs(i))/2+ic1(2)*(1+icfs(i))/2
              icol2=ic1(2)*(1-icfs(i))/2+ic1(1)*(1+icfs(i))/2
c              write(6,*) 'passed 4'
c q -> q g(cluster)
            elseif(if1+if2.eq.0.and.if1.ne.0) then
              if(ic1(1).ne.0) then
                icol1=ic1(1)
                icol2=ic2(1)
              elseif(ic1(2).ne.0) then
                icol1=ic2(2)
                icol2=ic1(2)
              endif
c              write(6,*) 'passed 5'
            else
              goto 9
            endif
          else
c final-final cluster then
            if((if1+if2)*if1*if2.ne.0) goto 9
            if(ic1(1).ne.0.and.ic1(1).eq.ic2(2)) then
              icol1=ic2(1)
              icol2=ic1(2)
c              write(6,*) 'passed 11'
            elseif(ic1(2).ne.0.and.ic1(2).eq.ic2(1)) then
              icol1=ic1(1)
              icol2=ic2(2)
c              write(6,*) 'passed 12'
            elseif(ic1(1).eq.ic1(2)) then
              icol1=ic2(1)
              icol2=ic2(2)
c              write(6,*) 'passed 13'
            elseif(ic2(1).eq.ic2(2)) then
              icol1=ic1(1)
              icol2=ic1(2)
c              write(6,*) 'passed 14'
c g-> q qbar clusters
            elseif(if1+if2.eq.0.and.if1.ne.0) then
              if(ic1(1).ne.0) then
                icol1=ic1(1)
                icol2=ic2(2)
              elseif(ic2(1).ne.0) then
                icol1=ic2(1)
                icol2=ic1(2)
              endif
c              write(6,*) 'passed 15'
            else
              goto 9
            endif
          endif
c
 8        tmp=kperp(clp(1,i),clp(1,j),icfs(i),icfs(j),cluopt,ktfac)
c          write(6,*) (clp(imin,i),imin=1,3)
c          write(6,*) (clp(imin,j),imin=1,3)
c          write(6,*) 'i,j:',i,j,' kt=',sqrt(tmp)
c          if(tmp.lt.kres) then
c below resolution, reject event
c            iflag=1
c            return
c          endif
c     resolve the symmetry between merging with either beam-like
C     particle. Make more probable the configuration with jet and beam
C     traveling in the same z direction
c          if(is.lt.0) then
c            tmp=tmp-sign(tiny,clp(3,i)*clp(3,j))
c          endif
          ktmin=min(ktmin,tmp)
          if(abs(ktmin-tmp).lt.1d-8*tmp) then
            imin=i
            jmin=j
            if(icol1.eq.icol2) then
              icmin1=0
              icmin2=0
            else
              icmin1=icol1
              icmin2=icol2
            endif
c            write(6,*) 'new min kt:',ktmin,' pair:',i,j
          endif
 9        continue
        enddo
 10     continue 
      enddo
c create new cluster
c      write(6,*) 'new cluster:',imin,jmin,' ktmin=',ktmin
c no new cluster found
      if(ktmin.eq.1d12) goto 100
      nmax=nmax+1
      icdau(1,nmax)=imin
      icdau(2,nmax)=jmin
      clkt(2,nmax)=ktmin
      is=icfs(imin)*icfs(jmin)
      icfl(nmax)=is*(icfs(imin)*icfl(imin)+icfs(jmin)*icfl(jmin))
      icfs(nmax)=is
c     merge pesudoparticles into cluster
      call clumrg(clp(1,imin),clp(1,jmin),icfs(imin),icfs(jmin),clp(1
     $     ,nmax),mrgopt)
c update i,j status
      icst(imin)=1
      icst(jmin)=1
      clkt(1,imin)=ktmin
      clkt(1,jmin)=ktmin
      icmot(imin)=nmax
      icmot(jmin)=nmax
c cluster colour
      iccol(1,nmax)=icmin1
      iccol(2,nmax)=icmin2
      if(nmax.lt.nclmax) goto 1
      return
 100  continue
c leave incomplete clusters with scale at ktres
c      write(6,*) 'incomplete cluster'
c completed clustering
      end

      subroutine cktmin
c determine ktmin for the event
      implicit none
      include 'alpgen.inc'
c inputs
      double precision ktmin
c local
      integer i,j,nmax,nclmax,is,imin,jmin,istart,ifs1,ifs2,ic1(2)
     $     ,ic2(2),ifl1,ifl2
      double precision kperp,tmp
c
      call alsclu
      nmax=nclext
      ktmin=qsq
      do i=1,nmax-1
        istart=i
        ifs1=icfs(i)
        ifl1=icfl(i)*ifs1
        ic1(1)=iccol(1,i)
        ic1(2)=iccol(2,i)
        if(ic1(1).eq.0.and.ic1(2).eq.0) goto 10
        do j=istart+1,nmax
          ifs2=icfs(j)
          ifl2=icfl(j)*ifs2
          ic2(1)=iccol(1,j)
          ic2(2)=iccol(2,j)
          if(ic2(1).eq.0.and.ic2(2).eq.0) goto 9
c          is=ifs1+ifs2
c     check clusterability of i+j:
c accept all  qg or gg clusters, check flavour for  qqbar clusters
          if(ifl1*ifl2*(ifl1+ifl2).ne.0) goto 9
c            if(ifl1+ifl2.eq.0) goto 8
c          else
c            goto 8
c            if(is.eq.0) then
c initial-final
c              if(ic1(1).eq.ic2(1).or.ic1(2).eq.ic2(2)) goto 8
c            else
c final-final or initial-initial
c              if(ic1(1).eq.ic2(2).or.ic1(2).eq.ic2(1)) goto 8
c            endif
c          endif
c          goto 9
c
 8        tmp=kperp(clp(1,i),clp(1,j),ifs1,ifs2,cluopt,ktfac)
c          write(6,*) 'i,j:',i,j,' kt=',tmp
          ktmin=min(ktmin,tmp)
 9        continue
        enddo
 10     continue 
      enddo
      kres=ktmin
      end



      function kperp(p1,p2,is1,is2,iopt,qscale)
      implicit none
      double precision kperp,p1(4),p2(4),qscale
      integer is1,is2,iopt
c local variables
      double precision pt1,pt2,y1,y2,dphi,mt1,mt2,q1,q2,m1,m2,tmp
c commons
c
      pt1=p1(1)**2+p1(2)**2
      pt2=p2(1)**2+p2(2)**2
      mt1=p1(4)**2-p1(3)**2
      mt2=p2(4)**2-p2(3)**2
      m1=mt1-pt1
      m2=mt2-pt2
      if(iopt.eq.1) then
        q1=pt1
        q2=pt2
      elseif(iopt.eq.2) then
        q1=mt1
        q2=mt2
      endif
      if(is1.gt.0.and.is2.gt.0) then
c     final-final
        y1=-0.5*log((p1(4)-p1(3))/(p1(4)+p1(3)))
        y2=-0.5*log((p2(4)-p2(3))/(p2(4)+p2(3)))
        dphi=acos( (p1(1)*p2(1)+p1(2)*p2(2) )/sqrt(pt1*pt2+1d-3))
        kperp=((y2-y1)**2+dphi**2)*min(q1,q2)
c     initial-final
      elseif(is1.lt.0.and.is2.gt.0) then
c     kperp=mt2
        kperp=q2+abs(q1-pt1)
c     final-initial
      elseif(is2.lt.0.and.is1.gt.0) then
c     kperp=mt1
        kperp=q1+abs(q2-pt2)
c     initial-initial
      else
        kperp=(p1(4)+p2(4))**2-(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3
     $       )+p2(3))**2
      endif
c     protect against sqrt(kperp) falling below the mass of final state
C     and below 2 gev
c
      kperp=qscale**2*max(kperp,abs(m1)+abs(m2))
      kperp=max(kperp,4d0)
      end

      subroutine alasrs(rewgt)
      implicit none
      include 'alpgen.inc'
c arguments
      double precision rewgt
c locals
      integer i
      double precision alfas_clu
c
      rewgt=1d0
      do i=nclext+1,ncltot
c        if(clkt(2,i).ne.kres) write(2,*) 'rewgt',i,alfas(clkt(2,i),xlclu
c     $       ,lpclu,-1)/asmax
        if(clkt(2,i).ne.kres) rewgt=rewgt*alfas_clu(clkt(2,i),xlclu
     $       ,lpclu,-1)/asmax
      enddo
      end


      subroutine rstcol(icol)
c set colour flags for and from clustering
      implicit none
      include 'alpgen.inc'
      integer i,j,icol
      if(icol.eq.1) then
        do j=1,nclext
          do i=1,2
            icu(i,j)=0
          enddo
        enddo
      elseif(icol.eq.2) then
        do j=1,nclext
          do i=1,2
            icu(i,j)=iccol(i,j)
          enddo
        enddo
      endif
      end

      subroutine alpclu(rewgt)
      implicit none
      include 'alpgen.inc'
      integer iflag
      double precision rewgt,w
      rewgt=1d0
      call alsclu
      call alcclu(iflag)
      if(iflag.eq.1) then
        rewgt=0d0
        return
      endif
c rescale alphas only during the weighted event generation
c      if(imode.le.1) then
      call alasrs(w)
      rewgt=rewgt*w
c      endif
      end

c function to setup MPI
#ifdef USE_MPI
      subroutine INITIALIZE_MPI
      include 'mpif.h'
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      integer ierr
      real timestamp

      call MPI_INIT(ierr)
      if(ierr.ne.0) then
         write(6,*) 'Error initializing mpi'
         call exit(-1)
      endif
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiworldsize, ierr)
      if(ierr.ne.0) then
         write(6,*) 'Error getting mpi comm world size'
         call exit(-1)
      endif
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, ierr)
      if(ierr.ne.0) then
         write(6,*) 'Error getting mpi comm rank'
         call exit(-1)
      endif

c     redirect stdout to dev/null for all ranks but 0
      if(mpirank.ne.0) then
         open(unit=6,file='/dev/null')
      endif


      
      write(rankstring,'(i8.8)')  mpirank

      call cpu_time(timestamp)
      write(6,*) ' INITIALIZE RANK ',mpirank,' AT ',timestamp
       
      end
c finalize MPI
      subroutine FINALIZE_MPI
      include 'mpif.h'
      integer mpirank,mpiworldsize
      character*8 rankstring
      common/mpi/mpirank,mpiworldsize,rankstring
      integer ierr
      real timestamp
      call cpu_time(timestamp) 
      write(6,*) ' FINALIZE RANK ',rankstring,' AT ',timestamp
      call MPI_FINALIZE(ierr)
      if(ierr.ne.0) then
         write(6,*) 'Error finalizing MPI'
         call exit(-1)
      endif
      end
#endif
