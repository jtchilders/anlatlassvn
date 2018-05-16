      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_par.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_physpar.h'
      logical verbose
      parameter(verbose=.true.)
      integer aemrun
      real *8 alfaem,pwhg_alphas,rdummy,totbr
      external alfaem,pwhg_alphas
      integer i,j,idummy
      real *8 alphaem_inv
      common/calphaem_inv/alphaem_inv
      real *8 powheginput
      external powheginput

      par_isrtinycsi = 1d-5
      par_isrtinyy =   1d-6
      par_fsrtinycsi = 1d-5
      par_fsrtinyy =   1d-6

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     number of light flavors
      st_nlight = 5

c     setting physical parameters
      write(*,*) 'POWHEG: loading-setting physical parameters'

c     branching ratio (used only in LH event file)
      call pickwdecay(-1000,rdummy,idummy,rdummy,idummy,rdummy,idummy
     $     ,rdummy,totbr)
      rad_branching=totbr

c     top mass
      topmass_pow=powheginput('topmass')
      if(topmass_pow.lt.0) then
         write(*,*) 'Input error: topmass ',topmass_pow
         call exit(1)
      endif

c     ew parameters
c     true inputs are wmass, sthw2, alphaem

c     W mass
      wmass_pow=powheginput('#wmass')
      if(wmass_pow.lt.0) wmass_pow=80.4d0

c     sthw2
      sthw2_pow=powheginput('#sthw2')
      if(sthw2_pow.lt.0) sthw2_pow=0.23113d0

c     alphaem
c     typical values for alphaem:
c     Thompson value:    1/137.0359895d0
c     at z mass (91.188) 1/127.934 (?)
c     at top mass (175)  1/127.011989

      aemrun=0
c     definition of alphaem_pow value, according to aemrun
      if(aemrun.eq.0) then
         alphaem_inv=powheginput('#alphaem_inv')
         if(alphaem_inv.lt.0) alphaem_inv=127.011989
         alphaem_pow=1d0/alphaem_inv
         zmass_pow=91.188d0     !Not relevant in POWHEG; needed only by set_madgraph_parameters
      elseif(aemrun.eq.1) then
         write(*,*) 'Invalid option for aemrun: program stops'
         call exit(1)
c$$$c     alphaem_pow is evaluated at the top mass value using the alfaem function.
c$$$c     In this case zmass is needed by the alfaem function to set a reference
c$$$c     point for the running of alphaem. 
c$$$c     zmass needed also by set_madgraph_parameters.
c$$$c     This reference value is read and used by the function alfaem itself
c$$$c     that will assume alfaem(zmass)=1/alphaem_inv.
c$$$         alphaem_inv=127.934
c$$$         zmass_pow=92d0
c$$$         alphaem_pow=alfaem(topmass_pow**2)
      else
         write(*,*) 'Error while setting aemrun'
         call exit(1)
      endif

cccccccccccccccccccccccccccccccccccc
c     !: to check with Chris code
c      alphaem_pow=1./127.011989
cccccccccccccccccccccccccccccccccccc

c     CKM matrix entries
      CKM_pow(1,1)= powheginput('#CKM_Vud')
      CKM_pow(1,2)= powheginput('#CKM_Vus') 
      CKM_pow(1,3)= powheginput('#CKM_Vub') 
      CKM_pow(2,1)= powheginput('#CKM_Vcd') 
      CKM_pow(2,2)= powheginput('#CKM_Vcs') 
      CKM_pow(2,3)= powheginput('#CKM_Vcb') 
      CKM_pow(3,1)= powheginput('#CKM_Vtd') 
      CKM_pow(3,2)= powheginput('#CKM_Vts') 
      CKM_pow(3,3)= powheginput('#CKM_Vtb') 

      if(CKM_pow(1,1).lt.0) CKM_pow(1,1)= 0.9740
      if(CKM_pow(1,2).lt.0) CKM_pow(1,2)= 0.2225
      if(CKM_pow(1,3).lt.0) CKM_pow(1,3)= 1d-6

      if(CKM_pow(2,1).lt.0) CKM_pow(2,1)= 0.2225
      if(CKM_pow(2,2).lt.0) CKM_pow(2,2)= 0.9740
      if(CKM_pow(2,3).lt.0) CKM_pow(2,3)= 1d-6

      if(CKM_pow(3,1).lt.0) CKM_pow(3,1)= 1d-6
      if(CKM_pow(3,2).lt.0) CKM_pow(3,2)= 1d-6
      if(CKM_pow(3,3).lt.0) CKM_pow(3,3)= 1d0

c     W width (only for decay)
      wwidth_pow=powheginput('#wwidth')
      if(wwidth_pow.lt.0) wwidth_pow=2.141d0

c     top width (only for decay)
      topwidth_pow=powheginput('#topwidth')
      if(topwidth_pow.lt.0) topwidth_pow=1.7d0

c     masses for reshuffling procedure of
c     outgoing particles
      physpar_ml(1)=powheginput('#lhfm/emass')
      if(physpar_ml(1).lt.0) physpar_ml(1)=0.000511
      physpar_ml(2)=powheginput('#lhfm/mumass')
      if(physpar_ml(2).lt.0) physpar_ml(2)=0.1056
      physpar_ml(3)=powheginput('#lhfm/taumass')
      if(physpar_ml(3).lt.0) physpar_ml(3)=1.777

      physpar_mq(4)=powheginput('#lhfm/cmass')
      if(physpar_mq(4).lt.0) physpar_mq(4)=1.5
      physpar_mq(5)=powheginput('#lhfm/bmass')
      if(physpar_mq(5).lt.0) physpar_mq(5)=5.0

      do i=1,3
         do j=1,3
            if(CKM_pow(i,j).lt.0d0) then
               write(*,*) 'Input error: CKM (i,j)= ',i,j
               call exit(1)
            endif
         enddo
      enddo

      do i=1,6
         do j=1,6
            CKM(i,j)=0d0
         enddo
      enddo
      CKM(2,1)=CKM_pow(1,1)
      CKM(1,2)=CKM_pow(1,1)
      CKM(2,3)=CKM_pow(1,2)
      CKM(3,2)=CKM_pow(1,2)
      CKM(2,5)=CKM_pow(1,3)
      CKM(5,2)=CKM_pow(1,3)

      CKM(4,1)=CKM_pow(2,1)
      CKM(1,4)=CKM_pow(2,1)
      CKM(4,3)=CKM_pow(2,2)
      CKM(3,4)=CKM_pow(2,2)
      CKM(4,5)=CKM_pow(2,3)
      CKM(5,4)=CKM_pow(2,3)

      CKM(6,1)=CKM_pow(3,1)
      CKM(1,6)=CKM_pow(3,1)
      CKM(6,3)=CKM_pow(3,2)
      CKM(3,6)=CKM_pow(3,2)
      CKM(6,5)=CKM_pow(3,3)
      CKM(5,6)=CKM_pow(3,3)

c$$$c     setting mcnlo parameters (needed for amplitudes subroutines)
c$$$      call set_mcnlo_parameters

c     setting madgraph parameters (needed for madgraph subroutines)
      call set_madgraph_parameters


      if(verbose) then
         write(*,*) '--------------------------------------'
         write(*,*) 'POWHEG: RELEVANT PARAMETERS'
         write(*,*) 'top mass       ',topmass_pow
         write(*,*) 'top width      ',topwidth_pow      
         write(*,*) '1/alphaem      ',1.d0/alphaem_pow
         write(*,*) 'W mass         ',wmass_pow
         write(*,*) 'W width        ',wwidth_pow
         write(*,*) 'sin2w          ',sthw2_pow
         write(*,*)'CKM matrix (rows:u,c,t columns:d,s,b )'
      write(*,'(a,3(f10.7))')' ',CKM_pow(1,1),CKM_pow(1,2),CKM_pow(1,3)
      write(*,'(a,3(f10.7))')' ',CKM_pow(2,1),CKM_pow(2,2),CKM_pow(2,3)
      write(*,'(a,3(f10.7))')' ',CKM_pow(3,1),CKM_pow(3,2),CKM_pow(3,3)
         write(*,*) 'lambda_QCD     ',st_lambda5MSB
         write(*,'(1X,A,f7.3,A,f15.7)') 'alpha_s(',91.2d0,')'
     $,pwhg_alphas(91.2d0**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.7)') 'alpha_s(',topmass_pow,')'
     $,pwhg_alphas(topmass_pow**2,st_lambda5MSB,st_nlight)
         write(*,*)
         write(*,*) '(top * primary W) branching ratio ',totbr
         write(*,*) '--------------------------------------'
      endif

      end

c-------------------------------------------------------------------------
      function alfaem(q2)
c Alpha_em(MSbar) at the scale q2 = q^2. 
c Uses alpha_Thomson below the electron mass, alpha(mass) below
c mu_mass and m_tau, and the evolution equation above m_tau, comnsidering the b threshold
c This function is taken from the MC@NLO and modified by SA&ER
c-------------------------------------------------------------------------
      implicit none
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer npoints,ideg
      parameter (npoints=3,ideg=3)
      real*8 ooa(npoints),xlogmu(npoints)
c 1/alpha_em at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data ooa     / 137.036, 135.95, 133.513 /
c logs of sqrt(q2) at m_e=0.000511,m_mu=0.1056,m_tau=1.777      
      data xlogmu  / -7.57914, -2.2481, 0.574927 /
      real *8 zm
      real*8 ooaz,xlq,b,q2
      real *8 alfaem

      real *8 alphaem_inv
      common/calphaem_inv/alphaem_inv

      zm=zmass_pow
      ooaz=alphaem_inv

      if(q2.lt.exp(2.*xlogmu(1))) then
         alfaem = 1.d0/ooa(1)	 
      elseif(q2.lt.exp(2.*xlogmu(2))) then
         xlq = log(q2)/2.d0
         alfaem = 1.d0/ooa(2)
      elseif(q2.lt.exp(2.*xlogmu(3))) then
         xlq = log(q2)/2.d0
         alfaem = 1.d0/ooa(3)
      elseif(q2.lt.5.**2) then
         b = 3 + 2*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2) - 2.*xlogmu(3)
         alfaem = 1d0/ooa(3)/(1.d0 - 1.d0/3.d0/pi/ooa(3)*b*xlq)
      else
         b = 3 + 3*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2/zm**2)
         alfaem = 1d0/ooaz/(1.d0 - 1.d0/3.d0/pi/ooaz*b*xlq)
      endif
      return
      end


c     setting of MADGRAPH inputs
      subroutine set_madgraph_parameters
      include 'PhysPars.h'
      include 'pwhg_math.h'

cccccccccccccccccccccccccccccccc    
c     common bl. originally present in lh_readin, needed
c     by my_setpara
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
ccccccccccccccccccccccccccccccccc
c     common bl. originally present in setpar.f, needed
c     by HELAS subroutine
C
C     BEAM POLARIZATION
C
      REAL*8 POL(2)
      common/to_polarization/ POL
      data POL/1d0,1d0/
ccccccccccccccccccccccccccccccc
     
      real *8 www

      write(*,*) 'POWHEG: set_madgraph_parameters called'

c     SM INPUTS
      alpha=alphaem_pow         !Mad
      zmass=zmass_pow           !Mad
c     alfas=alphas_pow          !can change event-by-event

c     YUKAWA (unuseful in single top)
      mcMS=0d0                  !Mad
      mbMS=0d0                  !Mad
c     mtMS=topmass_pow          !can change event-by-event

c     CKM
c     1st row
      Vud=CKM_pow(1,1)
      Vus=CKM_pow(1,2)
      Vub=CKM_pow(1,3)
c     2nd row
      Vcd=CKM_pow(2,1)
      Vcs=CKM_pow(2,2)
      Vcb=CKM_pow(2,3)
c     3rd row
      Vtd=CKM_pow(3,1)
      Vts=CKM_pow(3,2)
      Vtb=CKM_pow(3,3)

c     MASSES
      bmass=0d0                 !Mad
c     tmass=topmass_pow         !can change event-by-event

      wmass=wmass_pow           !Mad 

c     WIDTHS
c     twidth=0d0                !can change event-by-event
c     wwidth=0d0                !can change event-by-event

c     Setting of wm MadGraph parameter. This is used only to
c     calculate the g_w (weak coupling) used in HELAS subroutines.
c     To have the same coupling of POWHEG, the following ad-hoc definition
c     of Madgraph gfermi is mandatory.
c     The following is an inversion of the assignment formula
c     for wm (see my_setpara subroutine).
      www=zmass*sqrt(1-sthw2_pow)
      gfermi=pi*zmass**2*alpha/sqrt(2.)
      gfermi=gfermi/(zmass**2*www**2 - www**4)

c     setting of other remaining couplings is done by my_setpara on
c     an event by event basis
      end


ccccccccccccccccccccccccccccccccccccccccccccc
c     !: beginning of interface subroutines to madgraph     
ccccccccccccccccccccccccccccccccccccccccccccc
c     see subroutine my_setpara




c$$$c     setting of MCNLO inputs
c$$$      subroutine set_mcnlo_parameters
c$$$      include 'PhysPars.h'
c$$$      include 'stpcblks.h'
c$$$
c$$$      integer idrmax(1:3,3)
c$$$      common/cidrmax/idrmax
c$$$      integer i,j
c$$$
c$$$      write(*,*) 'POWHEG: set_mcnlo_parameters called'
c$$$
c$$$c Fills the array idrmax(jproc,ich) of maximum values for idr; depends 
c$$$c on process type and channel
c$$$      do i=2,3
c$$$        do j=1,2
c$$$          idrmax(i,j)=4
c$$$        enddo
c$$$      enddo
c$$$      idrmax(3,2)=8
c$$$      idrmax(1,3)=1
c$$$      idrmax(2,3)=7
c$$$      idrmax(3,3)=3
c$$$
c$$$      nl=5
c$$$
c$$$c     MCNLO fixed parameters.
c$$$c     In this case, both xm1 and xm2 could be set here, because
c$$$c     I will note change the values of the top and the w mass
c$$$c     during the run. However, to be as much similar as possible
c$$$c     to the t_ch code, I set these parameters before the
c$$$c     effective calls to mcnlo subroutines.
c$$$c$$$      xm1=topmass_pow        
c$$$c$$$      xm12=xm1**2   
c$$$c$$$
c$$$c$$$      xmw = wmass_pow
c$$$c$$$      xm2 = xmw
c$$$c$$$      xm22 = xm2**2
c$$$c$$$      xmW2 = xm22
c$$$
c$$$c     MCNLO parameters that can change during the run have to be assigned
c$$$c     in the proper place and not here.
c$$$c     In fact this subroutine is called only at the beginning.
c$$$c$$$         xmuf2h1 =st_mufact2
c$$$c$$$         xmuf2h2 =st_mufact2
c$$$c$$$         xmur2   =st_muren2
c$$$      end
