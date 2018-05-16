      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include '../include/pwhg_st.h'
      include '../include/pwhg_math.h'
      include '../include/pwhg_flg.h'
      integer i,j
      double precision CKM(1:6,1:6)
      common/cckm/CKM
      double precision masswindow_low,masswindow_high
      double precision powheginput,pwhg_alphas
      external powheginput,pwhg_alphas
      double precision alphaem_inv
      logical verbose
      parameter(verbose=.true.)
cccccccccccccccccccccc
c     needed to avoid numerical problems
c     in madgraph, when we go too close to the
c     singular limits. 
c     In particular, a problem in fvixxx was noticed,
c     in a particularly nasty point:
c     pg = (0.00043528840317960796   ,
c           1.6462894587476746e-08   ,
c           -5.3600767545535054e-09  ,
c           0.00043528840283528763 )
c     and the incoming quark has
c     pq = (6959.0131978776772, 0, 0, 6959.0131978776772)
c     (gluon being very soft/collinear to the z axis, 
c     with incoming emitting parton with essentially 7 TeV energy).
c     In fvixxx, one has
c     pf =  (-6959.012762589181, 1.646289504719789e-08, -5.3600759386540631e-09, -6959.012762589181)
c     and therefore it happens that
c     d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )
c     is not defined (pf2 rounded to exact 0).
c     With (isrcsi,isry,fsrcsi,fsry)=(1d-5,1d-6,1d-5,1d-6),
c     problem still present.
c     The following choice is therefore a compromise to avoid to end up
c     in NaN from madgraph routines because of roundings,
c     and it was found by doing some tests.
c     The drawback is that, when checking limits, some of the collinear
c     limits will not look like very accurate.
c     However the check with MCFM distributions and totals is OK.
      include '../include/pwhg_par.h'
      par_isrtinycsi = 1d-5
      par_isrtinyy =   1d-5
      par_fsrtinycsi = 1d-5
      par_fsrtinyy =   1d-5
c$$$      par_isrtinycsi = 1.d-8
c$$$      par_isrtinyy   = 1.d-8
c$$$      par_fsrtinycsi = 1.d-8
c$$$      par_fsrtinyy   = 1.d-8
cccccccccccccccccccccc

c     to do dedicated imp. sampling in collinear remnants
c     flg_collremnsamp=.true.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c     Do this if we are doing the computation in the 5f scheme
c$$$      print*, '********************************************'
c$$$      print*, '********************************************'
c$$$      print*, '***********      WARNING       *************'
c$$$      print*, '5-flavour running of alfa_s used'
c$$$      print*, 'FONLL correction factors included'
c$$$c     number of light flavors
c$$$      st_nlight = 5
c$$$c     putting this to false is potentially dangerous
c$$$      flg_lightpart_check=.false.
c$$$      print*, '***********      WARNING       *************'
c$$$      print*, '********************************************'
c$$$      print*, '********************************************'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Do this if we are doing the computation in the 4f scheme (4f PDF)
      print*, '********************************************'
      print*, '********************************************'
      print*, '***********      WARNING       *************'
      print*, '4-flavour running of alfa_s used'
      print*, 'No need for FONLL correction factors'
c     number of light flavors
      st_nlight = 4
c     putting this to false is potentially dangerous
      flg_lightpart_check=.true.
      print*, '***********      WARNING       *************'
      print*, '********************************************'
      print*, '********************************************'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     setting physical parameters
      write(*,*) 'POWHEG: loading-setting physical parameters'

c     top mass
      topmass_pow=powheginput('#topmass')
      if(topmass_pow.lt.0) then
         topmass_pow=175.
         write(*,*) 'Using default top mass: ',topmass_pow
      endif
c     top width
      topwidth_pow=0.

c     b mass
      bmass_pow=powheginput('#bmass')
      if(bmass_pow.lt.0) then
         bmass_pow=5.
         write(*,*) 'Using default b mass: ',bmass_pow
      endif

c     ew parameters
c     true inputs are wmass, sthw2, alphaem

c     W mass
      ph_Wmass=powheginput('#wmass')
      if(ph_Wmass.lt.0) then
         ph_Wmass=80.4
         write(*,*) 'Using default W mass: ',ph_Wmass
      endif
c     W width
      ph_Wwidth=0.

c     sthw2
      ph_sthw2=powheginput('#sthw2')
      if(ph_sthw2.lt.0) then
         ph_sthw2=0.23113
         write(*,*) 'Using default sthw2: ',ph_sthw2
      endif
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)

c     alphaem
      alphaem_inv=powheginput('#alphaem_inv')
      if(alphaem_inv.lt.0) then
         alphaem_inv=127.011989
         write(*,*) 'Using default 1/alphaem: ',alphaem_inv
      endif
      ph_alphaem=1d0/alphaem_inv

c     not really needed
      ph_Zmass=90.
      ph_Zwidth=2.

c     CKM matrix entries
      ph_CKM(1,1)= powheginput('#CKM_Vud')
      ph_CKM(1,2)= powheginput('#CKM_Vus') 
      ph_CKM(1,3)= powheginput('#CKM_Vub') 
      ph_CKM(2,1)= powheginput('#CKM_Vcd') 
      ph_CKM(2,2)= powheginput('#CKM_Vcs') 
      ph_CKM(2,3)= powheginput('#CKM_Vcb') 
      ph_CKM(3,1)= powheginput('#CKM_Vtd') 
      ph_CKM(3,2)= powheginput('#CKM_Vts') 
      ph_CKM(3,3)= powheginput('#CKM_Vtb') 

      if(ph_CKM(1,1).lt.0) ph_CKM(1,1)= 0.9740
      if(ph_CKM(1,2).lt.0) ph_CKM(1,2)= 0.2225
      if(ph_CKM(1,3).lt.0) ph_CKM(1,3)= 1d-6

      if(ph_CKM(2,1).lt.0) ph_CKM(2,1)= 0.2225
      if(ph_CKM(2,2).lt.0) ph_CKM(2,2)= 0.9740
      if(ph_CKM(2,3).lt.0) ph_CKM(2,3)= 1d-6

      if(ph_CKM(3,1).lt.0) ph_CKM(3,1)= 1d-6
      if(ph_CKM(3,2).lt.0) ph_CKM(3,2)= 1d-6
      if(ph_CKM(3,3).lt.0) ph_CKM(3,3)= 1d0

ccccccccccccccccccccccccccccccccc
c     useful to have CKM element avoiding
c     gymnastic with indexes
      do i=1,6
         do j=1,6
            CKM(i,j)=0d0
         enddo
      enddo
      CKM(2,1)=ph_CKM(1,1)
      CKM(1,2)=ph_CKM(1,1)
      CKM(2,3)=ph_CKM(1,2)
      CKM(3,2)=ph_CKM(1,2)
      CKM(2,5)=ph_CKM(1,3)
      CKM(5,2)=ph_CKM(1,3)

      CKM(4,1)=ph_CKM(2,1)
      CKM(1,4)=ph_CKM(2,1)
      CKM(4,3)=ph_CKM(2,2)
      CKM(3,4)=ph_CKM(2,2)
      CKM(4,5)=ph_CKM(2,3)
      CKM(5,4)=ph_CKM(2,3)

      CKM(6,1)=ph_CKM(3,1)
      CKM(1,6)=ph_CKM(3,1)
      CKM(6,3)=ph_CKM(3,2)
      CKM(3,6)=ph_CKM(3,2)
      CKM(6,5)=ph_CKM(3,3)
      CKM(5,6)=ph_CKM(3,3)
ccccccccccccccccccccccccccccccccccc



      if(verbose) then
         write(*,*) '--------------------------------------'
         write(*,*) 'POWHEG: RELEVANT PARAMETERS'
         write(*,*) 'top mass       ',topmass_pow
         write(*,*) 'top width      ',topwidth_pow
         write(*,*) 'b mass         ',bmass_pow
         write(*,*) 'W mass         ',ph_Wmass
         write(*,*) 'W width        ',ph_Wwidth
         write(*,*) '1/alphaem      ',1.d0/ph_alphaem
         write(*,*) 'sin2w          ',ph_sthw2
         write(*,*)'CKM matrix (rows:u,c,t columns:d,s,b )'
      write(*,'(a,3(f10.7))')' ',ph_CKM(1,1),ph_CKM(1,2),ph_CKM(1,3)
      write(*,'(a,3(f10.7))')' ',ph_CKM(2,1),ph_CKM(2,2),ph_CKM(2,3)
      write(*,'(a,3(f10.7))')' ',ph_CKM(3,1),ph_CKM(3,2),ph_CKM(3,3)
         write(*,*) 'lambda5MSB_QCD     ',st_lambda5MSB
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',91.1876d0,')'
     $,pwhg_alphas(91.1876d0**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',topmass_pow,')'
     $,pwhg_alphas(topmass_pow**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',topmass_pow/4d0,')'
     $,pwhg_alphas((topmass_pow/4d0)**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',4.75d0,')'
     $,pwhg_alphas((4.75d0)**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',3.00d0,')'
     $,pwhg_alphas((3.00d0)**2,st_lambda5MSB,st_nlight)
         write(*,'(1X,A,f7.3,A,f15.10)') 'alpha_s(',1d0,')'
     $,pwhg_alphas(1d0,st_lambda5MSB,st_nlight)
         write(*,*)
c         write(*,*) 'top branching ratio ',totbr
         write(*,*) '--------------------------------------'
      endif

cccccccccccccc
c     This can be safely deleted
c      st_nlight=5
c      flg_lightpart_check=.false.
c$$$         write(*,*) 'alpha_s(',1d0,')'
c$$$     $,pwhg_alphas(1.3d0,st_lambda5MSB,-1),st_lambda5MSB
cccccccccccccc

      end



