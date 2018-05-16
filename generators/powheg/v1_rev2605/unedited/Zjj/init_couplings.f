      subroutine init_couplings
      implicit none
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
cccccccccccccccccccccc
      include 'pwhg_par.h'
      include 'pwhg_flg.h'
      include 'pwhg_physpar.h'
cccccccccccccccccccccc
      double precision masswindow_low,masswindow_high
      double precision powheginput
      external powheginput
      logical verbose
      parameter(verbose=.true.)

      character*11 filename
      logical use_OLP_Interface
      data use_OLP_Interface/.false./
      save use_OLP_Interface


cccccccccccccccccccccccccccccccccccccccccc
c     QUADRUPLE (pt5>100, pt6>100, s56>100^2, nflav=2)
c     F77= /usr/local/intel/fc/11.0.74/bin/ia32/ifort \
c     -r16 -double_size 128 -extend-source -g -fpconstant -fpe0 -C -traceback -save

c     DOUBLE
c     Keep using the defaults, but here I list some
c     settings that were used to test the code...
cccccccccccccccccccccccccccccccccccccccccc
c     For very nasty phase space points selected (randomly) when
c     performing the limits check, numerical accuracy can be problematic
c     for some collinear limits.  To do a 'safe' test, use iseed 1234222
c     @ LHC7, and set tinies as follows

c$$$      par_isrtinycsi = 1.d-7
c$$$      par_isrtinyy   = 1.d-7
c$$$      par_fsrtinycsi = 1.d-7
c$$$      par_fsrtinyy   = 1.d-7

c$$$      par_isrtinycsi = 1.d-8
c$$$      par_isrtinyy   = 1.d-8
c$$$      par_fsrtinycsi = 1.d-8
c$$$      par_fsrtinyy   = 1.d-8


ccccccccccccccccc
c     OTHER PARAMETERS THAT CAN BE USEFUL
c$$$  flg_jacsing=.true.
c$$$  par_diexp = 2
c$$$  par_dijexp= 2
c$$$  par_2gsupp= 2

c$$$  flg_sijtheta=.true.
c$$$  flg_doublefsr=.true.
c$$$  flg_doublefsr=.false.

cccccccccccccccccccccccccccccccc
c     masses for reshuffling 
cccccccccccccccccccccccccccccccc

      physpar_ml(1)=0.511d-3
      physpar_ml(2)=0.1057d0
      physpar_ml(3)=1.777d0
c      physpar_mq(1)=0.33d0      ! up
c      physpar_mq(2)=0.33d0      ! down
c      physpar_mq(3)=0.50d0      ! strange
      physpar_mq(4)=1.50d0      ! charm
      physpar_mq(5)=4.80d0      ! bottom


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   INDEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Z mass
      ph_Zmass = powheginput('#Zmass')
      if(ph_Zmass.lt.0d0) ph_Zmass = 91.188d0

c     Z width
      ph_Zwidth = powheginput('#Zwidth')
      if(ph_Zwidth.lt.0d0) ph_Zwidth = 2.486d0

c     sinth_w
      ph_sthw2 = powheginput("#sthw2")
      if (ph_sthw2.le.0d0) ph_sthw2 = 0.2312d0

c     alpha em
      ph_alphaem = powheginput("#alphaem")
      if (ph_alphaem.le.0d0) ph_alphaem = 1d0/128.930d0

c     number of light flavors
      st_nlight = 5

c     lepton pair virtuality
      masswindow_low = powheginput("#masswindow_low")
      if (masswindow_low.le.0d0) masswindow_low=10d0
      masswindow_high = powheginput("#masswindow_high")
      if (masswindow_high.le.0d0) masswindow_high=10d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   DEPENDENT QUANTITIES       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ph_sthw = sqrt(ph_sthw2)
      ph_cthw = sqrt(1-ph_sthw2)
      ph_Zmass2 = ph_Zmass**2

c     set mass windows around Z-mass peak in unit of ph_Zwidth
c     It is used in the generation of the Born phase space
      ph_Zmass2low=max(0d0,ph_Zmass-masswindow_low*ph_Zwidth)
      if (ph_Zmass2low.lt.1d0) then
         write(*,*) '*************************************'
         write(*,*) 'WARNING: Z virtuality cutoff at 1 GeV'
         write(*,*) '         to avoid the photon pole    '
         write(*,*) '*************************************'
         ph_Zmass2low=1d0
      endif
      ph_Zmass2low=ph_Zmass2low**2
      ph_Zmass2high=ph_Zmass+masswindow_high*ph_Zwidth
      ph_Zmass2high=min(kn_sbeams,ph_Zmass2high**2)
      ph_ZmZw = ph_Zmass * ph_Zwidth
      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(verbose) then
      write(*,*) '*************************************'
      write(*,*) 'Z mass = ',ph_Zmass
      write(*,*) 'Z width = ',ph_Zwidth
c$$$      write(*,*) 'W mass = ',ph_Wmass
c$$$      write(*,*) 'W width = ',ph_Wwidth
      write(*,*) '1/alphaem = ',1d0/ph_alphaem
      write(*,*) 'sthw2 = ',ph_sthw2
      write(*,*) '*************************************'
      write(*,*)
      write(*,*) '*************************************'
      write(*,*) sqrt(ph_Zmass2low),'< M_Z <',sqrt(ph_Zmass2high)
      write(*,*) '*************************************'
      endif

c     old version
c$$$      ph_Zmass2low=max(0d0,ph_Zmass-masswindow_low*ph_Zwidth)
c$$$      ph_Zmass2low=ph_Zmass2low**2
c$$$      ph_Zmass2high=(ph_Zmass+masswindow_high*ph_Zwidth)**2
c$$$      ph_ZmZw = ph_Zmass * ph_Zwidth
c$$$      ph_unit_e = sqrt(4*pi*ph_alphaem)

      if(powheginput("#use-OLP-interface").eq.1) then
         use_OLP_Interface=.true.
      endif

      if(use_OLP_interface) then
         write(*,*) 
         write(*,*) ' Creating order file for external OLP '
         write(*,*) 
         filename="order.lh"
         call create_OLP_order(filename)
      endif

      if(powheginput("#check_ref_amp").eq.1) then
         if((ph_Zmass.ne.91.1876d0).or.
     $        (ph_Zwidth.ne.2.49d0).or.
     $        (ph_sthw2.ne.0.23d0).or.
     $        (ph_alphaem.ne.0.007763854599d0)) then
            write(*,*) ' WARNING: check_ref_amp activated '
            write(*,*) '          but physical parameters different'
            write(*,*) '          from reference values'
            write(*,*) ' Zmass_ref should be 91.1876          ',
     $           (ph_Zmass.eq.91.1876d0)
            write(*,*) ' Zwidth_ref should be 2.49            ',
     $           (ph_Zwidth.eq.2.49d0)
            write(*,*) ' sthw2_ref should be 0.23             ',
     $           (ph_sthw2.eq.0.23d0)
            write(*,*) ' alphaem_ref should be 0.007763854599 ',
     $           (ph_alphaem.eq.0.007763854599d0)
            call exit(-1)
         endif
      endif


      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   functions needed to build the contract file for OLP interface
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine create_OLP_order(filename)
      implicit none 
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      character *(*) filename
      character * 80 flavstring
      integer flav(nlegborn)
      integer iun,ios,i,j,strlen
      double precision sin_2th
      call newunit(iun)
      open(unit=iun,status='new',file=filename,iostat=ios)
      if(ios.eq.0) then
         write(*,*) " Writing ",filename
         write(iun,*) "# order file written by POWHEG-BOX"
         write(iun,*) 
c$$$         write(iun,*) "MatrixElementSquareType   CHSummed"
c$$$         write(iun,*) "CorrectionType            QCD"
c$$$         write(iun,*) "IRregularisation          CDR"
c$$$         write(iun,*) "MassiveParticleScheme     OnShell"
c$$$         write(iun,*) "OperationMode             CouplingsStrippedOff"
c$$$         write(iun,*) "#ModelFile                "
c$$$         write(iun,*) "SubdivideSubprocess       no"
c$$$         write(iun,*) "AlphasPower               2"
c$$$         write(iun,*) "AlphaPower                1"

c     squaring sin(2th) = 2 sin(th) cos(th)
         sin_2th=sqrt(4d0*ph_sthw2*(1d0-ph_sthw2))

         write(iun,*) "Z_mass    ",ph_Zmass
         write(iun,*) "Z_width   ",ph_Zwidth
         write(iun,*) "#sin_2th   ",sin_2th
         write(iun,*) "sin_th_2  ",ph_sthw2


         write(iun,*)
         do i=1,flst_nborn
            call from_numbers_to_string(nlegborn,flst_born(1,i)
     $           ,flavstring)
            write(iun,'(a)') flavstring(1:strlen(flavstring))
            call intassign(nlegborn,flst_born(1,i),flav)
            do j=1,nlegborn
               if(flav(j).eq.0) flav(j)=21
            enddo
            write(iun,'(6(a,i3))') " ",flav(1)," ",flav(2),"  -> "
     $           ,flav(3)," ",flav(4)," ",flav(5)," ",flav(6)
            write(iun,*)
         enddo
      else
         write(*,*) " File order.lh already present"
      endif
      close(iun)
      end
      
      subroutine from_numbers_to_string(n,flav,string)
      implicit none
      integer n,flav(n)
      include 'nlegborn.h'
      character * (*) string
      integer max_partnames
      parameter (max_partnames=16)
      character * 5 partnames(-max_partnames:max_partnames)
      data partnames/'vtbar','tau+','vmbar','mu+','vebar','e+' ,' ',' '
     $     ,' ',' ','tbar','bbar','cbar' ,'sbar','ubar','dbar','g','d'
     $     ,'u','s' ,'c','b','t' ,' ',' ',' ',' ','e-' ,'ve','mu-','vmu'
     $     ,'tau-','vtau'/
      integer j,nsp
      parameter (nsp=5)
      if(len(string).lt.nsp*(n+1)+5) then
         write(*,*)'from_numbers_to_string: string too short;'
         write(*,*)'Increase its size'
         call exit(-1)
      endif
      string='# '
      do j=1,n
         if (abs(flav(j)).le.max_partnames) then
            string(nsp*j:)=partnames(flav(j))
         else 
            string(nsp*j:)=' ??? '
         endif
      enddo
      do j=len(string)-5,2*nsp+1,-1
         string(j+5:j+5)=string(j:j)
      enddo
      string(2*nsp+3:2*nsp+8)=' ==> '
      end

      integer function strlen(st)
      integer i
      character	st*(*)
      i = len(st)
      do while (st(i:i) .eq. ' ')
        i = i - 1
      enddo
      strlen = i
      return
      end      

