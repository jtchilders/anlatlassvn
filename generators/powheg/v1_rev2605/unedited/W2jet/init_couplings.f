      subroutine init_couplings
      implicit none
      include 'pwhg_math.h'
      include 'pwhg_par.h'
      include 'pwhg_st.h'
      include "nlegborn.h"
      include "pwhg_kn.h"
      integer nf
      double precision  Two, Four, Rt2
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter(nf=5)
      include "zcouple.f"
      include "nwz.f"
      include "ckm.f"
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud

      double precision asmz,pwhg_alphas,esq
      double precision hmass, wmass, zmass, tmass,   bmass,
     &     hwidth,wwidth,zwidth,twidth
      integer nproc
      common/nproc/nproc      
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      data called/.false./
      save called
      real *8 powheginput 
      external powheginput 
      if(called) then
         return
      else
         called=.true.
      endif

C-----Fix basic parameters
      alpha=1/1.32506980d+02
      gfermi = 1.16639000d-05
      alfas = 0.119d0
      zmass = 9.11880000d+01
      mcMS = 0d0
      mbMS = 0d0
      mtMS = 174d0
      mtaMS = 1.777d0
      vud = 1d0
      bmass=0
      tmass=173.5d0
c Not used here
      hmass=125
      hwidth=0
      twidth=0
c

      rt2=sqrt(2d0)
      wmass=sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))

c      zwidth=2.44140351d+00
      zwidth=2.49d0
c      wwidth=2.04759951d+00
      wwidth=2.06d0

C     take Z mas from input, if present
      if (powheginput("#zmass") .gt. 0d0) then
         zmass = powheginput("#zmass")
         wmass=sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      endif
c if we insist that wmass is different:
      if (powheginput("#wmass") .gt. 0d0)  wmass = powheginput("#wmass")
      
      st_nlight=5

      esq=4*pi*alpha

C     now fill MCFM parameters 
      call fillMCFMCommon(esq,gfermi,
     &     hmass, wmass, zmass, tmass,   bmass,
     &     hwidth,wwidth,zwidth,twidth,
     &     st_alpha,st_muren2,st_mufact2)

C     set mur/muf 
      call setscalesbtilde 
      call st_mcfm
      call print_st_mcfm

C---setup specific MCFM variables for this process
      call ckmfill(nwz)

      kn_masses(1)=0
      kn_masses(2)=0
      kn_masses(3)=0
      kn_masses(4)=0
      kn_masses(5)=0
      kn_masses(6)=0

      end


      subroutine st_mcfm
      implicit none
      include 'pwhg_st.h' 
      include 'pwhg_math.h' 
      include 'qcdcouple.f' 
      include 'scale.f'
C     set strong coupling and scale for MCFM files 
      musq = st_muren2
      scale = sqrt(musq)
      as = st_alpha
      gsq = 4d0*pi*as
      ason2pi = as/(2d0*pi)
      ason4pi = as/(4d0*pi)
      end

      subroutine print_st_mcfm
      implicit none
      include 'qcdcouple.f' 
      include 'scale.f'
      write(6,*) 'musq=      ',musq
      write(6,*) 'alpha_s(musq)=',as
      end
