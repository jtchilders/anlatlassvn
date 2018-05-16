c -*- Fortran -*-

#include "SLHA.h"

c     SM parameters needed for slepton pair production
      real *8 ph_alphaem0,ph_dalphaemmz,
     &     ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth,ph_sthw,ph_cthw,ph_gf,
     &     ph_Zmass2,ph_Wmass2,ph_sthw2
      logical ph_runningemscale
      common/ph_common/
     &     ph_alphaem0,        ! fine structure constant
     &     ph_dalphaemmz,      ! universal Delta(alpha_em)(Mz^2) in ONS scheme
     &     ph_Zmass,ph_Zwidth, ! mass and width of Z boson
     &     ph_Wmass,ph_Wwidth, ! mass and width of W boson
     &     ph_sthw,ph_cthw,    ! sin(theta_w),cos(theta_w)
     &     ph_gf,              ! Fermi constant
     &     ph_Zmass2,ph_Wmass2,ph_sthw2, ! squared parameters
     &     ph_runningemscale   ! whether to use a running alpha_em

c     MSSM parameters needed for slepton pair production
c     slepton 3: momentum p3, negative charge, mass ph_slepton3mass [GeV]
c     slepton 4: momentum p4, positive charge, mass ph_slepton4mass [GeV]
c     Let Slm = (Sl1, Sl2) be mass eigenstates, Slg = (Sll, Slr) the
c     partners of the left and right handed lepton and Slm_i = U_ij Slg_j
c     Then we denote by U3l U_{type3,1} where type3=1 or 2 is the mass
c     eigenstate type of slepton 3. We use the described order of
c     the indices also for the Fortran "matrices".
      complex *16 ph_c_U3l,ph_c_U3r,ph_c_U4l,ph_c_U4r
      real *8     ph_slepton3mass,ph_slepton4mass
      common /ph_process/ ph_c_U3l,ph_c_U3r,ph_c_U4l,ph_c_U4r,
     &        ph_slepton3mass,ph_slepton4mass

c     MSSM parameters in SLHALib
      complex *16 slhadata(nslhadata)
      common /ph_slha/ slhadata

