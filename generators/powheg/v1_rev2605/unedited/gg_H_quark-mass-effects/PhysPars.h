c -*- Fortran -*-

c     Simulation parameters used in the calculations - Set to OS/DR(MS)bar from time to time according to
c     the value of massaren
      real * 8 ph_alphaem,ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth,ph_cthw,
     $     ph_cthw2,ph_sthw,ph_sthw2,ph_Hmass,ph_Hmass2, ph_Hwidth,
     $     ph_Hmass2low,ph_Hmass2high,ph_HmHw,ph_unit_e,ph_alphaemmz,
     $     ph_CKM(3,3),ph_GF,ph_topmass,ph_bottommass,ph_charmmass,
     $     ph_bottommass_mh,ph_cb_mh,ph_tanb_mh,ph_tanb_mz,ph_ma_mz,
     $     ph_MSUSY,ph_mbmb,ph_asmz_nnlo
      common/ph_common/ph_alphaem,ph_Zmass,ph_Zwidth,ph_Wmass,
     $     ph_Wwidth,ph_cthw,ph_cthw2,ph_sthw,ph_sthw2,
     $     ph_Hmass,ph_Hmass2,ph_Hwidth,ph_Hmass2low,ph_alphaemmz,
     $     ph_Hmass2high, ph_HmHw,
     $     ph_unit_e,ph_CKM,ph_GF,ph_topmass,ph_bottommass,
     $     ph_charmmass,ph_bottommass_mh,ph_tanb_mh,ph_cb_mh,
     $     ph_tanb_mz,ph_ma_mz,ph_MSUSY,ph_mbmb,ph_asmz_nnlo

c     Simulation parameters used in the calculations
      real * 8 ph_t1, ph_t2, ph_t1_2,ph_t2_2, ph_mbsb,
     $     ph_b1,ph_b2,ph_b1_2,ph_b2_2,ph_cb,
     $     ph_sb,ph_tanb,ph_ctt,ph_stt,ph_ctb,ph_stb,
     $     ph_alpha,q,ph_mumssm,ph_mg
      common /ph_mssm_common/ph_t1,ph_t2,ph_t1_2,ph_t2_2,
     $     ph_mbsb,ph_b1,ph_b2,ph_b1_2,ph_b2_2,ph_cb,
     $     ph_sb,ph_tanb,ph_ctt,ph_stt,ph_ctb,ph_stb,
     $     ph_alpha,q,ph_mumssm,ph_mg


c     Array with EW corrections for fast evaluation
      real * 8 cached_ew_corr(0:10000)
      integer ewbins
      common /ewcorrs/cached_ew_corr,ewbins

c     Only OS parameters
      real * 8 ph_mbpole
      common /ph_common/ph_mbpole

c     Maximumx numbers of fermions and scalars running in the loops
      integer nfer,nsca
      parameter(nfer=3)
      parameter(nsca=4)

c     Scalar renormalization scale
      real * 8 mumass2
      common /renfacscale/mumass2

      real * 8 trfer(nfer), trsca(nsca)
      common /representations/trfer, trsca

      real * 8 lambdafer(nfer), lambdasca(nsca)
      common /couplings/lambdafer, lambdasca

c     NNLO rescaling factor
      real * 8 nnlorescfactor
      common /nnlo/nnlorescfactor

c     Input customizable by the user
      real * 8 mfer(nfer), msca(nsca), mmaa
      common /particlesmasses/mmaa, mfer, msca

c     the number of active fermion and scalar fields running in the loops
      integer afer,asca
      common /active/afer,asca


c     log(mfer**2/q**2) used in the DRbar/MSbar calculation
      real*8 ferlogmratio(nfer)
      common /massratio/ferlogmratio

c     CP-even, neutral higgs type for MSSM
c     1 = light Higgs, 2 = heavy Higgs
      integer ih
      common /mssmhiggssel/ih
      
      real * 8 mh2,mh
      common /higgsvirtuality/mh2,mh
