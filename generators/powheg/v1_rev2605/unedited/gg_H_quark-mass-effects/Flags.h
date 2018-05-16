c -*- Fortran -*-
      integer flg_onshell,flg_ew,flg_noscalars,flg_lhscalars,
     $flg_fast_ew,flg_nnlo,flg_passarino,flg_mssm_q_mh,flg_hdecay
      common /gghflags/flg_onshell,flg_ew,flg_noscalars,flg_lhscalars
     $,flg_fast_ew,flg_nnlo,flg_passarino,flg_mssm_q_mh,flg_hdecay

c     Integer to choose between 0-OS/1-MS(DR)bar/2-Mixed
      integer massren
      common /renscheme/massren

c     Model parameter (0 = SM; 1 = MW; 2 = MSSM)
      integer model
      common /model/model
