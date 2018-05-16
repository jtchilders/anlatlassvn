c -*- Fortran -*-
      integer flg_onshell,flg_ew,flg_noscalars,flg_lhscalars,
     $flg_fast_ew,flg_nnlo,flg_passarino,flg_mssm_q_mh,flg_hdecay,
     $flg_fhdecay,flg_hHweightevents,flg_2HDMtype
      common /gghflags/flg_onshell,flg_ew,flg_noscalars,flg_lhscalars
     $,flg_fast_ew,flg_nnlo,flg_passarino,flg_mssm_q_mh,flg_hdecay,
     $flg_fhdecay,flg_hHweightevents,flg_2HDMtype

c     Integer to choose between 0-OS/1-MS(DR)bar/2-Mixed
c	  massrenb is for the OS=2 subcase for the bottom subroutine
      integer massren,massrenb
      common /renscheme/massren,massrenb
