c  precalculated lepton/quark tensors for WW --> 4 fermions 
c  with Real kinematics
c
c  lepton/quark spinors and W+- polarization vectors
      double complex wep_real, wve_real, wmu_real, wvm_real, 
     &               wp_real, wm_real, wpp_real, wmp_real, 
     &               xp_real, xm_real
      double precision qp_real, qm_real
      common /wpmspinor_real/ wep_real(6), wve_real(6), 
     &                   wmu_real(6), wvm_real(6),
     &                   wp_real(6), wm_real(6), 
     &                   qp_real(0:4), qm_real(0:4),
     &                   wpp_real(6), wmp_real(6), 
     &                   xp_real, xm_real 
c  fermionic tensors
      double complex aww_real, zww_real, aaww_real, azww_real, 
     &           zaww_real, zzww_real, 
     &           NCwpa_real,NCwpz_real,CCwpa_real,CCwpz_real,wwww5_real, 
     &           NCwma_real,NCwmz_real,CCwma_real,CCwmz_real,wwww6_real
      common /tensww_real/ aww_real(0:5), zww_real(0:5), 
     &                aaww_real(0:3,0:3,3), 
     &                azww_real(0:3,0:3,3),  
     &                zaww_real(0:3,0:3,3), 
     &                zzww_real(0:3,0:3,3),  
     &                wwww5_real(0:3,0:3,3), 
     &                wwww6_real(0:3,0:3,3),
     &                NCwpa_real(0:3,0:3,2,3),  
     &                NCwpz_real(0:3,0:3,2,3), 
     &                CCwpa_real(0:3,0:3,2,3),  
     &                CCwpz_real(0:3,0:3,2,3),
     &                NCwma_real(0:3,0:3,2,3),  
     &                NCwmz_real(0:3,0:3,2,3), 
     &                CCwma_real(0:3,0:3,2,3),  
     &                CCwmz_real(0:3,0:3,2,3)
