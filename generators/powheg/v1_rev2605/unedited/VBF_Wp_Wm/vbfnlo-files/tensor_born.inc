c  precalculated lepton/quark tensors for WW --> 4 fermions 
c  with Born kinematics
c
c  lepton/quark spinors and W+- polarization vectors
      double complex wep_born, wve_born, wmu_born, wvm_born, 
     &               wp_born, wm_born, wpp_born, wmp_born, 
     &               xp_born, xm_born
      double precision qp_born, qm_born
      common /wpmspinor_born/ wep_born(6), wve_born(6), 
     &                   wmu_born(6), wvm_born(6),
     &                   wp_born(6), wm_born(6), 
     &                   qp_born(0:4), qm_born(0:4),
     &                   wpp_born(6), wmp_born(6), 
     &                   xp_born, xm_born 
c  fermionic tensors
      double complex aww_born, zww_born, aaww_born, azww_born, 
     &           zaww_born, zzww_born, 
     &           NCwpa_born,NCwpz_born,CCwpa_born,CCwpz_born,wwww5_born, 
     &           NCwma_born,NCwmz_born,CCwma_born,CCwmz_born,wwww6_born
      common /tensww_born/ aww_born(0:5), zww_born(0:5), 
     &                aaww_born(0:3,0:3,3), 
     &                azww_born(0:3,0:3,3),  
     &                zaww_born(0:3,0:3,3), 
     &                zzww_born(0:3,0:3,3),  
     &                wwww5_born(0:3,0:3,3), 
     &                wwww6_born(0:3,0:3,3),
     &                NCwpa_born(0:3,0:3,2,3),  
     &                NCwpz_born(0:3,0:3,2,3), 
     &                CCwpa_born(0:3,0:3,2,3),  
     &                CCwpz_born(0:3,0:3,2,3),
     &                NCwma_born(0:3,0:3,2,3),  
     &                NCwmz_born(0:3,0:3,2,3), 
     &                CCwma_born(0:3,0:3,2,3),  
     &                CCwmz_born(0:3,0:3,2,3)
