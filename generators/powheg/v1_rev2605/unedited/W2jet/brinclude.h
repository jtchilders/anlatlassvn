c -*-Fortran-*-

      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer  br_nlegborn, br_nlegreal
      parameter (br_nlegborn=nlegborn-1,br_nlegreal=nlegreal-1)

      real * 8 brkn_pborn(0:3,br_nlegborn),brkn_preal(0:3,br_nlegreal),
     1   brkn_cmpborn(0:3,br_nlegborn),brkn_cmpreal(0:3,br_nlegreal),
     1   brkn_sborn,brkn_sreal,
     1   brkn_y,brkn_azi,brkn_csitilde,
     1   brkn_csimax,
     1   brkn_csi,
     1   brkn_jacborn,brkn_jacreal,brkn_xb1,brkn_xb2,
     1   brkn_x1,brkn_x2,brkn_masses(br_nlegreal),
     1   brkn_dijterm(0:br_nlegreal,br_nlegreal),
     1   brkn_ktmin,
     1   brkn_csimax_arr(0:br_nlegborn),brpar_diexp,brpar_dijexp

      integer brkn_emitter

      common/pwhg_brkn/brkn_pborn,brkn_cmpborn,brkn_preal,brkn_cmpreal,
     1   brkn_y,brkn_sborn,brkn_sreal,brkn_azi,brkn_csitilde,
     1   brkn_csimax,brkn_csi,brkn_jacborn,
     1   brkn_jacreal,brkn_xb1,brkn_xb2,brkn_x1,
     1   brkn_x2,brkn_masses,brkn_dijterm,brkn_ktmin,brkn_csimax_arr,
     1   brpar_diexp,brpar_dijexp,
c integers
     1   brkn_emitter

      save /pwhg_brkn/
      
