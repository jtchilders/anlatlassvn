      function pt2max_regular()
      real * 8 pt2max_regular
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'PhysPars.h'
      logical debug
      parameter (debug=.false.) 
      pt2max_regular=sqrt(kn_cmpreal(1,4)**2+kn_cmpreal(2,4)**2)
      end
