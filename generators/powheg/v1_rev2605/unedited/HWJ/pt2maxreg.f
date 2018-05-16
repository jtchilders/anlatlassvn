      function pt2max_regular()
      implicit none
      real * 8 pt2max_regular
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'PhysPars.h'
      real * 8 pt1sq,pt2sq
      pt1sq = kn_cmpreal(1,6)**2+kn_cmpreal(2,6)**2
      pt2sq = kn_cmpreal(1,7)**2+kn_cmpreal(2,7)**2
      pt2max_regular = max(pt1sq,pt2sq)

      end
