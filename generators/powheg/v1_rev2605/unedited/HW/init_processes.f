      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include 'pwhg_flg.h'
      include "coupl.inc"
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode
      integer i
      real * 8 powheginput
      external powheginput


c     in order to prevent a WRONG call to setlocalscales, set
c     MiNLO flags to false
      flg_minlo=.false.
      
      idvecbos=powheginput('idvecbos')

      call init_processes_born
      call init_processes_real

      st_nlight=5

      call init_couplings

      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo
      
      end



      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode
      integer jborn

      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           2
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=         -11
      flst_born(   5,   1)=          12

      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           4
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=         -11
      flst_born(   5,   2)=          12

      flst_born(   1,   3)=           2
      flst_born(   2,   3)=          -1
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=         -11
      flst_born(   5,   3)=          12

      flst_born(   1,   4)=           2
      flst_born(   2,   4)=          -3
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=         -11
      flst_born(   5,   4)=          12

      flst_born(   1,   5)=           2
      flst_born(   2,   5)=          -5
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=         -11
      flst_born(   5,   5)=          12

      flst_born(   1,   6)=           4
      flst_born(   2,   6)=          -1
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=         -11
      flst_born(   5,   6)=          12

      flst_born(   1,   7)=           4
      flst_born(   2,   7)=          -3
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=         -11
      flst_born(   5,   7)=          12

      flst_born(   1,   8)=           4
      flst_born(   2,   8)=          -5
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=         -11
      flst_born(   5,   8)=          12

      flst_born(   1,   9)=          -3
      flst_born(   2,   9)=           2
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=         -11
      flst_born(   5,   9)=          12

      flst_born(   1,  10)=          -3
      flst_born(   2,  10)=           4
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=         -11
      flst_born(   5,  10)=          12

      flst_born(   1,  11)=          -5
      flst_born(   2,  11)=           2
      flst_born(   3,  11)=          25
      flst_born(   4,  11)=         -11
      flst_born(   5,  11)=          12

      flst_born(   1,  12)=          -5
      flst_born(   2,  12)=           4
      flst_born(   3,  12)=          25
      flst_born(   4,  12)=         -11
      flst_born(   5,  12)=          12

      flst_nborn=          12

      if(idvecbos.eq.-24) then
         do jborn=1,flst_nborn
            call cconj(flst_born(:,jborn),nlegborn)
         enddo
      endif

      end

      subroutine cconj(flav,nflav)
      implicit none
      integer nflav,flav(nflav)
      integer k
      do k=1,nflav
         if(flav(k).ne.25) then
            flav(k)=-flav(k)
         endif
      enddo
      end


      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode
      integer jreal

      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=           2
      flst_real(   3,   1)=          25
      flst_real(   4,   1)=         -11
      flst_real(   5,   1)=          12
      flst_real(   6,   1)=           0

      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           4
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=         -11
      flst_real(   5,   2)=          12
      flst_real(   6,   2)=           0

      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=           0
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=         -11
      flst_real(   5,   3)=          12
      flst_real(   6,   3)=          -2

      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           0
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=         -11
      flst_real(   5,   4)=          12
      flst_real(   6,   4)=          -4

      flst_real(   1,   5)=           2
      flst_real(   2,   5)=          -1
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=         -11
      flst_real(   5,   5)=          12
      flst_real(   6,   5)=           0

      flst_real(   1,   6)=           2
      flst_real(   2,   6)=          -3
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=         -11
      flst_real(   5,   6)=          12
      flst_real(   6,   6)=           0

      flst_real(   1,   7)=           2
      flst_real(   2,   7)=          -5
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=         -11
      flst_real(   5,   7)=          12
      flst_real(   6,   7)=           0

      flst_real(   1,   8)=           2
      flst_real(   2,   8)=           0
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=         -11
      flst_real(   5,   8)=          12
      flst_real(   6,   8)=           1

      flst_real(   1,   9)=           2
      flst_real(   2,   9)=           0
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=         -11
      flst_real(   5,   9)=          12
      flst_real(   6,   9)=           3

      flst_real(   1,  10)=           2
      flst_real(   2,  10)=           0
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=         -11
      flst_real(   5,  10)=          12
      flst_real(   6,  10)=           5

      flst_real(   1,  11)=           4
      flst_real(   2,  11)=          -1
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=         -11
      flst_real(   5,  11)=          12
      flst_real(   6,  11)=           0

      flst_real(   1,  12)=           4
      flst_real(   2,  12)=          -3
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=         -11
      flst_real(   5,  12)=          12
      flst_real(   6,  12)=           0

      flst_real(   1,  13)=           4
      flst_real(   2,  13)=          -5
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=         -11
      flst_real(   5,  13)=          12
      flst_real(   6,  13)=           0

      flst_real(   1,  14)=           4
      flst_real(   2,  14)=           0
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=         -11
      flst_real(   5,  14)=          12
      flst_real(   6,  14)=           1

      flst_real(   1,  15)=           4
      flst_real(   2,  15)=           0
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=         -11
      flst_real(   5,  15)=          12
      flst_real(   6,  15)=           3

      flst_real(   1,  16)=           4
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=         -11
      flst_real(   5,  16)=          12
      flst_real(   6,  16)=           5

      flst_real(   1,  17)=          -3
      flst_real(   2,  17)=           2
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=         -11
      flst_real(   5,  17)=          12
      flst_real(   6,  17)=           0

      flst_real(   1,  18)=          -3
      flst_real(   2,  18)=           4
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=         -11
      flst_real(   5,  18)=          12
      flst_real(   6,  18)=           0

      flst_real(   1,  19)=          -3
      flst_real(   2,  19)=           0
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=         -11
      flst_real(   5,  19)=          12
      flst_real(   6,  19)=          -2

      flst_real(   1,  20)=          -3
      flst_real(   2,  20)=           0
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=         -11
      flst_real(   5,  20)=          12
      flst_real(   6,  20)=          -4

      flst_real(   1,  21)=          -5
      flst_real(   2,  21)=           2
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=         -11
      flst_real(   5,  21)=          12
      flst_real(   6,  21)=           0

      flst_real(   1,  22)=          -5
      flst_real(   2,  22)=           4
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=         -11
      flst_real(   5,  22)=          12
      flst_real(   6,  22)=           0

      flst_real(   1,  23)=          -5
      flst_real(   2,  23)=           0
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=         -11
      flst_real(   5,  23)=          12
      flst_real(   6,  23)=          -2

      flst_real(   1,  24)=          -5
      flst_real(   2,  24)=           0
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=         -11
      flst_real(   5,  24)=          12
      flst_real(   6,  24)=          -4

      flst_real(   1,  25)=           0
      flst_real(   2,  25)=          -1
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=         -11
      flst_real(   5,  25)=          12
      flst_real(   6,  25)=          -2

      flst_real(   1,  26)=           0
      flst_real(   2,  26)=          -1
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=         -11
      flst_real(   5,  26)=          12
      flst_real(   6,  26)=          -4

      flst_real(   1,  27)=           0
      flst_real(   2,  27)=           2
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=         -11
      flst_real(   5,  27)=          12
      flst_real(   6,  27)=           1

      flst_real(   1,  28)=           0
      flst_real(   2,  28)=           2
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=         -11
      flst_real(   5,  28)=          12
      flst_real(   6,  28)=           3

      flst_real(   1,  29)=           0
      flst_real(   2,  29)=           2
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=         -11
      flst_real(   5,  29)=          12
      flst_real(   6,  29)=           5

      flst_real(   1,  30)=           0
      flst_real(   2,  30)=           4
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=         -11
      flst_real(   5,  30)=          12
      flst_real(   6,  30)=           1

      flst_real(   1,  31)=           0
      flst_real(   2,  31)=           4
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=         -11
      flst_real(   5,  31)=          12
      flst_real(   6,  31)=           3

      flst_real(   1,  32)=           0
      flst_real(   2,  32)=           4
      flst_real(   3,  32)=          25
      flst_real(   4,  32)=         -11
      flst_real(   5,  32)=          12
      flst_real(   6,  32)=           5

      flst_real(   1,  33)=           0
      flst_real(   2,  33)=          -3
      flst_real(   3,  33)=          25
      flst_real(   4,  33)=         -11
      flst_real(   5,  33)=          12
      flst_real(   6,  33)=          -2

      flst_real(   1,  34)=           0
      flst_real(   2,  34)=          -3
      flst_real(   3,  34)=          25
      flst_real(   4,  34)=         -11
      flst_real(   5,  34)=          12
      flst_real(   6,  34)=          -4

      flst_real(   1,  35)=           0
      flst_real(   2,  35)=          -5
      flst_real(   3,  35)=          25
      flst_real(   4,  35)=         -11
      flst_real(   5,  35)=          12
      flst_real(   6,  35)=          -2

      flst_real(   1,  36)=           0
      flst_real(   2,  36)=          -5
      flst_real(   3,  36)=          25
      flst_real(   4,  36)=         -11
      flst_real(   5,  36)=          12
      flst_real(   6,  36)=          -4

      flst_nreal=          36

      if(idvecbos.eq.-24) then
         do jreal=1,flst_nreal
            call cconj(flst_real(:,jreal),nlegreal)
         enddo
      endif

      end

