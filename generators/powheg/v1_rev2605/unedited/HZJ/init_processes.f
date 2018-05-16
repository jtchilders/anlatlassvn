      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include "pwhg_flg.h"
      include "coupl.inc"
      integer i
      integer idvecbos,vdecaymode
      common/cvecbos/idvecbos,vdecaymode
      real * 8 powheginput
      external powheginput
 
 
      st_bornorder=1
      if(powheginput("#minlo").eq.1) then
         flg_minlo=.true.
         if(powheginput("#minlo_nnll").eq.1) then
            flg_minlo_nnll=.true.
         else
            flg_minlo_nnll=.false.
         endif
         flg_minlo_real=.false.
      else
         flg_minlo=.false.
         flg_minlo_nnll=.false.
         flg_minlo_real=.false.
      endif

     
      call init_processes_born
      call init_processes_real
      st_nlight=5
      call init_couplings
c$$$      if (tmass.eq.0d0) then
c$$$         st_nlight=6
c$$$      elseif(bmass.eq.0d0) then
c$$$         st_nlight=5
c$$$      elseif(cmass.eq.0d0) then
c$$$         st_nlight=4
c$$$      else
c$$$         st_nlight=3
c$$$      endif
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo
 
      return
      end
 
 
 
      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           1
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=         -11
      flst_born(   5,   1)=          11
      flst_born(   6,   1)=           0
 
      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           0
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=         -11
      flst_born(   5,   2)=          11
      flst_born(   6,   2)=          -1
 
      flst_born(   1,   3)=           1
      flst_born(   2,   3)=          -1
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=         -11
      flst_born(   5,   3)=          11
      flst_born(   6,   3)=           0
 
      flst_born(   1,   4)=           1
      flst_born(   2,   4)=           0
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=         -11
      flst_born(   5,   4)=          11
      flst_born(   6,   4)=           1
 
      flst_born(   1,   5)=          -2
      flst_born(   2,   5)=           2
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=         -11
      flst_born(   5,   5)=          11
      flst_born(   6,   5)=           0
 
      flst_born(   1,   6)=          -2
      flst_born(   2,   6)=           0
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=         -11
      flst_born(   5,   6)=          11
      flst_born(   6,   6)=          -2
 
      flst_born(   1,   7)=           2
      flst_born(   2,   7)=          -2
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=         -11
      flst_born(   5,   7)=          11
      flst_born(   6,   7)=           0
 
      flst_born(   1,   8)=           2
      flst_born(   2,   8)=           0
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=         -11
      flst_born(   5,   8)=          11
      flst_born(   6,   8)=           2
 
      flst_born(   1,   9)=          -4
      flst_born(   2,   9)=           4
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=         -11
      flst_born(   5,   9)=          11
      flst_born(   6,   9)=           0
 
      flst_born(   1,  10)=          -4
      flst_born(   2,  10)=           0
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=         -11
      flst_born(   5,  10)=          11
      flst_born(   6,  10)=          -4
 
      flst_born(   1,  11)=           4
      flst_born(   2,  11)=          -4
      flst_born(   3,  11)=          25
      flst_born(   4,  11)=         -11
      flst_born(   5,  11)=          11
      flst_born(   6,  11)=           0
 
      flst_born(   1,  12)=           4
      flst_born(   2,  12)=           0
      flst_born(   3,  12)=          25
      flst_born(   4,  12)=         -11
      flst_born(   5,  12)=          11
      flst_born(   6,  12)=           4
 
      flst_born(   1,  13)=          -3
      flst_born(   2,  13)=           3
      flst_born(   3,  13)=          25
      flst_born(   4,  13)=         -11
      flst_born(   5,  13)=          11
      flst_born(   6,  13)=           0
 
      flst_born(   1,  14)=          -3
      flst_born(   2,  14)=           0
      flst_born(   3,  14)=          25
      flst_born(   4,  14)=         -11
      flst_born(   5,  14)=          11
      flst_born(   6,  14)=          -3
 
      flst_born(   1,  15)=           3
      flst_born(   2,  15)=          -3
      flst_born(   3,  15)=          25
      flst_born(   4,  15)=         -11
      flst_born(   5,  15)=          11
      flst_born(   6,  15)=           0
 
      flst_born(   1,  16)=           3
      flst_born(   2,  16)=           0
      flst_born(   3,  16)=          25
      flst_born(   4,  16)=         -11
      flst_born(   5,  16)=          11
      flst_born(   6,  16)=           3
 
      flst_born(   1,  17)=          -5
      flst_born(   2,  17)=           5
      flst_born(   3,  17)=          25
      flst_born(   4,  17)=         -11
      flst_born(   5,  17)=          11
      flst_born(   6,  17)=           0
 
      flst_born(   1,  18)=          -5
      flst_born(   2,  18)=           0
      flst_born(   3,  18)=          25
      flst_born(   4,  18)=         -11
      flst_born(   5,  18)=          11
      flst_born(   6,  18)=          -5
 
      flst_born(   1,  19)=           5
      flst_born(   2,  19)=          -5
      flst_born(   3,  19)=          25
      flst_born(   4,  19)=         -11
      flst_born(   5,  19)=          11
      flst_born(   6,  19)=           0
 
      flst_born(   1,  20)=           5
      flst_born(   2,  20)=           0
      flst_born(   3,  20)=          25
      flst_born(   4,  20)=         -11
      flst_born(   5,  20)=          11
      flst_born(   6,  20)=           5
 
      flst_born(   1,  21)=           0
      flst_born(   2,  21)=          -1
      flst_born(   3,  21)=          25
      flst_born(   4,  21)=         -11
      flst_born(   5,  21)=          11
      flst_born(   6,  21)=          -1
 
      flst_born(   1,  22)=           0
      flst_born(   2,  22)=           1
      flst_born(   3,  22)=          25
      flst_born(   4,  22)=         -11
      flst_born(   5,  22)=          11
      flst_born(   6,  22)=           1
 
      flst_born(   1,  23)=           0
      flst_born(   2,  23)=          -2
      flst_born(   3,  23)=          25
      flst_born(   4,  23)=         -11
      flst_born(   5,  23)=          11
      flst_born(   6,  23)=          -2
 
      flst_born(   1,  24)=           0
      flst_born(   2,  24)=           2
      flst_born(   3,  24)=          25
      flst_born(   4,  24)=         -11
      flst_born(   5,  24)=          11
      flst_born(   6,  24)=           2
 
      flst_born(   1,  25)=           0
      flst_born(   2,  25)=          -4
      flst_born(   3,  25)=          25
      flst_born(   4,  25)=         -11
      flst_born(   5,  25)=          11
      flst_born(   6,  25)=          -4
 
      flst_born(   1,  26)=           0
      flst_born(   2,  26)=           4
      flst_born(   3,  26)=          25
      flst_born(   4,  26)=         -11
      flst_born(   5,  26)=          11
      flst_born(   6,  26)=           4
 
      flst_born(   1,  27)=           0
      flst_born(   2,  27)=          -3
      flst_born(   3,  27)=          25
      flst_born(   4,  27)=         -11
      flst_born(   5,  27)=          11
      flst_born(   6,  27)=          -3
 
      flst_born(   1,  28)=           0
      flst_born(   2,  28)=           3
      flst_born(   3,  28)=          25
      flst_born(   4,  28)=         -11
      flst_born(   5,  28)=          11
      flst_born(   6,  28)=           3
 
      flst_born(   1,  29)=           0
      flst_born(   2,  29)=          -5
      flst_born(   3,  29)=          25
      flst_born(   4,  29)=         -11
      flst_born(   5,  29)=          11
      flst_born(   6,  29)=          -5
 
      flst_born(   1,  30)=           0
      flst_born(   2,  30)=           5
      flst_born(   3,  30)=          25
      flst_born(   4,  30)=         -11
      flst_born(   5,  30)=          11
      flst_born(   6,  30)=           5
 
      flst_nborn=          30
 
      return
      end
 
 
 
      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=          -1
      flst_real(   3,   1)=          25
      flst_real(   4,   1)=         -11
      flst_real(   5,   1)=          11
      flst_real(   6,   1)=          -1
      flst_real(   7,   1)=          -1
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           1
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=         -11
      flst_real(   5,   2)=          11
      flst_real(   6,   2)=           1
      flst_real(   7,   2)=          -1
 
      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=           1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=         -11
      flst_real(   5,   3)=          11
      flst_real(   6,   3)=           2
      flst_real(   7,   3)=          -2
 
      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           1
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=         -11
      flst_real(   5,   4)=          11
      flst_real(   6,   4)=           4
      flst_real(   7,   4)=          -4
 
      flst_real(   1,   5)=          -1
      flst_real(   2,   5)=           1
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=         -11
      flst_real(   5,   5)=          11
      flst_real(   6,   5)=           3
      flst_real(   7,   5)=          -3
 
      flst_real(   1,   6)=          -1
      flst_real(   2,   6)=           1
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=         -11
      flst_real(   5,   6)=          11
      flst_real(   6,   6)=           5
      flst_real(   7,   6)=          -5
 
      flst_real(   1,   7)=          -1
      flst_real(   2,   7)=           1
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=         -11
      flst_real(   5,   7)=          11
      flst_real(   6,   7)=           0
      flst_real(   7,   7)=           0
 
      flst_real(   1,   8)=          -1
      flst_real(   2,   8)=          -2
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=         -11
      flst_real(   5,   8)=          11
      flst_real(   6,   8)=          -1
      flst_real(   7,   8)=          -2
 
      flst_real(   1,   9)=          -1
      flst_real(   2,   9)=           2
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=         -11
      flst_real(   5,   9)=          11
      flst_real(   6,   9)=          -1
      flst_real(   7,   9)=           2
 
      flst_real(   1,  10)=          -1
      flst_real(   2,  10)=          -4
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=         -11
      flst_real(   5,  10)=          11
      flst_real(   6,  10)=          -1
      flst_real(   7,  10)=          -4
 
      flst_real(   1,  11)=          -1
      flst_real(   2,  11)=           4
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=         -11
      flst_real(   5,  11)=          11
      flst_real(   6,  11)=          -1
      flst_real(   7,  11)=           4
 
      flst_real(   1,  12)=          -1
      flst_real(   2,  12)=          -3
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=         -11
      flst_real(   5,  12)=          11
      flst_real(   6,  12)=          -1
      flst_real(   7,  12)=          -3
 
      flst_real(   1,  13)=          -1
      flst_real(   2,  13)=           3
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=         -11
      flst_real(   5,  13)=          11
      flst_real(   6,  13)=          -1
      flst_real(   7,  13)=           3
 
      flst_real(   1,  14)=          -1
      flst_real(   2,  14)=          -5
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=         -11
      flst_real(   5,  14)=          11
      flst_real(   6,  14)=          -1
      flst_real(   7,  14)=          -5
 
      flst_real(   1,  15)=          -1
      flst_real(   2,  15)=           5
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=         -11
      flst_real(   5,  15)=          11
      flst_real(   6,  15)=          -1
      flst_real(   7,  15)=           5
 
      flst_real(   1,  16)=          -1
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=         -11
      flst_real(   5,  16)=          11
      flst_real(   6,  16)=          -1
      flst_real(   7,  16)=           0
 
      flst_real(   1,  17)=           1
      flst_real(   2,  17)=          -1
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=         -11
      flst_real(   5,  17)=          11
      flst_real(   6,  17)=           1
      flst_real(   7,  17)=          -1
 
      flst_real(   1,  18)=           1
      flst_real(   2,  18)=          -1
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=         -11
      flst_real(   5,  18)=          11
      flst_real(   6,  18)=           2
      flst_real(   7,  18)=          -2
 
      flst_real(   1,  19)=           1
      flst_real(   2,  19)=          -1
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=         -11
      flst_real(   5,  19)=          11
      flst_real(   6,  19)=           4
      flst_real(   7,  19)=          -4
 
      flst_real(   1,  20)=           1
      flst_real(   2,  20)=          -1
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=         -11
      flst_real(   5,  20)=          11
      flst_real(   6,  20)=           3
      flst_real(   7,  20)=          -3
 
      flst_real(   1,  21)=           1
      flst_real(   2,  21)=          -1
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=         -11
      flst_real(   5,  21)=          11
      flst_real(   6,  21)=           5
      flst_real(   7,  21)=          -5
 
      flst_real(   1,  22)=           1
      flst_real(   2,  22)=          -1
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=         -11
      flst_real(   5,  22)=          11
      flst_real(   6,  22)=           0
      flst_real(   7,  22)=           0
 
      flst_real(   1,  23)=           1
      flst_real(   2,  23)=           1
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=         -11
      flst_real(   5,  23)=          11
      flst_real(   6,  23)=           1
      flst_real(   7,  23)=           1
 
      flst_real(   1,  24)=           1
      flst_real(   2,  24)=          -2
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=         -11
      flst_real(   5,  24)=          11
      flst_real(   6,  24)=           1
      flst_real(   7,  24)=          -2
 
      flst_real(   1,  25)=           1
      flst_real(   2,  25)=           2
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=         -11
      flst_real(   5,  25)=          11
      flst_real(   6,  25)=           1
      flst_real(   7,  25)=           2
 
      flst_real(   1,  26)=           1
      flst_real(   2,  26)=          -4
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=         -11
      flst_real(   5,  26)=          11
      flst_real(   6,  26)=           1
      flst_real(   7,  26)=          -4
 
      flst_real(   1,  27)=           1
      flst_real(   2,  27)=           4
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=         -11
      flst_real(   5,  27)=          11
      flst_real(   6,  27)=           1
      flst_real(   7,  27)=           4
 
      flst_real(   1,  28)=           1
      flst_real(   2,  28)=          -3
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=         -11
      flst_real(   5,  28)=          11
      flst_real(   6,  28)=           1
      flst_real(   7,  28)=          -3
 
      flst_real(   1,  29)=           1
      flst_real(   2,  29)=           3
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=         -11
      flst_real(   5,  29)=          11
      flst_real(   6,  29)=           1
      flst_real(   7,  29)=           3
 
      flst_real(   1,  30)=           1
      flst_real(   2,  30)=          -5
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=         -11
      flst_real(   5,  30)=          11
      flst_real(   6,  30)=           1
      flst_real(   7,  30)=          -5
 
      flst_real(   1,  31)=           1
      flst_real(   2,  31)=           5
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=         -11
      flst_real(   5,  31)=          11
      flst_real(   6,  31)=           1
      flst_real(   7,  31)=           5
 
      flst_real(   1,  32)=           1
      flst_real(   2,  32)=           0
      flst_real(   3,  32)=          25
      flst_real(   4,  32)=         -11
      flst_real(   5,  32)=          11
      flst_real(   6,  32)=           1
      flst_real(   7,  32)=           0
 
      flst_real(   1,  33)=          -2
      flst_real(   2,  33)=          -1
      flst_real(   3,  33)=          25
      flst_real(   4,  33)=         -11
      flst_real(   5,  33)=          11
      flst_real(   6,  33)=          -1
      flst_real(   7,  33)=          -2
 
      flst_real(   1,  34)=          -2
      flst_real(   2,  34)=           1
      flst_real(   3,  34)=          25
      flst_real(   4,  34)=         -11
      flst_real(   5,  34)=          11
      flst_real(   6,  34)=           1
      flst_real(   7,  34)=          -2
 
      flst_real(   1,  35)=          -2
      flst_real(   2,  35)=          -2
      flst_real(   3,  35)=          25
      flst_real(   4,  35)=         -11
      flst_real(   5,  35)=          11
      flst_real(   6,  35)=          -2
      flst_real(   7,  35)=          -2
 
      flst_real(   1,  36)=          -2
      flst_real(   2,  36)=           2
      flst_real(   3,  36)=          25
      flst_real(   4,  36)=         -11
      flst_real(   5,  36)=          11
      flst_real(   6,  36)=           1
      flst_real(   7,  36)=          -1
 
      flst_real(   1,  37)=          -2
      flst_real(   2,  37)=           2
      flst_real(   3,  37)=          25
      flst_real(   4,  37)=         -11
      flst_real(   5,  37)=          11
      flst_real(   6,  37)=           2
      flst_real(   7,  37)=          -2
 
      flst_real(   1,  38)=          -2
      flst_real(   2,  38)=           2
      flst_real(   3,  38)=          25
      flst_real(   4,  38)=         -11
      flst_real(   5,  38)=          11
      flst_real(   6,  38)=           4
      flst_real(   7,  38)=          -4
 
      flst_real(   1,  39)=          -2
      flst_real(   2,  39)=           2
      flst_real(   3,  39)=          25
      flst_real(   4,  39)=         -11
      flst_real(   5,  39)=          11
      flst_real(   6,  39)=           3
      flst_real(   7,  39)=          -3
 
      flst_real(   1,  40)=          -2
      flst_real(   2,  40)=           2
      flst_real(   3,  40)=          25
      flst_real(   4,  40)=         -11
      flst_real(   5,  40)=          11
      flst_real(   6,  40)=           5
      flst_real(   7,  40)=          -5
 
      flst_real(   1,  41)=          -2
      flst_real(   2,  41)=           2
      flst_real(   3,  41)=          25
      flst_real(   4,  41)=         -11
      flst_real(   5,  41)=          11
      flst_real(   6,  41)=           0
      flst_real(   7,  41)=           0
 
      flst_real(   1,  42)=          -2
      flst_real(   2,  42)=          -4
      flst_real(   3,  42)=          25
      flst_real(   4,  42)=         -11
      flst_real(   5,  42)=          11
      flst_real(   6,  42)=          -2
      flst_real(   7,  42)=          -4
 
      flst_real(   1,  43)=          -2
      flst_real(   2,  43)=           4
      flst_real(   3,  43)=          25
      flst_real(   4,  43)=         -11
      flst_real(   5,  43)=          11
      flst_real(   6,  43)=          -2
      flst_real(   7,  43)=           4
 
      flst_real(   1,  44)=          -2
      flst_real(   2,  44)=          -3
      flst_real(   3,  44)=          25
      flst_real(   4,  44)=         -11
      flst_real(   5,  44)=          11
      flst_real(   6,  44)=          -2
      flst_real(   7,  44)=          -3
 
      flst_real(   1,  45)=          -2
      flst_real(   2,  45)=           3
      flst_real(   3,  45)=          25
      flst_real(   4,  45)=         -11
      flst_real(   5,  45)=          11
      flst_real(   6,  45)=          -2
      flst_real(   7,  45)=           3
 
      flst_real(   1,  46)=          -2
      flst_real(   2,  46)=          -5
      flst_real(   3,  46)=          25
      flst_real(   4,  46)=         -11
      flst_real(   5,  46)=          11
      flst_real(   6,  46)=          -2
      flst_real(   7,  46)=          -5
 
      flst_real(   1,  47)=          -2
      flst_real(   2,  47)=           5
      flst_real(   3,  47)=          25
      flst_real(   4,  47)=         -11
      flst_real(   5,  47)=          11
      flst_real(   6,  47)=          -2
      flst_real(   7,  47)=           5
 
      flst_real(   1,  48)=          -2
      flst_real(   2,  48)=           0
      flst_real(   3,  48)=          25
      flst_real(   4,  48)=         -11
      flst_real(   5,  48)=          11
      flst_real(   6,  48)=          -2
      flst_real(   7,  48)=           0
 
      flst_real(   1,  49)=           2
      flst_real(   2,  49)=          -1
      flst_real(   3,  49)=          25
      flst_real(   4,  49)=         -11
      flst_real(   5,  49)=          11
      flst_real(   6,  49)=          -1
      flst_real(   7,  49)=           2
 
      flst_real(   1,  50)=           2
      flst_real(   2,  50)=           1
      flst_real(   3,  50)=          25
      flst_real(   4,  50)=         -11
      flst_real(   5,  50)=          11
      flst_real(   6,  50)=           1
      flst_real(   7,  50)=           2
 
      flst_real(   1,  51)=           2
      flst_real(   2,  51)=          -2
      flst_real(   3,  51)=          25
      flst_real(   4,  51)=         -11
      flst_real(   5,  51)=          11
      flst_real(   6,  51)=           1
      flst_real(   7,  51)=          -1
 
      flst_real(   1,  52)=           2
      flst_real(   2,  52)=          -2
      flst_real(   3,  52)=          25
      flst_real(   4,  52)=         -11
      flst_real(   5,  52)=          11
      flst_real(   6,  52)=           2
      flst_real(   7,  52)=          -2
 
      flst_real(   1,  53)=           2
      flst_real(   2,  53)=          -2
      flst_real(   3,  53)=          25
      flst_real(   4,  53)=         -11
      flst_real(   5,  53)=          11
      flst_real(   6,  53)=           4
      flst_real(   7,  53)=          -4
 
      flst_real(   1,  54)=           2
      flst_real(   2,  54)=          -2
      flst_real(   3,  54)=          25
      flst_real(   4,  54)=         -11
      flst_real(   5,  54)=          11
      flst_real(   6,  54)=           3
      flst_real(   7,  54)=          -3
 
      flst_real(   1,  55)=           2
      flst_real(   2,  55)=          -2
      flst_real(   3,  55)=          25
      flst_real(   4,  55)=         -11
      flst_real(   5,  55)=          11
      flst_real(   6,  55)=           5
      flst_real(   7,  55)=          -5
 
      flst_real(   1,  56)=           2
      flst_real(   2,  56)=          -2
      flst_real(   3,  56)=          25
      flst_real(   4,  56)=         -11
      flst_real(   5,  56)=          11
      flst_real(   6,  56)=           0
      flst_real(   7,  56)=           0
 
      flst_real(   1,  57)=           2
      flst_real(   2,  57)=           2
      flst_real(   3,  57)=          25
      flst_real(   4,  57)=         -11
      flst_real(   5,  57)=          11
      flst_real(   6,  57)=           2
      flst_real(   7,  57)=           2
 
      flst_real(   1,  58)=           2
      flst_real(   2,  58)=          -4
      flst_real(   3,  58)=          25
      flst_real(   4,  58)=         -11
      flst_real(   5,  58)=          11
      flst_real(   6,  58)=           2
      flst_real(   7,  58)=          -4
 
      flst_real(   1,  59)=           2
      flst_real(   2,  59)=           4
      flst_real(   3,  59)=          25
      flst_real(   4,  59)=         -11
      flst_real(   5,  59)=          11
      flst_real(   6,  59)=           2
      flst_real(   7,  59)=           4
 
      flst_real(   1,  60)=           2
      flst_real(   2,  60)=          -3
      flst_real(   3,  60)=          25
      flst_real(   4,  60)=         -11
      flst_real(   5,  60)=          11
      flst_real(   6,  60)=           2
      flst_real(   7,  60)=          -3
 
      flst_real(   1,  61)=           2
      flst_real(   2,  61)=           3
      flst_real(   3,  61)=          25
      flst_real(   4,  61)=         -11
      flst_real(   5,  61)=          11
      flst_real(   6,  61)=           2
      flst_real(   7,  61)=           3
 
      flst_real(   1,  62)=           2
      flst_real(   2,  62)=          -5
      flst_real(   3,  62)=          25
      flst_real(   4,  62)=         -11
      flst_real(   5,  62)=          11
      flst_real(   6,  62)=           2
      flst_real(   7,  62)=          -5
 
      flst_real(   1,  63)=           2
      flst_real(   2,  63)=           5
      flst_real(   3,  63)=          25
      flst_real(   4,  63)=         -11
      flst_real(   5,  63)=          11
      flst_real(   6,  63)=           2
      flst_real(   7,  63)=           5
 
      flst_real(   1,  64)=           2
      flst_real(   2,  64)=           0
      flst_real(   3,  64)=          25
      flst_real(   4,  64)=         -11
      flst_real(   5,  64)=          11
      flst_real(   6,  64)=           2
      flst_real(   7,  64)=           0
 
      flst_real(   1,  65)=          -4
      flst_real(   2,  65)=          -1
      flst_real(   3,  65)=          25
      flst_real(   4,  65)=         -11
      flst_real(   5,  65)=          11
      flst_real(   6,  65)=          -1
      flst_real(   7,  65)=          -4
 
      flst_real(   1,  66)=          -4
      flst_real(   2,  66)=           1
      flst_real(   3,  66)=          25
      flst_real(   4,  66)=         -11
      flst_real(   5,  66)=          11
      flst_real(   6,  66)=           1
      flst_real(   7,  66)=          -4
 
      flst_real(   1,  67)=          -4
      flst_real(   2,  67)=          -2
      flst_real(   3,  67)=          25
      flst_real(   4,  67)=         -11
      flst_real(   5,  67)=          11
      flst_real(   6,  67)=          -2
      flst_real(   7,  67)=          -4
 
      flst_real(   1,  68)=          -4
      flst_real(   2,  68)=           2
      flst_real(   3,  68)=          25
      flst_real(   4,  68)=         -11
      flst_real(   5,  68)=          11
      flst_real(   6,  68)=           2
      flst_real(   7,  68)=          -4
 
      flst_real(   1,  69)=          -4
      flst_real(   2,  69)=          -4
      flst_real(   3,  69)=          25
      flst_real(   4,  69)=         -11
      flst_real(   5,  69)=          11
      flst_real(   6,  69)=          -4
      flst_real(   7,  69)=          -4
 
      flst_real(   1,  70)=          -4
      flst_real(   2,  70)=           4
      flst_real(   3,  70)=          25
      flst_real(   4,  70)=         -11
      flst_real(   5,  70)=          11
      flst_real(   6,  70)=           1
      flst_real(   7,  70)=          -1
 
      flst_real(   1,  71)=          -4
      flst_real(   2,  71)=           4
      flst_real(   3,  71)=          25
      flst_real(   4,  71)=         -11
      flst_real(   5,  71)=          11
      flst_real(   6,  71)=           2
      flst_real(   7,  71)=          -2
 
      flst_real(   1,  72)=          -4
      flst_real(   2,  72)=           4
      flst_real(   3,  72)=          25
      flst_real(   4,  72)=         -11
      flst_real(   5,  72)=          11
      flst_real(   6,  72)=           4
      flst_real(   7,  72)=          -4
 
      flst_real(   1,  73)=          -4
      flst_real(   2,  73)=           4
      flst_real(   3,  73)=          25
      flst_real(   4,  73)=         -11
      flst_real(   5,  73)=          11
      flst_real(   6,  73)=           3
      flst_real(   7,  73)=          -3
 
      flst_real(   1,  74)=          -4
      flst_real(   2,  74)=           4
      flst_real(   3,  74)=          25
      flst_real(   4,  74)=         -11
      flst_real(   5,  74)=          11
      flst_real(   6,  74)=           5
      flst_real(   7,  74)=          -5
 
      flst_real(   1,  75)=          -4
      flst_real(   2,  75)=           4
      flst_real(   3,  75)=          25
      flst_real(   4,  75)=         -11
      flst_real(   5,  75)=          11
      flst_real(   6,  75)=           0
      flst_real(   7,  75)=           0
 
      flst_real(   1,  76)=          -4
      flst_real(   2,  76)=          -3
      flst_real(   3,  76)=          25
      flst_real(   4,  76)=         -11
      flst_real(   5,  76)=          11
      flst_real(   6,  76)=          -4
      flst_real(   7,  76)=          -3
 
      flst_real(   1,  77)=          -4
      flst_real(   2,  77)=           3
      flst_real(   3,  77)=          25
      flst_real(   4,  77)=         -11
      flst_real(   5,  77)=          11
      flst_real(   6,  77)=          -4
      flst_real(   7,  77)=           3
 
      flst_real(   1,  78)=          -4
      flst_real(   2,  78)=          -5
      flst_real(   3,  78)=          25
      flst_real(   4,  78)=         -11
      flst_real(   5,  78)=          11
      flst_real(   6,  78)=          -4
      flst_real(   7,  78)=          -5
 
      flst_real(   1,  79)=          -4
      flst_real(   2,  79)=           5
      flst_real(   3,  79)=          25
      flst_real(   4,  79)=         -11
      flst_real(   5,  79)=          11
      flst_real(   6,  79)=          -4
      flst_real(   7,  79)=           5
 
      flst_real(   1,  80)=          -4
      flst_real(   2,  80)=           0
      flst_real(   3,  80)=          25
      flst_real(   4,  80)=         -11
      flst_real(   5,  80)=          11
      flst_real(   6,  80)=          -4
      flst_real(   7,  80)=           0
 
      flst_real(   1,  81)=           4
      flst_real(   2,  81)=          -1
      flst_real(   3,  81)=          25
      flst_real(   4,  81)=         -11
      flst_real(   5,  81)=          11
      flst_real(   6,  81)=          -1
      flst_real(   7,  81)=           4
 
      flst_real(   1,  82)=           4
      flst_real(   2,  82)=           1
      flst_real(   3,  82)=          25
      flst_real(   4,  82)=         -11
      flst_real(   5,  82)=          11
      flst_real(   6,  82)=           1
      flst_real(   7,  82)=           4
 
      flst_real(   1,  83)=           4
      flst_real(   2,  83)=          -2
      flst_real(   3,  83)=          25
      flst_real(   4,  83)=         -11
      flst_real(   5,  83)=          11
      flst_real(   6,  83)=          -2
      flst_real(   7,  83)=           4
 
      flst_real(   1,  84)=           4
      flst_real(   2,  84)=           2
      flst_real(   3,  84)=          25
      flst_real(   4,  84)=         -11
      flst_real(   5,  84)=          11
      flst_real(   6,  84)=           2
      flst_real(   7,  84)=           4
 
      flst_real(   1,  85)=           4
      flst_real(   2,  85)=          -4
      flst_real(   3,  85)=          25
      flst_real(   4,  85)=         -11
      flst_real(   5,  85)=          11
      flst_real(   6,  85)=           1
      flst_real(   7,  85)=          -1
 
      flst_real(   1,  86)=           4
      flst_real(   2,  86)=          -4
      flst_real(   3,  86)=          25
      flst_real(   4,  86)=         -11
      flst_real(   5,  86)=          11
      flst_real(   6,  86)=           2
      flst_real(   7,  86)=          -2
 
      flst_real(   1,  87)=           4
      flst_real(   2,  87)=          -4
      flst_real(   3,  87)=          25
      flst_real(   4,  87)=         -11
      flst_real(   5,  87)=          11
      flst_real(   6,  87)=           4
      flst_real(   7,  87)=          -4
 
      flst_real(   1,  88)=           4
      flst_real(   2,  88)=          -4
      flst_real(   3,  88)=          25
      flst_real(   4,  88)=         -11
      flst_real(   5,  88)=          11
      flst_real(   6,  88)=           3
      flst_real(   7,  88)=          -3
 
      flst_real(   1,  89)=           4
      flst_real(   2,  89)=          -4
      flst_real(   3,  89)=          25
      flst_real(   4,  89)=         -11
      flst_real(   5,  89)=          11
      flst_real(   6,  89)=           5
      flst_real(   7,  89)=          -5
 
      flst_real(   1,  90)=           4
      flst_real(   2,  90)=          -4
      flst_real(   3,  90)=          25
      flst_real(   4,  90)=         -11
      flst_real(   5,  90)=          11
      flst_real(   6,  90)=           0
      flst_real(   7,  90)=           0
 
      flst_real(   1,  91)=           4
      flst_real(   2,  91)=           4
      flst_real(   3,  91)=          25
      flst_real(   4,  91)=         -11
      flst_real(   5,  91)=          11
      flst_real(   6,  91)=           4
      flst_real(   7,  91)=           4
 
      flst_real(   1,  92)=           4
      flst_real(   2,  92)=          -3
      flst_real(   3,  92)=          25
      flst_real(   4,  92)=         -11
      flst_real(   5,  92)=          11
      flst_real(   6,  92)=           4
      flst_real(   7,  92)=          -3
 
      flst_real(   1,  93)=           4
      flst_real(   2,  93)=           3
      flst_real(   3,  93)=          25
      flst_real(   4,  93)=         -11
      flst_real(   5,  93)=          11
      flst_real(   6,  93)=           4
      flst_real(   7,  93)=           3
 
      flst_real(   1,  94)=           4
      flst_real(   2,  94)=          -5
      flst_real(   3,  94)=          25
      flst_real(   4,  94)=         -11
      flst_real(   5,  94)=          11
      flst_real(   6,  94)=           4
      flst_real(   7,  94)=          -5
 
      flst_real(   1,  95)=           4
      flst_real(   2,  95)=           5
      flst_real(   3,  95)=          25
      flst_real(   4,  95)=         -11
      flst_real(   5,  95)=          11
      flst_real(   6,  95)=           4
      flst_real(   7,  95)=           5
 
      flst_real(   1,  96)=           4
      flst_real(   2,  96)=           0
      flst_real(   3,  96)=          25
      flst_real(   4,  96)=         -11
      flst_real(   5,  96)=          11
      flst_real(   6,  96)=           4
      flst_real(   7,  96)=           0
 
      flst_real(   1,  97)=          -3
      flst_real(   2,  97)=          -1
      flst_real(   3,  97)=          25
      flst_real(   4,  97)=         -11
      flst_real(   5,  97)=          11
      flst_real(   6,  97)=          -1
      flst_real(   7,  97)=          -3
 
      flst_real(   1,  98)=          -3
      flst_real(   2,  98)=           1
      flst_real(   3,  98)=          25
      flst_real(   4,  98)=         -11
      flst_real(   5,  98)=          11
      flst_real(   6,  98)=           1
      flst_real(   7,  98)=          -3
 
      flst_real(   1,  99)=          -3
      flst_real(   2,  99)=          -2
      flst_real(   3,  99)=          25
      flst_real(   4,  99)=         -11
      flst_real(   5,  99)=          11
      flst_real(   6,  99)=          -2
      flst_real(   7,  99)=          -3
 
      flst_real(   1, 100)=          -3
      flst_real(   2, 100)=           2
      flst_real(   3, 100)=          25
      flst_real(   4, 100)=         -11
      flst_real(   5, 100)=          11
      flst_real(   6, 100)=           2
      flst_real(   7, 100)=          -3
 
      flst_real(   1, 101)=          -3
      flst_real(   2, 101)=          -4
      flst_real(   3, 101)=          25
      flst_real(   4, 101)=         -11
      flst_real(   5, 101)=          11
      flst_real(   6, 101)=          -4
      flst_real(   7, 101)=          -3
 
      flst_real(   1, 102)=          -3
      flst_real(   2, 102)=           4
      flst_real(   3, 102)=          25
      flst_real(   4, 102)=         -11
      flst_real(   5, 102)=          11
      flst_real(   6, 102)=           4
      flst_real(   7, 102)=          -3
 
      flst_real(   1, 103)=          -3
      flst_real(   2, 103)=          -3
      flst_real(   3, 103)=          25
      flst_real(   4, 103)=         -11
      flst_real(   5, 103)=          11
      flst_real(   6, 103)=          -3
      flst_real(   7, 103)=          -3
 
      flst_real(   1, 104)=          -3
      flst_real(   2, 104)=           3
      flst_real(   3, 104)=          25
      flst_real(   4, 104)=         -11
      flst_real(   5, 104)=          11
      flst_real(   6, 104)=           1
      flst_real(   7, 104)=          -1
 
      flst_real(   1, 105)=          -3
      flst_real(   2, 105)=           3
      flst_real(   3, 105)=          25
      flst_real(   4, 105)=         -11
      flst_real(   5, 105)=          11
      flst_real(   6, 105)=           2
      flst_real(   7, 105)=          -2
 
      flst_real(   1, 106)=          -3
      flst_real(   2, 106)=           3
      flst_real(   3, 106)=          25
      flst_real(   4, 106)=         -11
      flst_real(   5, 106)=          11
      flst_real(   6, 106)=           4
      flst_real(   7, 106)=          -4
 
      flst_real(   1, 107)=          -3
      flst_real(   2, 107)=           3
      flst_real(   3, 107)=          25
      flst_real(   4, 107)=         -11
      flst_real(   5, 107)=          11
      flst_real(   6, 107)=           3
      flst_real(   7, 107)=          -3
 
      flst_real(   1, 108)=          -3
      flst_real(   2, 108)=           3
      flst_real(   3, 108)=          25
      flst_real(   4, 108)=         -11
      flst_real(   5, 108)=          11
      flst_real(   6, 108)=           5
      flst_real(   7, 108)=          -5
 
      flst_real(   1, 109)=          -3
      flst_real(   2, 109)=           3
      flst_real(   3, 109)=          25
      flst_real(   4, 109)=         -11
      flst_real(   5, 109)=          11
      flst_real(   6, 109)=           0
      flst_real(   7, 109)=           0
 
      flst_real(   1, 110)=          -3
      flst_real(   2, 110)=          -5
      flst_real(   3, 110)=          25
      flst_real(   4, 110)=         -11
      flst_real(   5, 110)=          11
      flst_real(   6, 110)=          -3
      flst_real(   7, 110)=          -5
 
      flst_real(   1, 111)=          -3
      flst_real(   2, 111)=           5
      flst_real(   3, 111)=          25
      flst_real(   4, 111)=         -11
      flst_real(   5, 111)=          11
      flst_real(   6, 111)=          -3
      flst_real(   7, 111)=           5
 
      flst_real(   1, 112)=          -3
      flst_real(   2, 112)=           0
      flst_real(   3, 112)=          25
      flst_real(   4, 112)=         -11
      flst_real(   5, 112)=          11
      flst_real(   6, 112)=          -3
      flst_real(   7, 112)=           0
 
      flst_real(   1, 113)=           3
      flst_real(   2, 113)=          -1
      flst_real(   3, 113)=          25
      flst_real(   4, 113)=         -11
      flst_real(   5, 113)=          11
      flst_real(   6, 113)=          -1
      flst_real(   7, 113)=           3
 
      flst_real(   1, 114)=           3
      flst_real(   2, 114)=           1
      flst_real(   3, 114)=          25
      flst_real(   4, 114)=         -11
      flst_real(   5, 114)=          11
      flst_real(   6, 114)=           1
      flst_real(   7, 114)=           3
 
      flst_real(   1, 115)=           3
      flst_real(   2, 115)=          -2
      flst_real(   3, 115)=          25
      flst_real(   4, 115)=         -11
      flst_real(   5, 115)=          11
      flst_real(   6, 115)=          -2
      flst_real(   7, 115)=           3
 
      flst_real(   1, 116)=           3
      flst_real(   2, 116)=           2
      flst_real(   3, 116)=          25
      flst_real(   4, 116)=         -11
      flst_real(   5, 116)=          11
      flst_real(   6, 116)=           2
      flst_real(   7, 116)=           3
 
      flst_real(   1, 117)=           3
      flst_real(   2, 117)=          -4
      flst_real(   3, 117)=          25
      flst_real(   4, 117)=         -11
      flst_real(   5, 117)=          11
      flst_real(   6, 117)=          -4
      flst_real(   7, 117)=           3
 
      flst_real(   1, 118)=           3
      flst_real(   2, 118)=           4
      flst_real(   3, 118)=          25
      flst_real(   4, 118)=         -11
      flst_real(   5, 118)=          11
      flst_real(   6, 118)=           4
      flst_real(   7, 118)=           3
 
      flst_real(   1, 119)=           3
      flst_real(   2, 119)=          -3
      flst_real(   3, 119)=          25
      flst_real(   4, 119)=         -11
      flst_real(   5, 119)=          11
      flst_real(   6, 119)=           1
      flst_real(   7, 119)=          -1
 
      flst_real(   1, 120)=           3
      flst_real(   2, 120)=          -3
      flst_real(   3, 120)=          25
      flst_real(   4, 120)=         -11
      flst_real(   5, 120)=          11
      flst_real(   6, 120)=           2
      flst_real(   7, 120)=          -2
 
      flst_real(   1, 121)=           3
      flst_real(   2, 121)=          -3
      flst_real(   3, 121)=          25
      flst_real(   4, 121)=         -11
      flst_real(   5, 121)=          11
      flst_real(   6, 121)=           4
      flst_real(   7, 121)=          -4
 
      flst_real(   1, 122)=           3
      flst_real(   2, 122)=          -3
      flst_real(   3, 122)=          25
      flst_real(   4, 122)=         -11
      flst_real(   5, 122)=          11
      flst_real(   6, 122)=           3
      flst_real(   7, 122)=          -3
 
      flst_real(   1, 123)=           3
      flst_real(   2, 123)=          -3
      flst_real(   3, 123)=          25
      flst_real(   4, 123)=         -11
      flst_real(   5, 123)=          11
      flst_real(   6, 123)=           5
      flst_real(   7, 123)=          -5
 
      flst_real(   1, 124)=           3
      flst_real(   2, 124)=          -3
      flst_real(   3, 124)=          25
      flst_real(   4, 124)=         -11
      flst_real(   5, 124)=          11
      flst_real(   6, 124)=           0
      flst_real(   7, 124)=           0
 
      flst_real(   1, 125)=           3
      flst_real(   2, 125)=           3
      flst_real(   3, 125)=          25
      flst_real(   4, 125)=         -11
      flst_real(   5, 125)=          11
      flst_real(   6, 125)=           3
      flst_real(   7, 125)=           3
 
      flst_real(   1, 126)=           3
      flst_real(   2, 126)=          -5
      flst_real(   3, 126)=          25
      flst_real(   4, 126)=         -11
      flst_real(   5, 126)=          11
      flst_real(   6, 126)=           3
      flst_real(   7, 126)=          -5
 
      flst_real(   1, 127)=           3
      flst_real(   2, 127)=           5
      flst_real(   3, 127)=          25
      flst_real(   4, 127)=         -11
      flst_real(   5, 127)=          11
      flst_real(   6, 127)=           3
      flst_real(   7, 127)=           5
 
      flst_real(   1, 128)=           3
      flst_real(   2, 128)=           0
      flst_real(   3, 128)=          25
      flst_real(   4, 128)=         -11
      flst_real(   5, 128)=          11
      flst_real(   6, 128)=           3
      flst_real(   7, 128)=           0
 
      flst_real(   1, 129)=          -5
      flst_real(   2, 129)=          -1
      flst_real(   3, 129)=          25
      flst_real(   4, 129)=         -11
      flst_real(   5, 129)=          11
      flst_real(   6, 129)=          -1
      flst_real(   7, 129)=          -5
 
      flst_real(   1, 130)=          -5
      flst_real(   2, 130)=           1
      flst_real(   3, 130)=          25
      flst_real(   4, 130)=         -11
      flst_real(   5, 130)=          11
      flst_real(   6, 130)=           1
      flst_real(   7, 130)=          -5
 
      flst_real(   1, 131)=          -5
      flst_real(   2, 131)=          -2
      flst_real(   3, 131)=          25
      flst_real(   4, 131)=         -11
      flst_real(   5, 131)=          11
      flst_real(   6, 131)=          -2
      flst_real(   7, 131)=          -5
 
      flst_real(   1, 132)=          -5
      flst_real(   2, 132)=           2
      flst_real(   3, 132)=          25
      flst_real(   4, 132)=         -11
      flst_real(   5, 132)=          11
      flst_real(   6, 132)=           2
      flst_real(   7, 132)=          -5
 
      flst_real(   1, 133)=          -5
      flst_real(   2, 133)=          -4
      flst_real(   3, 133)=          25
      flst_real(   4, 133)=         -11
      flst_real(   5, 133)=          11
      flst_real(   6, 133)=          -4
      flst_real(   7, 133)=          -5
 
      flst_real(   1, 134)=          -5
      flst_real(   2, 134)=           4
      flst_real(   3, 134)=          25
      flst_real(   4, 134)=         -11
      flst_real(   5, 134)=          11
      flst_real(   6, 134)=           4
      flst_real(   7, 134)=          -5
 
      flst_real(   1, 135)=          -5
      flst_real(   2, 135)=          -3
      flst_real(   3, 135)=          25
      flst_real(   4, 135)=         -11
      flst_real(   5, 135)=          11
      flst_real(   6, 135)=          -3
      flst_real(   7, 135)=          -5
 
      flst_real(   1, 136)=          -5
      flst_real(   2, 136)=           3
      flst_real(   3, 136)=          25
      flst_real(   4, 136)=         -11
      flst_real(   5, 136)=          11
      flst_real(   6, 136)=           3
      flst_real(   7, 136)=          -5
 
      flst_real(   1, 137)=          -5
      flst_real(   2, 137)=          -5
      flst_real(   3, 137)=          25
      flst_real(   4, 137)=         -11
      flst_real(   5, 137)=          11
      flst_real(   6, 137)=          -5
      flst_real(   7, 137)=          -5
 
      flst_real(   1, 138)=          -5
      flst_real(   2, 138)=           5
      flst_real(   3, 138)=          25
      flst_real(   4, 138)=         -11
      flst_real(   5, 138)=          11
      flst_real(   6, 138)=           1
      flst_real(   7, 138)=          -1
 
      flst_real(   1, 139)=          -5
      flst_real(   2, 139)=           5
      flst_real(   3, 139)=          25
      flst_real(   4, 139)=         -11
      flst_real(   5, 139)=          11
      flst_real(   6, 139)=           2
      flst_real(   7, 139)=          -2
 
      flst_real(   1, 140)=          -5
      flst_real(   2, 140)=           5
      flst_real(   3, 140)=          25
      flst_real(   4, 140)=         -11
      flst_real(   5, 140)=          11
      flst_real(   6, 140)=           4
      flst_real(   7, 140)=          -4
 
      flst_real(   1, 141)=          -5
      flst_real(   2, 141)=           5
      flst_real(   3, 141)=          25
      flst_real(   4, 141)=         -11
      flst_real(   5, 141)=          11
      flst_real(   6, 141)=           3
      flst_real(   7, 141)=          -3
 
      flst_real(   1, 142)=          -5
      flst_real(   2, 142)=           5
      flst_real(   3, 142)=          25
      flst_real(   4, 142)=         -11
      flst_real(   5, 142)=          11
      flst_real(   6, 142)=           5
      flst_real(   7, 142)=          -5
 
      flst_real(   1, 143)=          -5
      flst_real(   2, 143)=           5
      flst_real(   3, 143)=          25
      flst_real(   4, 143)=         -11
      flst_real(   5, 143)=          11
      flst_real(   6, 143)=           0
      flst_real(   7, 143)=           0
 
      flst_real(   1, 144)=          -5
      flst_real(   2, 144)=           0
      flst_real(   3, 144)=          25
      flst_real(   4, 144)=         -11
      flst_real(   5, 144)=          11
      flst_real(   6, 144)=          -5
      flst_real(   7, 144)=           0
 
      flst_real(   1, 145)=           5
      flst_real(   2, 145)=          -1
      flst_real(   3, 145)=          25
      flst_real(   4, 145)=         -11
      flst_real(   5, 145)=          11
      flst_real(   6, 145)=          -1
      flst_real(   7, 145)=           5
 
      flst_real(   1, 146)=           5
      flst_real(   2, 146)=           1
      flst_real(   3, 146)=          25
      flst_real(   4, 146)=         -11
      flst_real(   5, 146)=          11
      flst_real(   6, 146)=           1
      flst_real(   7, 146)=           5
 
      flst_real(   1, 147)=           5
      flst_real(   2, 147)=          -2
      flst_real(   3, 147)=          25
      flst_real(   4, 147)=         -11
      flst_real(   5, 147)=          11
      flst_real(   6, 147)=          -2
      flst_real(   7, 147)=           5
 
      flst_real(   1, 148)=           5
      flst_real(   2, 148)=           2
      flst_real(   3, 148)=          25
      flst_real(   4, 148)=         -11
      flst_real(   5, 148)=          11
      flst_real(   6, 148)=           2
      flst_real(   7, 148)=           5
 
      flst_real(   1, 149)=           5
      flst_real(   2, 149)=          -4
      flst_real(   3, 149)=          25
      flst_real(   4, 149)=         -11
      flst_real(   5, 149)=          11
      flst_real(   6, 149)=          -4
      flst_real(   7, 149)=           5
 
      flst_real(   1, 150)=           5
      flst_real(   2, 150)=           4
      flst_real(   3, 150)=          25
      flst_real(   4, 150)=         -11
      flst_real(   5, 150)=          11
      flst_real(   6, 150)=           4
      flst_real(   7, 150)=           5
 
      flst_real(   1, 151)=           5
      flst_real(   2, 151)=          -3
      flst_real(   3, 151)=          25
      flst_real(   4, 151)=         -11
      flst_real(   5, 151)=          11
      flst_real(   6, 151)=          -3
      flst_real(   7, 151)=           5
 
      flst_real(   1, 152)=           5
      flst_real(   2, 152)=           3
      flst_real(   3, 152)=          25
      flst_real(   4, 152)=         -11
      flst_real(   5, 152)=          11
      flst_real(   6, 152)=           3
      flst_real(   7, 152)=           5
 
      flst_real(   1, 153)=           5
      flst_real(   2, 153)=          -5
      flst_real(   3, 153)=          25
      flst_real(   4, 153)=         -11
      flst_real(   5, 153)=          11
      flst_real(   6, 153)=           1
      flst_real(   7, 153)=          -1
 
      flst_real(   1, 154)=           5
      flst_real(   2, 154)=          -5
      flst_real(   3, 154)=          25
      flst_real(   4, 154)=         -11
      flst_real(   5, 154)=          11
      flst_real(   6, 154)=           2
      flst_real(   7, 154)=          -2
 
      flst_real(   1, 155)=           5
      flst_real(   2, 155)=          -5
      flst_real(   3, 155)=          25
      flst_real(   4, 155)=         -11
      flst_real(   5, 155)=          11
      flst_real(   6, 155)=           4
      flst_real(   7, 155)=          -4
 
      flst_real(   1, 156)=           5
      flst_real(   2, 156)=          -5
      flst_real(   3, 156)=          25
      flst_real(   4, 156)=         -11
      flst_real(   5, 156)=          11
      flst_real(   6, 156)=           3
      flst_real(   7, 156)=          -3
 
      flst_real(   1, 157)=           5
      flst_real(   2, 157)=          -5
      flst_real(   3, 157)=          25
      flst_real(   4, 157)=         -11
      flst_real(   5, 157)=          11
      flst_real(   6, 157)=           5
      flst_real(   7, 157)=          -5
 
      flst_real(   1, 158)=           5
      flst_real(   2, 158)=          -5
      flst_real(   3, 158)=          25
      flst_real(   4, 158)=         -11
      flst_real(   5, 158)=          11
      flst_real(   6, 158)=           0
      flst_real(   7, 158)=           0
 
      flst_real(   1, 159)=           5
      flst_real(   2, 159)=           5
      flst_real(   3, 159)=          25
      flst_real(   4, 159)=         -11
      flst_real(   5, 159)=          11
      flst_real(   6, 159)=           5
      flst_real(   7, 159)=           5
 
      flst_real(   1, 160)=           5
      flst_real(   2, 160)=           0
      flst_real(   3, 160)=          25
      flst_real(   4, 160)=         -11
      flst_real(   5, 160)=          11
      flst_real(   6, 160)=           5
      flst_real(   7, 160)=           0
 
      flst_real(   1, 161)=           0
      flst_real(   2, 161)=          -1
      flst_real(   3, 161)=          25
      flst_real(   4, 161)=         -11
      flst_real(   5, 161)=          11
      flst_real(   6, 161)=          -1
      flst_real(   7, 161)=           0
 
      flst_real(   1, 162)=           0
      flst_real(   2, 162)=           1
      flst_real(   3, 162)=          25
      flst_real(   4, 162)=         -11
      flst_real(   5, 162)=          11
      flst_real(   6, 162)=           1
      flst_real(   7, 162)=           0
 
      flst_real(   1, 163)=           0
      flst_real(   2, 163)=          -2
      flst_real(   3, 163)=          25
      flst_real(   4, 163)=         -11
      flst_real(   5, 163)=          11
      flst_real(   6, 163)=          -2
      flst_real(   7, 163)=           0
 
      flst_real(   1, 164)=           0
      flst_real(   2, 164)=           2
      flst_real(   3, 164)=          25
      flst_real(   4, 164)=         -11
      flst_real(   5, 164)=          11
      flst_real(   6, 164)=           2
      flst_real(   7, 164)=           0
 
      flst_real(   1, 165)=           0
      flst_real(   2, 165)=          -4
      flst_real(   3, 165)=          25
      flst_real(   4, 165)=         -11
      flst_real(   5, 165)=          11
      flst_real(   6, 165)=          -4
      flst_real(   7, 165)=           0
 
      flst_real(   1, 166)=           0
      flst_real(   2, 166)=           4
      flst_real(   3, 166)=          25
      flst_real(   4, 166)=         -11
      flst_real(   5, 166)=          11
      flst_real(   6, 166)=           4
      flst_real(   7, 166)=           0
 
      flst_real(   1, 167)=           0
      flst_real(   2, 167)=          -3
      flst_real(   3, 167)=          25
      flst_real(   4, 167)=         -11
      flst_real(   5, 167)=          11
      flst_real(   6, 167)=          -3
      flst_real(   7, 167)=           0
 
      flst_real(   1, 168)=           0
      flst_real(   2, 168)=           3
      flst_real(   3, 168)=          25
      flst_real(   4, 168)=         -11
      flst_real(   5, 168)=          11
      flst_real(   6, 168)=           3
      flst_real(   7, 168)=           0
 
      flst_real(   1, 169)=           0
      flst_real(   2, 169)=          -5
      flst_real(   3, 169)=          25
      flst_real(   4, 169)=         -11
      flst_real(   5, 169)=          11
      flst_real(   6, 169)=          -5
      flst_real(   7, 169)=           0
 
      flst_real(   1, 170)=           0
      flst_real(   2, 170)=           5
      flst_real(   3, 170)=          25
      flst_real(   4, 170)=         -11
      flst_real(   5, 170)=          11
      flst_real(   6, 170)=           5
      flst_real(   7, 170)=           0
 
      flst_real(   1, 171)=           0
      flst_real(   2, 171)=           0
      flst_real(   3, 171)=          25
      flst_real(   4, 171)=         -11
      flst_real(   5, 171)=          11
      flst_real(   6, 171)=           1
      flst_real(   7, 171)=          -1
 
      flst_real(   1, 172)=           0
      flst_real(   2, 172)=           0
      flst_real(   3, 172)=          25
      flst_real(   4, 172)=         -11
      flst_real(   5, 172)=          11
      flst_real(   6, 172)=           2
      flst_real(   7, 172)=          -2
 
      flst_real(   1, 173)=           0
      flst_real(   2, 173)=           0
      flst_real(   3, 173)=          25
      flst_real(   4, 173)=         -11
      flst_real(   5, 173)=          11
      flst_real(   6, 173)=           4
      flst_real(   7, 173)=          -4
 
      flst_real(   1, 174)=           0
      flst_real(   2, 174)=           0
      flst_real(   3, 174)=          25
      flst_real(   4, 174)=         -11
      flst_real(   5, 174)=          11
      flst_real(   6, 174)=           3
      flst_real(   7, 174)=          -3
 
      flst_real(   1, 175)=           0
      flst_real(   2, 175)=           0
      flst_real(   3, 175)=          25
      flst_real(   4, 175)=         -11
      flst_real(   5, 175)=          11
      flst_real(   6, 175)=           5
      flst_real(   7, 175)=          -5
 
      flst_nreal=         175
 
      return
      end
 
