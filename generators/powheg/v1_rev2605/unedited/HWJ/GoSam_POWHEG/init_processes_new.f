      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include "coupl.inc"
      integer i

      call init_processes_born
      call init_processes_real
      call init_couplings
      if (tmass.eq.0d0) then
         st_nlight=6
      elseif(bmass.eq.0d0) then
         st_nlight=5
      elseif(cmass.eq.0d0) then
         st_nlight=4
      else
         st_nlight=3
      endif
      do i=3,nlegreal
         if (abs(flst_real(i,1)).le.st_nlight) then
            flst_lightpart=i
            exit
         endif
      enddo

      write(*,*) 'st_nlight ',st_nlight
      stop
      return
      end



      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"

      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           2
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=         -11
      flst_born(   5,   1)=          12
      flst_born(   6,   1)=           0

      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           4
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=         -11
      flst_born(   5,   2)=          12
      flst_born(   6,   2)=           0

      flst_born(   1,   3)=          -1
      flst_born(   2,   3)=           0
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=         -11
      flst_born(   5,   3)=          12
      flst_born(   6,   3)=          -2

      flst_born(   1,   4)=          -1
      flst_born(   2,   4)=           0
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=         -11
      flst_born(   5,   4)=          12
      flst_born(   6,   4)=          -4

      flst_born(   1,   5)=           2
      flst_born(   2,   5)=          -1
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=         -11
      flst_born(   5,   5)=          12
      flst_born(   6,   5)=           0

      flst_born(   1,   6)=           2
      flst_born(   2,   6)=          -3
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=         -11
      flst_born(   5,   6)=          12
      flst_born(   6,   6)=           0

      flst_born(   1,   7)=           2
      flst_born(   2,   7)=          -5
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=         -11
      flst_born(   5,   7)=          12
      flst_born(   6,   7)=           0

      flst_born(   1,   8)=           2
      flst_born(   2,   8)=           0
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=         -11
      flst_born(   5,   8)=          12
      flst_born(   6,   8)=           1

      flst_born(   1,   9)=           2
      flst_born(   2,   9)=           0
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=         -11
      flst_born(   5,   9)=          12
      flst_born(   6,   9)=           3

      flst_born(   1,  10)=           2
      flst_born(   2,  10)=           0
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=         -11
      flst_born(   5,  10)=          12
      flst_born(   6,  10)=           5

      flst_born(   1,  11)=           4
      flst_born(   2,  11)=          -1
      flst_born(   3,  11)=          25
      flst_born(   4,  11)=         -11
      flst_born(   5,  11)=          12
      flst_born(   6,  11)=           0

      flst_born(   1,  12)=           4
      flst_born(   2,  12)=          -3
      flst_born(   3,  12)=          25
      flst_born(   4,  12)=         -11
      flst_born(   5,  12)=          12
      flst_born(   6,  12)=           0

      flst_born(   1,  13)=           4
      flst_born(   2,  13)=          -5
      flst_born(   3,  13)=          25
      flst_born(   4,  13)=         -11
      flst_born(   5,  13)=          12
      flst_born(   6,  13)=           0

      flst_born(   1,  14)=           4
      flst_born(   2,  14)=           0
      flst_born(   3,  14)=          25
      flst_born(   4,  14)=         -11
      flst_born(   5,  14)=          12
      flst_born(   6,  14)=           1

      flst_born(   1,  15)=           4
      flst_born(   2,  15)=           0
      flst_born(   3,  15)=          25
      flst_born(   4,  15)=         -11
      flst_born(   5,  15)=          12
      flst_born(   6,  15)=           3

      flst_born(   1,  16)=           4
      flst_born(   2,  16)=           0
      flst_born(   3,  16)=          25
      flst_born(   4,  16)=         -11
      flst_born(   5,  16)=          12
      flst_born(   6,  16)=           5

      flst_born(   1,  17)=          -3
      flst_born(   2,  17)=           2
      flst_born(   3,  17)=          25
      flst_born(   4,  17)=         -11
      flst_born(   5,  17)=          12
      flst_born(   6,  17)=           0

      flst_born(   1,  18)=          -3
      flst_born(   2,  18)=           4
      flst_born(   3,  18)=          25
      flst_born(   4,  18)=         -11
      flst_born(   5,  18)=          12
      flst_born(   6,  18)=           0

      flst_born(   1,  19)=          -3
      flst_born(   2,  19)=           0
      flst_born(   3,  19)=          25
      flst_born(   4,  19)=         -11
      flst_born(   5,  19)=          12
      flst_born(   6,  19)=          -2

      flst_born(   1,  20)=          -3
      flst_born(   2,  20)=           0
      flst_born(   3,  20)=          25
      flst_born(   4,  20)=         -11
      flst_born(   5,  20)=          12
      flst_born(   6,  20)=          -4

      flst_born(   1,  21)=          -5
      flst_born(   2,  21)=           2
      flst_born(   3,  21)=          25
      flst_born(   4,  21)=         -11
      flst_born(   5,  21)=          12
      flst_born(   6,  21)=           0

      flst_born(   1,  22)=          -5
      flst_born(   2,  22)=           4
      flst_born(   3,  22)=          25
      flst_born(   4,  22)=         -11
      flst_born(   5,  22)=          12
      flst_born(   6,  22)=           0

      flst_born(   1,  23)=          -5
      flst_born(   2,  23)=           0
      flst_born(   3,  23)=          25
      flst_born(   4,  23)=         -11
      flst_born(   5,  23)=          12
      flst_born(   6,  23)=          -2

      flst_born(   1,  24)=          -5
      flst_born(   2,  24)=           0
      flst_born(   3,  24)=          25
      flst_born(   4,  24)=         -11
      flst_born(   5,  24)=          12
      flst_born(   6,  24)=          -4

      flst_born(   1,  25)=           0
      flst_born(   2,  25)=          -1
      flst_born(   3,  25)=          25
      flst_born(   4,  25)=         -11
      flst_born(   5,  25)=          12
      flst_born(   6,  25)=          -2

      flst_born(   1,  26)=           0
      flst_born(   2,  26)=          -1
      flst_born(   3,  26)=          25
      flst_born(   4,  26)=         -11
      flst_born(   5,  26)=          12
      flst_born(   6,  26)=          -4

      flst_born(   1,  27)=           0
      flst_born(   2,  27)=           2
      flst_born(   3,  27)=          25
      flst_born(   4,  27)=         -11
      flst_born(   5,  27)=          12
      flst_born(   6,  27)=           1

      flst_born(   1,  28)=           0
      flst_born(   2,  28)=           2
      flst_born(   3,  28)=          25
      flst_born(   4,  28)=         -11
      flst_born(   5,  28)=          12
      flst_born(   6,  28)=           3

      flst_born(   1,  29)=           0
      flst_born(   2,  29)=           2
      flst_born(   3,  29)=          25
      flst_born(   4,  29)=         -11
      flst_born(   5,  29)=          12
      flst_born(   6,  29)=           5

      flst_born(   1,  30)=           0
      flst_born(   2,  30)=           4
      flst_born(   3,  30)=          25
      flst_born(   4,  30)=         -11
      flst_born(   5,  30)=          12
      flst_born(   6,  30)=           1

      flst_born(   1,  31)=           0
      flst_born(   2,  31)=           4
      flst_born(   3,  31)=          25
      flst_born(   4,  31)=         -11
      flst_born(   5,  31)=          12
      flst_born(   6,  31)=           3

      flst_born(   1,  32)=           0
      flst_born(   2,  32)=           4
      flst_born(   3,  32)=          25
      flst_born(   4,  32)=         -11
      flst_born(   5,  32)=          12
      flst_born(   6,  32)=           5

      flst_born(   1,  33)=           0
      flst_born(   2,  33)=          -3
      flst_born(   3,  33)=          25
      flst_born(   4,  33)=         -11
      flst_born(   5,  33)=          12
      flst_born(   6,  33)=          -2

      flst_born(   1,  34)=           0
      flst_born(   2,  34)=          -3
      flst_born(   3,  34)=          25
      flst_born(   4,  34)=         -11
      flst_born(   5,  34)=          12
      flst_born(   6,  34)=          -4

      flst_born(   1,  35)=           0
      flst_born(   2,  35)=          -5
      flst_born(   3,  35)=          25
      flst_born(   4,  35)=         -11
      flst_born(   5,  35)=          12
      flst_born(   6,  35)=          -2

      flst_born(   1,  36)=           0
      flst_born(   2,  36)=          -5
      flst_born(   3,  36)=          25
      flst_born(   4,  36)=         -11
      flst_born(   5,  36)=          12
      flst_born(   6,  36)=          -4

      flst_nborn=          36


c     call gosam_flst_born(1,3)

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
      flst_real(   5,   1)=          12
      flst_real(   6,   1)=          -1
      flst_real(   7,   1)=          -2

      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=          -1
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=         -11
      flst_real(   5,   2)=          12
      flst_real(   6,   2)=          -1
      flst_real(   7,   2)=          -4

      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=           1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=         -11
      flst_real(   5,   3)=          12
      flst_real(   6,   3)=           1
      flst_real(   7,   3)=          -2

      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           1
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=         -11
      flst_real(   5,   4)=          12
      flst_real(   6,   4)=           1
      flst_real(   7,   4)=          -4

      flst_real(   1,   5)=          -1
      flst_real(   2,   5)=           1
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=         -11
      flst_real(   5,   5)=          12
      flst_real(   6,   5)=          -2
      flst_real(   7,   5)=           3

      flst_real(   1,   6)=          -1
      flst_real(   2,   6)=           1
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=         -11
      flst_real(   5,   6)=          12
      flst_real(   6,   6)=          -2
      flst_real(   7,   6)=           5

      flst_real(   1,   7)=          -1
      flst_real(   2,   7)=           1
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=         -11
      flst_real(   5,   7)=          12
      flst_real(   6,   7)=          -4
      flst_real(   7,   7)=           3

      flst_real(   1,   8)=          -1
      flst_real(   2,   8)=           1
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=         -11
      flst_real(   5,   8)=          12
      flst_real(   6,   8)=          -4
      flst_real(   7,   8)=           5

      flst_real(   1,   9)=          -1
      flst_real(   2,   9)=          -2
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=         -11
      flst_real(   5,   9)=          12
      flst_real(   6,   9)=          -2
      flst_real(   7,   9)=          -2

      flst_real(   1,  10)=          -1
      flst_real(   2,  10)=          -2
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=         -11
      flst_real(   5,  10)=          12
      flst_real(   6,  10)=          -2
      flst_real(   7,  10)=          -4

      flst_real(   1,  11)=          -1
      flst_real(   2,  11)=           2
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=         -11
      flst_real(   5,  11)=          12
      flst_real(   6,  11)=           1
      flst_real(   7,  11)=          -1

      flst_real(   1,  12)=          -1
      flst_real(   2,  12)=           2
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=         -11
      flst_real(   5,  12)=          12
      flst_real(   6,  12)=          -1
      flst_real(   7,  12)=           3

      flst_real(   1,  13)=          -1
      flst_real(   2,  13)=           2
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=         -11
      flst_real(   5,  13)=          12
      flst_real(   6,  13)=          -1
      flst_real(   7,  13)=           5

      flst_real(   1,  14)=          -1
      flst_real(   2,  14)=           2
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=         -11
      flst_real(   5,  14)=          12
      flst_real(   6,  14)=           2
      flst_real(   7,  14)=          -2

      flst_real(   1,  15)=          -1
      flst_real(   2,  15)=           2
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=         -11
      flst_real(   5,  15)=          12
      flst_real(   6,  15)=           2
      flst_real(   7,  15)=          -4

      flst_real(   1,  16)=          -1
      flst_real(   2,  16)=           2
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=         -11
      flst_real(   5,  16)=          12
      flst_real(   6,  16)=           4
      flst_real(   7,  16)=          -4

      flst_real(   1,  17)=          -1
      flst_real(   2,  17)=           2
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=         -11
      flst_real(   5,  17)=          12
      flst_real(   6,  17)=           3
      flst_real(   7,  17)=          -3

      flst_real(   1,  18)=          -1
      flst_real(   2,  18)=           2
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=         -11
      flst_real(   5,  18)=          12
      flst_real(   6,  18)=           5
      flst_real(   7,  18)=          -5

      flst_real(   1,  19)=          -1
      flst_real(   2,  19)=           2
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=         -11
      flst_real(   5,  19)=          12
      flst_real(   6,  19)=           0
      flst_real(   7,  19)=           0

      flst_real(   1,  20)=          -1
      flst_real(   2,  20)=          -4
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=         -11
      flst_real(   5,  20)=          12
      flst_real(   6,  20)=          -2
      flst_real(   7,  20)=          -4

      flst_real(   1,  21)=          -1
      flst_real(   2,  21)=          -4
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=         -11
      flst_real(   5,  21)=          12
      flst_real(   6,  21)=          -4
      flst_real(   7,  21)=          -4

      flst_real(   1,  22)=          -1
      flst_real(   2,  22)=           4
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=         -11
      flst_real(   5,  22)=          12
      flst_real(   6,  22)=           1
      flst_real(   7,  22)=          -1

      flst_real(   1,  23)=          -1
      flst_real(   2,  23)=           4
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=         -11
      flst_real(   5,  23)=          12
      flst_real(   6,  23)=          -1
      flst_real(   7,  23)=           3

      flst_real(   1,  24)=          -1
      flst_real(   2,  24)=           4
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=         -11
      flst_real(   5,  24)=          12
      flst_real(   6,  24)=          -1
      flst_real(   7,  24)=           5

      flst_real(   1,  25)=          -1
      flst_real(   2,  25)=           4
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=         -11
      flst_real(   5,  25)=          12
      flst_real(   6,  25)=           2
      flst_real(   7,  25)=          -2

      flst_real(   1,  26)=          -1
      flst_real(   2,  26)=           4
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=         -11
      flst_real(   5,  26)=          12
      flst_real(   6,  26)=          -2
      flst_real(   7,  26)=           4

      flst_real(   1,  27)=          -1
      flst_real(   2,  27)=           4
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=         -11
      flst_real(   5,  27)=          12
      flst_real(   6,  27)=           4
      flst_real(   7,  27)=          -4

      flst_real(   1,  28)=          -1
      flst_real(   2,  28)=           4
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=         -11
      flst_real(   5,  28)=          12
      flst_real(   6,  28)=           3
      flst_real(   7,  28)=          -3

      flst_real(   1,  29)=          -1
      flst_real(   2,  29)=           4
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=         -11
      flst_real(   5,  29)=          12
      flst_real(   6,  29)=           5
      flst_real(   7,  29)=          -5

      flst_real(   1,  30)=          -1
      flst_real(   2,  30)=           4
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=         -11
      flst_real(   5,  30)=          12
      flst_real(   6,  30)=           0
      flst_real(   7,  30)=           0

      flst_real(   1,  31)=          -1
      flst_real(   2,  31)=          -3
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=         -11
      flst_real(   5,  31)=          12
      flst_real(   6,  31)=          -1
      flst_real(   7,  31)=          -2

      flst_real(   1,  32)=          -1
      flst_real(   2,  32)=          -3
      flst_real(   3,  32)=          25
      flst_real(   4,  32)=         -11
      flst_real(   5,  32)=          12
      flst_real(   6,  32)=          -1
      flst_real(   7,  32)=          -4

      flst_real(   1,  33)=          -1
      flst_real(   2,  33)=          -3
      flst_real(   3,  33)=          25
      flst_real(   4,  33)=         -11
      flst_real(   5,  33)=          12
      flst_real(   6,  33)=          -2
      flst_real(   7,  33)=          -3

      flst_real(   1,  34)=          -1
      flst_real(   2,  34)=          -3
      flst_real(   3,  34)=          25
      flst_real(   4,  34)=         -11
      flst_real(   5,  34)=          12
      flst_real(   6,  34)=          -4
      flst_real(   7,  34)=          -3

      flst_real(   1,  35)=          -1
      flst_real(   2,  35)=           3
      flst_real(   3,  35)=          25
      flst_real(   4,  35)=         -11
      flst_real(   5,  35)=          12
      flst_real(   6,  35)=          -2
      flst_real(   7,  35)=           3

      flst_real(   1,  36)=          -1
      flst_real(   2,  36)=           3
      flst_real(   3,  36)=          25
      flst_real(   4,  36)=         -11
      flst_real(   5,  36)=          12
      flst_real(   6,  36)=          -4
      flst_real(   7,  36)=           3

      flst_real(   1,  37)=          -1
      flst_real(   2,  37)=          -5
      flst_real(   3,  37)=          25
      flst_real(   4,  37)=         -11
      flst_real(   5,  37)=          12
      flst_real(   6,  37)=          -1
      flst_real(   7,  37)=          -2

      flst_real(   1,  38)=          -1
      flst_real(   2,  38)=          -5
      flst_real(   3,  38)=          25
      flst_real(   4,  38)=         -11
      flst_real(   5,  38)=          12
      flst_real(   6,  38)=          -1
      flst_real(   7,  38)=          -4

      flst_real(   1,  39)=          -1
      flst_real(   2,  39)=          -5
      flst_real(   3,  39)=          25
      flst_real(   4,  39)=         -11
      flst_real(   5,  39)=          12
      flst_real(   6,  39)=          -2
      flst_real(   7,  39)=          -5

      flst_real(   1,  40)=          -1
      flst_real(   2,  40)=          -5
      flst_real(   3,  40)=          25
      flst_real(   4,  40)=         -11
      flst_real(   5,  40)=          12
      flst_real(   6,  40)=          -4
      flst_real(   7,  40)=          -5

      flst_real(   1,  41)=          -1
      flst_real(   2,  41)=           5
      flst_real(   3,  41)=          25
      flst_real(   4,  41)=         -11
      flst_real(   5,  41)=          12
      flst_real(   6,  41)=          -2
      flst_real(   7,  41)=           5

      flst_real(   1,  42)=          -1
      flst_real(   2,  42)=           5
      flst_real(   3,  42)=          25
      flst_real(   4,  42)=         -11
      flst_real(   5,  42)=          12
      flst_real(   6,  42)=          -4
      flst_real(   7,  42)=           5

      flst_real(   1,  43)=          -1
      flst_real(   2,  43)=           0
      flst_real(   3,  43)=          25
      flst_real(   4,  43)=         -11
      flst_real(   5,  43)=          12
      flst_real(   6,  43)=          -2
      flst_real(   7,  43)=           0

      flst_real(   1,  44)=          -1
      flst_real(   2,  44)=           0
      flst_real(   3,  44)=          25
      flst_real(   4,  44)=         -11
      flst_real(   5,  44)=          12
      flst_real(   6,  44)=          -4
      flst_real(   7,  44)=           0

      flst_real(   1,  45)=           1
      flst_real(   2,  45)=          -1
      flst_real(   3,  45)=          25
      flst_real(   4,  45)=         -11
      flst_real(   5,  45)=          12
      flst_real(   6,  45)=           1
      flst_real(   7,  45)=          -2

      flst_real(   1,  46)=           1
      flst_real(   2,  46)=          -1
      flst_real(   3,  46)=          25
      flst_real(   4,  46)=         -11
      flst_real(   5,  46)=          12
      flst_real(   6,  46)=           1
      flst_real(   7,  46)=          -4

      flst_real(   1,  47)=           1
      flst_real(   2,  47)=          -1
      flst_real(   3,  47)=          25
      flst_real(   4,  47)=         -11
      flst_real(   5,  47)=          12
      flst_real(   6,  47)=          -2
      flst_real(   7,  47)=           3

      flst_real(   1,  48)=           1
      flst_real(   2,  48)=          -1
      flst_real(   3,  48)=          25
      flst_real(   4,  48)=         -11
      flst_real(   5,  48)=          12
      flst_real(   6,  48)=          -2
      flst_real(   7,  48)=           5

      flst_real(   1,  49)=           1
      flst_real(   2,  49)=          -1
      flst_real(   3,  49)=          25
      flst_real(   4,  49)=         -11
      flst_real(   5,  49)=          12
      flst_real(   6,  49)=          -4
      flst_real(   7,  49)=           3

      flst_real(   1,  50)=           1
      flst_real(   2,  50)=          -1
      flst_real(   3,  50)=          25
      flst_real(   4,  50)=         -11
      flst_real(   5,  50)=          12
      flst_real(   6,  50)=          -4
      flst_real(   7,  50)=           5

      flst_real(   1,  51)=           1
      flst_real(   2,  51)=           2
      flst_real(   3,  51)=          25
      flst_real(   4,  51)=         -11
      flst_real(   5,  51)=          12
      flst_real(   6,  51)=           1
      flst_real(   7,  51)=           1

      flst_real(   1,  52)=           1
      flst_real(   2,  52)=           2
      flst_real(   3,  52)=          25
      flst_real(   4,  52)=         -11
      flst_real(   5,  52)=          12
      flst_real(   6,  52)=           1
      flst_real(   7,  52)=           3

      flst_real(   1,  53)=           1
      flst_real(   2,  53)=           2
      flst_real(   3,  53)=          25
      flst_real(   4,  53)=         -11
      flst_real(   5,  53)=          12
      flst_real(   6,  53)=           1
      flst_real(   7,  53)=           5

      flst_real(   1,  54)=           1
      flst_real(   2,  54)=           4
      flst_real(   3,  54)=          25
      flst_real(   4,  54)=         -11
      flst_real(   5,  54)=          12
      flst_real(   6,  54)=           1
      flst_real(   7,  54)=           1

      flst_real(   1,  55)=           1
      flst_real(   2,  55)=           4
      flst_real(   3,  55)=          25
      flst_real(   4,  55)=         -11
      flst_real(   5,  55)=          12
      flst_real(   6,  55)=           1
      flst_real(   7,  55)=           3

      flst_real(   1,  56)=           1
      flst_real(   2,  56)=           4
      flst_real(   3,  56)=          25
      flst_real(   4,  56)=         -11
      flst_real(   5,  56)=          12
      flst_real(   6,  56)=           1
      flst_real(   7,  56)=           5

      flst_real(   1,  57)=           1
      flst_real(   2,  57)=          -3
      flst_real(   3,  57)=          25
      flst_real(   4,  57)=         -11
      flst_real(   5,  57)=          12
      flst_real(   6,  57)=           1
      flst_real(   7,  57)=          -2

      flst_real(   1,  58)=           1
      flst_real(   2,  58)=          -3
      flst_real(   3,  58)=          25
      flst_real(   4,  58)=         -11
      flst_real(   5,  58)=          12
      flst_real(   6,  58)=           1
      flst_real(   7,  58)=          -4

      flst_real(   1,  59)=           1
      flst_real(   2,  59)=          -5
      flst_real(   3,  59)=          25
      flst_real(   4,  59)=         -11
      flst_real(   5,  59)=          12
      flst_real(   6,  59)=           1
      flst_real(   7,  59)=          -2

      flst_real(   1,  60)=           1
      flst_real(   2,  60)=          -5
      flst_real(   3,  60)=          25
      flst_real(   4,  60)=         -11
      flst_real(   5,  60)=          12
      flst_real(   6,  60)=           1
      flst_real(   7,  60)=          -4

      flst_real(   1,  61)=          -2
      flst_real(   2,  61)=          -1
      flst_real(   3,  61)=          25
      flst_real(   4,  61)=         -11
      flst_real(   5,  61)=          12
      flst_real(   6,  61)=          -2
      flst_real(   7,  61)=          -2

      flst_real(   1,  62)=          -2
      flst_real(   2,  62)=          -1
      flst_real(   3,  62)=          25
      flst_real(   4,  62)=         -11
      flst_real(   5,  62)=          12
      flst_real(   6,  62)=          -2
      flst_real(   7,  62)=          -4

      flst_real(   1,  63)=          -2
      flst_real(   2,  63)=           2
      flst_real(   3,  63)=          25
      flst_real(   4,  63)=         -11
      flst_real(   5,  63)=          12
      flst_real(   6,  63)=           1
      flst_real(   7,  63)=          -2

      flst_real(   1,  64)=          -2
      flst_real(   2,  64)=           2
      flst_real(   3,  64)=          25
      flst_real(   4,  64)=         -11
      flst_real(   5,  64)=          12
      flst_real(   6,  64)=           1
      flst_real(   7,  64)=          -4

      flst_real(   1,  65)=          -2
      flst_real(   2,  65)=           2
      flst_real(   3,  65)=          25
      flst_real(   4,  65)=         -11
      flst_real(   5,  65)=          12
      flst_real(   6,  65)=          -2
      flst_real(   7,  65)=           3

      flst_real(   1,  66)=          -2
      flst_real(   2,  66)=           2
      flst_real(   3,  66)=          25
      flst_real(   4,  66)=         -11
      flst_real(   5,  66)=          12
      flst_real(   6,  66)=          -2
      flst_real(   7,  66)=           5

      flst_real(   1,  67)=          -2
      flst_real(   2,  67)=           2
      flst_real(   3,  67)=          25
      flst_real(   4,  67)=         -11
      flst_real(   5,  67)=          12
      flst_real(   6,  67)=          -4
      flst_real(   7,  67)=           3

      flst_real(   1,  68)=          -2
      flst_real(   2,  68)=           2
      flst_real(   3,  68)=          25
      flst_real(   4,  68)=         -11
      flst_real(   5,  68)=          12
      flst_real(   6,  68)=          -4
      flst_real(   7,  68)=           5

      flst_real(   1,  69)=          -2
      flst_real(   2,  69)=           4
      flst_real(   3,  69)=          25
      flst_real(   4,  69)=         -11
      flst_real(   5,  69)=          12
      flst_real(   6,  69)=           1
      flst_real(   7,  69)=          -2

      flst_real(   1,  70)=          -2
      flst_real(   2,  70)=           4
      flst_real(   3,  70)=          25
      flst_real(   4,  70)=         -11
      flst_real(   5,  70)=          12
      flst_real(   6,  70)=          -2
      flst_real(   7,  70)=           3

      flst_real(   1,  71)=          -2
      flst_real(   2,  71)=           4
      flst_real(   3,  71)=          25
      flst_real(   4,  71)=         -11
      flst_real(   5,  71)=          12
      flst_real(   6,  71)=          -2
      flst_real(   7,  71)=           5

      flst_real(   1,  72)=          -2
      flst_real(   2,  72)=          -3
      flst_real(   3,  72)=          25
      flst_real(   4,  72)=         -11
      flst_real(   5,  72)=          12
      flst_real(   6,  72)=          -2
      flst_real(   7,  72)=          -2

      flst_real(   1,  73)=          -2
      flst_real(   2,  73)=          -3
      flst_real(   3,  73)=          25
      flst_real(   4,  73)=         -11
      flst_real(   5,  73)=          12
      flst_real(   6,  73)=          -2
      flst_real(   7,  73)=          -4

      flst_real(   1,  74)=          -2
      flst_real(   2,  74)=          -5
      flst_real(   3,  74)=          25
      flst_real(   4,  74)=         -11
      flst_real(   5,  74)=          12
      flst_real(   6,  74)=          -2
      flst_real(   7,  74)=          -2

      flst_real(   1,  75)=          -2
      flst_real(   2,  75)=          -5
      flst_real(   3,  75)=          25
      flst_real(   4,  75)=         -11
      flst_real(   5,  75)=          12
      flst_real(   6,  75)=          -2
      flst_real(   7,  75)=          -4

      flst_real(   1,  76)=           2
      flst_real(   2,  76)=          -1
      flst_real(   3,  76)=          25
      flst_real(   4,  76)=         -11
      flst_real(   5,  76)=          12
      flst_real(   6,  76)=           1
      flst_real(   7,  76)=          -1

      flst_real(   1,  77)=           2
      flst_real(   2,  77)=          -1
      flst_real(   3,  77)=          25
      flst_real(   4,  77)=         -11
      flst_real(   5,  77)=          12
      flst_real(   6,  77)=          -1
      flst_real(   7,  77)=           3

      flst_real(   1,  78)=           2
      flst_real(   2,  78)=          -1
      flst_real(   3,  78)=          25
      flst_real(   4,  78)=         -11
      flst_real(   5,  78)=          12
      flst_real(   6,  78)=          -1
      flst_real(   7,  78)=           5

      flst_real(   1,  79)=           2
      flst_real(   2,  79)=          -1
      flst_real(   3,  79)=          25
      flst_real(   4,  79)=         -11
      flst_real(   5,  79)=          12
      flst_real(   6,  79)=           2
      flst_real(   7,  79)=          -2

      flst_real(   1,  80)=           2
      flst_real(   2,  80)=          -1
      flst_real(   3,  80)=          25
      flst_real(   4,  80)=         -11
      flst_real(   5,  80)=          12
      flst_real(   6,  80)=           2
      flst_real(   7,  80)=          -4

      flst_real(   1,  81)=           2
      flst_real(   2,  81)=          -1
      flst_real(   3,  81)=          25
      flst_real(   4,  81)=         -11
      flst_real(   5,  81)=          12
      flst_real(   6,  81)=           4
      flst_real(   7,  81)=          -4

      flst_real(   1,  82)=           2
      flst_real(   2,  82)=          -1
      flst_real(   3,  82)=          25
      flst_real(   4,  82)=         -11
      flst_real(   5,  82)=          12
      flst_real(   6,  82)=           3
      flst_real(   7,  82)=          -3

      flst_real(   1,  83)=           2
      flst_real(   2,  83)=          -1
      flst_real(   3,  83)=          25
      flst_real(   4,  83)=         -11
      flst_real(   5,  83)=          12
      flst_real(   6,  83)=           5
      flst_real(   7,  83)=          -5

      flst_real(   1,  84)=           2
      flst_real(   2,  84)=          -1
      flst_real(   3,  84)=          25
      flst_real(   4,  84)=         -11
      flst_real(   5,  84)=          12
      flst_real(   6,  84)=           0
      flst_real(   7,  84)=           0

      flst_real(   1,  85)=           2
      flst_real(   2,  85)=           1
      flst_real(   3,  85)=          25
      flst_real(   4,  85)=         -11
      flst_real(   5,  85)=          12
      flst_real(   6,  85)=           1
      flst_real(   7,  85)=           1

      flst_real(   1,  86)=           2
      flst_real(   2,  86)=           1
      flst_real(   3,  86)=          25
      flst_real(   4,  86)=         -11
      flst_real(   5,  86)=          12
      flst_real(   6,  86)=           1
      flst_real(   7,  86)=           3

      flst_real(   1,  87)=           2
      flst_real(   2,  87)=           1
      flst_real(   3,  87)=          25
      flst_real(   4,  87)=         -11
      flst_real(   5,  87)=          12
      flst_real(   6,  87)=           1
      flst_real(   7,  87)=           5

      flst_real(   1,  88)=           2
      flst_real(   2,  88)=          -2
      flst_real(   3,  88)=          25
      flst_real(   4,  88)=         -11
      flst_real(   5,  88)=          12
      flst_real(   6,  88)=           1
      flst_real(   7,  88)=          -2

      flst_real(   1,  89)=           2
      flst_real(   2,  89)=          -2
      flst_real(   3,  89)=          25
      flst_real(   4,  89)=         -11
      flst_real(   5,  89)=          12
      flst_real(   6,  89)=           1
      flst_real(   7,  89)=          -4

      flst_real(   1,  90)=           2
      flst_real(   2,  90)=          -2
      flst_real(   3,  90)=          25
      flst_real(   4,  90)=         -11
      flst_real(   5,  90)=          12
      flst_real(   6,  90)=          -2
      flst_real(   7,  90)=           3

      flst_real(   1,  91)=           2
      flst_real(   2,  91)=          -2
      flst_real(   3,  91)=          25
      flst_real(   4,  91)=         -11
      flst_real(   5,  91)=          12
      flst_real(   6,  91)=          -2
      flst_real(   7,  91)=           5

      flst_real(   1,  92)=           2
      flst_real(   2,  92)=          -2
      flst_real(   3,  92)=          25
      flst_real(   4,  92)=         -11
      flst_real(   5,  92)=          12
      flst_real(   6,  92)=          -4
      flst_real(   7,  92)=           3

      flst_real(   1,  93)=           2
      flst_real(   2,  93)=          -2
      flst_real(   3,  93)=          25
      flst_real(   4,  93)=         -11
      flst_real(   5,  93)=          12
      flst_real(   6,  93)=          -4
      flst_real(   7,  93)=           5

      flst_real(   1,  94)=           2
      flst_real(   2,  94)=           2
      flst_real(   3,  94)=          25
      flst_real(   4,  94)=         -11
      flst_real(   5,  94)=          12
      flst_real(   6,  94)=           1
      flst_real(   7,  94)=           2

      flst_real(   1,  95)=           2
      flst_real(   2,  95)=           2
      flst_real(   3,  95)=          25
      flst_real(   4,  95)=         -11
      flst_real(   5,  95)=          12
      flst_real(   6,  95)=           2
      flst_real(   7,  95)=           3

      flst_real(   1,  96)=           2
      flst_real(   2,  96)=           2
      flst_real(   3,  96)=          25
      flst_real(   4,  96)=         -11
      flst_real(   5,  96)=          12
      flst_real(   6,  96)=           2
      flst_real(   7,  96)=           5

      flst_real(   1,  97)=           2
      flst_real(   2,  97)=          -4
      flst_real(   3,  97)=          25
      flst_real(   4,  97)=         -11
      flst_real(   5,  97)=          12
      flst_real(   6,  97)=           1
      flst_real(   7,  97)=          -4

      flst_real(   1,  98)=           2
      flst_real(   2,  98)=          -4
      flst_real(   3,  98)=          25
      flst_real(   4,  98)=         -11
      flst_real(   5,  98)=          12
      flst_real(   6,  98)=          -4
      flst_real(   7,  98)=           3

      flst_real(   1,  99)=           2
      flst_real(   2,  99)=          -4
      flst_real(   3,  99)=          25
      flst_real(   4,  99)=         -11
      flst_real(   5,  99)=          12
      flst_real(   6,  99)=          -4
      flst_real(   7,  99)=           5

      flst_real(   1, 100)=           2
      flst_real(   2, 100)=           4
      flst_real(   3, 100)=          25
      flst_real(   4, 100)=         -11
      flst_real(   5, 100)=          12
      flst_real(   6, 100)=           1
      flst_real(   7, 100)=           2

      flst_real(   1, 101)=           2
      flst_real(   2, 101)=           4
      flst_real(   3, 101)=          25
      flst_real(   4, 101)=         -11
      flst_real(   5, 101)=          12
      flst_real(   6, 101)=           1
      flst_real(   7, 101)=           4

      flst_real(   1, 102)=           2
      flst_real(   2, 102)=           4
      flst_real(   3, 102)=          25
      flst_real(   4, 102)=         -11
      flst_real(   5, 102)=          12
      flst_real(   6, 102)=           2
      flst_real(   7, 102)=           3

      flst_real(   1, 103)=           2
      flst_real(   2, 103)=           4
      flst_real(   3, 103)=          25
      flst_real(   4, 103)=         -11
      flst_real(   5, 103)=          12
      flst_real(   6, 103)=           2
      flst_real(   7, 103)=           5

      flst_real(   1, 104)=           2
      flst_real(   2, 104)=           4
      flst_real(   3, 104)=          25
      flst_real(   4, 104)=         -11
      flst_real(   5, 104)=          12
      flst_real(   6, 104)=           4
      flst_real(   7, 104)=           3

      flst_real(   1, 105)=           2
      flst_real(   2, 105)=           4
      flst_real(   3, 105)=          25
      flst_real(   4, 105)=         -11
      flst_real(   5, 105)=          12
      flst_real(   6, 105)=           4
      flst_real(   7, 105)=           5

      flst_real(   1, 106)=           2
      flst_real(   2, 106)=          -3
      flst_real(   3, 106)=          25
      flst_real(   4, 106)=         -11
      flst_real(   5, 106)=          12
      flst_real(   6, 106)=           1
      flst_real(   7, 106)=          -1

      flst_real(   1, 107)=           2
      flst_real(   2, 107)=          -3
      flst_real(   3, 107)=          25
      flst_real(   4, 107)=         -11
      flst_real(   5, 107)=          12
      flst_real(   6, 107)=           1
      flst_real(   7, 107)=          -3

      flst_real(   1, 108)=           2
      flst_real(   2, 108)=          -3
      flst_real(   3, 108)=          25
      flst_real(   4, 108)=         -11
      flst_real(   5, 108)=          12
      flst_real(   6, 108)=           2
      flst_real(   7, 108)=          -2

      flst_real(   1, 109)=           2
      flst_real(   2, 109)=          -3
      flst_real(   3, 109)=          25
      flst_real(   4, 109)=         -11
      flst_real(   5, 109)=          12
      flst_real(   6, 109)=           2
      flst_real(   7, 109)=          -4

      flst_real(   1, 110)=           2
      flst_real(   2, 110)=          -3
      flst_real(   3, 110)=          25
      flst_real(   4, 110)=         -11
      flst_real(   5, 110)=          12
      flst_real(   6, 110)=           4
      flst_real(   7, 110)=          -4

      flst_real(   1, 111)=           2
      flst_real(   2, 111)=          -3
      flst_real(   3, 111)=          25
      flst_real(   4, 111)=         -11
      flst_real(   5, 111)=          12
      flst_real(   6, 111)=           3
      flst_real(   7, 111)=          -3

      flst_real(   1, 112)=           2
      flst_real(   2, 112)=          -3
      flst_real(   3, 112)=          25
      flst_real(   4, 112)=         -11
      flst_real(   5, 112)=          12
      flst_real(   6, 112)=          -3
      flst_real(   7, 112)=           5

      flst_real(   1, 113)=           2
      flst_real(   2, 113)=          -3
      flst_real(   3, 113)=          25
      flst_real(   4, 113)=         -11
      flst_real(   5, 113)=          12
      flst_real(   6, 113)=           5
      flst_real(   7, 113)=          -5

      flst_real(   1, 114)=           2
      flst_real(   2, 114)=          -3
      flst_real(   3, 114)=          25
      flst_real(   4, 114)=         -11
      flst_real(   5, 114)=          12
      flst_real(   6, 114)=           0
      flst_real(   7, 114)=           0

      flst_real(   1, 115)=           2
      flst_real(   2, 115)=           3
      flst_real(   3, 115)=          25
      flst_real(   4, 115)=         -11
      flst_real(   5, 115)=          12
      flst_real(   6, 115)=           1
      flst_real(   7, 115)=           3

      flst_real(   1, 116)=           2
      flst_real(   2, 116)=           3
      flst_real(   3, 116)=          25
      flst_real(   4, 116)=         -11
      flst_real(   5, 116)=          12
      flst_real(   6, 116)=           3
      flst_real(   7, 116)=           3

      flst_real(   1, 117)=           2
      flst_real(   2, 117)=           3
      flst_real(   3, 117)=          25
      flst_real(   4, 117)=         -11
      flst_real(   5, 117)=          12
      flst_real(   6, 117)=           3
      flst_real(   7, 117)=           5

      flst_real(   1, 118)=           2
      flst_real(   2, 118)=          -5
      flst_real(   3, 118)=          25
      flst_real(   4, 118)=         -11
      flst_real(   5, 118)=          12
      flst_real(   6, 118)=           1
      flst_real(   7, 118)=          -1

      flst_real(   1, 119)=           2
      flst_real(   2, 119)=          -5
      flst_real(   3, 119)=          25
      flst_real(   4, 119)=         -11
      flst_real(   5, 119)=          12
      flst_real(   6, 119)=           1
      flst_real(   7, 119)=          -5

      flst_real(   1, 120)=           2
      flst_real(   2, 120)=          -5
      flst_real(   3, 120)=          25
      flst_real(   4, 120)=         -11
      flst_real(   5, 120)=          12
      flst_real(   6, 120)=           2
      flst_real(   7, 120)=          -2

      flst_real(   1, 121)=           2
      flst_real(   2, 121)=          -5
      flst_real(   3, 121)=          25
      flst_real(   4, 121)=         -11
      flst_real(   5, 121)=          12
      flst_real(   6, 121)=           2
      flst_real(   7, 121)=          -4

      flst_real(   1, 122)=           2
      flst_real(   2, 122)=          -5
      flst_real(   3, 122)=          25
      flst_real(   4, 122)=         -11
      flst_real(   5, 122)=          12
      flst_real(   6, 122)=           4
      flst_real(   7, 122)=          -4

      flst_real(   1, 123)=           2
      flst_real(   2, 123)=          -5
      flst_real(   3, 123)=          25
      flst_real(   4, 123)=         -11
      flst_real(   5, 123)=          12
      flst_real(   6, 123)=           3
      flst_real(   7, 123)=          -3

      flst_real(   1, 124)=           2
      flst_real(   2, 124)=          -5
      flst_real(   3, 124)=          25
      flst_real(   4, 124)=         -11
      flst_real(   5, 124)=          12
      flst_real(   6, 124)=           3
      flst_real(   7, 124)=          -5

      flst_real(   1, 125)=           2
      flst_real(   2, 125)=          -5
      flst_real(   3, 125)=          25
      flst_real(   4, 125)=         -11
      flst_real(   5, 125)=          12
      flst_real(   6, 125)=           5
      flst_real(   7, 125)=          -5

      flst_real(   1, 126)=           2
      flst_real(   2, 126)=          -5
      flst_real(   3, 126)=          25
      flst_real(   4, 126)=         -11
      flst_real(   5, 126)=          12
      flst_real(   6, 126)=           0
      flst_real(   7, 126)=           0

      flst_real(   1, 127)=           2
      flst_real(   2, 127)=           5
      flst_real(   3, 127)=          25
      flst_real(   4, 127)=         -11
      flst_real(   5, 127)=          12
      flst_real(   6, 127)=           1
      flst_real(   7, 127)=           5

      flst_real(   1, 128)=           2
      flst_real(   2, 128)=           5
      flst_real(   3, 128)=          25
      flst_real(   4, 128)=         -11
      flst_real(   5, 128)=          12
      flst_real(   6, 128)=           3
      flst_real(   7, 128)=           5

      flst_real(   1, 129)=           2
      flst_real(   2, 129)=           5
      flst_real(   3, 129)=          25
      flst_real(   4, 129)=         -11
      flst_real(   5, 129)=          12
      flst_real(   6, 129)=           5
      flst_real(   7, 129)=           5

      flst_real(   1, 130)=           2
      flst_real(   2, 130)=           0
      flst_real(   3, 130)=          25
      flst_real(   4, 130)=         -11
      flst_real(   5, 130)=          12
      flst_real(   6, 130)=           1
      flst_real(   7, 130)=           0

      flst_real(   1, 131)=           2
      flst_real(   2, 131)=           0
      flst_real(   3, 131)=          25
      flst_real(   4, 131)=         -11
      flst_real(   5, 131)=          12
      flst_real(   6, 131)=           3
      flst_real(   7, 131)=           0

      flst_real(   1, 132)=           2
      flst_real(   2, 132)=           0
      flst_real(   3, 132)=          25
      flst_real(   4, 132)=         -11
      flst_real(   5, 132)=          12
      flst_real(   6, 132)=           5
      flst_real(   7, 132)=           0

      flst_real(   1, 133)=          -4
      flst_real(   2, 133)=          -1
      flst_real(   3, 133)=          25
      flst_real(   4, 133)=         -11
      flst_real(   5, 133)=          12
      flst_real(   6, 133)=          -2
      flst_real(   7, 133)=          -4

      flst_real(   1, 134)=          -4
      flst_real(   2, 134)=          -1
      flst_real(   3, 134)=          25
      flst_real(   4, 134)=         -11
      flst_real(   5, 134)=          12
      flst_real(   6, 134)=          -4
      flst_real(   7, 134)=          -4

      flst_real(   1, 135)=          -4
      flst_real(   2, 135)=           2
      flst_real(   3, 135)=          25
      flst_real(   4, 135)=         -11
      flst_real(   5, 135)=          12
      flst_real(   6, 135)=           1
      flst_real(   7, 135)=          -4

      flst_real(   1, 136)=          -4
      flst_real(   2, 136)=           2
      flst_real(   3, 136)=          25
      flst_real(   4, 136)=         -11
      flst_real(   5, 136)=          12
      flst_real(   6, 136)=          -4
      flst_real(   7, 136)=           3

      flst_real(   1, 137)=          -4
      flst_real(   2, 137)=           2
      flst_real(   3, 137)=          25
      flst_real(   4, 137)=         -11
      flst_real(   5, 137)=          12
      flst_real(   6, 137)=          -4
      flst_real(   7, 137)=           5

      flst_real(   1, 138)=          -4
      flst_real(   2, 138)=           4
      flst_real(   3, 138)=          25
      flst_real(   4, 138)=         -11
      flst_real(   5, 138)=          12
      flst_real(   6, 138)=           1
      flst_real(   7, 138)=          -2

      flst_real(   1, 139)=          -4
      flst_real(   2, 139)=           4
      flst_real(   3, 139)=          25
      flst_real(   4, 139)=         -11
      flst_real(   5, 139)=          12
      flst_real(   6, 139)=           1
      flst_real(   7, 139)=          -4

      flst_real(   1, 140)=          -4
      flst_real(   2, 140)=           4
      flst_real(   3, 140)=          25
      flst_real(   4, 140)=         -11
      flst_real(   5, 140)=          12
      flst_real(   6, 140)=          -2
      flst_real(   7, 140)=           3

      flst_real(   1, 141)=          -4
      flst_real(   2, 141)=           4
      flst_real(   3, 141)=          25
      flst_real(   4, 141)=         -11
      flst_real(   5, 141)=          12
      flst_real(   6, 141)=          -2
      flst_real(   7, 141)=           5

      flst_real(   1, 142)=          -4
      flst_real(   2, 142)=           4
      flst_real(   3, 142)=          25
      flst_real(   4, 142)=         -11
      flst_real(   5, 142)=          12
      flst_real(   6, 142)=          -4
      flst_real(   7, 142)=           3

      flst_real(   1, 143)=          -4
      flst_real(   2, 143)=           4
      flst_real(   3, 143)=          25
      flst_real(   4, 143)=         -11
      flst_real(   5, 143)=          12
      flst_real(   6, 143)=          -4
      flst_real(   7, 143)=           5

      flst_real(   1, 144)=          -4
      flst_real(   2, 144)=          -3
      flst_real(   3, 144)=          25
      flst_real(   4, 144)=         -11
      flst_real(   5, 144)=          12
      flst_real(   6, 144)=          -2
      flst_real(   7, 144)=          -4

      flst_real(   1, 145)=          -4
      flst_real(   2, 145)=          -3
      flst_real(   3, 145)=          25
      flst_real(   4, 145)=         -11
      flst_real(   5, 145)=          12
      flst_real(   6, 145)=          -4
      flst_real(   7, 145)=          -4

      flst_real(   1, 146)=          -4
      flst_real(   2, 146)=          -5
      flst_real(   3, 146)=          25
      flst_real(   4, 146)=         -11
      flst_real(   5, 146)=          12
      flst_real(   6, 146)=          -2
      flst_real(   7, 146)=          -4

      flst_real(   1, 147)=          -4
      flst_real(   2, 147)=          -5
      flst_real(   3, 147)=          25
      flst_real(   4, 147)=         -11
      flst_real(   5, 147)=          12
      flst_real(   6, 147)=          -4
      flst_real(   7, 147)=          -4

      flst_real(   1, 148)=           4
      flst_real(   2, 148)=          -1
      flst_real(   3, 148)=          25
      flst_real(   4, 148)=         -11
      flst_real(   5, 148)=          12
      flst_real(   6, 148)=           1
      flst_real(   7, 148)=          -1

      flst_real(   1, 149)=           4
      flst_real(   2, 149)=          -1
      flst_real(   3, 149)=          25
      flst_real(   4, 149)=         -11
      flst_real(   5, 149)=          12
      flst_real(   6, 149)=          -1
      flst_real(   7, 149)=           3

      flst_real(   1, 150)=           4
      flst_real(   2, 150)=          -1
      flst_real(   3, 150)=          25
      flst_real(   4, 150)=         -11
      flst_real(   5, 150)=          12
      flst_real(   6, 150)=          -1
      flst_real(   7, 150)=           5

      flst_real(   1, 151)=           4
      flst_real(   2, 151)=          -1
      flst_real(   3, 151)=          25
      flst_real(   4, 151)=         -11
      flst_real(   5, 151)=          12
      flst_real(   6, 151)=           2
      flst_real(   7, 151)=          -2

      flst_real(   1, 152)=           4
      flst_real(   2, 152)=          -1
      flst_real(   3, 152)=          25
      flst_real(   4, 152)=         -11
      flst_real(   5, 152)=          12
      flst_real(   6, 152)=          -2
      flst_real(   7, 152)=           4

      flst_real(   1, 153)=           4
      flst_real(   2, 153)=          -1
      flst_real(   3, 153)=          25
      flst_real(   4, 153)=         -11
      flst_real(   5, 153)=          12
      flst_real(   6, 153)=           4
      flst_real(   7, 153)=          -4

      flst_real(   1, 154)=           4
      flst_real(   2, 154)=          -1
      flst_real(   3, 154)=          25
      flst_real(   4, 154)=         -11
      flst_real(   5, 154)=          12
      flst_real(   6, 154)=           3
      flst_real(   7, 154)=          -3

      flst_real(   1, 155)=           4
      flst_real(   2, 155)=          -1
      flst_real(   3, 155)=          25
      flst_real(   4, 155)=         -11
      flst_real(   5, 155)=          12
      flst_real(   6, 155)=           5
      flst_real(   7, 155)=          -5

      flst_real(   1, 156)=           4
      flst_real(   2, 156)=          -1
      flst_real(   3, 156)=          25
      flst_real(   4, 156)=         -11
      flst_real(   5, 156)=          12
      flst_real(   6, 156)=           0
      flst_real(   7, 156)=           0

      flst_real(   1, 157)=           4
      flst_real(   2, 157)=           1
      flst_real(   3, 157)=          25
      flst_real(   4, 157)=         -11
      flst_real(   5, 157)=          12
      flst_real(   6, 157)=           1
      flst_real(   7, 157)=           1

      flst_real(   1, 158)=           4
      flst_real(   2, 158)=           1
      flst_real(   3, 158)=          25
      flst_real(   4, 158)=         -11
      flst_real(   5, 158)=          12
      flst_real(   6, 158)=           1
      flst_real(   7, 158)=           3

      flst_real(   1, 159)=           4
      flst_real(   2, 159)=           1
      flst_real(   3, 159)=          25
      flst_real(   4, 159)=         -11
      flst_real(   5, 159)=          12
      flst_real(   6, 159)=           1
      flst_real(   7, 159)=           5

      flst_real(   1, 160)=           4
      flst_real(   2, 160)=          -2
      flst_real(   3, 160)=          25
      flst_real(   4, 160)=         -11
      flst_real(   5, 160)=          12
      flst_real(   6, 160)=           1
      flst_real(   7, 160)=          -2

      flst_real(   1, 161)=           4
      flst_real(   2, 161)=          -2
      flst_real(   3, 161)=          25
      flst_real(   4, 161)=         -11
      flst_real(   5, 161)=          12
      flst_real(   6, 161)=          -2
      flst_real(   7, 161)=           3

      flst_real(   1, 162)=           4
      flst_real(   2, 162)=          -2
      flst_real(   3, 162)=          25
      flst_real(   4, 162)=         -11
      flst_real(   5, 162)=          12
      flst_real(   6, 162)=          -2
      flst_real(   7, 162)=           5

      flst_real(   1, 163)=           4
      flst_real(   2, 163)=           2
      flst_real(   3, 163)=          25
      flst_real(   4, 163)=         -11
      flst_real(   5, 163)=          12
      flst_real(   6, 163)=           1
      flst_real(   7, 163)=           2

      flst_real(   1, 164)=           4
      flst_real(   2, 164)=           2
      flst_real(   3, 164)=          25
      flst_real(   4, 164)=         -11
      flst_real(   5, 164)=          12
      flst_real(   6, 164)=           1
      flst_real(   7, 164)=           4

      flst_real(   1, 165)=           4
      flst_real(   2, 165)=           2
      flst_real(   3, 165)=          25
      flst_real(   4, 165)=         -11
      flst_real(   5, 165)=          12
      flst_real(   6, 165)=           2
      flst_real(   7, 165)=           3

      flst_real(   1, 166)=           4
      flst_real(   2, 166)=           2
      flst_real(   3, 166)=          25
      flst_real(   4, 166)=         -11
      flst_real(   5, 166)=          12
      flst_real(   6, 166)=           2
      flst_real(   7, 166)=           5

      flst_real(   1, 167)=           4
      flst_real(   2, 167)=           2
      flst_real(   3, 167)=          25
      flst_real(   4, 167)=         -11
      flst_real(   5, 167)=          12
      flst_real(   6, 167)=           4
      flst_real(   7, 167)=           3

      flst_real(   1, 168)=           4
      flst_real(   2, 168)=           2
      flst_real(   3, 168)=          25
      flst_real(   4, 168)=         -11
      flst_real(   5, 168)=          12
      flst_real(   6, 168)=           4
      flst_real(   7, 168)=           5

      flst_real(   1, 169)=           4
      flst_real(   2, 169)=          -4
      flst_real(   3, 169)=          25
      flst_real(   4, 169)=         -11
      flst_real(   5, 169)=          12
      flst_real(   6, 169)=           1
      flst_real(   7, 169)=          -2

      flst_real(   1, 170)=           4
      flst_real(   2, 170)=          -4
      flst_real(   3, 170)=          25
      flst_real(   4, 170)=         -11
      flst_real(   5, 170)=          12
      flst_real(   6, 170)=           1
      flst_real(   7, 170)=          -4

      flst_real(   1, 171)=           4
      flst_real(   2, 171)=          -4
      flst_real(   3, 171)=          25
      flst_real(   4, 171)=         -11
      flst_real(   5, 171)=          12
      flst_real(   6, 171)=          -2
      flst_real(   7, 171)=           3

      flst_real(   1, 172)=           4
      flst_real(   2, 172)=          -4
      flst_real(   3, 172)=          25
      flst_real(   4, 172)=         -11
      flst_real(   5, 172)=          12
      flst_real(   6, 172)=          -2
      flst_real(   7, 172)=           5

      flst_real(   1, 173)=           4
      flst_real(   2, 173)=          -4
      flst_real(   3, 173)=          25
      flst_real(   4, 173)=         -11
      flst_real(   5, 173)=          12
      flst_real(   6, 173)=          -4
      flst_real(   7, 173)=           3

      flst_real(   1, 174)=           4
      flst_real(   2, 174)=          -4
      flst_real(   3, 174)=          25
      flst_real(   4, 174)=         -11
      flst_real(   5, 174)=          12
      flst_real(   6, 174)=          -4
      flst_real(   7, 174)=           5

      flst_real(   1, 175)=           4
      flst_real(   2, 175)=           4
      flst_real(   3, 175)=          25
      flst_real(   4, 175)=         -11
      flst_real(   5, 175)=          12
      flst_real(   6, 175)=           1
      flst_real(   7, 175)=           4

      flst_real(   1, 176)=           4
      flst_real(   2, 176)=           4
      flst_real(   3, 176)=          25
      flst_real(   4, 176)=         -11
      flst_real(   5, 176)=          12
      flst_real(   6, 176)=           4
      flst_real(   7, 176)=           3

      flst_real(   1, 177)=           4
      flst_real(   2, 177)=           4
      flst_real(   3, 177)=          25
      flst_real(   4, 177)=         -11
      flst_real(   5, 177)=          12
      flst_real(   6, 177)=           4
      flst_real(   7, 177)=           5

      flst_real(   1, 178)=           4
      flst_real(   2, 178)=          -3
      flst_real(   3, 178)=          25
      flst_real(   4, 178)=         -11
      flst_real(   5, 178)=          12
      flst_real(   6, 178)=           1
      flst_real(   7, 178)=          -1

      flst_real(   1, 179)=           4
      flst_real(   2, 179)=          -3
      flst_real(   3, 179)=          25
      flst_real(   4, 179)=         -11
      flst_real(   5, 179)=          12
      flst_real(   6, 179)=           1
      flst_real(   7, 179)=          -3

      flst_real(   1, 180)=           4
      flst_real(   2, 180)=          -3
      flst_real(   3, 180)=          25
      flst_real(   4, 180)=         -11
      flst_real(   5, 180)=          12
      flst_real(   6, 180)=           2
      flst_real(   7, 180)=          -2

      flst_real(   1, 181)=           4
      flst_real(   2, 181)=          -3
      flst_real(   3, 181)=          25
      flst_real(   4, 181)=         -11
      flst_real(   5, 181)=          12
      flst_real(   6, 181)=          -2
      flst_real(   7, 181)=           4

      flst_real(   1, 182)=           4
      flst_real(   2, 182)=          -3
      flst_real(   3, 182)=          25
      flst_real(   4, 182)=         -11
      flst_real(   5, 182)=          12
      flst_real(   6, 182)=           4
      flst_real(   7, 182)=          -4

      flst_real(   1, 183)=           4
      flst_real(   2, 183)=          -3
      flst_real(   3, 183)=          25
      flst_real(   4, 183)=         -11
      flst_real(   5, 183)=          12
      flst_real(   6, 183)=           3
      flst_real(   7, 183)=          -3

      flst_real(   1, 184)=           4
      flst_real(   2, 184)=          -3
      flst_real(   3, 184)=          25
      flst_real(   4, 184)=         -11
      flst_real(   5, 184)=          12
      flst_real(   6, 184)=          -3
      flst_real(   7, 184)=           5

      flst_real(   1, 185)=           4
      flst_real(   2, 185)=          -3
      flst_real(   3, 185)=          25
      flst_real(   4, 185)=         -11
      flst_real(   5, 185)=          12
      flst_real(   6, 185)=           5
      flst_real(   7, 185)=          -5

      flst_real(   1, 186)=           4
      flst_real(   2, 186)=          -3
      flst_real(   3, 186)=          25
      flst_real(   4, 186)=         -11
      flst_real(   5, 186)=          12
      flst_real(   6, 186)=           0
      flst_real(   7, 186)=           0

      flst_real(   1, 187)=           4
      flst_real(   2, 187)=           3
      flst_real(   3, 187)=          25
      flst_real(   4, 187)=         -11
      flst_real(   5, 187)=          12
      flst_real(   6, 187)=           1
      flst_real(   7, 187)=           3

      flst_real(   1, 188)=           4
      flst_real(   2, 188)=           3
      flst_real(   3, 188)=          25
      flst_real(   4, 188)=         -11
      flst_real(   5, 188)=          12
      flst_real(   6, 188)=           3
      flst_real(   7, 188)=           3

      flst_real(   1, 189)=           4
      flst_real(   2, 189)=           3
      flst_real(   3, 189)=          25
      flst_real(   4, 189)=         -11
      flst_real(   5, 189)=          12
      flst_real(   6, 189)=           3
      flst_real(   7, 189)=           5

      flst_real(   1, 190)=           4
      flst_real(   2, 190)=          -5
      flst_real(   3, 190)=          25
      flst_real(   4, 190)=         -11
      flst_real(   5, 190)=          12
      flst_real(   6, 190)=           1
      flst_real(   7, 190)=          -1

      flst_real(   1, 191)=           4
      flst_real(   2, 191)=          -5
      flst_real(   3, 191)=          25
      flst_real(   4, 191)=         -11
      flst_real(   5, 191)=          12
      flst_real(   6, 191)=           1
      flst_real(   7, 191)=          -5

      flst_real(   1, 192)=           4
      flst_real(   2, 192)=          -5
      flst_real(   3, 192)=          25
      flst_real(   4, 192)=         -11
      flst_real(   5, 192)=          12
      flst_real(   6, 192)=           2
      flst_real(   7, 192)=          -2

      flst_real(   1, 193)=           4
      flst_real(   2, 193)=          -5
      flst_real(   3, 193)=          25
      flst_real(   4, 193)=         -11
      flst_real(   5, 193)=          12
      flst_real(   6, 193)=          -2
      flst_real(   7, 193)=           4

      flst_real(   1, 194)=           4
      flst_real(   2, 194)=          -5
      flst_real(   3, 194)=          25
      flst_real(   4, 194)=         -11
      flst_real(   5, 194)=          12
      flst_real(   6, 194)=           4
      flst_real(   7, 194)=          -4

      flst_real(   1, 195)=           4
      flst_real(   2, 195)=          -5
      flst_real(   3, 195)=          25
      flst_real(   4, 195)=         -11
      flst_real(   5, 195)=          12
      flst_real(   6, 195)=           3
      flst_real(   7, 195)=          -3

      flst_real(   1, 196)=           4
      flst_real(   2, 196)=          -5
      flst_real(   3, 196)=          25
      flst_real(   4, 196)=         -11
      flst_real(   5, 196)=          12
      flst_real(   6, 196)=           3
      flst_real(   7, 196)=          -5

      flst_real(   1, 197)=           4
      flst_real(   2, 197)=          -5
      flst_real(   3, 197)=          25
      flst_real(   4, 197)=         -11
      flst_real(   5, 197)=          12
      flst_real(   6, 197)=           5
      flst_real(   7, 197)=          -5

      flst_real(   1, 198)=           4
      flst_real(   2, 198)=          -5
      flst_real(   3, 198)=          25
      flst_real(   4, 198)=         -11
      flst_real(   5, 198)=          12
      flst_real(   6, 198)=           0
      flst_real(   7, 198)=           0

      flst_real(   1, 199)=           4
      flst_real(   2, 199)=           5
      flst_real(   3, 199)=          25
      flst_real(   4, 199)=         -11
      flst_real(   5, 199)=          12
      flst_real(   6, 199)=           1
      flst_real(   7, 199)=           5

      flst_real(   1, 200)=           4
      flst_real(   2, 200)=           5
      flst_real(   3, 200)=          25
      flst_real(   4, 200)=         -11
      flst_real(   5, 200)=          12
      flst_real(   6, 200)=           3
      flst_real(   7, 200)=           5

      flst_real(   1, 201)=           4
      flst_real(   2, 201)=           5
      flst_real(   3, 201)=          25
      flst_real(   4, 201)=         -11
      flst_real(   5, 201)=          12
      flst_real(   6, 201)=           5
      flst_real(   7, 201)=           5

      flst_real(   1, 202)=           4
      flst_real(   2, 202)=           0
      flst_real(   3, 202)=          25
      flst_real(   4, 202)=         -11
      flst_real(   5, 202)=          12
      flst_real(   6, 202)=           1
      flst_real(   7, 202)=           0

      flst_real(   1, 203)=           4
      flst_real(   2, 203)=           0
      flst_real(   3, 203)=          25
      flst_real(   4, 203)=         -11
      flst_real(   5, 203)=          12
      flst_real(   6, 203)=           3
      flst_real(   7, 203)=           0

      flst_real(   1, 204)=           4
      flst_real(   2, 204)=           0
      flst_real(   3, 204)=          25
      flst_real(   4, 204)=         -11
      flst_real(   5, 204)=          12
      flst_real(   6, 204)=           5
      flst_real(   7, 204)=           0

      flst_real(   1, 205)=          -3
      flst_real(   2, 205)=          -1
      flst_real(   3, 205)=          25
      flst_real(   4, 205)=         -11
      flst_real(   5, 205)=          12
      flst_real(   6, 205)=          -1
      flst_real(   7, 205)=          -2

      flst_real(   1, 206)=          -3
      flst_real(   2, 206)=          -1
      flst_real(   3, 206)=          25
      flst_real(   4, 206)=         -11
      flst_real(   5, 206)=          12
      flst_real(   6, 206)=          -1
      flst_real(   7, 206)=          -4

      flst_real(   1, 207)=          -3
      flst_real(   2, 207)=          -1
      flst_real(   3, 207)=          25
      flst_real(   4, 207)=         -11
      flst_real(   5, 207)=          12
      flst_real(   6, 207)=          -2
      flst_real(   7, 207)=          -3

      flst_real(   1, 208)=          -3
      flst_real(   2, 208)=          -1
      flst_real(   3, 208)=          25
      flst_real(   4, 208)=         -11
      flst_real(   5, 208)=          12
      flst_real(   6, 208)=          -4
      flst_real(   7, 208)=          -3

      flst_real(   1, 209)=          -3
      flst_real(   2, 209)=           1
      flst_real(   3, 209)=          25
      flst_real(   4, 209)=         -11
      flst_real(   5, 209)=          12
      flst_real(   6, 209)=           1
      flst_real(   7, 209)=          -2

      flst_real(   1, 210)=          -3
      flst_real(   2, 210)=           1
      flst_real(   3, 210)=          25
      flst_real(   4, 210)=         -11
      flst_real(   5, 210)=          12
      flst_real(   6, 210)=           1
      flst_real(   7, 210)=          -4

      flst_real(   1, 211)=          -3
      flst_real(   2, 211)=          -2
      flst_real(   3, 211)=          25
      flst_real(   4, 211)=         -11
      flst_real(   5, 211)=          12
      flst_real(   6, 211)=          -2
      flst_real(   7, 211)=          -2

      flst_real(   1, 212)=          -3
      flst_real(   2, 212)=          -2
      flst_real(   3, 212)=          25
      flst_real(   4, 212)=         -11
      flst_real(   5, 212)=          12
      flst_real(   6, 212)=          -2
      flst_real(   7, 212)=          -4

      flst_real(   1, 213)=          -3
      flst_real(   2, 213)=           2
      flst_real(   3, 213)=          25
      flst_real(   4, 213)=         -11
      flst_real(   5, 213)=          12
      flst_real(   6, 213)=           1
      flst_real(   7, 213)=          -1

      flst_real(   1, 214)=          -3
      flst_real(   2, 214)=           2
      flst_real(   3, 214)=          25
      flst_real(   4, 214)=         -11
      flst_real(   5, 214)=          12
      flst_real(   6, 214)=           1
      flst_real(   7, 214)=          -3

      flst_real(   1, 215)=          -3
      flst_real(   2, 215)=           2
      flst_real(   3, 215)=          25
      flst_real(   4, 215)=         -11
      flst_real(   5, 215)=          12
      flst_real(   6, 215)=           2
      flst_real(   7, 215)=          -2

      flst_real(   1, 216)=          -3
      flst_real(   2, 216)=           2
      flst_real(   3, 216)=          25
      flst_real(   4, 216)=         -11
      flst_real(   5, 216)=          12
      flst_real(   6, 216)=           2
      flst_real(   7, 216)=          -4

      flst_real(   1, 217)=          -3
      flst_real(   2, 217)=           2
      flst_real(   3, 217)=          25
      flst_real(   4, 217)=         -11
      flst_real(   5, 217)=          12
      flst_real(   6, 217)=           4
      flst_real(   7, 217)=          -4

      flst_real(   1, 218)=          -3
      flst_real(   2, 218)=           2
      flst_real(   3, 218)=          25
      flst_real(   4, 218)=         -11
      flst_real(   5, 218)=          12
      flst_real(   6, 218)=           3
      flst_real(   7, 218)=          -3

      flst_real(   1, 219)=          -3
      flst_real(   2, 219)=           2
      flst_real(   3, 219)=          25
      flst_real(   4, 219)=         -11
      flst_real(   5, 219)=          12
      flst_real(   6, 219)=          -3
      flst_real(   7, 219)=           5

      flst_real(   1, 220)=          -3
      flst_real(   2, 220)=           2
      flst_real(   3, 220)=          25
      flst_real(   4, 220)=         -11
      flst_real(   5, 220)=          12
      flst_real(   6, 220)=           5
      flst_real(   7, 220)=          -5

      flst_real(   1, 221)=          -3
      flst_real(   2, 221)=           2
      flst_real(   3, 221)=          25
      flst_real(   4, 221)=         -11
      flst_real(   5, 221)=          12
      flst_real(   6, 221)=           0
      flst_real(   7, 221)=           0

      flst_real(   1, 222)=          -3
      flst_real(   2, 222)=          -4
      flst_real(   3, 222)=          25
      flst_real(   4, 222)=         -11
      flst_real(   5, 222)=          12
      flst_real(   6, 222)=          -2
      flst_real(   7, 222)=          -4

      flst_real(   1, 223)=          -3
      flst_real(   2, 223)=          -4
      flst_real(   3, 223)=          25
      flst_real(   4, 223)=         -11
      flst_real(   5, 223)=          12
      flst_real(   6, 223)=          -4
      flst_real(   7, 223)=          -4

      flst_real(   1, 224)=          -3
      flst_real(   2, 224)=           4
      flst_real(   3, 224)=          25
      flst_real(   4, 224)=         -11
      flst_real(   5, 224)=          12
      flst_real(   6, 224)=           1
      flst_real(   7, 224)=          -1

      flst_real(   1, 225)=          -3
      flst_real(   2, 225)=           4
      flst_real(   3, 225)=          25
      flst_real(   4, 225)=         -11
      flst_real(   5, 225)=          12
      flst_real(   6, 225)=           1
      flst_real(   7, 225)=          -3

      flst_real(   1, 226)=          -3
      flst_real(   2, 226)=           4
      flst_real(   3, 226)=          25
      flst_real(   4, 226)=         -11
      flst_real(   5, 226)=          12
      flst_real(   6, 226)=           2
      flst_real(   7, 226)=          -2

      flst_real(   1, 227)=          -3
      flst_real(   2, 227)=           4
      flst_real(   3, 227)=          25
      flst_real(   4, 227)=         -11
      flst_real(   5, 227)=          12
      flst_real(   6, 227)=          -2
      flst_real(   7, 227)=           4

      flst_real(   1, 228)=          -3
      flst_real(   2, 228)=           4
      flst_real(   3, 228)=          25
      flst_real(   4, 228)=         -11
      flst_real(   5, 228)=          12
      flst_real(   6, 228)=           4
      flst_real(   7, 228)=          -4

      flst_real(   1, 229)=          -3
      flst_real(   2, 229)=           4
      flst_real(   3, 229)=          25
      flst_real(   4, 229)=         -11
      flst_real(   5, 229)=          12
      flst_real(   6, 229)=           3
      flst_real(   7, 229)=          -3

      flst_real(   1, 230)=          -3
      flst_real(   2, 230)=           4
      flst_real(   3, 230)=          25
      flst_real(   4, 230)=         -11
      flst_real(   5, 230)=          12
      flst_real(   6, 230)=          -3
      flst_real(   7, 230)=           5

      flst_real(   1, 231)=          -3
      flst_real(   2, 231)=           4
      flst_real(   3, 231)=          25
      flst_real(   4, 231)=         -11
      flst_real(   5, 231)=          12
      flst_real(   6, 231)=           5
      flst_real(   7, 231)=          -5

      flst_real(   1, 232)=          -3
      flst_real(   2, 232)=           4
      flst_real(   3, 232)=          25
      flst_real(   4, 232)=         -11
      flst_real(   5, 232)=          12
      flst_real(   6, 232)=           0
      flst_real(   7, 232)=           0

      flst_real(   1, 233)=          -3
      flst_real(   2, 233)=          -3
      flst_real(   3, 233)=          25
      flst_real(   4, 233)=         -11
      flst_real(   5, 233)=          12
      flst_real(   6, 233)=          -2
      flst_real(   7, 233)=          -3

      flst_real(   1, 234)=          -3
      flst_real(   2, 234)=          -3
      flst_real(   3, 234)=          25
      flst_real(   4, 234)=         -11
      flst_real(   5, 234)=          12
      flst_real(   6, 234)=          -4
      flst_real(   7, 234)=          -3

      flst_real(   1, 235)=          -3
      flst_real(   2, 235)=           3
      flst_real(   3, 235)=          25
      flst_real(   4, 235)=         -11
      flst_real(   5, 235)=          12
      flst_real(   6, 235)=           1
      flst_real(   7, 235)=          -2

      flst_real(   1, 236)=          -3
      flst_real(   2, 236)=           3
      flst_real(   3, 236)=          25
      flst_real(   4, 236)=         -11
      flst_real(   5, 236)=          12
      flst_real(   6, 236)=           1
      flst_real(   7, 236)=          -4

      flst_real(   1, 237)=          -3
      flst_real(   2, 237)=           3
      flst_real(   3, 237)=          25
      flst_real(   4, 237)=         -11
      flst_real(   5, 237)=          12
      flst_real(   6, 237)=          -2
      flst_real(   7, 237)=           3

      flst_real(   1, 238)=          -3
      flst_real(   2, 238)=           3
      flst_real(   3, 238)=          25
      flst_real(   4, 238)=         -11
      flst_real(   5, 238)=          12
      flst_real(   6, 238)=          -2
      flst_real(   7, 238)=           5

      flst_real(   1, 239)=          -3
      flst_real(   2, 239)=           3
      flst_real(   3, 239)=          25
      flst_real(   4, 239)=         -11
      flst_real(   5, 239)=          12
      flst_real(   6, 239)=          -4
      flst_real(   7, 239)=           3

      flst_real(   1, 240)=          -3
      flst_real(   2, 240)=           3
      flst_real(   3, 240)=          25
      flst_real(   4, 240)=         -11
      flst_real(   5, 240)=          12
      flst_real(   6, 240)=          -4
      flst_real(   7, 240)=           5

      flst_real(   1, 241)=          -3
      flst_real(   2, 241)=          -5
      flst_real(   3, 241)=          25
      flst_real(   4, 241)=         -11
      flst_real(   5, 241)=          12
      flst_real(   6, 241)=          -2
      flst_real(   7, 241)=          -3

      flst_real(   1, 242)=          -3
      flst_real(   2, 242)=          -5
      flst_real(   3, 242)=          25
      flst_real(   4, 242)=         -11
      flst_real(   5, 242)=          12
      flst_real(   6, 242)=          -2
      flst_real(   7, 242)=          -5

      flst_real(   1, 243)=          -3
      flst_real(   2, 243)=          -5
      flst_real(   3, 243)=          25
      flst_real(   4, 243)=         -11
      flst_real(   5, 243)=          12
      flst_real(   6, 243)=          -4
      flst_real(   7, 243)=          -3

      flst_real(   1, 244)=          -3
      flst_real(   2, 244)=          -5
      flst_real(   3, 244)=          25
      flst_real(   4, 244)=         -11
      flst_real(   5, 244)=          12
      flst_real(   6, 244)=          -4
      flst_real(   7, 244)=          -5

      flst_real(   1, 245)=          -3
      flst_real(   2, 245)=           5
      flst_real(   3, 245)=          25
      flst_real(   4, 245)=         -11
      flst_real(   5, 245)=          12
      flst_real(   6, 245)=          -2
      flst_real(   7, 245)=           5

      flst_real(   1, 246)=          -3
      flst_real(   2, 246)=           5
      flst_real(   3, 246)=          25
      flst_real(   4, 246)=         -11
      flst_real(   5, 246)=          12
      flst_real(   6, 246)=          -4
      flst_real(   7, 246)=           5

      flst_real(   1, 247)=          -3
      flst_real(   2, 247)=           0
      flst_real(   3, 247)=          25
      flst_real(   4, 247)=         -11
      flst_real(   5, 247)=          12
      flst_real(   6, 247)=          -2
      flst_real(   7, 247)=           0

      flst_real(   1, 248)=          -3
      flst_real(   2, 248)=           0
      flst_real(   3, 248)=          25
      flst_real(   4, 248)=         -11
      flst_real(   5, 248)=          12
      flst_real(   6, 248)=          -4
      flst_real(   7, 248)=           0

      flst_real(   1, 249)=           3
      flst_real(   2, 249)=          -1
      flst_real(   3, 249)=          25
      flst_real(   4, 249)=         -11
      flst_real(   5, 249)=          12
      flst_real(   6, 249)=          -2
      flst_real(   7, 249)=           3

      flst_real(   1, 250)=           3
      flst_real(   2, 250)=          -1
      flst_real(   3, 250)=          25
      flst_real(   4, 250)=         -11
      flst_real(   5, 250)=          12
      flst_real(   6, 250)=          -4
      flst_real(   7, 250)=           3

      flst_real(   1, 251)=           3
      flst_real(   2, 251)=           2
      flst_real(   3, 251)=          25
      flst_real(   4, 251)=         -11
      flst_real(   5, 251)=          12
      flst_real(   6, 251)=           1
      flst_real(   7, 251)=           3

      flst_real(   1, 252)=           3
      flst_real(   2, 252)=           2
      flst_real(   3, 252)=          25
      flst_real(   4, 252)=         -11
      flst_real(   5, 252)=          12
      flst_real(   6, 252)=           3
      flst_real(   7, 252)=           3

      flst_real(   1, 253)=           3
      flst_real(   2, 253)=           2
      flst_real(   3, 253)=          25
      flst_real(   4, 253)=         -11
      flst_real(   5, 253)=          12
      flst_real(   6, 253)=           3
      flst_real(   7, 253)=           5

      flst_real(   1, 254)=           3
      flst_real(   2, 254)=           4
      flst_real(   3, 254)=          25
      flst_real(   4, 254)=         -11
      flst_real(   5, 254)=          12
      flst_real(   6, 254)=           1
      flst_real(   7, 254)=           3

      flst_real(   1, 255)=           3
      flst_real(   2, 255)=           4
      flst_real(   3, 255)=          25
      flst_real(   4, 255)=         -11
      flst_real(   5, 255)=          12
      flst_real(   6, 255)=           3
      flst_real(   7, 255)=           3

      flst_real(   1, 256)=           3
      flst_real(   2, 256)=           4
      flst_real(   3, 256)=          25
      flst_real(   4, 256)=         -11
      flst_real(   5, 256)=          12
      flst_real(   6, 256)=           3
      flst_real(   7, 256)=           5

      flst_real(   1, 257)=           3
      flst_real(   2, 257)=          -3
      flst_real(   3, 257)=          25
      flst_real(   4, 257)=         -11
      flst_real(   5, 257)=          12
      flst_real(   6, 257)=           1
      flst_real(   7, 257)=          -2

      flst_real(   1, 258)=           3
      flst_real(   2, 258)=          -3
      flst_real(   3, 258)=          25
      flst_real(   4, 258)=         -11
      flst_real(   5, 258)=          12
      flst_real(   6, 258)=           1
      flst_real(   7, 258)=          -4

      flst_real(   1, 259)=           3
      flst_real(   2, 259)=          -3
      flst_real(   3, 259)=          25
      flst_real(   4, 259)=         -11
      flst_real(   5, 259)=          12
      flst_real(   6, 259)=          -2
      flst_real(   7, 259)=           3

      flst_real(   1, 260)=           3
      flst_real(   2, 260)=          -3
      flst_real(   3, 260)=          25
      flst_real(   4, 260)=         -11
      flst_real(   5, 260)=          12
      flst_real(   6, 260)=          -2
      flst_real(   7, 260)=           5

      flst_real(   1, 261)=           3
      flst_real(   2, 261)=          -3
      flst_real(   3, 261)=          25
      flst_real(   4, 261)=         -11
      flst_real(   5, 261)=          12
      flst_real(   6, 261)=          -4
      flst_real(   7, 261)=           3

      flst_real(   1, 262)=           3
      flst_real(   2, 262)=          -3
      flst_real(   3, 262)=          25
      flst_real(   4, 262)=         -11
      flst_real(   5, 262)=          12
      flst_real(   6, 262)=          -4
      flst_real(   7, 262)=           5

      flst_real(   1, 263)=           3
      flst_real(   2, 263)=          -5
      flst_real(   3, 263)=          25
      flst_real(   4, 263)=         -11
      flst_real(   5, 263)=          12
      flst_real(   6, 263)=          -2
      flst_real(   7, 263)=           3

      flst_real(   1, 264)=           3
      flst_real(   2, 264)=          -5
      flst_real(   3, 264)=          25
      flst_real(   4, 264)=         -11
      flst_real(   5, 264)=          12
      flst_real(   6, 264)=          -4
      flst_real(   7, 264)=           3

      flst_real(   1, 265)=          -5
      flst_real(   2, 265)=          -1
      flst_real(   3, 265)=          25
      flst_real(   4, 265)=         -11
      flst_real(   5, 265)=          12
      flst_real(   6, 265)=          -1
      flst_real(   7, 265)=          -2

      flst_real(   1, 266)=          -5
      flst_real(   2, 266)=          -1
      flst_real(   3, 266)=          25
      flst_real(   4, 266)=         -11
      flst_real(   5, 266)=          12
      flst_real(   6, 266)=          -1
      flst_real(   7, 266)=          -4

      flst_real(   1, 267)=          -5
      flst_real(   2, 267)=          -1
      flst_real(   3, 267)=          25
      flst_real(   4, 267)=         -11
      flst_real(   5, 267)=          12
      flst_real(   6, 267)=          -2
      flst_real(   7, 267)=          -5

      flst_real(   1, 268)=          -5
      flst_real(   2, 268)=          -1
      flst_real(   3, 268)=          25
      flst_real(   4, 268)=         -11
      flst_real(   5, 268)=          12
      flst_real(   6, 268)=          -4
      flst_real(   7, 268)=          -5

      flst_real(   1, 269)=          -5
      flst_real(   2, 269)=           1
      flst_real(   3, 269)=          25
      flst_real(   4, 269)=         -11
      flst_real(   5, 269)=          12
      flst_real(   6, 269)=           1
      flst_real(   7, 269)=          -2

      flst_real(   1, 270)=          -5
      flst_real(   2, 270)=           1
      flst_real(   3, 270)=          25
      flst_real(   4, 270)=         -11
      flst_real(   5, 270)=          12
      flst_real(   6, 270)=           1
      flst_real(   7, 270)=          -4

      flst_real(   1, 271)=          -5
      flst_real(   2, 271)=          -2
      flst_real(   3, 271)=          25
      flst_real(   4, 271)=         -11
      flst_real(   5, 271)=          12
      flst_real(   6, 271)=          -2
      flst_real(   7, 271)=          -2

      flst_real(   1, 272)=          -5
      flst_real(   2, 272)=          -2
      flst_real(   3, 272)=          25
      flst_real(   4, 272)=         -11
      flst_real(   5, 272)=          12
      flst_real(   6, 272)=          -2
      flst_real(   7, 272)=          -4

      flst_real(   1, 273)=          -5
      flst_real(   2, 273)=           2
      flst_real(   3, 273)=          25
      flst_real(   4, 273)=         -11
      flst_real(   5, 273)=          12
      flst_real(   6, 273)=           1
      flst_real(   7, 273)=          -1

      flst_real(   1, 274)=          -5
      flst_real(   2, 274)=           2
      flst_real(   3, 274)=          25
      flst_real(   4, 274)=         -11
      flst_real(   5, 274)=          12
      flst_real(   6, 274)=           1
      flst_real(   7, 274)=          -5

      flst_real(   1, 275)=          -5
      flst_real(   2, 275)=           2
      flst_real(   3, 275)=          25
      flst_real(   4, 275)=         -11
      flst_real(   5, 275)=          12
      flst_real(   6, 275)=           2
      flst_real(   7, 275)=          -2

      flst_real(   1, 276)=          -5
      flst_real(   2, 276)=           2
      flst_real(   3, 276)=          25
      flst_real(   4, 276)=         -11
      flst_real(   5, 276)=          12
      flst_real(   6, 276)=           2
      flst_real(   7, 276)=          -4

      flst_real(   1, 277)=          -5
      flst_real(   2, 277)=           2
      flst_real(   3, 277)=          25
      flst_real(   4, 277)=         -11
      flst_real(   5, 277)=          12
      flst_real(   6, 277)=           4
      flst_real(   7, 277)=          -4

      flst_real(   1, 278)=          -5
      flst_real(   2, 278)=           2
      flst_real(   3, 278)=          25
      flst_real(   4, 278)=         -11
      flst_real(   5, 278)=          12
      flst_real(   6, 278)=           3
      flst_real(   7, 278)=          -3

      flst_real(   1, 279)=          -5
      flst_real(   2, 279)=           2
      flst_real(   3, 279)=          25
      flst_real(   4, 279)=         -11
      flst_real(   5, 279)=          12
      flst_real(   6, 279)=           3
      flst_real(   7, 279)=          -5

      flst_real(   1, 280)=          -5
      flst_real(   2, 280)=           2
      flst_real(   3, 280)=          25
      flst_real(   4, 280)=         -11
      flst_real(   5, 280)=          12
      flst_real(   6, 280)=           5
      flst_real(   7, 280)=          -5

      flst_real(   1, 281)=          -5
      flst_real(   2, 281)=           2
      flst_real(   3, 281)=          25
      flst_real(   4, 281)=         -11
      flst_real(   5, 281)=          12
      flst_real(   6, 281)=           0
      flst_real(   7, 281)=           0

      flst_real(   1, 282)=          -5
      flst_real(   2, 282)=          -4
      flst_real(   3, 282)=          25
      flst_real(   4, 282)=         -11
      flst_real(   5, 282)=          12
      flst_real(   6, 282)=          -2
      flst_real(   7, 282)=          -4

      flst_real(   1, 283)=          -5
      flst_real(   2, 283)=          -4
      flst_real(   3, 283)=          25
      flst_real(   4, 283)=         -11
      flst_real(   5, 283)=          12
      flst_real(   6, 283)=          -4
      flst_real(   7, 283)=          -4

      flst_real(   1, 284)=          -5
      flst_real(   2, 284)=           4
      flst_real(   3, 284)=          25
      flst_real(   4, 284)=         -11
      flst_real(   5, 284)=          12
      flst_real(   6, 284)=           1
      flst_real(   7, 284)=          -1

      flst_real(   1, 285)=          -5
      flst_real(   2, 285)=           4
      flst_real(   3, 285)=          25
      flst_real(   4, 285)=         -11
      flst_real(   5, 285)=          12
      flst_real(   6, 285)=           1
      flst_real(   7, 285)=          -5

      flst_real(   1, 286)=          -5
      flst_real(   2, 286)=           4
      flst_real(   3, 286)=          25
      flst_real(   4, 286)=         -11
      flst_real(   5, 286)=          12
      flst_real(   6, 286)=           2
      flst_real(   7, 286)=          -2

      flst_real(   1, 287)=          -5
      flst_real(   2, 287)=           4
      flst_real(   3, 287)=          25
      flst_real(   4, 287)=         -11
      flst_real(   5, 287)=          12
      flst_real(   6, 287)=          -2
      flst_real(   7, 287)=           4

      flst_real(   1, 288)=          -5
      flst_real(   2, 288)=           4
      flst_real(   3, 288)=          25
      flst_real(   4, 288)=         -11
      flst_real(   5, 288)=          12
      flst_real(   6, 288)=           4
      flst_real(   7, 288)=          -4

      flst_real(   1, 289)=          -5
      flst_real(   2, 289)=           4
      flst_real(   3, 289)=          25
      flst_real(   4, 289)=         -11
      flst_real(   5, 289)=          12
      flst_real(   6, 289)=           3
      flst_real(   7, 289)=          -3

      flst_real(   1, 290)=          -5
      flst_real(   2, 290)=           4
      flst_real(   3, 290)=          25
      flst_real(   4, 290)=         -11
      flst_real(   5, 290)=          12
      flst_real(   6, 290)=           3
      flst_real(   7, 290)=          -5

      flst_real(   1, 291)=          -5
      flst_real(   2, 291)=           4
      flst_real(   3, 291)=          25
      flst_real(   4, 291)=         -11
      flst_real(   5, 291)=          12
      flst_real(   6, 291)=           5
      flst_real(   7, 291)=          -5

      flst_real(   1, 292)=          -5
      flst_real(   2, 292)=           4
      flst_real(   3, 292)=          25
      flst_real(   4, 292)=         -11
      flst_real(   5, 292)=          12
      flst_real(   6, 292)=           0
      flst_real(   7, 292)=           0

      flst_real(   1, 293)=          -5
      flst_real(   2, 293)=          -3
      flst_real(   3, 293)=          25
      flst_real(   4, 293)=         -11
      flst_real(   5, 293)=          12
      flst_real(   6, 293)=          -2
      flst_real(   7, 293)=          -3

      flst_real(   1, 294)=          -5
      flst_real(   2, 294)=          -3
      flst_real(   3, 294)=          25
      flst_real(   4, 294)=         -11
      flst_real(   5, 294)=          12
      flst_real(   6, 294)=          -2
      flst_real(   7, 294)=          -5

      flst_real(   1, 295)=          -5
      flst_real(   2, 295)=          -3
      flst_real(   3, 295)=          25
      flst_real(   4, 295)=         -11
      flst_real(   5, 295)=          12
      flst_real(   6, 295)=          -4
      flst_real(   7, 295)=          -3

      flst_real(   1, 296)=          -5
      flst_real(   2, 296)=          -3
      flst_real(   3, 296)=          25
      flst_real(   4, 296)=         -11
      flst_real(   5, 296)=          12
      flst_real(   6, 296)=          -4
      flst_real(   7, 296)=          -5

      flst_real(   1, 297)=          -5
      flst_real(   2, 297)=           3
      flst_real(   3, 297)=          25
      flst_real(   4, 297)=         -11
      flst_real(   5, 297)=          12
      flst_real(   6, 297)=          -2
      flst_real(   7, 297)=           3

      flst_real(   1, 298)=          -5
      flst_real(   2, 298)=           3
      flst_real(   3, 298)=          25
      flst_real(   4, 298)=         -11
      flst_real(   5, 298)=          12
      flst_real(   6, 298)=          -4
      flst_real(   7, 298)=           3

      flst_real(   1, 299)=          -5
      flst_real(   2, 299)=          -5
      flst_real(   3, 299)=          25
      flst_real(   4, 299)=         -11
      flst_real(   5, 299)=          12
      flst_real(   6, 299)=          -2
      flst_real(   7, 299)=          -5

      flst_real(   1, 300)=          -5
      flst_real(   2, 300)=          -5
      flst_real(   3, 300)=          25
      flst_real(   4, 300)=         -11
      flst_real(   5, 300)=          12
      flst_real(   6, 300)=          -4
      flst_real(   7, 300)=          -5

      flst_real(   1, 301)=          -5
      flst_real(   2, 301)=           5
      flst_real(   3, 301)=          25
      flst_real(   4, 301)=         -11
      flst_real(   5, 301)=          12
      flst_real(   6, 301)=           1
      flst_real(   7, 301)=          -2

      flst_real(   1, 302)=          -5
      flst_real(   2, 302)=           5
      flst_real(   3, 302)=          25
      flst_real(   4, 302)=         -11
      flst_real(   5, 302)=          12
      flst_real(   6, 302)=           1
      flst_real(   7, 302)=          -4

      flst_real(   1, 303)=          -5
      flst_real(   2, 303)=           5
      flst_real(   3, 303)=          25
      flst_real(   4, 303)=         -11
      flst_real(   5, 303)=          12
      flst_real(   6, 303)=          -2
      flst_real(   7, 303)=           3

      flst_real(   1, 304)=          -5
      flst_real(   2, 304)=           5
      flst_real(   3, 304)=          25
      flst_real(   4, 304)=         -11
      flst_real(   5, 304)=          12
      flst_real(   6, 304)=          -2
      flst_real(   7, 304)=           5

      flst_real(   1, 305)=          -5
      flst_real(   2, 305)=           5
      flst_real(   3, 305)=          25
      flst_real(   4, 305)=         -11
      flst_real(   5, 305)=          12
      flst_real(   6, 305)=          -4
      flst_real(   7, 305)=           3

      flst_real(   1, 306)=          -5
      flst_real(   2, 306)=           5
      flst_real(   3, 306)=          25
      flst_real(   4, 306)=         -11
      flst_real(   5, 306)=          12
      flst_real(   6, 306)=          -4
      flst_real(   7, 306)=           5

      flst_real(   1, 307)=          -5
      flst_real(   2, 307)=           0
      flst_real(   3, 307)=          25
      flst_real(   4, 307)=         -11
      flst_real(   5, 307)=          12
      flst_real(   6, 307)=          -2
      flst_real(   7, 307)=           0

      flst_real(   1, 308)=          -5
      flst_real(   2, 308)=           0
      flst_real(   3, 308)=          25
      flst_real(   4, 308)=         -11
      flst_real(   5, 308)=          12
      flst_real(   6, 308)=          -4
      flst_real(   7, 308)=           0

      flst_real(   1, 309)=           5
      flst_real(   2, 309)=          -1
      flst_real(   3, 309)=          25
      flst_real(   4, 309)=         -11
      flst_real(   5, 309)=          12
      flst_real(   6, 309)=          -2
      flst_real(   7, 309)=           5

      flst_real(   1, 310)=           5
      flst_real(   2, 310)=          -1
      flst_real(   3, 310)=          25
      flst_real(   4, 310)=         -11
      flst_real(   5, 310)=          12
      flst_real(   6, 310)=          -4
      flst_real(   7, 310)=           5

      flst_real(   1, 311)=           5
      flst_real(   2, 311)=           2
      flst_real(   3, 311)=          25
      flst_real(   4, 311)=         -11
      flst_real(   5, 311)=          12
      flst_real(   6, 311)=           1
      flst_real(   7, 311)=           5

      flst_real(   1, 312)=           5
      flst_real(   2, 312)=           2
      flst_real(   3, 312)=          25
      flst_real(   4, 312)=         -11
      flst_real(   5, 312)=          12
      flst_real(   6, 312)=           3
      flst_real(   7, 312)=           5

      flst_real(   1, 313)=           5
      flst_real(   2, 313)=           2
      flst_real(   3, 313)=          25
      flst_real(   4, 313)=         -11
      flst_real(   5, 313)=          12
      flst_real(   6, 313)=           5
      flst_real(   7, 313)=           5

      flst_real(   1, 314)=           5
      flst_real(   2, 314)=           4
      flst_real(   3, 314)=          25
      flst_real(   4, 314)=         -11
      flst_real(   5, 314)=          12
      flst_real(   6, 314)=           1
      flst_real(   7, 314)=           5

      flst_real(   1, 315)=           5
      flst_real(   2, 315)=           4
      flst_real(   3, 315)=          25
      flst_real(   4, 315)=         -11
      flst_real(   5, 315)=          12
      flst_real(   6, 315)=           3
      flst_real(   7, 315)=           5

      flst_real(   1, 316)=           5
      flst_real(   2, 316)=           4
      flst_real(   3, 316)=          25
      flst_real(   4, 316)=         -11
      flst_real(   5, 316)=          12
      flst_real(   6, 316)=           5
      flst_real(   7, 316)=           5

      flst_real(   1, 317)=           5
      flst_real(   2, 317)=          -3
      flst_real(   3, 317)=          25
      flst_real(   4, 317)=         -11
      flst_real(   5, 317)=          12
      flst_real(   6, 317)=          -2
      flst_real(   7, 317)=           5

      flst_real(   1, 318)=           5
      flst_real(   2, 318)=          -3
      flst_real(   3, 318)=          25
      flst_real(   4, 318)=         -11
      flst_real(   5, 318)=          12
      flst_real(   6, 318)=          -4
      flst_real(   7, 318)=           5

      flst_real(   1, 319)=           5
      flst_real(   2, 319)=          -5
      flst_real(   3, 319)=          25
      flst_real(   4, 319)=         -11
      flst_real(   5, 319)=          12
      flst_real(   6, 319)=           1
      flst_real(   7, 319)=          -2

      flst_real(   1, 320)=           5
      flst_real(   2, 320)=          -5
      flst_real(   3, 320)=          25
      flst_real(   4, 320)=         -11
      flst_real(   5, 320)=          12
      flst_real(   6, 320)=           1
      flst_real(   7, 320)=          -4

      flst_real(   1, 321)=           5
      flst_real(   2, 321)=          -5
      flst_real(   3, 321)=          25
      flst_real(   4, 321)=         -11
      flst_real(   5, 321)=          12
      flst_real(   6, 321)=          -2
      flst_real(   7, 321)=           3

      flst_real(   1, 322)=           5
      flst_real(   2, 322)=          -5
      flst_real(   3, 322)=          25
      flst_real(   4, 322)=         -11
      flst_real(   5, 322)=          12
      flst_real(   6, 322)=          -2
      flst_real(   7, 322)=           5

      flst_real(   1, 323)=           5
      flst_real(   2, 323)=          -5
      flst_real(   3, 323)=          25
      flst_real(   4, 323)=         -11
      flst_real(   5, 323)=          12
      flst_real(   6, 323)=          -4
      flst_real(   7, 323)=           3

      flst_real(   1, 324)=           5
      flst_real(   2, 324)=          -5
      flst_real(   3, 324)=          25
      flst_real(   4, 324)=         -11
      flst_real(   5, 324)=          12
      flst_real(   6, 324)=          -4
      flst_real(   7, 324)=           5

      flst_real(   1, 325)=           0
      flst_real(   2, 325)=          -1
      flst_real(   3, 325)=          25
      flst_real(   4, 325)=         -11
      flst_real(   5, 325)=          12
      flst_real(   6, 325)=          -2
      flst_real(   7, 325)=           0

      flst_real(   1, 326)=           0
      flst_real(   2, 326)=          -1
      flst_real(   3, 326)=          25
      flst_real(   4, 326)=         -11
      flst_real(   5, 326)=          12
      flst_real(   6, 326)=          -4
      flst_real(   7, 326)=           0

      flst_real(   1, 327)=           0
      flst_real(   2, 327)=           2
      flst_real(   3, 327)=          25
      flst_real(   4, 327)=         -11
      flst_real(   5, 327)=          12
      flst_real(   6, 327)=           1
      flst_real(   7, 327)=           0

      flst_real(   1, 328)=           0
      flst_real(   2, 328)=           2
      flst_real(   3, 328)=          25
      flst_real(   4, 328)=         -11
      flst_real(   5, 328)=          12
      flst_real(   6, 328)=           3
      flst_real(   7, 328)=           0

      flst_real(   1, 329)=           0
      flst_real(   2, 329)=           2
      flst_real(   3, 329)=          25
      flst_real(   4, 329)=         -11
      flst_real(   5, 329)=          12
      flst_real(   6, 329)=           5
      flst_real(   7, 329)=           0

      flst_real(   1, 330)=           0
      flst_real(   2, 330)=           4
      flst_real(   3, 330)=          25
      flst_real(   4, 330)=         -11
      flst_real(   5, 330)=          12
      flst_real(   6, 330)=           1
      flst_real(   7, 330)=           0

      flst_real(   1, 331)=           0
      flst_real(   2, 331)=           4
      flst_real(   3, 331)=          25
      flst_real(   4, 331)=         -11
      flst_real(   5, 331)=          12
      flst_real(   6, 331)=           3
      flst_real(   7, 331)=           0

      flst_real(   1, 332)=           0
      flst_real(   2, 332)=           4
      flst_real(   3, 332)=          25
      flst_real(   4, 332)=         -11
      flst_real(   5, 332)=          12
      flst_real(   6, 332)=           5
      flst_real(   7, 332)=           0

      flst_real(   1, 333)=           0
      flst_real(   2, 333)=          -3
      flst_real(   3, 333)=          25
      flst_real(   4, 333)=         -11
      flst_real(   5, 333)=          12
      flst_real(   6, 333)=          -2
      flst_real(   7, 333)=           0

      flst_real(   1, 334)=           0
      flst_real(   2, 334)=          -3
      flst_real(   3, 334)=          25
      flst_real(   4, 334)=         -11
      flst_real(   5, 334)=          12
      flst_real(   6, 334)=          -4
      flst_real(   7, 334)=           0

      flst_real(   1, 335)=           0
      flst_real(   2, 335)=          -5
      flst_real(   3, 335)=          25
      flst_real(   4, 335)=         -11
      flst_real(   5, 335)=          12
      flst_real(   6, 335)=          -2
      flst_real(   7, 335)=           0

      flst_real(   1, 336)=           0
      flst_real(   2, 336)=          -5
      flst_real(   3, 336)=          25
      flst_real(   4, 336)=         -11
      flst_real(   5, 336)=          12
      flst_real(   6, 336)=          -4
      flst_real(   7, 336)=           0

      flst_real(   1, 337)=           0
      flst_real(   2, 337)=           0
      flst_real(   3, 337)=          25
      flst_real(   4, 337)=         -11
      flst_real(   5, 337)=          12
      flst_real(   6, 337)=           1
      flst_real(   7, 337)=          -2

      flst_real(   1, 338)=           0
      flst_real(   2, 338)=           0
      flst_real(   3, 338)=          25
      flst_real(   4, 338)=         -11
      flst_real(   5, 338)=          12
      flst_real(   6, 338)=           1
      flst_real(   7, 338)=          -4

      flst_real(   1, 339)=           0
      flst_real(   2, 339)=           0
      flst_real(   3, 339)=          25
      flst_real(   4, 339)=         -11
      flst_real(   5, 339)=          12
      flst_real(   6, 339)=          -2
      flst_real(   7, 339)=           3

      flst_real(   1, 340)=           0
      flst_real(   2, 340)=           0
      flst_real(   3, 340)=          25
      flst_real(   4, 340)=         -11
      flst_real(   5, 340)=          12
      flst_real(   6, 340)=          -2
      flst_real(   7, 340)=           5

      flst_real(   1, 341)=           0
      flst_real(   2, 341)=           0
      flst_real(   3, 341)=          25
      flst_real(   4, 341)=         -11
      flst_real(   5, 341)=          12
      flst_real(   6, 341)=          -4
      flst_real(   7, 341)=           3

      flst_real(   1, 342)=           0
      flst_real(   2, 342)=           0
      flst_real(   3, 342)=          25
      flst_real(   4, 342)=         -11
      flst_real(   5, 342)=          12
      flst_real(   6, 342)=          -4
      flst_real(   7, 342)=           5

      flst_nreal=         342

      return
      end

