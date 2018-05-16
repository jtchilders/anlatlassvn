      subroutine init_processes
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
      include "pwhg_st.h"
      include "coupl.inc"
      integer i
 
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
 
      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           1
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=         -11
      flst_born(   5,   1)=          11
 
      flst_born(   1,   2)=           1
      flst_born(   2,   2)=          -1
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=         -11
      flst_born(   5,   2)=          11
 
      flst_born(   1,   3)=          -2
      flst_born(   2,   3)=           2
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=         -11
      flst_born(   5,   3)=          11
 
      flst_born(   1,   4)=           2
      flst_born(   2,   4)=          -2
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=         -11
      flst_born(   5,   4)=          11
 
      flst_born(   1,   5)=          -4
      flst_born(   2,   5)=           4
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=         -11
      flst_born(   5,   5)=          11
 
      flst_born(   1,   6)=           4
      flst_born(   2,   6)=          -4
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=         -11
      flst_born(   5,   6)=          11
 
      flst_born(   1,   7)=          -3
      flst_born(   2,   7)=           3
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=         -11
      flst_born(   5,   7)=          11
 
      flst_born(   1,   8)=           3
      flst_born(   2,   8)=          -3
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=         -11
      flst_born(   5,   8)=          11
 
      flst_born(   1,   9)=          -5
      flst_born(   2,   9)=           5
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=         -11
      flst_born(   5,   9)=          11
 
      flst_born(   1,  10)=           5
      flst_born(   2,  10)=          -5
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=         -11
      flst_born(   5,  10)=          11
 
      flst_nborn=          10
 
      return
      end
 
 
 
      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=           1
      flst_real(   3,   1)=          25
      flst_real(   4,   1)=         -11
      flst_real(   5,   1)=          11
      flst_real(   6,   1)=           0
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           0
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=         -11
      flst_real(   5,   2)=          11
      flst_real(   6,   2)=          -1
 
      flst_real(   1,   3)=           1
      flst_real(   2,   3)=          -1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=         -11
      flst_real(   5,   3)=          11
      flst_real(   6,   3)=           0
 
      flst_real(   1,   4)=           1
      flst_real(   2,   4)=           0
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=         -11
      flst_real(   5,   4)=          11
      flst_real(   6,   4)=           1
 
      flst_real(   1,   5)=          -2
      flst_real(   2,   5)=           2
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=         -11
      flst_real(   5,   5)=          11
      flst_real(   6,   5)=           0
 
      flst_real(   1,   6)=          -2
      flst_real(   2,   6)=           0
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=         -11
      flst_real(   5,   6)=          11
      flst_real(   6,   6)=          -2
 
      flst_real(   1,   7)=           2
      flst_real(   2,   7)=          -2
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=         -11
      flst_real(   5,   7)=          11
      flst_real(   6,   7)=           0
 
      flst_real(   1,   8)=           2
      flst_real(   2,   8)=           0
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=         -11
      flst_real(   5,   8)=          11
      flst_real(   6,   8)=           2
 
      flst_real(   1,   9)=          -4
      flst_real(   2,   9)=           4
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=         -11
      flst_real(   5,   9)=          11
      flst_real(   6,   9)=           0
 
      flst_real(   1,  10)=          -4
      flst_real(   2,  10)=           0
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=         -11
      flst_real(   5,  10)=          11
      flst_real(   6,  10)=          -4
 
      flst_real(   1,  11)=           4
      flst_real(   2,  11)=          -4
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=         -11
      flst_real(   5,  11)=          11
      flst_real(   6,  11)=           0
 
      flst_real(   1,  12)=           4
      flst_real(   2,  12)=           0
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=         -11
      flst_real(   5,  12)=          11
      flst_real(   6,  12)=           4
 
      flst_real(   1,  13)=          -3
      flst_real(   2,  13)=           3
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=         -11
      flst_real(   5,  13)=          11
      flst_real(   6,  13)=           0
 
      flst_real(   1,  14)=          -3
      flst_real(   2,  14)=           0
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=         -11
      flst_real(   5,  14)=          11
      flst_real(   6,  14)=          -3
 
      flst_real(   1,  15)=           3
      flst_real(   2,  15)=          -3
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=         -11
      flst_real(   5,  15)=          11
      flst_real(   6,  15)=           0
 
      flst_real(   1,  16)=           3
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=         -11
      flst_real(   5,  16)=          11
      flst_real(   6,  16)=           3
 
      flst_real(   1,  17)=          -5
      flst_real(   2,  17)=           5
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=         -11
      flst_real(   5,  17)=          11
      flst_real(   6,  17)=           0
 
      flst_real(   1,  18)=          -5
      flst_real(   2,  18)=           0
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=         -11
      flst_real(   5,  18)=          11
      flst_real(   6,  18)=          -5
 
      flst_real(   1,  19)=           5
      flst_real(   2,  19)=          -5
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=         -11
      flst_real(   5,  19)=          11
      flst_real(   6,  19)=           0
 
      flst_real(   1,  20)=           5
      flst_real(   2,  20)=           0
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=         -11
      flst_real(   5,  20)=          11
      flst_real(   6,  20)=           5
 
      flst_real(   1,  21)=           0
      flst_real(   2,  21)=          -1
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=         -11
      flst_real(   5,  21)=          11
      flst_real(   6,  21)=          -1
 
      flst_real(   1,  22)=           0
      flst_real(   2,  22)=           1
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=         -11
      flst_real(   5,  22)=          11
      flst_real(   6,  22)=           1
 
      flst_real(   1,  23)=           0
      flst_real(   2,  23)=          -2
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=         -11
      flst_real(   5,  23)=          11
      flst_real(   6,  23)=          -2
 
      flst_real(   1,  24)=           0
      flst_real(   2,  24)=           2
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=         -11
      flst_real(   5,  24)=          11
      flst_real(   6,  24)=           2
 
      flst_real(   1,  25)=           0
      flst_real(   2,  25)=          -4
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=         -11
      flst_real(   5,  25)=          11
      flst_real(   6,  25)=          -4
 
      flst_real(   1,  26)=           0
      flst_real(   2,  26)=           4
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=         -11
      flst_real(   5,  26)=          11
      flst_real(   6,  26)=           4
 
      flst_real(   1,  27)=           0
      flst_real(   2,  27)=          -3
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=         -11
      flst_real(   5,  27)=          11
      flst_real(   6,  27)=          -3
 
      flst_real(   1,  28)=           0
      flst_real(   2,  28)=           3
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=         -11
      flst_real(   5,  28)=          11
      flst_real(   6,  28)=           3
 
      flst_real(   1,  29)=           0
      flst_real(   2,  29)=          -5
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=         -11
      flst_real(   5,  29)=          11
      flst_real(   6,  29)=          -5
 
      flst_real(   1,  30)=           0
      flst_real(   2,  30)=           5
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=         -11
      flst_real(   5,  30)=          11
      flst_real(   6,  30)=           5
 
      flst_nreal=          30
 
      return
      end
 
