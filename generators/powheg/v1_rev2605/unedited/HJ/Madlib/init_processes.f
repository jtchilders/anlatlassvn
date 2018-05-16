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
 
      return
      end
 
 
 
      subroutine init_processes_born
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_born(   1,   1)=          -1
      flst_born(   2,   1)=           1
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=           0
 
      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           0
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=          -1
 
      flst_born(   1,   3)=           1
      flst_born(   2,   3)=          -1
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=           0
 
      flst_born(   1,   4)=           1
      flst_born(   2,   4)=           0
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=           1
 
      flst_born(   1,   5)=          -2
      flst_born(   2,   5)=           2
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=           0
 
      flst_born(   1,   6)=          -2
      flst_born(   2,   6)=           0
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=          -2
 
      flst_born(   1,   7)=           2
      flst_born(   2,   7)=          -2
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=           0
 
      flst_born(   1,   8)=           2
      flst_born(   2,   8)=           0
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=           2
 
      flst_born(   1,   9)=          -4
      flst_born(   2,   9)=           4
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=           0
 
      flst_born(   1,  10)=          -4
      flst_born(   2,  10)=           0
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=          -4
 
      flst_born(   1,  11)=           4
      flst_born(   2,  11)=          -4
      flst_born(   3,  11)=          25
      flst_born(   4,  11)=           0
 
      flst_born(   1,  12)=           4
      flst_born(   2,  12)=           0
      flst_born(   3,  12)=          25
      flst_born(   4,  12)=           4
 
      flst_born(   1,  13)=          -3
      flst_born(   2,  13)=           3
      flst_born(   3,  13)=          25
      flst_born(   4,  13)=           0
 
      flst_born(   1,  14)=          -3
      flst_born(   2,  14)=           0
      flst_born(   3,  14)=          25
      flst_born(   4,  14)=          -3
 
      flst_born(   1,  15)=           3
      flst_born(   2,  15)=          -3
      flst_born(   3,  15)=          25
      flst_born(   4,  15)=           0
 
      flst_born(   1,  16)=           3
      flst_born(   2,  16)=           0
      flst_born(   3,  16)=          25
      flst_born(   4,  16)=           3
 
      flst_born(   1,  17)=          -5
      flst_born(   2,  17)=           5
      flst_born(   3,  17)=          25
      flst_born(   4,  17)=           0
 
      flst_born(   1,  18)=          -5
      flst_born(   2,  18)=           0
      flst_born(   3,  18)=          25
      flst_born(   4,  18)=          -5
 
      flst_born(   1,  19)=           5
      flst_born(   2,  19)=          -5
      flst_born(   3,  19)=          25
      flst_born(   4,  19)=           0
 
      flst_born(   1,  20)=           5
      flst_born(   2,  20)=           0
      flst_born(   3,  20)=          25
      flst_born(   4,  20)=           5
 
      flst_born(   1,  21)=           0
      flst_born(   2,  21)=          -1
      flst_born(   3,  21)=          25
      flst_born(   4,  21)=          -1
 
      flst_born(   1,  22)=           0
      flst_born(   2,  22)=           1
      flst_born(   3,  22)=          25
      flst_born(   4,  22)=           1
 
      flst_born(   1,  23)=           0
      flst_born(   2,  23)=          -2
      flst_born(   3,  23)=          25
      flst_born(   4,  23)=          -2
 
      flst_born(   1,  24)=           0
      flst_born(   2,  24)=           2
      flst_born(   3,  24)=          25
      flst_born(   4,  24)=           2
 
      flst_born(   1,  25)=           0
      flst_born(   2,  25)=          -4
      flst_born(   3,  25)=          25
      flst_born(   4,  25)=          -4
 
      flst_born(   1,  26)=           0
      flst_born(   2,  26)=           4
      flst_born(   3,  26)=          25
      flst_born(   4,  26)=           4
 
      flst_born(   1,  27)=           0
      flst_born(   2,  27)=          -3
      flst_born(   3,  27)=          25
      flst_born(   4,  27)=          -3
 
      flst_born(   1,  28)=           0
      flst_born(   2,  28)=           3
      flst_born(   3,  28)=          25
      flst_born(   4,  28)=           3
 
      flst_born(   1,  29)=           0
      flst_born(   2,  29)=          -5
      flst_born(   3,  29)=          25
      flst_born(   4,  29)=          -5
 
      flst_born(   1,  30)=           0
      flst_born(   2,  30)=           5
      flst_born(   3,  30)=          25
      flst_born(   4,  30)=           5
 
      flst_born(   1,  31)=           0
      flst_born(   2,  31)=           0
      flst_born(   3,  31)=          25
      flst_born(   4,  31)=           0
 
      flst_nborn=          31
 
      return
      end
 
 
 
      subroutine init_processes_real
      implicit none
      include "nlegborn.h"
      include "pwhg_flst.h"
 
      flst_real(   1,   1)=          -1
      flst_real(   2,   1)=          -1
      flst_real(   3,   1)=          25
      flst_real(   4,   1)=          -1
      flst_real(   5,   1)=          -1
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           1
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=           1
      flst_real(   5,   2)=          -1
 
      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=           1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=           2
      flst_real(   5,   3)=          -2
 
      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           1
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=           4
      flst_real(   5,   4)=          -4
 
      flst_real(   1,   5)=          -1
      flst_real(   2,   5)=           1
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=           3
      flst_real(   5,   5)=          -3
 
      flst_real(   1,   6)=          -1
      flst_real(   2,   6)=           1
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=           5
      flst_real(   5,   6)=          -5
 
      flst_real(   1,   7)=          -1
      flst_real(   2,   7)=           1
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=           0
      flst_real(   5,   7)=           0
 
      flst_real(   1,   8)=          -1
      flst_real(   2,   8)=          -2
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=          -1
      flst_real(   5,   8)=          -2
 
      flst_real(   1,   9)=          -1
      flst_real(   2,   9)=           2
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=          -1
      flst_real(   5,   9)=           2
 
      flst_real(   1,  10)=          -1
      flst_real(   2,  10)=          -4
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=          -1
      flst_real(   5,  10)=          -4
 
      flst_real(   1,  11)=          -1
      flst_real(   2,  11)=           4
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=          -1
      flst_real(   5,  11)=           4
 
      flst_real(   1,  12)=          -1
      flst_real(   2,  12)=          -3
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=          -1
      flst_real(   5,  12)=          -3
 
      flst_real(   1,  13)=          -1
      flst_real(   2,  13)=           3
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=          -1
      flst_real(   5,  13)=           3
 
      flst_real(   1,  14)=          -1
      flst_real(   2,  14)=          -5
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=          -1
      flst_real(   5,  14)=          -5
 
      flst_real(   1,  15)=          -1
      flst_real(   2,  15)=           5
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=          -1
      flst_real(   5,  15)=           5
 
      flst_real(   1,  16)=          -1
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=          -1
      flst_real(   5,  16)=           0
 
      flst_real(   1,  17)=           1
      flst_real(   2,  17)=          -1
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=           1
      flst_real(   5,  17)=          -1
 
      flst_real(   1,  18)=           1
      flst_real(   2,  18)=          -1
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=           2
      flst_real(   5,  18)=          -2
 
      flst_real(   1,  19)=           1
      flst_real(   2,  19)=          -1
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=           4
      flst_real(   5,  19)=          -4
 
      flst_real(   1,  20)=           1
      flst_real(   2,  20)=          -1
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=           3
      flst_real(   5,  20)=          -3
 
      flst_real(   1,  21)=           1
      flst_real(   2,  21)=          -1
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=           5
      flst_real(   5,  21)=          -5
 
      flst_real(   1,  22)=           1
      flst_real(   2,  22)=          -1
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=           0
      flst_real(   5,  22)=           0
 
      flst_real(   1,  23)=           1
      flst_real(   2,  23)=           1
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=           1
      flst_real(   5,  23)=           1
 
      flst_real(   1,  24)=           1
      flst_real(   2,  24)=          -2
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=           1
      flst_real(   5,  24)=          -2
 
      flst_real(   1,  25)=           1
      flst_real(   2,  25)=           2
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=           1
      flst_real(   5,  25)=           2
 
      flst_real(   1,  26)=           1
      flst_real(   2,  26)=          -4
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=           1
      flst_real(   5,  26)=          -4
 
      flst_real(   1,  27)=           1
      flst_real(   2,  27)=           4
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=           1
      flst_real(   5,  27)=           4
 
      flst_real(   1,  28)=           1
      flst_real(   2,  28)=          -3
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=           1
      flst_real(   5,  28)=          -3
 
      flst_real(   1,  29)=           1
      flst_real(   2,  29)=           3
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=           1
      flst_real(   5,  29)=           3
 
      flst_real(   1,  30)=           1
      flst_real(   2,  30)=          -5
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=           1
      flst_real(   5,  30)=          -5
 
      flst_real(   1,  31)=           1
      flst_real(   2,  31)=           5
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=           1
      flst_real(   5,  31)=           5
 
      flst_real(   1,  32)=           1
      flst_real(   2,  32)=           0
      flst_real(   3,  32)=          25
      flst_real(   4,  32)=           1
      flst_real(   5,  32)=           0
 
      flst_real(   1,  33)=          -2
      flst_real(   2,  33)=          -1
      flst_real(   3,  33)=          25
      flst_real(   4,  33)=          -1
      flst_real(   5,  33)=          -2
 
      flst_real(   1,  34)=          -2
      flst_real(   2,  34)=           1
      flst_real(   3,  34)=          25
      flst_real(   4,  34)=           1
      flst_real(   5,  34)=          -2
 
      flst_real(   1,  35)=          -2
      flst_real(   2,  35)=          -2
      flst_real(   3,  35)=          25
      flst_real(   4,  35)=          -2
      flst_real(   5,  35)=          -2
 
      flst_real(   1,  36)=          -2
      flst_real(   2,  36)=           2
      flst_real(   3,  36)=          25
      flst_real(   4,  36)=           1
      flst_real(   5,  36)=          -1
 
      flst_real(   1,  37)=          -2
      flst_real(   2,  37)=           2
      flst_real(   3,  37)=          25
      flst_real(   4,  37)=           2
      flst_real(   5,  37)=          -2
 
      flst_real(   1,  38)=          -2
      flst_real(   2,  38)=           2
      flst_real(   3,  38)=          25
      flst_real(   4,  38)=           4
      flst_real(   5,  38)=          -4
 
      flst_real(   1,  39)=          -2
      flst_real(   2,  39)=           2
      flst_real(   3,  39)=          25
      flst_real(   4,  39)=           3
      flst_real(   5,  39)=          -3
 
      flst_real(   1,  40)=          -2
      flst_real(   2,  40)=           2
      flst_real(   3,  40)=          25
      flst_real(   4,  40)=           5
      flst_real(   5,  40)=          -5
 
      flst_real(   1,  41)=          -2
      flst_real(   2,  41)=           2
      flst_real(   3,  41)=          25
      flst_real(   4,  41)=           0
      flst_real(   5,  41)=           0
 
      flst_real(   1,  42)=          -2
      flst_real(   2,  42)=          -4
      flst_real(   3,  42)=          25
      flst_real(   4,  42)=          -2
      flst_real(   5,  42)=          -4
 
      flst_real(   1,  43)=          -2
      flst_real(   2,  43)=           4
      flst_real(   3,  43)=          25
      flst_real(   4,  43)=          -2
      flst_real(   5,  43)=           4
 
      flst_real(   1,  44)=          -2
      flst_real(   2,  44)=          -3
      flst_real(   3,  44)=          25
      flst_real(   4,  44)=          -2
      flst_real(   5,  44)=          -3
 
      flst_real(   1,  45)=          -2
      flst_real(   2,  45)=           3
      flst_real(   3,  45)=          25
      flst_real(   4,  45)=          -2
      flst_real(   5,  45)=           3
 
      flst_real(   1,  46)=          -2
      flst_real(   2,  46)=          -5
      flst_real(   3,  46)=          25
      flst_real(   4,  46)=          -2
      flst_real(   5,  46)=          -5
 
      flst_real(   1,  47)=          -2
      flst_real(   2,  47)=           5
      flst_real(   3,  47)=          25
      flst_real(   4,  47)=          -2
      flst_real(   5,  47)=           5
 
      flst_real(   1,  48)=          -2
      flst_real(   2,  48)=           0
      flst_real(   3,  48)=          25
      flst_real(   4,  48)=          -2
      flst_real(   5,  48)=           0
 
      flst_real(   1,  49)=           2
      flst_real(   2,  49)=          -1
      flst_real(   3,  49)=          25
      flst_real(   4,  49)=          -1
      flst_real(   5,  49)=           2
 
      flst_real(   1,  50)=           2
      flst_real(   2,  50)=           1
      flst_real(   3,  50)=          25
      flst_real(   4,  50)=           1
      flst_real(   5,  50)=           2
 
      flst_real(   1,  51)=           2
      flst_real(   2,  51)=          -2
      flst_real(   3,  51)=          25
      flst_real(   4,  51)=           1
      flst_real(   5,  51)=          -1
 
      flst_real(   1,  52)=           2
      flst_real(   2,  52)=          -2
      flst_real(   3,  52)=          25
      flst_real(   4,  52)=           2
      flst_real(   5,  52)=          -2
 
      flst_real(   1,  53)=           2
      flst_real(   2,  53)=          -2
      flst_real(   3,  53)=          25
      flst_real(   4,  53)=           4
      flst_real(   5,  53)=          -4
 
      flst_real(   1,  54)=           2
      flst_real(   2,  54)=          -2
      flst_real(   3,  54)=          25
      flst_real(   4,  54)=           3
      flst_real(   5,  54)=          -3
 
      flst_real(   1,  55)=           2
      flst_real(   2,  55)=          -2
      flst_real(   3,  55)=          25
      flst_real(   4,  55)=           5
      flst_real(   5,  55)=          -5
 
      flst_real(   1,  56)=           2
      flst_real(   2,  56)=          -2
      flst_real(   3,  56)=          25
      flst_real(   4,  56)=           0
      flst_real(   5,  56)=           0
 
      flst_real(   1,  57)=           2
      flst_real(   2,  57)=           2
      flst_real(   3,  57)=          25
      flst_real(   4,  57)=           2
      flst_real(   5,  57)=           2
 
      flst_real(   1,  58)=           2
      flst_real(   2,  58)=          -4
      flst_real(   3,  58)=          25
      flst_real(   4,  58)=           2
      flst_real(   5,  58)=          -4
 
      flst_real(   1,  59)=           2
      flst_real(   2,  59)=           4
      flst_real(   3,  59)=          25
      flst_real(   4,  59)=           2
      flst_real(   5,  59)=           4
 
      flst_real(   1,  60)=           2
      flst_real(   2,  60)=          -3
      flst_real(   3,  60)=          25
      flst_real(   4,  60)=           2
      flst_real(   5,  60)=          -3
 
      flst_real(   1,  61)=           2
      flst_real(   2,  61)=           3
      flst_real(   3,  61)=          25
      flst_real(   4,  61)=           2
      flst_real(   5,  61)=           3
 
      flst_real(   1,  62)=           2
      flst_real(   2,  62)=          -5
      flst_real(   3,  62)=          25
      flst_real(   4,  62)=           2
      flst_real(   5,  62)=          -5
 
      flst_real(   1,  63)=           2
      flst_real(   2,  63)=           5
      flst_real(   3,  63)=          25
      flst_real(   4,  63)=           2
      flst_real(   5,  63)=           5
 
      flst_real(   1,  64)=           2
      flst_real(   2,  64)=           0
      flst_real(   3,  64)=          25
      flst_real(   4,  64)=           2
      flst_real(   5,  64)=           0
 
      flst_real(   1,  65)=          -4
      flst_real(   2,  65)=          -1
      flst_real(   3,  65)=          25
      flst_real(   4,  65)=          -1
      flst_real(   5,  65)=          -4
 
      flst_real(   1,  66)=          -4
      flst_real(   2,  66)=           1
      flst_real(   3,  66)=          25
      flst_real(   4,  66)=           1
      flst_real(   5,  66)=          -4
 
      flst_real(   1,  67)=          -4
      flst_real(   2,  67)=          -2
      flst_real(   3,  67)=          25
      flst_real(   4,  67)=          -2
      flst_real(   5,  67)=          -4
 
      flst_real(   1,  68)=          -4
      flst_real(   2,  68)=           2
      flst_real(   3,  68)=          25
      flst_real(   4,  68)=           2
      flst_real(   5,  68)=          -4
 
      flst_real(   1,  69)=          -4
      flst_real(   2,  69)=          -4
      flst_real(   3,  69)=          25
      flst_real(   4,  69)=          -4
      flst_real(   5,  69)=          -4
 
      flst_real(   1,  70)=          -4
      flst_real(   2,  70)=           4
      flst_real(   3,  70)=          25
      flst_real(   4,  70)=           1
      flst_real(   5,  70)=          -1
 
      flst_real(   1,  71)=          -4
      flst_real(   2,  71)=           4
      flst_real(   3,  71)=          25
      flst_real(   4,  71)=           2
      flst_real(   5,  71)=          -2
 
      flst_real(   1,  72)=          -4
      flst_real(   2,  72)=           4
      flst_real(   3,  72)=          25
      flst_real(   4,  72)=           4
      flst_real(   5,  72)=          -4
 
      flst_real(   1,  73)=          -4
      flst_real(   2,  73)=           4
      flst_real(   3,  73)=          25
      flst_real(   4,  73)=           3
      flst_real(   5,  73)=          -3
 
      flst_real(   1,  74)=          -4
      flst_real(   2,  74)=           4
      flst_real(   3,  74)=          25
      flst_real(   4,  74)=           5
      flst_real(   5,  74)=          -5
 
      flst_real(   1,  75)=          -4
      flst_real(   2,  75)=           4
      flst_real(   3,  75)=          25
      flst_real(   4,  75)=           0
      flst_real(   5,  75)=           0
 
      flst_real(   1,  76)=          -4
      flst_real(   2,  76)=          -3
      flst_real(   3,  76)=          25
      flst_real(   4,  76)=          -4
      flst_real(   5,  76)=          -3
 
      flst_real(   1,  77)=          -4
      flst_real(   2,  77)=           3
      flst_real(   3,  77)=          25
      flst_real(   4,  77)=          -4
      flst_real(   5,  77)=           3
 
      flst_real(   1,  78)=          -4
      flst_real(   2,  78)=          -5
      flst_real(   3,  78)=          25
      flst_real(   4,  78)=          -4
      flst_real(   5,  78)=          -5
 
      flst_real(   1,  79)=          -4
      flst_real(   2,  79)=           5
      flst_real(   3,  79)=          25
      flst_real(   4,  79)=          -4
      flst_real(   5,  79)=           5
 
      flst_real(   1,  80)=          -4
      flst_real(   2,  80)=           0
      flst_real(   3,  80)=          25
      flst_real(   4,  80)=          -4
      flst_real(   5,  80)=           0
 
      flst_real(   1,  81)=           4
      flst_real(   2,  81)=          -1
      flst_real(   3,  81)=          25
      flst_real(   4,  81)=          -1
      flst_real(   5,  81)=           4
 
      flst_real(   1,  82)=           4
      flst_real(   2,  82)=           1
      flst_real(   3,  82)=          25
      flst_real(   4,  82)=           1
      flst_real(   5,  82)=           4
 
      flst_real(   1,  83)=           4
      flst_real(   2,  83)=          -2
      flst_real(   3,  83)=          25
      flst_real(   4,  83)=          -2
      flst_real(   5,  83)=           4
 
      flst_real(   1,  84)=           4
      flst_real(   2,  84)=           2
      flst_real(   3,  84)=          25
      flst_real(   4,  84)=           2
      flst_real(   5,  84)=           4
 
      flst_real(   1,  85)=           4
      flst_real(   2,  85)=          -4
      flst_real(   3,  85)=          25
      flst_real(   4,  85)=           1
      flst_real(   5,  85)=          -1
 
      flst_real(   1,  86)=           4
      flst_real(   2,  86)=          -4
      flst_real(   3,  86)=          25
      flst_real(   4,  86)=           2
      flst_real(   5,  86)=          -2
 
      flst_real(   1,  87)=           4
      flst_real(   2,  87)=          -4
      flst_real(   3,  87)=          25
      flst_real(   4,  87)=           4
      flst_real(   5,  87)=          -4
 
      flst_real(   1,  88)=           4
      flst_real(   2,  88)=          -4
      flst_real(   3,  88)=          25
      flst_real(   4,  88)=           3
      flst_real(   5,  88)=          -3
 
      flst_real(   1,  89)=           4
      flst_real(   2,  89)=          -4
      flst_real(   3,  89)=          25
      flst_real(   4,  89)=           5
      flst_real(   5,  89)=          -5
 
      flst_real(   1,  90)=           4
      flst_real(   2,  90)=          -4
      flst_real(   3,  90)=          25
      flst_real(   4,  90)=           0
      flst_real(   5,  90)=           0
 
      flst_real(   1,  91)=           4
      flst_real(   2,  91)=           4
      flst_real(   3,  91)=          25
      flst_real(   4,  91)=           4
      flst_real(   5,  91)=           4
 
      flst_real(   1,  92)=           4
      flst_real(   2,  92)=          -3
      flst_real(   3,  92)=          25
      flst_real(   4,  92)=           4
      flst_real(   5,  92)=          -3
 
      flst_real(   1,  93)=           4
      flst_real(   2,  93)=           3
      flst_real(   3,  93)=          25
      flst_real(   4,  93)=           4
      flst_real(   5,  93)=           3
 
      flst_real(   1,  94)=           4
      flst_real(   2,  94)=          -5
      flst_real(   3,  94)=          25
      flst_real(   4,  94)=           4
      flst_real(   5,  94)=          -5
 
      flst_real(   1,  95)=           4
      flst_real(   2,  95)=           5
      flst_real(   3,  95)=          25
      flst_real(   4,  95)=           4
      flst_real(   5,  95)=           5
 
      flst_real(   1,  96)=           4
      flst_real(   2,  96)=           0
      flst_real(   3,  96)=          25
      flst_real(   4,  96)=           4
      flst_real(   5,  96)=           0
 
      flst_real(   1,  97)=          -3
      flst_real(   2,  97)=          -1
      flst_real(   3,  97)=          25
      flst_real(   4,  97)=          -1
      flst_real(   5,  97)=          -3
 
      flst_real(   1,  98)=          -3
      flst_real(   2,  98)=           1
      flst_real(   3,  98)=          25
      flst_real(   4,  98)=           1
      flst_real(   5,  98)=          -3
 
      flst_real(   1,  99)=          -3
      flst_real(   2,  99)=          -2
      flst_real(   3,  99)=          25
      flst_real(   4,  99)=          -2
      flst_real(   5,  99)=          -3
 
      flst_real(   1, 100)=          -3
      flst_real(   2, 100)=           2
      flst_real(   3, 100)=          25
      flst_real(   4, 100)=           2
      flst_real(   5, 100)=          -3
 
      flst_real(   1, 101)=          -3
      flst_real(   2, 101)=          -4
      flst_real(   3, 101)=          25
      flst_real(   4, 101)=          -4
      flst_real(   5, 101)=          -3
 
      flst_real(   1, 102)=          -3
      flst_real(   2, 102)=           4
      flst_real(   3, 102)=          25
      flst_real(   4, 102)=           4
      flst_real(   5, 102)=          -3
 
      flst_real(   1, 103)=          -3
      flst_real(   2, 103)=          -3
      flst_real(   3, 103)=          25
      flst_real(   4, 103)=          -3
      flst_real(   5, 103)=          -3
 
      flst_real(   1, 104)=          -3
      flst_real(   2, 104)=           3
      flst_real(   3, 104)=          25
      flst_real(   4, 104)=           1
      flst_real(   5, 104)=          -1
 
      flst_real(   1, 105)=          -3
      flst_real(   2, 105)=           3
      flst_real(   3, 105)=          25
      flst_real(   4, 105)=           2
      flst_real(   5, 105)=          -2
 
      flst_real(   1, 106)=          -3
      flst_real(   2, 106)=           3
      flst_real(   3, 106)=          25
      flst_real(   4, 106)=           4
      flst_real(   5, 106)=          -4
 
      flst_real(   1, 107)=          -3
      flst_real(   2, 107)=           3
      flst_real(   3, 107)=          25
      flst_real(   4, 107)=           3
      flst_real(   5, 107)=          -3
 
      flst_real(   1, 108)=          -3
      flst_real(   2, 108)=           3
      flst_real(   3, 108)=          25
      flst_real(   4, 108)=           5
      flst_real(   5, 108)=          -5
 
      flst_real(   1, 109)=          -3
      flst_real(   2, 109)=           3
      flst_real(   3, 109)=          25
      flst_real(   4, 109)=           0
      flst_real(   5, 109)=           0
 
      flst_real(   1, 110)=          -3
      flst_real(   2, 110)=          -5
      flst_real(   3, 110)=          25
      flst_real(   4, 110)=          -3
      flst_real(   5, 110)=          -5
 
      flst_real(   1, 111)=          -3
      flst_real(   2, 111)=           5
      flst_real(   3, 111)=          25
      flst_real(   4, 111)=          -3
      flst_real(   5, 111)=           5
 
      flst_real(   1, 112)=          -3
      flst_real(   2, 112)=           0
      flst_real(   3, 112)=          25
      flst_real(   4, 112)=          -3
      flst_real(   5, 112)=           0
 
      flst_real(   1, 113)=           3
      flst_real(   2, 113)=          -1
      flst_real(   3, 113)=          25
      flst_real(   4, 113)=          -1
      flst_real(   5, 113)=           3
 
      flst_real(   1, 114)=           3
      flst_real(   2, 114)=           1
      flst_real(   3, 114)=          25
      flst_real(   4, 114)=           1
      flst_real(   5, 114)=           3
 
      flst_real(   1, 115)=           3
      flst_real(   2, 115)=          -2
      flst_real(   3, 115)=          25
      flst_real(   4, 115)=          -2
      flst_real(   5, 115)=           3
 
      flst_real(   1, 116)=           3
      flst_real(   2, 116)=           2
      flst_real(   3, 116)=          25
      flst_real(   4, 116)=           2
      flst_real(   5, 116)=           3
 
      flst_real(   1, 117)=           3
      flst_real(   2, 117)=          -4
      flst_real(   3, 117)=          25
      flst_real(   4, 117)=          -4
      flst_real(   5, 117)=           3
 
      flst_real(   1, 118)=           3
      flst_real(   2, 118)=           4
      flst_real(   3, 118)=          25
      flst_real(   4, 118)=           4
      flst_real(   5, 118)=           3
 
      flst_real(   1, 119)=           3
      flst_real(   2, 119)=          -3
      flst_real(   3, 119)=          25
      flst_real(   4, 119)=           1
      flst_real(   5, 119)=          -1
 
      flst_real(   1, 120)=           3
      flst_real(   2, 120)=          -3
      flst_real(   3, 120)=          25
      flst_real(   4, 120)=           2
      flst_real(   5, 120)=          -2
 
      flst_real(   1, 121)=           3
      flst_real(   2, 121)=          -3
      flst_real(   3, 121)=          25
      flst_real(   4, 121)=           4
      flst_real(   5, 121)=          -4
 
      flst_real(   1, 122)=           3
      flst_real(   2, 122)=          -3
      flst_real(   3, 122)=          25
      flst_real(   4, 122)=           3
      flst_real(   5, 122)=          -3
 
      flst_real(   1, 123)=           3
      flst_real(   2, 123)=          -3
      flst_real(   3, 123)=          25
      flst_real(   4, 123)=           5
      flst_real(   5, 123)=          -5
 
      flst_real(   1, 124)=           3
      flst_real(   2, 124)=          -3
      flst_real(   3, 124)=          25
      flst_real(   4, 124)=           0
      flst_real(   5, 124)=           0
 
      flst_real(   1, 125)=           3
      flst_real(   2, 125)=           3
      flst_real(   3, 125)=          25
      flst_real(   4, 125)=           3
      flst_real(   5, 125)=           3
 
      flst_real(   1, 126)=           3
      flst_real(   2, 126)=          -5
      flst_real(   3, 126)=          25
      flst_real(   4, 126)=           3
      flst_real(   5, 126)=          -5
 
      flst_real(   1, 127)=           3
      flst_real(   2, 127)=           5
      flst_real(   3, 127)=          25
      flst_real(   4, 127)=           3
      flst_real(   5, 127)=           5
 
      flst_real(   1, 128)=           3
      flst_real(   2, 128)=           0
      flst_real(   3, 128)=          25
      flst_real(   4, 128)=           3
      flst_real(   5, 128)=           0
 
      flst_real(   1, 129)=          -5
      flst_real(   2, 129)=          -1
      flst_real(   3, 129)=          25
      flst_real(   4, 129)=          -1
      flst_real(   5, 129)=          -5
 
      flst_real(   1, 130)=          -5
      flst_real(   2, 130)=           1
      flst_real(   3, 130)=          25
      flst_real(   4, 130)=           1
      flst_real(   5, 130)=          -5
 
      flst_real(   1, 131)=          -5
      flst_real(   2, 131)=          -2
      flst_real(   3, 131)=          25
      flst_real(   4, 131)=          -2
      flst_real(   5, 131)=          -5
 
      flst_real(   1, 132)=          -5
      flst_real(   2, 132)=           2
      flst_real(   3, 132)=          25
      flst_real(   4, 132)=           2
      flst_real(   5, 132)=          -5
 
      flst_real(   1, 133)=          -5
      flst_real(   2, 133)=          -4
      flst_real(   3, 133)=          25
      flst_real(   4, 133)=          -4
      flst_real(   5, 133)=          -5
 
      flst_real(   1, 134)=          -5
      flst_real(   2, 134)=           4
      flst_real(   3, 134)=          25
      flst_real(   4, 134)=           4
      flst_real(   5, 134)=          -5
 
      flst_real(   1, 135)=          -5
      flst_real(   2, 135)=          -3
      flst_real(   3, 135)=          25
      flst_real(   4, 135)=          -3
      flst_real(   5, 135)=          -5
 
      flst_real(   1, 136)=          -5
      flst_real(   2, 136)=           3
      flst_real(   3, 136)=          25
      flst_real(   4, 136)=           3
      flst_real(   5, 136)=          -5
 
      flst_real(   1, 137)=          -5
      flst_real(   2, 137)=          -5
      flst_real(   3, 137)=          25
      flst_real(   4, 137)=          -5
      flst_real(   5, 137)=          -5
 
      flst_real(   1, 138)=          -5
      flst_real(   2, 138)=           5
      flst_real(   3, 138)=          25
      flst_real(   4, 138)=           1
      flst_real(   5, 138)=          -1
 
      flst_real(   1, 139)=          -5
      flst_real(   2, 139)=           5
      flst_real(   3, 139)=          25
      flst_real(   4, 139)=           2
      flst_real(   5, 139)=          -2
 
      flst_real(   1, 140)=          -5
      flst_real(   2, 140)=           5
      flst_real(   3, 140)=          25
      flst_real(   4, 140)=           4
      flst_real(   5, 140)=          -4
 
      flst_real(   1, 141)=          -5
      flst_real(   2, 141)=           5
      flst_real(   3, 141)=          25
      flst_real(   4, 141)=           3
      flst_real(   5, 141)=          -3
 
      flst_real(   1, 142)=          -5
      flst_real(   2, 142)=           5
      flst_real(   3, 142)=          25
      flst_real(   4, 142)=           5
      flst_real(   5, 142)=          -5
 
      flst_real(   1, 143)=          -5
      flst_real(   2, 143)=           5
      flst_real(   3, 143)=          25
      flst_real(   4, 143)=           0
      flst_real(   5, 143)=           0
 
      flst_real(   1, 144)=          -5
      flst_real(   2, 144)=           0
      flst_real(   3, 144)=          25
      flst_real(   4, 144)=          -5
      flst_real(   5, 144)=           0
 
      flst_real(   1, 145)=           5
      flst_real(   2, 145)=          -1
      flst_real(   3, 145)=          25
      flst_real(   4, 145)=          -1
      flst_real(   5, 145)=           5
 
      flst_real(   1, 146)=           5
      flst_real(   2, 146)=           1
      flst_real(   3, 146)=          25
      flst_real(   4, 146)=           1
      flst_real(   5, 146)=           5
 
      flst_real(   1, 147)=           5
      flst_real(   2, 147)=          -2
      flst_real(   3, 147)=          25
      flst_real(   4, 147)=          -2
      flst_real(   5, 147)=           5
 
      flst_real(   1, 148)=           5
      flst_real(   2, 148)=           2
      flst_real(   3, 148)=          25
      flst_real(   4, 148)=           2
      flst_real(   5, 148)=           5
 
      flst_real(   1, 149)=           5
      flst_real(   2, 149)=          -4
      flst_real(   3, 149)=          25
      flst_real(   4, 149)=          -4
      flst_real(   5, 149)=           5
 
      flst_real(   1, 150)=           5
      flst_real(   2, 150)=           4
      flst_real(   3, 150)=          25
      flst_real(   4, 150)=           4
      flst_real(   5, 150)=           5
 
      flst_real(   1, 151)=           5
      flst_real(   2, 151)=          -3
      flst_real(   3, 151)=          25
      flst_real(   4, 151)=          -3
      flst_real(   5, 151)=           5
 
      flst_real(   1, 152)=           5
      flst_real(   2, 152)=           3
      flst_real(   3, 152)=          25
      flst_real(   4, 152)=           3
      flst_real(   5, 152)=           5
 
      flst_real(   1, 153)=           5
      flst_real(   2, 153)=          -5
      flst_real(   3, 153)=          25
      flst_real(   4, 153)=           1
      flst_real(   5, 153)=          -1
 
      flst_real(   1, 154)=           5
      flst_real(   2, 154)=          -5
      flst_real(   3, 154)=          25
      flst_real(   4, 154)=           2
      flst_real(   5, 154)=          -2
 
      flst_real(   1, 155)=           5
      flst_real(   2, 155)=          -5
      flst_real(   3, 155)=          25
      flst_real(   4, 155)=           4
      flst_real(   5, 155)=          -4
 
      flst_real(   1, 156)=           5
      flst_real(   2, 156)=          -5
      flst_real(   3, 156)=          25
      flst_real(   4, 156)=           3
      flst_real(   5, 156)=          -3
 
      flst_real(   1, 157)=           5
      flst_real(   2, 157)=          -5
      flst_real(   3, 157)=          25
      flst_real(   4, 157)=           5
      flst_real(   5, 157)=          -5
 
      flst_real(   1, 158)=           5
      flst_real(   2, 158)=          -5
      flst_real(   3, 158)=          25
      flst_real(   4, 158)=           0
      flst_real(   5, 158)=           0
 
      flst_real(   1, 159)=           5
      flst_real(   2, 159)=           5
      flst_real(   3, 159)=          25
      flst_real(   4, 159)=           5
      flst_real(   5, 159)=           5
 
      flst_real(   1, 160)=           5
      flst_real(   2, 160)=           0
      flst_real(   3, 160)=          25
      flst_real(   4, 160)=           5
      flst_real(   5, 160)=           0
 
      flst_real(   1, 161)=           0
      flst_real(   2, 161)=          -1
      flst_real(   3, 161)=          25
      flst_real(   4, 161)=          -1
      flst_real(   5, 161)=           0
 
      flst_real(   1, 162)=           0
      flst_real(   2, 162)=           1
      flst_real(   3, 162)=          25
      flst_real(   4, 162)=           1
      flst_real(   5, 162)=           0
 
      flst_real(   1, 163)=           0
      flst_real(   2, 163)=          -2
      flst_real(   3, 163)=          25
      flst_real(   4, 163)=          -2
      flst_real(   5, 163)=           0
 
      flst_real(   1, 164)=           0
      flst_real(   2, 164)=           2
      flst_real(   3, 164)=          25
      flst_real(   4, 164)=           2
      flst_real(   5, 164)=           0
 
      flst_real(   1, 165)=           0
      flst_real(   2, 165)=          -4
      flst_real(   3, 165)=          25
      flst_real(   4, 165)=          -4
      flst_real(   5, 165)=           0
 
      flst_real(   1, 166)=           0
      flst_real(   2, 166)=           4
      flst_real(   3, 166)=          25
      flst_real(   4, 166)=           4
      flst_real(   5, 166)=           0
 
      flst_real(   1, 167)=           0
      flst_real(   2, 167)=          -3
      flst_real(   3, 167)=          25
      flst_real(   4, 167)=          -3
      flst_real(   5, 167)=           0
 
      flst_real(   1, 168)=           0
      flst_real(   2, 168)=           3
      flst_real(   3, 168)=          25
      flst_real(   4, 168)=           3
      flst_real(   5, 168)=           0
 
      flst_real(   1, 169)=           0
      flst_real(   2, 169)=          -5
      flst_real(   3, 169)=          25
      flst_real(   4, 169)=          -5
      flst_real(   5, 169)=           0
 
      flst_real(   1, 170)=           0
      flst_real(   2, 170)=           5
      flst_real(   3, 170)=          25
      flst_real(   4, 170)=           5
      flst_real(   5, 170)=           0
 
      flst_real(   1, 171)=           0
      flst_real(   2, 171)=           0
      flst_real(   3, 171)=          25
      flst_real(   4, 171)=           1
      flst_real(   5, 171)=          -1
 
      flst_real(   1, 172)=           0
      flst_real(   2, 172)=           0
      flst_real(   3, 172)=          25
      flst_real(   4, 172)=           2
      flst_real(   5, 172)=          -2
 
      flst_real(   1, 173)=           0
      flst_real(   2, 173)=           0
      flst_real(   3, 173)=          25
      flst_real(   4, 173)=           4
      flst_real(   5, 173)=          -4
 
      flst_real(   1, 174)=           0
      flst_real(   2, 174)=           0
      flst_real(   3, 174)=          25
      flst_real(   4, 174)=           3
      flst_real(   5, 174)=          -3
 
      flst_real(   1, 175)=           0
      flst_real(   2, 175)=           0
      flst_real(   3, 175)=          25
      flst_real(   4, 175)=           5
      flst_real(   5, 175)=          -5
 
      flst_real(   1, 176)=           0
      flst_real(   2, 176)=           0
      flst_real(   3, 176)=          25
      flst_real(   4, 176)=           0
      flst_real(   5, 176)=           0
 
      flst_nreal=         176
 
      return
      end
 
