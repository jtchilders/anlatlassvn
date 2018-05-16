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
      flst_born(   2,   1)=          -1
      flst_born(   3,   1)=          25
      flst_born(   4,   1)=          -1
      flst_born(   5,   1)=          -1
 
      flst_born(   1,   2)=          -1
      flst_born(   2,   2)=           1
      flst_born(   3,   2)=          25
      flst_born(   4,   2)=           1
      flst_born(   5,   2)=          -1
 
      flst_born(   1,   3)=          -1
      flst_born(   2,   3)=           1
      flst_born(   3,   3)=          25
      flst_born(   4,   3)=           2
      flst_born(   5,   3)=          -2
 
      flst_born(   1,   4)=          -1
      flst_born(   2,   4)=           1
      flst_born(   3,   4)=          25
      flst_born(   4,   4)=           4
      flst_born(   5,   4)=          -4
 
      flst_born(   1,   5)=          -1
      flst_born(   2,   5)=           1
      flst_born(   3,   5)=          25
      flst_born(   4,   5)=           3
      flst_born(   5,   5)=          -3
 
      flst_born(   1,   6)=          -1
      flst_born(   2,   6)=           1
      flst_born(   3,   6)=          25
      flst_born(   4,   6)=           5
      flst_born(   5,   6)=          -5
 
      flst_born(   1,   7)=          -1
      flst_born(   2,   7)=           1
      flst_born(   3,   7)=          25
      flst_born(   4,   7)=           0
      flst_born(   5,   7)=           0
 
      flst_born(   1,   8)=          -1
      flst_born(   2,   8)=          -2
      flst_born(   3,   8)=          25
      flst_born(   4,   8)=          -1
      flst_born(   5,   8)=          -2
 
      flst_born(   1,   9)=          -1
      flst_born(   2,   9)=           2
      flst_born(   3,   9)=          25
      flst_born(   4,   9)=          -1
      flst_born(   5,   9)=           2
 
      flst_born(   1,  10)=          -1
      flst_born(   2,  10)=          -4
      flst_born(   3,  10)=          25
      flst_born(   4,  10)=          -1
      flst_born(   5,  10)=          -4
 
      flst_born(   1,  11)=          -1
      flst_born(   2,  11)=           4
      flst_born(   3,  11)=          25
      flst_born(   4,  11)=          -1
      flst_born(   5,  11)=           4
 
      flst_born(   1,  12)=          -1
      flst_born(   2,  12)=          -3
      flst_born(   3,  12)=          25
      flst_born(   4,  12)=          -1
      flst_born(   5,  12)=          -3
 
      flst_born(   1,  13)=          -1
      flst_born(   2,  13)=           3
      flst_born(   3,  13)=          25
      flst_born(   4,  13)=          -1
      flst_born(   5,  13)=           3
 
      flst_born(   1,  14)=          -1
      flst_born(   2,  14)=          -5
      flst_born(   3,  14)=          25
      flst_born(   4,  14)=          -1
      flst_born(   5,  14)=          -5
 
      flst_born(   1,  15)=          -1
      flst_born(   2,  15)=           5
      flst_born(   3,  15)=          25
      flst_born(   4,  15)=          -1
      flst_born(   5,  15)=           5
 
      flst_born(   1,  16)=          -1
      flst_born(   2,  16)=           0
      flst_born(   3,  16)=          25
      flst_born(   4,  16)=          -1
      flst_born(   5,  16)=           0
 
      flst_born(   1,  17)=           1
      flst_born(   2,  17)=          -1
      flst_born(   3,  17)=          25
      flst_born(   4,  17)=           1
      flst_born(   5,  17)=          -1
 
      flst_born(   1,  18)=           1
      flst_born(   2,  18)=          -1
      flst_born(   3,  18)=          25
      flst_born(   4,  18)=           2
      flst_born(   5,  18)=          -2
 
      flst_born(   1,  19)=           1
      flst_born(   2,  19)=          -1
      flst_born(   3,  19)=          25
      flst_born(   4,  19)=           4
      flst_born(   5,  19)=          -4
 
      flst_born(   1,  20)=           1
      flst_born(   2,  20)=          -1
      flst_born(   3,  20)=          25
      flst_born(   4,  20)=           3
      flst_born(   5,  20)=          -3
 
      flst_born(   1,  21)=           1
      flst_born(   2,  21)=          -1
      flst_born(   3,  21)=          25
      flst_born(   4,  21)=           5
      flst_born(   5,  21)=          -5
 
      flst_born(   1,  22)=           1
      flst_born(   2,  22)=          -1
      flst_born(   3,  22)=          25
      flst_born(   4,  22)=           0
      flst_born(   5,  22)=           0
 
      flst_born(   1,  23)=           1
      flst_born(   2,  23)=           1
      flst_born(   3,  23)=          25
      flst_born(   4,  23)=           1
      flst_born(   5,  23)=           1
 
      flst_born(   1,  24)=           1
      flst_born(   2,  24)=          -2
      flst_born(   3,  24)=          25
      flst_born(   4,  24)=           1
      flst_born(   5,  24)=          -2
 
      flst_born(   1,  25)=           1
      flst_born(   2,  25)=           2
      flst_born(   3,  25)=          25
      flst_born(   4,  25)=           1
      flst_born(   5,  25)=           2
 
      flst_born(   1,  26)=           1
      flst_born(   2,  26)=          -4
      flst_born(   3,  26)=          25
      flst_born(   4,  26)=           1
      flst_born(   5,  26)=          -4
 
      flst_born(   1,  27)=           1
      flst_born(   2,  27)=           4
      flst_born(   3,  27)=          25
      flst_born(   4,  27)=           1
      flst_born(   5,  27)=           4
 
      flst_born(   1,  28)=           1
      flst_born(   2,  28)=          -3
      flst_born(   3,  28)=          25
      flst_born(   4,  28)=           1
      flst_born(   5,  28)=          -3
 
      flst_born(   1,  29)=           1
      flst_born(   2,  29)=           3
      flst_born(   3,  29)=          25
      flst_born(   4,  29)=           1
      flst_born(   5,  29)=           3
 
      flst_born(   1,  30)=           1
      flst_born(   2,  30)=          -5
      flst_born(   3,  30)=          25
      flst_born(   4,  30)=           1
      flst_born(   5,  30)=          -5
 
      flst_born(   1,  31)=           1
      flst_born(   2,  31)=           5
      flst_born(   3,  31)=          25
      flst_born(   4,  31)=           1
      flst_born(   5,  31)=           5
 
      flst_born(   1,  32)=           1
      flst_born(   2,  32)=           0
      flst_born(   3,  32)=          25
      flst_born(   4,  32)=           1
      flst_born(   5,  32)=           0
 
      flst_born(   1,  33)=          -2
      flst_born(   2,  33)=          -1
      flst_born(   3,  33)=          25
      flst_born(   4,  33)=          -1
      flst_born(   5,  33)=          -2
 
      flst_born(   1,  34)=          -2
      flst_born(   2,  34)=           1
      flst_born(   3,  34)=          25
      flst_born(   4,  34)=           1
      flst_born(   5,  34)=          -2
 
      flst_born(   1,  35)=          -2
      flst_born(   2,  35)=          -2
      flst_born(   3,  35)=          25
      flst_born(   4,  35)=          -2
      flst_born(   5,  35)=          -2
 
      flst_born(   1,  36)=          -2
      flst_born(   2,  36)=           2
      flst_born(   3,  36)=          25
      flst_born(   4,  36)=           1
      flst_born(   5,  36)=          -1
 
      flst_born(   1,  37)=          -2
      flst_born(   2,  37)=           2
      flst_born(   3,  37)=          25
      flst_born(   4,  37)=           2
      flst_born(   5,  37)=          -2
 
      flst_born(   1,  38)=          -2
      flst_born(   2,  38)=           2
      flst_born(   3,  38)=          25
      flst_born(   4,  38)=           4
      flst_born(   5,  38)=          -4
 
      flst_born(   1,  39)=          -2
      flst_born(   2,  39)=           2
      flst_born(   3,  39)=          25
      flst_born(   4,  39)=           3
      flst_born(   5,  39)=          -3
 
      flst_born(   1,  40)=          -2
      flst_born(   2,  40)=           2
      flst_born(   3,  40)=          25
      flst_born(   4,  40)=           5
      flst_born(   5,  40)=          -5
 
      flst_born(   1,  41)=          -2
      flst_born(   2,  41)=           2
      flst_born(   3,  41)=          25
      flst_born(   4,  41)=           0
      flst_born(   5,  41)=           0
 
      flst_born(   1,  42)=          -2
      flst_born(   2,  42)=          -4
      flst_born(   3,  42)=          25
      flst_born(   4,  42)=          -2
      flst_born(   5,  42)=          -4
 
      flst_born(   1,  43)=          -2
      flst_born(   2,  43)=           4
      flst_born(   3,  43)=          25
      flst_born(   4,  43)=          -2
      flst_born(   5,  43)=           4
 
      flst_born(   1,  44)=          -2
      flst_born(   2,  44)=          -3
      flst_born(   3,  44)=          25
      flst_born(   4,  44)=          -2
      flst_born(   5,  44)=          -3
 
      flst_born(   1,  45)=          -2
      flst_born(   2,  45)=           3
      flst_born(   3,  45)=          25
      flst_born(   4,  45)=          -2
      flst_born(   5,  45)=           3
 
      flst_born(   1,  46)=          -2
      flst_born(   2,  46)=          -5
      flst_born(   3,  46)=          25
      flst_born(   4,  46)=          -2
      flst_born(   5,  46)=          -5
 
      flst_born(   1,  47)=          -2
      flst_born(   2,  47)=           5
      flst_born(   3,  47)=          25
      flst_born(   4,  47)=          -2
      flst_born(   5,  47)=           5
 
      flst_born(   1,  48)=          -2
      flst_born(   2,  48)=           0
      flst_born(   3,  48)=          25
      flst_born(   4,  48)=          -2
      flst_born(   5,  48)=           0
 
      flst_born(   1,  49)=           2
      flst_born(   2,  49)=          -1
      flst_born(   3,  49)=          25
      flst_born(   4,  49)=          -1
      flst_born(   5,  49)=           2
 
      flst_born(   1,  50)=           2
      flst_born(   2,  50)=           1
      flst_born(   3,  50)=          25
      flst_born(   4,  50)=           1
      flst_born(   5,  50)=           2
 
      flst_born(   1,  51)=           2
      flst_born(   2,  51)=          -2
      flst_born(   3,  51)=          25
      flst_born(   4,  51)=           1
      flst_born(   5,  51)=          -1
 
      flst_born(   1,  52)=           2
      flst_born(   2,  52)=          -2
      flst_born(   3,  52)=          25
      flst_born(   4,  52)=           2
      flst_born(   5,  52)=          -2
 
      flst_born(   1,  53)=           2
      flst_born(   2,  53)=          -2
      flst_born(   3,  53)=          25
      flst_born(   4,  53)=           4
      flst_born(   5,  53)=          -4
 
      flst_born(   1,  54)=           2
      flst_born(   2,  54)=          -2
      flst_born(   3,  54)=          25
      flst_born(   4,  54)=           3
      flst_born(   5,  54)=          -3
 
      flst_born(   1,  55)=           2
      flst_born(   2,  55)=          -2
      flst_born(   3,  55)=          25
      flst_born(   4,  55)=           5
      flst_born(   5,  55)=          -5
 
      flst_born(   1,  56)=           2
      flst_born(   2,  56)=          -2
      flst_born(   3,  56)=          25
      flst_born(   4,  56)=           0
      flst_born(   5,  56)=           0
 
      flst_born(   1,  57)=           2
      flst_born(   2,  57)=           2
      flst_born(   3,  57)=          25
      flst_born(   4,  57)=           2
      flst_born(   5,  57)=           2
 
      flst_born(   1,  58)=           2
      flst_born(   2,  58)=          -4
      flst_born(   3,  58)=          25
      flst_born(   4,  58)=           2
      flst_born(   5,  58)=          -4
 
      flst_born(   1,  59)=           2
      flst_born(   2,  59)=           4
      flst_born(   3,  59)=          25
      flst_born(   4,  59)=           2
      flst_born(   5,  59)=           4
 
      flst_born(   1,  60)=           2
      flst_born(   2,  60)=          -3
      flst_born(   3,  60)=          25
      flst_born(   4,  60)=           2
      flst_born(   5,  60)=          -3
 
      flst_born(   1,  61)=           2
      flst_born(   2,  61)=           3
      flst_born(   3,  61)=          25
      flst_born(   4,  61)=           2
      flst_born(   5,  61)=           3
 
      flst_born(   1,  62)=           2
      flst_born(   2,  62)=          -5
      flst_born(   3,  62)=          25
      flst_born(   4,  62)=           2
      flst_born(   5,  62)=          -5
 
      flst_born(   1,  63)=           2
      flst_born(   2,  63)=           5
      flst_born(   3,  63)=          25
      flst_born(   4,  63)=           2
      flst_born(   5,  63)=           5
 
      flst_born(   1,  64)=           2
      flst_born(   2,  64)=           0
      flst_born(   3,  64)=          25
      flst_born(   4,  64)=           2
      flst_born(   5,  64)=           0
 
      flst_born(   1,  65)=          -4
      flst_born(   2,  65)=          -1
      flst_born(   3,  65)=          25
      flst_born(   4,  65)=          -1
      flst_born(   5,  65)=          -4
 
      flst_born(   1,  66)=          -4
      flst_born(   2,  66)=           1
      flst_born(   3,  66)=          25
      flst_born(   4,  66)=           1
      flst_born(   5,  66)=          -4
 
      flst_born(   1,  67)=          -4
      flst_born(   2,  67)=          -2
      flst_born(   3,  67)=          25
      flst_born(   4,  67)=          -2
      flst_born(   5,  67)=          -4
 
      flst_born(   1,  68)=          -4
      flst_born(   2,  68)=           2
      flst_born(   3,  68)=          25
      flst_born(   4,  68)=           2
      flst_born(   5,  68)=          -4
 
      flst_born(   1,  69)=          -4
      flst_born(   2,  69)=          -4
      flst_born(   3,  69)=          25
      flst_born(   4,  69)=          -4
      flst_born(   5,  69)=          -4
 
      flst_born(   1,  70)=          -4
      flst_born(   2,  70)=           4
      flst_born(   3,  70)=          25
      flst_born(   4,  70)=           1
      flst_born(   5,  70)=          -1
 
      flst_born(   1,  71)=          -4
      flst_born(   2,  71)=           4
      flst_born(   3,  71)=          25
      flst_born(   4,  71)=           2
      flst_born(   5,  71)=          -2
 
      flst_born(   1,  72)=          -4
      flst_born(   2,  72)=           4
      flst_born(   3,  72)=          25
      flst_born(   4,  72)=           4
      flst_born(   5,  72)=          -4
 
      flst_born(   1,  73)=          -4
      flst_born(   2,  73)=           4
      flst_born(   3,  73)=          25
      flst_born(   4,  73)=           3
      flst_born(   5,  73)=          -3
 
      flst_born(   1,  74)=          -4
      flst_born(   2,  74)=           4
      flst_born(   3,  74)=          25
      flst_born(   4,  74)=           5
      flst_born(   5,  74)=          -5
 
      flst_born(   1,  75)=          -4
      flst_born(   2,  75)=           4
      flst_born(   3,  75)=          25
      flst_born(   4,  75)=           0
      flst_born(   5,  75)=           0
 
      flst_born(   1,  76)=          -4
      flst_born(   2,  76)=          -3
      flst_born(   3,  76)=          25
      flst_born(   4,  76)=          -4
      flst_born(   5,  76)=          -3
 
      flst_born(   1,  77)=          -4
      flst_born(   2,  77)=           3
      flst_born(   3,  77)=          25
      flst_born(   4,  77)=          -4
      flst_born(   5,  77)=           3
 
      flst_born(   1,  78)=          -4
      flst_born(   2,  78)=          -5
      flst_born(   3,  78)=          25
      flst_born(   4,  78)=          -4
      flst_born(   5,  78)=          -5
 
      flst_born(   1,  79)=          -4
      flst_born(   2,  79)=           5
      flst_born(   3,  79)=          25
      flst_born(   4,  79)=          -4
      flst_born(   5,  79)=           5
 
      flst_born(   1,  80)=          -4
      flst_born(   2,  80)=           0
      flst_born(   3,  80)=          25
      flst_born(   4,  80)=          -4
      flst_born(   5,  80)=           0
 
      flst_born(   1,  81)=           4
      flst_born(   2,  81)=          -1
      flst_born(   3,  81)=          25
      flst_born(   4,  81)=          -1
      flst_born(   5,  81)=           4
 
      flst_born(   1,  82)=           4
      flst_born(   2,  82)=           1
      flst_born(   3,  82)=          25
      flst_born(   4,  82)=           1
      flst_born(   5,  82)=           4
 
      flst_born(   1,  83)=           4
      flst_born(   2,  83)=          -2
      flst_born(   3,  83)=          25
      flst_born(   4,  83)=          -2
      flst_born(   5,  83)=           4
 
      flst_born(   1,  84)=           4
      flst_born(   2,  84)=           2
      flst_born(   3,  84)=          25
      flst_born(   4,  84)=           2
      flst_born(   5,  84)=           4
 
      flst_born(   1,  85)=           4
      flst_born(   2,  85)=          -4
      flst_born(   3,  85)=          25
      flst_born(   4,  85)=           1
      flst_born(   5,  85)=          -1
 
      flst_born(   1,  86)=           4
      flst_born(   2,  86)=          -4
      flst_born(   3,  86)=          25
      flst_born(   4,  86)=           2
      flst_born(   5,  86)=          -2
 
      flst_born(   1,  87)=           4
      flst_born(   2,  87)=          -4
      flst_born(   3,  87)=          25
      flst_born(   4,  87)=           4
      flst_born(   5,  87)=          -4
 
      flst_born(   1,  88)=           4
      flst_born(   2,  88)=          -4
      flst_born(   3,  88)=          25
      flst_born(   4,  88)=           3
      flst_born(   5,  88)=          -3
 
      flst_born(   1,  89)=           4
      flst_born(   2,  89)=          -4
      flst_born(   3,  89)=          25
      flst_born(   4,  89)=           5
      flst_born(   5,  89)=          -5
 
      flst_born(   1,  90)=           4
      flst_born(   2,  90)=          -4
      flst_born(   3,  90)=          25
      flst_born(   4,  90)=           0
      flst_born(   5,  90)=           0
 
      flst_born(   1,  91)=           4
      flst_born(   2,  91)=           4
      flst_born(   3,  91)=          25
      flst_born(   4,  91)=           4
      flst_born(   5,  91)=           4
 
      flst_born(   1,  92)=           4
      flst_born(   2,  92)=          -3
      flst_born(   3,  92)=          25
      flst_born(   4,  92)=           4
      flst_born(   5,  92)=          -3
 
      flst_born(   1,  93)=           4
      flst_born(   2,  93)=           3
      flst_born(   3,  93)=          25
      flst_born(   4,  93)=           4
      flst_born(   5,  93)=           3
 
      flst_born(   1,  94)=           4
      flst_born(   2,  94)=          -5
      flst_born(   3,  94)=          25
      flst_born(   4,  94)=           4
      flst_born(   5,  94)=          -5
 
      flst_born(   1,  95)=           4
      flst_born(   2,  95)=           5
      flst_born(   3,  95)=          25
      flst_born(   4,  95)=           4
      flst_born(   5,  95)=           5
 
      flst_born(   1,  96)=           4
      flst_born(   2,  96)=           0
      flst_born(   3,  96)=          25
      flst_born(   4,  96)=           4
      flst_born(   5,  96)=           0
 
      flst_born(   1,  97)=          -3
      flst_born(   2,  97)=          -1
      flst_born(   3,  97)=          25
      flst_born(   4,  97)=          -1
      flst_born(   5,  97)=          -3
 
      flst_born(   1,  98)=          -3
      flst_born(   2,  98)=           1
      flst_born(   3,  98)=          25
      flst_born(   4,  98)=           1
      flst_born(   5,  98)=          -3
 
      flst_born(   1,  99)=          -3
      flst_born(   2,  99)=          -2
      flst_born(   3,  99)=          25
      flst_born(   4,  99)=          -2
      flst_born(   5,  99)=          -3
 
      flst_born(   1, 100)=          -3
      flst_born(   2, 100)=           2
      flst_born(   3, 100)=          25
      flst_born(   4, 100)=           2
      flst_born(   5, 100)=          -3
 
      flst_born(   1, 101)=          -3
      flst_born(   2, 101)=          -4
      flst_born(   3, 101)=          25
      flst_born(   4, 101)=          -4
      flst_born(   5, 101)=          -3
 
      flst_born(   1, 102)=          -3
      flst_born(   2, 102)=           4
      flst_born(   3, 102)=          25
      flst_born(   4, 102)=           4
      flst_born(   5, 102)=          -3
 
      flst_born(   1, 103)=          -3
      flst_born(   2, 103)=          -3
      flst_born(   3, 103)=          25
      flst_born(   4, 103)=          -3
      flst_born(   5, 103)=          -3
 
      flst_born(   1, 104)=          -3
      flst_born(   2, 104)=           3
      flst_born(   3, 104)=          25
      flst_born(   4, 104)=           1
      flst_born(   5, 104)=          -1
 
      flst_born(   1, 105)=          -3
      flst_born(   2, 105)=           3
      flst_born(   3, 105)=          25
      flst_born(   4, 105)=           2
      flst_born(   5, 105)=          -2
 
      flst_born(   1, 106)=          -3
      flst_born(   2, 106)=           3
      flst_born(   3, 106)=          25
      flst_born(   4, 106)=           4
      flst_born(   5, 106)=          -4
 
      flst_born(   1, 107)=          -3
      flst_born(   2, 107)=           3
      flst_born(   3, 107)=          25
      flst_born(   4, 107)=           3
      flst_born(   5, 107)=          -3
 
      flst_born(   1, 108)=          -3
      flst_born(   2, 108)=           3
      flst_born(   3, 108)=          25
      flst_born(   4, 108)=           5
      flst_born(   5, 108)=          -5
 
      flst_born(   1, 109)=          -3
      flst_born(   2, 109)=           3
      flst_born(   3, 109)=          25
      flst_born(   4, 109)=           0
      flst_born(   5, 109)=           0
 
      flst_born(   1, 110)=          -3
      flst_born(   2, 110)=          -5
      flst_born(   3, 110)=          25
      flst_born(   4, 110)=          -3
      flst_born(   5, 110)=          -5
 
      flst_born(   1, 111)=          -3
      flst_born(   2, 111)=           5
      flst_born(   3, 111)=          25
      flst_born(   4, 111)=          -3
      flst_born(   5, 111)=           5
 
      flst_born(   1, 112)=          -3
      flst_born(   2, 112)=           0
      flst_born(   3, 112)=          25
      flst_born(   4, 112)=          -3
      flst_born(   5, 112)=           0
 
      flst_born(   1, 113)=           3
      flst_born(   2, 113)=          -1
      flst_born(   3, 113)=          25
      flst_born(   4, 113)=          -1
      flst_born(   5, 113)=           3
 
      flst_born(   1, 114)=           3
      flst_born(   2, 114)=           1
      flst_born(   3, 114)=          25
      flst_born(   4, 114)=           1
      flst_born(   5, 114)=           3
 
      flst_born(   1, 115)=           3
      flst_born(   2, 115)=          -2
      flst_born(   3, 115)=          25
      flst_born(   4, 115)=          -2
      flst_born(   5, 115)=           3
 
      flst_born(   1, 116)=           3
      flst_born(   2, 116)=           2
      flst_born(   3, 116)=          25
      flst_born(   4, 116)=           2
      flst_born(   5, 116)=           3
 
      flst_born(   1, 117)=           3
      flst_born(   2, 117)=          -4
      flst_born(   3, 117)=          25
      flst_born(   4, 117)=          -4
      flst_born(   5, 117)=           3
 
      flst_born(   1, 118)=           3
      flst_born(   2, 118)=           4
      flst_born(   3, 118)=          25
      flst_born(   4, 118)=           4
      flst_born(   5, 118)=           3
 
      flst_born(   1, 119)=           3
      flst_born(   2, 119)=          -3
      flst_born(   3, 119)=          25
      flst_born(   4, 119)=           1
      flst_born(   5, 119)=          -1
 
      flst_born(   1, 120)=           3
      flst_born(   2, 120)=          -3
      flst_born(   3, 120)=          25
      flst_born(   4, 120)=           2
      flst_born(   5, 120)=          -2
 
      flst_born(   1, 121)=           3
      flst_born(   2, 121)=          -3
      flst_born(   3, 121)=          25
      flst_born(   4, 121)=           4
      flst_born(   5, 121)=          -4
 
      flst_born(   1, 122)=           3
      flst_born(   2, 122)=          -3
      flst_born(   3, 122)=          25
      flst_born(   4, 122)=           3
      flst_born(   5, 122)=          -3
 
      flst_born(   1, 123)=           3
      flst_born(   2, 123)=          -3
      flst_born(   3, 123)=          25
      flst_born(   4, 123)=           5
      flst_born(   5, 123)=          -5
 
      flst_born(   1, 124)=           3
      flst_born(   2, 124)=          -3
      flst_born(   3, 124)=          25
      flst_born(   4, 124)=           0
      flst_born(   5, 124)=           0
 
      flst_born(   1, 125)=           3
      flst_born(   2, 125)=           3
      flst_born(   3, 125)=          25
      flst_born(   4, 125)=           3
      flst_born(   5, 125)=           3
 
      flst_born(   1, 126)=           3
      flst_born(   2, 126)=          -5
      flst_born(   3, 126)=          25
      flst_born(   4, 126)=           3
      flst_born(   5, 126)=          -5
 
      flst_born(   1, 127)=           3
      flst_born(   2, 127)=           5
      flst_born(   3, 127)=          25
      flst_born(   4, 127)=           3
      flst_born(   5, 127)=           5
 
      flst_born(   1, 128)=           3
      flst_born(   2, 128)=           0
      flst_born(   3, 128)=          25
      flst_born(   4, 128)=           3
      flst_born(   5, 128)=           0
 
      flst_born(   1, 129)=          -5
      flst_born(   2, 129)=          -1
      flst_born(   3, 129)=          25
      flst_born(   4, 129)=          -1
      flst_born(   5, 129)=          -5
 
      flst_born(   1, 130)=          -5
      flst_born(   2, 130)=           1
      flst_born(   3, 130)=          25
      flst_born(   4, 130)=           1
      flst_born(   5, 130)=          -5
 
      flst_born(   1, 131)=          -5
      flst_born(   2, 131)=          -2
      flst_born(   3, 131)=          25
      flst_born(   4, 131)=          -2
      flst_born(   5, 131)=          -5
 
      flst_born(   1, 132)=          -5
      flst_born(   2, 132)=           2
      flst_born(   3, 132)=          25
      flst_born(   4, 132)=           2
      flst_born(   5, 132)=          -5
 
      flst_born(   1, 133)=          -5
      flst_born(   2, 133)=          -4
      flst_born(   3, 133)=          25
      flst_born(   4, 133)=          -4
      flst_born(   5, 133)=          -5
 
      flst_born(   1, 134)=          -5
      flst_born(   2, 134)=           4
      flst_born(   3, 134)=          25
      flst_born(   4, 134)=           4
      flst_born(   5, 134)=          -5
 
      flst_born(   1, 135)=          -5
      flst_born(   2, 135)=          -3
      flst_born(   3, 135)=          25
      flst_born(   4, 135)=          -3
      flst_born(   5, 135)=          -5
 
      flst_born(   1, 136)=          -5
      flst_born(   2, 136)=           3
      flst_born(   3, 136)=          25
      flst_born(   4, 136)=           3
      flst_born(   5, 136)=          -5
 
      flst_born(   1, 137)=          -5
      flst_born(   2, 137)=          -5
      flst_born(   3, 137)=          25
      flst_born(   4, 137)=          -5
      flst_born(   5, 137)=          -5
 
      flst_born(   1, 138)=          -5
      flst_born(   2, 138)=           5
      flst_born(   3, 138)=          25
      flst_born(   4, 138)=           1
      flst_born(   5, 138)=          -1
 
      flst_born(   1, 139)=          -5
      flst_born(   2, 139)=           5
      flst_born(   3, 139)=          25
      flst_born(   4, 139)=           2
      flst_born(   5, 139)=          -2
 
      flst_born(   1, 140)=          -5
      flst_born(   2, 140)=           5
      flst_born(   3, 140)=          25
      flst_born(   4, 140)=           4
      flst_born(   5, 140)=          -4
 
      flst_born(   1, 141)=          -5
      flst_born(   2, 141)=           5
      flst_born(   3, 141)=          25
      flst_born(   4, 141)=           3
      flst_born(   5, 141)=          -3
 
      flst_born(   1, 142)=          -5
      flst_born(   2, 142)=           5
      flst_born(   3, 142)=          25
      flst_born(   4, 142)=           5
      flst_born(   5, 142)=          -5
 
      flst_born(   1, 143)=          -5
      flst_born(   2, 143)=           5
      flst_born(   3, 143)=          25
      flst_born(   4, 143)=           0
      flst_born(   5, 143)=           0
 
      flst_born(   1, 144)=          -5
      flst_born(   2, 144)=           0
      flst_born(   3, 144)=          25
      flst_born(   4, 144)=          -5
      flst_born(   5, 144)=           0
 
      flst_born(   1, 145)=           5
      flst_born(   2, 145)=          -1
      flst_born(   3, 145)=          25
      flst_born(   4, 145)=          -1
      flst_born(   5, 145)=           5
 
      flst_born(   1, 146)=           5
      flst_born(   2, 146)=           1
      flst_born(   3, 146)=          25
      flst_born(   4, 146)=           1
      flst_born(   5, 146)=           5
 
      flst_born(   1, 147)=           5
      flst_born(   2, 147)=          -2
      flst_born(   3, 147)=          25
      flst_born(   4, 147)=          -2
      flst_born(   5, 147)=           5
 
      flst_born(   1, 148)=           5
      flst_born(   2, 148)=           2
      flst_born(   3, 148)=          25
      flst_born(   4, 148)=           2
      flst_born(   5, 148)=           5
 
      flst_born(   1, 149)=           5
      flst_born(   2, 149)=          -4
      flst_born(   3, 149)=          25
      flst_born(   4, 149)=          -4
      flst_born(   5, 149)=           5
 
      flst_born(   1, 150)=           5
      flst_born(   2, 150)=           4
      flst_born(   3, 150)=          25
      flst_born(   4, 150)=           4
      flst_born(   5, 150)=           5
 
      flst_born(   1, 151)=           5
      flst_born(   2, 151)=          -3
      flst_born(   3, 151)=          25
      flst_born(   4, 151)=          -3
      flst_born(   5, 151)=           5
 
      flst_born(   1, 152)=           5
      flst_born(   2, 152)=           3
      flst_born(   3, 152)=          25
      flst_born(   4, 152)=           3
      flst_born(   5, 152)=           5
 
      flst_born(   1, 153)=           5
      flst_born(   2, 153)=          -5
      flst_born(   3, 153)=          25
      flst_born(   4, 153)=           1
      flst_born(   5, 153)=          -1
 
      flst_born(   1, 154)=           5
      flst_born(   2, 154)=          -5
      flst_born(   3, 154)=          25
      flst_born(   4, 154)=           2
      flst_born(   5, 154)=          -2
 
      flst_born(   1, 155)=           5
      flst_born(   2, 155)=          -5
      flst_born(   3, 155)=          25
      flst_born(   4, 155)=           4
      flst_born(   5, 155)=          -4
 
      flst_born(   1, 156)=           5
      flst_born(   2, 156)=          -5
      flst_born(   3, 156)=          25
      flst_born(   4, 156)=           3
      flst_born(   5, 156)=          -3
 
      flst_born(   1, 157)=           5
      flst_born(   2, 157)=          -5
      flst_born(   3, 157)=          25
      flst_born(   4, 157)=           5
      flst_born(   5, 157)=          -5
 
      flst_born(   1, 158)=           5
      flst_born(   2, 158)=          -5
      flst_born(   3, 158)=          25
      flst_born(   4, 158)=           0
      flst_born(   5, 158)=           0
 
      flst_born(   1, 159)=           5
      flst_born(   2, 159)=           5
      flst_born(   3, 159)=          25
      flst_born(   4, 159)=           5
      flst_born(   5, 159)=           5
 
      flst_born(   1, 160)=           5
      flst_born(   2, 160)=           0
      flst_born(   3, 160)=          25
      flst_born(   4, 160)=           5
      flst_born(   5, 160)=           0
 
      flst_born(   1, 161)=           0
      flst_born(   2, 161)=          -1
      flst_born(   3, 161)=          25
      flst_born(   4, 161)=          -1
      flst_born(   5, 161)=           0
 
      flst_born(   1, 162)=           0
      flst_born(   2, 162)=           1
      flst_born(   3, 162)=          25
      flst_born(   4, 162)=           1
      flst_born(   5, 162)=           0
 
      flst_born(   1, 163)=           0
      flst_born(   2, 163)=          -2
      flst_born(   3, 163)=          25
      flst_born(   4, 163)=          -2
      flst_born(   5, 163)=           0
 
      flst_born(   1, 164)=           0
      flst_born(   2, 164)=           2
      flst_born(   3, 164)=          25
      flst_born(   4, 164)=           2
      flst_born(   5, 164)=           0
 
      flst_born(   1, 165)=           0
      flst_born(   2, 165)=          -4
      flst_born(   3, 165)=          25
      flst_born(   4, 165)=          -4
      flst_born(   5, 165)=           0
 
      flst_born(   1, 166)=           0
      flst_born(   2, 166)=           4
      flst_born(   3, 166)=          25
      flst_born(   4, 166)=           4
      flst_born(   5, 166)=           0
 
      flst_born(   1, 167)=           0
      flst_born(   2, 167)=          -3
      flst_born(   3, 167)=          25
      flst_born(   4, 167)=          -3
      flst_born(   5, 167)=           0
 
      flst_born(   1, 168)=           0
      flst_born(   2, 168)=           3
      flst_born(   3, 168)=          25
      flst_born(   4, 168)=           3
      flst_born(   5, 168)=           0
 
      flst_born(   1, 169)=           0
      flst_born(   2, 169)=          -5
      flst_born(   3, 169)=          25
      flst_born(   4, 169)=          -5
      flst_born(   5, 169)=           0
 
      flst_born(   1, 170)=           0
      flst_born(   2, 170)=           5
      flst_born(   3, 170)=          25
      flst_born(   4, 170)=           5
      flst_born(   5, 170)=           0
 
      flst_born(   1, 171)=           0
      flst_born(   2, 171)=           0
      flst_born(   3, 171)=          25
      flst_born(   4, 171)=           1
      flst_born(   5, 171)=          -1
 
      flst_born(   1, 172)=           0
      flst_born(   2, 172)=           0
      flst_born(   3, 172)=          25
      flst_born(   4, 172)=           2
      flst_born(   5, 172)=          -2
 
      flst_born(   1, 173)=           0
      flst_born(   2, 173)=           0
      flst_born(   3, 173)=          25
      flst_born(   4, 173)=           4
      flst_born(   5, 173)=          -4
 
      flst_born(   1, 174)=           0
      flst_born(   2, 174)=           0
      flst_born(   3, 174)=          25
      flst_born(   4, 174)=           3
      flst_born(   5, 174)=          -3
 
      flst_born(   1, 175)=           0
      flst_born(   2, 175)=           0
      flst_born(   3, 175)=          25
      flst_born(   4, 175)=           5
      flst_born(   5, 175)=          -5
 
      flst_born(   1, 176)=           0
      flst_born(   2, 176)=           0
      flst_born(   3, 176)=          25
      flst_born(   4, 176)=           0
      flst_born(   5, 176)=           0
 
      flst_nborn=         176
 
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
      flst_real(   6,   1)=           0
 
      flst_real(   1,   2)=          -1
      flst_real(   2,   2)=           1
      flst_real(   3,   2)=          25
      flst_real(   4,   2)=           1
      flst_real(   5,   2)=          -1
      flst_real(   6,   2)=           0
 
      flst_real(   1,   3)=          -1
      flst_real(   2,   3)=           1
      flst_real(   3,   3)=          25
      flst_real(   4,   3)=           2
      flst_real(   5,   3)=          -2
      flst_real(   6,   3)=           0
 
      flst_real(   1,   4)=          -1
      flst_real(   2,   4)=           1
      flst_real(   3,   4)=          25
      flst_real(   4,   4)=           4
      flst_real(   5,   4)=          -4
      flst_real(   6,   4)=           0
 
      flst_real(   1,   5)=          -1
      flst_real(   2,   5)=           1
      flst_real(   3,   5)=          25
      flst_real(   4,   5)=           3
      flst_real(   5,   5)=          -3
      flst_real(   6,   5)=           0
 
      flst_real(   1,   6)=          -1
      flst_real(   2,   6)=           1
      flst_real(   3,   6)=          25
      flst_real(   4,   6)=           5
      flst_real(   5,   6)=          -5
      flst_real(   6,   6)=           0
 
      flst_real(   1,   7)=          -1
      flst_real(   2,   7)=           1
      flst_real(   3,   7)=          25
      flst_real(   4,   7)=           0
      flst_real(   5,   7)=           0
      flst_real(   6,   7)=           0
 
      flst_real(   1,   8)=          -1
      flst_real(   2,   8)=          -2
      flst_real(   3,   8)=          25
      flst_real(   4,   8)=          -1
      flst_real(   5,   8)=          -2
      flst_real(   6,   8)=           0
 
      flst_real(   1,   9)=          -1
      flst_real(   2,   9)=           2
      flst_real(   3,   9)=          25
      flst_real(   4,   9)=          -1
      flst_real(   5,   9)=           2
      flst_real(   6,   9)=           0
 
      flst_real(   1,  10)=          -1
      flst_real(   2,  10)=          -4
      flst_real(   3,  10)=          25
      flst_real(   4,  10)=          -1
      flst_real(   5,  10)=          -4
      flst_real(   6,  10)=           0
 
      flst_real(   1,  11)=          -1
      flst_real(   2,  11)=           4
      flst_real(   3,  11)=          25
      flst_real(   4,  11)=          -1
      flst_real(   5,  11)=           4
      flst_real(   6,  11)=           0
 
      flst_real(   1,  12)=          -1
      flst_real(   2,  12)=          -3
      flst_real(   3,  12)=          25
      flst_real(   4,  12)=          -1
      flst_real(   5,  12)=          -3
      flst_real(   6,  12)=           0
 
      flst_real(   1,  13)=          -1
      flst_real(   2,  13)=           3
      flst_real(   3,  13)=          25
      flst_real(   4,  13)=          -1
      flst_real(   5,  13)=           3
      flst_real(   6,  13)=           0
 
      flst_real(   1,  14)=          -1
      flst_real(   2,  14)=          -5
      flst_real(   3,  14)=          25
      flst_real(   4,  14)=          -1
      flst_real(   5,  14)=          -5
      flst_real(   6,  14)=           0
 
      flst_real(   1,  15)=          -1
      flst_real(   2,  15)=           5
      flst_real(   3,  15)=          25
      flst_real(   4,  15)=          -1
      flst_real(   5,  15)=           5
      flst_real(   6,  15)=           0
 
      flst_real(   1,  16)=          -1
      flst_real(   2,  16)=           0
      flst_real(   3,  16)=          25
      flst_real(   4,  16)=           1
      flst_real(   5,  16)=          -1
      flst_real(   6,  16)=          -1
 
      flst_real(   1,  17)=          -1
      flst_real(   2,  17)=           0
      flst_real(   3,  17)=          25
      flst_real(   4,  17)=          -1
      flst_real(   5,  17)=           2
      flst_real(   6,  17)=          -2
 
      flst_real(   1,  18)=          -1
      flst_real(   2,  18)=           0
      flst_real(   3,  18)=          25
      flst_real(   4,  18)=          -1
      flst_real(   5,  18)=           4
      flst_real(   6,  18)=          -4
 
      flst_real(   1,  19)=          -1
      flst_real(   2,  19)=           0
      flst_real(   3,  19)=          25
      flst_real(   4,  19)=          -1
      flst_real(   5,  19)=           3
      flst_real(   6,  19)=          -3
 
      flst_real(   1,  20)=          -1
      flst_real(   2,  20)=           0
      flst_real(   3,  20)=          25
      flst_real(   4,  20)=          -1
      flst_real(   5,  20)=           5
      flst_real(   6,  20)=          -5
 
      flst_real(   1,  21)=          -1
      flst_real(   2,  21)=           0
      flst_real(   3,  21)=          25
      flst_real(   4,  21)=          -1
      flst_real(   5,  21)=           0
      flst_real(   6,  21)=           0
 
      flst_real(   1,  22)=           1
      flst_real(   2,  22)=          -1
      flst_real(   3,  22)=          25
      flst_real(   4,  22)=           1
      flst_real(   5,  22)=          -1
      flst_real(   6,  22)=           0
 
      flst_real(   1,  23)=           1
      flst_real(   2,  23)=          -1
      flst_real(   3,  23)=          25
      flst_real(   4,  23)=           2
      flst_real(   5,  23)=          -2
      flst_real(   6,  23)=           0
 
      flst_real(   1,  24)=           1
      flst_real(   2,  24)=          -1
      flst_real(   3,  24)=          25
      flst_real(   4,  24)=           4
      flst_real(   5,  24)=          -4
      flst_real(   6,  24)=           0
 
      flst_real(   1,  25)=           1
      flst_real(   2,  25)=          -1
      flst_real(   3,  25)=          25
      flst_real(   4,  25)=           3
      flst_real(   5,  25)=          -3
      flst_real(   6,  25)=           0
 
      flst_real(   1,  26)=           1
      flst_real(   2,  26)=          -1
      flst_real(   3,  26)=          25
      flst_real(   4,  26)=           5
      flst_real(   5,  26)=          -5
      flst_real(   6,  26)=           0
 
      flst_real(   1,  27)=           1
      flst_real(   2,  27)=          -1
      flst_real(   3,  27)=          25
      flst_real(   4,  27)=           0
      flst_real(   5,  27)=           0
      flst_real(   6,  27)=           0
 
      flst_real(   1,  28)=           1
      flst_real(   2,  28)=           1
      flst_real(   3,  28)=          25
      flst_real(   4,  28)=           1
      flst_real(   5,  28)=           1
      flst_real(   6,  28)=           0
 
      flst_real(   1,  29)=           1
      flst_real(   2,  29)=          -2
      flst_real(   3,  29)=          25
      flst_real(   4,  29)=           1
      flst_real(   5,  29)=          -2
      flst_real(   6,  29)=           0
 
      flst_real(   1,  30)=           1
      flst_real(   2,  30)=           2
      flst_real(   3,  30)=          25
      flst_real(   4,  30)=           1
      flst_real(   5,  30)=           2
      flst_real(   6,  30)=           0
 
      flst_real(   1,  31)=           1
      flst_real(   2,  31)=          -4
      flst_real(   3,  31)=          25
      flst_real(   4,  31)=           1
      flst_real(   5,  31)=          -4
      flst_real(   6,  31)=           0
 
      flst_real(   1,  32)=           1
      flst_real(   2,  32)=           4
      flst_real(   3,  32)=          25
      flst_real(   4,  32)=           1
      flst_real(   5,  32)=           4
      flst_real(   6,  32)=           0
 
      flst_real(   1,  33)=           1
      flst_real(   2,  33)=          -3
      flst_real(   3,  33)=          25
      flst_real(   4,  33)=           1
      flst_real(   5,  33)=          -3
      flst_real(   6,  33)=           0
 
      flst_real(   1,  34)=           1
      flst_real(   2,  34)=           3
      flst_real(   3,  34)=          25
      flst_real(   4,  34)=           1
      flst_real(   5,  34)=           3
      flst_real(   6,  34)=           0
 
      flst_real(   1,  35)=           1
      flst_real(   2,  35)=          -5
      flst_real(   3,  35)=          25
      flst_real(   4,  35)=           1
      flst_real(   5,  35)=          -5
      flst_real(   6,  35)=           0
 
      flst_real(   1,  36)=           1
      flst_real(   2,  36)=           5
      flst_real(   3,  36)=          25
      flst_real(   4,  36)=           1
      flst_real(   5,  36)=           5
      flst_real(   6,  36)=           0
 
      flst_real(   1,  37)=           1
      flst_real(   2,  37)=           0
      flst_real(   3,  37)=          25
      flst_real(   4,  37)=           1
      flst_real(   5,  37)=           1
      flst_real(   6,  37)=          -1
 
      flst_real(   1,  38)=           1
      flst_real(   2,  38)=           0
      flst_real(   3,  38)=          25
      flst_real(   4,  38)=           1
      flst_real(   5,  38)=           2
      flst_real(   6,  38)=          -2
 
      flst_real(   1,  39)=           1
      flst_real(   2,  39)=           0
      flst_real(   3,  39)=          25
      flst_real(   4,  39)=           1
      flst_real(   5,  39)=           4
      flst_real(   6,  39)=          -4
 
      flst_real(   1,  40)=           1
      flst_real(   2,  40)=           0
      flst_real(   3,  40)=          25
      flst_real(   4,  40)=           1
      flst_real(   5,  40)=           3
      flst_real(   6,  40)=          -3
 
      flst_real(   1,  41)=           1
      flst_real(   2,  41)=           0
      flst_real(   3,  41)=          25
      flst_real(   4,  41)=           1
      flst_real(   5,  41)=           5
      flst_real(   6,  41)=          -5
 
      flst_real(   1,  42)=           1
      flst_real(   2,  42)=           0
      flst_real(   3,  42)=          25
      flst_real(   4,  42)=           1
      flst_real(   5,  42)=           0
      flst_real(   6,  42)=           0
 
      flst_real(   1,  43)=          -2
      flst_real(   2,  43)=          -1
      flst_real(   3,  43)=          25
      flst_real(   4,  43)=          -1
      flst_real(   5,  43)=          -2
      flst_real(   6,  43)=           0
 
      flst_real(   1,  44)=          -2
      flst_real(   2,  44)=           1
      flst_real(   3,  44)=          25
      flst_real(   4,  44)=           1
      flst_real(   5,  44)=          -2
      flst_real(   6,  44)=           0
 
      flst_real(   1,  45)=          -2
      flst_real(   2,  45)=          -2
      flst_real(   3,  45)=          25
      flst_real(   4,  45)=          -2
      flst_real(   5,  45)=          -2
      flst_real(   6,  45)=           0
 
      flst_real(   1,  46)=          -2
      flst_real(   2,  46)=           2
      flst_real(   3,  46)=          25
      flst_real(   4,  46)=           1
      flst_real(   5,  46)=          -1
      flst_real(   6,  46)=           0
 
      flst_real(   1,  47)=          -2
      flst_real(   2,  47)=           2
      flst_real(   3,  47)=          25
      flst_real(   4,  47)=           2
      flst_real(   5,  47)=          -2
      flst_real(   6,  47)=           0
 
      flst_real(   1,  48)=          -2
      flst_real(   2,  48)=           2
      flst_real(   3,  48)=          25
      flst_real(   4,  48)=           4
      flst_real(   5,  48)=          -4
      flst_real(   6,  48)=           0
 
      flst_real(   1,  49)=          -2
      flst_real(   2,  49)=           2
      flst_real(   3,  49)=          25
      flst_real(   4,  49)=           3
      flst_real(   5,  49)=          -3
      flst_real(   6,  49)=           0
 
      flst_real(   1,  50)=          -2
      flst_real(   2,  50)=           2
      flst_real(   3,  50)=          25
      flst_real(   4,  50)=           5
      flst_real(   5,  50)=          -5
      flst_real(   6,  50)=           0
 
      flst_real(   1,  51)=          -2
      flst_real(   2,  51)=           2
      flst_real(   3,  51)=          25
      flst_real(   4,  51)=           0
      flst_real(   5,  51)=           0
      flst_real(   6,  51)=           0
 
      flst_real(   1,  52)=          -2
      flst_real(   2,  52)=          -4
      flst_real(   3,  52)=          25
      flst_real(   4,  52)=          -2
      flst_real(   5,  52)=          -4
      flst_real(   6,  52)=           0
 
      flst_real(   1,  53)=          -2
      flst_real(   2,  53)=           4
      flst_real(   3,  53)=          25
      flst_real(   4,  53)=          -2
      flst_real(   5,  53)=           4
      flst_real(   6,  53)=           0
 
      flst_real(   1,  54)=          -2
      flst_real(   2,  54)=          -3
      flst_real(   3,  54)=          25
      flst_real(   4,  54)=          -2
      flst_real(   5,  54)=          -3
      flst_real(   6,  54)=           0
 
      flst_real(   1,  55)=          -2
      flst_real(   2,  55)=           3
      flst_real(   3,  55)=          25
      flst_real(   4,  55)=          -2
      flst_real(   5,  55)=           3
      flst_real(   6,  55)=           0
 
      flst_real(   1,  56)=          -2
      flst_real(   2,  56)=          -5
      flst_real(   3,  56)=          25
      flst_real(   4,  56)=          -2
      flst_real(   5,  56)=          -5
      flst_real(   6,  56)=           0
 
      flst_real(   1,  57)=          -2
      flst_real(   2,  57)=           5
      flst_real(   3,  57)=          25
      flst_real(   4,  57)=          -2
      flst_real(   5,  57)=           5
      flst_real(   6,  57)=           0
 
      flst_real(   1,  58)=          -2
      flst_real(   2,  58)=           0
      flst_real(   3,  58)=          25
      flst_real(   4,  58)=           1
      flst_real(   5,  58)=          -1
      flst_real(   6,  58)=          -2
 
      flst_real(   1,  59)=          -2
      flst_real(   2,  59)=           0
      flst_real(   3,  59)=          25
      flst_real(   4,  59)=           2
      flst_real(   5,  59)=          -2
      flst_real(   6,  59)=          -2
 
      flst_real(   1,  60)=          -2
      flst_real(   2,  60)=           0
      flst_real(   3,  60)=          25
      flst_real(   4,  60)=          -2
      flst_real(   5,  60)=           4
      flst_real(   6,  60)=          -4
 
      flst_real(   1,  61)=          -2
      flst_real(   2,  61)=           0
      flst_real(   3,  61)=          25
      flst_real(   4,  61)=          -2
      flst_real(   5,  61)=           3
      flst_real(   6,  61)=          -3
 
      flst_real(   1,  62)=          -2
      flst_real(   2,  62)=           0
      flst_real(   3,  62)=          25
      flst_real(   4,  62)=          -2
      flst_real(   5,  62)=           5
      flst_real(   6,  62)=          -5
 
      flst_real(   1,  63)=          -2
      flst_real(   2,  63)=           0
      flst_real(   3,  63)=          25
      flst_real(   4,  63)=          -2
      flst_real(   5,  63)=           0
      flst_real(   6,  63)=           0
 
      flst_real(   1,  64)=           2
      flst_real(   2,  64)=          -1
      flst_real(   3,  64)=          25
      flst_real(   4,  64)=          -1
      flst_real(   5,  64)=           2
      flst_real(   6,  64)=           0
 
      flst_real(   1,  65)=           2
      flst_real(   2,  65)=           1
      flst_real(   3,  65)=          25
      flst_real(   4,  65)=           1
      flst_real(   5,  65)=           2
      flst_real(   6,  65)=           0
 
      flst_real(   1,  66)=           2
      flst_real(   2,  66)=          -2
      flst_real(   3,  66)=          25
      flst_real(   4,  66)=           1
      flst_real(   5,  66)=          -1
      flst_real(   6,  66)=           0
 
      flst_real(   1,  67)=           2
      flst_real(   2,  67)=          -2
      flst_real(   3,  67)=          25
      flst_real(   4,  67)=           2
      flst_real(   5,  67)=          -2
      flst_real(   6,  67)=           0
 
      flst_real(   1,  68)=           2
      flst_real(   2,  68)=          -2
      flst_real(   3,  68)=          25
      flst_real(   4,  68)=           4
      flst_real(   5,  68)=          -4
      flst_real(   6,  68)=           0
 
      flst_real(   1,  69)=           2
      flst_real(   2,  69)=          -2
      flst_real(   3,  69)=          25
      flst_real(   4,  69)=           3
      flst_real(   5,  69)=          -3
      flst_real(   6,  69)=           0
 
      flst_real(   1,  70)=           2
      flst_real(   2,  70)=          -2
      flst_real(   3,  70)=          25
      flst_real(   4,  70)=           5
      flst_real(   5,  70)=          -5
      flst_real(   6,  70)=           0
 
      flst_real(   1,  71)=           2
      flst_real(   2,  71)=          -2
      flst_real(   3,  71)=          25
      flst_real(   4,  71)=           0
      flst_real(   5,  71)=           0
      flst_real(   6,  71)=           0
 
      flst_real(   1,  72)=           2
      flst_real(   2,  72)=           2
      flst_real(   3,  72)=          25
      flst_real(   4,  72)=           2
      flst_real(   5,  72)=           2
      flst_real(   6,  72)=           0
 
      flst_real(   1,  73)=           2
      flst_real(   2,  73)=          -4
      flst_real(   3,  73)=          25
      flst_real(   4,  73)=           2
      flst_real(   5,  73)=          -4
      flst_real(   6,  73)=           0
 
      flst_real(   1,  74)=           2
      flst_real(   2,  74)=           4
      flst_real(   3,  74)=          25
      flst_real(   4,  74)=           2
      flst_real(   5,  74)=           4
      flst_real(   6,  74)=           0
 
      flst_real(   1,  75)=           2
      flst_real(   2,  75)=          -3
      flst_real(   3,  75)=          25
      flst_real(   4,  75)=           2
      flst_real(   5,  75)=          -3
      flst_real(   6,  75)=           0
 
      flst_real(   1,  76)=           2
      flst_real(   2,  76)=           3
      flst_real(   3,  76)=          25
      flst_real(   4,  76)=           2
      flst_real(   5,  76)=           3
      flst_real(   6,  76)=           0
 
      flst_real(   1,  77)=           2
      flst_real(   2,  77)=          -5
      flst_real(   3,  77)=          25
      flst_real(   4,  77)=           2
      flst_real(   5,  77)=          -5
      flst_real(   6,  77)=           0
 
      flst_real(   1,  78)=           2
      flst_real(   2,  78)=           5
      flst_real(   3,  78)=          25
      flst_real(   4,  78)=           2
      flst_real(   5,  78)=           5
      flst_real(   6,  78)=           0
 
      flst_real(   1,  79)=           2
      flst_real(   2,  79)=           0
      flst_real(   3,  79)=          25
      flst_real(   4,  79)=           1
      flst_real(   5,  79)=          -1
      flst_real(   6,  79)=           2
 
      flst_real(   1,  80)=           2
      flst_real(   2,  80)=           0
      flst_real(   3,  80)=          25
      flst_real(   4,  80)=           2
      flst_real(   5,  80)=           2
      flst_real(   6,  80)=          -2
 
      flst_real(   1,  81)=           2
      flst_real(   2,  81)=           0
      flst_real(   3,  81)=          25
      flst_real(   4,  81)=           2
      flst_real(   5,  81)=           4
      flst_real(   6,  81)=          -4
 
      flst_real(   1,  82)=           2
      flst_real(   2,  82)=           0
      flst_real(   3,  82)=          25
      flst_real(   4,  82)=           2
      flst_real(   5,  82)=           3
      flst_real(   6,  82)=          -3
 
      flst_real(   1,  83)=           2
      flst_real(   2,  83)=           0
      flst_real(   3,  83)=          25
      flst_real(   4,  83)=           2
      flst_real(   5,  83)=           5
      flst_real(   6,  83)=          -5
 
      flst_real(   1,  84)=           2
      flst_real(   2,  84)=           0
      flst_real(   3,  84)=          25
      flst_real(   4,  84)=           2
      flst_real(   5,  84)=           0
      flst_real(   6,  84)=           0
 
      flst_real(   1,  85)=          -4
      flst_real(   2,  85)=          -1
      flst_real(   3,  85)=          25
      flst_real(   4,  85)=          -1
      flst_real(   5,  85)=          -4
      flst_real(   6,  85)=           0
 
      flst_real(   1,  86)=          -4
      flst_real(   2,  86)=           1
      flst_real(   3,  86)=          25
      flst_real(   4,  86)=           1
      flst_real(   5,  86)=          -4
      flst_real(   6,  86)=           0
 
      flst_real(   1,  87)=          -4
      flst_real(   2,  87)=          -2
      flst_real(   3,  87)=          25
      flst_real(   4,  87)=          -2
      flst_real(   5,  87)=          -4
      flst_real(   6,  87)=           0
 
      flst_real(   1,  88)=          -4
      flst_real(   2,  88)=           2
      flst_real(   3,  88)=          25
      flst_real(   4,  88)=           2
      flst_real(   5,  88)=          -4
      flst_real(   6,  88)=           0
 
      flst_real(   1,  89)=          -4
      flst_real(   2,  89)=          -4
      flst_real(   3,  89)=          25
      flst_real(   4,  89)=          -4
      flst_real(   5,  89)=          -4
      flst_real(   6,  89)=           0
 
      flst_real(   1,  90)=          -4
      flst_real(   2,  90)=           4
      flst_real(   3,  90)=          25
      flst_real(   4,  90)=           1
      flst_real(   5,  90)=          -1
      flst_real(   6,  90)=           0
 
      flst_real(   1,  91)=          -4
      flst_real(   2,  91)=           4
      flst_real(   3,  91)=          25
      flst_real(   4,  91)=           2
      flst_real(   5,  91)=          -2
      flst_real(   6,  91)=           0
 
      flst_real(   1,  92)=          -4
      flst_real(   2,  92)=           4
      flst_real(   3,  92)=          25
      flst_real(   4,  92)=           4
      flst_real(   5,  92)=          -4
      flst_real(   6,  92)=           0
 
      flst_real(   1,  93)=          -4
      flst_real(   2,  93)=           4
      flst_real(   3,  93)=          25
      flst_real(   4,  93)=           3
      flst_real(   5,  93)=          -3
      flst_real(   6,  93)=           0
 
      flst_real(   1,  94)=          -4
      flst_real(   2,  94)=           4
      flst_real(   3,  94)=          25
      flst_real(   4,  94)=           5
      flst_real(   5,  94)=          -5
      flst_real(   6,  94)=           0
 
      flst_real(   1,  95)=          -4
      flst_real(   2,  95)=           4
      flst_real(   3,  95)=          25
      flst_real(   4,  95)=           0
      flst_real(   5,  95)=           0
      flst_real(   6,  95)=           0
 
      flst_real(   1,  96)=          -4
      flst_real(   2,  96)=          -3
      flst_real(   3,  96)=          25
      flst_real(   4,  96)=          -4
      flst_real(   5,  96)=          -3
      flst_real(   6,  96)=           0
 
      flst_real(   1,  97)=          -4
      flst_real(   2,  97)=           3
      flst_real(   3,  97)=          25
      flst_real(   4,  97)=          -4
      flst_real(   5,  97)=           3
      flst_real(   6,  97)=           0
 
      flst_real(   1,  98)=          -4
      flst_real(   2,  98)=          -5
      flst_real(   3,  98)=          25
      flst_real(   4,  98)=          -4
      flst_real(   5,  98)=          -5
      flst_real(   6,  98)=           0
 
      flst_real(   1,  99)=          -4
      flst_real(   2,  99)=           5
      flst_real(   3,  99)=          25
      flst_real(   4,  99)=          -4
      flst_real(   5,  99)=           5
      flst_real(   6,  99)=           0
 
      flst_real(   1, 100)=          -4
      flst_real(   2, 100)=           0
      flst_real(   3, 100)=          25
      flst_real(   4, 100)=           1
      flst_real(   5, 100)=          -1
      flst_real(   6, 100)=          -4
 
      flst_real(   1, 101)=          -4
      flst_real(   2, 101)=           0
      flst_real(   3, 101)=          25
      flst_real(   4, 101)=           2
      flst_real(   5, 101)=          -2
      flst_real(   6, 101)=          -4
 
      flst_real(   1, 102)=          -4
      flst_real(   2, 102)=           0
      flst_real(   3, 102)=          25
      flst_real(   4, 102)=           4
      flst_real(   5, 102)=          -4
      flst_real(   6, 102)=          -4
 
      flst_real(   1, 103)=          -4
      flst_real(   2, 103)=           0
      flst_real(   3, 103)=          25
      flst_real(   4, 103)=          -4
      flst_real(   5, 103)=           3
      flst_real(   6, 103)=          -3
 
      flst_real(   1, 104)=          -4
      flst_real(   2, 104)=           0
      flst_real(   3, 104)=          25
      flst_real(   4, 104)=          -4
      flst_real(   5, 104)=           5
      flst_real(   6, 104)=          -5
 
      flst_real(   1, 105)=          -4
      flst_real(   2, 105)=           0
      flst_real(   3, 105)=          25
      flst_real(   4, 105)=          -4
      flst_real(   5, 105)=           0
      flst_real(   6, 105)=           0
 
      flst_real(   1, 106)=           4
      flst_real(   2, 106)=          -1
      flst_real(   3, 106)=          25
      flst_real(   4, 106)=          -1
      flst_real(   5, 106)=           4
      flst_real(   6, 106)=           0
 
      flst_real(   1, 107)=           4
      flst_real(   2, 107)=           1
      flst_real(   3, 107)=          25
      flst_real(   4, 107)=           1
      flst_real(   5, 107)=           4
      flst_real(   6, 107)=           0
 
      flst_real(   1, 108)=           4
      flst_real(   2, 108)=          -2
      flst_real(   3, 108)=          25
      flst_real(   4, 108)=          -2
      flst_real(   5, 108)=           4
      flst_real(   6, 108)=           0
 
      flst_real(   1, 109)=           4
      flst_real(   2, 109)=           2
      flst_real(   3, 109)=          25
      flst_real(   4, 109)=           2
      flst_real(   5, 109)=           4
      flst_real(   6, 109)=           0
 
      flst_real(   1, 110)=           4
      flst_real(   2, 110)=          -4
      flst_real(   3, 110)=          25
      flst_real(   4, 110)=           1
      flst_real(   5, 110)=          -1
      flst_real(   6, 110)=           0
 
      flst_real(   1, 111)=           4
      flst_real(   2, 111)=          -4
      flst_real(   3, 111)=          25
      flst_real(   4, 111)=           2
      flst_real(   5, 111)=          -2
      flst_real(   6, 111)=           0
 
      flst_real(   1, 112)=           4
      flst_real(   2, 112)=          -4
      flst_real(   3, 112)=          25
      flst_real(   4, 112)=           4
      flst_real(   5, 112)=          -4
      flst_real(   6, 112)=           0
 
      flst_real(   1, 113)=           4
      flst_real(   2, 113)=          -4
      flst_real(   3, 113)=          25
      flst_real(   4, 113)=           3
      flst_real(   5, 113)=          -3
      flst_real(   6, 113)=           0
 
      flst_real(   1, 114)=           4
      flst_real(   2, 114)=          -4
      flst_real(   3, 114)=          25
      flst_real(   4, 114)=           5
      flst_real(   5, 114)=          -5
      flst_real(   6, 114)=           0
 
      flst_real(   1, 115)=           4
      flst_real(   2, 115)=          -4
      flst_real(   3, 115)=          25
      flst_real(   4, 115)=           0
      flst_real(   5, 115)=           0
      flst_real(   6, 115)=           0
 
      flst_real(   1, 116)=           4
      flst_real(   2, 116)=           4
      flst_real(   3, 116)=          25
      flst_real(   4, 116)=           4
      flst_real(   5, 116)=           4
      flst_real(   6, 116)=           0
 
      flst_real(   1, 117)=           4
      flst_real(   2, 117)=          -3
      flst_real(   3, 117)=          25
      flst_real(   4, 117)=           4
      flst_real(   5, 117)=          -3
      flst_real(   6, 117)=           0
 
      flst_real(   1, 118)=           4
      flst_real(   2, 118)=           3
      flst_real(   3, 118)=          25
      flst_real(   4, 118)=           4
      flst_real(   5, 118)=           3
      flst_real(   6, 118)=           0
 
      flst_real(   1, 119)=           4
      flst_real(   2, 119)=          -5
      flst_real(   3, 119)=          25
      flst_real(   4, 119)=           4
      flst_real(   5, 119)=          -5
      flst_real(   6, 119)=           0
 
      flst_real(   1, 120)=           4
      flst_real(   2, 120)=           5
      flst_real(   3, 120)=          25
      flst_real(   4, 120)=           4
      flst_real(   5, 120)=           5
      flst_real(   6, 120)=           0
 
      flst_real(   1, 121)=           4
      flst_real(   2, 121)=           0
      flst_real(   3, 121)=          25
      flst_real(   4, 121)=           1
      flst_real(   5, 121)=          -1
      flst_real(   6, 121)=           4
 
      flst_real(   1, 122)=           4
      flst_real(   2, 122)=           0
      flst_real(   3, 122)=          25
      flst_real(   4, 122)=           2
      flst_real(   5, 122)=          -2
      flst_real(   6, 122)=           4
 
      flst_real(   1, 123)=           4
      flst_real(   2, 123)=           0
      flst_real(   3, 123)=          25
      flst_real(   4, 123)=           4
      flst_real(   5, 123)=           4
      flst_real(   6, 123)=          -4
 
      flst_real(   1, 124)=           4
      flst_real(   2, 124)=           0
      flst_real(   3, 124)=          25
      flst_real(   4, 124)=           4
      flst_real(   5, 124)=           3
      flst_real(   6, 124)=          -3
 
      flst_real(   1, 125)=           4
      flst_real(   2, 125)=           0
      flst_real(   3, 125)=          25
      flst_real(   4, 125)=           4
      flst_real(   5, 125)=           5
      flst_real(   6, 125)=          -5
 
      flst_real(   1, 126)=           4
      flst_real(   2, 126)=           0
      flst_real(   3, 126)=          25
      flst_real(   4, 126)=           4
      flst_real(   5, 126)=           0
      flst_real(   6, 126)=           0
 
      flst_real(   1, 127)=          -3
      flst_real(   2, 127)=          -1
      flst_real(   3, 127)=          25
      flst_real(   4, 127)=          -1
      flst_real(   5, 127)=          -3
      flst_real(   6, 127)=           0
 
      flst_real(   1, 128)=          -3
      flst_real(   2, 128)=           1
      flst_real(   3, 128)=          25
      flst_real(   4, 128)=           1
      flst_real(   5, 128)=          -3
      flst_real(   6, 128)=           0
 
      flst_real(   1, 129)=          -3
      flst_real(   2, 129)=          -2
      flst_real(   3, 129)=          25
      flst_real(   4, 129)=          -2
      flst_real(   5, 129)=          -3
      flst_real(   6, 129)=           0
 
      flst_real(   1, 130)=          -3
      flst_real(   2, 130)=           2
      flst_real(   3, 130)=          25
      flst_real(   4, 130)=           2
      flst_real(   5, 130)=          -3
      flst_real(   6, 130)=           0
 
      flst_real(   1, 131)=          -3
      flst_real(   2, 131)=          -4
      flst_real(   3, 131)=          25
      flst_real(   4, 131)=          -4
      flst_real(   5, 131)=          -3
      flst_real(   6, 131)=           0
 
      flst_real(   1, 132)=          -3
      flst_real(   2, 132)=           4
      flst_real(   3, 132)=          25
      flst_real(   4, 132)=           4
      flst_real(   5, 132)=          -3
      flst_real(   6, 132)=           0
 
      flst_real(   1, 133)=          -3
      flst_real(   2, 133)=          -3
      flst_real(   3, 133)=          25
      flst_real(   4, 133)=          -3
      flst_real(   5, 133)=          -3
      flst_real(   6, 133)=           0
 
      flst_real(   1, 134)=          -3
      flst_real(   2, 134)=           3
      flst_real(   3, 134)=          25
      flst_real(   4, 134)=           1
      flst_real(   5, 134)=          -1
      flst_real(   6, 134)=           0
 
      flst_real(   1, 135)=          -3
      flst_real(   2, 135)=           3
      flst_real(   3, 135)=          25
      flst_real(   4, 135)=           2
      flst_real(   5, 135)=          -2
      flst_real(   6, 135)=           0
 
      flst_real(   1, 136)=          -3
      flst_real(   2, 136)=           3
      flst_real(   3, 136)=          25
      flst_real(   4, 136)=           4
      flst_real(   5, 136)=          -4
      flst_real(   6, 136)=           0
 
      flst_real(   1, 137)=          -3
      flst_real(   2, 137)=           3
      flst_real(   3, 137)=          25
      flst_real(   4, 137)=           3
      flst_real(   5, 137)=          -3
      flst_real(   6, 137)=           0
 
      flst_real(   1, 138)=          -3
      flst_real(   2, 138)=           3
      flst_real(   3, 138)=          25
      flst_real(   4, 138)=           5
      flst_real(   5, 138)=          -5
      flst_real(   6, 138)=           0
 
      flst_real(   1, 139)=          -3
      flst_real(   2, 139)=           3
      flst_real(   3, 139)=          25
      flst_real(   4, 139)=           0
      flst_real(   5, 139)=           0
      flst_real(   6, 139)=           0
 
      flst_real(   1, 140)=          -3
      flst_real(   2, 140)=          -5
      flst_real(   3, 140)=          25
      flst_real(   4, 140)=          -3
      flst_real(   5, 140)=          -5
      flst_real(   6, 140)=           0
 
      flst_real(   1, 141)=          -3
      flst_real(   2, 141)=           5
      flst_real(   3, 141)=          25
      flst_real(   4, 141)=          -3
      flst_real(   5, 141)=           5
      flst_real(   6, 141)=           0
 
      flst_real(   1, 142)=          -3
      flst_real(   2, 142)=           0
      flst_real(   3, 142)=          25
      flst_real(   4, 142)=           1
      flst_real(   5, 142)=          -1
      flst_real(   6, 142)=          -3
 
      flst_real(   1, 143)=          -3
      flst_real(   2, 143)=           0
      flst_real(   3, 143)=          25
      flst_real(   4, 143)=           2
      flst_real(   5, 143)=          -2
      flst_real(   6, 143)=          -3
 
      flst_real(   1, 144)=          -3
      flst_real(   2, 144)=           0
      flst_real(   3, 144)=          25
      flst_real(   4, 144)=           4
      flst_real(   5, 144)=          -4
      flst_real(   6, 144)=          -3
 
      flst_real(   1, 145)=          -3
      flst_real(   2, 145)=           0
      flst_real(   3, 145)=          25
      flst_real(   4, 145)=           3
      flst_real(   5, 145)=          -3
      flst_real(   6, 145)=          -3
 
      flst_real(   1, 146)=          -3
      flst_real(   2, 146)=           0
      flst_real(   3, 146)=          25
      flst_real(   4, 146)=          -3
      flst_real(   5, 146)=           5
      flst_real(   6, 146)=          -5
 
      flst_real(   1, 147)=          -3
      flst_real(   2, 147)=           0
      flst_real(   3, 147)=          25
      flst_real(   4, 147)=          -3
      flst_real(   5, 147)=           0
      flst_real(   6, 147)=           0
 
      flst_real(   1, 148)=           3
      flst_real(   2, 148)=          -1
      flst_real(   3, 148)=          25
      flst_real(   4, 148)=          -1
      flst_real(   5, 148)=           3
      flst_real(   6, 148)=           0
 
      flst_real(   1, 149)=           3
      flst_real(   2, 149)=           1
      flst_real(   3, 149)=          25
      flst_real(   4, 149)=           1
      flst_real(   5, 149)=           3
      flst_real(   6, 149)=           0
 
      flst_real(   1, 150)=           3
      flst_real(   2, 150)=          -2
      flst_real(   3, 150)=          25
      flst_real(   4, 150)=          -2
      flst_real(   5, 150)=           3
      flst_real(   6, 150)=           0
 
      flst_real(   1, 151)=           3
      flst_real(   2, 151)=           2
      flst_real(   3, 151)=          25
      flst_real(   4, 151)=           2
      flst_real(   5, 151)=           3
      flst_real(   6, 151)=           0
 
      flst_real(   1, 152)=           3
      flst_real(   2, 152)=          -4
      flst_real(   3, 152)=          25
      flst_real(   4, 152)=          -4
      flst_real(   5, 152)=           3
      flst_real(   6, 152)=           0
 
      flst_real(   1, 153)=           3
      flst_real(   2, 153)=           4
      flst_real(   3, 153)=          25
      flst_real(   4, 153)=           4
      flst_real(   5, 153)=           3
      flst_real(   6, 153)=           0
 
      flst_real(   1, 154)=           3
      flst_real(   2, 154)=          -3
      flst_real(   3, 154)=          25
      flst_real(   4, 154)=           1
      flst_real(   5, 154)=          -1
      flst_real(   6, 154)=           0
 
      flst_real(   1, 155)=           3
      flst_real(   2, 155)=          -3
      flst_real(   3, 155)=          25
      flst_real(   4, 155)=           2
      flst_real(   5, 155)=          -2
      flst_real(   6, 155)=           0
 
      flst_real(   1, 156)=           3
      flst_real(   2, 156)=          -3
      flst_real(   3, 156)=          25
      flst_real(   4, 156)=           4
      flst_real(   5, 156)=          -4
      flst_real(   6, 156)=           0
 
      flst_real(   1, 157)=           3
      flst_real(   2, 157)=          -3
      flst_real(   3, 157)=          25
      flst_real(   4, 157)=           3
      flst_real(   5, 157)=          -3
      flst_real(   6, 157)=           0
 
      flst_real(   1, 158)=           3
      flst_real(   2, 158)=          -3
      flst_real(   3, 158)=          25
      flst_real(   4, 158)=           5
      flst_real(   5, 158)=          -5
      flst_real(   6, 158)=           0
 
      flst_real(   1, 159)=           3
      flst_real(   2, 159)=          -3
      flst_real(   3, 159)=          25
      flst_real(   4, 159)=           0
      flst_real(   5, 159)=           0
      flst_real(   6, 159)=           0
 
      flst_real(   1, 160)=           3
      flst_real(   2, 160)=           3
      flst_real(   3, 160)=          25
      flst_real(   4, 160)=           3
      flst_real(   5, 160)=           3
      flst_real(   6, 160)=           0
 
      flst_real(   1, 161)=           3
      flst_real(   2, 161)=          -5
      flst_real(   3, 161)=          25
      flst_real(   4, 161)=           3
      flst_real(   5, 161)=          -5
      flst_real(   6, 161)=           0
 
      flst_real(   1, 162)=           3
      flst_real(   2, 162)=           5
      flst_real(   3, 162)=          25
      flst_real(   4, 162)=           3
      flst_real(   5, 162)=           5
      flst_real(   6, 162)=           0
 
      flst_real(   1, 163)=           3
      flst_real(   2, 163)=           0
      flst_real(   3, 163)=          25
      flst_real(   4, 163)=           1
      flst_real(   5, 163)=          -1
      flst_real(   6, 163)=           3
 
      flst_real(   1, 164)=           3
      flst_real(   2, 164)=           0
      flst_real(   3, 164)=          25
      flst_real(   4, 164)=           2
      flst_real(   5, 164)=          -2
      flst_real(   6, 164)=           3
 
      flst_real(   1, 165)=           3
      flst_real(   2, 165)=           0
      flst_real(   3, 165)=          25
      flst_real(   4, 165)=           4
      flst_real(   5, 165)=          -4
      flst_real(   6, 165)=           3
 
      flst_real(   1, 166)=           3
      flst_real(   2, 166)=           0
      flst_real(   3, 166)=          25
      flst_real(   4, 166)=           3
      flst_real(   5, 166)=           3
      flst_real(   6, 166)=          -3
 
      flst_real(   1, 167)=           3
      flst_real(   2, 167)=           0
      flst_real(   3, 167)=          25
      flst_real(   4, 167)=           3
      flst_real(   5, 167)=           5
      flst_real(   6, 167)=          -5
 
      flst_real(   1, 168)=           3
      flst_real(   2, 168)=           0
      flst_real(   3, 168)=          25
      flst_real(   4, 168)=           3
      flst_real(   5, 168)=           0
      flst_real(   6, 168)=           0
 
      flst_real(   1, 169)=          -5
      flst_real(   2, 169)=          -1
      flst_real(   3, 169)=          25
      flst_real(   4, 169)=          -1
      flst_real(   5, 169)=          -5
      flst_real(   6, 169)=           0
 
      flst_real(   1, 170)=          -5
      flst_real(   2, 170)=           1
      flst_real(   3, 170)=          25
      flst_real(   4, 170)=           1
      flst_real(   5, 170)=          -5
      flst_real(   6, 170)=           0
 
      flst_real(   1, 171)=          -5
      flst_real(   2, 171)=          -2
      flst_real(   3, 171)=          25
      flst_real(   4, 171)=          -2
      flst_real(   5, 171)=          -5
      flst_real(   6, 171)=           0
 
      flst_real(   1, 172)=          -5
      flst_real(   2, 172)=           2
      flst_real(   3, 172)=          25
      flst_real(   4, 172)=           2
      flst_real(   5, 172)=          -5
      flst_real(   6, 172)=           0
 
      flst_real(   1, 173)=          -5
      flst_real(   2, 173)=          -4
      flst_real(   3, 173)=          25
      flst_real(   4, 173)=          -4
      flst_real(   5, 173)=          -5
      flst_real(   6, 173)=           0
 
      flst_real(   1, 174)=          -5
      flst_real(   2, 174)=           4
      flst_real(   3, 174)=          25
      flst_real(   4, 174)=           4
      flst_real(   5, 174)=          -5
      flst_real(   6, 174)=           0
 
      flst_real(   1, 175)=          -5
      flst_real(   2, 175)=          -3
      flst_real(   3, 175)=          25
      flst_real(   4, 175)=          -3
      flst_real(   5, 175)=          -5
      flst_real(   6, 175)=           0
 
      flst_real(   1, 176)=          -5
      flst_real(   2, 176)=           3
      flst_real(   3, 176)=          25
      flst_real(   4, 176)=           3
      flst_real(   5, 176)=          -5
      flst_real(   6, 176)=           0
 
      flst_real(   1, 177)=          -5
      flst_real(   2, 177)=          -5
      flst_real(   3, 177)=          25
      flst_real(   4, 177)=          -5
      flst_real(   5, 177)=          -5
      flst_real(   6, 177)=           0
 
      flst_real(   1, 178)=          -5
      flst_real(   2, 178)=           5
      flst_real(   3, 178)=          25
      flst_real(   4, 178)=           1
      flst_real(   5, 178)=          -1
      flst_real(   6, 178)=           0
 
      flst_real(   1, 179)=          -5
      flst_real(   2, 179)=           5
      flst_real(   3, 179)=          25
      flst_real(   4, 179)=           2
      flst_real(   5, 179)=          -2
      flst_real(   6, 179)=           0
 
      flst_real(   1, 180)=          -5
      flst_real(   2, 180)=           5
      flst_real(   3, 180)=          25
      flst_real(   4, 180)=           4
      flst_real(   5, 180)=          -4
      flst_real(   6, 180)=           0
 
      flst_real(   1, 181)=          -5
      flst_real(   2, 181)=           5
      flst_real(   3, 181)=          25
      flst_real(   4, 181)=           3
      flst_real(   5, 181)=          -3
      flst_real(   6, 181)=           0
 
      flst_real(   1, 182)=          -5
      flst_real(   2, 182)=           5
      flst_real(   3, 182)=          25
      flst_real(   4, 182)=           5
      flst_real(   5, 182)=          -5
      flst_real(   6, 182)=           0
 
      flst_real(   1, 183)=          -5
      flst_real(   2, 183)=           5
      flst_real(   3, 183)=          25
      flst_real(   4, 183)=           0
      flst_real(   5, 183)=           0
      flst_real(   6, 183)=           0
 
      flst_real(   1, 184)=          -5
      flst_real(   2, 184)=           0
      flst_real(   3, 184)=          25
      flst_real(   4, 184)=           1
      flst_real(   5, 184)=          -1
      flst_real(   6, 184)=          -5
 
      flst_real(   1, 185)=          -5
      flst_real(   2, 185)=           0
      flst_real(   3, 185)=          25
      flst_real(   4, 185)=           2
      flst_real(   5, 185)=          -2
      flst_real(   6, 185)=          -5
 
      flst_real(   1, 186)=          -5
      flst_real(   2, 186)=           0
      flst_real(   3, 186)=          25
      flst_real(   4, 186)=           4
      flst_real(   5, 186)=          -4
      flst_real(   6, 186)=          -5
 
      flst_real(   1, 187)=          -5
      flst_real(   2, 187)=           0
      flst_real(   3, 187)=          25
      flst_real(   4, 187)=           3
      flst_real(   5, 187)=          -3
      flst_real(   6, 187)=          -5
 
      flst_real(   1, 188)=          -5
      flst_real(   2, 188)=           0
      flst_real(   3, 188)=          25
      flst_real(   4, 188)=           5
      flst_real(   5, 188)=          -5
      flst_real(   6, 188)=          -5
 
      flst_real(   1, 189)=          -5
      flst_real(   2, 189)=           0
      flst_real(   3, 189)=          25
      flst_real(   4, 189)=          -5
      flst_real(   5, 189)=           0
      flst_real(   6, 189)=           0
 
      flst_real(   1, 190)=           5
      flst_real(   2, 190)=          -1
      flst_real(   3, 190)=          25
      flst_real(   4, 190)=          -1
      flst_real(   5, 190)=           5
      flst_real(   6, 190)=           0
 
      flst_real(   1, 191)=           5
      flst_real(   2, 191)=           1
      flst_real(   3, 191)=          25
      flst_real(   4, 191)=           1
      flst_real(   5, 191)=           5
      flst_real(   6, 191)=           0
 
      flst_real(   1, 192)=           5
      flst_real(   2, 192)=          -2
      flst_real(   3, 192)=          25
      flst_real(   4, 192)=          -2
      flst_real(   5, 192)=           5
      flst_real(   6, 192)=           0
 
      flst_real(   1, 193)=           5
      flst_real(   2, 193)=           2
      flst_real(   3, 193)=          25
      flst_real(   4, 193)=           2
      flst_real(   5, 193)=           5
      flst_real(   6, 193)=           0
 
      flst_real(   1, 194)=           5
      flst_real(   2, 194)=          -4
      flst_real(   3, 194)=          25
      flst_real(   4, 194)=          -4
      flst_real(   5, 194)=           5
      flst_real(   6, 194)=           0
 
      flst_real(   1, 195)=           5
      flst_real(   2, 195)=           4
      flst_real(   3, 195)=          25
      flst_real(   4, 195)=           4
      flst_real(   5, 195)=           5
      flst_real(   6, 195)=           0
 
      flst_real(   1, 196)=           5
      flst_real(   2, 196)=          -3
      flst_real(   3, 196)=          25
      flst_real(   4, 196)=          -3
      flst_real(   5, 196)=           5
      flst_real(   6, 196)=           0
 
      flst_real(   1, 197)=           5
      flst_real(   2, 197)=           3
      flst_real(   3, 197)=          25
      flst_real(   4, 197)=           3
      flst_real(   5, 197)=           5
      flst_real(   6, 197)=           0
 
      flst_real(   1, 198)=           5
      flst_real(   2, 198)=          -5
      flst_real(   3, 198)=          25
      flst_real(   4, 198)=           1
      flst_real(   5, 198)=          -1
      flst_real(   6, 198)=           0
 
      flst_real(   1, 199)=           5
      flst_real(   2, 199)=          -5
      flst_real(   3, 199)=          25
      flst_real(   4, 199)=           2
      flst_real(   5, 199)=          -2
      flst_real(   6, 199)=           0
 
      flst_real(   1, 200)=           5
      flst_real(   2, 200)=          -5
      flst_real(   3, 200)=          25
      flst_real(   4, 200)=           4
      flst_real(   5, 200)=          -4
      flst_real(   6, 200)=           0
 
      flst_real(   1, 201)=           5
      flst_real(   2, 201)=          -5
      flst_real(   3, 201)=          25
      flst_real(   4, 201)=           3
      flst_real(   5, 201)=          -3
      flst_real(   6, 201)=           0
 
      flst_real(   1, 202)=           5
      flst_real(   2, 202)=          -5
      flst_real(   3, 202)=          25
      flst_real(   4, 202)=           5
      flst_real(   5, 202)=          -5
      flst_real(   6, 202)=           0
 
      flst_real(   1, 203)=           5
      flst_real(   2, 203)=          -5
      flst_real(   3, 203)=          25
      flst_real(   4, 203)=           0
      flst_real(   5, 203)=           0
      flst_real(   6, 203)=           0
 
      flst_real(   1, 204)=           5
      flst_real(   2, 204)=           5
      flst_real(   3, 204)=          25
      flst_real(   4, 204)=           5
      flst_real(   5, 204)=           5
      flst_real(   6, 204)=           0
 
      flst_real(   1, 205)=           5
      flst_real(   2, 205)=           0
      flst_real(   3, 205)=          25
      flst_real(   4, 205)=           1
      flst_real(   5, 205)=          -1
      flst_real(   6, 205)=           5
 
      flst_real(   1, 206)=           5
      flst_real(   2, 206)=           0
      flst_real(   3, 206)=          25
      flst_real(   4, 206)=           2
      flst_real(   5, 206)=          -2
      flst_real(   6, 206)=           5
 
      flst_real(   1, 207)=           5
      flst_real(   2, 207)=           0
      flst_real(   3, 207)=          25
      flst_real(   4, 207)=           4
      flst_real(   5, 207)=          -4
      flst_real(   6, 207)=           5
 
      flst_real(   1, 208)=           5
      flst_real(   2, 208)=           0
      flst_real(   3, 208)=          25
      flst_real(   4, 208)=           3
      flst_real(   5, 208)=          -3
      flst_real(   6, 208)=           5
 
      flst_real(   1, 209)=           5
      flst_real(   2, 209)=           0
      flst_real(   3, 209)=          25
      flst_real(   4, 209)=           5
      flst_real(   5, 209)=           5
      flst_real(   6, 209)=          -5
 
      flst_real(   1, 210)=           5
      flst_real(   2, 210)=           0
      flst_real(   3, 210)=          25
      flst_real(   4, 210)=           5
      flst_real(   5, 210)=           0
      flst_real(   6, 210)=           0
 
      flst_real(   1, 211)=           0
      flst_real(   2, 211)=          -1
      flst_real(   3, 211)=          25
      flst_real(   4, 211)=           1
      flst_real(   5, 211)=          -1
      flst_real(   6, 211)=          -1
 
      flst_real(   1, 212)=           0
      flst_real(   2, 212)=          -1
      flst_real(   3, 212)=          25
      flst_real(   4, 212)=          -1
      flst_real(   5, 212)=           2
      flst_real(   6, 212)=          -2
 
      flst_real(   1, 213)=           0
      flst_real(   2, 213)=          -1
      flst_real(   3, 213)=          25
      flst_real(   4, 213)=          -1
      flst_real(   5, 213)=           4
      flst_real(   6, 213)=          -4
 
      flst_real(   1, 214)=           0
      flst_real(   2, 214)=          -1
      flst_real(   3, 214)=          25
      flst_real(   4, 214)=          -1
      flst_real(   5, 214)=           3
      flst_real(   6, 214)=          -3
 
      flst_real(   1, 215)=           0
      flst_real(   2, 215)=          -1
      flst_real(   3, 215)=          25
      flst_real(   4, 215)=          -1
      flst_real(   5, 215)=           5
      flst_real(   6, 215)=          -5
 
      flst_real(   1, 216)=           0
      flst_real(   2, 216)=          -1
      flst_real(   3, 216)=          25
      flst_real(   4, 216)=          -1
      flst_real(   5, 216)=           0
      flst_real(   6, 216)=           0
 
      flst_real(   1, 217)=           0
      flst_real(   2, 217)=           1
      flst_real(   3, 217)=          25
      flst_real(   4, 217)=           1
      flst_real(   5, 217)=           1
      flst_real(   6, 217)=          -1
 
      flst_real(   1, 218)=           0
      flst_real(   2, 218)=           1
      flst_real(   3, 218)=          25
      flst_real(   4, 218)=           1
      flst_real(   5, 218)=           2
      flst_real(   6, 218)=          -2
 
      flst_real(   1, 219)=           0
      flst_real(   2, 219)=           1
      flst_real(   3, 219)=          25
      flst_real(   4, 219)=           1
      flst_real(   5, 219)=           4
      flst_real(   6, 219)=          -4
 
      flst_real(   1, 220)=           0
      flst_real(   2, 220)=           1
      flst_real(   3, 220)=          25
      flst_real(   4, 220)=           1
      flst_real(   5, 220)=           3
      flst_real(   6, 220)=          -3
 
      flst_real(   1, 221)=           0
      flst_real(   2, 221)=           1
      flst_real(   3, 221)=          25
      flst_real(   4, 221)=           1
      flst_real(   5, 221)=           5
      flst_real(   6, 221)=          -5
 
      flst_real(   1, 222)=           0
      flst_real(   2, 222)=           1
      flst_real(   3, 222)=          25
      flst_real(   4, 222)=           1
      flst_real(   5, 222)=           0
      flst_real(   6, 222)=           0
 
      flst_real(   1, 223)=           0
      flst_real(   2, 223)=          -2
      flst_real(   3, 223)=          25
      flst_real(   4, 223)=           1
      flst_real(   5, 223)=          -1
      flst_real(   6, 223)=          -2
 
      flst_real(   1, 224)=           0
      flst_real(   2, 224)=          -2
      flst_real(   3, 224)=          25
      flst_real(   4, 224)=           2
      flst_real(   5, 224)=          -2
      flst_real(   6, 224)=          -2
 
      flst_real(   1, 225)=           0
      flst_real(   2, 225)=          -2
      flst_real(   3, 225)=          25
      flst_real(   4, 225)=          -2
      flst_real(   5, 225)=           4
      flst_real(   6, 225)=          -4
 
      flst_real(   1, 226)=           0
      flst_real(   2, 226)=          -2
      flst_real(   3, 226)=          25
      flst_real(   4, 226)=          -2
      flst_real(   5, 226)=           3
      flst_real(   6, 226)=          -3
 
      flst_real(   1, 227)=           0
      flst_real(   2, 227)=          -2
      flst_real(   3, 227)=          25
      flst_real(   4, 227)=          -2
      flst_real(   5, 227)=           5
      flst_real(   6, 227)=          -5
 
      flst_real(   1, 228)=           0
      flst_real(   2, 228)=          -2
      flst_real(   3, 228)=          25
      flst_real(   4, 228)=          -2
      flst_real(   5, 228)=           0
      flst_real(   6, 228)=           0
 
      flst_real(   1, 229)=           0
      flst_real(   2, 229)=           2
      flst_real(   3, 229)=          25
      flst_real(   4, 229)=           1
      flst_real(   5, 229)=          -1
      flst_real(   6, 229)=           2
 
      flst_real(   1, 230)=           0
      flst_real(   2, 230)=           2
      flst_real(   3, 230)=          25
      flst_real(   4, 230)=           2
      flst_real(   5, 230)=           2
      flst_real(   6, 230)=          -2
 
      flst_real(   1, 231)=           0
      flst_real(   2, 231)=           2
      flst_real(   3, 231)=          25
      flst_real(   4, 231)=           2
      flst_real(   5, 231)=           4
      flst_real(   6, 231)=          -4
 
      flst_real(   1, 232)=           0
      flst_real(   2, 232)=           2
      flst_real(   3, 232)=          25
      flst_real(   4, 232)=           2
      flst_real(   5, 232)=           3
      flst_real(   6, 232)=          -3
 
      flst_real(   1, 233)=           0
      flst_real(   2, 233)=           2
      flst_real(   3, 233)=          25
      flst_real(   4, 233)=           2
      flst_real(   5, 233)=           5
      flst_real(   6, 233)=          -5
 
      flst_real(   1, 234)=           0
      flst_real(   2, 234)=           2
      flst_real(   3, 234)=          25
      flst_real(   4, 234)=           2
      flst_real(   5, 234)=           0
      flst_real(   6, 234)=           0
 
      flst_real(   1, 235)=           0
      flst_real(   2, 235)=          -4
      flst_real(   3, 235)=          25
      flst_real(   4, 235)=           1
      flst_real(   5, 235)=          -1
      flst_real(   6, 235)=          -4
 
      flst_real(   1, 236)=           0
      flst_real(   2, 236)=          -4
      flst_real(   3, 236)=          25
      flst_real(   4, 236)=           2
      flst_real(   5, 236)=          -2
      flst_real(   6, 236)=          -4
 
      flst_real(   1, 237)=           0
      flst_real(   2, 237)=          -4
      flst_real(   3, 237)=          25
      flst_real(   4, 237)=           4
      flst_real(   5, 237)=          -4
      flst_real(   6, 237)=          -4
 
      flst_real(   1, 238)=           0
      flst_real(   2, 238)=          -4
      flst_real(   3, 238)=          25
      flst_real(   4, 238)=          -4
      flst_real(   5, 238)=           3
      flst_real(   6, 238)=          -3
 
      flst_real(   1, 239)=           0
      flst_real(   2, 239)=          -4
      flst_real(   3, 239)=          25
      flst_real(   4, 239)=          -4
      flst_real(   5, 239)=           5
      flst_real(   6, 239)=          -5
 
      flst_real(   1, 240)=           0
      flst_real(   2, 240)=          -4
      flst_real(   3, 240)=          25
      flst_real(   4, 240)=          -4
      flst_real(   5, 240)=           0
      flst_real(   6, 240)=           0
 
      flst_real(   1, 241)=           0
      flst_real(   2, 241)=           4
      flst_real(   3, 241)=          25
      flst_real(   4, 241)=           1
      flst_real(   5, 241)=          -1
      flst_real(   6, 241)=           4
 
      flst_real(   1, 242)=           0
      flst_real(   2, 242)=           4
      flst_real(   3, 242)=          25
      flst_real(   4, 242)=           2
      flst_real(   5, 242)=          -2
      flst_real(   6, 242)=           4
 
      flst_real(   1, 243)=           0
      flst_real(   2, 243)=           4
      flst_real(   3, 243)=          25
      flst_real(   4, 243)=           4
      flst_real(   5, 243)=           4
      flst_real(   6, 243)=          -4
 
      flst_real(   1, 244)=           0
      flst_real(   2, 244)=           4
      flst_real(   3, 244)=          25
      flst_real(   4, 244)=           4
      flst_real(   5, 244)=           3
      flst_real(   6, 244)=          -3
 
      flst_real(   1, 245)=           0
      flst_real(   2, 245)=           4
      flst_real(   3, 245)=          25
      flst_real(   4, 245)=           4
      flst_real(   5, 245)=           5
      flst_real(   6, 245)=          -5
 
      flst_real(   1, 246)=           0
      flst_real(   2, 246)=           4
      flst_real(   3, 246)=          25
      flst_real(   4, 246)=           4
      flst_real(   5, 246)=           0
      flst_real(   6, 246)=           0
 
      flst_real(   1, 247)=           0
      flst_real(   2, 247)=          -3
      flst_real(   3, 247)=          25
      flst_real(   4, 247)=           1
      flst_real(   5, 247)=          -1
      flst_real(   6, 247)=          -3
 
      flst_real(   1, 248)=           0
      flst_real(   2, 248)=          -3
      flst_real(   3, 248)=          25
      flst_real(   4, 248)=           2
      flst_real(   5, 248)=          -2
      flst_real(   6, 248)=          -3
 
      flst_real(   1, 249)=           0
      flst_real(   2, 249)=          -3
      flst_real(   3, 249)=          25
      flst_real(   4, 249)=           4
      flst_real(   5, 249)=          -4
      flst_real(   6, 249)=          -3
 
      flst_real(   1, 250)=           0
      flst_real(   2, 250)=          -3
      flst_real(   3, 250)=          25
      flst_real(   4, 250)=           3
      flst_real(   5, 250)=          -3
      flst_real(   6, 250)=          -3
 
      flst_real(   1, 251)=           0
      flst_real(   2, 251)=          -3
      flst_real(   3, 251)=          25
      flst_real(   4, 251)=          -3
      flst_real(   5, 251)=           5
      flst_real(   6, 251)=          -5
 
      flst_real(   1, 252)=           0
      flst_real(   2, 252)=          -3
      flst_real(   3, 252)=          25
      flst_real(   4, 252)=          -3
      flst_real(   5, 252)=           0
      flst_real(   6, 252)=           0
 
      flst_real(   1, 253)=           0
      flst_real(   2, 253)=           3
      flst_real(   3, 253)=          25
      flst_real(   4, 253)=           1
      flst_real(   5, 253)=          -1
      flst_real(   6, 253)=           3
 
      flst_real(   1, 254)=           0
      flst_real(   2, 254)=           3
      flst_real(   3, 254)=          25
      flst_real(   4, 254)=           2
      flst_real(   5, 254)=          -2
      flst_real(   6, 254)=           3
 
      flst_real(   1, 255)=           0
      flst_real(   2, 255)=           3
      flst_real(   3, 255)=          25
      flst_real(   4, 255)=           4
      flst_real(   5, 255)=          -4
      flst_real(   6, 255)=           3
 
      flst_real(   1, 256)=           0
      flst_real(   2, 256)=           3
      flst_real(   3, 256)=          25
      flst_real(   4, 256)=           3
      flst_real(   5, 256)=           3
      flst_real(   6, 256)=          -3
 
      flst_real(   1, 257)=           0
      flst_real(   2, 257)=           3
      flst_real(   3, 257)=          25
      flst_real(   4, 257)=           3
      flst_real(   5, 257)=           5
      flst_real(   6, 257)=          -5
 
      flst_real(   1, 258)=           0
      flst_real(   2, 258)=           3
      flst_real(   3, 258)=          25
      flst_real(   4, 258)=           3
      flst_real(   5, 258)=           0
      flst_real(   6, 258)=           0
 
      flst_real(   1, 259)=           0
      flst_real(   2, 259)=          -5
      flst_real(   3, 259)=          25
      flst_real(   4, 259)=           1
      flst_real(   5, 259)=          -1
      flst_real(   6, 259)=          -5
 
      flst_real(   1, 260)=           0
      flst_real(   2, 260)=          -5
      flst_real(   3, 260)=          25
      flst_real(   4, 260)=           2
      flst_real(   5, 260)=          -2
      flst_real(   6, 260)=          -5
 
      flst_real(   1, 261)=           0
      flst_real(   2, 261)=          -5
      flst_real(   3, 261)=          25
      flst_real(   4, 261)=           4
      flst_real(   5, 261)=          -4
      flst_real(   6, 261)=          -5
 
      flst_real(   1, 262)=           0
      flst_real(   2, 262)=          -5
      flst_real(   3, 262)=          25
      flst_real(   4, 262)=           3
      flst_real(   5, 262)=          -3
      flst_real(   6, 262)=          -5
 
      flst_real(   1, 263)=           0
      flst_real(   2, 263)=          -5
      flst_real(   3, 263)=          25
      flst_real(   4, 263)=           5
      flst_real(   5, 263)=          -5
      flst_real(   6, 263)=          -5
 
      flst_real(   1, 264)=           0
      flst_real(   2, 264)=          -5
      flst_real(   3, 264)=          25
      flst_real(   4, 264)=          -5
      flst_real(   5, 264)=           0
      flst_real(   6, 264)=           0
 
      flst_real(   1, 265)=           0
      flst_real(   2, 265)=           5
      flst_real(   3, 265)=          25
      flst_real(   4, 265)=           1
      flst_real(   5, 265)=          -1
      flst_real(   6, 265)=           5
 
      flst_real(   1, 266)=           0
      flst_real(   2, 266)=           5
      flst_real(   3, 266)=          25
      flst_real(   4, 266)=           2
      flst_real(   5, 266)=          -2
      flst_real(   6, 266)=           5
 
      flst_real(   1, 267)=           0
      flst_real(   2, 267)=           5
      flst_real(   3, 267)=          25
      flst_real(   4, 267)=           4
      flst_real(   5, 267)=          -4
      flst_real(   6, 267)=           5
 
      flst_real(   1, 268)=           0
      flst_real(   2, 268)=           5
      flst_real(   3, 268)=          25
      flst_real(   4, 268)=           3
      flst_real(   5, 268)=          -3
      flst_real(   6, 268)=           5
 
      flst_real(   1, 269)=           0
      flst_real(   2, 269)=           5
      flst_real(   3, 269)=          25
      flst_real(   4, 269)=           5
      flst_real(   5, 269)=           5
      flst_real(   6, 269)=          -5
 
      flst_real(   1, 270)=           0
      flst_real(   2, 270)=           5
      flst_real(   3, 270)=          25
      flst_real(   4, 270)=           5
      flst_real(   5, 270)=           0
      flst_real(   6, 270)=           0
 
      flst_real(   1, 271)=           0
      flst_real(   2, 271)=           0
      flst_real(   3, 271)=          25
      flst_real(   4, 271)=           1
      flst_real(   5, 271)=          -1
      flst_real(   6, 271)=           0
 
      flst_real(   1, 272)=           0
      flst_real(   2, 272)=           0
      flst_real(   3, 272)=          25
      flst_real(   4, 272)=           2
      flst_real(   5, 272)=          -2
      flst_real(   6, 272)=           0
 
      flst_real(   1, 273)=           0
      flst_real(   2, 273)=           0
      flst_real(   3, 273)=          25
      flst_real(   4, 273)=           4
      flst_real(   5, 273)=          -4
      flst_real(   6, 273)=           0
 
      flst_real(   1, 274)=           0
      flst_real(   2, 274)=           0
      flst_real(   3, 274)=          25
      flst_real(   4, 274)=           3
      flst_real(   5, 274)=          -3
      flst_real(   6, 274)=           0
 
      flst_real(   1, 275)=           0
      flst_real(   2, 275)=           0
      flst_real(   3, 275)=          25
      flst_real(   4, 275)=           5
      flst_real(   5, 275)=          -5
      flst_real(   6, 275)=           0
 
      flst_real(   1, 276)=           0
      flst_real(   2, 276)=           0
      flst_real(   3, 276)=          25
      flst_real(   4, 276)=           0
      flst_real(   5, 276)=           0
      flst_real(   6, 276)=           0
 
      flst_nreal=         276
 
      return
      end
 
