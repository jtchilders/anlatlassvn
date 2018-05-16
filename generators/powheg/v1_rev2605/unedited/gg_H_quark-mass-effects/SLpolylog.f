***********************************************************************
*
*     MISCELLANEOUS ROUTINES FOR POLYLOGARITHMS
*
*
*
********** real dilogarithm used in the two-loop SUSY contribs  *******
*

      double precision function p_Li2(x) 

      implicit none

      double complex CLI2,z
      double precision x

      z = DCMPLX(x,0d0)
      p_Li2 = DBLE(CLI2(z))     ! call the (slower) complex version

      return
      end
