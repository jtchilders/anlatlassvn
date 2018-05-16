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

      double complex CLI2,z, dilog
      double precision x

      z = DCMPLX(x,0d0)
      p_Li2 = DBLE(Dilog(z))     ! call the (slower) complex version

      return
      end

*
********** complex dilogarithm lifted from FeynHiggs ******************
*
* Dilog.F
* complex dilogarithm
* this file is part of FeynHiggs
* last modified 20 Oct 05 th

      double complex function Dilog(z)
      implicit none
      double complex z
      
      double complex Dilogsum
      external Dilogsum
      
      double precision absz, abs1z
      double complex t, mlogz
      
      double precision pi, zeta2
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (zeta2 = pi*pi/6D0)
      
      absz = abs(z)
      if( absz .lt. 1D-20 ) then
         Dilog = -log(1 - z)
         return
      endif
      
      abs1z = abs(1 - z)
	if( abs1z .lt. 1D-20 ) then
           Dilog = zeta2
           return
	endif
        
	if( DBLE(z) .gt. .5D0 ) then
           mlogz = -log(z)
           t = zeta2 + mlogz*log(1 - z)
           if( abs1z .gt. 1 ) then
              Dilog = Dilogsum(log(1 - 1/z)) + zeta2 +
     $             .5D0*log(z - 1)**2 + t
           else
	    Dilog = -Dilogsum(mlogz) + t
         endif
      else
         if( absz .gt. 1 ) then
	    Dilog = -Dilogsum(-log(1 - 1/z)) - zeta2 - .5D0*log(-z)**2
         else
	    Dilog = Dilogsum(-log(1 - z))
         endif
      endif
      end
      
      double complex function Dilogsum(w)
      implicit none
      double complex w
      
      double complex u, t
      integer k
      
      double precision b2, b4, b6, b8, b10, b12, b14
      double precision b16, b18, b20, b22
      parameter (b2 = 1/6D0)
      parameter (b4 = -1/30D0)
      parameter (b6 = 1/42D0)
      parameter (b8 = -1/30D0)
      parameter (b10 = 5/66D0)
      parameter (b12 = -691/2730D0)
      parameter (b14 = 7/6D0)
      parameter (b16 = -3617/510D0)
      parameter (b18 = 43867/798D0)
      parameter (b20 = -174611/330D0)
      parameter (b22 = 854513/138D0)
      
      double precision bernoulliB(11)
      data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
     &     b16, b18, b20, b22/
      
      Dilogsum = w*(1 - .25D0*w)
      if( abs(w) .lt. 1D-10 ) return
      
      u = w
      do k = 1, 11
         u = u*w**2/DBLE(2*k*(2*k + 1))
         t = u*bernoulliB(k)
         Dilogsum = Dilogsum + t
         if( abs(t) .lt. 1D-16*abs(Dilogsum) ) return
      enddo

      end

