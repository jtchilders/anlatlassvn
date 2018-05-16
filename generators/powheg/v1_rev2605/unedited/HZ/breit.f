      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,wt)       
      implicit none
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
c
c     msq is distributed with a probability proportional to
c     d msq/((rwidth*rmass)^2+(msq-rmass**2)**2)), The weight is
c     proportional to ((rwidth*rmass)^2+(msq-rmass**2)**2)),
c     with a coefficient (almax-almin)/(rmass*rwidth)
c
      double precision x1,mminsq,mmaxsq,rmass,rwidth,msq,wt
      double precision almin,almax,al,tanal
      include 'pwhg_math.h'
      logical zerowidth 

      zerowidth = .false. 
c--- in case the maximum msq is very small, just generate linearly for safety
      if (mmaxsq .lt. rmass*1d-3) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

      if (zerowidth) then
          tanal=0d0
          almax=+pi/2d0 
          almin=-pi/2d0 
      else
          almin=atan((mminsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*x1+almin
          tanal=tan(al)
      endif

      msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
      wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
      return
      end
