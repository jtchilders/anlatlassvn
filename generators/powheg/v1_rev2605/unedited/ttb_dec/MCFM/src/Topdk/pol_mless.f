      subroutine pol_mless(p,i,f)
c-----massless vector polarization  subroutine
      implicit none
      integer i,pol
      double complex p(4)
      double precision p0,px,py,pz
      double precision pv,ct,st,cphi,sphi
      double complex f(4)

      p0=dreal(p(1))
      px=dreal(p(2))
      py=dreal(p(3))
      pz=dreal(p(4))

      pv=dsqrt(dabs(p0**2))
      ct=pz/pv
      st=dsqrt(dabs(1d0-ct**2))

      if (st.lt.1d-8) then 
         cphi=1d0
         sphi=0d0
      else
         cphi= px/pv/st
         sphi= py/pv/st
      endif

c      the following ifstatement distinguishes between 
c      positive and negative energies
      if ( p0.gt.0d0) then  
      pol=i
      else
      pol=-i
      endif

      f(1)=dcmplx(0d0,0d0)
      f(2)=dcmplx(ct*cphi/dsqrt(2d0),+pol*sphi/dsqrt(2d0))
      f(3)=dcmplx(ct*sphi/dsqrt(2d0),-pol*cphi/dsqrt(2d0))
      f(4)=dcmplx(-st/dsqrt(2d0),0d0)

      return
      end

