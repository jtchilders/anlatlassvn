      subroutine setreal(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 p(0:3,nlegs)
      integer rflav(nlegs)
      real * 8 amp2
      integer i,j
      
      amp2 = 0d0
      i=rflav(1)
      j=rflav(2)
      if((i.eq.j).and.(i.eq.0)) then
c     g g -> h g
         call M2_gg_hg(p,amp2)
      elseif(((i*j.eq.0).and.(i.ne.j)).and.(i.ne.0)) then
c     q g -> h q
         call M2_qg_hq(p,amp2)
      elseif(((i*j.eq.0).and.(i.ne.j)).and.(j.ne.0)) then
c     g q -> h q
         call M2_gq_hq(p,amp2)
      elseif((i+j.eq.0).and.(i.ne.j)) then
c     q aq -> h g
         call M2_qaq_hg(p,amp2)
      else
         write(*,*) 'ERROR setreal: unammissible flavour!'
         write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4)
         call exit(1)
      endif
      
      if(amp2.eq.0d0) then
         write(*,*) 'WARNING setreal: returning 0 amplitude!'
         write(*,*) rflav(1),' ',rflav(2),'->',rflav(3),' ',rflav(4)
      endif

c     cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
      amp2 = amp2/(st_alpha/(2*pi))
      end
      
      
     
      subroutine M2_gg_hg(pphy,amp2)
c     Real matrix element times normalizations and averages.
c     IMPORTANT the flux factor 1/2s is intentionally missing
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'Flags.h'
      include 'pwhg_st.h'
      include 'pwhg_br.h'
      real * 8 s,v2,t,u,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      real* 8 amp2
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 pphy(0:3,nlegs)
      real * 8 m2
      real * 8 agg2
      real * 8 upt1,upt2
      real * 8 tpt1,tpt2
      real * 8 ewcorr
      common/ew/ewcorr
      
      s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
     #-(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2
      
c     CAVEAT to avoid rounding errors causing t to become 0d0.
c      t=(pphy(0,1)-pphy(0,3))**2-(pphy(1,1)-pphy(1,3))**2
c     #-(pphy(2,1)-pphy(2,3))**2 -(pphy(3,1)-pphy(3,3))**2
      upt1 = (pphy(0,1)-pphy(0,3))**2-(pphy(3,1)-pphy(3,3))**2
      upt2 = -(pphy(1,1)-pphy(1,3))**2-(pphy(2,1)-pphy(2,3))**2
      u = upt1+upt2
      
c     CAVEAT to avoid rounding errors causing u/t to become 0d0.
c      u=(pphy(0,1)-pphy(0,4))**2-(pphy(1,1)-pphy(1,4))**2
c     #-(pphy(2,1)-pphy(2,4))**2 -(pphy(3,1)-pphy(3,4))**2
      tpt1 = (pphy(0,1)-pphy(0,4))**2-(pphy(3,1)-pphy(3,4))**2
      tpt2 = -(pphy(2,1)-pphy(2,4))**2 -(pphy(1,1)-pphy(1,4))**2
      t = tpt1+tpt2

      m2=(pphy(0,3))**2-(pphy(1,3))**2
     #-(pphy(2,3))**2 -(pphy(3,3))**2
      
      xnorm=32d0/(256.d0*3.d0*pi*v2)
      call agg(s,t,u,agg2)
      tmp=(agg2/(s*t*u))*9*(m2**4)*4


      amp2=xnorm*tmp*(st_alpha**3)*ph_GF
      if (flg_ew.eq.1) then
         amp2 = amp2*(1d0+ewcorr)
      endif

      end
      
      subroutine M2_qaq_hg(pphy,amp2)
c Real matrix element times normalizations and averages.
c     IMPORTANT the flux factor 1/2s is missing
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      real * 8 s,v2,t,u,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      real* 8 amp2
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 pphy(0:3,nlegs)
      real * 8 aqqbar2
      real * 8 m2
      real * 8 tpt1,tpt2
      real * 8 upt1,upt2
      
       s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
     #    -(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2
       
c     CAVEAT to avoid rounding errors causing t to become 0d0.
c      t=(pphy(0,1)-pphy(0,3))**2-(pphy(1,1)-pphy(1,3))**2
c     #-(pphy(2,1)-pphy(2,3))**2 -(pphy(3,1)-pphy(3,3))**2
      upt1 = (pphy(0,1)-pphy(0,3))**2-(pphy(3,1)-pphy(3,3))**2
      upt2 = -(pphy(1,1)-pphy(1,3))**2-(pphy(2,1)-pphy(2,3))**2
      u = upt1+upt2

c     CAVEAT to avoid rounding errors causing u/t to become 0d0.
c      u=(pphy(0,1)-pphy(0,4))**2-(pphy(1,1)-pphy(1,4))**2
c     #-(pphy(2,1)-pphy(2,4))**2 -(pphy(3,1)-pphy(3,4))**2
      tpt1 = (pphy(0,1)-pphy(0,4))**2-(pphy(3,1)-pphy(3,4))**2
      tpt2 = -(pphy(2,1)-pphy(2,4))**2 -(pphy(1,1)-pphy(1,4))**2
      t = tpt1+tpt2

       m2=(pphy(0,3))**2-(pphy(1,3))**2
     #   -(pphy(2,3))**2 -(pphy(3,3))**2
       
c From eq.(3.1) of NPB359(91)283

      xnorm=16/(36.d0*9.d0*pi*v2)
      call aqqbar(s,aqqbar2)
      tmp = (aqqbar2*9*(m2**2)*4*(u**2+t**2))/(4*s*((u+t)**2))
      amp2=xnorm*tmp*(st_alpha**3)*ph_GF
      end
      
      
      subroutine M2_qg_hq(pphy,amp2)
c Real matrix element times normalizations and averages.
c IMPORTANT the flux factor 1/2s is missing
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'PhysPars.h'
      include 'Flags.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 s,v2,t,u,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      real* 8 amp2
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 pphy(0:3,nlegs)
      real * 8 m2
      real * 8 aqqbar2
      real * 8 tpt1,tpt2
      real * 8 upt1,upt2
      real * 8 ewcorr
      common/ew/ewcorr

      
      s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
     #-(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2

c     CAVEAT to avoid rounding errors causing t to become 0d0.
c      t=(pphy(0,1)-pphy(0,3))**2-(pphy(1,1)-pphy(1,3))**2
c     #-(pphy(2,1)-pphy(2,3))**2 -(pphy(3,1)-pphy(3,3))**2
      upt1 = (pphy(0,1)-pphy(0,3))**2-(pphy(3,1)-pphy(3,3))**2
      upt2 = -(pphy(1,1)-pphy(1,3))**2-(pphy(2,1)-pphy(2,3))**2
      u = upt1+upt2

c     CAVEAT to avoid rounding errors causing t to become 0d0.
c     u=(pphy(0,1)-pphy(0,4))**2-(pphy(1,1)-pphy(1,4))**2
c     #-(pphy(2,1)-pphy(2,4))**2 -(pphy(3,1)-pphy(3,4))**2
      tpt1 = (pphy(0,1)-pphy(0,4))**2-(pphy(3,1)-pphy(3,4))**2
      tpt2 = -(pphy(2,1)-pphy(2,4))**2 -(pphy(1,1)-pphy(1,4))**2
      t = tpt1+tpt2

      m2=(pphy(0,3))**2-(pphy(1,3))**2
     #-(pphy(2,3))**2 -(pphy(3,3))**2

      xnorm=-16/(96.d0*9.d0*pi*v2)
      call aqqbar(t,aqqbar2)
      tmp = (aqqbar2*9*(m2**2)*4*(u**2+s**2))/(4*t*((u+s)**2))
      amp2=xnorm*tmp*(st_alpha**3)*ph_GF
      if (flg_ew.eq.1) then
         amp2 = amp2*(1d0+ewcorr)
      endif
      end



      subroutine M2_gq_hq(pphy,amp2)
c Real matrix element times normalizations and averages.
c IMPORTANT the flux factor 1/2s is missing
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'Flags.h'
      include 'PhysPars.h' 
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 s,v2,t,u,tmp,xnorm
      parameter (v2=0.70710678118654757d0)
      real* 8 amp2
      integer nlegs
      parameter (nlegs=nlegreal)
      real * 8 pphy(0:3,nlegs)
      real * 8 m2
      real * 8 aqqbar2
      real * 8 tpt1,tpt2
      real * 8 upt1,upt2
      real * 8 ewcorr
      common/ew/ewcorr

       s=(pphy(0,1)+pphy(0,2))**2-(pphy(1,1)+pphy(1,2))**2
     #    -(pphy(2,1)+pphy(2,2))**2 -(pphy(3,1)+pphy(3,2))**2

c     CAVEAT to avoid rounding errors causing t to become 0d0.
c      t=(pphy(0,1)-pphy(0,3))**2-(pphy(1,1)-pphy(1,3))**2
c     #-(pphy(2,1)-pphy(2,3))**2 -(pphy(3,1)-pphy(3,3))**2
      upt1 = (pphy(0,1)-pphy(0,3))**2-(pphy(3,1)-pphy(3,3))**2
      upt2 = -(pphy(1,1)-pphy(1,3))**2-(pphy(2,1)-pphy(2,3))**2
      u = upt1+upt2

c     CAVEAT to avoid rounding errors causing u/t to become 0d0.
c     u=(pphy(0,1)-pphy(0,4))**2-(pphy(1,1)-pphy(1,4))**2
c     #-(pphy(2,1)-pphy(2,4))**2 -(pphy(3,1)-pphy(3,4))**2
      tpt1 = (pphy(0,1)-pphy(0,4))**2-(pphy(3,1)-pphy(3,4))**2
      tpt2 = -(pphy(2,1)-pphy(2,4))**2 -(pphy(1,1)-pphy(1,4))**2
      t = tpt1+tpt2

       m2=(pphy(0,3))**2-(pphy(1,3))**2
     #   -(pphy(2,3))**2 -(pphy(3,3))**2

      xnorm=-16/(96.d0*9.d0*pi*v2)
      call aqqbar(u,aqqbar2)
      tmp = (aqqbar2*9*(m2**2)*4*(t**2+s**2))/(4*u*((t+s)**2))
      amp2=xnorm*tmp*(st_alpha**3)*ph_GF
      if (flg_ew.eq.1) then
         amp2 = amp2*(1d0+ewcorr)
      endif
      end
