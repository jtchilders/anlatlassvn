c     returns 2 Re(M_B * M_V)/(as/(2pi)), 
c     where M_B is the Born amplitude and 
c     M_V is the finite part of the virtual amplitude
c     The as/(2pi) factor is attached at a later point
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'Flags.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)      
      real * 8 virtual,dummyjk(nlegborn,nlegborn)
      real * 8 born,dummymunu(0:3,0:3,nlegborn)
      real * 8 s,dotp
      external dotp
      real * 8 g2l,prevg2l
      save prevg2l
      data prevg2l /-1d0/
      real * 8 gluongamma
      real * 8 renormcorr
      real * 8 bornqcd
      complex * 16 ampl1loop
      common /bornampl/ampl1loop,bornqcd

      s=2d0*dotp(p(0,1),p(0,2))

      call compborn(p,vflav,born,dummyjk,dummymunu)

c     Trick to avoid calculating every time the virtual correction factor if the Higgs boson is to be produced on-shell.
c$$$      if (flg_onshell.eq.1) then
c$$$         if (prevg2l.ne.-1d0) then
c$$$            g2l = prevg2l
c$$$         else
c$$$            call twoloopfactors(g2l)
c$$$            prevg2l = g2l
c$$$         endif
c$$$      else
c$$$         call twoloopfactors(g2l)
c$$$      endif

      call twoloopfactors(g2l)

      virtual=(2d0*g2l + pi*pi*CA - CA*(log(st_muren2/s))**2)*born
c      write(*,*) 'g2l', g2l
c      stop
      end

