C     returns 2 Re(M_B * M_V)/(as/(2pi)), 
C     where M_B is the Born amplitude and 
C     M_V is the finite parte of the virtual amplitude.
C     The as/(2pi) factor is attached at a later point.
      subroutine setvirtual(p,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'PhysPars.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      real *8 p(0:3,nlegborn),virtual
      integer vflav(nlegborn)
      real * 8 c0,q2,llog
      real *8 born,bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn)
      real *8 dotp
      external dotp

      virtual = 0d0      
      q2=2*dotp(p(0,1),p(0,2))
      
      call setborn(p,vflav,born,bornjk,bmunu)

c     2*Re(finite coefficient of the triangle V_mu)
      llog = log(st_muren2/q2)
      c0 = (pi**2-8-3*llog-llog**2)

c     finite part of the virtual diagram
      virtual=born*c0*CF
      end
