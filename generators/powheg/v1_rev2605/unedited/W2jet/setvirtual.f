      subroutine setvirtual(pin,vflav,virtual)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'PhysPars.h'

      include 'constants.f'
      include 'flags.f'
      include 'mcfmtopwhg.f'

      integer nlegs,ro
      parameter (nlegs=nlegborn)
      double precision pin(0:3,nlegborn)
      integer vflav(nlegborn),myin(nlegs)
      double precision virtual
      double precision bornjk(nlegborn,nlegborn)
      double precision bmunu(0:3,0:3,nlegborn),born,tmp
      integer j,i 
      double precision mcfmp(mxpart,4),born2
      integer nq 
      double precision powheginput
c set mcfm strong coupling and scales
      call st_mcfm
c ---
C     vector for redirection of powheg vector onto mcfm
         myin(:)=-1
         myin(1)=1
         myin(2)=2
         myin(3)=3
         myin(4)=4
         myin(5)=5
         myin(6)=6

*  MCFM notation
*     q(-p1) +qbar(-p2) -> nu(p3)+e+(p4)+j(p5)+j(p6)
*     q(-p1) +qbar(-p2) -> e-(p3)+nu~(p4)+j(p5)+j(p6)

      do j=1,nlegs
      do ro=1,4
      if (myin(j) .gt. 0) then
      if (j .le. 2) then
c incoming partons
      mcfmp(myin(j),ro)=-pin(pwhg(ro),j)
      else
      mcfmp(myin(j),ro)=+pin(pwhg(ro),j)
      endif
      endif
      enddo
      enddo

      call qqb_w2jet_v_pwhg(mcfmp,vflav,virtual)

       

      end
