      subroutine setvirtual(pin,vflav,virtual)
      implicit none
      include 'nlegborn.h'

C     MCFM include files
      include 'constants.f'
      include 'mcfmtopwhg.f'

      integer nlegs
      parameter (nlegs=nlegborn)
      integer ro,j,vflav(nlegborn),myin(nlegs)
      double precision pin(0:3,nlegborn),virtual,mcfmp(mxpart,4)
      real *8 powheginput 
      external powheginput 
      integer, save :: fakevirt 
      logical, save :: firsttime = .true. 
c set mcfm strong coupling and scales
      call st_mcfm
c ---

      if (firsttime) then 
         fakevirt=powheginput("#fakevirt")
         if (fakevirt == 1) write(*,*) 'WARNING: Using fakevirt !'         
         firsttime = .false. 
      endif

C     vector for redirection of powheg vector onto mcfm
         myin(:)=-1
         myin(1)=1
         myin(2)=2
         myin(3)=3
         myin(4)=4
         myin(5)=5
         myin(6)=6


*  MCFM notation
*     q(-p1) +qbar(-p2) -> e-(p3)+e+(p4)+j(p5)+j(p6)

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

C----return value in MSBR scheme divided by ason2pi
      if (fakevirt == 1) then 
      call qqb_z2jet_pwhg(mcfmp,vflav,virtual) ! this is born 
      virtual = virtual*0.5d0 
      else
      call qqb_z2jet_v_pwhg(mcfmp,vflav,virtual)
      endif
       
      end
