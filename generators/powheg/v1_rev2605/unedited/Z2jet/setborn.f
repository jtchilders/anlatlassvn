      subroutine setborn(p,bflav,bres,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'constants.f'
      include 'mcfmtopwhg.f'
      include 'Bmcfm.f'
      integer nlegs
      parameter (nlegs=nlegborn)
      integer bflav(nlegs),myin(nlegs),j,ro,ro1,ro2
      double precision p(0:3,nlegs),bres,bornjk(nlegs,nlegs),
     & bmunu(0:3,0:3,nlegs),btemp(4,4,2),bsum1,bsum2,mcfmp(mxpart,4),
     & bmnmad(0:3,0:3,nlegs),born,res1
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
*     q(-p1) +qbar(-p2) -> e-(p3)+e+(p4)+j(p5)+j(p6)

      do j=1,nlegs
      do ro=1,4
      if (myin(j) .gt. 0) then
      if (j .le. 2) then
c incoming partons
      mcfmp(myin(j),ro)=-p(pwhg(ro),j)
      else
      mcfmp(myin(j),ro)=+p(pwhg(ro),j)
      endif
      endif
      enddo
      enddo


      Bmnmad(:,:,:)=zip
      Bmunu(:,:,:)=zip
      bornjk(:,:)=zip
      call qqb_z2jet_pwhg(mcfmp,bflav,bres)
      call qqb_z2jet_cs_pwhg(bflav,bornjk)

C     New Bmunu
      do ro1=1,4
      do ro2=1,4
      Bmunu(pwhg(ro1),pwhg(ro2),:)=Bmcfm(ro1,ro2,:)
      enddo
      enddo
      return
      end

      subroutine finalize_lh
      implicit none
      include 'LesHouches.h'
c     Set up the resonances whose mass must be preserved
c     on the Les Houches interface.
c e- is 11, anti-nue is -12; so
c idup(3)+idup(4) is 0 for Z, 1 for W+, -1 for W-
      if(idup(3)+idup(4).eq.0) then
         call add_resonance(23,3,4)
      else
         call add_resonance(24*(idup(3)+idup(4)),3,4)
      endif
      end


      
