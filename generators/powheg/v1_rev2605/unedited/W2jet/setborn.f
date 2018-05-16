      subroutine setborn(p,bflav,bres,bornjk,bmunu)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'constants.f'
      include 'plabel.f'
      include 'process.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'mcfmtopwhg.f'
      include 'ckmallowed.f'
      include 'mmsq_cs.f'
      include 'Bmcfm.f'
      integer nlegs,ro1,ro2
      parameter (nlegs=nlegborn)
      integer bflav(nlegs),myin(nlegs)
      integer i7,i4,i3,i8,i5,i6,fl,j,k,mu,nu,mu1,nu1,ro
      double precision p(0:3,nlegs),bres,bornjk(nlegs,nlegs),
     &  bmunu(0:3,0:3,nlegs),btemp(4,4,2),bsum1,bsum2,mcfmp(mxpart,4),
     &  msq(-nf:nf,-nf:nf),B(0:3,0:3,nlegs),
     &  bornmad,bornjkmad(nlegs,nlegs),bmnmad(0:3,0:3,nlegs)
      logical iseven,isodd,isaquark
      external qqb_w2jet_gvec
c set mcfm strong coupling and scales
      call st_mcfm
c ---
      i3=3
      i4=4
      i5=5
      i6=6
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
      mcfmp(myin(j),ro)=-p(pwhg(ro),j)
      else
      mcfmp(myin(j),ro)=+p(pwhg(ro),j)
      endif
      endif
      enddo
      enddo

      B(:,:,:)=zip
      Bmunu(:,:,:)=zip
      bornjk(:,:)=zip
      call qqb_w2jet_pwhg(mcfmp,bflav,bres)
      call qqb_w2jet_cs_pwhg(bflav,bornjk)

C     Old Bmunu
c      call setupBmunu(mcfmp,myin,qqb_w2jet_gvec,nlegs,bflav,B)
c      if ((bflav(1) == 0) .and. (bflav(2) == 0)) then
C-----remove sum over final state quarks included by qqb_w2jet_gvec
c      B(:,:,:)=0.5d0*B(:,:,:) 
c      endif

C     New Bmunu
      do ro1=1,4
      do ro2=1,4
      Bmunu(pwhg(ro1),pwhg(ro2),:)=Bmcfm(ro1,ro2,:)
      enddo
      enddo

c      call sborn_proc(p,bflav,bornmad,bornjkmad,bmnmad)
c      write(6,*) 'bflav',bflav
c      write(6,*) 'bres   ',bres
c      write(6,*) 'bornmad',bornmad
c      pause

c      do mu=0,3
c      do nu=0,3
c      do j=1,2
c      If (abs((B(mu,nu,j)-Bmunu(mu,nu,j))/Bres) .gt. 1d-7) then
c      write(6,*) 'j',j
c      write(6,*) 'gvec'
c      write(6,*) B(0,0,j),B(0,1,j),B(0,2,j),B(0,3,j)
c      write(6,*) B(1,0,j),B(1,1,j),B(1,2,j),B(1,3,j)
c      write(6,*) B(2,0,j),B(2,1,j),B(2,2,j),B(2,3,j)
c      write(6,*) B(3,0,j),B(3,1,j),B(3,2,j),B(3,3,j)
c      write(6,*) 'mcfm'
c      write(6,*) Bmunu(0,0,j),Bmunu(0,1,j),Bmunu(0,2,j),Bmunu(0,3,j)
c      write(6,*) Bmunu(1,0,j),Bmunu(1,1,j),Bmunu(1,2,j),Bmunu(1,3,j)
c      write(6,*) Bmunu(2,0,j),Bmunu(2,1,j),Bmunu(2,2,j),Bmunu(2,3,j)
c      write(6,*) Bmunu(3,0,j),Bmunu(3,1,j),Bmunu(3,2,j),Bmunu(3,3,j)
c      pause
c      enddo
c      enddo
c      enddo
 

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

      
