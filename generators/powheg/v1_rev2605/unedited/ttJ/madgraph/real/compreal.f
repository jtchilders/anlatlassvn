      subroutine compreal(p,flav,amp2)
      implicit none
      real * 8 p(0:3,1: 6)
      integer flav( 6)
      real * 8 amp2
      real * 8 madp(0:3,1: 6)
      integer madflav( 6),perm( 6)
      logical flavequiv_perm
      external flavequiv_perm
      integer i, mu

      madflav(1)=5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbb_ttxbb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbg_ttxgb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_ttxub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_ttxub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_ttxub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_ttxub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxbx_ttxbxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxg_ttxgbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scg_ttxgc(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxg_ttxgcx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_ttxud(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgb_ttxgb(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgbx_ttxgbx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgc_ttxgc(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgcx_ttxgcx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_ttxub(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_ttxub(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_ttxub(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_ttxub(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_ttxubx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_ttxud(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_ttxudx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxug(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxug(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxug(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_ttxuu(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_ttxuu(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_ttxuu(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_ttxuu(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      madflav(6)=5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_ttxuxb(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_ttxuxbx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_ttxuxd(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_ttxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxuxg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      madflav(6)=-5
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxbbx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxddx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      madflav(6)=0
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxgg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      madflav(6)=-2
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_ttxuxux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      madflav(6)=-1
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_ttxuxux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      madflav(6)=-3
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_ttxuxux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      madflav(6)=-4
      if (flavequiv_perm( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_ttxuxux(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

