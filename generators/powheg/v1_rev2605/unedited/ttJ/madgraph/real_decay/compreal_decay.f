      subroutine compreal_decay(p,flav,amp2)
      implicit none
      real * 8 p(0:3,1: 10)
      integer flav( 10)
      real * 8 amp2
      real * 8 madp(0:3,1: 10)
      integer madflav( 10),perm( 10)
      logical flavequiv_perm
      external flavequiv_perm
      integer i, mu

      madflav(1)=5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbb_bepvebxemvexbb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbg_bepvebxemvexbg(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbu_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=5
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbux_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxbx_bepvebxemvexbxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-5
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxg_bepvebxemvexbxg(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxu_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-5
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxux_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgb_bepvebxemvexbg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-5
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgbx_bepvebxemvexbxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sub_bepvebxemvexub(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Subx_bepvebxemvexubx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_bepvebxemvexuc(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_bepvebxemvexucx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexug(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_bepvebxemvexuu(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_bepvebxemvexuu(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_bepvebxemvexuu(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_bepvebxemvexuu(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxb_bepvebxemvexuxb(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxbx_bepvebxemvexuxbx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_bepvebxemvexuxc(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_bepvebxemvexuxcx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexuxg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=5
      madflav(10)=-5
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexbbx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexccx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      madflav(10)=0
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexgg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-2
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-2
      madflav(10)=-2
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_bepvebxemvexuxux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-4
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-4
      madflav(10)=-4
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_bepvebxemvexuxux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-1
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-1
      madflav(10)=-1
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_bepvebxemvexuxux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-3
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=-3
      madflav(10)=-3
      if (flavequiv_perm( 10,flav,madflav,perm)) then
      do i=1, 10 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_bepvebxemvexuxux(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

