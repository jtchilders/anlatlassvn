      subroutine compborn_decay(p,flav,amp2)
      implicit none
      real * 8 p(0:3,1: 9)
      integer flav( 9)
      real * 8 amp2
      real * 8 madp(0:3,1: 9)
      integer madflav( 9),perm( 9)
      logical flavequiv_perm
      external flavequiv_perm
      integer i, mu

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=5
      madflav(4)=-11
      madflav(5)=12
      madflav(6)=-5
      madflav(7)=11
      madflav(8)=-12
      madflav(9)=0
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbg_bepvebxemvexb(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxg_bepvebxemvexbx(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgb_bepvebxemvexb(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgbx_bepvebxemvexbx(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_bepvebxemvexu(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_bepvebxemvexux(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexg(madp,amp2)
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
      if (flavequiv_perm( 9,flav,madflav,perm)) then
      do i=1, 9 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_bepvebxemvexg(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

