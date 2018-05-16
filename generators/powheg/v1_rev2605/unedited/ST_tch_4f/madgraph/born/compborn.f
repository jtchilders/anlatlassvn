      subroutine compborn(p,flav,amp2)
      implicit none
      real * 8 p(0:3,1: 5)
      integer flav( 5)
      real * 8 amp2
      real * 8 madp(0:3,1: 5)
      integer madflav( 5),perm( 5)
      logical flavequiv_perm
      external flavequiv_perm
      integer i, mu

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scg_tbxd(madp,amp2)

c$$$ccccccccccccccccccccccccccccccc
c$$$      print*, amp2
c$$$      stop
c$$$ccccccccccccccccccccccccccccccc

      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scg_tbxs(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxg_tbxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxg_tbxux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgc_tbxd(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgc_tbxs(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgdx_tbxcx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgdx_tbxux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgsx_tbxcx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgsx_tbxux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_tbxd(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_tbxs(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxg_tbxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxg_tbxux(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_tbxd(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_tbxs(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

