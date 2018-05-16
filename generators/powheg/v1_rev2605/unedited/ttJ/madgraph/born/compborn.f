      subroutine compborn(p,flav,amp2,amp2jk,amp2munu)
      implicit none
      real * 8 p(0:3,1: 5)
      integer flav( 5)
      real * 8 amp2,amp2jk( 5, 5),amp2munu(0:3,0:3, 5)
      real * 8 madp(0:3,1: 5)
      integer madflav( 5),perm( 5)
      logical flavequiv_perm
      external flavequiv_perm
      integer i, mu

      madflav(1)=5
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbbx_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=5
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbg_ttxb(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-5
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxb_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-5
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-5
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sbxg_ttxbx(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=5
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgb_ttxb(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=-5
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-5
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgbx_ttxbx(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=0
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_ttxu(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-2
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-1
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-3
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=-4
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_ttxux(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-6
      madflav(5)=0
      if (flavequiv_perm( 5,flav,madflav,perm)) then
      do i=1, 5 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_ttxg(madp,amp2,amp2jk,amp2munu)
      call reorder_born_perm(perm,amp2jk,amp2munu)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

  
      subroutine reorder_born_perm(perm,bjk,bmunu)
      implicit none
      integer i,j,mu,nu,perm(  5)
      real * 8 bjk(  5,   5),bmunu(0:3,0:3,  5) 
      real * 8 bjkx(  5,   5),bmunux(0:3,0:3,  5)
c Local copy
      do i=1,  5
         do j=i+1,  5
            bjkx(i,j)=bjk(i,j)
            bjkx(j,i)=bjkx(i,j)
         enddo 
         do mu=0,3
            do nu=0,3
               bmunux(mu,nu,i)=bmunu(mu,nu,i)
            enddo
         enddo
      enddo 
c
      do i=1,  5
         do j=i+1,  5
            bjk(i,j)=bjkx(perm(i),perm(j))
            bjk(j,i)=bjk(i,j)
        enddo 
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,i)=bmunux(mu,nu,perm(i))
            enddo
         enddo
      enddo
      end 

