      subroutine compreal_tb(p,flav,amp2)
      implicit none
      real * 8 p(0:3,1: 6)
      integer flav( 6)
      real * 8 amp2
      real * 8 madp(0:3,1: 6)
      integer madflav( 6),perm( 6)
      logical flavequiv_perm_real
      external flavequiv_perm_real
      integer i, mu

      madflav(1)=4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scd_txbcc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scd_txbuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scd_txbuc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scs_txbcc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scs_txbuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scs_txbuc(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxcx_txbcxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxcx_txbcxsx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxd_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxd_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxd_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxdx_txbdxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxg_txbdxg(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxg_txbsxg(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxs_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxs_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxs_txbccx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxs_txbdxs(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxsx_txbdxsx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxsx_txbsxsx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbcxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbcxsx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbcxsx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbuxdx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbuxsx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxux_txbuxsx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdc_txbcc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdc_txbuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdc_txbuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdcx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdcx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdcx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdd_txbcd(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdd_txbud(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdg_txbcg(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdg_txbug(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sds_txbcd(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sds_txbcs(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_txbuc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_txbuu(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdux_txbuxc(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxcx_txbdxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_txbdxdx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgcx_txbdxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgcx_txbsxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgd_txbcg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgd_txbug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_txbcdx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_txbcsx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_txbudx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_txbusx(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgs_txbcg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgs_txbug(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_txbdxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgux_txbsxg(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssc_txbcc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssc_txbuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssc_txbuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sscx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sscx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sscx_txbccx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sscx_txbdxs(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssd_txbcd(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssd_txbcs(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssg_txbcg(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssg_txbug(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sss_txbcs(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sss_txbus(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssu_txbuc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssu_txbuu(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbdxs(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbuux(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssux_txbuxc(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxcx_txbdxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxcx_txbsxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxux_txbdxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxux_txbsxsx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbudx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbudx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbudx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbusx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbusx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbusx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sucx_txbusx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_txbuc(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_txbuu(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sus_txbuc(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sus_txbuu(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbudx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbudx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbudx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbusx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbusx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbusx(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbcxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbcxsx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbcxsx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbuxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbuxsx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxcx_txbuxsx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxd_txbuxc(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_txbdxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_txbdxg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=0
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxg_txbsxg(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbdxs(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbuux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxs_txbuxc(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxsx_txbdxsx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxsx_txbsxsx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbudx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=1
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=3
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_txbusx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_txbuxdx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-2
      madflav(3)=-6
      madflav(4)=5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxux_txbuxsx(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

