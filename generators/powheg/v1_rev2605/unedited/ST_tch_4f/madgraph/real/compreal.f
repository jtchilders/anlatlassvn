      subroutine compreal(p,flav,amp2)
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
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scc_tbxcd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scc_tbxcs(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scd_tbxdd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scdx_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scdx_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scdx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scdx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scdx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scg_tbxdg(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scg_tbxsg(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scs_tbxds(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scs_tbxss(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scsx_tbxssx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scsx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scsx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scsx_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxcd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxcs(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxcs(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxud(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxus(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scu_tbxus(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxdx_tbxcxcx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxdx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxdx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxsx_tbxcxcx(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxsx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Scxsx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdc_tbxdd(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sddx_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sddx_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdsx_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdsx_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdu_tbxdd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxc_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxc_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxcx_tbxcxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxcx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxcx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxd_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxd_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxdx_tbxcxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxdx_tbxuxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxg_tbxcxg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxg_tbxuxg(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxsx_tbxcxdx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxsx_tbxcxsx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxu_tbxuux(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sdxux_tbxuxux(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgc_tbxdg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgc_tbxsg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgdx_tbxcxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgdx_tbxuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_tbxcxs(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgg_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgsx_tbxcxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgsx_tbxuxg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_tbxdg(madp,amp2)
      return
      endif

      madflav(1)=0
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sgu_tbxsg(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssc_tbxds(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssc_tbxss(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssu_tbxds(madp,amp2)
      return
      endif

      madflav(1)=3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssu_tbxss(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxc_tbxssx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxc_tbxuxc(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxcx_tbxcxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxcx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxcx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxd_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxd_tbxcxd(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxdx_tbxcxdx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxdx_tbxcxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxg_tbxcxg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxg_tbxuxg(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxsx_tbxcxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxsx_tbxuxsx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxu_tbxssx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxu_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxu_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxu_tbxuux(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxux_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-3
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Ssxux_tbxuxux(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxcd(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxcs(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxcs(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxud(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxus(madp,amp2)
      return
      endif

      madflav(1)=1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suc_tbxus(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sud_tbxdd(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_tbxddx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sudx_tbxuux(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_tbxdg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=0
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=0
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sug_tbxsg(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sus_tbxds(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Sus_tbxss(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=3
      madflav(6)=-3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Susx_tbxssx(madp,amp2)

c$$$      print*, '--------------------'
c$$$      print*, madflav
c$$$      print*, flav
c$$$      print*, perm


      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Susx_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=4
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=4
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Susx_tbxucx(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Susx_tbxuux(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_tbxud(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suu_tbxus(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suux_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=4
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxc_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-1
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxdx_tbxuxux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-4
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxsx_tbxuxcx(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=-3
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=-2
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxsx_tbxuxux(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=1
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_tbxuxd(madp,amp2)
      return
      endif

      madflav(1)=-2
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-2
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-4
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-4
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_tbxuxs(madp,amp2)
      return
      endif

      madflav(1)=-1
      madflav(2)=2
      madflav(3)=6
      madflav(4)=-5
      madflav(5)=-1
      madflav(6)=3
      if (flavequiv_perm_real( 6,flav,madflav,perm)) then
      do i=1, 6 
         do mu=0,3
            madp(mu,perm(i))=p(mu,i)
         enddo
      enddo

      call Suxu_tbxuxs(madp,amp2)
      return
      endif

      write(*,*) 'ERROR: the flavour list', flav
      write(*,*) 'is not in the list of MADGRAPH routines'
      call exit(1)
      end

