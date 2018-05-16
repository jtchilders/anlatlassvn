      subroutine choose_real_process_full(p,flav,amp2)
      real * 8 p(0:3,1: 5)
      integer flav( 5)
      real * 8 amp2

c     !ER: I want that if no process is called, then amp2 is null
      amp2=0d0

      if ((flav(1).eq.5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call Sbb_twmb(p,amp2) !not resonant

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sbbx_twmbx(p,amp2) !3!

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sbbx_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sbbx_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call Sbd_twmb(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Sbg_twmg(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call Sbs_twmb(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Sbu_twmu(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Sbu_twmu(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call Sbu_twmu(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call Sbu_twmu(p,amp2)









      elseif ((flav(1).eq.5).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Sbux_twmux(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Sbux_twmux(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sbux_twmux(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sbux_twmux(p,amp2)











      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sbxb_twmbx(p,amp2) !3!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sbxb_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sbxb_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sbxd_twmbx(p,amp2) !not

      elseif ((flav(1).eq.-5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sbxs_twmbx(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call Sdb_twmb(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sdbx_twmbx(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call Sdd_twmd(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sddx_twmdx(p,amp2) !3!

      elseif ((flav(1).eq.1).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Sdg_twmg(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Sdu_twmu(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Sdu_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call Sdu_twmu(p,amp2)










      elseif ((flav(1).eq.1).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Sdux_twmux(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Sdux_twmux(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sdux_twmux(p,amp2)






      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sdxd_twmdx(p,amp2) !3!

      elseif ((flav(1).eq.0).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Sgb_twmg(p,amp2) !not

      elseif ((flav(1).eq.0).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Sgd_twmg(p,amp2) !not

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Sgg_twmbx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Sgg_twmdx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sgg_twmsx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Sgs_twmg(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call Ssb_twmb(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Ssbx_twmbx(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call Ssg_twmg(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call Sss_twms(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Sssx_twmsx(p,amp2) !3!

      elseif ((flav(1).eq.3).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Ssu_twmu(p,amp2) !not
      elseif ((flav(1).eq.3).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Ssu_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call Ssu_twmu(p,amp2)







      elseif ((flav(1).eq.3).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Ssux_twmux(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Ssux_twmux(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Ssux_twmux(p,amp2)







      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Ssxs_twmsx(p,amp2) !3!

      elseif ((flav(1).eq.2).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Sub_twmu(p,amp2) !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Sub_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call Sub_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call Sub_twmu(p,amp2)










      elseif ((flav(1).eq.2).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Sud_twmu(p,amp2)  !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Sud_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call Sud_twmu(p,amp2)







      elseif ((flav(1).eq.2).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call Sus_twmu(p,amp2) !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call Sus_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call Sus_twmu(p,amp2)









      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suux_twmbx(p,amp2) !1!

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suux_twmbx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suux_twmbx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suux_twmbx(p,amp2)









      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suux_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suux_twmdx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suux_twmdx(p,amp2)





      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suux_twmsx(p,amp2) !1!
      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suux_twmsx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suux_twmsx(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Suxb_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Suxb_twmux(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suxb_twmux(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suxb_twmux(p,amp2)






      elseif ((flav(1).eq.-2).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Suxd_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Suxd_twmux(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suxd_twmux(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call Suxs_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call Suxs_twmux(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suxs_twmux(p,amp2)









      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suxu_twmbx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call Suxu_twmbx(p,amp2)







      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suxu_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suxu_twmdx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call Suxu_twmdx(p,amp2)






      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suxu_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suxu_twmsx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call Suxu_twmsx(p,amp2)






      else
         write(*,*) 'Error in MAD_real_routines'
         call exit(1)

      endif
      end

      SUBROUTINE Sbb_twmb(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b -> t w- b  
C  
C Crossing   1 is b b -> t w- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bb_twmb
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bb_twmb(maxamps), jamp2bb_twmb(0:maxamps)
      common/to_ampsbb_twmb/  amp2bb_twmb,       jamp2bb_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbb_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bb_twmb(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bb_twmb(ihel)=0d0
              jamp2bb_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bb_twmb(0))
              jamp2bb_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bb_twmb(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bb_twmb(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bb_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bb_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bb_twmb(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b -> t w- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bb_twmb(maxamps), jamp2bb_twmb(0:maxamps)
      common/to_ampsbb_twmb/  amp2bb_twmb,       jamp2bb_twmb
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 5, 1]T[ 3, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     9   )) 
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      bb_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bb_twmb =bb_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bb_twmb(i)=amp2bb_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bb_twmb(i)=Jamp2bb_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbbx_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- b~  
C  
C Crossing   1 is b b~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bbx_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmbx(maxamps), jamp2bbx_twmbx(0:maxamps)
      common/to_ampsbbx_twmbx/  amp2bbx_twmbx,       jamp2bbx_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bbx_twmbx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bbx_twmbx(ihel)=0d0
              jamp2bbx_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_twmbx(0))
              jamp2bbx_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bbx_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmbx(maxamps), jamp2bbx_twmbx(0:maxamps)
      common/to_ampsbbx_twmbx/  amp2bbx_twmbx,       jamp2bbx_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 2, 1]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTB ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      bbx_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_twmbx =bbx_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_twmbx(i)=amp2bbx_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_twmbx(i)=Jamp2bbx_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbbx_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- d~  
C  
C Crossing   1 is b b~ -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bbx_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmdx(maxamps), jamp2bbx_twmdx(0:maxamps)
      common/to_ampsbbx_twmdx/  amp2bbx_twmdx,       jamp2bbx_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bbx_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bbx_twmdx(ihel)=0d0
              jamp2bbx_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_twmdx(0))
              jamp2bbx_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bbx_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bbx_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmdx(maxamps), jamp2bbx_twmdx(0:maxamps)
      common/to_ampsbbx_twmdx/  amp2bbx_twmdx,       jamp2bbx_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bbx_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_twmdx =bbx_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_twmdx(i)=amp2bbx_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_twmdx(i)=Jamp2bbx_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbbx_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- s~  
C  
C Crossing   1 is b b~ -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bbx_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmsx(maxamps), jamp2bbx_twmsx(0:maxamps)
      common/to_ampsbbx_twmsx/  amp2bbx_twmsx,       jamp2bbx_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bbx_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bbx_twmsx(ihel)=0d0
              jamp2bbx_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_twmsx(0))
              jamp2bbx_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bbx_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bbx_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b b~ -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bbx_twmsx(maxamps), jamp2bbx_twmsx(0:maxamps)
      common/to_ampsbbx_twmsx/  amp2bbx_twmsx,       jamp2bbx_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bbx_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_twmsx =bbx_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_twmsx(i)=amp2bbx_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_twmsx(i)=Jamp2bbx_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbd_twmb(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b d -> t w- b  
C  
C Crossing   1 is b d -> t w- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bd_twmb
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bd_twmb(maxamps), jamp2bd_twmb(0:maxamps)
      common/to_ampsbd_twmb/  amp2bd_twmb,       jamp2bd_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbd_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bd_twmb(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bd_twmb(ihel)=0d0
              jamp2bd_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bd_twmb(0))
              jamp2bd_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bd_twmb(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bd_twmb(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bd_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bd_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bd_twmb(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b d -> t w- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bd_twmb(maxamps), jamp2bd_twmb(0:maxamps)
      common/to_ampsbd_twmb/  amp2bd_twmb,       jamp2bd_twmb
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bd_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bd_twmb =bd_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bd_twmb(i)=amp2bd_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bd_twmb(i)=Jamp2bd_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbg_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b g -> t w- g  
C  
C Crossing   1 is b g -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bg_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bg_twmg(maxamps), jamp2bg_twmg(0:maxamps)
      common/to_ampsbg_twmg/  amp2bg_twmg,       jamp2bg_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbg_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2bg_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bg_twmg(ihel)=0d0
              jamp2bg_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bg_twmg(0))
              jamp2bg_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bg_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bg_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bg_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bg_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bg_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b g -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bg_twmg(maxamps), jamp2bg_twmg(0:maxamps)
      common/to_ampsbg_twmg/  amp2bg_twmg,       jamp2bg_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 1, 2, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 1, 5, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     8   ))
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))    
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                        
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      bg_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bg_twmg =bg_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bg_twmg(i)=amp2bg_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bg_twmg(i)=Jamp2bg_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbs_twmb(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b s -> t w- b  
C  
C Crossing   1 is b s -> t w- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bs_twmb
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bs_twmb(maxamps), jamp2bs_twmb(0:maxamps)
      common/to_ampsbs_twmb/  amp2bs_twmb,       jamp2bs_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbs_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bs_twmb(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bs_twmb(ihel)=0d0
              jamp2bs_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bs_twmb(0))
              jamp2bs_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bs_twmb(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bs_twmb(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bs_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bs_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bs_twmb(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b s -> t w- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bs_twmb(maxamps), jamp2bs_twmb(0:maxamps)
      common/to_ampsbs_twmb/  amp2bs_twmb,       jamp2bs_twmb
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bs_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bs_twmb =bs_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bs_twmb(i)=amp2bs_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bs_twmb(i)=Jamp2bs_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbu_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b u -> t w- u  
C  
C Crossing   1 is b u -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bu_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bu_twmu(maxamps), jamp2bu_twmu(0:maxamps)
      common/to_ampsbu_twmu/  amp2bu_twmu,       jamp2bu_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbu_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bu_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bu_twmu(ihel)=0d0
              jamp2bu_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bu_twmu(0))
              jamp2bu_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bu_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bu_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bu_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bu_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bu_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b u -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bu_twmu(maxamps), jamp2bu_twmu(0:maxamps)
      common/to_ampsbu_twmu/  amp2bu_twmu,       jamp2bu_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 3, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))   
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bu_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bu_twmu =bu_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bu_twmu(i)=amp2bu_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bu_twmu(i)=Jamp2bu_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbux_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b u~ -> t w- u~  
C  
C Crossing   1 is b u~ -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bux_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bux_twmux(maxamps), jamp2bux_twmux(0:maxamps)
      common/to_ampsbux_twmux/  amp2bux_twmux,       jamp2bux_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbux_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bux_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bux_twmux(ihel)=0d0
              jamp2bux_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bux_twmux(0))
              jamp2bux_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bux_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bux_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bux_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bux_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bux_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b u~ -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bux_twmux(maxamps), jamp2bux_twmux(0:maxamps)
      common/to_ampsbux_twmux/  amp2bux_twmux,       jamp2bux_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 2, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,1   ))       
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bux_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bux_twmux =bux_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bux_twmux(i)=amp2bux_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bux_twmux(i)=Jamp2bux_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbxb_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- b~  
C  
C Crossing   1 is b~ b -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bxb_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmbx(maxamps), jamp2bxb_twmbx(0:maxamps)
      common/to_ampsbxb_twmbx/  amp2bxb_twmbx,       jamp2bxb_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bxb_twmbx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bxb_twmbx(ihel)=0d0
              jamp2bxb_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_twmbx(0))
              jamp2bxb_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bxb_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bxb_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmbx(maxamps), jamp2bxb_twmbx(0:maxamps)
      common/to_ampsbxb_twmbx/  amp2bxb_twmbx,       jamp2bxb_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 1, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))   
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTB ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      bxb_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_twmbx =bxb_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_twmbx(i)=amp2bxb_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_twmbx(i)=Jamp2bxb_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbxb_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- d~  
C  
C Crossing   1 is b~ b -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bxb_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmdx(maxamps), jamp2bxb_twmdx(0:maxamps)
      common/to_ampsbxb_twmdx/  amp2bxb_twmdx,       jamp2bxb_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxb_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bxb_twmdx(ihel)=0d0
              jamp2bxb_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_twmdx(0))
              jamp2bxb_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bxb_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bxb_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmdx(maxamps), jamp2bxb_twmdx(0:maxamps)
      common/to_ampsbxb_twmdx/  amp2bxb_twmdx,       jamp2bxb_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bxb_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_twmdx =bxb_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_twmdx(i)=amp2bxb_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_twmdx(i)=Jamp2bxb_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbxb_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- s~  
C  
C Crossing   1 is b~ b -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bxb_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmsx(maxamps), jamp2bxb_twmsx(0:maxamps)
      common/to_ampsbxb_twmsx/  amp2bxb_twmsx,       jamp2bxb_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxb_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bxb_twmsx(ihel)=0d0
              jamp2bxb_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_twmsx(0))
              jamp2bxb_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bxb_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bxb_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxb_twmsx(maxamps), jamp2bxb_twmsx(0:maxamps)
      common/to_ampsbxb_twmsx/  amp2bxb_twmsx,       jamp2bxb_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bxb_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_twmsx =bxb_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_twmsx(i)=amp2bxb_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_twmsx(i)=Jamp2bxb_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbxd_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ d -> t w- b~  
C  
C Crossing   1 is b~ d -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bxd_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxd_twmbx(maxamps), jamp2bxd_twmbx(0:maxamps)
      common/to_ampsbxd_twmbx/  amp2bxd_twmbx,       jamp2bxd_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxd_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxd_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bxd_twmbx(ihel)=0d0
              jamp2bxd_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxd_twmbx(0))
              jamp2bxd_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxd_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bxd_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bxd_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxd_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxd_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ d -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxd_twmbx(maxamps), jamp2bxd_twmbx(0:maxamps)
      common/to_ampsbxd_twmbx/  amp2bxd_twmbx,       jamp2bxd_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 1, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bxd_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxd_twmbx =bxd_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxd_twmbx(i)=amp2bxd_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxd_twmbx(i)=Jamp2bxd_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sbxs_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ s -> t w- b~  
C  
C Crossing   1 is b~ s -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 bxs_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxs_twmbx(maxamps), jamp2bxs_twmbx(0:maxamps)
      common/to_ampsbxs_twmbx/  amp2bxs_twmbx,       jamp2bxs_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxs_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxs_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2bxs_twmbx(ihel)=0d0
              jamp2bxs_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxs_twmbx(0))
              jamp2bxs_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxs_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=bxs_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2bxs_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxs_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxs_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ s -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2bxs_twmbx(maxamps), jamp2bxs_twmbx(0:maxamps)
      common/to_ampsbxs_twmbx/  amp2bxs_twmbx,       jamp2bxs_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 1, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bxs_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxs_twmbx =bxs_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxs_twmbx(i)=amp2bxs_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxs_twmbx(i)=Jamp2bxs_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdb_twmb(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d b -> t w- b  
C  
C Crossing   1 is d b -> t w- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 db_twmb
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2db_twmb(maxamps), jamp2db_twmb(0:maxamps)
      common/to_ampsdb_twmb/  amp2db_twmb,       jamp2db_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdb_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2db_twmb(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2db_twmb(ihel)=0d0
              jamp2db_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2db_twmb(0))
              jamp2db_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=db_twmb(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=db_twmb(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2db_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2db_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION db_twmb(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d b -> t w- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2db_twmb(maxamps), jamp2db_twmb(0:maxamps)
      common/to_ampsdb_twmb/  amp2db_twmb,       jamp2db_twmb
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 3, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      db_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          db_twmb =db_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2db_twmb(i)=amp2db_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2db_twmb(i)=Jamp2db_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdbx_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d b~ -> t w- b~  
C  
C Crossing   1 is d b~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 dbx_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dbx_twmbx(maxamps), jamp2dbx_twmbx(0:maxamps)
      common/to_ampsdbx_twmbx/  amp2dbx_twmbx,       jamp2dbx_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdbx_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2dbx_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2dbx_twmbx(ihel)=0d0
              jamp2dbx_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dbx_twmbx(0))
              jamp2dbx_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=dbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2dbx_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dbx_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dbx_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d b~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dbx_twmbx(maxamps), jamp2dbx_twmbx(0:maxamps)
      common/to_ampsdbx_twmbx/  amp2dbx_twmbx,       jamp2dbx_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 2, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                        
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      dbx_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dbx_twmbx =dbx_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dbx_twmbx(i)=amp2dbx_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dbx_twmbx(i)=Jamp2dbx_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdd_twmd(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d d -> t w- d  
C  
C Crossing   1 is d d -> t w- d  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 dd_twmd
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dd_twmd(maxamps), jamp2dd_twmd(0:maxamps)
      common/to_ampsdd_twmd/  amp2dd_twmd,       jamp2dd_twmd

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdd_twmd/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2dd_twmd(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2dd_twmd(ihel)=0d0
              jamp2dd_twmd(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dd_twmd(0))
              jamp2dd_twmd(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dd_twmd(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=dd_twmd(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2dd_twmd(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dd_twmd(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dd_twmd(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d d -> t w- d  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dd_twmd(maxamps), jamp2dd_twmd(0:maxamps)
      common/to_ampsdd_twmd/  amp2dd_twmd,       jamp2dd_twmd
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 5, 1]T[ 3, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     9   ))
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      dd_twmd = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dd_twmd =dd_twmd+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dd_twmd(i)=amp2dd_twmd(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dd_twmd(i)=Jamp2dd_twmd(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sddx_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d d~ -> t w- d~  
C  
C Crossing   1 is d d~ -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ddx_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ddx_twmdx(maxamps), jamp2ddx_twmdx(0:maxamps)
      common/to_ampsddx_twmdx/  amp2ddx_twmdx,       jamp2ddx_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixddx_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ddx_twmdx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ddx_twmdx(ihel)=0d0
              jamp2ddx_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ddx_twmdx(0))
              jamp2ddx_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ddx_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ddx_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ddx_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ddx_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ddx_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d d~ -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ddx_twmdx(maxamps), jamp2ddx_twmdx(0:maxamps)
      common/to_ampsddx_twmdx/  amp2ddx_twmdx,       jamp2ddx_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 2, 1]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTD ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      ddx_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ddx_twmdx =ddx_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ddx_twmdx(i)=amp2ddx_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ddx_twmdx(i)=Jamp2ddx_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdg_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d g -> t w- g  
C  
C Crossing   1 is d g -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 dg_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dg_twmg(maxamps), jamp2dg_twmg(0:maxamps)
      common/to_ampsdg_twmg/  amp2dg_twmg,       jamp2dg_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdg_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2dg_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2dg_twmg(ihel)=0d0
              jamp2dg_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dg_twmg(0))
              jamp2dg_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dg_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=dg_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2dg_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dg_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dg_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d g -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dg_twmg(maxamps), jamp2dg_twmg(0:maxamps)
      common/to_ampsdg_twmg/  amp2dg_twmg,       jamp2dg_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 1, 2, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 1, 5, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))    
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     8   ))
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                         
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      dg_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dg_twmg =dg_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dg_twmg(i)=amp2dg_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dg_twmg(i)=Jamp2dg_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdu_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d u -> t w- u  
C  
C Crossing   1 is d u -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 du_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2du_twmu(maxamps), jamp2du_twmu(0:maxamps)
      common/to_ampsdu_twmu/  amp2du_twmu,       jamp2du_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdu_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2du_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2du_twmu(ihel)=0d0
              jamp2du_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2du_twmu(0))
              jamp2du_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=du_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=du_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2du_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2du_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION du_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d u -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2du_twmu(maxamps), jamp2du_twmu(0:maxamps)
      common/to_ampsdu_twmu/  amp2du_twmu,       jamp2du_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 3, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      du_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          du_twmu =du_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2du_twmu(i)=amp2du_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2du_twmu(i)=Jamp2du_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdux_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d u~ -> t w- u~  
C  
C Crossing   1 is d u~ -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 dux_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dux_twmux(maxamps), jamp2dux_twmux(0:maxamps)
      common/to_ampsdux_twmux/  amp2dux_twmux,       jamp2dux_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdux_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2dux_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2dux_twmux(ihel)=0d0
              jamp2dux_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dux_twmux(0))
              jamp2dux_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dux_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=dux_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2dux_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dux_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dux_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d u~ -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dux_twmux(maxamps), jamp2dux_twmux(0:maxamps)
      common/to_ampsdux_twmux/  amp2dux_twmux,       jamp2dux_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 2, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))    
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      dux_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dux_twmux =dux_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dux_twmux(i)=amp2dux_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dux_twmux(i)=Jamp2dux_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sdxd_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d~ d -> t w- d~  
C  
C Crossing   1 is d~ d -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 dxd_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dxd_twmdx(maxamps), jamp2dxd_twmdx(0:maxamps)
      common/to_ampsdxd_twmdx/  amp2dxd_twmdx,       jamp2dxd_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdxd_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2dxd_twmdx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2dxd_twmdx(ihel)=0d0
              jamp2dxd_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dxd_twmdx(0))
              jamp2dxd_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dxd_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=dxd_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2dxd_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dxd_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dxd_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d~ d -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2dxd_twmdx(maxamps), jamp2dxd_twmdx(0:maxamps)
      common/to_ampsdxd_twmdx/  amp2dxd_twmdx,       jamp2dxd_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 1, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTD ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      dxd_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dxd_twmdx =dxd_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dxd_twmdx(i)=amp2dxd_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dxd_twmdx(i)=Jamp2dxd_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgb_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g b -> t w- g  
C  
C Crossing   1 is g b -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gb_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gb_twmg(maxamps), jamp2gb_twmg(0:maxamps)
      common/to_ampsgb_twmg/  amp2gb_twmg,       jamp2gb_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgb_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gb_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gb_twmg(ihel)=0d0
              jamp2gb_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gb_twmg(0))
              jamp2gb_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gb_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gb_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gb_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gb_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gb_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g b -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gb_twmg(maxamps), jamp2gb_twmg(0:maxamps)
      common/to_ampsgb_twmg/  amp2gb_twmg,       jamp2gb_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 2, 1, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 2, 5, 1]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                  
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                   
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gb_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gb_twmg =gb_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gb_twmg(i)=amp2gb_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gb_twmg(i)=Jamp2gb_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgd_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g d -> t w- g  
C  
C Crossing   1 is g d -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gd_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gd_twmg(maxamps), jamp2gd_twmg(0:maxamps)
      common/to_ampsgd_twmg/  amp2gd_twmg,       jamp2gd_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgd_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gd_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gd_twmg(ihel)=0d0
              jamp2gd_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gd_twmg(0))
              jamp2gd_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gd_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gd_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gd_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gd_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gd_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g d -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gd_twmg(maxamps), jamp2gd_twmg(0:maxamps)
      common/to_ampsgd_twmg/  amp2gd_twmg,       jamp2gd_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 2, 1, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 2, 5, 1]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                          
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gd_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gd_twmg =gd_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gd_twmg(i)=amp2gd_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gd_twmg(i)=Jamp2gd_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgg_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- b~  
C  
C Crossing   1 is g g -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gg_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmbx(maxamps), jamp2gg_twmbx(0:maxamps)
      common/to_ampsgg_twmbx/  amp2gg_twmbx,       jamp2gg_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_twmbx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gg_twmbx(ihel)=0d0
              jamp2gg_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_twmbx(0))
              jamp2gg_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gg_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gg_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  16, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !: ausiliary wave functions
      complex*16 w6tmp(6),w9tmp(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmbx(maxamps), jamp2gg_twmbx(0:maxamps)
      common/to_ampsgg_twmbx/  amp2gg_twmbx,       jamp2gg_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c     !: create amp1 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
c$$$      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
c$$$     &     7   ))   
c$$$      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     7   ))   
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
c     !: create amp2 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTB ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W6tmp(1    ))     
      CALL FVOXXX(W6tmp(1    ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))   
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTB ,AMP(2   ))     
c     !: create amp3 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
c$$$      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
c$$$     &     10  ))                                                          
c$$$      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
c     !: create amp4 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTB ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W9tmp(1    ))     
      CALL FVOXXX(W9tmp(1    ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTB ,AMP(4   ))      
c     !: create amp5
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
c     !: create amp6
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
c     !: create amp7 with NONZERO twidth
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTB ,AMP(7   ))          
c     !: create amp8 
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_twmbx =gg_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_twmbx(i)=amp2gg_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_twmbx(i)=Jamp2gg_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgg_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- d~  
C  
C Crossing   1 is g g -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gg_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmdx(maxamps), jamp2gg_twmdx(0:maxamps)
      common/to_ampsgg_twmdx/  amp2gg_twmdx,       jamp2gg_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_twmdx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gg_twmdx(ihel)=0d0
              jamp2gg_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_twmdx(0))
              jamp2gg_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gg_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gg_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  16, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !: ausiliary wave functions
      complex*16 w6tmp(6),w9tmp(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmdx(maxamps), jamp2gg_twmdx(0:maxamps)
      common/to_ampsgg_twmdx/  amp2gg_twmdx,       jamp2gg_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c     !: create amp1 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
c$$$      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTD ,BMASS   ,ZERO    ,W(1,        
c$$$     &     7   ))   
c$$$      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTD ,BMASS   ,ZERO    ,W(1,        
     &     7   ))   
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
c     !: create amp2 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTD ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W6tmp(1    ))     
      CALL FVOXXX(W6tmp(1    ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))   
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTD ,AMP(2   ))     
c     !: create amp3 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
c$$$      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTD ,BMASS   ,ZERO    ,W(1,        
c$$$     &     10  ))                                                          
c$$$      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTD ,BMASS   ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
c     !: create amp4 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTD ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W9tmp(1    ))     
      CALL FVOXXX(W9tmp(1    ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTD ,AMP(4   ))      
c     !: create amp5
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,BMASS   ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
c     !: create amp6
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
c     !: create amp7 with NONZERO twidth
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTD ,AMP(7   ))          
c     !: create amp8 
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_twmdx =gg_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_twmdx(i)=amp2gg_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_twmdx(i)=Jamp2gg_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgg_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- s~  
C  
C Crossing   1 is g g -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gg_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmsx(maxamps), jamp2gg_twmsx(0:maxamps)
      common/to_ampsgg_twmsx/  amp2gg_twmsx,       jamp2gg_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_twmsx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) / 256/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gg_twmsx(ihel)=0d0
              jamp2gg_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_twmsx(0))
              jamp2gg_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gg_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gg_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g g -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  16, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !: ausiliary wave functions
      complex*16 w6tmp(6),w9tmp(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_twmsx(maxamps), jamp2gg_twmsx(0:maxamps)
      common/to_ampsgg_twmsx/  amp2gg_twmsx,       jamp2gg_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 5, 2, 1]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 5, 1, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c     !: create amp1 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
c$$$      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTS ,BMASS   ,ZERO    ,W(1,        
c$$$     &     7   ))   
c$$$      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTS ,BMASS   ,ZERO    ,W(1,        
     &     7   ))   
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))     
c     !: create amp2 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTS ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W6tmp(1    ))     
      CALL FVOXXX(W6tmp(1    ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))   
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTS ,AMP(2   ))     
c     !: create amp3 with twidth=0
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
c$$$      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTS ,BMASS   ,ZERO    ,W(1,        
c$$$     &     10  ))                                                          
c$$$      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTS ,BMASS   ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
c     !: create amp4 with NONZERO twidth
c$$$      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
c$$$      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTS ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W9tmp(1    ))     
      CALL FVOXXX(W9tmp(1    ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTS ,AMP(4   ))      
c     !: create amp5
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,BMASS   ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
c     !: create amp6
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
c     !: create amp7 with NONZERO twidth
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTS ,AMP(7   ))          
c     !: create amp8 
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_twmsx =gg_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_twmsx(i)=amp2gg_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_twmsx(i)=Jamp2gg_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sgs_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g s -> t w- g  
C  
C Crossing   1 is g s -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 gs_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gs_twmg(maxamps), jamp2gs_twmg(0:maxamps)
      common/to_ampsgs_twmg/  amp2gs_twmg,       jamp2gs_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgs_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gs_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2gs_twmg(ihel)=0d0
              jamp2gs_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gs_twmg(0))
              jamp2gs_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gs_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=gs_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2gs_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gs_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gs_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g s -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gs_twmg(maxamps), jamp2gs_twmg(0:maxamps)
      common/to_ampsgs_twmg/  amp2gs_twmg,       jamp2gs_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 2, 1, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 2, 5, 1]                                             
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,ZERO  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                          
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gs_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gs_twmg =gs_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gs_twmg(i)=amp2gs_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gs_twmg(i)=Jamp2gs_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssb_twmb(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s b -> t w- b  
C  
C Crossing   1 is s b -> t w- b  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 sb_twmb
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sb_twmb(maxamps), jamp2sb_twmb(0:maxamps)
      common/to_ampssb_twmb/  amp2sb_twmb,       jamp2sb_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsb_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sb_twmb(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2sb_twmb(ihel)=0d0
              jamp2sb_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sb_twmb(0))
              jamp2sb_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sb_twmb(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=sb_twmb(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2sb_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sb_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sb_twmb(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s b -> t w- b  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sb_twmb(maxamps), jamp2sb_twmb(0:maxamps)
      common/to_ampssb_twmb/  amp2sb_twmb,       jamp2sb_twmb
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 3, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),BMASS ,NHEL(5   ),+1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      sb_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sb_twmb =sb_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sb_twmb(i)=amp2sb_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sb_twmb(i)=Jamp2sb_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssbx_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s b~ -> t w- b~  
C  
C Crossing   1 is s b~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 sbx_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sbx_twmbx(maxamps), jamp2sbx_twmbx(0:maxamps)
      common/to_ampssbx_twmbx/  amp2sbx_twmbx,       jamp2sbx_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsbx_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sbx_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2sbx_twmbx(ihel)=0d0
              jamp2sbx_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sbx_twmbx(0))
              jamp2sbx_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=sbx_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2sbx_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sbx_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sbx_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s b~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sbx_twmbx(maxamps), jamp2sbx_twmbx(0:maxamps)
      common/to_ampssbx_twmbx/  amp2sbx_twmbx,       jamp2sbx_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 2, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),BMASS ,NHEL(2   ),-1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))   
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      sbx_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sbx_twmbx =sbx_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sbx_twmbx(i)=amp2sbx_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sbx_twmbx(i)=Jamp2sbx_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssg_twmg(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s g -> t w- g  
C  
C Crossing   1 is s g -> t w- g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 sg_twmg
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sg_twmg(maxamps), jamp2sg_twmg(0:maxamps)
      common/to_ampssg_twmg/  amp2sg_twmg,       jamp2sg_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsg_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2sg_twmg(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2sg_twmg(ihel)=0d0
              jamp2sg_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sg_twmg(0))
              jamp2sg_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sg_twmg(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=sg_twmg(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2sg_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sg_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sg_twmg(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s g -> t w- g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   8,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  14, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sg_twmg(maxamps), jamp2sg_twmg(0:maxamps)
      common/to_ampssg_twmg/  amp2sg_twmg,       jamp2sg_twmg
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /    16,   -2/                            
C               T[ 3, 1, 2, 5]                                             
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,   16/                            
C               T[ 3, 1, 5, 2]                                             
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,ZERO  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     8   ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
c$$$      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,ZERO  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
c$$$      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     14  ))                                                          
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      sg_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sg_twmg =sg_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sg_twmg(i)=amp2sg_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sg_twmg(i)=Jamp2sg_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sss_twms(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s s -> t w- s  
C  
C Crossing   1 is s s -> t w- s  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ss_twms
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ss_twms(maxamps), jamp2ss_twms(0:maxamps)
      common/to_ampsss_twms/  amp2ss_twms,       jamp2ss_twms

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixss_twms/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ss_twms(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ss_twms(ihel)=0d0
              jamp2ss_twms(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ss_twms(0))
              jamp2ss_twms(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ss_twms(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ss_twms(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ss_twms(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ss_twms(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ss_twms(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s s -> t w- s  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ss_twms(maxamps), jamp2ss_twms(0:maxamps)
      common/to_ampsss_twms/  amp2ss_twms,       jamp2ss_twms
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 5, 1]T[ 3, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     9   ))                                                          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      ss_twms = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ss_twms =ss_twms+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ss_twms(i)=amp2ss_twms(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ss_twms(i)=Jamp2ss_twms(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sssx_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s s~ -> t w- s~  
C  
C Crossing   1 is s s~ -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ssx_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ssx_twmsx(maxamps), jamp2ssx_twmsx(0:maxamps)
      common/to_ampsssx_twmsx/  amp2ssx_twmsx,       jamp2ssx_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixssx_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ssx_twmsx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ssx_twmsx(ihel)=0d0
              jamp2ssx_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ssx_twmsx(0))
              jamp2ssx_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ssx_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ssx_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ssx_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ssx_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ssx_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s s~ -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ssx_twmsx(maxamps), jamp2ssx_twmsx(0:maxamps)
      common/to_ampsssx_twmsx/  amp2ssx_twmsx,       jamp2ssx_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 2, 1]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTS ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      ssx_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ssx_twmsx =ssx_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ssx_twmsx(i)=amp2ssx_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ssx_twmsx(i)=Jamp2ssx_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssu_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s u -> t w- u  
C  
C Crossing   1 is s u -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 su_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2su_twmu(maxamps), jamp2su_twmu(0:maxamps)
      common/to_ampssu_twmu/  amp2su_twmu,       jamp2su_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsu_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2su_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2su_twmu(ihel)=0d0
              jamp2su_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2su_twmu(0))
              jamp2su_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=su_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=su_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2su_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2su_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION su_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s u -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2su_twmu(maxamps), jamp2su_twmu(0:maxamps)
      common/to_ampssu_twmu/  amp2su_twmu,       jamp2su_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 1]T[ 3, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      su_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          su_twmu =su_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2su_twmu(i)=amp2su_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2su_twmu(i)=Jamp2su_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssux_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s u~ -> t w- u~  
C  
C Crossing   1 is s u~ -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 sux_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sux_twmux(maxamps), jamp2sux_twmux(0:maxamps)
      common/to_ampssux_twmux/  amp2sux_twmux,       jamp2sux_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsux_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sux_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2sux_twmux(ihel)=0d0
              jamp2sux_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sux_twmux(0))
              jamp2sux_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sux_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=sux_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2sux_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sux_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sux_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s u~ -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sux_twmux(maxamps), jamp2sux_twmux(0:maxamps)
      common/to_ampssux_twmux/  amp2sux_twmux,       jamp2sux_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 2, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                        
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      sux_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sux_twmux =sux_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sux_twmux(i)=amp2sux_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sux_twmux(i)=Jamp2sux_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Ssxs_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : s~ s -> t w- s~  
C  
C Crossing   1 is s~ s -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 sxs_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sxs_twmsx(maxamps), jamp2sxs_twmsx(0:maxamps)
      common/to_ampssxs_twmsx/  amp2sxs_twmsx,       jamp2sxs_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsxs_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2sxs_twmsx(0) /   2/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2sxs_twmsx(ihel)=0d0
              jamp2sxs_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sxs_twmsx(0))
              jamp2sxs_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sxs_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=sxs_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2sxs_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sxs_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sxs_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : s~ s -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   4,NEIGEN=  2) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  10, NCOLOR=   2) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2sxs_twmsx(maxamps), jamp2sxs_twmsx(0:maxamps)
      common/to_ampssxs_twmsx/  amp2sxs_twmsx,       jamp2sxs_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            3/                                       
      DATA (CF(i,1  ),i=1  ,2  ) /     6,   -2/                            
C               T[ 3, 5]T[ 1, 2]                                           
      DATA Denom(2  )/            3/                                       
      DATA (CF(i,2  ),i=1  ,2  ) /    -2,    6/                            
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL FVOXXX(W(1,3   ),W(1,9   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,4   ),GWFTS ,AMP(3   ))          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,9   ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      sxs_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sxs_twmsx =sxs_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sxs_twmsx(i)=amp2sxs_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sxs_twmsx(i)=Jamp2sxs_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sub_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u b -> t w- u  
C  
C Crossing   1 is u b -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ub_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ub_twmu(maxamps), jamp2ub_twmu(0:maxamps)
      common/to_ampsub_twmu/  amp2ub_twmu,       jamp2ub_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixub_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2ub_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ub_twmu(ihel)=0d0
              jamp2ub_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ub_twmu(0))
              jamp2ub_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ub_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ub_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ub_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ub_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ub_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u b -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ub_twmu(maxamps), jamp2ub_twmu(0:maxamps)
      common/to_ampsub_twmu/  amp2ub_twmu,       jamp2ub_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      ub_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ub_twmu =ub_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ub_twmu(i)=amp2ub_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ub_twmu(i)=Jamp2ub_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sud_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u d -> t w- u  
C  
C Crossing   1 is u d -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 ud_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ud_twmu(maxamps), jamp2ud_twmu(0:maxamps)
      common/to_ampsud_twmu/  amp2ud_twmu,       jamp2ud_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixud_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2ud_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2ud_twmu(ihel)=0d0
              jamp2ud_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ud_twmu(0))
              jamp2ud_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ud_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=ud_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2ud_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ud_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ud_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u d -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2ud_twmu(maxamps), jamp2ud_twmu(0:maxamps)
      common/to_ampsud_twmu/  amp2ud_twmu,       jamp2ud_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                          
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      ud_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ud_twmu =ud_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ud_twmu(i)=amp2ud_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ud_twmu(i)=Jamp2ud_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Sus_twmu(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u s -> t w- u  
C  
C Crossing   1 is u s -> t w- u  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 us_twmu
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2us_twmu(maxamps), jamp2us_twmu(0:maxamps)
      common/to_ampsus_twmu/  amp2us_twmu,       jamp2us_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixus_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2us_twmu(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2us_twmu(ihel)=0d0
              jamp2us_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2us_twmu(0))
              jamp2us_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=us_twmu(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=us_twmu(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2us_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2us_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION us_twmu(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u s -> t w- u  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2us_twmu(maxamps), jamp2us_twmu(0:maxamps)
      common/to_ampsus_twmu/  amp2us_twmu,       jamp2us_twmu
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 5, 2]T[ 3, 1]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))                                                        
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      us_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          us_twmu =us_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2us_twmu(i)=amp2us_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2us_twmu(i)=Jamp2us_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suux_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- b~  
C  
C Crossing   1 is u u~ -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uux_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmbx(maxamps), jamp2uux_twmbx(0:maxamps)
      common/to_ampsuux_twmbx/  amp2uux_twmbx,       jamp2uux_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uux_twmbx(ihel)=0d0
              jamp2uux_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_twmbx(0))
              jamp2uux_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uux_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uux_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmbx(maxamps), jamp2uux_twmbx(0:maxamps)
      common/to_ampsuux_twmbx/  amp2uux_twmbx,       jamp2uux_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_twmbx =uux_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_twmbx(i)=amp2uux_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_twmbx(i)=Jamp2uux_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suux_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- d~  
C  
C Crossing   1 is u u~ -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uux_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmdx(maxamps), jamp2uux_twmdx(0:maxamps)
      common/to_ampsuux_twmdx/  amp2uux_twmdx,       jamp2uux_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uux_twmdx(ihel)=0d0
              jamp2uux_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_twmdx(0))
              jamp2uux_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uux_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uux_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmdx(maxamps), jamp2uux_twmdx(0:maxamps)
      common/to_ampsuux_twmdx/  amp2uux_twmdx,       jamp2uux_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_twmdx =uux_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_twmdx(i)=amp2uux_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_twmdx(i)=Jamp2uux_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suux_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- s~  
C  
C Crossing   1 is u u~ -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uux_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmsx(maxamps), jamp2uux_twmsx(0:maxamps)
      common/to_ampsuux_twmsx/  amp2uux_twmsx,       jamp2uux_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uux_twmsx(ihel)=0d0
              jamp2uux_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_twmsx(0))
              jamp2uux_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uux_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uux_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uux_twmsx(maxamps), jamp2uux_twmsx(0:maxamps)
      common/to_ampsuux_twmsx/  amp2uux_twmsx,       jamp2uux_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 1]T[ 2, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_twmsx =uux_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_twmsx(i)=amp2uux_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_twmsx(i)=Jamp2uux_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxb_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ b -> t w- u~  
C  
C Crossing   1 is u~ b -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxb_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxb_twmux(maxamps), jamp2uxb_twmux(0:maxamps)
      common/to_ampsuxb_twmux/  amp2uxb_twmux,       jamp2uxb_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxb_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxb_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxb_twmux(ihel)=0d0
              jamp2uxb_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxb_twmux(0))
              jamp2uxb_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxb_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxb_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxb_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxb_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxb_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ b -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxb_twmux(maxamps), jamp2uxb_twmux(0:maxamps)
      common/to_ampsuxb_twmux/  amp2uxb_twmux,       jamp2uxb_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 1, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,ZERO  ,W(1,        
     &     6   )) 
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxb_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxb_twmux =uxb_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxb_twmux(i)=amp2uxb_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxb_twmux(i)=Jamp2uxb_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxd_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> t w- u~  
C  
C Crossing   1 is u~ d -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxd_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxd_twmux(maxamps), jamp2uxd_twmux(0:maxamps)
      common/to_ampsuxd_twmux/  amp2uxd_twmux,       jamp2uxd_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxd_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxd_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxd_twmux(ihel)=0d0
              jamp2uxd_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxd_twmux(0))
              jamp2uxd_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxd_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxd_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxd_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxd_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxd_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ d -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxd_twmux(maxamps), jamp2uxd_twmux(0:maxamps)
      common/to_ampsuxd_twmux/  amp2uxd_twmux,       jamp2uxd_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 1, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   ))  
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxd_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxd_twmux =uxd_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxd_twmux(i)=amp2uxd_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxd_twmux(i)=Jamp2uxd_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxs_twmux(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ s -> t w- u~  
C  
C Crossing   1 is u~ s -> t w- u~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxs_twmux
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxs_twmux(maxamps), jamp2uxs_twmux(0:maxamps)
      common/to_ampsuxs_twmux/  amp2uxs_twmux,       jamp2uxs_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxs_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxs_twmux(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxs_twmux(ihel)=0d0
              jamp2uxs_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxs_twmux(0))
              jamp2uxs_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxs_twmux(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxs_twmux(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxs_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxs_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxs_twmux(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ s -> t w- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxs_twmux(maxamps), jamp2uxs_twmux(0:maxamps)
      common/to_ampsuxs_twmux/  amp2uxs_twmux,       jamp2uxs_twmux
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 5]T[ 1, 2]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
c$$$      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
c$$$     &     6   )) 
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,ZERO  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxs_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxs_twmux =uxs_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxs_twmux(i)=amp2uxs_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxs_twmux(i)=Jamp2uxs_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxu_twmbx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- b~  
C  
C Crossing   1 is u~ u -> t w- b~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxu_twmbx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmbx(maxamps), jamp2uxu_twmbx(0:maxamps)
      common/to_ampsuxu_twmbx/  amp2uxu_twmbx,       jamp2uxu_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_twmbx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxu_twmbx(ihel)=0d0
              jamp2uxu_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_twmbx(0))
              jamp2uxu_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_twmbx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxu_twmbx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxu_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_twmbx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- b~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmbx(maxamps), jamp2uxu_twmbx(0:maxamps)
      common/to_ampsuxu_twmbx/  amp2uxu_twmbx,       jamp2uxu_twmbx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),BMASS ,NHEL(5   ),-1*IC(5   ),W(1,5   ))       
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTB ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_twmbx =uxu_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_twmbx(i)=amp2uxu_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_twmbx(i)=Jamp2uxu_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxu_twmdx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- d~  
C  
C Crossing   1 is u~ u -> t w- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxu_twmdx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmdx(maxamps), jamp2uxu_twmdx(0:maxamps)
      common/to_ampsuxu_twmdx/  amp2uxu_twmdx,       jamp2uxu_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_twmdx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxu_twmdx(ihel)=0d0
              jamp2uxu_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_twmdx(0))
              jamp2uxu_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_twmdx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxu_twmdx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxu_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_twmdx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmdx(maxamps), jamp2uxu_twmdx(0:maxamps)
      common/to_ampsuxu_twmdx/  amp2uxu_twmdx,       jamp2uxu_twmdx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTD ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_twmdx =uxu_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_twmdx(i)=amp2uxu_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_twmdx(i)=Jamp2uxu_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE Suxu_twmsx(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- s~  
C  
C Crossing   1 is u~ u -> t w- s~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps_r.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  48, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 uxu_twmsx
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
      REAL     xran1
      EXTERNAL xran1
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmsx(maxamps), jamp2uxu_twmsx(0:maxamps)
      common/to_ampsuxu_twmsx/  amp2uxu_twmsx,       jamp2uxu_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_twmsx(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) /-1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) /-1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1,-1,-1, 0,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1,-1,-1, 0, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  33),IHEL=1, 5) / 1,-1, 1, 0,-1/
      DATA (NHEL(IHEL,  34),IHEL=1, 5) / 1,-1, 1, 0, 1/
      DATA (NHEL(IHEL,  35),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  36),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  37),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  38),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  39),IHEL=1, 5) / 1, 1,-1, 0,-1/
      DATA (NHEL(IHEL,  40),IHEL=1, 5) / 1, 1,-1, 0, 1/
      DATA (NHEL(IHEL,  41),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  42),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  43),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  44),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  45),IHEL=1, 5) / 1, 1, 1, 0,-1/
      DATA (NHEL(IHEL,  46),IHEL=1, 5) / 1, 1, 1, 0, 1/
      DATA (NHEL(IHEL,  47),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  48),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
       
      IF (multi_channel) THEN
          DO IHEL=1,NGRAPHS
              amp2uxu_twmsx(ihel)=0d0
              jamp2uxu_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_twmsx(0))
              jamp2uxu_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_twmsx(P ,NHEL(1,IHEL),JC(1))            
               DO JJ=1,nincoming
                 IF(POL(JJ).NE.1d0.AND.
     &              NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                   T=T*ABS(POL(JJ))
                 ELSE IF(POL(JJ).NE.1d0)THEN
                   T=T*(2d0-ABS(POL(JJ)))
                 ENDIF
               ENDDO
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
          JHEL = 1
          ISUM_HEL=MIN(ISUM_HEL,NGOOD)
      ELSE              !RANDOM HELICITY
          DO J=1,ISUM_HEL
              JHEL=JHEL+1
              IF (JHEL .GT. NGOOD) JHEL=1
              HWGT = REAL(NGOOD)/REAL(ISUM_HEL)
              IHEL = IGOOD(JHEL)
              T=uxu_twmsx(P ,NHEL(1,IHEL),JC(1))            
              DO JJ=1,nincoming
                IF(POL(JJ).NE.1d0.AND.
     &             NHEL(JJ,IHEL).EQ.INT(SIGN(1d0,POL(JJ)))) THEN
                  T=T*ABS(POL(JJ))
                ELSE IF(POL(JJ).NE.1d0)THEN
                  T=T*(2d0-ABS(POL(JJ)))
                ENDIF
              ENDDO
              ANS(IPROC)=ANS(IPROC)+T*HWGT
          ENDDO
          IF (ISUM_HEL .EQ. 1) THEN
              WRITE(HEL_BUFF,'(16i5)')(NHEL(i,IHEL),i=1,nexternal)
          ENDIF
      ENDIF
      IF (MULTI_CHANNEL) THEN
          XTOT=0D0
          DO IHEL=1,MAPCONFIG(0)
              XTOT=XTOT+AMP2uxu_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_twmsx(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ u -> t w- s~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=   2,NEIGEN=  1) 
      include "genps_r.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   8, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2uxu_twmsx(maxamps), jamp2uxu_twmsx(0:maxamps)
      common/to_ampsuxu_twmsx/  amp2uxu_twmsx,       jamp2uxu_twmsx
      include "../coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     2/                                  
C               T[ 3, 2]T[ 1, 5]                                           
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),ZERO ,NHEL(1   ),-1*IC(1   ),W(1,1   ))        
      CALL IXXXXX(P(0,2   ),ZERO ,NHEL(2   ),+1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),WMASS ,NHEL(4   ),+1*IC(4   ),W(1,4   ))       
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL JIOXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,6   ))     
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,4   ),GWFTS ,AMP(1   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,6   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_twmsx =uxu_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_twmsx(i)=amp2uxu_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_twmsx(i)=Jamp2uxu_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
