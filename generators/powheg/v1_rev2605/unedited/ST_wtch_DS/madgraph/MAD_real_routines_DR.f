      subroutine choose_real_process_DR(p,flav,amp2)
      real * 8 p(0:3,1: 5)
      integer flav( 5)
      real * 8 amp2

      if ((flav(1).eq.5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call DR_Sbb_twmb(p,amp2) !not resonant

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sbbx_twmbx(p,amp2) !3!

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sbbx_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.5).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sbbx_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call DR_Sbd_twmb(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Sbg_twmg(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call DR_Sbs_twmb(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Sbu_twmu(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Sbu_twmu(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call DR_Sbu_twmu(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call DR_Sbu_twmu(p,amp2)









      elseif ((flav(1).eq.5).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Sbux_twmux(p,amp2) !not

      elseif ((flav(1).eq.5).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Sbux_twmux(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sbux_twmux(p,amp2)

      elseif ((flav(1).eq.5).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sbux_twmux(p,amp2)











      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sbxb_twmbx(p,amp2) !3!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sbxb_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sbxb_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.-5).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sbxd_twmbx(p,amp2) !not

      elseif ((flav(1).eq.-5).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sbxs_twmbx(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call DR_Sdb_twmb(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sdbx_twmbx(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call DR_Sdd_twmd(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sddx_twmdx(p,amp2) !3!

      elseif ((flav(1).eq.1).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Sdg_twmg(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Sdu_twmu(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Sdu_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call DR_Sdu_twmu(p,amp2)










      elseif ((flav(1).eq.1).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Sdux_twmux(p,amp2) !not

      elseif ((flav(1).eq.1).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Sdux_twmux(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sdux_twmux(p,amp2)






      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sdxd_twmdx(p,amp2) !3!

      elseif ((flav(1).eq.0).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Sgb_twmg(p,amp2) !not

      elseif ((flav(1).eq.0).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Sgd_twmg(p,amp2) !not

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Sgg_twmbx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Sgg_twmdx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sgg_twmsx(p,amp2) !2,4,7!

      elseif ((flav(1).eq.0).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Sgs_twmg(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.5)
     $)then
      call DR_Ssb_twmb(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Ssbx_twmbx(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.0)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.0)
     $)then
      call DR_Ssg_twmg(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call DR_Sss_twms(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Sssx_twmsx(p,amp2) !3!

      elseif ((flav(1).eq.3).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Ssu_twmu(p,amp2) !not
      elseif ((flav(1).eq.3).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Ssu_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call DR_Ssu_twmu(p,amp2)







      elseif ((flav(1).eq.3).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Ssux_twmux(p,amp2) !not

      elseif ((flav(1).eq.3).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Ssux_twmux(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Ssux_twmux(p,amp2)







      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Ssxs_twmsx(p,amp2) !3!

      elseif ((flav(1).eq.2).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Sub_twmu(p,amp2) !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Sub_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call DR_Sub_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call DR_Sub_twmu(p,amp2)










      elseif ((flav(1).eq.2).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Sud_twmu(p,amp2)  !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Sud_twmu(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.3)
     $)then
      call DR_Sud_twmu(p,amp2)







      elseif ((flav(1).eq.2).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.2)
     $)then
      call DR_Sus_twmu(p,amp2) !not

      elseif ((flav(1).eq.4).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.4)
     $)then
      call DR_Sus_twmu(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.1)
     $)then
      call DR_Sus_twmu(p,amp2)









      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suux_twmbx(p,amp2) !1!

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suux_twmbx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suux_twmbx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suux_twmbx(p,amp2)









      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suux_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suux_twmdx(p,amp2)

      elseif ((flav(1).eq.3).and.(flav(2).eq.-3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suux_twmdx(p,amp2)





      elseif ((flav(1).eq.2).and.(flav(2).eq.-2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suux_twmsx(p,amp2) !1!
      elseif ((flav(1).eq.4).and.(flav(2).eq.-4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suux_twmsx(p,amp2)

      elseif ((flav(1).eq.1).and.(flav(2).eq.-1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suux_twmsx(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Suxb_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Suxb_twmux(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suxb_twmux(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.5)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suxb_twmux(p,amp2)






      elseif ((flav(1).eq.-2).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Suxd_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Suxd_twmux(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suxd_twmux(p,amp2)








      elseif ((flav(1).eq.-2).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-2)
     $)then
      call DR_Suxs_twmux(p,amp2) !not

      elseif ((flav(1).eq.-4).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-4)
     $)then
      call DR_Suxs_twmux(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suxs_twmux(p,amp2)









      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suxu_twmbx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suxu_twmbx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-5)
     $)then
      call DR_Suxu_twmbx(p,amp2)







      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suxu_twmdx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suxu_twmdx(p,amp2)

      elseif ((flav(1).eq.-3).and.(flav(2).eq.3)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-1)
     $)then
      call DR_Suxu_twmdx(p,amp2)






      elseif ((flav(1).eq.-2).and.(flav(2).eq.2)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suxu_twmsx(p,amp2) !1!

      elseif ((flav(1).eq.-4).and.(flav(2).eq.4)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suxu_twmsx(p,amp2)

      elseif ((flav(1).eq.-1).and.(flav(2).eq.1)
     $.and.(flav(3).eq.6)
     $.and.(flav(4).eq.-24)
     $.and.(flav(5).eq.-3)
     $)then
      call DR_Suxu_twmsx(p,amp2)






      else
         write(*,*) 'Error in MAD_real_routines'
         call exit(1)

      endif
      end

      SUBROUTINE DR_Sbb_twmb(P1,ANS)
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
      REAL*8 bb_DR_twmb
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
      Double Precision amp2bb_DR_twmb(maxamps), jamp2bb_DR_twmb(0:maxamps)
      common/to_ampsbb_DR_twmb/  amp2bb_DR_twmb,       jamp2bb_DR_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbb_DR_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bb_DR_twmb(0) /   2/          
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
              amp2bb_DR_twmb(ihel)=0d0
              jamp2bb_DR_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bb_DR_twmb(0))
              jamp2bb_DR_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bb_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              T=bb_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bb_DR_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bb_DR_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bb_DR_twmb(P,NHEL,IC)
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
      Double Precision amp2bb_DR_twmb(maxamps), jamp2bb_DR_twmb(0:maxamps)
      common/to_ampsbb_DR_twmb/  amp2bb_DR_twmb,       jamp2bb_DR_twmb
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      bb_DR_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bb_DR_twmb =bb_DR_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bb_DR_twmb(i)=amp2bb_DR_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bb_DR_twmb(i)=Jamp2bb_DR_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbbx_twmbx(P1,ANS)
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
      REAL*8 bbx_DR_twmbx
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
      Double Precision amp2bbx_DR_twmbx(maxamps), jamp2bbx_DR_twmbx(0:maxamps)
      common/to_ampsbbx_DR_twmbx/  amp2bbx_DR_twmbx,       jamp2bbx_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bbx_DR_twmbx(0) /   2/          
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
              amp2bbx_DR_twmbx(ihel)=0d0
              jamp2bbx_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_DR_twmbx(0))
              jamp2bbx_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bbx_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2bbx_DR_twmbx(maxamps), jamp2bbx_DR_twmbx(0:maxamps)
      common/to_ampsbbx_DR_twmbx/  amp2bbx_DR_twmbx,       jamp2bbx_DR_twmbx
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      bbx_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_DR_twmbx =bbx_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_DR_twmbx(i)=amp2bbx_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_DR_twmbx(i)=Jamp2bbx_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbbx_twmdx(P1,ANS)
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
      REAL*8 bbx_DR_twmdx
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
      Double Precision amp2bbx_DR_twmdx(maxamps), jamp2bbx_DR_twmdx(0:maxamps)
      common/to_ampsbbx_DR_twmdx/  amp2bbx_DR_twmdx,       jamp2bbx_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bbx_DR_twmdx(0) /   1/          
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
              amp2bbx_DR_twmdx(ihel)=0d0
              jamp2bbx_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_DR_twmdx(0))
              jamp2bbx_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bbx_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bbx_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2bbx_DR_twmdx(maxamps), jamp2bbx_DR_twmdx(0:maxamps)
      common/to_ampsbbx_DR_twmdx/  amp2bbx_DR_twmdx,       jamp2bbx_DR_twmdx
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
      amp(1)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bbx_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_DR_twmdx =bbx_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_DR_twmdx(i)=amp2bbx_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_DR_twmdx(i)=Jamp2bbx_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbbx_twmsx(P1,ANS)
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
      REAL*8 bbx_DR_twmsx
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
      Double Precision amp2bbx_DR_twmsx(maxamps), jamp2bbx_DR_twmsx(0:maxamps)
      common/to_ampsbbx_DR_twmsx/  amp2bbx_DR_twmsx,       jamp2bbx_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbbx_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bbx_DR_twmsx(0) /   1/          
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
              amp2bbx_DR_twmsx(ihel)=0d0
              jamp2bbx_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bbx_DR_twmsx(0))
              jamp2bbx_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bbx_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bbx_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bbx_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bbx_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bbx_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2bbx_DR_twmsx(maxamps), jamp2bbx_DR_twmsx(0:maxamps)
      common/to_ampsbbx_DR_twmsx/  amp2bbx_DR_twmsx,       jamp2bbx_DR_twmsx
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
      amp(1)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bbx_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bbx_DR_twmsx =bbx_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bbx_DR_twmsx(i)=amp2bbx_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bbx_DR_twmsx(i)=Jamp2bbx_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbd_twmb(P1,ANS)
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
      REAL*8 bd_DR_twmb
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
      Double Precision amp2bd_DR_twmb(maxamps), jamp2bd_DR_twmb(0:maxamps)
      common/to_ampsbd_DR_twmb/  amp2bd_DR_twmb,       jamp2bd_DR_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbd_DR_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bd_DR_twmb(0) /   1/          
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
              amp2bd_DR_twmb(ihel)=0d0
              jamp2bd_DR_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bd_DR_twmb(0))
              jamp2bd_DR_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bd_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              T=bd_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bd_DR_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bd_DR_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bd_DR_twmb(P,NHEL,IC)
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
      Double Precision amp2bd_DR_twmb(maxamps), jamp2bd_DR_twmb(0:maxamps)
      common/to_ampsbd_DR_twmb/  amp2bd_DR_twmb,       jamp2bd_DR_twmb
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bd_DR_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bd_DR_twmb =bd_DR_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bd_DR_twmb(i)=amp2bd_DR_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bd_DR_twmb(i)=Jamp2bd_DR_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbg_twmg(P1,ANS)
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
      REAL*8 bg_DR_twmg
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
      Double Precision amp2bg_DR_twmg(maxamps), jamp2bg_DR_twmg(0:maxamps)
      common/to_ampsbg_DR_twmg/  amp2bg_DR_twmg,       jamp2bg_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbg_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2bg_DR_twmg(0) /   2/          
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
              amp2bg_DR_twmg(ihel)=0d0
              jamp2bg_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bg_DR_twmg(0))
              jamp2bg_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=bg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bg_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bg_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bg_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2bg_DR_twmg(maxamps), jamp2bg_DR_twmg(0:maxamps)
      common/to_ampsbg_DR_twmg/  amp2bg_DR_twmg,       jamp2bg_DR_twmg
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
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      bg_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bg_DR_twmg =bg_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bg_DR_twmg(i)=amp2bg_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bg_DR_twmg(i)=Jamp2bg_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbs_twmb(P1,ANS)
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
      REAL*8 bs_DR_twmb
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
      Double Precision amp2bs_DR_twmb(maxamps), jamp2bs_DR_twmb(0:maxamps)
      common/to_ampsbs_DR_twmb/  amp2bs_DR_twmb,       jamp2bs_DR_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbs_DR_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bs_DR_twmb(0) /   1/          
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
              amp2bs_DR_twmb(ihel)=0d0
              jamp2bs_DR_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bs_DR_twmb(0))
              jamp2bs_DR_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bs_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              T=bs_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bs_DR_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bs_DR_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bs_DR_twmb(P,NHEL,IC)
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
      Double Precision amp2bs_DR_twmb(maxamps), jamp2bs_DR_twmb(0:maxamps)
      common/to_ampsbs_DR_twmb/  amp2bs_DR_twmb,       jamp2bs_DR_twmb
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bs_DR_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bs_DR_twmb =bs_DR_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bs_DR_twmb(i)=amp2bs_DR_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bs_DR_twmb(i)=Jamp2bs_DR_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbu_twmu(P1,ANS)
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
      REAL*8 bu_DR_twmu
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
      Double Precision amp2bu_DR_twmu(maxamps), jamp2bu_DR_twmu(0:maxamps)
      common/to_ampsbu_DR_twmu/  amp2bu_DR_twmu,       jamp2bu_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbu_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bu_DR_twmu(0) /   1/          
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
              amp2bu_DR_twmu(ihel)=0d0
              jamp2bu_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bu_DR_twmu(0))
              jamp2bu_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bu_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=bu_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bu_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bu_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bu_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2bu_DR_twmu(maxamps), jamp2bu_DR_twmu(0:maxamps)
      common/to_ampsbu_DR_twmu/  amp2bu_DR_twmu,       jamp2bu_DR_twmu
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bu_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bu_DR_twmu =bu_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bu_DR_twmu(i)=amp2bu_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bu_DR_twmu(i)=Jamp2bu_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbux_twmux(P1,ANS)
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
      REAL*8 bux_DR_twmux
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
      Double Precision amp2bux_DR_twmux(maxamps), jamp2bux_DR_twmux(0:maxamps)
      common/to_ampsbux_DR_twmux/  amp2bux_DR_twmux,       jamp2bux_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbux_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bux_DR_twmux(0) /   1/          
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
              amp2bux_DR_twmux(ihel)=0d0
              jamp2bux_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bux_DR_twmux(0))
              jamp2bux_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=bux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bux_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bux_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bux_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2bux_DR_twmux(maxamps), jamp2bux_DR_twmux(0:maxamps)
      common/to_ampsbux_DR_twmux/  amp2bux_DR_twmux,       jamp2bux_DR_twmux
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bux_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bux_DR_twmux =bux_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bux_DR_twmux(i)=amp2bux_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bux_DR_twmux(i)=Jamp2bux_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbxb_twmbx(P1,ANS)
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
      REAL*8 bxb_DR_twmbx
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
      Double Precision amp2bxb_DR_twmbx(maxamps), jamp2bxb_DR_twmbx(0:maxamps)
      common/to_ampsbxb_DR_twmbx/  amp2bxb_DR_twmbx,       jamp2bxb_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2bxb_DR_twmbx(0) /   2/          
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
              amp2bxb_DR_twmbx(ihel)=0d0
              jamp2bxb_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_DR_twmbx(0))
              jamp2bxb_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bxb_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bxb_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2bxb_DR_twmbx(maxamps), jamp2bxb_DR_twmbx(0:maxamps)
      common/to_ampsbxb_DR_twmbx/  amp2bxb_DR_twmbx,       jamp2bxb_DR_twmbx
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      bxb_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_DR_twmbx =bxb_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_DR_twmbx(i)=amp2bxb_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_DR_twmbx(i)=Jamp2bxb_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbxb_twmdx(P1,ANS)
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
      REAL*8 bxb_DR_twmdx
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
      Double Precision amp2bxb_DR_twmdx(maxamps), jamp2bxb_DR_twmdx(0:maxamps)
      common/to_ampsbxb_DR_twmdx/  amp2bxb_DR_twmdx,       jamp2bxb_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxb_DR_twmdx(0) /   1/          
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
              amp2bxb_DR_twmdx(ihel)=0d0
              jamp2bxb_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_DR_twmdx(0))
              jamp2bxb_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bxb_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bxb_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2bxb_DR_twmdx(maxamps), jamp2bxb_DR_twmdx(0:maxamps)
      common/to_ampsbxb_DR_twmdx/  amp2bxb_DR_twmdx,       jamp2bxb_DR_twmdx
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
      amp(1)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bxb_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_DR_twmdx =bxb_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_DR_twmdx(i)=amp2bxb_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_DR_twmdx(i)=Jamp2bxb_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbxb_twmsx(P1,ANS)
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
      REAL*8 bxb_DR_twmsx
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
      Double Precision amp2bxb_DR_twmsx(maxamps), jamp2bxb_DR_twmsx(0:maxamps)
      common/to_ampsbxb_DR_twmsx/  amp2bxb_DR_twmsx,       jamp2bxb_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxb_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxb_DR_twmsx(0) /   1/          
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
              amp2bxb_DR_twmsx(ihel)=0d0
              jamp2bxb_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxb_DR_twmsx(0))
              jamp2bxb_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxb_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bxb_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bxb_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxb_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxb_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2bxb_DR_twmsx(maxamps), jamp2bxb_DR_twmsx(0:maxamps)
      common/to_ampsbxb_DR_twmsx/  amp2bxb_DR_twmsx,       jamp2bxb_DR_twmsx
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
      amp(1)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      bxb_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxb_DR_twmsx =bxb_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxb_DR_twmsx(i)=amp2bxb_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxb_DR_twmsx(i)=Jamp2bxb_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbxd_twmbx(P1,ANS)
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
      REAL*8 bxd_DR_twmbx
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
      Double Precision amp2bxd_DR_twmbx(maxamps), jamp2bxd_DR_twmbx(0:maxamps)
      common/to_ampsbxd_DR_twmbx/  amp2bxd_DR_twmbx,       jamp2bxd_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxd_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxd_DR_twmbx(0) /   1/          
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
              amp2bxd_DR_twmbx(ihel)=0d0
              jamp2bxd_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxd_DR_twmbx(0))
              jamp2bxd_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxd_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bxd_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bxd_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxd_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxd_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2bxd_DR_twmbx(maxamps), jamp2bxd_DR_twmbx(0:maxamps)
      common/to_ampsbxd_DR_twmbx/  amp2bxd_DR_twmbx,       jamp2bxd_DR_twmbx
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bxd_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxd_DR_twmbx =bxd_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxd_DR_twmbx(i)=amp2bxd_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxd_DR_twmbx(i)=Jamp2bxd_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sbxs_twmbx(P1,ANS)
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
      REAL*8 bxs_DR_twmbx
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
      Double Precision amp2bxs_DR_twmbx(maxamps), jamp2bxs_DR_twmbx(0:maxamps)
      common/to_ampsbxs_DR_twmbx/  amp2bxs_DR_twmbx,       jamp2bxs_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixbxs_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2bxs_DR_twmbx(0) /   1/          
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
              amp2bxs_DR_twmbx(ihel)=0d0
              jamp2bxs_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2bxs_DR_twmbx(0))
              jamp2bxs_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=bxs_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=bxs_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2bxs_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2bxs_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION bxs_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2bxs_DR_twmbx(maxamps), jamp2bxs_DR_twmbx(0:maxamps)
      common/to_ampsbxs_DR_twmbx/  amp2bxs_DR_twmbx,       jamp2bxs_DR_twmbx
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      bxs_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          bxs_DR_twmbx =bxs_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2bxs_DR_twmbx(i)=amp2bxs_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2bxs_DR_twmbx(i)=Jamp2bxs_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdb_twmb(P1,ANS)
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
      REAL*8 db_DR_twmb
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
      Double Precision amp2db_DR_twmb(maxamps), jamp2db_DR_twmb(0:maxamps)
      common/to_ampsdb_DR_twmb/  amp2db_DR_twmb,       jamp2db_DR_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdb_DR_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2db_DR_twmb(0) /   1/          
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
              amp2db_DR_twmb(ihel)=0d0
              jamp2db_DR_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2db_DR_twmb(0))
              jamp2db_DR_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=db_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              T=db_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2db_DR_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2db_DR_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION db_DR_twmb(P,NHEL,IC)
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
      Double Precision amp2db_DR_twmb(maxamps), jamp2db_DR_twmb(0:maxamps)
      common/to_ampsdb_DR_twmb/  amp2db_DR_twmb,       jamp2db_DR_twmb
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      db_DR_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          db_DR_twmb =db_DR_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2db_DR_twmb(i)=amp2db_DR_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2db_DR_twmb(i)=Jamp2db_DR_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdbx_twmbx(P1,ANS)
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
      REAL*8 dbx_DR_twmbx
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
      Double Precision amp2dbx_DR_twmbx(maxamps), jamp2dbx_DR_twmbx(0:maxamps)
      common/to_ampsdbx_DR_twmbx/  amp2dbx_DR_twmbx,       jamp2dbx_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdbx_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2dbx_DR_twmbx(0) /   1/          
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
              amp2dbx_DR_twmbx(ihel)=0d0
              jamp2dbx_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dbx_DR_twmbx(0))
              jamp2dbx_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=dbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2dbx_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dbx_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dbx_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2dbx_DR_twmbx(maxamps), jamp2dbx_DR_twmbx(0:maxamps)
      common/to_ampsdbx_DR_twmbx/  amp2dbx_DR_twmbx,       jamp2dbx_DR_twmbx
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      dbx_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dbx_DR_twmbx =dbx_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dbx_DR_twmbx(i)=amp2dbx_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dbx_DR_twmbx(i)=Jamp2dbx_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdd_twmd(P1,ANS)
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
      REAL*8 dd_DR_twmd
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
      Double Precision amp2dd_DR_twmd(maxamps), jamp2dd_DR_twmd(0:maxamps)
      common/to_ampsdd_DR_twmd/  amp2dd_DR_twmd,       jamp2dd_DR_twmd

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdd_DR_twmd/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2dd_DR_twmd(0) /   2/          
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
              amp2dd_DR_twmd(ihel)=0d0
              jamp2dd_DR_twmd(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dd_DR_twmd(0))
              jamp2dd_DR_twmd(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dd_DR_twmd(P ,NHEL(1,IHEL),JC(1))            
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
              T=dd_DR_twmd(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2dd_DR_twmd(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dd_DR_twmd(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dd_DR_twmd(P,NHEL,IC)
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
      Double Precision amp2dd_DR_twmd(maxamps), jamp2dd_DR_twmd(0:maxamps)
      common/to_ampsdd_DR_twmd/  amp2dd_DR_twmd,       jamp2dd_DR_twmd
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      dd_DR_twmd = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dd_DR_twmd =dd_DR_twmd+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dd_DR_twmd(i)=amp2dd_DR_twmd(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dd_DR_twmd(i)=Jamp2dd_DR_twmd(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sddx_twmdx(P1,ANS)
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
      REAL*8 ddx_DR_twmdx
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
      Double Precision amp2ddx_DR_twmdx(maxamps), jamp2ddx_DR_twmdx(0:maxamps)
      common/to_ampsddx_DR_twmdx/  amp2ddx_DR_twmdx,       jamp2ddx_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixddx_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ddx_DR_twmdx(0) /   2/          
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
              amp2ddx_DR_twmdx(ihel)=0d0
              jamp2ddx_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ddx_DR_twmdx(0))
              jamp2ddx_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ddx_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=ddx_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2ddx_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ddx_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ddx_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2ddx_DR_twmdx(maxamps), jamp2ddx_DR_twmdx(0:maxamps)
      common/to_ampsddx_DR_twmdx/  amp2ddx_DR_twmdx,       jamp2ddx_DR_twmdx
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      ddx_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ddx_DR_twmdx =ddx_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ddx_DR_twmdx(i)=amp2ddx_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ddx_DR_twmdx(i)=Jamp2ddx_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdg_twmg(P1,ANS)
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
      REAL*8 dg_DR_twmg
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
      Double Precision amp2dg_DR_twmg(maxamps), jamp2dg_DR_twmg(0:maxamps)
      common/to_ampsdg_DR_twmg/  amp2dg_DR_twmg,       jamp2dg_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdg_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2dg_DR_twmg(0) /   2/          
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
              amp2dg_DR_twmg(ihel)=0d0
              jamp2dg_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dg_DR_twmg(0))
              jamp2dg_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=dg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2dg_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dg_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dg_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2dg_DR_twmg(maxamps), jamp2dg_DR_twmg(0:maxamps)
      common/to_ampsdg_DR_twmg/  amp2dg_DR_twmg,       jamp2dg_DR_twmg
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
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      dg_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dg_DR_twmg =dg_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dg_DR_twmg(i)=amp2dg_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dg_DR_twmg(i)=Jamp2dg_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdu_twmu(P1,ANS)
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
      REAL*8 du_DR_twmu
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
      Double Precision amp2du_DR_twmu(maxamps), jamp2du_DR_twmu(0:maxamps)
      common/to_ampsdu_DR_twmu/  amp2du_DR_twmu,       jamp2du_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdu_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2du_DR_twmu(0) /   1/          
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
              amp2du_DR_twmu(ihel)=0d0
              jamp2du_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2du_DR_twmu(0))
              jamp2du_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=du_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=du_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2du_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2du_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION du_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2du_DR_twmu(maxamps), jamp2du_DR_twmu(0:maxamps)
      common/to_ampsdu_DR_twmu/  amp2du_DR_twmu,       jamp2du_DR_twmu
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      du_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          du_DR_twmu =du_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2du_DR_twmu(i)=amp2du_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2du_DR_twmu(i)=Jamp2du_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdux_twmux(P1,ANS)
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
      REAL*8 dux_DR_twmux
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
      Double Precision amp2dux_DR_twmux(maxamps), jamp2dux_DR_twmux(0:maxamps)
      common/to_ampsdux_DR_twmux/  amp2dux_DR_twmux,       jamp2dux_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdux_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2dux_DR_twmux(0) /   1/          
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
              amp2dux_DR_twmux(ihel)=0d0
              jamp2dux_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dux_DR_twmux(0))
              jamp2dux_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=dux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2dux_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dux_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dux_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2dux_DR_twmux(maxamps), jamp2dux_DR_twmux(0:maxamps)
      common/to_ampsdux_DR_twmux/  amp2dux_DR_twmux,       jamp2dux_DR_twmux
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      dux_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dux_DR_twmux =dux_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dux_DR_twmux(i)=amp2dux_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dux_DR_twmux(i)=Jamp2dux_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sdxd_twmdx(P1,ANS)
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
      REAL*8 dxd_DR_twmdx
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
      Double Precision amp2dxd_DR_twmdx(maxamps), jamp2dxd_DR_twmdx(0:maxamps)
      common/to_ampsdxd_DR_twmdx/  amp2dxd_DR_twmdx,       jamp2dxd_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixdxd_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2dxd_DR_twmdx(0) /   2/          
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
              amp2dxd_DR_twmdx(ihel)=0d0
              jamp2dxd_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2dxd_DR_twmdx(0))
              jamp2dxd_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=dxd_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=dxd_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2dxd_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2dxd_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION dxd_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2dxd_DR_twmdx(maxamps), jamp2dxd_DR_twmdx(0:maxamps)
      common/to_ampsdxd_DR_twmdx/  amp2dxd_DR_twmdx,       jamp2dxd_DR_twmdx
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      dxd_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          dxd_DR_twmdx =dxd_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2dxd_DR_twmdx(i)=amp2dxd_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2dxd_DR_twmdx(i)=Jamp2dxd_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgb_twmg(P1,ANS)
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
      REAL*8 gb_DR_twmg
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
      Double Precision amp2gb_DR_twmg(maxamps), jamp2gb_DR_twmg(0:maxamps)
      common/to_ampsgb_DR_twmg/  amp2gb_DR_twmg,       jamp2gb_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgb_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gb_DR_twmg(0) /   2/          
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
              amp2gb_DR_twmg(ihel)=0d0
              jamp2gb_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gb_DR_twmg(0))
              jamp2gb_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gb_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=gb_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gb_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gb_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gb_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2gb_DR_twmg(maxamps), jamp2gb_DR_twmg(0:maxamps)
      common/to_ampsgb_DR_twmg/  amp2gb_DR_twmg,       jamp2gb_DR_twmg
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gb_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gb_DR_twmg =gb_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gb_DR_twmg(i)=amp2gb_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gb_DR_twmg(i)=Jamp2gb_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgd_twmg(P1,ANS)
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
      REAL*8 gd_DR_twmg
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
      Double Precision amp2gd_DR_twmg(maxamps), jamp2gd_DR_twmg(0:maxamps)
      common/to_ampsgd_DR_twmg/  amp2gd_DR_twmg,       jamp2gd_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgd_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gd_DR_twmg(0) /   2/          
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
              amp2gd_DR_twmg(ihel)=0d0
              jamp2gd_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gd_DR_twmg(0))
              jamp2gd_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gd_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=gd_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gd_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gd_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gd_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2gd_DR_twmg(maxamps), jamp2gd_DR_twmg(0:maxamps)
      common/to_ampsgd_DR_twmg/  amp2gd_DR_twmg,       jamp2gd_DR_twmg
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gd_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gd_DR_twmg =gd_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gd_DR_twmg(i)=amp2gd_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gd_DR_twmg(i)=Jamp2gd_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgg_twmbx(P1,ANS)
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
      REAL*8 gg_DR_twmbx
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
      Double Precision amp2gg_DR_twmbx(maxamps), jamp2gg_DR_twmbx(0:maxamps)
      common/to_ampsgg_DR_twmbx/  amp2gg_DR_twmbx,       jamp2gg_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_DR_twmbx(0) /   2/          
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
              amp2gg_DR_twmbx(ihel)=0d0
              jamp2gg_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_DR_twmbx(0))
              jamp2gg_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=gg_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gg_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_DR_twmbx(P,NHEL,IC)
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
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_DR_twmbx(maxamps), jamp2gg_DR_twmbx(0:maxamps)
      common/to_ampsgg_DR_twmbx/  amp2gg_DR_twmbx,       jamp2gg_DR_twmbx
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
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTB ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTB ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,BMASS   ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,BMASS   ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTB ,AMP(7   ))          
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      amp(2)=(0,0) !:
      amp(4)=(0,0) !:
      amp(7)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_DR_twmbx =gg_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_DR_twmbx(i)=amp2gg_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_DR_twmbx(i)=Jamp2gg_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgg_twmdx(P1,ANS)
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
      REAL*8 gg_DR_twmdx
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
      Double Precision amp2gg_DR_twmdx(maxamps), jamp2gg_DR_twmdx(0:maxamps)
      common/to_ampsgg_DR_twmdx/  amp2gg_DR_twmdx,       jamp2gg_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_DR_twmdx(0) /   2/          
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
              amp2gg_DR_twmdx(ihel)=0d0
              jamp2gg_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_DR_twmdx(0))
              jamp2gg_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=gg_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gg_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_DR_twmdx(P,NHEL,IC)
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
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_DR_twmdx(maxamps), jamp2gg_DR_twmdx(0:maxamps)
      common/to_ampsgg_DR_twmdx/  amp2gg_DR_twmdx,       jamp2gg_DR_twmdx
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
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTD ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTD ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTD ,AMP(7   ))          
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      amp(2)=(0,0) !:
      amp(4)=(0,0) !:
      amp(7)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_DR_twmdx =gg_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_DR_twmdx(i)=amp2gg_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_DR_twmdx(i)=Jamp2gg_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgg_twmsx(P1,ANS)
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
      REAL*8 gg_DR_twmsx
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
      Double Precision amp2gg_DR_twmsx(maxamps), jamp2gg_DR_twmsx(0:maxamps)
      common/to_ampsgg_DR_twmsx/  amp2gg_DR_twmsx,       jamp2gg_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgg_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gg_DR_twmsx(0) /   2/          
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
              amp2gg_DR_twmsx(ihel)=0d0
              jamp2gg_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gg_DR_twmsx(0))
              jamp2gg_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gg_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=gg_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gg_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gg_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gg_DR_twmsx(P,NHEL,IC)
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
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2gg_DR_twmsx(maxamps), jamp2gg_DR_twmsx(0:maxamps)
      common/to_ampsgg_DR_twmsx/  amp2gg_DR_twmsx,       jamp2gg_DR_twmsx
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
      CALL IXXXXX(P(0,5   ),ZERO ,NHEL(5   ),-1*IC(5   ),W(1,5   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,7   ),W(1,1   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,6   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL IOVXXX(W(1,5   ),W(1,8   ),W(1,4   ),GWFTS ,AMP(2   ))          
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,9   ))     
      CALL FVOXXX(W(1,9   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     10  ))                                                          
      CALL IOVXXX(W(1,5   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL FVOXXX(W(1,9   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,11  ))     
      CALL IOVXXX(W(1,5   ),W(1,11  ),W(1,4   ),GWFTS ,AMP(4   ))          
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     12  ))                                                          
      CALL FVIXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL IOVXXX(W(1,13  ),W(1,12  ),W(1,2   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,14  ))     
      CALL IOVXXX(W(1,14  ),W(1,12  ),W(1,1   ),GG ,AMP(6   ))             
      CALL JVVXXX(W(1,1   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,15  ))      
      CALL FVOXXX(W(1,3   ),W(1,15  ),GG ,TMASS   ,TWIDTH  ,W(1,16  ))     
      CALL IOVXXX(W(1,5   ),W(1,16  ),W(1,4   ),GWFTS ,AMP(7   ))          
      CALL IOVXXX(W(1,5   ),W(1,12  ),W(1,15  ),GG ,AMP(8   ))             
      amp(2)=(0,0) !:
      amp(4)=(0,0) !:
      amp(7)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   5)+AMP(   7)+AMP(   8)
      JAMP(   2) = -AMP(   3)-AMP(   4)-AMP(   6)-AMP(   7)-AMP(   8)
      gg_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gg_DR_twmsx =gg_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gg_DR_twmsx(i)=amp2gg_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gg_DR_twmsx(i)=Jamp2gg_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sgs_twmg(P1,ANS)
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
      REAL*8 gs_DR_twmg
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
      Double Precision amp2gs_DR_twmg(maxamps), jamp2gs_DR_twmg(0:maxamps)
      common/to_ampsgs_DR_twmg/  amp2gs_DR_twmg,       jamp2gs_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixgs_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2gs_DR_twmg(0) /   2/          
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
              amp2gs_DR_twmg(ihel)=0d0
              jamp2gs_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2gs_DR_twmg(0))
              jamp2gs_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=gs_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=gs_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2gs_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2gs_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION gs_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2gs_DR_twmg(maxamps), jamp2gs_DR_twmg(0:maxamps)
      common/to_ampsgs_DR_twmg/  amp2gs_DR_twmg,       jamp2gs_DR_twmg
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,1   ),G ,ZERO    ,ZERO    ,W(1,7   ))      
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,1   ),GG ,TMASS   ,TWIDTH  ,W(1,8   ))     
      CALL FVOXXX(W(1,8   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     9   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,9   ),W(1,5   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,6   ),W(1,10  ),W(1,1   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,6   ),W(1,8   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,11  ),W(1,7   ),GG ,AMP(5   ))             
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,1   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,2   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = -AMP(   1)+AMP(   2)+AMP(   4)-AMP(   5)+AMP(   6)
      JAMP(   2) = +AMP(   1)+AMP(   3)+AMP(   5)+AMP(   7)+AMP(   8)
      gs_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          gs_DR_twmg =gs_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2gs_DR_twmg(i)=amp2gs_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2gs_DR_twmg(i)=Jamp2gs_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssb_twmb(P1,ANS)
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
      REAL*8 sb_DR_twmb
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
      Double Precision amp2sb_DR_twmb(maxamps), jamp2sb_DR_twmb(0:maxamps)
      common/to_ampssb_DR_twmb/  amp2sb_DR_twmb,       jamp2sb_DR_twmb

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsb_DR_twmb/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sb_DR_twmb(0) /   1/          
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
              amp2sb_DR_twmb(ihel)=0d0
              jamp2sb_DR_twmb(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sb_DR_twmb(0))
              jamp2sb_DR_twmb(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sb_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              T=sb_DR_twmb(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2sb_DR_twmb(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sb_DR_twmb(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sb_DR_twmb(P,NHEL,IC)
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
      Double Precision amp2sb_DR_twmb(maxamps), jamp2sb_DR_twmb(0:maxamps)
      common/to_ampssb_DR_twmb/  amp2sb_DR_twmb,       jamp2sb_DR_twmb
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      sb_DR_twmb = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sb_DR_twmb =sb_DR_twmb+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sb_DR_twmb(i)=amp2sb_DR_twmb(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sb_DR_twmb(i)=Jamp2sb_DR_twmb(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssbx_twmbx(P1,ANS)
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
      REAL*8 sbx_DR_twmbx
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
      Double Precision amp2sbx_DR_twmbx(maxamps), jamp2sbx_DR_twmbx(0:maxamps)
      common/to_ampssbx_DR_twmbx/  amp2sbx_DR_twmbx,       jamp2sbx_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsbx_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sbx_DR_twmbx(0) /   1/          
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
              amp2sbx_DR_twmbx(ihel)=0d0
              jamp2sbx_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sbx_DR_twmbx(0))
              jamp2sbx_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=sbx_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2sbx_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sbx_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sbx_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2sbx_DR_twmbx(maxamps), jamp2sbx_DR_twmbx(0:maxamps)
      common/to_ampssbx_DR_twmbx/  amp2sbx_DR_twmbx,       jamp2sbx_DR_twmbx
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      sbx_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sbx_DR_twmbx =sbx_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sbx_DR_twmbx(i)=amp2sbx_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sbx_DR_twmbx(i)=Jamp2sbx_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssg_twmg(P1,ANS)
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
      REAL*8 sg_DR_twmg
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
      Double Precision amp2sg_DR_twmg(maxamps), jamp2sg_DR_twmg(0:maxamps)
      common/to_ampssg_DR_twmg/  amp2sg_DR_twmg,       jamp2sg_DR_twmg

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsg_DR_twmg/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    8/          
      DATA jamp2sg_DR_twmg(0) /   2/          
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
              amp2sg_DR_twmg(ihel)=0d0
              jamp2sg_DR_twmg(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sg_DR_twmg(0))
              jamp2sg_DR_twmg(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              T=sg_DR_twmg(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2sg_DR_twmg(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sg_DR_twmg(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sg_DR_twmg(P,NHEL,IC)
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
      Double Precision amp2sg_DR_twmg(maxamps), jamp2sg_DR_twmg(0:maxamps)
      common/to_ampssg_DR_twmg/  amp2sg_DR_twmg,       jamp2sg_DR_twmg
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
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,6   ))     
      CALL FVOXXX(W(1,6   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     7   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,7   ),W(1,5   ),GG ,AMP(1   ))             
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     8   ))                                                          
      CALL JVVXXX(W(1,5   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,9   ))      
      CALL IOVXXX(W(1,8   ),W(1,3   ),W(1,9   ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,5   ),GG ,TMASS   ,TWIDTH  ,W(1,10  ))     
      CALL IOVXXX(W(1,8   ),W(1,10  ),W(1,2   ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,8   ),W(1,6   ),W(1,5   ),GG ,AMP(4   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     11  ))                                                          
      CALL FVIXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,12  ))     
      CALL IOVXXX(W(1,12  ),W(1,11  ),W(1,2   ),GG ,AMP(5   ))             
      CALL IOVXXX(W(1,1   ),W(1,11  ),W(1,9   ),GG ,AMP(6   ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,13  ))     
      CALL FVIXXX(W(1,13  ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     14  ))                                                          
      CALL IOVXXX(W(1,14  ),W(1,3   ),W(1,5   ),GG ,AMP(7   ))             
      CALL IOVXXX(W(1,13  ),W(1,11  ),W(1,5   ),GG ,AMP(8   ))             
      JAMP(   1) = +AMP(   1)-AMP(   2)+AMP(   4)+AMP(   5)-AMP(   6)
      JAMP(   2) = +AMP(   2)+AMP(   3)+AMP(   6)+AMP(   7)+AMP(   8)
      sg_DR_twmg = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sg_DR_twmg =sg_DR_twmg+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sg_DR_twmg(i)=amp2sg_DR_twmg(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sg_DR_twmg(i)=Jamp2sg_DR_twmg(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sss_twms(P1,ANS)
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
      REAL*8 ss_DR_twms
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
      Double Precision amp2ss_DR_twms(maxamps), jamp2ss_DR_twms(0:maxamps)
      common/to_ampsss_DR_twms/  amp2ss_DR_twms,       jamp2ss_DR_twms

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixss_DR_twms/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ss_DR_twms(0) /   2/          
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
              amp2ss_DR_twms(ihel)=0d0
              jamp2ss_DR_twms(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ss_DR_twms(0))
              jamp2ss_DR_twms(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ss_DR_twms(P ,NHEL(1,IHEL),JC(1))            
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
              T=ss_DR_twms(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2ss_DR_twms(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ss_DR_twms(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ss_DR_twms(P,NHEL,IC)
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
      Double Precision amp2ss_DR_twms(maxamps), jamp2ss_DR_twms(0:maxamps)
      common/to_ampsss_DR_twms/  amp2ss_DR_twms,       jamp2ss_DR_twms
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     9   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL IOVXXX(W(1,9   ),W(1,3   ),W(1,10  ),GG ,AMP(3   ))             
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,10  ),GG ,AMP(4   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      ss_DR_twms = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ss_DR_twms =ss_DR_twms+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ss_DR_twms(i)=amp2ss_DR_twms(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ss_DR_twms(i)=Jamp2ss_DR_twms(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sssx_twmsx(P1,ANS)
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
      REAL*8 ssx_DR_twmsx
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
      Double Precision amp2ssx_DR_twmsx(maxamps), jamp2ssx_DR_twmsx(0:maxamps)
      common/to_ampsssx_DR_twmsx/  amp2ssx_DR_twmsx,       jamp2ssx_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixssx_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2ssx_DR_twmsx(0) /   2/          
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
              amp2ssx_DR_twmsx(ihel)=0d0
              jamp2ssx_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ssx_DR_twmsx(0))
              jamp2ssx_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ssx_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=ssx_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2ssx_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ssx_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ssx_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2ssx_DR_twmsx(maxamps), jamp2ssx_DR_twmsx(0:maxamps)
      common/to_ampsssx_DR_twmsx/  amp2ssx_DR_twmsx,       jamp2ssx_DR_twmsx
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      JAMP(   2) = -AMP(   3)-AMP(   4)
      ssx_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ssx_DR_twmsx =ssx_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ssx_DR_twmsx(i)=amp2ssx_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ssx_DR_twmsx(i)=Jamp2ssx_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssu_twmu(P1,ANS)
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
      REAL*8 su_DR_twmu
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
      Double Precision amp2su_DR_twmu(maxamps), jamp2su_DR_twmu(0:maxamps)
      common/to_ampssu_DR_twmu/  amp2su_DR_twmu,       jamp2su_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsu_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2su_DR_twmu(0) /   1/          
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
              amp2su_DR_twmu(ihel)=0d0
              jamp2su_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2su_DR_twmu(0))
              jamp2su_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=su_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=su_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2su_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2su_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION su_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2su_DR_twmu(maxamps), jamp2su_DR_twmu(0:maxamps)
      common/to_ampssu_DR_twmu/  amp2su_DR_twmu,       jamp2su_DR_twmu
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,2   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      su_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          su_DR_twmu =su_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2su_DR_twmu(i)=amp2su_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2su_DR_twmu(i)=Jamp2su_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssux_twmux(P1,ANS)
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
      REAL*8 sux_DR_twmux
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
      Double Precision amp2sux_DR_twmux(maxamps), jamp2sux_DR_twmux(0:maxamps)
      common/to_ampssux_DR_twmux/  amp2sux_DR_twmux,       jamp2sux_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsux_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2sux_DR_twmux(0) /   1/          
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
              amp2sux_DR_twmux(ihel)=0d0
              jamp2sux_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sux_DR_twmux(0))
              jamp2sux_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=sux_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2sux_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sux_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sux_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2sux_DR_twmux(maxamps), jamp2sux_DR_twmux(0:maxamps)
      common/to_ampssux_DR_twmux/  amp2sux_DR_twmux,       jamp2sux_DR_twmux
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
      CALL FVIXXX(W(1,1   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,1   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      sux_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sux_DR_twmux =sux_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sux_DR_twmux(i)=amp2sux_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sux_DR_twmux(i)=Jamp2sux_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Ssxs_twmsx(P1,ANS)
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
      REAL*8 sxs_DR_twmsx
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
      Double Precision amp2sxs_DR_twmsx(maxamps), jamp2sxs_DR_twmsx(0:maxamps)
      common/to_ampssxs_DR_twmsx/  amp2sxs_DR_twmsx,       jamp2sxs_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixsxs_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    4/          
      DATA jamp2sxs_DR_twmsx(0) /   2/          
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
              amp2sxs_DR_twmsx(ihel)=0d0
              jamp2sxs_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2sxs_DR_twmsx(0))
              jamp2sxs_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=sxs_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=sxs_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2sxs_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2sxs_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION sxs_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2sxs_DR_twmsx(maxamps), jamp2sxs_DR_twmsx(0:maxamps)
      common/to_ampssxs_DR_twmsx/  amp2sxs_DR_twmsx,       jamp2sxs_DR_twmsx
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
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
      amp(3)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      JAMP(   2) = +AMP(   3)+AMP(   4)
      sxs_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          sxs_DR_twmsx =sxs_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2sxs_DR_twmsx(i)=amp2sxs_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2sxs_DR_twmsx(i)=Jamp2sxs_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sub_twmu(P1,ANS)
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
      REAL*8 ub_DR_twmu
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
      Double Precision amp2ub_DR_twmu(maxamps), jamp2ub_DR_twmu(0:maxamps)
      common/to_ampsub_DR_twmu/  amp2ub_DR_twmu,       jamp2ub_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixub_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2ub_DR_twmu(0) /   1/          
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
              amp2ub_DR_twmu(ihel)=0d0
              jamp2ub_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ub_DR_twmu(0))
              jamp2ub_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ub_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=ub_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2ub_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ub_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ub_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2ub_DR_twmu(maxamps), jamp2ub_DR_twmu(0:maxamps)
      common/to_ampsub_DR_twmu/  amp2ub_DR_twmu,       jamp2ub_DR_twmu
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      ub_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ub_DR_twmu =ub_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ub_DR_twmu(i)=amp2ub_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ub_DR_twmu(i)=Jamp2ub_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sud_twmu(P1,ANS)
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
      REAL*8 ud_DR_twmu
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
      Double Precision amp2ud_DR_twmu(maxamps), jamp2ud_DR_twmu(0:maxamps)
      common/to_ampsud_DR_twmu/  amp2ud_DR_twmu,       jamp2ud_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixud_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2ud_DR_twmu(0) /   1/          
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
              amp2ud_DR_twmu(ihel)=0d0
              jamp2ud_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2ud_DR_twmu(0))
              jamp2ud_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=ud_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=ud_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2ud_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2ud_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION ud_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2ud_DR_twmu(maxamps), jamp2ud_DR_twmu(0:maxamps)
      common/to_ampsud_DR_twmu/  amp2ud_DR_twmu,       jamp2ud_DR_twmu
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      ud_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          ud_DR_twmu =ud_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2ud_DR_twmu(i)=amp2ud_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2ud_DR_twmu(i)=Jamp2ud_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Sus_twmu(P1,ANS)
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
      REAL*8 us_DR_twmu
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
      Double Precision amp2us_DR_twmu(maxamps), jamp2us_DR_twmu(0:maxamps)
      common/to_ampsus_DR_twmu/  amp2us_DR_twmu,       jamp2us_DR_twmu

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixus_DR_twmu/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2us_DR_twmu(0) /   1/          
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
              amp2us_DR_twmu(ihel)=0d0
              jamp2us_DR_twmu(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2us_DR_twmu(0))
              jamp2us_DR_twmu(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=us_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              T=us_DR_twmu(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2us_DR_twmu(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2us_DR_twmu(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION us_DR_twmu(P,NHEL,IC)
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
      Double Precision amp2us_DR_twmu(maxamps), jamp2us_DR_twmu(0:maxamps)
      common/to_ampsus_DR_twmu/  amp2us_DR_twmu,       jamp2us_DR_twmu
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      us_DR_twmu = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          us_DR_twmu =us_DR_twmu+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2us_DR_twmu(i)=amp2us_DR_twmu(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2us_DR_twmu(i)=Jamp2us_DR_twmu(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suux_twmbx(P1,ANS)
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
      REAL*8 uux_DR_twmbx
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
      Double Precision amp2uux_DR_twmbx(maxamps), jamp2uux_DR_twmbx(0:maxamps)
      common/to_ampsuux_DR_twmbx/  amp2uux_DR_twmbx,       jamp2uux_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_DR_twmbx(0) /   1/          
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
              amp2uux_DR_twmbx(ihel)=0d0
              jamp2uux_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_DR_twmbx(0))
              jamp2uux_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uux_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uux_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2uux_DR_twmbx(maxamps), jamp2uux_DR_twmbx(0:maxamps)
      common/to_ampsuux_DR_twmbx/  amp2uux_DR_twmbx,       jamp2uux_DR_twmbx
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
      amp(1)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_DR_twmbx =uux_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_DR_twmbx(i)=amp2uux_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_DR_twmbx(i)=Jamp2uux_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suux_twmdx(P1,ANS)
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
      REAL*8 uux_DR_twmdx
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
      Double Precision amp2uux_DR_twmdx(maxamps), jamp2uux_DR_twmdx(0:maxamps)
      common/to_ampsuux_DR_twmdx/  amp2uux_DR_twmdx,       jamp2uux_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_DR_twmdx(0) /   1/          
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
              amp2uux_DR_twmdx(ihel)=0d0
              jamp2uux_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_DR_twmdx(0))
              jamp2uux_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uux_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uux_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2uux_DR_twmdx(maxamps), jamp2uux_DR_twmdx(0:maxamps)
      common/to_ampsuux_DR_twmdx/  amp2uux_DR_twmdx,       jamp2uux_DR_twmdx
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
      amp(1)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_DR_twmdx =uux_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_DR_twmdx(i)=amp2uux_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_DR_twmdx(i)=Jamp2uux_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suux_twmsx(P1,ANS)
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
      REAL*8 uux_DR_twmsx
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
      Double Precision amp2uux_DR_twmsx(maxamps), jamp2uux_DR_twmsx(0:maxamps)
      common/to_ampsuux_DR_twmsx/  amp2uux_DR_twmsx,       jamp2uux_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuux_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uux_DR_twmsx(0) /   1/          
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
              amp2uux_DR_twmsx(ihel)=0d0
              jamp2uux_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uux_DR_twmsx(0))
              jamp2uux_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uux_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uux_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uux_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uux_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uux_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2uux_DR_twmsx(maxamps), jamp2uux_DR_twmsx(0:maxamps)
      common/to_ampsuux_DR_twmsx/  amp2uux_DR_twmsx,       jamp2uux_DR_twmsx
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
      amp(1)=(0,0) !:
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uux_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uux_DR_twmsx =uux_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uux_DR_twmsx(i)=amp2uux_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uux_DR_twmsx(i)=Jamp2uux_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxb_twmux(P1,ANS)
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
      REAL*8 uxb_DR_twmux
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
      Double Precision amp2uxb_DR_twmux(maxamps), jamp2uxb_DR_twmux(0:maxamps)
      common/to_ampsuxb_DR_twmux/  amp2uxb_DR_twmux,       jamp2uxb_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxb_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxb_DR_twmux(0) /   1/          
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
              amp2uxb_DR_twmux(ihel)=0d0
              jamp2uxb_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxb_DR_twmux(0))
              jamp2uxb_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxb_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxb_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxb_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxb_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxb_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2uxb_DR_twmux(maxamps), jamp2uxb_DR_twmux(0:maxamps)
      common/to_ampsuxb_DR_twmux/  amp2uxb_DR_twmux,       jamp2uxb_DR_twmux
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTB ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTB ,BMASS   ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxb_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxb_DR_twmux =uxb_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxb_DR_twmux(i)=amp2uxb_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxb_DR_twmux(i)=Jamp2uxb_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxd_twmux(P1,ANS)
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
      REAL*8 uxd_DR_twmux
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
      Double Precision amp2uxd_DR_twmux(maxamps), jamp2uxd_DR_twmux(0:maxamps)
      common/to_ampsuxd_DR_twmux/  amp2uxd_DR_twmux,       jamp2uxd_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxd_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxd_DR_twmux(0) /   1/          
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
              amp2uxd_DR_twmux(ihel)=0d0
              jamp2uxd_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxd_DR_twmux(0))
              jamp2uxd_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxd_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxd_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxd_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxd_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxd_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2uxd_DR_twmux(maxamps), jamp2uxd_DR_twmux(0:maxamps)
      common/to_ampsuxd_DR_twmux/  amp2uxd_DR_twmux,       jamp2uxd_DR_twmux
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTD ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTD ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxd_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxd_DR_twmux =uxd_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxd_DR_twmux(i)=amp2uxd_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxd_DR_twmux(i)=Jamp2uxd_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxs_twmux(P1,ANS)
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
      REAL*8 uxs_DR_twmux
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
      Double Precision amp2uxs_DR_twmux(maxamps), jamp2uxs_DR_twmux(0:maxamps)
      common/to_ampsuxs_DR_twmux/  amp2uxs_DR_twmux,       jamp2uxs_DR_twmux

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxs_DR_twmux/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxs_DR_twmux(0) /   1/          
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
              amp2uxs_DR_twmux(ihel)=0d0
              jamp2uxs_DR_twmux(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxs_DR_twmux(0))
              jamp2uxs_DR_twmux(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxs_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxs_DR_twmux(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxs_DR_twmux(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxs_DR_twmux(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxs_DR_twmux(P,NHEL,IC)
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
      Double Precision amp2uxs_DR_twmux(maxamps), jamp2uxs_DR_twmux(0:maxamps)
      common/to_ampsuxs_DR_twmux/  amp2uxs_DR_twmux,       jamp2uxs_DR_twmux
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
      CALL FVIXXX(W(1,2   ),W(1,4   ),GWFTS ,TMASS   ,TWIDTH  ,W(1,        
     &     6   ))                                                          
      CALL JIOXXX(W(1,5   ),W(1,1   ),GG ,ZERO    ,ZERO    ,W(1,7   ))     
      CALL IOVXXX(W(1,6   ),W(1,3   ),W(1,7   ),GG ,AMP(1   ))             
      CALL FVOXXX(W(1,3   ),W(1,4   ),GWFTS ,ZERO    ,ZERO    ,W(1,        
     &     8   ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,8   ),W(1,7   ),GG ,AMP(2   ))             
      JAMP(   1) = -AMP(   1)-AMP(   2)
      uxs_DR_twmux = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxs_DR_twmux =uxs_DR_twmux+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxs_DR_twmux(i)=amp2uxs_DR_twmux(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxs_DR_twmux(i)=Jamp2uxs_DR_twmux(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxu_twmbx(P1,ANS)
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
      REAL*8 uxu_DR_twmbx
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
      Double Precision amp2uxu_DR_twmbx(maxamps), jamp2uxu_DR_twmbx(0:maxamps)
      common/to_ampsuxu_DR_twmbx/  amp2uxu_DR_twmbx,       jamp2uxu_DR_twmbx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_DR_twmbx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_DR_twmbx(0) /   1/          
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
              amp2uxu_DR_twmbx(ihel)=0d0
              jamp2uxu_DR_twmbx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_DR_twmbx(0))
              jamp2uxu_DR_twmbx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxu_DR_twmbx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxu_DR_twmbx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_DR_twmbx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_DR_twmbx(P,NHEL,IC)
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
      Double Precision amp2uxu_DR_twmbx(maxamps), jamp2uxu_DR_twmbx(0:maxamps)
      common/to_ampsuxu_DR_twmbx/  amp2uxu_DR_twmbx,       jamp2uxu_DR_twmbx
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
      amp(1)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_DR_twmbx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_DR_twmbx =uxu_DR_twmbx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_DR_twmbx(i)=amp2uxu_DR_twmbx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_DR_twmbx(i)=Jamp2uxu_DR_twmbx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxu_twmdx(P1,ANS)
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
      REAL*8 uxu_DR_twmdx
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
      Double Precision amp2uxu_DR_twmdx(maxamps), jamp2uxu_DR_twmdx(0:maxamps)
      common/to_ampsuxu_DR_twmdx/  amp2uxu_DR_twmdx,       jamp2uxu_DR_twmdx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_DR_twmdx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_DR_twmdx(0) /   1/          
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
              amp2uxu_DR_twmdx(ihel)=0d0
              jamp2uxu_DR_twmdx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_DR_twmdx(0))
              jamp2uxu_DR_twmdx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxu_DR_twmdx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxu_DR_twmdx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_DR_twmdx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_DR_twmdx(P,NHEL,IC)
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
      Double Precision amp2uxu_DR_twmdx(maxamps), jamp2uxu_DR_twmdx(0:maxamps)
      common/to_ampsuxu_DR_twmdx/  amp2uxu_DR_twmdx,       jamp2uxu_DR_twmdx
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
      amp(1)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_DR_twmdx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_DR_twmdx =uxu_DR_twmdx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_DR_twmdx(i)=amp2uxu_DR_twmdx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_DR_twmdx(i)=Jamp2uxu_DR_twmdx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
      SUBROUTINE DR_Suxu_twmsx(P1,ANS)
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
      REAL*8 uxu_DR_twmsx
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
      Double Precision amp2uxu_DR_twmsx(maxamps), jamp2uxu_DR_twmsx(0:maxamps)
      common/to_ampsuxu_DR_twmsx/  amp2uxu_DR_twmsx,       jamp2uxu_DR_twmsx

      character*79         hel_buff
      common/to_helicity/  hel_buff

      REAL*8 POL(2)
      common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      common/to_matrixuxu_DR_twmsx/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.false./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    2/          
      DATA jamp2uxu_DR_twmsx(0) /   1/          
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
              amp2uxu_DR_twmsx(ihel)=0d0
              jamp2uxu_DR_twmsx(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2uxu_DR_twmsx(0))
              jamp2uxu_DR_twmsx(ihel)=0d0
          ENDDO
      ENDIF
      ANS(IPROC) = 0D0
      write(hel_buff,'(16i5)') (0,i=1,nexternal)
      IF (ISUM_HEL .EQ. 0 .OR. NTRY .LT. 10) THEN
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=uxu_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              T=uxu_DR_twmsx(P ,NHEL(1,IHEL),JC(1))            
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
              XTOT=XTOT+AMP2uxu_DR_twmsx(MAPCONFIG(IHEL))
          ENDDO
          IF (XTOT.NE.0D0) THEN
              ANS(IPROC)=ANS(IPROC)*AMP2uxu_DR_twmsx(MAPCONFIG(ICONFIG))/XTOT
          ELSE
              ANS(IPROC)=0D0
          ENDIF
      ENDIF
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION uxu_DR_twmsx(P,NHEL,IC)
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
      Double Precision amp2uxu_DR_twmsx(maxamps), jamp2uxu_DR_twmsx(0:maxamps)
      common/to_ampsuxu_DR_twmsx/  amp2uxu_DR_twmsx,       jamp2uxu_DR_twmsx
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
      amp(1)=(0,0) !:
      JAMP(   1) = +AMP(   1)+AMP(   2)
      uxu_DR_twmsx = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          uxu_DR_twmsx =uxu_DR_twmsx+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2uxu_DR_twmsx(i)=amp2uxu_DR_twmsx(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2uxu_DR_twmsx(i)=Jamp2uxu_DR_twmsx(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       
