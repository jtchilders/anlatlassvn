      double precision ptmin_def, ymax_def, rll_min, rll_max, ymax_part
      double precision jetgenktp
      integer min_jets, njet_min 	
      double precision rjjmin, rjl_min, ymax_l, ptlmin,ptlmax, etajjmin
      double precision dely_jlep,dely_jveto
      double precision mdijmin,delmll,mdijmax
      double precision ptmin_def_default(4) /20,20,20,20/

      logical ysign, lrapidgap,grapidgap,lveto
      double precision ptmin_veto,ymax_veto
      double precision rgg_min, rgg_max, rjg_min, rlg_min, ymax_g
      double precision ptgmin,ptgmax, dely_jphot
      double precision isol_phopart, eff_phopart,ptinmin
      double precision mllmin,mllmax,mlAmin, mlAmax, MTmax,MTmin

      common/jetdef/rjjmin,ymax_part,njet_min
      common/jetalgo/jetgenktp	
      
      common/basiccuts/ptmin_def(4), ymax_def, min_jets

      common/photoncuts/isol_phopart, eff_phopart

      common/vbfcuts/rgg_min,rgg_max,rjg_min,rlg_min,ymax_g,ptgmin,dely_jphot,
     &	             rll_min,rll_max,rjl_min,ymax_l,ptlmax,ptlmin,etajjmin,
     &               mdijmin,mdijmax,ymax_veto,ptmin_veto,delmll,dely_jlep,
     &               dely_jveto,grapidgap,lveto,ysign,lrapidgap,ptgmax

      common/tribosoncuts/mllmin,mllmax,mlAmin,mlAmax,MTmax,MTmin,ptinmin

