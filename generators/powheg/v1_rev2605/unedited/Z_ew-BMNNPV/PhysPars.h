c -*- Fortran -*-
      real * 8
     $     ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth,ph_Hmass,ph_Hwidth,
     $     ph_cthw,ph_sthw,ph_sthw2,
     $     ph_Zmass2,ph_Zmass2low,ph_Zmass2high,
     $     ph_Wmass2,ph_Wmass2low,ph_Wmass2high,ph_ZmZw,ph_WmWw,
     $     ph_Hmass2,ph_HmHw,ph_mTop,ph_mBot,ph_mCha,ph_mStr,
     $     ph_mUp,ph_mDown,ph_mEl,ph_mMuon,ph_mTau,
     $     ph_alphaem,ph_unit_e,ph_CKM(3,3),ph_gmu,ph_alpha_mz
      
      common/ph_common/
     $     ph_Zmass,ph_Zwidth,ph_Wmass,ph_Wwidth,ph_Hmass,ph_Hwidth,
     $     ph_cthw,ph_sthw,ph_sthw2,
     $     ph_Zmass2,ph_Zmass2low,ph_Zmass2high,
     $     ph_Wmass2,ph_Wmass2low,ph_Wmass2high,ph_ZmZw,ph_WmWw,
     $     ph_Hmass2,ph_HmHw,ph_mTop,ph_mBot,ph_mCha,ph_mStr,
     $     ph_mUp,ph_mDown,ph_mEl,ph_mMuon,ph_mTau,
     $     ph_alphaem,ph_unit_e,ph_CKM,ph_gmu,ph_alpha_mz
      
	  logical complexmasses
      common/complexm/complexmasses

	  logical iftopinloop
      common/topinloop/iftopinloop

      integer scheme
      common/sch/scheme
