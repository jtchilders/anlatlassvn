c -*- Fortran -*-

c flg_nlotest: perform calls to outfun for testing NLO output
c flg_withsubtr: counterterms are included in NLO test
c flg_withdamp: Born zero procedure
c flg_withreg: if regular regions exist or not
c flg_smartsig: remember or not equal squared amplitude
c flg_bornonly: do the Born contribution only

c flg_lightpart_check: in genflavreglist perform or not perform the check
c                     that there are no coloured light partons before flst_lightpart

c flg_debug: outputs extra infos for radiation region and randon numbers on the LHEF   
c flg_withnegweights: allows negative weighted events
c flg_jacsing:  importance sampling for singular jacobians in case of massless FSR
c flg_weightedev:  allows weighted events
c flg_pdfreweight: outputs extra infos useful for pdf's reweighting on the LHEF
c flg_collremnsamp:  importance sampling for collinear remnants
c flg_reweight: outputs extra infos for reweighting LH events

      logical flg_nlotest,flg_withsubtr,flg_withdamp,flg_withreg,
     1     flg_smartsig,flg_bornonly,flg_debug,flg_withnegweights,
     2     flg_jacsing,flg_weightedev,flg_pdfreweight,flg_collremnsamp,
     3     flg_lightpart_check,flg_btlscalereal,flg_btlscalect,
     4     flg_bornzerodamp,flg_ckkwscalup,
     5     flg_minlo,flg_reweight,flg_newweight,flg_fastbtlbound,
     6     flg_storemintupb,flg_doublefsr
      character * 1 flg_btildepart
      character * 20 flg_processid
      common/pwhg_flg/flg_nlotest,flg_withsubtr,flg_withdamp,
     2     flg_withreg,flg_smartsig,flg_bornonly,flg_debug,
     3     flg_withnegweights,flg_jacsing,flg_weightedev,
     4     flg_pdfreweight,flg_collremnsamp,flg_lightpart_check,
     5     flg_btlscalereal,flg_btlscalect,
     6     flg_bornzerodamp,flg_ckkwscalup,
     7     flg_minlo,flg_reweight,flg_newweight,flg_fastbtlbound,
     8     flg_storemintupb,flg_doublefsr,
     9     flg_btildepart,flg_processid
      save /pwhg_flg/
