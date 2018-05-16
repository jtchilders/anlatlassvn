c -*- Fortran -*-

c flg_nlotest: perform calls to outfun for testing NLO output
c flg_withsubtr: counterterms are included in NLO test
c flg_withdamp: Born zero procedure
c flg_withreg: if regular regions exist or not
c flg_smartsig: remember or not equal suqred amplitude
c flg_bornonly: do the Born contribution only
      logical flg_fastbtlbound
      common/pwhg_flg_add/flg_fastbtlbound

