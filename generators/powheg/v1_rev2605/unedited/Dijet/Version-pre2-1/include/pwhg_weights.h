      integer weights_max,weights_num
      parameter (weights_max=20)
      real * 8 weights_val(weights_max),
     1         weights_renfac(weights_max),
     2         weights_facfac(weights_max)
      integer weights_npdf1(weights_max),
     1        weights_npdf2(weights_max)
      character * 3 weights_whichpdf(weights_max)
      common/pwhg_weights/weights_val,weights_renfac,weights_facfac,
     1    weights_npdf1,weights_npdf2,weights_num,weights_whichpdf
