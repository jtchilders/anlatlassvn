********************************************************************************
*** Storing the values of the anomalous vector boson couplings               ***
********************************************************************************

** Parameterisation flag, universal formfactor, switch for individual/universal
** formfactor
      integer trianom, ffexponent
      double precision ffmassscale2
      logical formfact, formfac_ind

      common /anomFF/ ffmassscale2, trianom, ffexponent
      common /anomFF/ formfact, formfac_ind


** Trianom = 1 parameterisation (_0 = initial values, before formfactors)
      double precision fwww_0, fw_0, fb_0
      double precision fwww, fw, fb
** Trianom = 1 formfactors
      double precision massscale2FWWW, massscale2FW, massscale2FB
      integer ffexpFWWW, ffexpFW, ffexpFB

      common /trianom1/ fwww_0, fw_0, fb_0, fwww, fw, fb
      common /trianom1/ massscale2FWWW, massscale2FW, massscale2FB
      common /trianom1/ ffexpFWWW, ffexpFW, ffexpFB

** Trianom = 2 parameterisation (_0 = initial values, before formfactors)
      double precision lambda0_0, zDg0_0, aDkappa0_0, zDkappa0_0
      double precision lambda0, zDg0, aDkappa0, zDkappa0
** Trianom = 2 formfactors
      double precision massscale2KZ, massscale2KA, massscale2L, massscale2G
      integer ffexpKZ, ffexpKA, ffexpL, ffexpG

      common /trianom2/ lambda0_0, zDg0_0, aDkappa0_0, zDkappa0_0
      common /trianom2/ lambda0, zDg0, aDkappa0, zDkappa0
      common /trianom2/ massscale2KZ, massscale2KA, massscale2L, massscale2G
      common /trianom2/ ffexpKZ, ffexpKA, ffexpL, ffexpG


** Dimension-8 operators (for triboson processes only)
      double precision ft0_0, ft1_0, ft2_0, ft5_0, ft6_0, ft7_0
      double precision fm0_0, fm1_0, fm2_0, fm3_0, fm4_0, fm5_0
      double precision fm6_0, fm7_0, fs0_0, fs1_0
      double precision ft0, ft1, ft2, ft5, ft6, ft7
      double precision fm0, fm1, fm2, fm3, fm4, fm5
      double precision fm6, fm7, fs0, fs1

      common /anomDim8/ ft0_0, ft1_0, ft2_0, ft5_0, ft6_0, ft7_0
      common /anomDim8/ fm0_0, fm1_0, fm2_0, fm3_0, fm4_0, fm5_0
      common /anomDim8/ fm6_0, fm7_0, fs0_0, fs1_0
      common /anomDim8/ ft0, ft1, ft2, ft5, ft6, ft7 
      common /anomDim8/ fm0, fm1, fm2, fm3, fm4, fm5
      common /anomDim8/ fm6, fm7, fs0, fs1


** Dimension-6 operators used only for WWjj production
      double precision fbw_0, fdw_0, fww_0, fbb_0
      double precision fbwtilde_0, fdwtilde_0, fwwtilde_0, fbbtilde_0
      double precision fwtilde_0, fbtilde_0, fwwwtilde_0
      double precision fbw, fdw, fww, fbb
      double precision fbwtilde, fdwtilde, fwwtilde, fbbtilde
      double precision fwtilde, fbtilde, fwwwtilde
      logical ofs, processAnom

      common /anomDim6/ fbw_0, fdw_0, fww_0, fbb_0
      common /anomDim6/ fbwtilde_0, fdwtilde_0, fwwtilde_0, fbbtilde_0
      common /anomDim6/ fwtilde_0, fbtilde_0, fwwwtilde_0
      common /anomDim6/ fbw, fdw, fww, fbb
      common /anomDim6/ fbwtilde, fdwtilde, fwwtilde, fbbtilde
      common /anomDim6/ fwtilde, fbtilde, fwwwtilde, ofs, processAnom

