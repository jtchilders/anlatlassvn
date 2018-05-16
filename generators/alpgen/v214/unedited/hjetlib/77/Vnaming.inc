C *** MODULE naming ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      integer nmax
      parameter (nmax = 10)
      integer nlor
      parameter (nlor = 6)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
      integer gluon, wm, wp, z, photon, qdbr, qbtbr, qupbrl, qdbrl, ebrl
     1   , nuebr, qupbrr, qdbrr, qcbr, qcbrl, qcbrr, ebrr, qsbrl, qsbrr
     2   , qtopbr, qupbr, qup, qd, qbt, qupl, qdl, qupr, qdr, el, nue, 
     3   qc, qcr, er, qsl, qsr, qtop, qcl, xpp, ypp, zzww, higgs, xglu, 
     4   yglu
      parameter (gluon = 1, wm = 4, wp = 3, z = 2, photon = 5, qdbr = 6
     1   , qbtbr = 7, qupbrl = 8, qdbrl = 9, ebrl = 10, nuebr = 11, 
     2   qupbrr = 12, qdbrr = 13, qcbr = 14, qcbrl = 15, qcbrr = 16, 
     3   ebrr = 17, qsbrl = 18, qsbrr = 19, qtopbr = 20, qupbr = 21, qup
     4    = 22, qd = 23, qbt = 24, qupl = 25, qdl = 26, qupr = 27, qdr
     5    = 28, el = 29, nue = 30, qc = 31, qcr = 32, er = 33, qsl = 34
     6   , qsr = 35, qtop = 36, qcl = 37, xpp = 38, ypp = 39, zzww = 40
     7   , higgs = 41, xglu = 42, yglu = 43)
      integer prtcl__chcg(43)
      integer prtcl__cnschg(43)
      integer prtcl__lortyp(43)
      integer prtcl__nlordof(43)
      integer prtcl__colrep(43)
      integer prtcl__zcoup(43)
      integer chcg(43), cnschg(43)!chcg(j) is the charged conjugate of j
      integer glu, frm, frmbr, msgbs, frmhl, frmbrhl, frmhr, frmbrhr, 
     1   axww, axzzww, pht, scal, xgl
      parameter (glu = 1, frm = 3, frmbr = 4, msgbs = 2, frmhl = 5, 
     1   frmbrhl = 6, frmhr = 7, frmbrhr = 8, axww = 9, axzzww = 10, pht
     2    = 11, scal = 12, xgl = 13)
      integer up, dn, lpt
      parameter (up = 1, dn = 2, lpt = 3)
      integer conv(43)
      common /naming/prtcl__chcg, prtcl__cnschg, prtcl__lortyp, 
     1   prtcl__nlordof, prtcl__colrep, prtcl__zcoup, chcg, cnschg, conv
C *** END MODULE naming ***
