C *** MODULE prcint ***
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
      common /naming/prtcl__chcg, prtcl__cnschg, prtcl__lortyp, 
     1   prtcl__nlordof, prtcl__colrep, prtcl__zcoup, chcg, cnschg, conv
      integer gluon
      parameter (gluon = 1)
      integer wm
      parameter (wm = 4)
      integer wp
      parameter (wp = 3)
      integer z
      parameter (z = 2)
      integer photon
      parameter (photon = 5)
      integer qdbr
      parameter (qdbr = 6)
      integer qbtbr
      parameter (qbtbr = 7)
      integer qupbrl
      parameter (qupbrl = 8)
      integer qdbrl
      parameter (qdbrl = 9)
      integer ebrl
      parameter (ebrl = 10)
      integer nuebr
      parameter (nuebr = 11)
      integer qupbrr
      parameter (qupbrr = 12)
      integer qdbrr
      parameter (qdbrr = 13)
      integer qcbr
      parameter (qcbr = 14)
      integer qcbrl
      parameter (qcbrl = 15)
      integer qcbrr
      parameter (qcbrr = 16)
      integer ebrr
      parameter (ebrr = 17)
      integer qsbrl
      parameter (qsbrl = 18)
      integer qsbrr
      parameter (qsbrr = 19)
      integer qtopbr
      parameter (qtopbr = 20)
      integer qupbr
      parameter (qupbr = 21)
      integer qup
      parameter (qup = 22)
      integer qd
      parameter (qd = 23)
      integer qbt
      parameter (qbt = 24)
      integer qupl
      parameter (qupl = 25)
      integer qdl
      parameter (qdl = 26)
      integer qupr
      parameter (qupr = 27)
      integer qdr
      parameter (qdr = 28)
      integer el
      parameter (el = 29)
      integer nue
      parameter (nue = 30)
      integer qc
      parameter (qc = 31)
      integer qcr
      parameter (qcr = 32)
      integer er
      parameter (er = 33)
      integer qsl
      parameter (qsl = 34)
      integer qsr
      parameter (qsr = 35)
      integer qtop
      parameter (qtop = 36)
      integer qcl
      parameter (qcl = 37)
      integer xpp
      parameter (xpp = 38)
      integer ypp
      parameter (ypp = 39)
      integer zzww
      parameter (zzww = 40)
      integer higgs
      parameter (higgs = 41)
      integer xglu
      parameter (xglu = 42)
      integer yglu
      parameter (yglu = 43)
      integer prtcl__chcg(43)
      integer prtcl__cnschg(43)
      integer prtcl__lortyp(43)
      integer prtcl__nlordof(43)
      integer prtcl__colrep(43)
      integer prtcl__zcoup(43)
      integer chcg(43)
      integer cnschg(43)
      integer glu
      parameter (glu = 1)
      integer frm
      parameter (frm = 3)
      integer frmbr
      parameter (frmbr = 4)
      integer msgbs
      parameter (msgbs = 2)
      integer frmhl
      parameter (frmhl = 5)
      integer frmbrhl
      parameter (frmbrhl = 6)
      integer frmhr
      parameter (frmhr = 7)
      integer frmbrhr
      parameter (frmbrhr = 8)
      integer axww
      parameter (axww = 9)
      integer axzzww
      parameter (axzzww = 10)
      integer pht
      parameter (pht = 11)
      integer scal
      parameter (scal = 12)
      integer xgl
      parameter (xgl = 13)
      integer up
      parameter (up = 1)
      integer dn
      parameter (dn = 2)
      integer lpt
      parameter (lpt = 3)
      integer conv(43)
c                        ! NPRCMX is the number of processes which needs
      integer nprcmx, nintmax
      parameter (nprcmx = 1, nintmax = 50)
!                      to be computed simultaneously and NINTMAX the mazimal number of interaction terms
      integer tabint(151,1)      ! interactions used for the computation
      integer tabintmrg(451,1)   ! interactions used for the computation
      integer opam(50,1), opfs(150,1), nvasfs(151,1), nvasam(51,1)
      double precision cpam(50,1), cpfs(150,1)
      integer nprc
      common /prcint/cpam, cpfs, tabint, tabintmrg, opam, opfs, nvasfs, 
     1   nvasam, nprc
C *** END MODULE prcint ***
