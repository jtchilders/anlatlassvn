!
      subroutine inital
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 20:13:17 2/17/05
C...Switches:                     
      integer nmax
      parameter (nmax = 10)
      integer nlor
      parameter (nlor = 6)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
      common /couplings/masses, width, coup, gstrong, g2weak, gbar, ctw
     1   , stw, glup, gldn, gllp, grup, grdn, grlp, emch, trihcp, yuk, 
     2   qrtc, qrtca, qrtc2, qrtc3, qrtchh, qrtchg, ggh, oper
      double precision masses(43)
      double precision width(43)
      integer oper(43,43,43)
      double precision coup(43,43,43)
      double precision gstrong
      double precision g2weak
      double precision gbar
      double precision ctw
      double precision stw
      double precision glup
      double precision gldn
      double precision gllp
      double precision grup
      double precision grdn
      double precision grlp
      double precision emch
      double precision trihcp
      double precision yuk
      double precision qrtc(3,4)
      double precision qrtca(3,4)
      double precision qrtc2(2,4)
      double precision qrtc3(2,4)
      double precision qrtchh(3)
      double precision qrtchg(2)
      double precision ggh
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
      common /prcint/cpam, cpfs, tabint, tabintmrg, opam, opfs, nvasfs, 
     1   nvasam, nprc
      integer nprcmx
      parameter (nprcmx = 1)
      integer nintmax
      parameter (nintmax = 40)
      integer tabint(121,1)
      integer tabintmrg(361,1)
      integer opam(40,1)
      integer opfs(120,1)
      integer nvasfs(121,1)
      integer nvasam(41,1)
      double precision cpam(40,1)
      double precision cpfs(120,1)
      integer nprc
      double precision apar(100)
      common /alphapar/apar
      integer j1
      tabint(1,1) = 13
      tabint(2,1) = 1
      tabint(3,1) = 1
      tabint(4,1) = 41
      tabint(5,1) = 1
      tabint(6,1) = 1
      tabint(7,1) = 1
      tabint(8,1) = 1
      tabint(9,1) = 42
      tabint(10,1) = 41
      tabint(11,1) = 1
      tabint(12,1) = 1
      tabint(13,1) = 43
      tabint(14,1) = 42
      tabint(15,1) = 42
      tabint(16,1) = 41
      tabint(17,1) = 25
      tabint(18,1) = 8
      tabint(19,1) = 1
      tabint(20,1) = 27
      tabint(21,1) = 12
      tabint(22,1) = 1
      tabint(23,1) = 28
      tabint(24,1) = 13
      tabint(25,1) = 1
      tabint(26,1) = 26
      tabint(27,1) = 9
      tabint(28,1) = 1
      tabint(29,1) = 32
      tabint(30,1) = 16
      tabint(31,1) = 1
      tabint(32,1) = 37
      tabint(33,1) = 15
      tabint(34,1) = 1
      tabint(35,1) = 35
      tabint(36,1) = 19
      tabint(37,1) = 1
      tabint(38,1) = 34
      tabint(39,1) = 18
      tabint(40,1) = 1
 
! $  tabint(1:121,1)=(/40, qupbrl,qdl,wp,qupl,qdbrl,wm,qcbrl,qsl,wp,qcl,qsbrl,wm,wp,wm,higgs, &
! $                      xpp,wm,wm,ypp,wp,wp, zzww,z,z,wp,wm,zzww,wp,wm,z,qupl,qupbrl,z,qupr,qupbrr,z, &
! $                      qdl,qdbrl,z,qdr,qdbrr,z ,qcl,qcbrl,z,qcr,qcbrr,z,qsl,qsbrl,z,qsr,qsbrr,z, &
! $                      higgs,z,z,higgs,higgs,higgs, &
! $                      qupl,qupbrl,gluon, qupr,qupbrr,gluon, qdr,qdbrr,gluon,qdl,qdbrl,gluon,  &
! $                      qcr,qcbrr,gluon,qcl,qcbrl,gluon,qsr,qsbrr,gluon,qsl,qsbrl,gluon,         &
! $                      qupl,qupbrl,photon, qupr,qupbrr,photon, qdr,qdbrr,photon,qdl,qdbrl,photon,  &
! $                      qcr,qcbrr,photon,qcl,qcbrl,photon,qsr,qsbrr,photon,qsl,qsbrl,photon,         &
! $                      wp,wm,photon, zzww,z,photon,  zzww,photon,photon ,gluon,gluon,gluon         /)
! $ ,gluon,qbt,qbtbr,gluon,       &
! $                      qtop,qtopbr,gluon, gluon,gluon,gluon,higgs,qtop,qtopbr,higgs,qbt,qbtbr /)       ! the first entry is
 
!    the number of relevant interaction, the following entries in group of three label the particles
!    entering the lagrangian interaction term (the ordering is immaterial)
!
!
! initializing coupling constants
      gstrong = 0.D0
      g2weak = 0.D0
      do j1 = 1, 43
         masses(j1) = 0.D0
         width(j1) = 0.D0
      end do
      gstrong = 1.D0
      g2weak = apar(51)*sqrt(2.D0)
      gbar = apar(55)*sqrt(2.D0)
      ctw = apar(53)
      stw = apar(54)
      masses(3) = apar(3)
      masses(4) = masses(3)
      width(3) = apar(4)
      width(4) = width(3)
      masses(2) = apar(1)
      width(2) = apar(2)
      masses(41) = apar(5)
      width(41) = apar(6)
      masses(24) = apar(15)
      masses(7) = masses(24)
      masses(36) = apar(16)
      masses(20) = masses(36)
      ggh = apar(61)
!
      end 
!
!
      subroutine matrix0(flvmlm, posi, imp, hel, result)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 20:13:17 2/17/05
C...Switches:                     
!
      integer nmax
      parameter (nmax = 10)
      integer nlor
      parameter (nlor = 6)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
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
      integer prtcl__chcg(43)
      integer prtcl__cnschg(43)
      integer prtcl__lortyp(43)
      integer prtcl__nlordof(43)
      integer prtcl__colrep(43)
      integer prtcl__zcoup(43)
      integer chcg(43)
      integer cnschg(43)
      integer conv(43)
      common /naming/prtcl__chcg, prtcl__cnschg, prtcl__lortyp, 
     1   prtcl__nlordof, prtcl__colrep, prtcl__zcoup, chcg, cnschg, conv
      INCLUDE 'Vdimen.inc'
      common /extsrc/src, exflg
      double complex src(6)
      logical exflg
      common /prcint/cpam, cpfs, tabint, tabintmrg, opam, opfs, nvasfs, 
     1   nvasam, nprc
      integer prtcl1x__chcg(43)
      integer prtcl1x__cnschg(43)
      integer prtcl1x__lortyp(43)
      integer prtcl1x__nlordof(43)
      integer prtcl1x__colrep(43)
      integer prtcl1x__zcoup(43)
      integer chcg1x(43)
      integer cnschg1x(43)
      integer conv1x(43)
      integer tabint(121,1)
      integer tabintmrg(361,1)
      integer opam(40,1)
      integer opfs(120,1)
      integer nvasfs(121,1)
      integer nvasam(41,1)
      double precision cpam(40,1)
      double precision cpfs(120,1)
      integer nprc
      integer flvmlm(10), posi(2), hel(10)
      double precision imp(40)
      double complex result
!
      integer j1, nfr, flv(10), npr, nlp, nlepbr, spn, hl(10)
      double precision impul(40), plep(4), plepbr(4), mw, ww
      double complex lp(4), lpbr(4), w(4), wpol(4)
      integer clr(20)
      integer color(20)                          !coulor string
      common /colore/color
      integer prcss__flv(10)
      integer prcss__hel(10)
      integer prcss__col(20)
      integer prcss__nfrm
      double precision prcss__mom(40)
      integer convert(-100:100),cnvmlsl(-16:16),cnvmlsr(-16:16),dim(7)
      save convert, cnvmlsl, cnvmlsr, dim
      logical massless(-16:16)
      save massless
      double precision apar(100)
      common /alphapar/apar
      integer flgdual                    !dual (0) or su3 (1) amplitudes
      common /dual/flgdual
      integer rep(10), nprt, nglu, nqrk, nlep, ngb, nphot
      integer j2
      data convert/ 75*-101, higgs, wp, z, photon, -101, 8*-101, nuebr, 
     1   5*-101, qtopbr, qbtbr, qcbr, -101, qupbr, qdbr, gluon, qd, qup
     2   , -101, qc, qbt, qtop, 5*-101, nue, 8*-101, -101, photon, z, wm
     3   , higgs, 75*-101/ 
      data cnvmlsl/ 4*-101, nuebr, ebrl, 6*-101, qcbrl, qsbrl, qupbrl, 
     1   qdbrl, gluon, qdl, qupl, qsl, qcl, 6*-101, el, nue, 4*-101/ 
      data cnvmlsr/ 5*-101, ebrr, 6*-101, qcbrr, qsbrr, qupbrr, qdbrr, 
     1   gluon, qdr, qupr, qsr, qcr, 6*-101, er, 5*-101/ 
      data massless/ 4*.FALSE., 2*.TRUE., 6*.FALSE., 4*.TRUE., .FALSE., 
     1   4*.TRUE., 6*.FALSE., 2*.TRUE., 4*.FALSE./ 
!
      nprc = 1
!
      do j2 = 1, 40
         impul(j2) = -imp(j2)
      end do
      do j2 = 1, 4
         impul(j2+4*(posi(1)-1)) = -impul(j2+4*(posi(1)-1))
      end do
      do j2 = 1, 4
         impul(j2+4*(posi(2)-1)) = -impul(j2+4*(posi(2)-1))
      end do
      do j2 = 1, 20
         clr(j2) = color(j2)
      end do
      do j2 = 1, 10
         hl(j2) = hel(j2)
         flv(j2) = -flvmlm(j2)
      end do
      flv(posi(1)) = -flv(posi(1))
      flv(posi(2)) = -flv(posi(2))
      nfr = 0
      npr = 0
      do j1 = 1, 10
         if (abs(flvmlm(j1)) .lt. 17) then
            if (flvmlm(j1) .ne. 0) then
               nfr = nfr + 1
               if (flvmlm(j1) .lt. 0) hl(j1) = 2*hl(j1)
            endif
         endif
         if (flvmlm(j1) .ne. 1001) npr = npr + 1
      end do
!
      do j2 = 1, 10 - npr
         flv(j2+npr) = -1001
      end do
! $  write(*,*) 'fl',flv
! $  write(*,*) 'hl',hl
!
      exflg = .FALSE.
      do j1 = 1, 10
         if (flv(j1) .eq. (-1001)) then
            flv(j1) = 1001
         else
            if (abs(flv(j1)) .lt. 17) then
               if (massless(flv(j1))) then
                  if (hl(j1) .gt. 0) then
                     flv(j1) = cnvmlsl(flv(j1))
                  else
                     flv(j1) = cnvmlsr(flv(j1))
                  endif
               else
                  flv(j1) = convert(flv(j1))
               endif
            else
               flv(j1) = convert(flv(j1))
            endif
         endif
      end do
!
      do j2 = 1, 10
         prcss__flv(j2) = flv(j2)
         prcss__hel(j2) = hl(j2)
      end do
      do j2 = 1, 20
         prcss__col(j2) = clr(j2)
      end do
      prcss__nfrm = nfr
      do j2 = 1, 40
         prcss__mom(j2) = impul(j2)
      end do
      if (npr.eq.3 .or. npr.eq.4 .or. npr.eq.5 .or. npr.eq.6 .or. npr
     1   .eq.7 .or. npr.eq.8 .or. npr.eq.9 .or. npr.eq.10) then
         dim(1) = 5370
         dim(2) = 5500
         dim(3) = 5550
         dim(4) = 5780
! $    case (9)
! $      dim(d_offsh)=2120; dim(d_offshcol)=1850; dim(d_ampfld)=790; dim(d_ampcol)=400
! $    case (8)
! $      dim(d_offsh)=1060; dim(d_offshcol)=1000; dim(d_ampfld)=400; dim(d_ampcol)=200
! $    case (7)
! $      dim(d_offsh)=500; dim(d_offshcol)=430; dim(d_ampfld)=210; dim(d_ampcol)=110
! $    case (6)
! $      dim(d_offsh)=240; dim(d_offshcol)=230; dim(d_ampfld)=110; dim(d_ampcol)=60
! $    case (5)
! $      dim(d_offsh)=110; dim(d_offshcol)=90; dim(d_ampfld)=50; dim(d_ampcol)=30
! $    case (4)
! $      dim(d_offsh)=50; dim(d_offshcol)=40; dim(d_ampfld)=20; dim(d_ampcol)=10
! $    case (3)
! $      dim(d_offsh)=10; dim(d_offshcol)=20; dim(d_ampfld)=10; dim(d_ampcol)=10
! $    case (10)
! $      dim(d_offsh)=4370; dim(d_offshcol)=4150; dim(d_ampfld)=1550; dim(d_ampcol)=780
! $    case (3,4,5,6,7,8,9,10)
! $      dim(d_offsh)=5370; dim(d_offshcol)=5500; dim(d_ampfld)=5550; dim(d_ampcol)=5780
      else
         write (*, *) 'too large/small number of particles', npr
      endif
      dim(5) = 2**npr
      dim(6) = 2**npr
      dim(7) = 2**npr
      if (flgdual .eq. 1) then
         call amp (prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1      prcss__mom, dim, result)
      else if (flgdual.eq.0 .or. flgdual.eq.2) then
         call ampd (prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1      prcss__mom, dim, result)
      else
         write (*, *) 'failure in INITALP'
         stop 
      endif
!
      end 
 
 
 
 
 
 
