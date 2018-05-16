      block data extsrc_bd 
      INCLUDE 'Vextsrc.inc'
      data exflg/ .FALSE./ 
!
      end 
      block data couplings_bd 
      INCLUDE 'Vcouplings.inc'
      data width/ 43*-99.D0/ 
      data masses/ 43*-99.D0/ 
!
      end 
!
      subroutine initprc
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
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
      integer nmax
      parameter (nmax = 10)
      integer nlor
      parameter (nlor = 6)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
      common /prcint/cpam, cpfs, tabint, tabintmrg, opam, opfs, nvasfs, 
     1   nvasam, nprc
      integer nprcmx
      parameter (nprcmx = 1)
      integer nintmax
      parameter (nintmax = 50)
      integer tabint(151,1)
      integer tabintmrg(451,1)
      integer opam(50,1)
      integer opfs(150,1)
      integer nvasfs(151,1)
      integer nvasam(51,1)
      double precision cpam(50,1)
      double precision cpfs(150,1)
      integer nprc
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
c                                                !auxiliary variables
      integer j1, j2, pos, j3, j4, tmp, nit, nitm, x1, x2, x3
      integer tmpa(3)                            !auxiliary variables
      integer j5
!
! reordering fields in ascendig order, according to alpha naming scheme, TABINT is used for the amplitude
!
      do j1 = 1, 1
         nit = tabint(1,j1)
         do j2 = 1, nit
            pos = 3*(j2 - 1) + 2
            do j3 = pos, pos + 2
               do j4 = j3 + 1, pos + 2
                  if (tabint(j3,j1) .gt. tabint(j4,j1)) then
                     tmp = tabint(j3,j1)
                     tabint(j3,j1) = tabint(j4,j1)
                     tabint(j4,j1) = tmp
                  endif
               end do
            end do
         end do
      end do
!
!  generating TABINTMRG which is used for the iterative steps before amplitude evaluation: to each
!  three-particle vertex do correspond three equation of motion (< three if equal particles are involved
!  again TABINTMRG is ordered according to ALPHA NAMING scheme
!
      do j1 = 1, 1
         nitm = 0
         nit = tabint(1,j1)
         do j2 = 1, nit
            do j5 = 1, 3
               tmpa(j5) = tabint(j5-2+3*j2,j1)
               tabintmrg(j5+1+3*nitm,j1) = tmpa(j5)
            end do
            nitm = nitm + 1
            if (tmpa(2) .ne. tmpa(1)) then
               tabintmrg(1+1+3*nitm,j1) = tmpa(2)
               tabintmrg(2+1+3*nitm,j1) = tmpa(1)
               tabintmrg(3+1+3*nitm,j1) = tmpa(3)
               nitm = nitm + 1
            endif
            if (tmpa(3).ne.tmpa(1) .and. tmpa(3).ne.tmpa(2)) then
               tabintmrg(1+1+3*nitm,j1) = tmpa(3)
               tabintmrg(2+1+3*nitm,j1) = tmpa(1)
               tabintmrg(3+1+3*nitm,j1) = tmpa(2)
               nitm = nitm + 1
            endif
         end do
         tabintmrg(1,j1) = nitm
      end do
!
      do j1 = 1, 1
         nitm = tabintmrg(1,j1)
         do j2 = 1, nitm
            do j3 = j2 + 1, nitm
               if (tabintmrg(3*(j2-1)+2,j1) .gt. tabintmrg(3*(j3-1)+2,j1
     1            )) then
                  do j5 = 1, 3
                     tmpa(j5) = tabintmrg(j5-2+3*j2,j1)
                     tabintmrg(j5-2+3*j2,j1) = tabintmrg(j5-2+3*j3,j1)
                  end do
                  do j5 = 1, 3
                     tabintmrg(j5-2+3*j3,j1) = tmpa(j5)
                  end do
               endif
            end do
         end do
         nvasfs(1,j1) = 0
         do j2 = 1, nitm
            x1 = tabintmrg(3*(j2-1)+2,j1)
            x2 = tabintmrg(3*(j2-1)+3,j1)
            x3 = tabintmrg(3*(j2-1)+4,j1)
            opfs(j2,j1) = oper(x1,x2,x3)
            cpfs(j2,j1) = coup(x1,x2,x3)
            if (x1.eq.1 .or. x2.eq.1 .or. x3.eq.1) then
               if (.not.(x1.eq.41 .or. x2.eq.41 .or. x3.eq.41)) then
                  nvasfs(1,j1) = nvasfs(1,j1) + 1
                  nvasfs(nvasfs(1,j1)+1,j1) = j2
               endif
            endif
         end do
      end do
      do j1 = 1, 1
         nitm = tabint(1,j1)
         do j2 = 1, nitm
            do j3 = j2 + 1, nitm
               if (tabint(3*(j2-1)+2,j1) .gt. tabint(3*(j3-1)+2,j1)) 
     1            then
                  do j5 = 1, 3
                     tmpa(j5) = tabint(j5-2+3*j2,j1)
                     tabint(j5-2+3*j2,j1) = tabint(j5-2+3*j3,j1)
                  end do
                  do j5 = 1, 3
                     tabint(j5-2+3*j3,j1) = tmpa(j5)
                  end do
               endif
            end do
         end do
         nvasam(1,j1) = 0
         do j2 = 1, nitm
            x1 = tabint(3*(j2-1)+2,j1)
            x2 = tabint(3*(j2-1)+3,j1)
            x3 = tabint(3*(j2-1)+4,j1)
            opam(j2,j1) = oper(x1,x2,x3)
            cpam(j2,j1) = coup(x1,x2,x3)
            if (x1.eq.1 .or. x2.eq.1 .or. x3.eq.1) then
               if (.not.(x1.eq.41 .or. x2.eq.41 .or. x3.eq.41)) then
                  nvasam(1,j1) = nvasam(1,j1) + 1
                  nvasam(nvasam(1,j1)+1,j1) = j2
               endif
            endif
         end do
      end do
!
      end 
      block data strfld_bd 
      INCLUDE 'Vstrfld.inc'
      data fldptr/ 430*-9999/ 
      data colrep/ 43*-9999/ 
      data nlordof/ 43*-9999/ 
      data lortyp/ 43*-9999/ 
!
      end 
!
!
      subroutine initmom
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
! auxiliary variables
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer nmax
      parameter (nmax = 10)
      integer nlor
      parameter (nlor = 6)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexcmax
      parameter (nexcmax = 13)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
      integer j1, i1, j2, aux, zr, j3, nf, nf0
      integer s(10)
      double precision pm
      integer popcnt0
      integer j4
      do j4 = 1, 2047
         momd(j4) = 0
      end do
      do j1 = 1, 1024
         momd(j1) = popcnt0(j1)
      end do
!
      do j1 = 2, 5
         aux = 0        ! number of possible J2+J3=J1 non simmetric case
         do j2 = 1, 5
            do j3 = 1, 5
               if (j3 + j2 .eq. j1) then
                  aux = aux + 1
                  nexc(1-1+2*aux,j1) = j2
                  nexc(2-1+2*aux,j1) = j3
               endif
            end do
         end do
         nexc(1,j1) = aux
      end do
!
      do j1 = 2, 5
         aux = 0           ! number of possible J2+J3=J1  simmetric case
         do j2 = 1, 5
            do j3 = j2, 5
               if (j3 + j2 .eq. j1) then
                  aux = aux + 1
                  nexcsm(1-1+2*aux,j1) = j2
                  nexcsm(2-1+2*aux,j1) = j3
               endif
            end do
         end do
         nexcsm(1,j1) = aux
      end do
!
      end 
!
      subroutine fillpr(i1, i2, i3, i4, i5, i6, tmp__chcg, tmp__cnschg, 
     1   tmp__lortyp, tmp__nlordof, tmp__colrep, tmp__zcoup)
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
      integer i1, i2, i3, i4, i5, i6
      integer tmp__chcg
      integer tmp__cnschg
      integer tmp__lortyp
      integer tmp__nlordof
      integer tmp__colrep
      integer tmp__zcoup
!
      tmp__chcg = i1
      tmp__cnschg = i2
      tmp__lortyp = i3
      tmp__nlordof = i4
      tmp__colrep = i5
      tmp__zcoup = i6
!
      end 
!
      subroutine inittab
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
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
      integer j1
      call fillpr (1, 0, 1, 6, 3, -99, prtcl__chcg(1), prtcl__cnschg(1)
     1   , prtcl__lortyp(1), prtcl__nlordof(1), prtcl__colrep(1), 
     2   prtcl__zcoup(1))
      call fillpr (5, 0, 11, 4, 1, -99, prtcl__chcg(5), prtcl__cnschg(5)
     1   , prtcl__lortyp(5), prtcl__nlordof(5), prtcl__colrep(5), 
     2   prtcl__zcoup(5))
      call fillpr (21, 1, 3, 4, 2, 1, prtcl__chcg(22), prtcl__cnschg(22)
     1   , prtcl__lortyp(22), prtcl__nlordof(22), prtcl__colrep(22), 
     2   prtcl__zcoup(22))
      call fillpr (14, 2, 3, 4, 2, 1, prtcl__chcg(31), prtcl__cnschg(31)
     1   , prtcl__lortyp(31), prtcl__nlordof(31), prtcl__colrep(31), 
     2   prtcl__zcoup(31))
      call fillpr (6, 3, 3, 4, 2, 2, prtcl__chcg(23), prtcl__cnschg(23)
     1   , prtcl__lortyp(23), prtcl__nlordof(23), prtcl__colrep(23), 
     2   prtcl__zcoup(23))
      call fillpr (7, 4, 3, 4, 2, 2, prtcl__chcg(24), prtcl__cnschg(24)
     1   , prtcl__lortyp(24), prtcl__nlordof(24), prtcl__colrep(24), 
     2   prtcl__zcoup(24))
      call fillpr (22, 5, 4, 4, 2, 1, prtcl__chcg(21), prtcl__cnschg(21)
     1   , prtcl__lortyp(21), prtcl__nlordof(21), prtcl__colrep(21), 
     2   prtcl__zcoup(21))
      call fillpr (31, 6, 4, 4, 2, 1, prtcl__chcg(14), prtcl__cnschg(14)
     1   , prtcl__lortyp(14), prtcl__nlordof(14), prtcl__colrep(14), 
     2   prtcl__zcoup(14))
      call fillpr (23, 7, 4, 4, 2, 2, prtcl__chcg(6), prtcl__cnschg(6), 
     1   prtcl__lortyp(6), prtcl__nlordof(6), prtcl__colrep(6), 
     2   prtcl__zcoup(6))
      call fillpr (24, 8, 4, 4, 2, 2, prtcl__chcg(7), prtcl__cnschg(7), 
     1   prtcl__lortyp(7), prtcl__nlordof(7), prtcl__colrep(7), 
     2   prtcl__zcoup(7))
      call fillpr (8, 9, 5, 2, 2, -99, prtcl__chcg(25), prtcl__cnschg(25
     1   ), prtcl__lortyp(25), prtcl__nlordof(25), prtcl__colrep(25), 
     2   prtcl__zcoup(25))
      call fillpr (15, 10, 5, 2, 2, -99, prtcl__chcg(37), prtcl__cnschg(
     1   37), prtcl__lortyp(37), prtcl__nlordof(37), prtcl__colrep(37), 
     2   prtcl__zcoup(37))
      call fillpr (9, 11, 5, 2, 2, -99, prtcl__chcg(26), prtcl__cnschg(
     1   26), prtcl__lortyp(26), prtcl__nlordof(26), prtcl__colrep(26), 
     2   prtcl__zcoup(26))
      call fillpr (12, 12, 7, 2, 2, -99, prtcl__chcg(27), prtcl__cnschg(
     1   27), prtcl__lortyp(27), prtcl__nlordof(27), prtcl__colrep(27), 
     2   prtcl__zcoup(27))
      call fillpr (16, 13, 7, 2, 2, -99, prtcl__chcg(32), prtcl__cnschg(
     1   32), prtcl__lortyp(32), prtcl__nlordof(32), prtcl__colrep(32), 
     2   prtcl__zcoup(32))
      call fillpr (13, 14, 7, 2, 2, -99, prtcl__chcg(28), prtcl__cnschg(
     1   28), prtcl__lortyp(28), prtcl__nlordof(28), prtcl__colrep(28), 
     2   prtcl__zcoup(28))
      call fillpr (25, 15, 6, 2, 2, -99, prtcl__chcg(8), prtcl__cnschg(8
     1   ), prtcl__lortyp(8), prtcl__nlordof(8), prtcl__colrep(8), 
     2   prtcl__zcoup(8))
      call fillpr (37, 16, 6, 2, 2, -99, prtcl__chcg(15), prtcl__cnschg(
     1   15), prtcl__lortyp(15), prtcl__nlordof(15), prtcl__colrep(15), 
     2   prtcl__zcoup(15))
      call fillpr (26, 17, 6, 2, 2, -99, prtcl__chcg(9), prtcl__cnschg(9
     1   ), prtcl__lortyp(9), prtcl__nlordof(9), prtcl__colrep(9), 
     2   prtcl__zcoup(9))
      call fillpr (27, 18, 8, 2, 2, 1, prtcl__chcg(12), prtcl__cnschg(12
     1   ), prtcl__lortyp(12), prtcl__nlordof(12), prtcl__colrep(12), 
     2   prtcl__zcoup(12))
      call fillpr (32, 19, 8, 2, 2, 1, prtcl__chcg(16), prtcl__cnschg(16
     1   ), prtcl__lortyp(16), prtcl__nlordof(16), prtcl__colrep(16), 
     2   prtcl__zcoup(16))
      call fillpr (28, 20, 8, 2, 2, 2, prtcl__chcg(13), prtcl__cnschg(13
     1   ), prtcl__lortyp(13), prtcl__nlordof(13), prtcl__colrep(13), 
     2   prtcl__zcoup(13))
      call fillpr (4, 21, 2, 5, 1, -99, prtcl__chcg(3), prtcl__cnschg(3)
     1   , prtcl__lortyp(3), prtcl__nlordof(3), prtcl__colrep(3), 
     2   prtcl__zcoup(3))
      call fillpr (3, 22, 2, 5, 1, -99, prtcl__chcg(4), prtcl__cnschg(4)
     1   , prtcl__lortyp(4), prtcl__nlordof(4), prtcl__colrep(4), 
     2   prtcl__zcoup(4))
      call fillpr (2, 0, 2, 5, 1, -99, prtcl__chcg(2), prtcl__cnschg(2)
     1   , prtcl__lortyp(2), prtcl__nlordof(2), prtcl__colrep(2), 
     2   prtcl__zcoup(2))
      call fillpr (29, 23, 6, 2, 1, -99, prtcl__chcg(10), prtcl__cnschg(
     1   10), prtcl__lortyp(10), prtcl__nlordof(10), prtcl__colrep(10), 
     2   prtcl__zcoup(10))
      call fillpr (10, 24, 5, 2, 1, -99, prtcl__chcg(29), prtcl__cnschg(
     1   29), prtcl__lortyp(29), prtcl__nlordof(29), prtcl__colrep(29), 
     2   prtcl__zcoup(29))
      call fillpr (33, 25, 8, 2, 1, -99, prtcl__chcg(17), prtcl__cnschg(
     1   17), prtcl__lortyp(17), prtcl__nlordof(17), prtcl__colrep(17), 
     2   prtcl__zcoup(17))
      call fillpr (17, 26, 7, 2, 1, -99, prtcl__chcg(33), prtcl__cnschg(
     1   33), prtcl__lortyp(33), prtcl__nlordof(33), prtcl__colrep(33), 
     2   prtcl__zcoup(33))
      call fillpr (30, 27, 6, 2, 1, -99, prtcl__chcg(11), prtcl__cnschg(
     1   11), prtcl__lortyp(11), prtcl__nlordof(11), prtcl__colrep(11), 
     2   prtcl__zcoup(11))
      call fillpr (11, 28, 5, 2, 1, -99, prtcl__chcg(30), prtcl__cnschg(
     1   30), prtcl__lortyp(30), prtcl__nlordof(30), prtcl__colrep(30), 
     2   prtcl__zcoup(30))
      call fillpr (18, 29, 5, 2, 2, -99, prtcl__chcg(34), prtcl__cnschg(
     1   34), prtcl__lortyp(34), prtcl__nlordof(34), prtcl__colrep(34), 
     2   prtcl__zcoup(34))
      call fillpr (19, 30, 7, 2, 2, -99, prtcl__chcg(35), prtcl__cnschg(
     1   35), prtcl__lortyp(35), prtcl__nlordof(35), prtcl__colrep(35), 
     2   prtcl__zcoup(35))
      call fillpr (34, 31, 6, 2, 2, -99, prtcl__chcg(18), prtcl__cnschg(
     1   18), prtcl__lortyp(18), prtcl__nlordof(18), prtcl__colrep(18), 
     2   prtcl__zcoup(18))
      call fillpr (35, 32, 8, 2, 2, -99, prtcl__chcg(19), prtcl__cnschg(
     1   19), prtcl__lortyp(19), prtcl__nlordof(19), prtcl__colrep(19), 
     2   prtcl__zcoup(19))
      call fillpr (20, 33, 3, 4, 2, 1, prtcl__chcg(36), prtcl__cnschg(36
     1   ), prtcl__lortyp(36), prtcl__nlordof(36), prtcl__colrep(36), 
     2   prtcl__zcoup(36))
      call fillpr (36, 34, 4, 4, 2, 1, prtcl__chcg(20), prtcl__cnschg(20
     1   ), prtcl__lortyp(20), prtcl__nlordof(20), prtcl__colrep(20), 
     2   prtcl__zcoup(20))
      call fillpr (39, 35, 9, 1, 1, -99, prtcl__chcg(38), prtcl__cnschg(
     1   38), prtcl__lortyp(38), prtcl__nlordof(38), prtcl__colrep(38), 
     2   prtcl__zcoup(38))
      call fillpr (38, 36, 9, 1, 1, -99, prtcl__chcg(39), prtcl__cnschg(
     1   39), prtcl__lortyp(39), prtcl__nlordof(39), prtcl__colrep(39), 
     2   prtcl__zcoup(39))
      call fillpr (40, 0, 10, 2, 1, -99, prtcl__chcg(40), prtcl__cnschg(
     1   40), prtcl__lortyp(40), prtcl__nlordof(40), prtcl__colrep(40), 
     2   prtcl__zcoup(40))
      call fillpr (41, 0, 12, 3, 1, -99, prtcl__chcg(41), prtcl__cnschg(
     1   41), prtcl__lortyp(41), prtcl__nlordof(41), prtcl__colrep(41), 
     2   prtcl__zcoup(41))
      call fillpr (42, 0, 13, 3, 3, -99, prtcl__chcg(43), prtcl__cnschg(
     1   43), prtcl__lortyp(43), prtcl__nlordof(43), prtcl__colrep(43), 
     2   prtcl__zcoup(43))
      call fillpr (43, 0, 13, 3, 3, -99, prtcl__chcg(42), prtcl__cnschg(
     1   42), prtcl__lortyp(42), prtcl__nlordof(42), prtcl__colrep(42), 
     2   prtcl__zcoup(42))
!
      do j1 = 1, 43
         chcg(j1) = prtcl__chcg(j1)
         cnschg(j1) = prtcl__cnschg(j1)
         lortyp(j1) = prtcl__lortyp(j1)
         nlordof(j1) = prtcl__nlordof(j1)
         colrep(j1) = prtcl__colrep(j1)
         conv(j1) = prtcl__zcoup(j1)
      end do
!
      end 
!
      subroutine initcoup
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      double precision stw2, sqd2, ch(43)
      integer tmp1(43), tmp2(43)
      integer j1, nt
      doubleprecision d1(5), d2(5), d3(5)
      integer j2, j3
!
      call inital
!
      sqd2 = sqrt(2.D0)
      oper(1,1,1) = 1
      coup(1,1,1) = gstrong/sqd2
!
      oper(42,41,42) = 47
      coup(42,41,42) = ggh
      oper(41,42,42) = 48
      coup(41,42,42) = ggh
      oper(1,1,43) = 49
      coup(1,1,43) = gstrong/sqd2
      oper(43,1,1) = 50
      coup(43,1,1) = gstrong/sqd2
      oper(1,41,42) = 51
      coup(1,41,42) = ggh
      oper(41,1,42) = 52
      coup(41,1,42) = ggh
      oper(42,1,41) = 53
      coup(42,1,41) = ggh
      oper(1,1,41) = 45
      coup(1,1,41) = ggh
      oper(41,1,1) = 46
      coup(41,1,1) = ggh
      tmp1(1) = 6
      tmp1(2) = 7
      tmp1(3) = 14
      tmp1(4) = 20
      tmp1(5) = 21
      tmp2(1) = 23
      tmp2(2) = 24
      tmp2(3) = 31
      tmp2(4) = 36
      tmp2(5) = 22
      do j1 = 1, 5
         oper(41,tmp1(j1),tmp2(j1)) = 40
         oper(tmp1(j1),tmp2(j1),41) = 41
         oper(tmp2(j1),tmp1(j1),41) = 41
         yuk = -masses(tmp1(j1))*gbar/(sqd2*2.D0*masses(2))
         coup(41,tmp1(j1),tmp2(j1)) = yuk
         coup(tmp1(j1),tmp2(j1),41) = yuk
         coup(tmp2(j1),tmp1(j1),41) = yuk
      end do
      trihcp = -masses(41)**2/masses(2)*0.25*gbar/sqd2
      qrtchh(1) = 6.D0
      qrtchh(2) = gbar/(8.D0*masses(2)*sqd2)*2.D0
      qrtchh(3) = -2.D0/trihcp
      oper(41,41,41) = 39
      coup(41,41,41) = trihcp
      qrtchg(1) = g2weak/(4.D0*masses(3)*sqd2)
      qrtchg(2) = gbar/(4.D0*masses(2)*sqd2)
      oper(3,4,41) = 38
      coup(3,4,41) = masses(3)*g2weak/sqd2
      oper(4,3,41) = 38
      coup(4,3,41) = masses(3)*g2weak/sqd2
      oper(41,3,4) = 37
      coup(41,3,4) = masses(3)*g2weak/sqd2
      oper(2,2,41) = 38
      coup(2,2,41) = masses(2)*gbar/sqd2
      oper(41,2,2) = 37
      coup(41,2,2) = masses(2)*gbar/sqd2
      do j3 = 1, 4
         do j2 = 1, 3
            qrtc(j2,j3) = 0.D0
         end do
      end do
      do j3 = 1, 4
         do j2 = 1, 2
            qrtc2(j2,j3) = 0.D0
            qrtc3(j2,j3) = 0.D0
         end do
      end do
      qrtc2(1,1) = 0.D0
      qrtc2(2,1) = 2.D0
      qrtc2(1,2) = 1.D0
      qrtc2(2,2) = 0.D0
      qrtc2(1,3) = 1.D0
      qrtc2(2,3) = 0.D0
      qrtc2(1,4) = 0.D0
      qrtc2(2,4) = 2.D0
      qrtc3(1,1) = 2.D0
      qrtc3(2,1) = 0.D0
      qrtc3(1,2) = 0.D0
      qrtc3(2,2) = 1.D0
      qrtc3(1,4) = 2.D0
      qrtc3(2,4) = 0.D0
      oper(3,4,40) = 28
      oper(4,3,40) = 28
      oper(40,3,4) = 29
      coup(3,4,40) = -g2weak/sqd2
      coup(4,3,40) = -g2weak/sqd2
      coup(40,3,4) = -g2weak/sqd2
      oper(2,2,40) = 28
      oper(40,2,2) = 29
      coup(2,2,40) = -g2weak/sqd2*ctw*ctw
      coup(40,2,2) = -g2weak/sqd2*ctw*ctw
      oper(2,5,40) = 33
      oper(5,2,40) = 32
      oper(40,2,5) = 35
      coup(2,5,40) = -g2weak/sqd2*ctw*stw
      coup(40,2,5) = -g2weak/sqd2*ctw*stw
      coup(5,2,40) = -g2weak/sqd2*ctw*stw
      oper(5,5,40) = 34
      oper(40,5,5) = 36
      coup(5,5,40) = -g2weak/sqd2*stw*stw
      coup(40,5,5) = -g2weak/sqd2*stw*stw
      oper(3,3,39) = 27
      oper(39,3,3) = 26
      coup(3,3,39) = -g2weak*sqd2
      coup(39,3,3) = -g2weak*sqd2
      oper(4,4,38) = 27
      oper(38,4,4) = 26
      coup(4,4,38) = g2weak/sqd2
      coup(38,4,4) = g2weak/sqd2
      oper(2,3,4) = 25
      oper(3,2,4) = 25
      oper(4,2,3) = 25
      coup(2,3,4) = -g2weak*ctw/sqd2
      coup(3,2,4) = g2weak*ctw/sqd2
      coup(4,2,3) = -g2weak*ctw/sqd2
      oper(5,3,4) = 30
      oper(3,4,5) = 31
      oper(4,3,5) = 31
      coup(5,3,4) = -g2weak*stw/sqd2
      coup(3,4,5) = -g2weak*stw/sqd2
      coup(4,3,5) = g2weak*stw/sqd2
      do j3 = 1, 4
         do j2 = 1, 3
            qrtc(j2,j3) = 0.D0
         end do
      end do
      qrtc(1,1) = -1.D0/ctw
      qrtc(2,2) = -1.D0/ctw
      qrtc(2,3) = 1.D0/ctw                       !Z5 W+ W-
      qrtc(2,1) = 1.D0
      qrtc(1,2) = 1.D0
      qrtc(3,3) = 1.D0                           !W+5 Z W-
      qrtc(3,1) = -1.D0
      qrtc(3,2) = 1.D0
      qrtc(1,3) = 1.D0                           !W-5 Z W+
      do j3 = 1, 4
         do j2 = 1, 3
c! qrtca(1,4)=-1.d0/stw; qrtca(3,2)=1.d0/stw; qrtca(3,3)=-1.d0/stw;  !Z5
            qrtca(j2,j3) = 0.D0
         end do
      end do
      qrtca(2,4) = 1.D0
      qrtca(1,2) = -1.D0
      qrtca(2,3) = -1.D0                         !W+5 Z W-
      qrtca(3,4) = -1.D0
      qrtca(2,2) = -1.D0
      qrtca(1,3) = -1.D0                         !W-5 Z W+
      tmp1(1) = 21
      tmp1(2) = 6
      tmp1(3) = 14
      tmp1(4) = 7
      tmp1(5) = 20
      tmp2(1) = 22
      tmp2(2) = 23
      tmp2(3) = 31
      tmp2(4) = 24
      tmp2(5) = 36
      do j1 = 1, 5
         oper(1,tmp1(j1),tmp2(j1)) = 4
         oper(tmp2(j1),1,tmp1(j1)) = 5
         oper(tmp1(j1),1,tmp2(j1)) = 6
         coup(1,tmp1(j1),tmp2(j1)) = gstrong/sqd2
         coup(tmp2(j1),1,tmp1(j1)) = gstrong/sqd2
         coup(tmp1(j1),1,tmp2(j1)) = gstrong/sqd2
      end do
      tmp1(1) = 8
      tmp1(2) = 9
      tmp1(3) = 15
      tmp1(4) = 18
      tmp2(1) = 25
      tmp2(2) = 26
      tmp2(3) = 37
      tmp2(4) = 34
      do j1 = 1, 4
         oper(1,tmp1(j1),tmp2(j1)) = 10
         oper(tmp2(j1),1,tmp1(j1)) = 11
         oper(tmp1(j1),1,tmp2(j1)) = 12
         coup(1,tmp1(j1),tmp2(j1)) = gstrong/sqd2
         coup(tmp2(j1),1,tmp1(j1)) = gstrong/sqd2
         coup(tmp1(j1),1,tmp2(j1)) = gstrong/sqd2
      end do
      tmp1(1) = 12
      tmp1(2) = 13
      tmp1(3) = 16
      tmp1(4) = 19
      tmp2(1) = 27
      tmp2(2) = 28
      tmp2(3) = 32
      tmp2(4) = 35
      do j1 = 1, 4
         oper(1,tmp1(j1),tmp2(j1)) = 16
         oper(tmp2(j1),1,tmp1(j1)) = 17
         oper(tmp1(j1),1,tmp2(j1)) = 18
         coup(1,tmp1(j1),tmp2(j1)) = gstrong/sqd2
         coup(tmp2(j1),1,tmp1(j1)) = gstrong/sqd2
         coup(tmp1(j1),1,tmp2(j1)) = gstrong/sqd2
      end do
      tmp1(1) = 6
      tmp1(2) = 7
      tmp2(1) = 22
      tmp2(2) = 36
      do j1 = 1, 2
         oper(4,tmp1(j1),tmp2(j1)) = 7
         oper(tmp2(j1),4,tmp1(j1)) = 8
         oper(tmp1(j1),4,tmp2(j1)) = 9
         coup(4,tmp1(j1),tmp2(j1)) = g2weak/2.D0
         coup(tmp2(j1),4,tmp1(j1)) = g2weak/2.D0
         coup(tmp1(j1),4,tmp2(j1)) = g2weak/2.D0
      end do
      tmp1(1) = 21
      tmp1(2) = 20
      tmp2(1) = 23
      tmp2(2) = 24
      do j1 = 1, 2
         oper(3,tmp1(j1),tmp2(j1)) = 7
         oper(tmp2(j1),3,tmp1(j1)) = 8
         oper(tmp1(j1),3,tmp2(j1)) = 9
         coup(3,tmp1(j1),tmp2(j1)) = g2weak/2.D0
         coup(tmp2(j1),3,tmp1(j1)) = g2weak/2.D0
         coup(tmp1(j1),3,tmp2(j1)) = g2weak/2.D0
      end do
      tmp1(1) = 9
      tmp1(2) = 10
      tmp1(3) = 18
      tmp2(1) = 25
      tmp2(2) = 30
      tmp2(3) = 37
      do j1 = 1, 3
         oper(4,tmp1(j1),tmp2(j1)) = 13
         oper(tmp2(j1),4,tmp1(j1)) = 14
         oper(tmp1(j1),4,tmp2(j1)) = 15
         coup(4,tmp1(j1),tmp2(j1)) = g2weak/2.D0
         coup(tmp2(j1),4,tmp1(j1)) = g2weak/2.D0
         coup(tmp1(j1),4,tmp2(j1)) = g2weak/2.D0
      end do
      tmp1(1) = 8
      tmp1(2) = 11
      tmp1(3) = 15
      tmp2(1) = 26
      tmp2(2) = 29
      tmp2(3) = 34
      do j1 = 1, 3
         oper(3,tmp1(j1),tmp2(j1)) = 13
         oper(tmp2(j1),3,tmp1(j1)) = 14
         oper(tmp1(j1),3,tmp2(j1)) = 15
         coup(3,tmp1(j1),tmp2(j1)) = g2weak/2.D0
         coup(tmp2(j1),3,tmp1(j1)) = g2weak/2.D0
         coup(tmp1(j1),3,tmp2(j1)) = g2weak/2.D0
      end do
      stw2 = stw**2
      glup = -gbar/sqd2*(1.D0 - 4.D0/3.D0*stw2)
      gldn = gbar/sqd2*(1.D0 - 2.D0/3.D0*stw2)
      gllp = gbar/sqd2*(1.D0 - 2.D0*stw2)
      grup = gbar*2.D0*sqd2/3.D0*stw2
      grdn = -gbar*sqd2/3.D0*stw2
      grlp = -gbar*stw2*sqd2
      emch = -g2weak*stw/sqd2
      tmp1(1) = 8
      tmp1(2) = 15
      tmp1(3) = 9
      tmp1(4) = 10
      tmp1(5) = 18
      tmp1(6) = 11
      tmp2(1) = 25
      tmp2(2) = 37
      tmp2(3) = 26
      tmp2(4) = 29
      tmp2(5) = 34
      tmp2(6) = 30
      ch(1) = glup/2.D0
      ch(2) = glup/2.D0
      ch(3) = gldn/2.D0
      ch(4) = gllp/2.D0
      ch(5) = gldn/2.D0
      ch(6) = -gbar/sqd2/2.D0
      do j1 = 1, 6
         oper(2,tmp1(j1),tmp2(j1)) = 13
         oper(tmp2(j1),2,tmp1(j1)) = 14
         oper(tmp1(j1),2,tmp2(j1)) = 15
         coup(2,tmp1(j1),tmp2(j1)) = ch(j1)
         coup(tmp2(j1),2,tmp1(j1)) = ch(j1)
         coup(tmp1(j1),2,tmp2(j1)) = ch(j1)
      end do
      tmp1(1) = 12
      tmp1(2) = 16
      tmp1(3) = 13
      tmp1(4) = 17
      tmp1(5) = 19
      tmp2(1) = 27
      tmp2(2) = 32
      tmp2(3) = 28
      tmp2(4) = 33
      tmp2(5) = 35
      ch(1) = grup/2.D0
      ch(2) = grup/2.D0
      ch(3) = grdn/2.D0
      ch(4) = grlp/2.D0
      ch(5) = grdn/2.D0
      do j1 = 1, 5
         oper(2,tmp1(j1),tmp2(j1)) = 19
         oper(tmp2(j1),2,tmp1(j1)) = 20
         oper(tmp1(j1),2,tmp2(j1)) = 21
         coup(2,tmp1(j1),tmp2(j1)) = ch(j1)
         coup(tmp2(j1),2,tmp1(j1)) = ch(j1)
         coup(tmp1(j1),2,tmp2(j1)) = ch(j1)
      end do
      tmp1(1) = 21
      tmp1(2) = 6
      tmp1(3) = 7
      tmp1(4) = 20
      tmp2(1) = 22
      tmp2(2) = 23
      tmp2(3) = 24
      tmp2(4) = 36
      do j1 = 1, 4
         oper(2,tmp1(j1),tmp2(j1)) = 22
         oper(tmp2(j1),2,tmp1(j1)) = 23
         oper(tmp1(j1),2,tmp2(j1)) = 24
         coup(2,tmp1(j1),tmp2(j1)) = 1.D0/2.D0
         coup(tmp2(j1),2,tmp1(j1)) = 1.D0/2.D0
         coup(tmp1(j1),2,tmp2(j1)) = 1.D0/2.D0
      end do
      tmp1(1) = 21
      tmp1(2) = 6
      tmp1(3) = 14
      tmp1(4) = 7
      tmp1(5) = 20
      tmp2(1) = 22
      tmp2(2) = 23
      tmp2(3) = 31
      tmp2(4) = 24
      tmp2(5) = 36
      d1(1) = 2.D0/3.D0*emch
      d1(2) = -1.D0/3.D0*emch
      d1(3) = 2.D0/3.D0*emch
      d1(4) = -1.D0/3.D0*emch
      d1(5) = 2.D0/3.D0*emch
      do j2 = 1, 5
         ch(j2) = -d1(j2)
      end do
      do j1 = 1, 5
         oper(5,tmp1(j1),tmp2(j1)) = 4
         oper(tmp2(j1),5,tmp1(j1)) = 5
         oper(tmp1(j1),5,tmp2(j1)) = 6
         coup(5,tmp1(j1),tmp2(j1)) = ch(j1)
         coup(tmp2(j1),5,tmp1(j1)) = ch(j1)
         coup(tmp1(j1),5,tmp2(j1)) = ch(j1)
      end do
      tmp1(1) = 8
      tmp1(2) = 9
      tmp1(3) = 15
      tmp1(4) = 10
      tmp1(5) = 18
      tmp2(1) = 25
      tmp2(2) = 26
      tmp2(3) = 37
      tmp2(4) = 29
      tmp2(5) = 34
      d2(1) = 2.D0/3.D0*emch
      d2(2) = -1.D0/3.D0*emch
      d2(3) = 2.D0/3.D0*emch
      d2(4) = -emch
      d2(5) = -1.D0/3.D0*emch
      do j2 = 1, 5
         ch(j2) = -d2(j2)
      end do
      do j1 = 1, 5
         oper(5,tmp1(j1),tmp2(j1)) = 10
         oper(tmp2(j1),5,tmp1(j1)) = 11
         oper(tmp1(j1),5,tmp2(j1)) = 12
         coup(5,tmp1(j1),tmp2(j1)) = ch(j1)
         coup(tmp2(j1),5,tmp1(j1)) = ch(j1)
         coup(tmp1(j1),5,tmp2(j1)) = ch(j1)
      end do
      tmp1(1) = 12
      tmp1(2) = 13
      tmp1(3) = 16
      tmp1(4) = 17
      tmp1(5) = 19
      tmp2(1) = 27
      tmp2(2) = 28
      tmp2(3) = 32
      tmp2(4) = 33
      tmp2(5) = 35
      d3(1) = 2.D0/3.D0*emch
      d3(2) = -1.D0/3.D0*emch
      d3(3) = 2.D0/3.D0*emch
      d3(4) = -emch
      d3(5) = -1.D0/3.0*emch
      do j2 = 1, 5
         ch(j2) = -d3(j2)
      end do
      do j1 = 1, 5
         oper(5,tmp1(j1),tmp2(j1)) = 16
         oper(tmp2(j1),5,tmp1(j1)) = 17
         oper(tmp1(j1),5,tmp2(j1)) = 18
         coup(5,tmp1(j1),tmp2(j1)) = ch(j1)
         coup(tmp2(j1),5,tmp1(j1)) = ch(j1)
         coup(tmp1(j1),5,tmp2(j1)) = ch(j1)
      end do
!
      end 
!
      subroutine rncpl
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vrunning.inc'
      double precision apar(100)
      common /alphapar/apar
!
c                   !alpha strong, to be changed on event by event basis
      rnas = apar(56)
!
      end 
!
      subroutine amp(prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1   prcss__mom, dim0, mtel)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexcmax
      parameter (nexcmax = 13)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
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
      parameter (nintmax = 50)
      integer tabint(151,1)
      integer tabintmrg(451,1)
      integer opam(50,1)
      integer opfs(150,1)
      integer nvasfs(151,1)
      integer nvasam(51,1)
      double precision cpam(50,1)
      double precision cpfs(150,1)
      integer nprc
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
      INCLUDE 'Vincol.inc'
      INCLUDE 'Vdimen.inc'
      INCLUDE 'Vutilities.inc'
      INCLUDE 'Vrunning.inc'
      integer prcss__flv(10)
      integer prcss__hel(10)
      integer prcss__col(20)
      integer prcss__nfrm
      double precision prcss__mom(40)
      integer dim0(7)
      double complex mtel                        !amplitude
!
      integer nca, dim(7)
      parameter (nca = 3)
      integer nprt, nfld, flvold, padd, nconf, j5, padd0
      integer j1,j2,iter,steptm,stepex,ex2,ex3,nf2,nf3,psfl2,psfl3
      integer pscfl2,pscfl3,nl1,nl2,nl3,j3,j4,beg3,nc,psfl3i,pscfl3i
      integer p1, p2, p3, posfld, poscol, nlold
      integer nf1, psfl1, pscfl1, pref, beg
      integer cl(2), clnw(2)
      integer psnew(1050)
      integer fststp
      save fststp
      integer tabmrgpr(451)
      integer exc(13)
      integer offshcol(5500), ampcol(5780), ordtmp(1050)
      double complex offsh(5370), ampfld(5550)
      double complex fusion(6), polar(6)
      double complex amplitude, tmp
      double precision colcff(2)
      double precision impul(4,1050)
      double precision imp(4)
      double precision perm
      integer colconv(0:3,0:3)
      save colconv
      integer cmbcol(-11:11,-11:11)
      save cmbcol
!
      integer nphf, even, nmom, ot, npmit, tmpi, ci, cf, j5b, ptr, nfr2
      integer operot, jop, op(150), sim
      logical choose, choose2, choose3, choose4, cnsold(0:43)
      double precision imp3(4), imp2(4)
      double precision coupot, cpot1, cp(150)
      integer psfit, pscit, psfa, psca
      save psfit, pscit, psfa, psca
      logical dbg
      parameter (dbg = .FALSE.)
      integer j9, j10
      doubleprecision d1(4), d2(4)
      doublecomplex x1
      integer j6, j7, j8
!
! once for all inizialization
      dim(1) = 5370
      dim(2) = 5500
      dim(3) = 5550
      dim(4) = 5780
      dim(5) = 1050
      dim(6) = 1050
      dim(7) = 1050
      data psca/ 0/ 
      data psfa/ 0/ 
      data pscit/ 0/ 
      data psfit/ 0/ 
      data fststp/ 0/ 
      if (fststp .eq. 0) then
         fststp = 1
! inizialization of arrays
         call initcoup
         call initmom
         call inittab
         call colprdsu3
!
! inizialization of merging tables
         call initprc
!
!   converting from (qb,q) coulor represenation to internal coulor representation
!   1,8 gluon; (-) 9,11 (anti) quarks; 0 coulorless
         do j7 = 1, 4
            do j6 = 1, 4
               colconv(j6-1,j7-1) = 1000
            end do
         end do
         colconv(1,1) = 1
         colconv(2,2) = 2
         colconv(1,2) = 3
         colconv(1,3) = 4
         colconv(2,3) = 5
         colconv(2,1) = 6
         colconv(3,1) = 7
         colconv(3,2) = 8
         colconv(0,1) = 9
         colconv(0,2) = 10
         colconv(0,3) = 11
         colconv(1,0) = -9
         colconv(2,0) = -10
         colconv(3,0) = -11
         colconv(0,0) = 0
         do j7 = 1, 23
            do j6 = 1, 23
               cmbcol(j6-12,j7-12) = 0
            end do
         end do
         cmbcol(1,1) = 1
         cmbcol(2,2) = 1
         cmbcol(3,6) = 1
         cmbcol(6,3) = 1
         cmbcol(4,7) = 1
         cmbcol(7,4) = 1
         cmbcol(5,8) = 1
         cmbcol(8,5) = 1
         cmbcol(9,-9) = 1
         cmbcol(10,-10) = 1
         cmbcol(11,-11) = 1
         cmbcol(-9,9) = 1
         cmbcol(-10,10) = 1
         cmbcol(-11,11) = 1
         cmbcol(0,0) = 1
      endif
      do j8 = 1, 5
         do j7 = 1, 43
            do j6 = 1, 2
               fldptr(j6,j7,j8) = -9999
            end do
         end do
      end do
      nprt = 0                             ! number of external particle
      do j1 = 1, 10
         if (prcss__flv(j1) .eq. 1001) go to 2
         nprt = nprt + 1
      end do
    2 continue
      nphf = nprt/2
      even = mod(nprt,2)                         ! for optmization
c                     ! maximum number of different combinations of nprt
      nmom = 2**nprt - 1
      call rncpl
      posfld = 1
c           ! array OFFSH, POSFLD and POSCOL allows to find the position
      poscol = 1
      flvold = -99
      do j1 = 1, nprt
         if (prcss__flv(j1) .ne. flvold) then
            if (flvold .gt. 0) fldptr(1,flvold,1) = nfld
            nfld = 1
            flvold = prcss__flv(j1)
            fldptr(2,flvold,1) = poscol
         else
            nfld = nfld + 1
         endif
         p1 = 2**(j1 - 1)
         if (colconv(prcss__col(2*j1-1),prcss__col(2*j1)) .eq. 2) then
            offshcol(1-1+poscol) = p1
            offshcol(2-1+poscol) = posfld
            offshcol(3-1+poscol) = 1
            do j6 = 1, nlordof(prcss__flv(j1))
               offsh(j6-1+posfld) = 0.D0
            end do
            posfld = posfld + nlordof(prcss__flv(j1))
            poscol = poscol + 3
            nfld = nfld + 1
         endif
c                                       !OFFSH(posfld,....) contains the
         call initpol (prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1      prcss__mom, j1, offsh(posfld))
         offshcol(1-1+poscol) = p1
         offshcol(2-1+poscol) = posfld
         offshcol(3-1+poscol) = colconv(prcss__col(2*j1-1),prcss__col(2*
     1      j1))
         do j6 = 1, 4
            impul(j6,offshcol(poscol)) = prcss__mom(j6+4*(j1-1))
         end do
         posfld = posfld + nlordof(prcss__flv(j1))
         poscol = poscol + 3
         if (colconv(prcss__col(2*j1-1),prcss__col(2*j1)) .eq. 1) then
            offshcol(1-1+poscol) = p1
            offshcol(2-1+poscol) = posfld
            offshcol(3-1+poscol) = 2
            do j6 = 1, nlordof(prcss__flv(j1))
               offsh(j6-1+posfld) = 0
            end do
            posfld = posfld + nlordof(prcss__flv(j1))
            poscol = poscol + 3
            nfld = nfld + 1
         endif
      end do
      fldptr(1,flvold,1) = nfld
      do j6 = 1, 451
         tabmrgpr(j6) = tabintmrg(j6,nprc)
      end do
      do j6 = 1, tabintmrg(1,nprc)
         cp(j6) = cpfs(j6,nprc)
         op(j6) = opfs(j6,nprc)
      end do
      do j1 = 1, nvasfs(1,nprc)
         cp(nvasfs(j1+1,nprc)) = cp(nvasfs(j1+1,nprc))*rnas
      end do
      do j6 = 1, nmom
         psnew(j6) = 0
      end do
      do iter = 2, nphf
         sim = 0
         do j6 = 1, 44
            cnsold(j6-1) = .FALSE.
         end do
         steptm = 2                              !for optimization
         flvold = -99
         do j1 = 1, tabmrgpr(1)
c                  !labels the flavour of the daughter offshell particle
            fl1 = tabmrgpr(steptm)
c                 !fl2 and fl3 label the flavour of the parents offshell
            fl2 = tabmrgpr(steptm+1)
            fl3 = tabmrgpr(steptm+2)             !particles to be merged
            coupot = cp(j1)
            operot = op(j1)                      !oper(fl1,fl2,fl3)
            nl1 = nlordof(fl1)
            nl2 = nlordof(fl2)                   !degrees of freedom
            nl3 = nlordof(fl3)
            if (fl1 .ne. flvold) then
               if (flvold .gt. 0) fldptr(1,chcg(flvold),iter) = nconf
               flvold = fl1
               if (cnsold(cnschg(fl1))) then
                  do j6 = 1, nmom
                     psnew(j6) = 0               !for optimization
                  end do
                  do j6 = 1, 44
                     cnsold(j6-1) = .FALSE.
                  end do
               endif
               cnsold(cnschg(fl1)) = .TRUE.
               nconf = 0
               fldptr(2,chcg(fl1),iter) = poscol
            endif
            steptm = steptm + 3
            if (fl2 .ne. fl3) then
               if (sim .ne. 1) then
                  do j6 = 1, 13
                     exc(j6) = nexc(j6,iter)
                  end do
               endif
               sim = 1
            else
               if (sim .ne. 2) then
                  do j6 = 1, 13
                     exc(j6) = nexcsm(j6,iter)
                  end do
               endif
               sim = 2
            endif
            stepex = 2
            do j2 = 1, exc(1)
               ex2 = exc(stepex)
               ex3 = exc(stepex+1)
               stepex = stepex + 2
               nf2 = fldptr(1,fl2,ex2)   !number of excitation of field2
               if (nf2 .gt. 0) then
                  nf3 = fldptr(1,fl3,ex3)
                  if (nf3 .gt. 0) then
                     pscfl2 = fldptr(2,fl2,ex2)
c                         !starting position of field2 in array OFFSHCOL
                     pscfl3 = fldptr(2,fl3,ex3)
c                            !starting position of field2 in array OFFSH
                     psfl2 = offshcol(pscfl2+1)
c                            !starting position of field3 in array OFFSH
                     psfl3 = offshcol(pscfl3+1)
                     do j3 = 1, nf2
c                    !if the flavour of field 2 and 3 are the same don't
                        if (ex2.eq.ex3 .and. fl2.eq.fl3) then
c                                     !repeat twice the same computation
                           beg3 = j3 + 1
                           psfl3i = psfl3 + nl3*(beg3 - 2)
                           pscfl3i = pscfl3 + 3*(beg3 - 2)
                        else
                           beg3 = 1
                           psfl3i = psfl3 - nl3
                           pscfl3i = pscfl3 - 3
                        endif
                        p2 = offshcol(pscfl2)
                        cl(1) = offshcol(pscfl2+2)
                        do j6 = 1, 4
                           imp2(j6) = impul(j6,p2)
                        end do
                        choose = iter.eq.nphf .and. even.eq.0
                        if (prcss__nfrm .gt. 2) then
                           nfr2 = 0
                           do j4 = 0, prcss__nfrm - 1
                              if (btest(p2,j4)) nfr2 = nfr2 + 1
                           end do
                        else
                           nfr2 = 0
                        endif
                        do j4 = beg3, nf3
                           pscfl3i = pscfl3i + 3
c                            !starting position of field3 in array OFFSH
                           psfl3i = offshcol(pscfl3i+1)
                           p3 = offshcol(pscfl3i)
c           !momentum of field1: resulting out of the merger of field2,3
                           p1 = p2 + p3
                           if (momd(p1) .eq. iter) then
c               !for even number of particle we accept merger of npart/2
                              if (.not.(choose .and. btest(p1+1,0))) 
     1                           then
                                 cl(2) = offshcol(pscfl3i+2)
                                 if (colrep(fl1).eq.2 .or. colrep(fl1)
     1                              .eq.3) then  !quarks and gluons
c  !pointer to deal with a one dimensional array as two dimensional, opt
                                    ptr = (cl(1)+12) + 23*(cl(2)+11)
c          !number of new coulored object: 1,2 2 only for neutral gluons
                                    nc = tbnew(ptr)
                                    if (nc .eq. 0) go to 4
                                    ptr = 2*ptr - 1
                                    colcff(1) = tbcoeff(ptr)
c              ! colcff=(/tbcoeff(ptr),tbcoeff(ptr+1)/)  !coulor clebsch
                                    colcff(2) = tbcoeff(ptr+1)
                                    clnw(1) = tbmom(ptr)
c!      clnw=tbmom(ptr:ptr+1)                   !coulor of the new objec
                                    clnw(2) = tbmom(ptr+1)
c                                                !coulorless object
                                 else if (colrep(fl1) .eq. 1) then
                                    if (cmbcol(cl(1),cl(2)) .eq. 0) 
     1                                 go to 4
                                    nc = 1
                                    colcff(1) = 1.D0
c                                                ! colcff=(/1.,0./)
                                    colcff(2) = 0.D0
                                    clnw(1) = 0
                                    clnw(2) = 0  ! clnw=0
                                 else
                                    write (*, *) 
     1'something wrong in color assignment'
                                 endif
                                 do j6 = 1, 4
                                    imp3(j6) = impul(j6,p3)
c                                            !momentum of the new fields
                                    imp(j6) = imp2(j6) + imp3(j6)
c                                   !array to store momenta combinations
                                    impul(j6,p1) = imp(j6)
                                 end do
                                 j9 = 0
                                 do j6 = 1, 4
                                    j9 = j9 + 1
                                    d1(j9) = -imp(j6)
                                 end do
c                           !via the proper trilinear interaction (oper)
                                 call fuse (d1, imp2, offsh(psfl2), imp3
     1                              , offsh(psfl3i), operot, fusion)
c                                !multiplying by the propagator -> polar
                                 call prpgt (chcg(fl1), imp, fusion, 
     1                              polar)
                                 if (nfr2 .gt. 0) then
c                             !returning proper fermi statistics in perm
                                    call permtn (prcss__nfrm, p2, p3, 
     1                                 perm)
                                 else
                                    perm = 1.D0
                                 endif
c                        !returning the position, in OFFSHCOL, of field1
                                 padd = psnew(p1)
c                                                !new configuration
                                 if (padd .eq. 0) then
                                    psnew(p1) = poscol
                                    if (diagonal) then
                                    if (clnw(1) .eq. 1) then
                                    nc = 2
                                    colcff(2) = 0.D0
                                    clnw(2) = 2
                                    endif
                                    endif
                                    do j5 = 1, nc
                                    cpot1 = coupot*colcff(j5)*perm
                                    do j6 = 1, nl1
c                                                !storing field1
                                    offsh(j6-1+posfld) = cpot1*polar(j6)
                                    end do
                                    offshcol(1-1+poscol) = p1
                                    offshcol(2-1+poscol) = posfld
                                    offshcol(3-1+poscol) = clnw(j5)
c    !number of new excitation of the same flavour and number of momenta
                                    nconf = nconf + 1
                                    poscol = poscol + 3
                                    posfld = posfld + nl1
                                    end do
c                              !new contribution to an old configuration
                                 else
                                    j5b = 0
                                    if (diagonal) then
                                    if (clnw(1) .eq. 2) j5b = nl1
                                    endif
                                    do j5 = 1, nc
                                    cpot1 = coupot*colcff(j5)*perm
c                                           !position of field1 in OFFSH
                                    padd0 = offshcol(padd+1) + j5b
                                    do j6 = 1, nl1
c                                           !adding the new contribution
                                    offsh(j6-1+padd0) = offsh(j6-1+padd0
     1                                 ) + cpot1*polar(j6)
                                    end do
                                    j5b = j5b + nl1
                                    end do
                                 endif
                              endif
                           endif
    4                      continue
                        end do
                        pscfl2 = pscfl2 + 3
                        psfl2 = offshcol(pscfl2+1)
                     end do
                  endif
               endif
            end do
         end do
         fldptr(1,chcg(fl1),iter) = nconf
      end do
      if (fststp.eq.1 .and. dbg) then
         fststp = 2
         if (poscol .gt. pscit) then
            pscit = poscol
            write (*, *) 'iteration,poscol,posfld', poscol, posfld
         endif
         if (posfld .gt. psfit) then
            psfit = posfld
            write (*, *) 'iteration,poscol,posfld', poscol, posfld
         endif
      endif
      if (poscol.gt.dim(2) .or. posfld.gt.dim(1)) then
         write (*, *) 
     1      'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
         write (*, *) 'iteration,poscol,posfld', poscol, posfld, dim(2)
     1      , dim(1)
         stop 
      endif
      amplitude = (0.,0.)
      do j6 = 1, 151
         tabmrgpr(j6) = tabint(j6,nprc)
      end do
      do j6 = 1, tabint(1,nprc)
         cp(j6) = cpam(j6,nprc)
         op(j6) = opam(j6,nprc)
      end do
      do j1 = 1, nvasam(1,nprc)
         cp(nvasam(j1+1,nprc)) = cp(nvasam(j1+1,nprc))*rnas
      end do
      do j6 = 1, 44
         cnsold(j6-1) = .FALSE.
      end do
      do j6 = 1, nmom
         psnew(j6) = 0    !none of these mergin has been computed before
c!reporting the order of field1 momenta to avoid repeating the same comp
         ordtmp(j6) = 0
      end do
      steptm = 2
      flvold = -99
      do j1 = 1, tabmrgpr(1) + 1
         if (j1 .eq. tabmrgpr(1)+1) then
c    !last step compute the final contribution to the amplitude and exit
            flvold = 9999
         else
            fl1 = tabmrgpr(steptm)
            fl2 = tabmrgpr(steptm+1)
            fl3 = tabmrgpr(steptm+2)
            steptm = steptm + 3
            nl1 = nlordof(fl1)
            nl2 = nlordof(fl2)
            nl3 = nlordof(fl3)
            coupot = cp(j1)
            operot = op(j1)                      !oper(fl1,fl2,fl3)
         endif
         if (fl1 .ne. flvold) then
            if (fststp.eq.2 .and. j1.ne.1) then
               if (poscol.gt.dim(4) .or. posfld.gt.dim(3)) then
                  write (*, *) 'dim', dim
                  write (*, *) 
     1               'ill dimensioned array AMPCOL and/or AMPFLD in AMP'
                  write (*, *) 'amplitude,poscol,posfld', poscol, posfld
     1               , dim(4), dim(3)
                  stop 
               endif
               fststp = 2
               if (dbg) then
                  if (poscol .gt. psca) then
                     psca = poscol
                     write (*, *) 'amplitude,poscol,posfld', poscol, 
     1                  posfld
                  endif
                  if (posfld .gt. psfa) then
                     psfa = posfld
                     write (*, *) 'amplitude,poscol,posfld', poscol, 
     1                  posfld
                  endif
               endif
            endif
            if (flvold .gt. 0) then
c                                         !flvold=9999 last contribution
               if (flvold .eq. 9999) flvold = fl1
               cnsold(cnschg(flvold)) = .TRUE.
c                                  !number of lorentz degrees of freedom
               nlold = nlordof(flvold)
               do iter = 1, nphf
                  nf1 = fldptr(1,flvold,iter)    !same as in iteration
                  if (nf1 .gt. 0) then
                     pscfl1 = fldptr(2,flvold,iter)
                     psfl1 = offshcol(pscfl1+1)
                     pref = 0                    !for optimization
                     do j2 = 1, nf1
                        psfl1 = offshcol(pscfl1+1)
                        p1 = offshcol(pscfl1)
                        p1 = nmom - p1
                        padd = psnew(p1)    !locating position of field1
c                                 !no 2,3 merging compatible with field1
                        if (padd .eq. 0) then
                           pscfl1 = pscfl1 + 3
                           go to 6
                        endif
                        if (p1 .ne. pref) then
                           pref = p1
                        else
                           padd = padd + 3 !there are two neutral gluons
                        endif
                        cl(1) = offshcol(pscfl1+2)
                        cl(2) = ampcol(padd+2)
                        if (cmbcol(cl(1),cl(2)) .eq. 1) then
                           padd = ampcol(padd+1)
                           x1 = 0
                           do j6 = 1, nlold
                              x1 = x1 + offsh(j6-1+psfl1)*ampfld(j6-1+
     1                           padd)
                           end do
                           tmp = x1
                           if (prcss__nfrm .gt. 2) then
c                             !returning proper fermi statistics in perm
                              call permtn(prcss__nfrm,nmom-p1,p1,perm)
                           else
                              perm = 1.D0
                           endif
                           amplitude = amplitude + tmp*perm
                           pscfl1 = pscfl1 + 3
                        endif
    6                   continue
                     end do
                  endif
               end do
               if (j1 .eq. tabmrgpr(1)+1) go to 7!end of computation
               if (cnsold(cnschg(fl1))) then
                  do j6 = 1, nmom
c                         !none of these mergin has been computed before
                     psnew(j6) = 0
c!reporting the order of field1 momenta to avoid repeating the same comp
                     ordtmp(j6) = 0
                  end do
               endif
            endif
            flvold = fl1
            posfld = 1                           !pointer to offsh
            poscol = 1                           !pointer to offshcol
         endif
         do iter = 1, nphf
            npmit = nprt - iter
            choose2 = iter.eq.nphf .and. even.eq.0
            nf1 = fldptr(1,fl1,iter)
            if (nf1 .gt. 0) then
               pscfl1 = fldptr(2,fl1,iter)
               if (fl1 .eq. fl2) then
                  beg = iter
               else
                  beg = 1
               endif
               do ex2 = beg, nphf
                  choose = fl1.eq.fl2 .and. iter.eq.ex2
                  ex3 = nprt - iter - ex2
                  if (ex3.ge.1 .and. ex3.le.nphf) then
c                                      !again to deal with equal flavour
                     if (fl2.ne.fl3 .or. ex3.ge.ex2) then
                        nf2 = fldptr(1,fl2,ex2)
                        if (nf2 .gt. 0) then
                           nf3 = fldptr(1,fl3,ex3)
                           if (nf3 .gt. 0) then
                              pscfl2 = fldptr(2,fl2,ex2)
                              pscfl3 = fldptr(2,fl3,ex3)
                              choose3 = ex2.eq.ex3 .and. fl2.eq.fl3
                              j5 = pscfl1
                              do j4 = 1, nf1
c                 !now ordtmp register the ordering of momenta of field1
                                 ordtmp(offshcol(j5)) = j4
                                 j5 = j5 + 3
                              end do
                              do j3 = 1, nf2
                                 cl(1) = offshcol(pscfl2+2)
                                 psfl2 = offshcol(pscfl2+1)
                                 if (choose3) then
                                    beg3 = j3 + 1
                                    pscfl3i = pscfl3 + 3*(beg3 - 2)
                                 else
                                    beg3 = 1
                                    pscfl3i = pscfl3 - 3
                                 endif
                                 p2 = offshcol(pscfl2)
                                 do j6 = 1, 4
                                    imp2(j6) = impul(j6,p2)
                                 end do
                                 if (prcss__nfrm .gt. 2) then
                                    nfr2 = 0
                                    do j4 = 0, prcss__nfrm - 1
                                    if (btest(p2,j4)) nfr2 = nfr2 + 1
                                    end do
                                 else
                                    nfr2 = 0
                                 endif
                                 do j4 = beg3, nf3
                                    pscfl3i = pscfl3i + 3
                                    p3 = offshcol(pscfl3i)
                                    p1 = p2 + p3
                                    if (momd(p1) .eq. npmit) then
                                    ot = ordtmp(nmom-p1)
c!field1 does not contain the frequency to match the proposed 2,3 fusion
                                    if (ot .ne. 0) then
c                         !to avoid repeating the same computation twice
                                    if (.not.(choose .and. ordtmp(p2)
     1                                 .le.ot)) then
                                    if (.not.(choose2 .and. btest(p1,0))
     1                                 ) then
c                                        ! cl=(/ci,offshcol(pscfl3i+2)/)
                                    cl(2) = offshcol(pscfl3i+2)
                                    if (colrep(fl1).eq.2 .or. colrep(fl1
     1                                 ).eq.3) then
                                    ptr = (cl(1)+12) + 23*(cl(2)+11)
                                    nc = tbnew(ptr)
                                    if (nc .eq. 0) go to 10
                                    ptr = 2*ptr - 1
                                    do j6 = 1, 2
                                    colcff(j6) = tbcoeff(j6-1+ptr)
                                    clnw(j6) = tbmom(j6-1+ptr)
                                    end do
                                    else if (colrep(fl1) .eq. 1) then
                                    if (cmbcol(cl(1),cl(2)) .eq. 0) 
     1                                 go to 10
                                    nc = 1
                                    colcff(1) = 1.D0
c                                                !  colcff=(/1.,0./)
                                    colcff(2) = 0.D0
                                    clnw(1) = 0
                                    clnw(2) = 0  !  clnw=0
                                    else
                                    write (*, *) 
     1'something wrong in color assignment'
                                    endif
                                    do j6 = 1, 4
                                    imp3(j6) = impul(j6,p3)
                                    imp(j6) = imp2(j6) + imp3(j6)
                                    end do
                                    psfl3i = offshcol(pscfl3i+1)
                                    j10 = 0
                                    do j6 = 1, 4
                                    j10 = j10 + 1
                                    d2(j10) = -imp(j6)
                                    end do
                                    call fuse (d2, imp2, offsh(psfl2), 
     1                                 imp3, offsh(psfl3i), operot, 
     2                                 fusion)
                                    if (nfr2 .gt. 0) then
                                    call permtn (prcss__nfrm, p2, p3, 
     1                                 perm)
                                    else
                                    perm = 1.D0
                                    endif
                                    padd = psnew(p1)
c                                                !new configuration
                                    if (padd .eq. 0) then
                                    psnew(p1) = poscol
                                    if (diagonal) then
                                    if (clnw(1) .eq. 1) then
                                    nc = 2
                                    colcff(2) = 0.D0
                                    clnw(2) = 2
                                    endif
                                    endif
                                    do j5 = 1, nc
                                    cpot1 = coupot*colcff(j5)*perm
                                    do j6 = 1, nl1
                                    ampfld(j6-1+posfld) = cpot1*fusion(
     1                                 j6)
                                    end do
                                    ampcol(1-1+poscol) = p1
                                    ampcol(2-1+poscol) = posfld
                                    ampcol(3-1+poscol) = clnw(j5)
                                    nconf = nconf + 1
                                    posfld = posfld + nl1
                                    poscol = poscol + 3
                                    end do
c                              !new contribution to an old configuration
                                    else
                                    j5b = 0
                                    if (diagonal) then
                                    if (clnw(1) .eq. 2) j5b = nl1
                                    endif
                                    do j5 = 1, nc
                                    padd0 = ampcol(padd+1) + j5b
                                    cpot1 = coupot*colcff(j5)*perm
                                    do j6 = 1, nl1
                                    ampfld(j6-1+padd0) = ampfld(j6-1+
     1                                 padd0) + cpot1*fusion(j6)
                                    end do
                                    j5b = j5b + nl1
                                    end do
                                    endif
                                    endif
                                    endif
                                    endif
                                    endif
   10                               continue
                                 end do
                                 pscfl2 = pscfl2 + 3
                              end do
                           endif
                        endif
                     endif
                  endif
               end do
            endif
         end do
    7    continue
      end do
      mtel = amplitude
      fststp = 1
!
      end 
!
      subroutine initpol(prcss__flv, prcss__hel, prcss__col, prcss__nfrm
     1   , prcss__mom, prt, fld)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexcmax
      parameter (nexcmax = 13)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
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
      common /extsrc/src, exflg
      double complex src(6)
      logical exflg
      integer prcss__flv(10)                     !input   process
      integer prcss__hel(10)
      integer prcss__col(20)
      integer prcss__nfrm
      double precision prcss__mom(40)
      integer prt                              !prt-th to be initialized
      double complex fld(6)  !polarization vector of particle number prt
!
      integer tmp
      double precision im(4)
      double complex tmpc(4)
!
      integer maxpar
      parameter (maxpar = 10)
      integer vcnt
      character idecay
      common /tdecflag/idecay
      double complex vtp(4), utpb(4)
      common /tspinors/vtp, utpb
      double complex eps(4,10)
      common /vdec/eps
      integer j1
!
      if (prt .eq. 1) vcnt = 0
!
      tmp = prcss__hel(prt)
      do j1 = 1, 4
         im(j1) = prcss__mom(j1+4*(prt-1))
      end do
      if (im(1) .lt. 0) then
         do j1 = 1, 4
            im(j1) = -im(j1)
         end do
      endif
      if (lortyp(prcss__flv(prt)).eq.1 .or. lortyp(prcss__flv(prt)).eq.
     1   11) then                                !gluon source
!      if (prt.eq.1) then; fld(1:3)=im(2:4); fld(4:6)=0.d0; return; endif
         call sourcemassboson (tmp, tmpc, im, masses(prcss__flv(prt)))
         do j1 = 1, 3
            fld(j1) = tmpc(j1+1)
            fld(j1+3) = 0.D0
         end do
c              !fermion source (incoming fermion outcoming anti fermion)
      else if (lortyp(prcss__flv(prt)) .eq. 3) then
         if (prcss__flv(prt).eq.20 .and. idecay.eq.'y') then
            do j1 = 1, 4
               fld(j1) = vtp(j1)
            end do
         else
            call fermionsources (tmp, fld, im, masses(prcss__flv(prt)))
         endif
c                                                !fermion bar source
      else if (lortyp(prcss__flv(prt)) .eq. 4) then
         if (prcss__flv(prt).eq.36 .and. idecay.eq.'y') then
            do j1 = 1, 4
               fld(j1) = utpb(j1)
            end do
         else
            call fermionbarsources(tmp,fld,im,masses(prcss__flv(prt)))
         endif
c              !fermion source (incoming fermion outcoming anti fermion)
      else if (lortyp(prcss__flv(prt)) .eq. 5) then
         if (tmp .lt. 0) then
            write (*, *) 'wrong helicity assignment', tmp
            stop 
         endif
         call fermionsourceshl (tmp, fld, im)
c                                                !fermion bar source
      else if (lortyp(prcss__flv(prt)) .eq. 6) then
         if (tmp .lt. 0) then
            write (*, *) 'wrong helicity assignment', tmp
            stop 
         endif
         call fermionbarsourceshl (tmp, fld, im)
c              !fermion source (incoming fermion outcoming anti fermion)
      else if (lortyp(prcss__flv(prt)) .eq. 7) then
         if (tmp .gt. 0) then
            write (*, *) 'wrong helicity assignment', tmp
            stop 
         endif
         call fermionsourceshl (tmp, fld, im)
c                                                !fermion bar source
      else if (lortyp(prcss__flv(prt)) .eq. 8) then
         if (tmp .gt. 0) then
            write (*, *) 'wrong helicity assignment', tmp
            stop 
         endif
         call fermionbarsourceshl (tmp, fld, im)
      else if (lortyp(prcss__flv(prt)) .eq. 2) then
         if (exflg) then                         !w,z source
            do j1 = 1, 4
               fld(j1) = src(j1)
            end do
            fld(5) = 0.D0
         else
            if (idecay .eq. 'y') then
               vcnt = vcnt + 1
               do j1 = 1, 4
                  tmpc(j1) = eps(j1,vcnt)
               end do
            else
               call sourcemassboson (tmp, tmpc, im, masses(prcss__flv(
     1            prt)))
            endif
            do j1 = 1, 4
               fld(j1) = tmpc(j1)
            end do
            fld(5) = 0.D0
         endif
      else if (lortyp(prcss__flv(prt)) .eq. 12) then
         fld(1) = 1.D0
         do j1 = 1, 2
            fld(j1+1) = 0.D0
         end do
      else
         write (6, *) 'wrong lorenz type in INITPOL', lortyp(prcss__flv(
     1      prt)), prcss__flv(prt)
         stop 
      endif
!
      end 
!
      subroutine fuse(p1, p2, amp2, p3, amp3, oper, nwam)
      INCLUDE 'Vincnst.inc'
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      common /naming/prtcl__chcg, prtcl__cnschg, prtcl__lortyp, 
     1   prtcl__nlordof, prtcl__colrep, prtcl__zcoup, chcg, cnschg, conv
      integer wm
      parameter (wm = 4)
      integer wp
      parameter (wp = 3)
      integer z
      parameter (z = 2)
      integer prtcl__chcg(43)
      integer prtcl__cnschg(43)
      integer prtcl__lortyp(43)
      integer prtcl__nlordof(43)
      integer prtcl__colrep(43)
      integer prtcl__zcoup(43)
      integer chcg(43)
      integer cnschg(43)
      integer conv(43)
      common /couplings/masses, width, coup, gstrong, g2weak, gbar, ctw
     1   , stw, glup, gldn, gllp, grup, grdn, grlp, emch, trihcp, yuk, 
     2   qrtc, qrtca, qrtc2, qrtc3, qrtchh, qrtchg, ggh, oper1x
      double precision masses(43)
      double precision width(43)
      integer oper1x(43,43,43)
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
      common /utilities/fl1, fl2, fl3, diagonal
      integer fl1
      integer fl2
      integer fl3
      logical diagonal
c                                  !on/off-shell amplitudes to be merged
      double complex amp2(6), amp3(6)
      double precision p1(4), p2(4), p3(4)       !particle momenta
c       !labelling the interaction type and the flavour of the 3rd field
      integer oper
      double complex nwam(6)                     !result of the merging
!
      double precision cp
!
      diagonal = .FALSE.
!
      if (oper .eq. 1) then
c                                           !trilinear gluon interaction
         call triglu (p1(2), p2(2), amp2, p3(2), amp3, nwam)
      else if (oper .eq. 4) then
         call fbftog (amp2, amp3, nwam)   !fusion of fermions into gluon
      else if (oper .eq. 5) then
c                                     !fusion of f-bar  gluon into f-bar
         call gfbtofb (amp2, amp3, nwam)
      else if (oper .eq. 6) then
         call gftof (amp2, amp3, nwam)!fusion of f-bar  gluon into f-bar
      else if (oper .eq. 7) then
         call fbftow (amp2, amp3, nwam)   !fusion of fermions into gluon
      else if (oper .eq. 8) then
c                                     !fusion of f-bar  gluon into f-bar
         call wfbtofb (amp2, amp3, nwam)
      else if (oper .eq. 9) then
         call wftof (amp2, amp3, nwam)!fusion of f-bar  gluon into f-bar
      else if (oper .eq. 10) then
         call fbftogl (amp2, amp3, nwam)  !fusion of fermions into gluon
      else if (oper .eq. 11) then
c                                     !fusion of f-bar  gluon into f-bar
         call gfbtofbl (amp2, amp3, nwam)
      else if (oper .eq. 12) then
c                                     !fusion of f-bar  gluon into f-bar
         call gftofl (amp2, amp3, nwam)
      else if (oper .eq. 13) then
         call fbftowl (amp2, amp3, nwam)  !fusion of fermions into gluon
      else if (oper .eq. 14) then
c                                     !fusion of f-bar  gluon into f-bar
         call wfbtofbl (amp2, amp3, nwam)
      else if (oper .eq. 15) then
c                                     !fusion of f-bar  gluon into f-bar
         call wftofl (amp2, amp3, nwam)
      else if (oper .eq. 16) then
         call fbftogr (amp2, amp3, nwam)  !fusion of fermions into gluon
      else if (oper .eq. 17) then
c                                     !fusion of f-bar  gluon into f-bar
         call gfbtofbr (amp2, amp3, nwam)
      else if (oper .eq. 18) then
c                                     !fusion of f-bar  gluon into f-bar
         call gftofr (amp2, amp3, nwam)
      else if (oper .eq. 19) then
         call fbftowr (amp2, amp3, nwam)  !fusion of fermions into gluon
      else if (oper .eq. 20) then
c                                     !fusion of f-bar  gluon into f-bar
         call wfbtofbr (amp2, amp3, nwam)
      else if (oper .eq. 21) then
c                                     !fusion of f-bar  gluon into f-bar
         call wftofr (amp2, amp3, nwam)
      else if (oper .eq. 22) then
c                                         !fusion of fermions into gluon
         call fbftoz (fl3, amp2, amp3, nwam)
      else if (oper .eq. 23) then
c                                     !fusion of f-bar  gluon into f-bar
         call zfbtofb (fl3, amp2, amp3, nwam)
      else if (oper .eq. 24) then
c                                     !fusion of f-bar  gluon into f-bar
         call zftof (fl3, amp2, amp3, nwam)
      else if (oper .eq. 25) then
c                                           !trilinear gluon interaction
         call triw (qrtc(1,fl1-1), p1, p2, amp2, p3, amp3, nwam)
      else if (oper .eq. 26) then
         call wwtoax (amp2, amp3, nwam)          !fusion W+ W+ into X++
      else if (oper .eq. 27) then
         call waxtow (amp2, amp3, nwam)          !fusion W+ Y++ into W+
      else if (oper .eq. 28) then
c                                                !fusion W+ Y++ into W+
         call wax2tow (qrtc2(1,fl1-1), amp2, amp3, nwam)
      else if (oper .eq. 29) then
c                                                !fusion W+ Y++ into W+
         call wwtoax2 (qrtc3(1,fl2-1), amp2, amp3, nwam)
      else if (oper .eq. 30) then
c                                             ! fusion W+ W- into photon
         call triphw (qrtca(1,fl1-1), p1, p2, amp2, p3, amp3, nwam)
      else if (oper .eq. 31) then
c                                               ! fusion W photot into W
         call triwph (qrtca(1,fl1-1), p1, p2, amp2, p3, amp3, nwam)
      else if (oper .eq. 32) then
c                                                !fusion W+ Y++ into W+
         call wax2top (qrtc2(1,fl1-1), amp2, amp3, nwam)
      else if (oper .eq. 33) then
c                                                !fusion W+ Y++ into W+
         call pax2tow (qrtc2(1,fl1-1), amp2, amp3, nwam)
      else if (oper .eq. 34) then
c                                                !fusion W+ Y++ into W+
         call pax2top (qrtc2(1,fl1-1), amp2, amp3, nwam)
      else if (oper .eq. 35) then
c                                                !fusion W+ Y++ into W+
         call zptoax2 (qrtc3(1,fl2-1), amp2, amp3, nwam)
      else if (oper .eq. 36) then
c                                                !fusion W+ Y++ into W+
         call pptoax2 (qrtc3(1,fl2-1), amp2, amp3, nwam)
      else if (oper .eq. 37) then
         if (fl2.eq.3 .or. fl2.eq.4) then
            cp = qrtchg(1)
         else if (fl2 .eq. 2) then
            cp = qrtchg(2)
         else
            write (*, *) 'wrong flavour in fuse'
            stop 
         endif
         call wwtoh (cp, amp2, amp3, nwam)       !fusion W+ W+ into X++
      else if (oper .eq. 38) then
         if (fl2.eq.3 .or. fl2.eq.4) then
            cp = qrtchg(1)
         else if (fl2 .eq. 2) then
            cp = qrtchg(2)
         else
            write (*, *) 'wrong flavour in fuse'
            stop 
         endif
         call whtow (cp, amp2, amp3, nwam)       !fusion W+ W+ into X++
      else if (oper .eq. 39) then
         call trih (amp2, amp3, nwam)
      else if (oper .eq. 40) then
         call fbftoh (amp2, amp3, nwam)
      else if (oper .eq. 41) then
         call fhtof (amp2, amp3, nwam)
      else if (oper .eq. 45) then
         call gluhtoglu (p1, p2, amp2, amp3, nwam)
         diagonal = .TRUE.
      else if (oper .eq. 46) then
         call gluglutoh (p2, amp2, p3, amp3, nwam)
      else if (oper .eq. 47) then
         call hxglutoxglu (amp2, amp3, nwam)
         diagonal = .TRUE.
      else if (oper .eq. 48) then
         call xgluxglutoh (amp2, amp3, nwam)
      else if (oper .eq. 49) then
         call gluyglutoglu (amp2, amp3, nwam)
      else if (oper .eq. 50) then
         call gluglutoyglu (amp2, amp3, nwam)
      else if (oper .eq. 51) then
         call hxglutoglu (p1, amp2, amp3, nwam)
         diagonal = .TRUE.
      else if (oper .eq. 52) then
         call gluxglutoh (p2, amp2, amp3, nwam)
      else if (oper .eq. 53) then
         call gluhtoxglu (p2, amp2, amp3, nwam)
         diagonal = .TRUE.
      else
         write (*, *) 'WRONG lorentz operator selection', oper
         stop 
      endif
!
      end 
!
      subroutine gluyglutoglu(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = -(amp2(2)*amp3(3)-amp2(3)*amp3(2))
      nwam(2) = -(amp2(3)*amp3(1)-amp2(1)*amp3(3))
      nwam(3) = -(amp2(1)*amp3(2)-amp2(2)*amp3(1))
      do j1 = 1, 3
         nwam(j1+3) = 0.D0
      end do
!
      end 
!
      subroutine gluglutoyglu(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = amp3(2)*amp2(3) - amp3(3)*amp2(2)
      nwam(2) = amp3(3)*amp2(1) - amp3(1)*amp2(3)
      nwam(3) = amp3(1)*amp2(2) - amp3(2)*amp2(1)
!
      end 
!
      subroutine hxglutoglu(p1, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision p1(4)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
      do j1 = 1, 3
         nwam(j1+3) = 0.D0
      end do
      nwam(1) = -(amp3(2)*p1(4)-amp3(3)*p1(3))
      nwam(2) = -(amp3(3)*p1(2)-amp3(1)*p1(4))
      nwam(3) = -(amp3(1)*p1(3)-amp3(2)*p1(2))
      do j1 = 1, 3
         nwam(j1) = nwam(j1)*amp2(1)
      end do
!
      end 
!
      subroutine gluxglutoh(p2, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision p2(4)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = -amp3(1)*(p2(3)*amp2(3)-p2(4)*amp2(2))
      nwam(1) = nwam(1) - amp3(2)*(p2(4)*amp2(1)-p2(2)*amp2(3))
      nwam(1) = nwam(1) - amp3(3)*(p2(2)*amp2(2)-p2(3)*amp2(1))
      do j1 = 1, 2
         nwam(j1+1) = (0.D0,0.D0)
      end do
!
      end 
!
      subroutine gluhtoxglu(p2, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision p2(4)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = (-p2(3)*amp2(3)) + p2(4)*amp2(2)
      nwam(2) = (-p2(4)*amp2(1)) + p2(2)*amp2(3)
      nwam(3) = (-p2(2)*amp2(2)) + p2(3)*amp2(1)
      do j1 = 1, 3
         nwam(j1) = nwam(j1)*amp3(1)
      end do
!
      end 
!
      subroutine xgluxglutoh(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1)*amp3(j1)
      end do
      nwam(1) = x1
      do j1 = 1, 2
         nwam(j1+1) = (0.D0,0.D0)
      end do
!
      end 
!
      subroutine hxglutoxglu(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
      do j1 = 1, 3
         nwam(j1) = amp3(j1)
         nwam(j1) = amp2(1)*nwam(j1)
      end do
!
      end 
!
      subroutine gluhtoglu(p1, p2, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision p1(4), p2(4)
      double complex nwam(6)
!
      double complex tmp
      doubleprecision d1
      doublecomplex x1
      integer j1
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + p1(j1+1)*p2(j1+1)
      end do
      tmp = p1(1)*p2(1) - d1
      do j1 = 1, 3
         nwam(j1) = -tmp*amp2(j1)
      end do
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + p1(j1+1)*amp2(j1)
      end do
      tmp = -x1
      do j1 = 1, 3
         nwam(j1) = nwam(j1) + tmp*p2(j1+1)
         nwam(j1) = nwam(j1)*amp3(1)
         nwam(j1+3) = 0.D0
      end do
!
      end 
!
      subroutine gluglutoh(p2, amp2, p3, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision p3(4), p2(4)
      double complex nwam(6)
      doubleprecision d1
      doublecomplex x1, x2, x3
      integer j1
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + p2(j1+1)*p3(j1+1)
      end do
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1)*amp3(j1)
      end do
      nwam(1) = (p2(1)*p3(1)-d1)*(-x1)
      x2 = 0
      do j1 = 1, 3
         x2 = x2 + p2(j1+1)*amp3(j1)
      end do
      x3 = 0
      do j1 = 1, 3
         x3 = x3 + amp2(j1)*p3(j1+1)
      end do
      nwam(1) = nwam(1) - x2*x3
      do j1 = 1, 2
         nwam(j1+1) = (0.D0,0.D0)
      end do
!
      end 
!
c                                       !trilinear gluon interaction   !
      subroutine triglu(p1, p2, amp2, p3, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision p1(3), p2(3), p3(3)
      double complex nwam(6)
!
      double complex a2a3, a2p, a3p
!
      integer j1
      doublecomplex x1, x2, x3
      integer j2
      x1 = 0
      do j2 = 1, 3
         x1 = x1 + amp2(j2)*amp3(j2)
      end do
      a2a3 = x1
      x2 = 0
      do j2 = 1, 3
         x2 = x2 + amp2(j2)*(p3(j2)-p1(j2))
      end do
      a2p = x2
      x3 = 0
      do j2 = 1, 3
         x3 = x3 + amp3(j2)*(p1(j2)-p2(j2))
      end do
      a3p = x3
      do j2 = 1, 3
         nwam(j2) = a2a3*(p2(j2)-p3(j2)) + a2p*amp3(j2) + a3p*amp2(j2)
      end do
      nwam(1) = nwam(1) + amp3(5)*amp2(3) - amp3(6)*amp2(2) - amp2(5)*
     1   amp3(3) + amp2(6)*amp3(2)
      nwam(2) = nwam(2) + amp3(6)*amp2(1) - amp3(4)*amp2(3) - amp2(6)*
     1   amp3(1) + amp2(4)*amp3(3)
      nwam(3) = nwam(3) + amp3(4)*amp2(2) - amp3(5)*amp2(1) - amp2(4)*
     1   amp3(2) + amp2(5)*amp3(1)
      nwam(4) = (-amp2(2)*amp3(3)) + amp2(3)*amp3(2)
      nwam(5) = (-amp2(3)*amp3(1)) + amp2(1)*amp3(3)
      nwam(6) = (-amp2(1)*amp3(2)) + amp2(2)*amp3(1)
!
      end 
!
c                                             !trilinear WWZ interaction
      subroutine triphw(qrt, p1, p2, amp2, p3, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(5), amp3(5)
      double precision p1(4), p2(4), p3(4)
      double precision qrt(3)
      double complex nwam(4)
!
      double complex a2a3, a2p, a3p, ax2, ax3
!
      integer j1
      doublecomplex x1, x2, x3
      integer j2
      x1 = 0
      do j2 = 1, 3
         x1 = x1 + amp2(j2+1)*amp3(j2+1)
      end do
      a2a3 = amp2(1)*amp3(1) - x1
      x2 = 0
      do j2 = 1, 3
         x2 = x2 + amp2(j2+1)*(p3(j2+1)-p1(j2+1))
      end do
      a2p = amp2(1)*(p3(1)-p1(1)) - x2
      x3 = 0
      do j2 = 1, 3
         x3 = x3 + amp3(j2+1)*(p1(j2+1)-p2(j2+1))
      end do
      a3p = amp3(1)*(p1(1)-p2(1)) - x3
      ax3 = qrt(3)*amp3(5)
      ax2 = qrt(2)*amp2(5)
      do j2 = 1, 3
         nwam(j2) = a2a3*(p2(j2+1)-p3(j2+1)) + (a2p + ax2)*amp3(j2+1) + 
     1      (a3p + ax3)*amp2(j2+1)
      end do
      nwam(4) = qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
      end 
!
c                                             !trilinear WWZ interaction
      subroutine triw(qrt, p1, p2, amp2, p3, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(5), amp3(5)
      double precision p1(4), p2(4), p3(4)
      double precision qrt(3)
      double complex nwam(5)
!
      double complex a2a3, a2p, a3p, ax2, ax3
!
      integer j1
      doublecomplex x1, x2, x3
      integer j2
      x1 = 0
      do j2 = 1, 3
         x1 = x1 + amp2(j2+1)*amp3(j2+1)
      end do
      a2a3 = amp2(1)*amp3(1) - x1
      x2 = 0
      do j2 = 1, 3
         x2 = x2 + amp2(j2+1)*(p3(j2+1)-p1(j2+1))
      end do
      a2p = amp2(1)*(p3(1)-p1(1)) - x2
      x3 = 0
      do j2 = 1, 3
         x3 = x3 + amp3(j2+1)*(p1(j2+1)-p2(j2+1))
      end do
      a3p = amp3(1)*(p1(1)-p2(1)) - x3
      ax3 = qrt(3)*amp3(5)
      ax2 = qrt(2)*amp2(5)
      do j2 = 1, 4
         nwam(j2) = a2a3*(p2(j2)-p3(j2)) + (a2p + ax2)*amp3(j2) + (a3p
     1       + ax3)*amp2(j2)
      end do
      nwam(1) = -nwam(1)
      nwam(5) = qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
      end 
!
c                                             !trilinear WWZ interaction
      subroutine triwph(qrt, p1, p2, amp2, p3, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(5)
      double complex amp3(4)
      double precision p1(4), p2(4), p3(4)
      double precision qrt(3)
      double complex nwam(5)
!
      double complex a2a3, a2p, a3p, ax2, ax3
!
      integer j1
      doublecomplex x1, x2, x3
      integer j2
      x1 = 0
      do j2 = 1, 3
         x1 = x1 + amp2(j2+1)*amp3(j2)
      end do
      a2a3 = -x1
      x2 = 0
      do j2 = 1, 3
         x2 = x2 + amp2(j2+1)*(p3(j2+1)-p1(j2+1))
      end do
      a2p = amp2(1)*(p3(1)-p1(1)) - x2
      x3 = 0
      do j2 = 1, 3
         x3 = x3 + amp3(j2)*(p1(j2+1)-p2(j2+1))
      end do
      a3p = -x3
      ax3 = qrt(3)*amp3(4)
      ax2 = qrt(2)*amp2(5)
      nwam(1) = a2a3*(p2(1)-p3(1)) + (a3p + ax3)*amp2(1)
      do j2 = 1, 3
         nwam(j2+1) = a2a3*(p2(j2+1)-p3(j2+1)) + (a2p + ax2)*amp3(j2) + 
     1      (a3p + ax3)*amp2(j2+1)
      end do
      nwam(1) = -nwam(1)
      nwam(5) = qrt(1)*a2a3
!  nwam(1:4)=(ax2)*amp3(1:4)+(ax3)*amp2(1:4)
!  nwam(1)=-nwam(1)
!  nwam(5)=qrt(1)*a2a3
!
      end 
!
      subroutine waxtow(amp2, amp3, nwam)      !fusion of W+ Y++ into W+
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
      do j1 = 1, 4
         nwam(j1) = amp3(1)*amp2(j1)
      end do
      do j1 = 1, 3
         nwam(j1+1) = -nwam(j1+1)
      end do
      nwam(5) = 0.D0
!
      end 
!
      subroutine pax2tow(qrt, amp2, amp3, nwam)!fusion of W+ Y++ into W+
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(6)
!
      double complex tmp
      integer j1
!
      tmp = qrt(1)*amp3(1) + qrt(2)*amp3(2)
      do j1 = 1, 3
         nwam(j1+1) = -tmp*amp2(j1)
      end do
!  nwam(1)=-nwam(1)
      nwam(1) = 0.D0
      nwam(5) = 0.D0
!
      end 
!
      subroutine pax2top(qrt, amp2, amp3, nwam)!fusion of W+ Y++ into W+
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(6)
!
      double complex tmp
      integer j1
!
      tmp = qrt(1)*amp3(1) + qrt(2)*amp3(2)
      do j1 = 1, 3
         nwam(j1) = -tmp*amp2(j1)
      end do
!  nwam(1)=-nwam(1)
      nwam(4) = 0.D0
!
      end 
!
      subroutine wax2top(qrt, amp2, amp3, nwam)!fusion of W+ Y++ into W+
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(6)
!
      double complex tmp
      integer j1
!
      tmp = qrt(1)*amp3(1) + qrt(2)*amp3(2)
      do j1 = 1, 3
         nwam(j1) = -tmp*amp2(j1+1)
      end do
      nwam(4) = 0.D0
!
      end 
!
      subroutine wax2tow(qrt, amp2, amp3, nwam)!fusion of W+ Y++ into W+
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(6)
!
      double complex tmp
      integer j1
!
      tmp = qrt(1)*amp3(1) + qrt(2)*amp3(2)
      do j1 = 1, 4
         nwam(j1) = tmp*amp2(j1)
      end do
!  nwam(1)=-nwam(1)
      do j1 = 1, 3
         nwam(j1+1) = -nwam(j1+1)
      end do
      nwam(5) = 0.D0
!
      end 
!
      subroutine wwtoax(amp2, amp3, nwam)      !fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1+1)*amp3(j1+1)
      end do
      nwam(1) = amp3(1)*amp2(1) - x1
!
      end 
!
      subroutine wwtoh(cp, amp2, amp3, nwam)   !fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double precision cp
      double complex amp2(4), amp3(4)
      double complex nwam(3)
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1+1)*amp3(j1+1)
      end do
      nwam(1) = amp3(1)*amp2(1) - x1
      nwam(2) = cp*nwam(1)
      nwam(3) = 0.D0
!
      end 
!
      subroutine whtow(cp, amp2, amp3, nwam)   !fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double precision cp
      double complex amp2(4), amp3(3)
      double complex nwam(5)
!
      double complex ax
      integer j1
!
      ax = amp3(1) + cp*amp3(2)
      nwam(1) = ax*amp2(1)
      do j1 = 1, 3
         nwam(j1+1) = -ax*amp2(j1+1)
      end do
      nwam(5) = 0.D0
!
      end 
!
      subroutine trih(amp2, amp3, nwam)        !fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      double complex amp2(3), amp3(3)
      double complex nwam(3)
!
      double complex ax
      integer j1
!
      ax = amp3(1)*amp2(1)
      do j1 = 1, 3
         nwam(j1) = ax*qrtchh(j1)
      end do
! nwam(1)=ax*qrtchh(1); nwam(2)=ax*qrtchh(3);  nwam(3)=ax*qrtchh(2);
      nwam(1) = nwam(1) + amp2(1)*(qrtchh(2)*amp3(2)+qrtchh(3)*amp3(3))
     1    + amp3(1)*(qrtchh(2)*amp2(2)+qrtchh(3)*amp2(3))
!
      end 
!
      subroutine fbftoh(amp2, amp3, nwam)      !fusion of W+ W+ into X++
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double complex amp2(4), amp3(4)
      double complex nwam(3)
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 4
         x1 = x1 + amp2(j1)*amp3(j1)
      end do
      nwam(1) = x1
      do j1 = 1, 2
         nwam(j1+1) = (0.D0,0.D0)
      end do
!
      end 
!
      subroutine fhtof(amp2, amp3, nwam)       !fusion of W+ W+ into X++
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double complex amp2(4), amp3(3)
      double complex nwam(4)
      integer j1
      do j1 = 1, 4
         nwam(j1) = amp3(1)*amp2(j1)
      end do
!
      end 
!
      subroutine wwtoax2(qrt, amp2, amp3, nwam)!fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(2)
!
      double complex tmp
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1+1)*amp3(j1+1)
      end do
      tmp = amp3(1)*amp2(1) - x1
      do j1 = 1, 2
         nwam(j1) = tmp*qrt(j1)
      end do
!
      end 
!
      subroutine zptoax2(qrt, amp2, amp3, nwam)!fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(2)
!
      double complex tmp
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1+1)*amp3(j1)
      end do
      tmp = -x1
      do j1 = 1, 2
         nwam(j1) = tmp*qrt(j1)
      end do
!
      end 
!
      subroutine pptoax2(qrt, amp2, amp3, nwam)!fusion of W+ W+ into X++
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double precision qrt(2)
      double complex nwam(2)
!
      double complex tmp
      doublecomplex x1
      integer j1
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + amp2(j1)*amp3(j1)
      end do
      tmp = -x1
      do j1 = 1, 2
         nwam(j1) = tmp*qrt(j1)
      end do
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftog(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = amp2(1)*amp3(4) + amp2(2)*amp3(3) - amp2(3)*amp3(2) - 
     1   amp2(4)*amp3(1)
      nwam(2) = (-amp2(1)*amp3(4)) + amp2(2)*amp3(3) + amp2(3)*amp3(2)
     1    - amp2(4)*amp3(1)
      nwam(2) = (0.D0,1.D0)*nwam(2)
      nwam(3) = amp2(1)*amp3(3) - amp2(2)*amp3(4) - amp2(3)*amp3(1) + 
     1   amp2(4)*amp3(2)
      do j1 = 1, 3
         nwam(j1+3) = 0
      end do
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftow(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double complex c11, c12, c21, c22
!
      c11 = amp2(3)*amp3(1)
      c12 = amp2(3)*amp3(2)
      c21 = amp2(4)*amp3(1)
      c22 = amp2(4)*amp3(2)
!
      nwam(1) = c11 + c22
      nwam(2) = c12 + c21
      nwam(3) = (0.D0,1.D0)*(c21 - c12)
      nwam(4) = c11 - c22
      nwam(5) = (0.D0,0.D0)
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftogr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = amp2(1)*amp3(2) + amp2(2)*amp3(1)
      nwam(2) = (-amp2(1)*amp3(2)) + amp2(2)*amp3(1)
      nwam(2) = (0.D0,1.D0)*nwam(2)
      nwam(3) = amp2(1)*amp3(1) - amp2(2)*amp3(2)
      do j1 = 1, 3
         nwam(j1+3) = 0
      end do
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftogl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      integer j1
!
      nwam(1) = (-amp2(1)*amp3(2)) - amp2(2)*amp3(1)
      nwam(2) = amp2(1)*amp3(2) - amp2(2)*amp3(1)
      nwam(2) = (0.D0,1.D0)*nwam(2)
      nwam(3) = (-amp2(1)*amp3(1)) + amp2(2)*amp3(2)
      do j1 = 1, 3
         nwam(j1+3) = 0
      end do
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftowl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double complex c11, c12, c21, c22
!
      c11 = amp2(1)*amp3(1)
      c12 = amp2(1)*amp3(2)
      c21 = amp2(2)*amp3(1)
      c22 = amp2(2)*amp3(2)
!
      nwam(1) = c11 + c22
      nwam(2) = c12 + c21
      nwam(3) = (0.D0,1.D0)*(c21 - c12)
      nwam(4) = c11 - c22
      nwam(5) = (0.D0,0.D0)
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftowr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double complex c11, c12, c21, c22
!
      c11 = amp2(1)*amp3(1)
      c12 = amp2(1)*amp3(2)
      c21 = amp2(2)*amp3(1)
      c22 = amp2(2)*amp3(2)
!
      nwam(1) = c11 + c22
      nwam(2) = -(c12 + c21)
      nwam(3) = (0.D0,-1.D0)*(c21 - c12)
      nwam(4) = c22 - c11
      nwam(5) = (0.D0,0.D0)
!
      end 
!
      subroutine glgr(fl, gl, gr)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      integer fl
      double precision gl, gr
!
      if (conv(fl) .eq. 1) then
         gl = glup
         gr = grup
      else if (conv(fl) .eq. 2) then
         gl = gldn
         gr = grdn
      else if (conv(fl) .eq. 3) then
         gl = gllp
         gr = grlp
      else
         write (*, *) 'something wrong in vma neutral interaction'
         stop 
      endif
!
      end 
!
c                          !fusion of fermion and fermion-bar into gluon
      subroutine fbftoz(flv, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      integer flv
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double complex c11, c12, c21, c22
      double precision gl, gr
      integer j1
!
      call glgr (flv, gl, gr)
!
      c11 = amp2(3)*amp3(1)
      c12 = amp2(3)*amp3(2)
      c21 = amp2(4)*amp3(1)
      c22 = amp2(4)*amp3(2)
!
      nwam(1) = c11 + c22
      nwam(2) = c12 + c21
      nwam(3) = (0.D0,1.D0)*(c21 - c12)
      nwam(4) = c11 - c22
      do j1 = 1, 6
         nwam(j1) = gl*nwam(j1)
      end do
!
      c11 = amp2(1)*amp3(3)
      c12 = amp2(1)*amp3(4)
      c21 = amp2(2)*amp3(3)
      c22 = amp2(2)*amp3(4)
!
      nwam(1) = nwam(1) + gr*(c11 + c22)
      nwam(2) = nwam(2) - gr*(c12 + c21)
      nwam(3) = nwam(3) - gr*(0.D0,1.D0)*(c21 - c12)
      nwam(4) = nwam(4) - gr*(c11 - c22)
      nwam(5) = 0.D0
!
      end 
!
c                        !fusion of fermionbar and gluon into fermionbar
      subroutine gfbtofb(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1)=(-amp3(3)*amp2(3))-amp3(4)*(amp2(1)+(0.D0,1.D0)*amp2(2))
      nwam(2)=(-amp3(3)*(amp2(1)-(0.D0,1.D0)*amp2(2)))+amp3(4)*amp2(3)
      nwam(3) = amp3(1)*amp2(3) + amp3(2)*(amp2(1)+(0.D0,1.D0)*amp2(2))
      nwam(4) = amp3(1)*(amp2(1)-(0.D0,1.D0)*amp2(2)) - amp3(2)*amp2(3)
!
      end 
!
      subroutine wfbtofb(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = (amp2(1)+amp2(4))*amp3(3) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(4)
      nwam(2) = (amp2(1)-amp2(4))*amp3(4) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(3)
      nwam(3) = 0.D0
      nwam(4) = 0.D0
!
      end 
!
c                        !fusion of fermionbar and gluon into fermionbar
      subroutine gfbtofbl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1)=(-amp3(1)*amp2(3))-amp3(2)*(amp2(1)+(0.D0,1.D0)*amp2(2))
      nwam(2)=(-amp3(1)*(amp2(1)-(0.D0,1.D0)*amp2(2)))+amp3(2)*amp2(3)
!
      end 
!
c                        !fusion of fermionbar and gluon into fermionbar
      subroutine gfbtofbr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = amp3(1)*amp2(3) + amp3(2)*(amp2(1)+(0.D0,1.D0)*amp2(2))
      nwam(2) = amp3(1)*(amp2(1)-(0.D0,1.D0)*amp2(2)) - amp3(2)*amp2(3)
!
      end 
!
      subroutine wfbtofbl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = (amp2(1)+amp2(4))*amp3(1) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(2) = (amp2(1)-amp2(4))*amp3(2) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
!
      end 
!
      subroutine wfbtofbr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = (amp2(1)-amp2(4))*amp3(1) - (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(2) = (amp2(1)+amp2(4))*amp3(2) - (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
!
      end 
!
      subroutine zfbtofb(flv, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      integer flv
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double precision gl, gr
      integer j1
!
      call glgr (flv, gl, gr)
!
      nwam(3) = (amp2(1)-amp2(4))*amp3(1) - (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(4) = (amp2(1)+amp2(4))*amp3(2) - (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
      do j1 = 1, 2
         nwam(j1+2) = gr*nwam(j1+2)
      end do
!
      nwam(1) = (amp2(1)+amp2(4))*amp3(3) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(4)
      nwam(2) = (amp2(1)-amp2(4))*amp3(4) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(3)
      do j1 = 1, 2
         nwam(j1) = gl*nwam(j1)
      end do
!
      end 
!
      subroutine gftof(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = amp3(3)*amp2(3) + amp3(4)*(amp2(1)-(0.D0,1.D0)*amp2(2))
      nwam(2) = amp3(3)*(amp2(1)+(0.D0,1.D0)*amp2(2)) - amp3(4)*amp2(3)
      nwam(3)=(-amp3(1)*amp2(3))-amp3(2)*(amp2(1)-(0.D0,1.D0)*amp2(2))
      nwam(4)=(-amp3(1)*(amp2(1)+(0.D0,1.D0)*amp2(2)))+amp3(2)*amp2(3)
!
      end 
!
      subroutine wftof(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = 0.D0
      nwam(2) = 0.D0
      nwam(3) = (amp2(1)+amp2(4))*amp3(1) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(4) = (amp2(1)-amp2(4))*amp3(2) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
!
      end 
!
      subroutine gftofl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1)=(-amp3(1)*amp2(3))-amp3(2)*(amp2(1)-(0.D0,1.D0)*amp2(2))
      nwam(2)=(-amp3(1)*(amp2(1)+(0.D0,1.D0)*amp2(2)))+amp3(2)*amp2(3)
!
      end 
!
      subroutine gftofr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = amp3(1)*amp2(3) + amp3(2)*(amp2(1)-(0.D0,1.D0)*amp2(2))
      nwam(2) = amp3(1)*(amp2(1)+(0.D0,1.D0)*amp2(2)) - amp3(2)*amp2(3)
!
      end 
!
      subroutine wftofl(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = (amp2(1)+amp2(4))*amp3(1) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(2) = (amp2(1)-amp2(4))*amp3(2) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
!
      end 
!
      subroutine wftofr(amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      double complex amp2(6), amp3(6)
      double complex nwam(6)
!
      nwam(1) = (amp2(1)-amp2(4))*amp3(1) - (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(2) = (amp2(1)+amp2(4))*amp3(2) - (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
!
      end 
!
      subroutine zftof(flv, amp2, amp3, nwam)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
      integer flv
      double complex amp2(6), amp3(6)
      double complex nwam(6)
      double precision gl, gr
      integer j1
!
      call glgr (flv, gl, gr)
!
      nwam(3) = (amp2(1)+amp2(4))*amp3(1) + (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(2)
      nwam(4) = (amp2(1)-amp2(4))*amp3(2) + (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(1)
      do j1 = 1, 2
         nwam(j1+2) = gl*nwam(j1+2)
      end do
!
      nwam(1) = (amp2(1)-amp2(4))*amp3(3) - (amp2(2)-(0.D0,1.D0)*amp2(3)
     1   )*amp3(4)
      nwam(2) = (amp2(1)+amp2(4))*amp3(4) - (amp2(2)+(0.D0,1.D0)*amp2(3)
     1   )*amp3(3)
      do j1 = 1, 2
         nwam(j1) = gr*nwam(j1)
      end do
!
      end 
!
      subroutine prpgt(fl, imp, fus, pol)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nmax
      parameter (nmax = 10)
      integer npmax
      parameter (npmax = 43)
      integer nmrgmx
      parameter (nmrgmx = 5001)
      integer nlor
      parameter (nlor = 6)
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
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
      integer fl                                 !flavour
      double precision imp(4)                    !momenta
c                       !amplitude prior to matching with the propagator
      double complex fus(6)
      double complex pol(6)                     !propagator times fusion
      integer j1
!
      if (lortyp(fl) .eq. 1) then          !massless gauge boson, F_munu
         call gbzms (imp, fus, pol)
      else if (lortyp(fl) .eq. 11) then    !massless gauge boson, F_munu
         call gbph (imp, fus, pol)
      else if (lortyp(fl) .eq. 2) then           !massless gauge boson
         call gbms (imp, fus, masses(fl), width(fl), pol)
      else if (lortyp(fl) .eq. 3) then           !fermion
         call psl (imp, fus, masses(fl), pol)
      else if (lortyp(fl) .eq. 4) then           !fermion-bar
         call psltr (imp, fus, masses(fl), pol)
      else if (lortyp(fl) .eq. 5) then           !fermion
         call psll (imp, fus, pol)
      else if (lortyp(fl) .eq. 6) then           !fermion-bar
         call psltrl (imp, fus, pol)
      else if (lortyp(fl) .eq. 7) then           !fermion
         call pslr (imp, fus, pol)
      else if (lortyp(fl) .eq. 8) then           !fermion-bar
         call psltrr (imp, fus, pol)
      else if (lortyp(fl) .eq. 9) then           !aux X++
         pol(1) = fus(1)
      else if (lortyp(fl) .eq. 10) then          !aux X++
         do j1 = 1, 2
            pol(j1) = fus(j1)
         end do
      else if (lortyp(fl) .eq. 12) then
         call scl (imp, fus, masses(fl), width(fl), pol)
      else if (lortyp(fl) .eq. 13) then
         do j1 = 1, 6
            pol(j1) = fus(j1)
         end do
      else
         write (*, *) 'Wrong choice of propagator'
         stop 
      endif
!
      end 
!
      subroutine gbzms(imp, fus, pol)   !massless gauge boson propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq, ip0
      double complex pa
!
      logical flg90
      common /f90/flg90! to avoid spurious divergences in temporal gauge
      doubleprecision d1, d2
      doublecomplex x1
      integer j1
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      ip0 = imp(1)*imp(1)
      if (ip0 .eq. 0.D0) then
         write (*, *) 'singular coulomb propagator'
         stop 
      endif
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + imp(j1+1)*fus(j1)
      end do
      pa = x1
      pa = pa/ip0
      do j1 = 1, 3
         pol(j1) = (fus(j1)-imp(j1+1)*pa)*pq
         pol(j1+3) = -fus(j1+3)
      end do
!
      if (.not.flg90) then
         d2 = 0
         do j1 = 1, 3
            d2 = d2 + abs(imp(j1+1))
         end do
         ip0 = d2
         if (ip0.ne.0.D0 .and. abs(imp(1)/ip0).lt.1.D-4) flg90 = .TRUE.
      endif
!
      end 
!
      subroutine gbph(imp, fus, pol)    !massless gauge boson propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq, ip0
      double complex pa
!
      logical flg90
      common /f90/flg90! to avoid spurious divergences in temporal gauge
      doubleprecision d1, d2
      doublecomplex x1
      integer j1
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      ip0 = imp(1)*imp(1)
      if (ip0 .eq. 0.D0) then
         write (*, *) 'singular coulomb propagator'
         stop 
      endif
      x1 = 0
      do j1 = 1, 3
         x1 = x1 + imp(j1+1)*fus(j1)
      end do
      pa = x1
      pa = pa/ip0
      do j1 = 1, 3
         pol(j1) = (fus(j1)-imp(j1+1)*pa)*pq
      end do
      pol(4) = fus(4)
!
      if (.not.flg90) then
         d2 = 0
         do j1 = 1, 3
            d2 = d2 + abs(imp(j1+1))
         end do
         ip0 = d2
         if (ip0.ne.0.D0 .and. abs(imp(1)/ip0).lt.1.D-4) flg90 = .TRUE.
      endif
!
      end 
!
c                                       !massless gauge boson propagator
      subroutine gbms(imp, fus, mass, width, pol)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4), mass, width
      double complex fus(6)
      double complex pol(6)
!
      double precision mq, theta
      double complex pa, pq
!
      character wmode*2
      character resonance
      double precision winsize
      common /gauinv/winsize, resonance, wmode
      doubleprecision d1
      doublecomplex x1
      integer j1
!
      if (resonance .eq. 'y') then
         do j1 = 1, 6
            pol(j1) = (0.D0,0.D0)
         end do
         return 
      endif
      mq = mass*mass
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1 - mq             ! +(0.d0,1.d0)*mass*width
      if (wmode .eq. 'yy') then
         pq = 1.D0/(pq + (0.D0,1.D0)*mass*width)
         resonance = 'n'
      else if (wmode .eq. 'yn') then
         if (abs(pq) .lt. winsize*mass*width) then
            resonance = 'y'
         else
            resonance = 'n'
         endif
         if (pq .eq. 0.D0) pq = 1.D0
         pq = 1/pq
      else if (wmode .eq. 'nn') then
         pq = 1/pq
      endif
      x1 = 0
      do j1 = 1, 4
         x1 = x1 + imp(j1)*fus(j1)
      end do
      pa = x1
      pa = pa/mq
      do j1 = 1, 4
         pol(j1) = (-fus(1)) + imp(1)*pa
      end do
      do j1 = 1, 3
         pol(j1+1) = fus(j1+1) + imp(j1+1)*pa
      end do
      do j1 = 1, 4
         pol(j1) = pol(j1)*pq
      end do
      pol(5) = fus(5)
!  pol(1)=-pol(1)
!
      end 
!
c                                       !massless gauge boson propagator
      subroutine scl(imp, fus, mass, width, pol)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4), mass, width
      double complex fus(3)
      double complex pol(3)
!
      double precision mq, theta
      double complex pa, pq
!
      character wmode*2
      character resonance
      double precision winsize
      common /gauinv/winsize, resonance, wmode
      doubleprecision d1
      integer j1
!
      if (resonance .eq. 'y') then
         do j1 = 1, 3
            pol(j1) = (0.D0,0.D0)
         end do
         return 
      endif
      mq = mass*mass
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1 - mq
      if (wmode .eq. 'yy') then
         pq = 1.D0/(pq + (0.D0,1.D0)*mass*width)
         resonance = 'n'
      else if (wmode .eq. 'yn') then
         if (abs(pq) .lt. winsize*mass*width) then
            resonance = 'y'
         else
            resonance = 'n'
         endif
         if (pq .eq. 0.D0) pq = 1.D0
         pq = 1/pq
      else if (wmode .eq. 'nn') then
         pq = 1/pq
      endif
      pol(1) = pq*fus(1)
      pol(2) = fus(3)
      pol(3) = fus(2)
!
      end 
!
      subroutine psl(imp, fus, mass, pol)        !fermion propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double precision mass
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1) = fus(1)*mass + fus(3)*(imp(1)-imp(4)) - fus(4)*(imp(2)-
     1   (0.D0,1.D0)*imp(3))
      pol(2) = fus(2)*mass - fus(3)*(imp(2)+(0.D0,1.D0)*imp(3)) + fus(4)
     1   *(imp(1)+imp(4))
      pol(3) = fus(1)*(imp(4)+imp(1)) + fus(2)*(imp(2)-(0.D0,1.D0)*imp(3
     1   )) + fus(3)*mass
      pol(4) = fus(1)*(imp(2)+(0.D0,1.D0)*imp(3)) + fus(2)*(imp(1)-imp(4
     1   )) + fus(4)*mass
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = pq - mass**2
      pq = 1/pq
      do j1 = 1, 4
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine pslr(imp, fus, pol)             !fermion propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1)=fus(1)*(imp(4)+imp(1))+fus(2)*(imp(2)-(0.D0,1.D0)*imp(3))
      pol(2)=fus(1)*(imp(2)+(0.D0,1.D0)*imp(3))+fus(2)*(imp(1)-imp(4))
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      do j1 = 1, 2
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine psll(imp, fus, pol)             !fermion propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1)=fus(1)*(imp(1)-imp(4))-fus(2)*(imp(2)-(0.D0,1.D0)*imp(3))
      pol(2) = (-fus(1)*(imp(2)+(0.D0,1.D0)*imp(3))) + fus(2)*(imp(1)+
     1   imp(4))
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      do j1 = 1, 2
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine psltr(imp, fus, mass, pol)      !fermion bar propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double precision mass
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1) = fus(1)*mass - fus(3)*(imp(4)+imp(1)) - fus(4)*(imp(2)+
     1   (0.D0,1.D0)*imp(3))
      pol(2) = fus(2)*mass - fus(3)*(imp(2)-(0.D0,1.D0)*imp(3)) + fus(4)
     1   *(imp(4)-imp(1))
      pol(3) = fus(1)*(imp(4)-imp(1)) + fus(2)*(imp(2)+(0.D0,1.D0)*imp(3
     1   )) + fus(3)*mass
      pol(4) = fus(1)*(imp(2)-(0.D0,1.D0)*imp(3)) - fus(2)*(imp(1)+imp(4
     1   )) + fus(4)*mass
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = pq - mass**2
      pq = 1/pq
      do j1 = 1, 4
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine psltrr(imp, fus, pol)           !fermion bar propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1) = (-fus(1)*(imp(4)+imp(1))) - fus(2)*(imp(2)+(0.D0,1.D0)*
     1   imp(3))
      pol(2) = (-fus(1)*(imp(2)-(0.D0,1.D0)*imp(3))) + fus(2)*(imp(4)-
     1   imp(1))
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      do j1 = 1, 2
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine psltrl(imp, fus, pol)           !fermion bar propagator
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nlor
      parameter (nlor = 6)
      double precision imp(4)
      double complex fus(6)
      double complex pol(6)
!
      double precision pq
      doubleprecision d1
      integer j1
!
      pol(1)=fus(1)*(imp(4)-imp(1))+fus(2)*(imp(2)+(0.D0,1.D0)*imp(3))
      pol(2)=fus(1)*(imp(2)-(0.D0,1.D0)*imp(3))-fus(2)*(imp(1)+imp(4))
      d1 = 0
      do j1 = 1, 3
         d1 = d1 + imp(j1+1)*imp(j1+1)
      end do
      pq = imp(1)*imp(1) - d1
      pq = 1/pq
      do j1 = 1, 2
         pol(j1) = pol(j1)*pq
      end do
!
      end 
!
      subroutine colprdsu3
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincol.inc'
      integer j1, j2, j3, i2, k1, k2, i1
c               !tabmom(i,j,k) = coulor(j) * coulor(k) (i=1,2 because of
      integer tabmom(2,-11:11,-11:11)
!                                                      neutral gluons
c                        ! number of coulors of the new object, == 0,1,2
      integer tabnew(-11:11,-11:11)
      double precision tabcoeff(2,-11:11,-11:11) !coulor clebsch
      integer j4, j5, j6
      do j5 = 1, 23
         do j4 = 1, 23
            tabnew(j4-12,j5-12) = 0
         end do
      end do
      do j6 = 1, 23
         do j5 = 1, 23
            do j4 = 1, 2
               tabmom(j4,j5-12,j6-12) = 0
               tabcoeff(j4,j5-12,j6-12) = 0
            end do
         end do
      end do
!
      tabnew(1,1) = 0
!
      tabnew(1,3) = 1
      tabcoeff(1,1,3) = 1./sqrt(2.D0)
      tabmom(1,1,3) = 3
      tabnew(1,4) = 1
      tabcoeff(1,1,4) = sqrt(2.D0)
      tabmom(1,1,4) = 4
      tabnew(1,6) = 1
      tabcoeff(1,1,6) = -1./sqrt(2.D0)
      tabmom(1,1,6) = 6
      tabnew(1,2) = 0
      tabnew(1,5) = 1
      tabcoeff(1,1,5) = 1./sqrt(2.D0)
      tabmom(1,1,5) = 5
      tabnew(1,7) = 1
      tabcoeff(1,1,7) = -sqrt(2.D0)
      tabmom(1,1,7) = 7
      tabnew(1,8) = 1
      tabcoeff(1,1,8) = -1./sqrt(2.D0)
      tabmom(1,1,8) = 8
      tabnew(3,1) = 1
      tabcoeff(1,3,1) = -1./sqrt(2.D0)
      tabmom(1,3,1) = 3
      tabnew(3,3) = 0
      tabnew(3,4) = 0
      tabnew(3,6) = 2
      tabcoeff(1,3,6) = 1./sqrt(2.D0)
      tabcoeff(2,3,6) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,3,6) = 1
      tabmom(2,3,6) = 2
      tabnew(3,2) = 1
      tabcoeff(1,3,2) = sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,3,2) = 3
      tabnew(3,5) = 1
      tabcoeff(1,3,5) = 1.
      tabmom(1,3,5) = 4
      tabnew(3,7) = 1
      tabcoeff(1,3,7) = -1.
      tabmom(1,3,7) = 8
      tabnew(3,8) = 0
      tabnew(4,1) = 1
      tabcoeff(1,4,1) = -sqrt(2.D0)
      tabmom(1,4,1) = 4
      tabnew(4,3) = 0
      tabnew(4,4) = 0
      tabnew(4,6) = 1
      tabcoeff(1,4,6) = -1.
      tabmom(1,4,6) = 5
      tabnew(4,2) = 0
      tabnew(4,5) = 0
      tabnew(4,7) = 2
      tabcoeff(1,4,7) = sqrt(2.D0)
      tabcoeff(2,4,7) = 0.D0
      tabmom(1,4,7) = 1
      tabmom(2,4,7) = 2
      tabnew(4,8) = 1
      tabcoeff(1,4,8) = 1.
      tabmom(1,4,8) = 3
      tabnew(6,1) = 1
      tabcoeff(1,6,1) = 1./sqrt(2.D0)
      tabmom(1,6,1) = 6
      tabnew(6,3) = 2
      tabcoeff(1,6,3) = -1./sqrt(2.D0)
      tabcoeff(2,6,3) = sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,6,3) = 1
      tabmom(2,6,3) = 2
      tabnew(6,4) = 1
      tabcoeff(1,6,4) = 1.
      tabmom(1,6,4) = 5
      tabnew(6,6) = 0
      tabnew(6,2) = 1
      tabcoeff(1,6,2) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,6,2) = 6
      tabnew(6,5) = 0
      tabnew(6,7) = 0
      tabnew(6,8) = 1
      tabcoeff(1,6,8) = -1.
      tabmom(1,6,8) = 7
      tabnew(2,1) = 0
      tabnew(2,3) = 1
      tabcoeff(1,2,3) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,2,3) = 3
      tabnew(2,4) = 0
      tabnew(2,6) = 1
      tabcoeff(1,2,6) = sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,2,6) = 6
      tabnew(2,2) = 0
      tabnew(2,5) = 1
      tabcoeff(1,2,5) = sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,2,5) = 5
      tabnew(2,7) = 0
      tabnew(2,8) = 1
      tabcoeff(1,2,8) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,2,8) = 8
      tabnew(5,1) = 1
      tabcoeff(1,5,1) = -1./sqrt(2.D0)
      tabmom(1,5,1) = 5
      tabnew(5,3) = 1
      tabcoeff(1,5,3) = -1
      tabmom(1,5,3) = 4
      tabnew(5,4) = 0
      tabnew(5,6) = 0
      tabnew(5,2) = 1
      tabcoeff(1,5,2) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,5,2) = 5
      tabnew(5,5) = 0
      tabnew(5,7) = 1
      tabcoeff(1,5,7) = 1.
      tabmom(1,5,7) = 6
      tabnew(5,8) = 2
      tabcoeff(1,5,8) = 1./sqrt(2.D0)
      tabcoeff(2,5,8) = sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,5,8) = 1
      tabmom(2,5,8) = 2
      tabnew(7,1) = 1
      tabcoeff(1,7,1) = sqrt(2.D0)
      tabmom(1,7,1) = 7
      tabnew(7,3) = 1
      tabcoeff(1,7,3) = 1.
      tabmom(1,7,3) = 8
      tabnew(7,4) = 2
      tabcoeff(1,7,4) = -sqrt(2.D0)
      tabcoeff(2,7,4) = 0.D0
      tabmom(1,7,4) = 1
      tabmom(2,7,4) = 2
      tabnew(7,6) = 0
      tabnew(7,2) = 0
      tabnew(7,5) = 1
      tabcoeff(1,7,5) = -1.
      tabmom(1,7,5) = 6
      tabnew(7,7) = 0
      tabnew(7,8) = 0
      tabnew(8,1) = 1
      tabcoeff(1,8,1) = 1./sqrt(2.D0)
      tabmom(1,8,1) = 8
      tabnew(8,3) = 0
      tabnew(8,4) = 1
      tabcoeff(1,8,4) = -1.
      tabmom(1,8,4) = 3
      tabnew(8,6) = 1
      tabcoeff(1,8,6) = 1.
      tabmom(1,8,6) = 7
      tabnew(8,2) = 1
      tabcoeff(1,8,2) = sqrt(3.D0/2.D0)
      tabmom(1,8,2) = 8
      tabnew(8,5) = 2
      tabcoeff(1,8,5) = -1./sqrt(2.D0)
      tabcoeff(2,8,5) = -sqrt(3.D0)/sqrt(2.D0)
      tabmom(1,8,5) = 1
      tabmom(2,8,5) = 2
      tabnew(8,7) = 0
      tabnew(8,8) = 0
      tabnew(9,-9) = 2
      tabcoeff(1,9,-9) = 1./sqrt(2.D0)
      tabcoeff(2,9,-9) = -1./sqrt(6.D0)
      tabmom(1,9,-9) = 1
      tabmom(2,9,-9) = 2
      tabnew(9,-10) = 1
      tabcoeff(1,9,-10) = 1.
      tabmom(1,9,-10) = 6
      tabnew(9,-11) = 1
      tabcoeff(1,9,-11) = 1.
      tabmom(1,9,-11) = 7
      tabnew(10,-9) = 1
      tabcoeff(1,10,-9) = 1.
      tabmom(1,10,-9) = 3
      tabnew(10,-10) = 2
      tabcoeff(1,10,-10) = 0.D0
      tabcoeff(2,10,-10) = sqrt(2.D0/3.D0)
      tabmom(1,10,-10) = 1
      tabmom(2,10,-10) = 2
      tabnew(10,-11) = 1
      tabcoeff(1,10,-11) = 1.
      tabmom(1,10,-11) = 8
      tabnew(11,-9) = 1
      tabcoeff(1,11,-9) = 1.
      tabmom(1,11,-9) = 4
      tabnew(11,-10) = 1
      tabcoeff(1,11,-10) = 1.
      tabmom(1,11,-10) = 5
      tabnew(11,-11) = 2
      tabcoeff(1,11,-11) = -1./sqrt(2.D0)
      tabcoeff(2,11,-11) = -1./sqrt(6.D0)
      tabmom(1,11,-11) = 1
      tabmom(2,11,-11) = 2
      tabnew(1,9) = 1
      tabcoeff(1,1,9) = 1./sqrt(2.D0)
      tabmom(1,1,9) = 9
      tabnew(1,11) = 1
      tabcoeff(1,1,11) = -1./sqrt(2.D0)
      tabmom(1,1,11) = 11
      tabnew(3,9) = 1
      tabcoeff(1,3,9) = 1.
      tabmom(1,3,9) = 10
      tabnew(4,9) = 1
      tabcoeff(1,4,9) = 1.
      tabmom(1,4,9) = 11
      tabnew(6,10) = 1
      tabcoeff(1,6,10) = 1.
      tabmom(1,6,10) = 9
      tabnew(2,9) = 1
      tabcoeff(1,2,9) = -1./sqrt(6.D0)
      tabmom(1,2,9) = 9
      tabnew(2,10) = 1
      tabcoeff(1,2,10) = sqrt(2.D0/3.D0)
      tabmom(1,2,10) = 10
      tabnew(2,11) = 1
      tabcoeff(1,2,11) = -1./sqrt(6.D0)
      tabmom(1,2,11) = 11
      tabnew(5,10) = 1
      tabcoeff(1,5,10) = 1.
      tabmom(1,5,10) = 11
      tabnew(7,11) = 1
      tabcoeff(1,7,11) = 1.
      tabmom(1,7,11) = 9
      tabnew(8,11) = 1
      tabcoeff(1,8,11) = 1.
      tabmom(1,8,11) = 10
      tabnew(1,-9) = 1
      tabcoeff(1,1,-9) = 1./sqrt(2.D0)
      tabmom(1,1,-9) = -9
      tabnew(1,-11) = 1
      tabcoeff(1,1,-11) = -1./sqrt(2.D0)
      tabmom(1,1,-11) = -11
      tabnew(3,-10) = 1
      tabcoeff(1,3,-10) = 1.
      tabmom(1,3,-10) = -9
      tabnew(4,-11) = 1
      tabcoeff(1,4,-11) = 1.
      tabmom(1,4,-11) = -9
      tabnew(6,-9) = 1
      tabcoeff(1,6,-9) = 1.
      tabmom(1,6,-9) = -10
      tabnew(2,-9) = 1
      tabcoeff(1,2,-9) = -1./sqrt(6.D0)
      tabmom(1,2,-9) = -9
      tabnew(2,-10) = 1
      tabcoeff(1,2,-10) = sqrt(2.D0/3.D0)
      tabmom(1,2,-10) = -10
      tabnew(2,-11) = 1
      tabcoeff(1,2,-11) = -1./sqrt(6.D0)
      tabmom(1,2,-11) = -11
      tabnew(5,-11) = 1
      tabcoeff(1,5,-11) = 1.
      tabmom(1,5,-11) = -10
      tabnew(7,-9) = 1
      tabcoeff(1,7,-9) = 1.
      tabmom(1,7,-9) = -11
      tabnew(8,-10) = 1
      tabcoeff(1,8,-10) = 1.
      tabmom(1,8,-10) = -11
      do j1 = 9, 11
         tabnew(0,j1) = 1
         tabcoeff(1,0,j1) = 1.
         tabmom(1,0,j1) = j1
         tabnew(j1,0) = 1
         tabcoeff(1,j1,0) = 1.
         tabmom(1,j1,0) = j1
      end do
      do j1 = -11, -9
         tabnew(0,j1) = 1
         tabcoeff(1,0,j1) = 1.
         tabmom(1,0,j1) = j1
         tabnew(j1,0) = 1
         tabcoeff(1,j1,0) = 1.
         tabmom(1,j1,0) = j1
      end do
      k1 = 0
      k2 = 0
      do j1 = 1, 23
         do j2 = 1, 23
            i1 = j1 - 12
            i2 = j2 - 12
            k1 = k1 + 1
            tbnew(k1) = tabnew(i2,i1)
            do j3 = 1, 2
               k2 = k2 + 1
               tbmom(k2) = tabmom(j3,i2,i1)
               tbcoeff(k2) = tabcoeff(j3,i2,i1)
            end do
         end do
      end do
      do j4 = 1, 8
         tbnew(j4+265) = 1
         tbcoeff(j4*2+529) = 1.D0
      end do
      tbmom(1*2+529) = 1
      tbmom(2*2+529) = 2
      tbmom(3*2+529) = 3
      tbmom(4*2+529) = 4
      tbmom(5*2+529) = 5
      tbmom(6*2+529) = 6
      tbmom(7*2+529) = 7
      tbmom(8*2+529) = 8
      do j4 = 1, 8
         tbnew(j4*23+265) = 1
         tbcoeff(j4*46+529) = 1.D0
      end do
      tbmom(1*46+529) = 1
      tbmom(2*46+529) = 2
      tbmom(3*46+529) = 3
      tbmom(4*46+529) = 4
      tbmom(5*46+529) = 5
      tbmom(6*46+529) = 6
      tbmom(7*46+529) = 7
      tbmom(8*46+529) = 8
!
      end 
!
      subroutine permtn(nfr, p1, p2, perm)
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer nfr, p1, p2
      double precision perm
!
      integer odd, j3, j1
!
      j3 = 0
      odd = 0
      do j1 = 0, nfr - 1
         if (btest(p2,j1)) then
            odd = odd + 1
         else
            if (btest(odd,0)) then
               if (btest(p1,j1)) j3 = j3 + 1
            endif
         endif
      end do
!
      if (btest(j3,0)) then
         perm = -1.D0
      else
         perm = 1.D0
      endif
!
      end 
!
      subroutine opt(j)
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      integer j
!
      j = j + 1
!
      end 
!
      function popcnt0 (i)
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      integer i, popcnt0
      integer j1
!
      popcnt0 = 0
      do j1 = 0, 31
         if (btest(i,j1)) popcnt0 = popcnt0 + 1
      end do
!
      return 
      end 
!
      subroutine ampd(prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1   prcss__mom, dim0, mtel)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      common /strfld/momd, lortyp, nlordof, colrep, fldptr, nexc, nexcsm
      integer momd(2047)
      integer lortyp(43)
      integer nlordof(43)
      integer colrep(43)
      integer fldptr(2,43,5)
      integer nexcmax
      parameter (nexcmax = 13)
      integer nexc(13,2:5)
      integer nexcsm(13,2:5)
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
      parameter (nintmax = 50)
      integer tabint(151,1)
      integer tabintmrg(451,1)
      integer opam(50,1)
      integer opfs(150,1)
      integer nvasfs(151,1)
      integer nvasam(51,1)
      double precision cpam(50,1)
      double precision cpfs(150,1)
      integer nprc
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
      INCLUDE 'Vincol.inc'
      INCLUDE 'Vdimen.inc'
      INCLUDE 'Vutilities.inc'
      INCLUDE 'Vrunning.inc'
      integer prcss__flv(10)
      integer prcss__hel(10)
      integer prcss__col(20)
      integer prcss__nfrm
      double precision prcss__mom(40)
      integer dim0(7)
      double complex mtel                        !amplitude
!
      integer nca, dim(7)
      parameter (nca = 4)
      integer nprt, nfld, flvold, padd, nconf, j5, padd0
      integer j1,j2,iter,steptm,stepex,ex2,ex3,nf2,nf3,psfl2,psfl3
      integer pscfl2,pscfl3,nl1,nl2,nl3,j3,j4,beg3,nc,psfl3i,pscfl3i
      integer p1, p2, p3, posfld, poscol, nlold
      integer nf1, psfl1, pscfl1, pref, beg
      integer cl(4), clnw(2)
      integer psnew(1050)
      integer fststp
      save fststp
      integer tabmrgpr(451)
      integer exc(13), excsm(13), excns(13)
      integer offshcol(5500), ampcol(5780), ordtmp(1050)
      double complex offsh(5370), ampfld(5550)
      double complex fusion(6), polar(6)
      double complex amplitude, tmp
      double precision colcff
      double precision impul(4,1050)
      double precision imp(4)
      double precision perm
!
      integer nphf, even, nmom, ot, npmit, tmpi, ci, cf, ptr, nfr2
      integer operot, jop, op(150), smold
      logical choose, choose2, choose3, choose4, cnsold(0:43)
      double precision imp3(4), imp2(4)
      double precision coupot, cpot1, cp(150)
      integer psfit, pscit, psfa, psca
      save psfit, pscit, psfa, psca
      logical dbg
      parameter (dbg = .FALSE.)
      integer j9, j10
      doubleprecision d1(4), d2(4)
      doublecomplex x1
      integer j6, j7, j8
!
! once for all inizialization
      dim(1) = 5370
      dim(2) = 5500
      dim(3) = 5550
      dim(4) = 5780
      dim(5) = 1050
      dim(6) = 1050
      dim(7) = 1050
      data psca/ 0/ 
      data psfa/ 0/ 
      data pscit/ 0/ 
      data psfit/ 0/ 
      data fststp/ 0/ 
      if (fststp .eq. 0) then
         fststp = 1
! inizialization of arrays
         call initcoup
         call initmom
         call inittab
         call colprdsu3
!
! inizialization of merging tables
         call initprc
!
      endif
!
      do j8 = 1, 5
         do j7 = 1, 43
            do j6 = 1, 2
               fldptr(j6,j7,j8) = -9999
            end do
         end do
      end do
      nprt = 0                             ! number of external particle
      do j1 = 1, 10
         if (prcss__flv(j1) .eq. 1001) go to 2
         nprt = nprt + 1
      end do
    2 continue
!
      nphf = nprt/2                              ! for optmization
      even = mod(nprt,2)                         ! for optmization
c                     ! maximum number of different combinations of nprt
      nmom = 2**nprt - 1
!                            momenta, each momenta enters at most once
! inizialization step
!
      call rncpl               !loading running couplings (alpha-strong)
!
c             ! each off shell amplitude will be collected into a single
      posfld = 1
c           ! array OFFSH, POSFLD and POSCOL allows to find the position
      poscol = 1
      flvold = -99
      do j1 = 1, nprt
         if (prcss__flv(j1) .ne. flvold) then
            if (flvold .gt. 0) fldptr(1,flvold,1) = nfld
            nfld = 1
            flvold = prcss__flv(j1)
            fldptr(2,flvold,1) = poscol
         else
            nfld = nfld + 1
         endif
         p1 = 2**(j1 - 1)
c                                       !OFFSH(posfld,....) contains the
         call initpol (prcss__flv, prcss__hel, prcss__col, prcss__nfrm, 
     1      prcss__mom, j1, offsh(posfld))
!                                        polarization of the jth particle
         offshcol(1-1+poscol) = p1
         offshcol(2-1+poscol) = posfld
         offshcol(3-1+poscol) = prcss__col(2*j1-1)
         offshcol(4-1+poscol) = prcss__col(2*j1)
         do j6 = 1, 4
            impul(j6,offshcol(poscol)) = prcss__mom(j6+4*(j1-1))
         end do
         posfld = posfld + nlordof(prcss__flv(j1))
         poscol = poscol + 4
      end do
      fldptr(1,flvold,1) = nfld
!
! iteration
! we merge on/off shell field2 and field3 into field1 via equation of motion up to the step required by the algorithm
      do j6 = 1, 451
         tabmrgpr(j6) = tabintmrg(j6,nprc)
      end do
      do j6 = 1, tabintmrg(1,nprc)
         cp(j6) = cpfs(j6,nprc)
         op(j6) = opfs(j6,nprc)
      end do
      do j1 = 1, nvasfs(1,nprc)
         cp(nvasfs(j1+1,nprc)) = cp(nvasfs(j1+1,nprc))*rnas
      end do
      do j6 = 1, nmom
         psnew(j6) = 0
      end do
      do iter = 2, nphf
         do j6 = 1, 44
            cnsold(j6-1) = .FALSE.
         end do
         steptm = 2                              !for optimization
         flvold = -99
         do j6 = 1, 13
            excsm(j6) = nexcsm(j6,iter)
            excns(j6) = nexc(j6,iter)
         end do
         do j1 = 1, tabmrgpr(1)
c                 !labels rthe flavour of the daughter offshell particle
            fl1 = tabmrgpr(steptm)
c                 !fl2 and fl3 label the flavour of the parents offshell
            fl2 = tabmrgpr(steptm+1)
            fl3 = tabmrgpr(steptm+2)             !particles to be merged
!
            coupot = cp(j1)                      !coup(fl1,fl2,fl3)
            operot = op(j1)                      !oper(fl1,fl2,fl3)
!
c                  !nl1,nl2,nl3 are the corresponding numbers of lorentz
            nl1 = nlordof(fl1)
            nl2 = nlordof(fl2)                   !degrees of freedom
            nl3 = nlordof(fl3)
            if (fl1 .ne. flvold) then
               if (flvold .gt. 0) then
                  fldptr(1,chcg(flvold),iter) = nconf
               else
                  smold = 0
                  if (fl2 .eq. fl3) smold = 1
               endif
               flvold = fl1
               if (cnsold(cnschg(fl1))) then
                  do j6 = 1, nmom
                     psnew(j6) = 0               !for optimization
                  end do
                  do j6 = 1, 44
                     cnsold(j6-1) = .FALSE.
                  end do
               endif
               cnsold(cnschg(fl1)) = .TRUE.
               nconf = 0
               fldptr(2,chcg(fl1),iter) = poscol
            endif
            steptm = steptm + 3
            if (fl2 .ne. fl3) then
               if (smold .eq. 0) then
                  do j6 = 1, 13
                     exc(j6) = excns(j6)
                  end do
               endif
               smold = 1
            else
               if (smold .eq. 1) then
                  do j6 = 1, 13
                     exc(j6) = excsm(j6)
                  end do
               endif
               smold = 0
            endif
            stepex = 2
            do j2 = 1, exc(1)
               ex2 = exc(stepex)
               ex3 = exc(stepex+1)
               stepex = stepex + 2
               nf2 = fldptr(1,fl2,ex2)   !number of excitation of field2
               nf3 = fldptr(1,fl3,ex3)   !number of excitation of field3
!
               if (nf2.gt.0 .and. nf3.gt.0) then
!
c                         !starting position of field2 in array OFFSHCOL
                  pscfl2 = fldptr(2,fl2,ex2)
c                         !starting position of field2 in array OFFSHCOL
                  pscfl3 = fldptr(2,fl3,ex3)
c                            !starting position of field2 in array OFFSH
                  psfl2 = offshcol(pscfl2+1)
c                            !starting position of field3 in array OFFSH
                  psfl3 = offshcol(pscfl3+1)
                  do j3 = 1, nf2
c                    !if the flavour of field 2 and 3 are the same don't
                     if (ex2.eq.ex3 .and. fl2.eq.fl3) then
                        beg3 = j3 + 1 !repeat twice the same computation
                        psfl3i = psfl3 + nl3*(beg3 - 2)
                        pscfl3i = pscfl3 + 4*(beg3 - 2)
                     else
                        beg3 = 1
                        psfl3i = psfl3 - nl3
                        pscfl3i = pscfl3 - 4
                     endif
!
c                                !momentum of field2 off-shell amplitude
                     p2 = offshcol(pscfl2)
!          cl(1:2)=(/offshcol(pscfl2+2:pscfl2+3)/)
                     cl(1) = offshcol(pscfl2+2)
                     cl(2) = offshcol(pscfl2+3)
                     do j6 = 1, 4
                        imp2(j6) = impul(j6,p2)
                     end do
!          imp2(1)=impul(1,p2); imp2(2)=impul(2,p2); imp2(3)=impul(3,p2); imp2(4)=impul(4,p2);
!
c       !for optimization: true for even number of external particle and
                     choose = iter.eq.nphf .and. even.eq.0
!                                                    !the resulting off shell field is the merging of half this number
!
                     if (prcss__nfrm .gt. 2) then
                        nfr2 = 0
                        do j4 = 0, prcss__nfrm - 1
                           if (btest(p2,j4)) nfr2 = nfr2 + 1
                        end do
                     else
                        nfr2 = 0
                     endif
                     do j4 = beg3, nf3
                        pscfl3i = pscfl3i + 4
c                            !starting position of field3 in array OFFSH
                        psfl3i = offshcol(pscfl3i+1)
 
!
                        p3 = offshcol(pscfl3i)   !momentum of field3
c           !momentum of field1: resulting out of the merger of field2,3
                        p1 = p2 + p3
                        if (momd(p1) .eq. iter) then
c               !for even number of particle we accept merger of npart/2
                           if (.not.(choose .and. btest(p1+1,0))) then
!                                                        !momenta only if momenta 1 contribute
!            cl(3:4)=(/offshcol(pscfl3i+2:pscfl3i+3)/)  !contains the coulors of field2 and 3
                              cl(3) = offshcol(pscfl3i+2)
c                                  !contains the coulors of field2 and 3
                              cl(4) = offshcol(pscfl3i+3)
c                                                !quarks and gluons
                              if (colrep(fl1) .eq. 3) then
                                 if (cl(2) .eq. 0) go to 4
                                 if (cl(2) .eq. cl(3)) then
                                    if (cl(1) .ne. 0) then
                                    if (cl(1) .eq. cl(4)) go to 4
                                    colcff = 1.D0
                                    else
                                    colcff = -0.333333333333333D0
                                    endif
                                    clnw(1) = cl(1)
                                    clnw(2) = cl(4)
                                 else if (cl(1) .eq. cl(4)) then
                                    if (cl(1) .ne. 0) then
                                    colcff = -1.D0
                                    else
                                    colcff = 1.D0
                                    endif
                                    clnw(1) = cl(3)
                                    clnw(2) = cl(2)
                                 else if (colrep(fl2).eq.1 .and. colrep(
     1                                 fl3).eq.3) then
                                    colcff = 1.D0
                                    do j6 = 1, 2
                                    clnw(j6) = cl(j6+2)
                                    end do
                                 else if (colrep(fl2).eq.3 .and. colrep(
     1                                 fl3).eq.1) then
                                    colcff = 1.D0
                                    do j6 = 1, 2
                                    clnw(j6) = cl(j6)
                                    end do
                                 else
                                    go to 4
                                 endif
c                                                !coulorless object
                              else if (colrep(fl1) .eq. 2) then
                                 if (cl(2) .eq. cl(3)) then
                                    clnw(1) = cl(1)
                                    clnw(2) = cl(4)
                                 else if (cl(1) .eq. cl(4)) then
                                    clnw(1) = cl(3)
                                    clnw(2) = cl(2)
                                 else
                                    go to 4
                                 endif
                                 colcff = 1.D0
c                                                !coulorless object
                              else if (colrep(fl1) .eq. 1) then
                                 if (cl(2).eq.cl(3) .and. cl(1).eq.cl(4)
     1                              ) then
                                    clnw(1) = 0
                                    clnw(2) = 0
                                 else
                                    go to 4
                                 endif
                                 colcff = 1.D0
                              else
                                 write (*, *) 
     1'something wrong in color assignment'
                              endif
                              do j6 = 1, 4
                                 imp3(j6) = impul(j6,p3)
c                                            !momentum of the new fields
                                 imp(j6) = imp2(j6) + imp3(j6)
c                                   !array to store momenta combinations
                                 impul(j6,p1) = imp(j6)
                              end do
                              j9 = 0
                              do j6 = 1, 4
                                 j9 = j9 + 1
                                 d1(j9) = -imp(j6)
                              end do
c                           !via the proper trilinear interaction (oper)
                              call fuse (d1, imp2, offsh(psfl2), imp3, 
     1                           offsh(psfl3i), operot, fusion)
c                                !multiplying by the propagator -> polar
                              call prpgt (chcg(fl1), imp, fusion, polar)
                              if (nfr2 .gt. 0) then
c                             !returning proper fermi statistics in perm
                                 call permtn (prcss__nfrm, p2, p3, perm)
                              else
                                 perm = 1.D0
                              endif
c                        !returning the position, in OFFSHCOL, of field1
                              padd = psnew(p1)
                              if (padd .eq. 0) then
                                 psnew(p1) = poscol
                                 cpot1 = coupot*colcff*perm
                                 do j6 = 1, nl1
c                                                !storing field1
                                    offsh(j6-1+posfld) = cpot1*polar(j6)
                                 end do
                                 offshcol(1-1+poscol) = p1
                                 offshcol(2-1+poscol) = posfld
                                 offshcol(3-1+poscol) = clnw(1)
                                 offshcol(4-1+poscol) = clnw(2)
c    !number of new excitation of the same flavour and number of momenta
                                 nconf = nconf + 1
                                 poscol = poscol + 4
                                 posfld = posfld + nl1
c                              !new contribution to an old configuration
                              else
                                 cpot1 = coupot*colcff*perm
c                                           !position of field1 in OFFSH
                                 padd0 = offshcol(padd+1)
                                 do j6 = 1, nl1
c                                           !adding the new contribution
                                    offsh(j6-1+padd0) = offsh(j6-1+padd0
     1                                 ) + cpot1*polar(j6)
                                 end do
                              endif
                           endif
                        endif
    4                   continue
                     end do
                     pscfl2 = pscfl2 + 4
                     psfl2 = offshcol(pscfl2+1)
                  end do
               endif
            end do
         end do
         fldptr(1,chcg(fl1),iter) = nconf
      end do
      if (fststp.eq.1 .and. dbg) then
         fststp = 2
         if (poscol .gt. pscit) then
            pscit = poscol
            write (*, *) 'iteration,poscol,posfld', poscol, posfld
         endif
         if (posfld .gt. psfit) then
            psfit = posfld
            write (*, *) 'iteration,poscol,posfld', poscol, posfld
         endif
      endif
      if (poscol.gt.dim(2) .or. psfit.gt.dim(1)) then
         write (*, *) 
     1      'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
         write (*, *) 'iteration,poscol,posfld', poscol, posfld, dim(2)
     1      , dim(1)
         write (*, *) 'prcss', prcss__flv, prcss__hel, prcss__col, 
     1      prcss__nfrm, prcss__mom
         stop 
      endif
      amplitude = (0.,0.)
      do j6 = 1, 151
         tabmrgpr(j6) = tabint(j6,nprc)
      end do
      do j6 = 1, tabint(1,nprc)
         cp(j6) = cpam(j6,nprc)
         op(j6) = opam(j6,nprc)
      end do
      do j1 = 1, nvasam(1,nprc)
         cp(nvasam(j1+1,nprc)) = cp(nvasam(j1+1,nprc))*rnas
      end do
      do j6 = 1, nmom
         psnew(j6) = 0
c      !reporting the order of field1 momenta to avoid repeating the sam
         ordtmp(j6) = 0
      end do
      do j6 = 1, 44
         cnsold(j6-1) = .FALSE.
      end do
      steptm = 2
      flvold = -99
      do j1 = 1, tabmrgpr(1) + 1
         if (j1 .eq. tabmrgpr(1)+1) then
c    !last step compute the final contribution to the amplitude and exit
            flvold = 9999
         else
            fl1 = tabmrgpr(steptm)
            fl2 = tabmrgpr(steptm+1)
            fl3 = tabmrgpr(steptm+2)
            steptm = steptm + 3
            nl1 = nlordof(fl1)
            nl2 = nlordof(fl2)
            nl3 = nlordof(fl3)
            coupot = cp(j1)
            operot = op(j1)                      !oper(fl1,fl2,fl3)
         endif
         if (fl1 .ne. flvold) then
            if (fststp.eq.2 .and. j1.ne.1) then
               if (poscol.gt.dim(4) .or. posfld.gt.dim(3)) then
                  write (*, *) 
     1'ill dimensioned array OFFSH and/or OFFSHCOL in AMP'
                  write (*, *) 'amplitude,poscol,posfld', poscol, posfld
     1               , dim(4), dim(3), j1
                  stop 
               endif
               if (dbg) then
                  fststp = 2
                  if (poscol .gt. psca) then
                     psca = poscol
                     write (*, *) 'amplitude,poscol,posfld', poscol, 
     1                  posfld
                  endif
                  if (posfld .gt. psfa) then
                     psfa = posfld
                     write (*, *) 'amplitude,poscol,posfld', poscol, 
     1                  posfld
                  endif
               endif
            endif
            if (flvold .gt. 0) then
c                                         !flvold=9999 last contribution
               if (flvold .eq. 9999) flvold = fl1
c                                  !number of lorentz degrees of freedom
               nlold = nlordof(flvold)
               do iter = 1, nphf
                  nf1 = fldptr(1,flvold,iter)    !same as in iteration
                  if (nf1 .gt. 0) then
                     pscfl1 = fldptr(2,flvold,iter)
                     psfl1 = offshcol(pscfl1+1)
                     pref = 0                    !for optimization
                     do j2 = 1, nf1
                        psfl1 = offshcol(pscfl1+1)
                        p1 = offshcol(pscfl1)
                        p1 = nmom - p1
                        padd = psnew(p1)    !locating position of field1
c                                 !no 2,3 merging compatible with field1
                        if (padd .eq. 0) then
                           pscfl1 = pscfl1 + 4
                           go to 6
                        endif
                        if (p1 .ne. pref) then
                           pref = p1
                        else
                           padd = padd + 4 !there are two neutral gluons
                        endif
                        cl(1) = offshcol(pscfl1+2)
                        cl(2) = offshcol(pscfl1+3)
                        cl(3) = ampcol(padd+2)
                        cl(4) = ampcol(padd+3)
                        if(abs(cl(1)-cl(4))+abs(cl(2)-cl(3)).eq.0)then
                           padd = ampcol(padd+1)
                           x1 = 0
                           do j6 = 1, nlold
                              x1 = x1 + offsh(j6-1+psfl1)*ampfld(j6-1+
     1                           padd)
                           end do
                           tmp = x1
                           if (prcss__nfrm .ge. 2) then
c                             !returning proper fermi statistics in perm
                              call permtn(prcss__nfrm,nmom-p1,p1,perm)
                           else
                              perm = 1.D0
                           endif
                           amplitude = amplitude + tmp*perm
                           pscfl1 = pscfl1 + 4
                        endif
    6                   continue
                     end do
                  endif
               end do
               if (j1 .eq. tabmrgpr(1)+1) go to 7!end of computation
            endif
            flvold = fl1
            posfld = 1                           !pointer to offsh
            poscol = 1                           !pointer to offshcol
            if (cnsold(cnschg(fl1))) then
               do j6 = 1, nmom
c                         !none of these mergin has been computed before
                  psnew(j6) = 0
c!reporting the order of field1 momenta to avoid repeating the same comp
                  ordtmp(j6) = 0
               end do
            endif
            cnsold(cnschg(flvold)) = .TRUE.
         endif
         do iter = 1, nphf
            npmit = nprt - iter
            choose2 = iter.eq.nphf .and. even.eq.0
            nf1 = fldptr(1,fl1,iter)
            if (nf1 .gt. 0) then
               pscfl1 = fldptr(2,fl1,iter)
               if (fl1 .eq. fl2) then
                  beg = iter
               else
                  beg = 1
               endif
               do ex2 = beg, nphf
                  choose = fl1.eq.fl2 .and. iter.eq.ex2
                  ex3 = nprt - iter - ex2
                  if (ex3.ge.1 .and. ex3.le.nphf) then
c                                      !again to deal with equal flavour
                     if (fl2.ne.fl3 .or. ex3.ge.ex2) then
                        nf2 = fldptr(1,fl2,ex2)
                        nf3 = fldptr(1,fl3,ex3)
                        if (nf2.gt.0 .and. nf3.gt.0) then
                           pscfl2 = fldptr(2,fl2,ex2)
                           pscfl3 = fldptr(2,fl3,ex3)
                           choose3 = ex2.eq.ex3 .and. fl2.eq.fl3
                           j5 = pscfl1
                           do j4 = 1, nf1
c                 !now ordtmp register the ordering of momenta of field1
                              ordtmp(offshcol(j5)) = j4
                              j5 = j5 + 4
                           end do
                           do j3 = 1, nf2
                              cl(1) = offshcol(pscfl2+2)
                              cl(2) = offshcol(pscfl2+3)
                              psfl2 = offshcol(pscfl2+1)
                              if (choose3) then
                                 beg3 = j3 + 1
                                 pscfl3i = pscfl3 + 4*(beg3 - 2)
                              else
                                 beg3 = 1
                                 pscfl3i = pscfl3 - 4
                              endif
                              p2 = offshcol(pscfl2)
                              do j6 = 1, 4
                                 imp2(j6) = impul(j6,p2)
                              end do
                              if (prcss__nfrm .gt. 2) then
                                 nfr2 = 0
                                 do j4 = 0, prcss__nfrm - 1
                                    if (btest(p2,j4)) nfr2 = nfr2 + 1
                                 end do
                              else
                                 nfr2 = 0
                              endif
                              do j4 = beg3, nf3
                                 pscfl3i = pscfl3i + 4
                                 p3 = offshcol(pscfl3i)
                                 p1 = p2 + p3
                                 if (momd(p1) .eq. npmit) then
                                    ot = ordtmp(nmom-p1)
c!field1 does not contain the frequency to match the proposed 2,3 fusion
                                    if (ot .ne. 0) then
c                         !to avoid repeating the same computation twice
                                    if (.not.(choose .and. ordtmp(p2)
     1                                 .le.ot)) then
                                    if (.not.(choose2 .and. btest(p1,0))
     1                                 ) then
                                    cl(3) = offshcol(pscfl3i+2)
                                    cl(4) = offshcol(pscfl3i+3)
c                                                !quarks and gluons
                                    if (colrep(fl1) .eq. 3) then
                                    if (cl(2) .eq. 0) go to 10
                                    if (cl(2) .eq. cl(3)) then
                                    if (cl(1) .ne. 0) then
                                    if (cl(1) .eq. cl(4)) go to 10
                                    colcff = 1.D0
                                    else
                                    colcff = 1.D0
                                    endif
                                    clnw(1) = cl(1)
                                    clnw(2) = cl(4)
                                    else if (cl(1) .eq. cl(4)) then
                                    if (cl(1) .ne. 0) then
                                    colcff = -1.D0
                                    else
                                    colcff = 1.D0
                                    endif
                                    clnw(1) = cl(3)
                                    clnw(2) = cl(2)
                                    else if (colrep(fl2).eq.1 .and. 
     1                                 colrep(fl3).eq.3) then
                                    colcff = 1.D0
                                    do j6 = 1, 2
                                    clnw(j6) = cl(j6+2)
                                    end do
                                    else if (colrep(fl2).eq.3 .and. 
     1                                 colrep(fl3).eq.1) then
                                    colcff = 1.D0
                                    do j6 = 1, 2
                                    clnw(j6) = cl(j6)
                                    end do
                                    else
                                    go to 10
                                    endif
c                                                !coulorless object
                                    else if (colrep(fl1) .eq. 2) then
                                    if (cl(2) .eq. cl(3)) then
                                    clnw(1) = cl(1)
                                    clnw(2) = cl(4)
                                    else if (cl(1) .eq. cl(4)) then
                                    clnw(1) = cl(3)
                                    clnw(2) = cl(2)
                                    else
                                    go to 10
                                    endif
                                    colcff = 1.D0
c                                                !coulorless object
                                    else if (colrep(fl1) .eq. 1) then
                                    if (cl(2).eq.cl(3) .and. cl(1).eq.cl
     1                                 (4)) then
                                    clnw(1) = 0
                                    clnw(2) = 0
                                    else
                                    go to 10
                                    endif
                                    colcff = 1.D0
                                    else
                                    write (*, *) 
     1'something wrong in color assignment'
                                    endif
                                    do j6 = 1, 4
                                    imp3(j6) = impul(j6,p3)
                                    imp(j6) = imp2(j6) + imp3(j6)
                                    end do
                                    psfl3i = offshcol(pscfl3i+1)
                                    j10 = 0
                                    do j6 = 1, 4
                                    j10 = j10 + 1
                                    d2(j10) = -imp(j6)
                                    end do
                                    call fuse (d2, imp2, offsh(psfl2), 
     1                                 imp3, offsh(psfl3i), operot, 
     2                                 fusion)
                                    if (nfr2 .gt. 0) then
                                    call permtn (prcss__nfrm, p2, p3, 
     1                                 perm)
                                    else
                                    perm = 1.D0
                                    endif
                                    padd = psnew(p1)
                                    if (padd .eq. 0) then
                                    psnew(p1) = poscol
                                    cpot1 = coupot*colcff*perm
                                    do j6 = 1, nl1
                                    ampfld(j6-1+posfld) = cpot1*fusion(
     1                                 j6)
                                    end do
                                    ampcol(1-1+poscol) = p1
                                    ampcol(2-1+poscol) = posfld
                                    ampcol(3-1+poscol) = clnw(1)
                                    ampcol(4-1+poscol) = clnw(2)
                                    nconf = nconf + 1
                                    posfld = posfld + nl1
                                    poscol = poscol + 4
c                              !new contribution to an old configuration
                                    else
                                    padd0 = ampcol(padd+1)
                                    cpot1 = coupot*colcff*perm
                                    do j6 = 1, nl1
                                    ampfld(j6-1+padd0) = ampfld(j6-1+
     1                                 padd0) + cpot1*fusion(j6)
                                    end do
                                    endif
                                    endif
                                    endif
                                    endif
                                 endif
   10                            continue
                              end do
                              pscfl2 = pscfl2 + 4
                           end do
                        endif
                     endif
                  endif
               end do
            endif
         end do
    7    continue
      end do
      mtel = amplitude
!
      end 
!***************************************************************************
      subroutine sourcemassboson(spinsour, fieldaux, mom, mass)
!***************************************************************************
!
! This subroutine given the four momentum of a massive boson MOM the required
! source polarization SPINSOUR and the boson mass MASS returns in
! FIELDAUX the source term
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      integer n                                  !loop index
!
      double precision knorm2     !squared modulus of the three momentum
      double precision mass                      !massive boson mass
      double precision massa2                    !squared mass
      double precision mom(4)                    !four momentum
      double precision nz(4), nx(4)
      save nz, nx
c                                                !longitudinal (1)
      double precision sour1(4), sour2(4), sour3(4)
!                                and transverse (2,3) polarizations
      integer spinsour                           !required polarization
      double precision scalar3
      double precision xnorm                   !to normalize the sources
      double complex fieldaux(4)    !to return the relevant polarization
      data nz/ 0., 0., 0., 1./ 
      data nx/ 0., 1., 0., 0./ 
      data sour1/ 4*0./ 
      data sour2/ 4*0./ 
      data sour3/ 4*0./ 
      doubleprecision d1, d2, d3
!
!
!
!
      knorm2 = mom(2)**2 + mom(3)**2 + mom(4)**2
!
      if (knorm2 .eq. 0) then
         sour1(2) = 1.
         sour2(3) = 1.
         sour3(4) = 1.
         return 
      endif
      massa2 = mass**2
!
!    pol. longitudinal
!
      do n = 2, 4
         sour1(n) = mom(n)
      end do
      sour1(1) = knorm2/mom(1)
      d1 = scalar3(sour1(1),sour1(1))
      xnorm = sqrt(d1)
      if (xnorm .ne. 0.) then
         do n = 1, 4
            sour1(n) = sour1(n)/xnorm
         end do
      endif
!
!    pol. transverse +
!
      if (mom(4)**2 .gt. mom(2)**2) then
         call vector3prod (mom, nx, sour2)
      else
         call vector3prod (mom, nz, sour2)
      endif
      d2 = scalar3(sour2(1),sour2(1))
      xnorm = sqrt(d2)
      do n = 1, 4
         sour2(n) = sour2(n)/xnorm
      end do
!
!    pol. transverse -
!
      call vector3prod (mom, sour2, sour3)
      d3 = scalar3(sour3(1),sour3(1))
      xnorm = sqrt(d3)
      do n = 1, 4
         sour3(n) = sour3(n)/xnorm
      end do
!
      do n = 1, 4
         fieldaux(n) = 0.
         if (spinsour .eq. 0) fieldaux(n) = sour1(n)
         if (spinsour .eq. 1) fieldaux(n) = sour2(n)
         if (spinsour .eq. (-1)) fieldaux(n) = sour3(n)
      end do
!
      return 
      end 
!*********************************************************************
      subroutine vector3prod(v1, v2, vfin)
!*********************************************************************
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision v1(4)
      double precision v2(4)
      double precision vfin(4)
!
      vfin(2) = v1(3)*v2(4) - v1(4)*v2(3)
      vfin(3) = v1(4)*v2(2) - v1(2)*v2(4)
      vfin(4) = v1(2)*v2(3) - v1(3)*v2(2)
      vfin(1) = 0
!
      return 
      end 
!***********************************************************************
      function scalar3 (v1, v2)
!***********************************************************************
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision scalar3
      double precision v1(4)
      double precision v2(4)
!
      scalar3 = abs(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
!
      return 
      end 
!***********************************************************************
      subroutine fermionsourceshl(spinsour, fieldaux, mom)
!***********************************************************************
!
! This subroutine return the initialized fermion field.
! Elicity heigenstates. Massless fermions
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mom(4)          !momenta of the external particle
      double precision p, p3p, p3m, coeffp, coeffm
!
c                         !array returning the fermion field configurati
      double complex fieldaux(2)
      double complex im      !immaginary unity in complex representation
      save im
      double complex p1p
      data im/ (0.,1.)/ 
      integer spinsour
      double precision usd
!
! Static variables
!
!
!
!
      p = sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
      p3p = p + mom(4)
      p3m = p - mom(4)
      p1p = mom(2) + mom(3)*im
!
      if (abs(spinsour) .eq. 1) spinsour = -spinsour
      if (abs(p3m).lt.1.D-10 .or. abs(p3p).lt.1.D-10) then
         call fshl (spinsour, fieldaux, mom)
         return 
      endif
!
      coeffp = 1./sqrt(p3p)
      coeffm = 1./sqrt(p3m)
!
! "spin up" ingoing fermion
!
      if (spinsour .eq. 1) then
         fieldaux(1) = coeffp*p3p
         fieldaux(2) = coeffp*p1p
      endif
!
! "spin down" ingoing fermion
!
      if (spinsour .eq. (-1)) then
         fieldaux(1) = coeffm*p3m
         fieldaux(2) = -coeffm*p1p
      endif
!
! "spin up" outgoing antifermion
!
      if (spinsour .eq. 2) then
         fieldaux(1) = -coeffm*p3m
         fieldaux(2) = coeffm*p1p
      endif
!
! "spin down" outgoing antifermion
!
      if (spinsour .eq. (-2)) then
         fieldaux(1) = coeffp*p3p
         fieldaux(2) = coeffp*p1p
      endif
!
      return 
      end 
 
!***********************************************************************
      subroutine fermionbarsourceshl(spinsour, fieldaux, mom)
!***********************************************************************
!
! This subroutine return the initialized fermionbar field.
! Elicity heigenstates. Massless fermion
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
c                          !array containing the momenta of the external
      double precision mom(4)
      double precision p, p3p, p3m, coeffp, coeffm
!
c                            !array returning the fermion field configur
      double complex fieldaux(2)
      double complex im      !immaginary unity in complex representation
      save im
      double complex p1p
      data im/ (0.,1.)/ 
      integer spinsour
!
! Static variables
!
!
!
      p = sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
      p3p = p + mom(4)
      p3m = p - mom(4)
      p1p = mom(2) - mom(3)*im
!
      if (abs(spinsour) .eq. 1) spinsour = -spinsour
      if (abs(p3m).lt.1.D-10 .or. abs(p3p).lt.1.D-10) then
         call fbshl (spinsour, fieldaux, mom)
         return 
      endif
!
      coeffp = 1./sqrt(p3p)
      coeffm = 1./sqrt(p3m)
!
! "spin up" ingoing fermion
!
      if (spinsour .eq. 1) then
         fieldaux(1) = coeffp*p3p
         fieldaux(2) = coeffp*p1p
      endif
!
! "spin down" ingoing fermion
!
      if (spinsour .eq. (-1)) then
         fieldaux(1) = coeffm*p3m
         fieldaux(2) = -coeffm*p1p
      endif
!
! "spin up" outgoing antifermion
!
      if (spinsour .eq. 2) then
         fieldaux(1) = -coeffm*p3m
         fieldaux(2) = coeffm*p1p
      endif
!
! "spin down" outgoing antifermion
!
      if (spinsour .eq. (-2)) then
         fieldaux(1) = coeffp*p3p
         fieldaux(2) = coeffp*p1p
      endif
!
      return 
      end 
!***********************************************************************
      subroutine fshl(spinsour, fieldaux, mom)
!***********************************************************************
!
! This subroutine return the initialized fermion field. Massless fermion
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
c                      !array containing the momenta of the external par
      double precision mom(4)
      integer spinsour          !array containing the spin of the source
!
c                            !array returning the fermion field configur
      double complex fieldaux(2)
      double complex im      !immaginary unity in complex representation
      save im
      data im/ (0.,1.)/ 
!
! Static variables
!
!
!
! "spin up" ingoing fermion
!
      if (mom(4) .lt. 0) spinsour = -spinsour
      if (spinsour .eq. 1) then
         fieldaux(1) = sqrt(2.D0*mom(1))
         fieldaux(2) = 0.D0
!
! "spin down" ingoing fermion
!
      else if (spinsour .eq. (-1)) then
         fieldaux(1) = 0.D0
         fieldaux(2) = sqrt(2.D0*mom(1))
!
! "spin up" outgoing antifermion
!
      else if (spinsour .eq. 2) then
         fieldaux(1) = sqrt(2.D0*mom(1))*mom(4)/abs(mom(4))
         fieldaux(2) = 0.D0
!
! "spin down" outgoing antifermion
!
      else if (spinsour .eq. (-2)) then
         fieldaux(1) = 0.D0
         fieldaux(2) = -sqrt(2.D0*mom(1))*mom(4)/abs(mom(4))
      endif
!
      return 
      end 
!***********************************************************************
      subroutine fbshl(spinsour, fieldaux, mom)
!***********************************************************************
!
! This subroutine return the initialized fermion field.
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mass!containing the mass of the external particle
c                     !array containing the momenta of the external part
      double precision mom(4)
      integer spinsour          !array containing the spin of the source
!
c                          !array returning the fermionbar field configu
      double complex fieldaux(2)
      double complex im      !immaginary unity in complex representation
      save im
      data im/ (0.,1.)/ 
!
! Static variables
!
!
!
! "spin up" outgoing fermion
!
      if (mom(4) .gt. 0) spinsour = -spinsour
      if (spinsour .eq. 1) then
         fieldaux(1) = sqrt(2.D0*mom(1))
         fieldaux(2) = 0.D0
!
! "spin down" outgoing fermion
!
      else if (spinsour .eq. (-1)) then
         fieldaux(2) = sqrt(2.D0*mom(1))
         fieldaux(1) = 0.D0
!
! "spin up" ingoing antifermion
!
      else if (spinsour .eq. 2) then
         fieldaux(1) = -sqrt(2.D0*mom(1))*mom(4)/abs(mom(4))
         fieldaux(2) = 0.D0
!
! "spin down" ingoing antifermion
!
      else if (spinsour .eq. (-2)) then
         fieldaux(2) = sqrt(2.D0*mom(1))*mom(4)/abs(mom(4))
         fieldaux(1) = 0.D0
      endif
!
      return 
      end 
!***********************************************************************
      subroutine fermionsources(spinsour, fieldaux, mom, mass)
!***********************************************************************
!
! This subroutine return the initialized fermion field.
! Elicity heigenstates
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mass!containing the mass of the external particle
c                          !array containing the momenta of the external
      double precision mom(4)
      double precision p, p3p, p3m, mp, coeffp, coeffm
!
c                         !array returning the fermion field configurati
      double complex fieldaux(4)
      double complex im      !immaginary unity in complex representation
      save im
      double complex p1p
      data im/ (0.,1.)/ 
      integer spinsour
      double precision usd
!
! Static variables
!
!
!
!
      usd = 1.D0/sqrt(2.D0)
      p = sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
      p3p = p + mom(4)
      p3m = p - mom(4)
      mp = mass + mom(1)
      p1p = mom(2) + mom(3)*im
!
      if (abs(spinsour) .eq. 1) spinsour = -spinsour
      if (abs(p3m).lt.1.D-10 .or. abs(p3p).lt.1.D-10) then
         call fs (spinsour, fieldaux, mom, mass)
         return 
      endif
!
      coeffp = 1./sqrt(2.*p*mp*p3p)
      coeffm = 1./sqrt(2.*p*mp*p3m)
!
! "spin up" ingoing fermion
!
      if (spinsour .eq. 1) then
         fieldaux(1) = usd*coeffp*p3p*(mp - p)
         fieldaux(2) = usd*coeffp*p1p*(mp - p)
         fieldaux(3) = usd*coeffp*p3p*(p + mp)
         fieldaux(4) = usd*coeffp*p1p*(p + mp)
      endif
!
! "spin down" ingoing fermion
!
      if (spinsour .eq. (-1)) then
         fieldaux(1) = usd*coeffm*p3m*(mp + p)
         fieldaux(2) = -usd*coeffm*p1p*(mp + p)
         fieldaux(3) = -usd*coeffm*p3m*(p - mp)
         fieldaux(4) = usd*coeffm*p1p*(p - mp)
      endif
!
! "spin up" outgoing antifermion
!
      if (spinsour .eq. 2) then
         fieldaux(1) = -usd*coeffm*p3m*(p + mp)
         fieldaux(2) = usd*coeffm*p1p*(p + mp)
         fieldaux(3) = usd*coeffm*p3m*(mp - p)
         fieldaux(4) = -usd*coeffm*p1p*(mp - p)
      endif
!
! "spin down" outgoing antifermion
!
      if (spinsour .eq. (-2)) then
         fieldaux(1) = usd*coeffp*p3p*(p - mp)
         fieldaux(2) = usd*coeffp*p1p*(p - mp)
         fieldaux(3) = usd*coeffp*p3p*(mp + p)
         fieldaux(4) = usd*coeffp*p1p*(mp + p)
      endif
!
      return 
      end 
 
!***********************************************************************
      subroutine fermionbarsources(spinsour, fieldaux, mom, mass)
!***********************************************************************
!
! This subroutine return the initialized fermionbar field.
! Elicity heigenstates
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mass!containing the mass of the external particle
c                          !array containing the momenta of the external
      double precision mom(4)
      double precision p, p3p, p3m, mp, coeffp, coeffm
!
c                            !array returning the fermion field configur
      double complex fieldaux(4)
      double complex im      !immaginary unity in complex representation
      save im
      double complex p1p
      data im/ (0.,1.)/ 
      integer spinsour
!
! Static variables
!
!
!
      p = sqrt(mom(2)**2+mom(3)**2+mom(4)**2)
      p3p = p + mom(4)
      p3m = p - mom(4)
      mp = mass + mom(1)
      p1p = mom(2) - mom(3)*im
!
      if (abs(spinsour) .eq. 1) spinsour = -spinsour
      if (abs(p3m).lt.1.D-10 .or. abs(p3p).lt.1.D-10) then
         call fbs (spinsour, fieldaux, mom, mass)
         return 
      endif
!
      coeffp = 1./sqrt(p*mp*p3p)/2.D0
      coeffm = 1./sqrt(p*mp*p3m)/2.D0
!
! "spin up" ingoing fermion
!
      if (spinsour .eq. 1) then
         fieldaux(1) = coeffp*p3p*(mp + p)
         fieldaux(2) = coeffp*p1p*(mp + p)
         fieldaux(3) = -coeffp*p3p*(p - mp)
         fieldaux(4) = -coeffp*p1p*(p - mp)
      endif
!
! "spin down" ingoing fermion
!
      if (spinsour .eq. (-1)) then
         fieldaux(1) = coeffm*p3m*(mp - p)
         fieldaux(2) = -coeffm*p1p*(mp - p)
         fieldaux(3) = coeffm*p3m*(p + mp)
         fieldaux(4) = -coeffm*p1p*(p + mp)
      endif
!
! "spin up" outgoing antifermion
!
      if (spinsour .eq. 2) then
         fieldaux(1) = -coeffm*p3m*(p - mp)
         fieldaux(2) = coeffm*p1p*(p - mp)
         fieldaux(3) = -coeffm*p3m*(mp + p)
         fieldaux(4) = coeffm*p1p*(mp + p)
      endif
!
! "spin down" outgoing antifermion
!
      if (spinsour .eq. (-2)) then
         fieldaux(1) = coeffp*p3p*(p + mp)
         fieldaux(2) = coeffp*p1p*(p + mp)
         fieldaux(3) = -coeffp*p3p*(mp - p)
         fieldaux(4) = -coeffp*p1p*(mp - p)
      endif
!
      return 
      end 
!***********************************************************************
      subroutine fs(spinsour, fieldaux, mom, mass)
!***********************************************************************
!
! This subroutine return the initialized fermion field.
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mass!containing the mass of the external particle
c                      !array containing the momenta of the external par
      double precision mom(4)
      integer spinsour          !array containing the spin of the source
!
c                            !array returning the fermion field configur
      double complex fieldaux(4)
      double complex im
      double complex f1
      double complex f2
      double complex f3
      double complex f4      !immaginary unity in complex representation
      save im
      data im/ (0.,1.)/ 
!
! Static variables
!
!
!
! "spin up" ingoing fermion
!
      if (mom(4) .lt. 0) spinsour = -spinsour
      if (spinsour .eq. 1) then
         f1 = sqrt(mom(1)+mass)
         f2 = 0.
         f3 = mom(4)/sqrt(mass + mom(1))
         f4 = mom(2)/sqrt(mass + mom(1)) + mom(3)/sqrt(mass + mom(1))*im
!
! "spin down" ingoing fermion
!
      else if (spinsour .eq. (-1)) then
         f1 = 0.
         f2 = sqrt(mom(1)+mass)
         f3 = mom(2)/sqrt(mass + mom(1)) - mom(3)/sqrt(mass + mom(1))*im
         f4 = -mom(4)/sqrt(mass + mom(1))
!
! "spin up" outgoing antifermion
!
      else if (spinsour .eq. 2) then
         f1 = mom(4)/sqrt(mass + mom(1))
         f2 = mom(2)/sqrt(mass + mom(1)) + mom(3)/sqrt(mass + mom(1))*im
         f3 = sqrt(mom(1)+mass)
         f4 = 0.
!
! "spin down" outgoing antifermion
!
      else if (spinsour .eq. (-2)) then
         f1 = mom(2)/sqrt(mass + mom(1)) - mom(3)/sqrt(mass + mom(1))*im
         f2 = -mom(4)/sqrt(mass + mom(1))
         f3 = 0.
         f4 = sqrt(mom(1)+mass)
      endif
!
      fieldaux(1) = sqrt(1.D0/2.D0)*(f1 - f3)
      fieldaux(2) = sqrt(1.D0/2.D0)*(f2 - f4)
      fieldaux(3) = sqrt(1.D0/2.D0)*(f1 + f3)
      fieldaux(4) = sqrt(1.D0/2.D0)*(f2 + f4)
!
      return 
      end 
!***********************************************************************
      subroutine fbs(spinsour, fieldaux, mom, mass)
!***********************************************************************
!
! This subroutine return the initialized fermion field.
!
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
!
      double precision mass!containing the mass of the external particle
c                     !array containing the momenta of the external part
      double precision mom(4)
      integer spinsour          !array containing the spin of the source
!
c                          !array returning the fermionbar field configu
      double complex fieldaux(4)
      double complex im
      double complex f1
      double complex f2
      double complex f3
      double complex f4      !immaginary unity in complex representation
      save im
      data im/ (0.,1.)/ 
!
! Static variables
!
!
!
! "spin up" outgoing fermion
!
      if (mom(4) .gt. 0) spinsour = -spinsour
      if (spinsour .eq. 1) then
         f1 = sqrt(mom(1)+mass)
         f2 = 0.
         f3 = -mom(4)/sqrt(mass + mom(1))
         f4=(-mom(2)/sqrt(mass+mom(1)))+mom(3)/sqrt(mass+mom(1))*im
!
! "spin down" outgoing fermion
!
      else if (spinsour .eq. (-1)) then
         f1 = 0.
         f2 = sqrt(mom(1)+mass)
         f3=(-mom(2)/sqrt(mass+mom(1)))-mom(3)/sqrt(mass+mom(1))*im
         f4 = mom(4)/sqrt(mass + mom(1))
!
! "spin up" ingoing antifermion
!
      else if (spinsour .eq. 2) then
         f1 = -mom(4)/sqrt(mass + mom(1))
         f2=(-mom(2)/sqrt(mass+mom(1)))+mom(3)/sqrt(mass+mom(1))*im
         f3 = sqrt(mom(1)+mass)
         f4 = 0.
!
! "spin down" ingoing antifermion
!
      else if (spinsour .eq. (-2)) then
         f1=(-mom(2)/sqrt(mass+mom(1)))-mom(3)/sqrt(mass+mom(1))*im
         f2 = mom(4)/sqrt(mass + mom(1))
         f3 = 0.
         f4 = sqrt(mom(1)+mass)
      endif
!
      fieldaux(1) = sqrt(1.D0/2.D0)*(f1 - f3)
      fieldaux(2) = sqrt(1.D0/2.D0)*(f2 - f4)
      fieldaux(3) = sqrt(1.D0/2.D0)*(f1 + f3)
      fieldaux(4) = sqrt(1.D0/2.D0)*(f2 + f4)
!
      return 
      end 
 
 
 
!
      subroutine inital
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
      parameter (nintmax = 50)
      integer tabint(151,1)
      integer tabintmrg(451,1)
      integer opam(50,1)
      integer opfs(150,1)
      integer nvasfs(151,1)
      integer nvasam(51,1)
      double precision cpam(50,1)
      double precision cpfs(150,1)
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
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
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
      integer tabint(151,1)
      integer tabintmrg(451,1)
      integer opam(50,1)
      integer opfs(150,1)
      integer nvasfs(151,1)
      integer nvasam(51,1)
      double precision cpam(50,1)
      double precision cpfs(150,1)
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
 
 
 
 
 
 
