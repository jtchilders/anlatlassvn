C *** MODULE couplings ***
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 12:11:27 7/12/05
C...Switches:                     
      INCLUDE 'Vincnst.inc'
!
c                                      ! MASSES(j) mass of j-th particle
      double precision masses(43), width(43)
!
c                       ! OPER(j,k,l) labels the operator to be used for
      integer oper(43,43,43)
!                                                        particles j,k and l interaction
c                        ! COUP(j,k,l) coupling constant for OPER(j,k,l)
      double precision coup(43,43,43)
!
      double precision gstrong, g2weak, gbar, ctw, stw, glup, gldn, gllp
     1   , grup, grdn, grlp, emch, trihcp, yuk! strong coupling constant
c                 !for quartic gauge couplings: 1-3==wp  4-6==wm  7-9==z
      double precision qrtc(3,4), qrtca(3,4), qrtc2(2,4), qrtc3(2,4)
c                           !for quartic higgs and higgs-gauge couplings
      double precision qrtchh(3), qrtchg(2)
!
      double precision ggh            ! effective higgs-glu-glu coupling
      common /couplings/masses, width, coup, gstrong, g2weak, gbar, ctw
     1   , stw, glup, gldn, gllp, grup, grdn, grlp, emch, trihcp, yuk, 
     2   qrtc, qrtca, qrtc2, qrtc3, qrtchh, qrtchg, ggh, oper
C *** END MODULE couplings ***
