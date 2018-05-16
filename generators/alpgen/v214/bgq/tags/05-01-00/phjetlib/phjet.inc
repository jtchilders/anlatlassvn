c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      real*8 pin,pout,pjet,pphot,ptj,etaj,ptp,etap,drjj,drpj,drpp
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar-2),
     +              pphot(4,maxpar-2),ptj(maxpar-2),
     +              etaj(maxpar-2),ptp(maxpar-2),etap(maxpar-2),
     +              drjj(maxpar-2,maxpar-2),
     +              drpj(maxpar-2,maxpar-2),drpp(maxpar-2,maxpar-2)
