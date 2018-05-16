c
c event kinematics, ordered according to generation ordering
      double precision pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      double precision pin,pout,pjet,ptj,etaj,drjj
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar),
     + ptj(maxpar),etaj(maxpar),drjj(maxpar,maxpar)
