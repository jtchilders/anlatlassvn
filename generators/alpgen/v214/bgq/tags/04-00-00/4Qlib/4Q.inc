c event kinematics, ordered according to generation ordering
      double precision pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      double precision pin,pout,pjet,
     + pbott,pbbar,ptop,ptopbar,ptj,ptb,ptbb,etaj,etab,etabb,
     + drjj,drbj,drbbj,drbb
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar),
     + pbott(4,3),pbbar(4,3),ptop(4,3),ptopbar(4,3),
     + ptj(maxpar),ptb(3),ptbb(3),     
     + etab(3),etabb(3),etaj(maxpar),
     + drjj(maxpar,maxpar),
     + drbj(maxpar,maxpar),drbbj(maxpar,maxpar),drbb(15)
