c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      integer nfspart
      real*8 pin,pout,pjet,
     + pbott,pbbar,plp,plm,pnub,pnu,pz,ptj,ptb,ptbb,etaj,etab,etabb,
     + etalp,etalm,etanub,etanu,ptlp,ptlm,ptnub,ptnu,drjj,
     + drbj,drbbj,drbb,drlj,drlb,drlbb,pbjet,ptbjet,etabjet
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar-2),
     + pbott(4),pbbar(4),plp(4),plm(4),pnu(4),pnub(4),pz(4),
     + ptj(maxpar-2),ptb,ptbb,ptlp,ptlm,ptnu,ptnub,     
     + etab,etabb,etaj(maxpar-2),etalp,etalm,etanub,etanu,
     + drjj(maxpar-2,maxpar-2),
     + drbj(maxpar),drbbj(maxpar),drbb,drlj(maxpar-2),drlb,drlbb,
     + pbjet(4),ptbjet,etabjet,nfspart
