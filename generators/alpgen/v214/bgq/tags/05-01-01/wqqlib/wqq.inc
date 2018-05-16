c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      real*8 pin,pout,pjet,
     + pbott,pbbar,plep,pnu,pw,ptj,ptb,ptbb,etaj,etab,etabb,
     + etalep,ptlep,ptmiss,drjj,
     + drbj,drbbj,drbb,drlj,drlb,drlbb,pbjet,ptbjet,etabjet,
     + wchrg
      integer nfspart
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar-2),
     + pbott(4),pbbar(4),plep(4),pnu(4),pw(4),ptj(maxpar-2),ptb,
     + ptbb,ptlep,ptmiss,     
     + etab,etabb,etaj(maxpar-2),etalep,
     + drjj(maxpar-2,maxpar-2),
     + drbj(maxpar),drbbj(maxpar),drbb,drlj(maxpar-2),drlb,drlbb,
     + pbjet(4),ptbjet,etabjet,wchrg,nfspart
