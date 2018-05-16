c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      real*8 pin,pout,pjet,
     + plp,plm,pnu,pnub,pz,ptj,etaj,
     + etalp,etalm,ptlp,ptlm,drjj,
     + drlj
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar-2),
     + plp(4),plm(4),pnu(4),pnub(4),pz(4),ptj(maxpar-2),
     + ptlp,ptlm,     
     + etaj(maxpar-2),etalp,etalm,
     + drjj(maxpar-2,maxpar-2),
     + drlj(maxpar-2)
