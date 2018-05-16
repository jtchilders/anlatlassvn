c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      real*8 pin,pout,pjet,
     + plep,pnu,pw,ptj,etaj,
     + etalep,ptlep,ptmiss,drjj,drlj,
     + pphot,ptphot,etaphot,drjphot,drphotphot,drlphot
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar-2),
     + plep(4),pnu(4),pw(4),ptj(maxpar-2),
     + ptlep,ptmiss,etaj(maxpar-2),etalep,
     + drjj(maxpar-2,maxpar-2),drlj(maxpar-2),
     + pphot(4,maxpar-2),ptphot(maxpar-2),etaphot(maxpar-2),
     + drjphot(maxpar-2,maxpar-2),drphotphot(maxpar-2,maxpar-2),
     + drlphot(maxpar-2)



