c event kinematics, ordered according to generation ordering
      double precision pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      double precision pin,pout,pjet,phiggs,ptop,ptbar,pbott,
     + pbbar,ptb,ptbb,etab,etabb,drbb,ptj,etaj,drjj
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar),
     + phiggs(4),ptop(4),ptbar(4),pbott(4),pbbar(4),
     + ptb,ptbb,etab,etabb,drbb,
     + ptj(maxpar),etaj(maxpar),drjj(maxpar,maxpar)
c
c   for t tbar decay
c
      integer jtl,jtbl
      character*1 idecay
      common/tdecflag/idecay
      common/topind/jtl,jtbl
      real*8 idec(4,3,2)
      common/decay/idec     ! idec(1:4,i,j) four momentum of the I-th 
C            decay product of the J-th particle (jtl=top, jtbl=topbar); 
C            the four momentum is (p_x,p_y,p_z,E)
C            for the top(topbr) I=1 is b (bbar) I=2 is fermion 
C            I=3 is fermionbar (W-> f fb) 
