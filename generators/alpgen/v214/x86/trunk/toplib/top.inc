c
c event kinematics, ordered according to generation ordering
      double precision pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c
c event kinematics, explicit labelling of flavours, leptons, etc
      double precision pin,pout,pjet,
     + pbott,pbbar,ptop,ptbar,ptj,ptb,ptbb,etaj,etab,etabb,drjj,
     + drbj,drbbj,drbb,pbjet,ptbjet,etabjet
      common/usrevt/pin(4,2),pout(4,maxpar-2),pjet(4,maxpar),
     + pbott(4),pbbar(4),ptop(4),ptbar(4),ptj(maxpar),ptb,ptbb,
     + etab,etabb,etaj(maxpar),drjj(maxpar,maxpar),
     + drbj(maxpar),drbbj(maxpar),drbb,
     + pbjet(4),ptbjet,etabjet
c
c   for t tbar decay
c
      integer jtl,jtbl
      character*1 idecay
      common/tdecflag/idecay
      real*8 idec(4,3,2)
      common/decay/idec     ! idec(1:4,i,j) four momentum of the I-th 
C            decay product of the J-th particle (jtl=top, jtbl=topbar); 
C            the four momentum is (p_x,p_y,p_z,E)
C            for the top(topbr) I=1 is b (bbar) I=2 is fermion 
C            I=3 is fermionbar (W-> f fb) 
c
c  for W decay
c
      real*8 wdec(4,2)
      common/decay/wdec     ! wdec(1:4,i) four momentum of the I-th 
C            decay product of the W; 
C            the four momentum is (p_x,p_y,p_z,E)
C            I=1 is fermion I=2 is fermionbar 
C

