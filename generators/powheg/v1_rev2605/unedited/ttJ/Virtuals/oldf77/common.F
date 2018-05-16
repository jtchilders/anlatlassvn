C$Modified: Sun Aug 13 02:00:59 2006 by puwer $
      double precision N,Ca,CF,TR
      double precision kappa
      double precision alphas
      double precision g2
      parameter (kappa=2.d0/3.d0) ! Value used by Stefan D.
C Stefan D's alpha_s(mt,1-loop)
*      parameter (alphas=0.11773621487454d0)
C Stefan D's alpha_s(mt,2-loop)
C      parameter (alphas=0.10752054927347d0)
C Stefan W's alpha_s
**      parameter (alphas=0.1082647127128654d0)
*      parameter (g2=4.d0*pi*alphas)
      double precision mu_ren,mu_fac,sqrt_s,s,mt,shad,sqrt_shad
     #    
      common/settings/g2,alphas,mu_ren,mu_fac,sqrt_s,s,sqrt_shad
     #     ,shad,mt
      parameter (N=3.d0, CA=N, CF = 0.5d0/N*(N**2-1.d0), TR =0.5d0)
**********************************************************************
      double precision etavalues(100), rhovalues(100)
      integer nval 
      common/points/etavalues,rhovalues,nval
**********************************************************************
      integer PDFSET,CTEQ6L,CTEQ6L1,CTEQ6M,MRST2003cNLO,MRST2003cNNLO
      integer CTEQ6M1
      parameter (CTEQ6L=3,CTEQ6L1=4,CTEQ6M=1,CTEQ6M1=200,
     &     MRST2003cNLO=10,MRST2003cNNLO=11)
      common/pdfstuff/pdfset
**********************************************************************
*      parameter (mt=174.d0)
Color matrices:
      double precision cg1g2(6,6),cg1g3(6,6),cg2g3(6,6),cg1q(6,6),
     >     cg2q(6,6),cg3q(6,6),cg1qb(6,6),cg2qb(6,6),cg3qb(6,6)
     >     ,cqqb(6,6),clo(6,6),cxy(6,6)
      common/Color/cg1g2,cg1g3,cg2g3,cg1q,
     >     cg2q,cg3q,cg1qb,cg2qb,cg3qb,cqqb,clo,cxy
c     run-control:
      double precision ktmin
      integer iwhattodo,iproc,iobs,ifast
      integer logunit,dataunit
      parameter (logunit=11,dataunit=10)
      character*(80) jobfile,outfile,logfile,valuefile
      character*(70) runstring
      common/control/iwhattodo,iobs,iproc,ifast,ktmin,      
     #     jobfile,outfile,logfile,valuefile,runstring
      integer its,pts
      common/vegasparam/its,pts
C
C     runtime errors
C     ierr = -1 precision loss, 
C     ierr = -2 n.n>0 in dggg, 
C     ierr = -3 p.p<0 in boost
      integer ierr
      common/errors/ierr
#ifdef HBOOK
      integer nwpawc
      parameter (nwpawc=1000000)
      real hmem(nwpawc)
      common /pawc/hmem
#endif

