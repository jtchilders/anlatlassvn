c-------------------------------------------------------------------
      subroutine usrcut(lnot,wusr)
c-------------------------------------------------------------------
c PRIMARY CUTS ALREADY APPLIED TO PHASE-SPACE GENERATION:
c     ptjmin < pt(jet) < ptjmax for all light jets
c     -etajmax < eta(jet) < etajmax for all light jets
c     delta R(jj) > drjmin for all (light jet, light jet) pairs
c     delta R(bj) > drjmin for all (heavy quark, light jet) pairs
c     ptbmin < pt(b) < ptbmax for all heavy quarks
c     -etabmax < eta(b) < etabmax for all heavy quarks 
c     delta R(bb) > drjmin for all (heavy quark, heavy quark) pairs
c
c USE THIS ROUTINE TO ENFORCE LEPTON CUTS, OR OTHER CUTS    
      implicit none
      include 'alpgen.inc'
      include 'hjet.inc'
      double precision wusr
      integer lnot,i
c initialize output parameters
      lnot=0
      wusr=1d0
c     
 5    return
c 10   lnot= 1
      end

c-------------------------------------------------------------------
      subroutine alshis
c-------------------------------------------------------------------
      include 'alpgen.inc'
      include 'hjet.inc'
      real ptbin,ptmax,xmbin,xmmax
      ptbin=4e0
      ptmax=400e0
      xmbin=50e0
      xmmax=5000e0
*
c      do i=1,njets
c        call mbook(i,'pt_jet',ptbin,0e0,ptmax)
c      enddo
      end

c-------------------------------------------------------------------
      subroutine alfhis
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'hjet.inc'
      integer i
      real  xnorm
      character *1 jet(9)
      data jet/'1','2','3','4','5','6','7','8','9'/
c
      double precision fcount
c
      open(unit=99,file=topfile,err=101,status='unknown')
      if(imode.le.1) then
         xnorm=sngl(avgwgt/totwgt)
      elseif(imode.eq.2) then
         xnorm=1e0/real(unwev)
      else
         write(*,*) 'imode type not allowed, stop'
         stop
      endif
c      do i=1,200
c         if(i.ne.61) call mopera(i,'F',i,i,xnorm,1.)
c         call mfinal(i)
c      enddo 
c      do i=1,njets
c         call mtop(i,99,'pt'//jet(i),' ','LOG')
c      enddo
 100  close(99)
 101  return
      end


      subroutine monitor(n,mon_fname)
c     This routine is called by default every 100K events.
c     The user can use it to get regular updates on the run
c     while this is progressing. Textual output can be written to file
c     fname, where partial cross-sections and and generation
c     efficiencies have already been printed by default
      implicit none
      include 'alpgen.inc'
      include 'hjet.inc'
      integer n
      character *50 mon_fname
c
      if(evgen) then
         if(mod(n,1000000).eq.0) then
c     save histograms' contents
            call msave
c     print out histograms
            call alfhis
c     restore original contents, to proceed with analysis
            call mrestore
         endif 
      endif 
      end
c-------------------------------------------------------------------
      subroutine aleana(jproc,wgt)
c     analyse event, fill histograms
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'hjet.inc' 
      real*8 wgt
      real rwgt,mjets
      real*8 ptc(maxpar),etac(maxpar),mjj
      common/ptanal/mjj,ptc,etac
      real*8 ptlep,ptj
      integer i,iw,jproc
c
      rwgt=real(wgt)
      if(rwgt.lt.0e0) then
         write(*,*) 'negative wgt=',wgt
         return
      elseif (rwgt.eq.0e0) then
         return
      endif 
      end
*
* INPUT:  QJS_0 ,QJS_X, QJS_Y, QJS_Z, QI_0, QI_1, QI_2, QI_3
* OUTPUT: QJ_0, QJ_X, QJ_Y, QJ_Z, DJ2
*
      SUBROUTINE BOOSTPIC(QJS,QI,QJ)
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QJS(0:3),QI(0:3),QJ(0:3)
      DIMENSION QJT(3),AN(3),QJSL(3)
*
      QIL= SQRT(QI(1)*QI(1)+QI(2)*QI(2)+QI(3)*QI(3))
*
      ZERO = 1.D-13
*
      IF(QI(0).GT.ZERO.AND.QIL.GT.ZERO) THEN
         BI= QIL/QI(0)
         IF(BI.GE.1.D0) THEN
            GML = 1.D13
         ELSE
            GML= 1.D0/SQRT(1.D0-BI*BI)
         ENDIF
      ELSE
         BI = 0.D0
         GML = 1.D0
      ENDIF
*
      IF(ABS(QIL).LT.ZERO) THEN
         AN(1)= 0.D0
         AN(2)= 0.D0
         AN(3)= 0.D0
      ELSE
         AN(1)= -QI(1)/QIL
         AN(2)= -QI(2)/QIL
         AN(3)= -QI(3)/QIL
      ENDIF
*
      QJSLM= QJS(1)*AN(1)+QJS(2)*AN(2)+QJS(3)*AN(3)
      QJSL(1)= QJSLM*AN(1)
      QJSL(2)= QJSLM*AN(2)
      QJSL(3)= QJSLM*AN(3)
      QJT(1)= QJS(1)-QJSL(1)
      QJT(2)= QJS(2)-QJSL(2)
      QJT(3)= QJS(3)-QJSL(3)
*
      QJ(0)= GML*(QJS(0)-BI*QJSLM)
      QJLM= GML*(-BI*QJS(0)+QJSLM)
*
      QJ(1)= QJLM*AN(1)+QJT(1)
      QJ(2)= QJLM*AN(2)+QJT(2)
      QJ(3)= QJLM*AN(3)+QJT(3)
*
      RETURN
      END

      SUBROUTINE PDECAY(P,PM,Q1,Q1M,Q2,Q2M,WT)
C--- Decays a particle with momentum P into two particles with momenta
C--- Q(1) and Q(2). WT is the phase space density beta_cm
C--- The decay is spherically symmetric in the decay C_of_M frame
      implicit none
      double precision pi,twopi
      PARAMETER(PI=3.14159,TWOPI=2.*PI)
      double precision q2e,qp,ctheta,stheta,phi,qplab,qplon,qptr,pmod
     $     ,ptr,wt,bet,gam
      double precision P(4),Q1(4),Q2(4),V(3),U(3),pm,q1m,q2m
      double precision x1,x2
      integer i
      CALL RANDA(X1)
      CALL RANDA(X2)
      Q2E=(PM**2-Q1M**2+Q2M**2)/(2.*PM)
      QP=SQRT(MAX(Q2E**2-Q2M**2,0.d0))
      CTHETA=2.*real(x1)-1.
      STHETA=SQRT(1.-CTHETA**2)
      PHI=TWOPI*real(x2)
      QPLON=QP*CTHETA
      QPTR=QP*STHETA
      PMOD=SQRT(P(1)**2+P(2)**2+P(3)**2)
      PTR=SQRT(P(2)**2+P(3)**2)                              

C--- if the decaying particle moves along the X axis:
      IF(PTR.LT.1.E-4) THEN
        V(1)=0.
        V(2)=1.
        V(3)=0.
        U(1)=0.
        U(2)=0.
        U(3)=1.
      ELSE
C--- 
        V(1)=0.
        V(2)=P(3)/PTR
        V(3)=-P(2)/PTR
        U(1)=PTR/PMOD
        U(2)=-P(1)*P(2)/PTR/PMOD
        U(3)=-P(1)*P(3)/PTR/PMOD
      ENDIF
      GAM=P(4)/PM
      BET=PMOD/P(4)
      QPLAB=GAM*(QPLON+BET*Q2E)
      DO I=1,3
      Q2(I)=QPLAB*P(I)/PMOD+QPTR*(V(I)*SIN(PHI)+U(I)*COS(PHI))
      Q1(I)=P(I)-Q2(I)
      END DO
      Q2(4)=GAM*(Q2E+BET*QPLON)
      Q1(4)=P(4)-Q2(4)
      WT=2.*QP/PM
      END            


