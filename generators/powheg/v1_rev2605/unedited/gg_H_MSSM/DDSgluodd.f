      subroutine DDSgluodd(it,massren,SM,mA,mt,mbpole,mbsb,mbrun,mg,
     $     T1,T2,s2t,B1,B2,s2b,mu,tanb,At,Ab,q,
     $     a1tf,a1bf,a1t,a1b,a2t,a2b)

c     Two-loop top-stop-gluino and bottom-sbottom-gluino contributions 
c     to the effective vertex for gg->A in the MSSM. 
c     Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
c     Based on G. Degrassi, S. Di Vita and P. Slavich, arXiv:1107.0914
c
c     20/12/2012: added modes OS=3 and OS=4 for "resummed" mb^pole
c     18/06/2012: first release
c
c     I/O PARAMETERS:
c     it = expansion for the two-loop top-stop-gluino contributions:   
c          (0) limit mA = 0; (1) asymptotic exp. O(Msusy^-2)  
c     OS = renormalization scheme for the quark masses in 1-loop part
c           0 = all DRbar (at the scale q), 1 = all on-shell, 
c           2 = top on-shell, bottom mixed [with mb^MSSM,DRbar(q)],
c           3 = top on-shell, bottom mixed [with mb^pole/(1+Db)], 
c           4 = top on-shell, bottom mixed [with mb^pole/(1+Db)(1-Db/tanb^2)] 
c     SM = include (1) or exclude (0) the 2-loop quark-gluon contributions 
c     mA = pseudoscalar Higgs mass, 
c     mt = top mass (mt^MSSM,DRbar(q) if OS = 0, otherwise mt^pole),
c     mbpole = pole bottom mass,
c     mbsb = effective bottom "mass" used in the sbottom sector,
c     mbrun = MSSM DRbar bottom mass at the scale q if OS = 0,
c             effective mass used only in the Yukawa coupling if OS = 2,3,4
c     mg = m_gluino, 
c     T1 = m_stop1^2, T2 = m_stop2^2, s2t = sin(2*theta_stop)
c     B1 = m_sbot1^2, B2 = m_sbot2^2, s2b = sin(2*theta_sbot)
c     mu = Higgs mixing parameter, tanb = tan(beta), 
c     At = Higgs-stop-stop trilinear, Ab = Higgs-sbot-sbot trilinear,
c     q = renormalization scale (relevant only for OS = 0 or OS = 2),
c     a1tf = one-loop top contribution to the ggA vertex, full
c     a1bf = one-loop bot contribution to the ggA vertex, full
c     a1t = one-loop top contribution to the ggA vertex, limit mA = 0
c     a1b = one-loop bot contribution to the ggA vertex, O(mb^2/mA^2) 
c     a2t = two-loop top/stop contribution to the ggA vertex
c     a2b = two-loop bot/sbot contribution to the ggA vertex

      implicit none

      integer it,OS,SM,massren

      double precision mA,mt,mbpole,mbsb,mbrun,mg,T1,T2,s2t,B1,B2,s2b,
     $     mu,tanb,At,Ab,q
      double complex a1tf,a1bf,a1t,a1b,a2b,a2t,da2b

      double precision TF,CF,Kttg,K1t0,glutop,mb,dmb,dmbsusy,fac
      double complex K1tf,K1bf,K1b,glubot,Kttg2,Kbbg

c     We translate here from the POWHEG internal flags to the Slavich ones
      if(massren.eq.2) then
        OS=3
      elseif(massren.eq.1) then
        OS=0
      else
        write(*,*) 'Error: unrecognized massaren in DDSgluodd.f'
        stop
      endif

      TF = 1d0/2d0
      CF = 4d0/3d0

      if(OS.eq.0) then          ! bottom mass in different schemes
         mb = mbrun         
         fac = 1d0
      elseif(OS.eq.1) then
         mb = mbpole
         fac = 1d0
      elseif(OS.ge.2.and.OS.le.4) then
         mb = mbpole         
         fac = mbrun/mbpole
      endif

c     1-loop quark contributions

      call oneloop(mA,mt,mb,K1tf,K1t0,K1bf,K1b)

      a1tf = TF*K1tf/tanb
      a1bf = TF*K1bf*tanb*fac

      a1t = TF*K1t0/tanb
      a1b = TF*K1b*tanb*fac

c     2-loop quark-gluon contributions
c     (NOTE: the top contribution is in the limit mA=0)

      call getglu(OS,mA,mt,mb,q,glutop,glubot)

      a2t = SM*TF*glutop/tanb
      
c     build the 2-loop top contribution

      if(it.eq.0) then          ! limit mA=0

         call getKttg(mt,mg,T1,T2,mu,tanb,Kttg)

         a2t = a2t + TF*Kttg/tanb

      elseif(it.eq.1) then      ! asymptotic expansion up to O(Msusy^-2)

         call getKttg2(mA,mt,mg,T1,T2,s2t,At,mu,tanb,Kttg2)
         a2t = a2t + TF*Kttg2/tanb            
         
      else
         write(*,*) 'wrong entry for it'
      endif

c     build the 2-loop bottom contribution

      call getKbbg(OS,mA,mb,mbsb,fac,mg,B1,B2,s2b,mu,Ab,tanb,q,Kbbg)

      a2b = SM*TF*glubot*tanb*fac + TF*Kbbg*tanb 
      
      if(OS.eq.2) then          ! shift to the "mixed" scheme
         
         dmb = SM*(3d0/2d0*log(mb/q) - 5d0/4d0)
     $        + dmbsusy(mbsb,mg,B1,B2,s2b,q)
         
         a2b = a2b - TF*tanb*K1b*CF*dmb*fac
         
      elseif(OS.eq.3.or.OS.eq.4) then

         call effshiftpseudo(OS,fac,K1b,tanb,mu,mg,B1,B2,TF,CF,da2b)

         a2b = a2b + da2b
         
      endif

      return
      end

*
***********************************************************************
*

      subroutine oneloop(mA,mt,mb,K1tf,K1t0,K1bf,K1b)

      implicit none

      double precision mA,mt,mb,K1t0,taub,taut,pi
      double complex K1tf,K1bf,Khalf,K1b,loghb

      parameter (pi = 3.1415926535897932384626433832795029D0)

      taut = 4*mt**2/mA**2
      taub = 4*mb**2/mA**2

      K1tf = Khalf(taut)
      K1t0 = -2d0

      K1bf = Khalf(taub)
      loghb = dcmplx(log(4/taub),-pi)      
      K1b = taub/2d0*loghb**2

      return
      end

*
***********************************************************************
*

      subroutine getglu(OS,mA,mt,mb,q,glutop,glubot)

c     two-loop quark-gluon contributions to the ggA vertex
c     eqs.(19)-(26) of 1107.0914
c
c     NOTE: the top-gluon contribution is in the limit mA=0

      implicit none
      
      integer OS
      double precision mA,mt,mb,q,glutop
      double complex glubot
      
      double precision CF,CA,taub,z2,z3,pi,F1t,F2t,F3t,logt,logb

      double complex loghb,F1b,F2b,F3b

      parameter (pi = 3.1415926535897932384626433832795029D0)

      CF = 4d0/3d0
      CA = 3d0
      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0

      taub = 4*mb**2/mA**2

      logt = 2*log(mt/q)
      logb = 2*log(mb/q)
      loghb = dcmplx(log(4/taub),-pi)

c     top contributions 

      F1t = 0d0
      F2t = 0d0
      F3t = -2d0
      
c     bottom contributions

      F1b = -taub*(9d0/5d0*z2**2-z3+(2-z2-4*z3)*loghb
     $     -(1-z2)*loghb**2+loghb**3/4d0+loghb**4/48d0)

      F2b = 3d0*taub/4d0*(2*loghb-loghb**2)

      F3b = taub*(8d0/5d0*z2**2+3*z3-3*z3*loghb+(1+2*z2)/4d0*loghb**2
     $     +loghb**4/48d0)

c     build results depending on the scheme

      if(OS.eq.0) then
         glutop = CF*(F1t+F2t*(logt-1d0/3d0))+CA*F3t
         glubot = CF*(F1b+F2b*(logb-1d0/3d0))+CA*F3b
      else
         glutop = CF*(F1t+4d0/3d0*F2t)+CA*F3t
         glubot = CF*(F1b+4d0/3d0*F2b)+CA*F3b
      endif

      return
      end

*
***********************************************************************
*

      subroutine getKttg(mt,mg,T1,T2,mu,tanb,Kttg)

c     two-loop top-stop-gluino contributions in the limit mA=0
c     eq.(29) of 1107.0914

      implicit none

      double precision mt,mg,T1,T2,mu,tanb,cotb,ft,Kttg

      cotb = 1d0/tanb
      
      Kttg = mu*mt/(T1-T2)*(cotb+tanb)*(ft(mt,mg,T1)-ft(mt,mg,T2)) 

      return
      end

*
***********************************************************************
*

      subroutine getKttg2(mA,mt,mg,T1,T2,s2t,At,mu,tanb,Kttg)

c     two-loop top-stop-gluino contributions 
c     in the asymptotic expansion up to O(1/Msusy^2)
c     eqs.(30)-(34) of 1107.0914

      implicit none

      double precision ma,mt,mg,T1,T2,s2t,At,mu,tanb
      double complex Kttg

      double precision CF,CA,Li2,tau,Yt,mt2,mg2,ma2,x1,x2,logx1,logx2
      double complex K1l,Khalf,B0att,myB0c,
     $     domCF,reststopCF,fact,topCF,topCA,dueCF,dueCA,treCF,treCA

      CF = 4d0/3d0
      CA = 3d0

      ma2 = ma**2
      mt2 = mt**2
      mg2 = mg**2
      
      x1 = T1/mg2
      x2 = T2/mg2
      logx1 = log(x1)
      logx2 = log(x2)

      Yt = At - mu*tanb

      tau = 4*mt2/ma2
      
      K1l = Khalf(tau)

      B0att = myB0c(ma2,mt2,mt2,mt2)

      domCF = -CF/2d0*K1l*mg/mt*(s2t/2d0-mt*Yt/(T1-T2))*
     $     (x1/(1-x1)*logx1-x2/(1-x2)*logx2)

      fact = 2*Log(mg2/mt2) - 3 - 3d0/2d0*K1l + 2*B0att

      topCF = CF/4d0/(1-x1)**3*((1-x1**2+2*x1*logx1)*fact
     $     -8*x1*Li2(1-x1)-2*x1*(3+x1)*logx1)
     $      - CF/4d0/(1-x2)**3*((1-x2**2+2*x2*logx2)*fact
     $     -8*x2*Li2(1-x2)-2*x2*(3+x2)*logx2)

      fact = Log(mt2/mg2) + 1 + K1l/2d0 - B0att

      topCA = CA/2d0/(1-x1)**2*((1-x1+x1*logx1)*fact
     $     +2*x1*Li2(1-x1)+x1*(1+x1)*logx1)
     $      - CA/2d0/(1-x2)**2*((1-x2+x2*logx2)*fact
     $     +2*x2*Li2(1-x2)+x2*(1+x2)*logx2)     

      reststopCF = CF/(x1-x2)**2*Yt/mg*(
     $     (x1**2*(1-2*x2)/2d0/(1-x1)/(1-x2)
     $     +x1/2/(1-x1)**2*(x1**2-2*x2+x1*x2)*logx1)-
     $     (x2**2*(1-2*x1)/2d0/(1-x2)/(1-x1)
     $     +x2/2/(1-x2)**2*(x2**2-2*x1+x2*x1)*logx2))
      
      dueCF = CF/4d0/(1-x1)**3*(2*(1-x1**2+2*x1*logx1)*log(mg2/mt2)
     $     -8*x1*Li2(1-x1)+(1-x1**2)*(1+K1l/2d0)
     $     -2*x1*(2+x1-K1l/2d0)*logx1)
     $      - CF/4d0/(1-x2)**3*(2*(1-x2**2+2*x2*logx2)*log(mg2/mt2)
     $     -8*x2*Li2(1-x2)+(1-x2**2)*(1+K1l/2d0)
     $     -2*x2*(2+x2-K1l/2d0)*logx2)

      dueCA = CA/2d0/(1-x1)**2*((1-x1+x1*logx1)*log(mt2/mg2)
     $     +2*x1*Li2(1-x1)+x1*(1+x1)*logx1)
     $      - CA/2d0/(1-x2)**2*((1-x2+x2*logx2)*log(mt2/mg2)
     $     +2*x2*Li2(1-x2)+x2*(1+x2)*logx2)

      fact = 2 + K1l - B0att

      treCF = CF/6d0/(1-x1)**4*(-2-3*x1+6*x1**2-x1**3-6*x1*logx1)*fact
     $      + CF/6d0/(1-x2)**4*(-2-3*x2+6*x2**2-x2**3-6*x2*logx2)*fact

      fact = 2 + K1l - 2*B0att
      
      treCA = CA/8d0/(1-x1)**3*(1-x1**2+2*x1*logx1)*fact
     $      + CA/8d0/(1-x2)**3*(1-x2**2+2*x2*logx2)*fact

      Kttg = domCF - mt/mg*s2t*(topCF+topCA+(1+K1l/2d0)*reststopCF)
     $     +2*mt2*Yt/mg/(T1-T2)*(dueCF+dueCA)+mt2/mg2*(treCF+treCA)
     $     -K1l/2d0*mA2/(T1-T2)*reststopCF

      return
      end

*
***********************************************************************
*

      subroutine getKbbg(OS,mA,mb,mbsb,fac,mg,B1,B2,s2b,mu,Ab,tanb,q,
     $     Kbbg)

c     two-loop bottom-sbottom-gluino contributions 
c     in the asymptotic expansion up to O(mb/Msusy) and O(mb^2/mA^2)
c     eqs.(35) and (41) of 1107.0914

      implicit none

      integer OS
      double precision mA,mb,mbsb,fac,mg,B1,B2,s2b,mu,Ab,tanb,q
      double complex Kbbg

      double precision cotb,g,x1,x2,CF,Yb,taub,pi,dmbsusy
      double complex loghb,K1l,F2b,dom,rest,Rb

      parameter (pi = 3.1415926535897932384626433832795029D0)

      CF = 4d0/3d0

      cotb = 1d0/tanb
      Yb = Ab - mu*cotb

      g = mg**2
      x1 = B1/g
      x2 = B2/g

      taub = 4*mb**2/mA**2
      loghb = dcmplx(log(4/taub),-pi)

      K1l = taub/2d0*loghb**2

      F2b = 3d0*taub/4d0*(2*loghb-loghb**2)

c     NOTE: we rescale by "fac" only the terms due to the A-bot-bot coupling

      dom = -CF/2d0*K1l*mg/mbsb*(fac*s2b/2d0-mbsb*Yb/(B1-B2))
     $     *(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2))
     
      rest= -mbsb/mg*s2b*(Rb(Yb,mA,mg,x1,x2,fac)-Rb(Yb,mA,mg,x2,x1,fac))

      Kbbg = dom + rest

      if(OS.eq.0) then          ! shift to DRbar 
         Kbbg = Kbbg + fac*4d0/3d0*CF*F2b*dmbsusy(mbsb,mg,B1,B2,s2b,q)
      endif

      return
      end

*     
***********************************************************************
*

      subroutine effshiftpseudo(OS,fac,K1b,tanb,mu,mg,B1,B2,TF,CF,dH)

c     additional shift for OS = 3 or 4, i.e. when the one-loop OS bottom
c     contribution is rescaled as in the effective Lagrangian approach

      implicit none

      integer OS
      double precision fac,tanb,mu,mg,B1,B2,TF,CF
      double complex K1b,H1l,dH

      double precision x1,x2,Deltab

      x1 = B1/mg**2
      x2 = B2/mg**2

      Deltab = CF/2d0*mg*mu*tanb/(B1-B2)
     $     *(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2))
      
      H1l = TF*tanb*K1b
      
      if(OS.eq.3) then
         dH = H1l*Deltab*fac
      elseif(OS.eq.4) then
         dH = H1l*Deltab*(fac+1d0/tanb**2)
      else
         write(*,*) 'WRONG OS in subroutine effshiftpseudo'
      endif
      
      return
      end

*
***********************************************************************
*
      
      double complex function Rb(Yb,mA,mg,x1,x2,fac)

c     eq.(36) of 1107.0914

      implicit none

      double precision Yb,mA,mg,x1,x2,fac,CF,CA,logx1,Li2,pi
      double complex loghg
      parameter (pi = 3.1415926535897932384626433832795029D0)

      CF = 4d0/3d0
      CA = 3d0

      logx1 = log(x1)      
      loghg = dcmplx(log(mA**2/mg**2),-pi)
      
      Rb = fac*CF/4d0/(1-x1)**3*(
     $     (1-x1**2+2*x1*logx1)*(1-2*loghg)
     $     -8*x1*Li2(1-x1)-2*x1*(3+x1)*logx1)
     $     +fac*CA/2d0/(1-x1)**2*(
     $     (1-x1+x1*logx1)*(loghg-1)
     $     +2*x1*Li2(1-x1)+x1*(1+x1)*logx1)
     $     +CF/(x1-x2)**2*Yb/mg*(
     $     x1**2*(1-2*x2)/2d0/(1-x1)/(1-x2)
     $     +x1/2d0/(1-x1)**2*(x1**2-2*x2+x1*x2)*logx1)

      return
      end

*
***********************************************************************
*

      double precision function ft(mt,mg,T1)

c     eq.(28) of 1107.0914

      implicit none

      double precision mt,mg,T1
      double precision CF,Nc,t,g,phi,phigtT1,del,delt

      CF = 4d0/3d0
      Nc = 3d0

      t = mt**2
      g = mg**2

      del = delt(g,t,T1)
      phigtT1 = phi(g,t,T1)

      ft = CF*mg/mt*((g-t+T1)*t*Log(t/g)+(g+t-T1)*T1*Log(T1/g)
     $     + 2*g*t*phigtT1)/del
     $     + Nc*mt/mg*((T1-t-g)*T1*Log(t/g)+(t-T1-g)*T1*Log(T1/g)
     $     +(t+T1-g)*g*phigtT1)/del
            
      return
      end

*     
***********************************************************************
*

      double precision function dmbsusy(mb,mg,B1,B2,s2b,q)

c     SUSY contribution to bottom self-energy, eqs.(39)-(40) of 1107.0914
 
      implicit none

      double precision mb,mg,B1,B2,s2b,q

      double precision x1,x2,f1,f2

      x1 = B1/mg**2
      x2 = B2/mg**2

      f1 = (x1-3)/4d0/(1-x1)+x1*(x1-2)/2d0/(1-x1)**2*log(x1)
      f2 = (x2-3)/4d0/(1-x2)+x2*(x2-2)/2d0/(1-x2)**2*log(x2)

      dmbsusy = -1/4d0*(2*log(mg/q)+f1+f2
     $     +mg/mb*s2b*(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2)))

      return
      end

*     
***********************************************************************
*

      double precision function dmtsusypseudo(mt,mg,T1,T2,s2t,q)

c     SUSY contribution to top self-energy, from eq.(B2) of hep-ph/0105096

      implicit none

      double precision mt,mg,T1,T2,s2t,q,q2,t,g,
     $     logG,logT1,logT2,B0tgT1,B0tgT2
      double complex myB0c

      t = mt**2
      g = mg**2
      q2 = q**2

      logG = log(g/q2)
      logT1 = log(T1/q2)
      logT2 = log(T2/q2)

      B0tgT1 = dreal(myB0c(t,g,T1,q2))
      B0tgT2 = dreal(myB0c(t,g,T2,q2))

      dmtsusypseudo = g/t*(logG-1) - 0.5d0*(
     $     T1/t*(LogT1-1)-(g+t-T1-2*s2t*mg*mt)/t*B0tgT1+
     $     T2/t*(LogT2-1)-(g+t-T2+2*s2t*mg*mt)/t*B0tgT2)

      return
      end

*
***********************************************************************
*

      function Khalf(x)

      double precision x
      double complex Khalf,myff

      Khalf = -2*x*myff(x)
      
      return
      end


*
***********************************************************************
*

      double complex function myB0c(q,m1,m2,mu2) 

c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
      
      implicit none
      double precision q,m1,m2,Omega,mu2,myB0,ImB0,sm1,sm2,pi
      parameter (pi = 3.1415926535897932384626433832795029D0)

      if(q.eq.0d0) then

         if(m1.eq.0d0.and.m2.ne.0d0) then
            myB0 = 1d0-Log(m2/mu2)
         elseif(m1.ne.0d0.and.m2.eq.0d0) then
            myB0 = 1d0-Log(m1/mu2)
         elseif(abs(m1-m2).le.1d-8) then
            myB0 = -Log(m1/mu2)
         else
            myB0 = 1d0 - Log(m2/mu2) + m1/(m1-m2)*Log(m2/m1)
         endif
         
      else

         if(m1.eq.0d0.and.m2.ne.0d0) then
            
            if(m2.ne.q) then
               myB0 = -(Log(m2/mu2)-2-(m2/q-1d0)*Log(abs(1d0-q/m2))) 
            else 
               myB0 = -(Log(m2/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.ne.0d0) then
            
            if(m1.ne.q) then
               myB0 = -(Log(m1/mu2)-2-(m1/q-1d0)*Log(abs(1d0-q/m1))) 
            else
               myB0 = -(Log(m1/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.eq.0d0) then
            
            myB0 = -(Log(q/mu2) - 2)
            
         else
            
            myB0 = -( log(q/mu2)-2.d0 + 
     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*log(m1/q) +
     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*log(m2/q) +
     3           2.d0*Omega(m1/q,m2/q))
            
         endif
         
      endif

      sm1 = sqrt(m1)
      sm2 = sqrt(m2)

      if(q.le.(sm1+sm2)**2) then
         ImB0 = 0d0
      else
         ImB0 = pi*sqrt(q-(sm1+sm2)**2)*sqrt(q-(sm1-sm2)**2)/q
      endif

      myB0c = dcmplx(myB0,ImB0)

      return
      end
      
c     function Omega(a,b) contained in myB0
c     optimized on 11/01/2010
!
!      double precision function Omega(a,b)
!      double precision a,b,cbig,sqCbig
!      Cbig = 0.5d0*(a+b) - 0.25d0*(a-b)*(a-b) -0.25d0
!      if(Cbig.gt.0d0) then
!         sqCbig = sqrt(Cbig)
!         Omega = sqCbig*
!     1        (atan((1 + a - b)/(2*sqCbig)) +
!     2        atan((1 - a + b)/(2*sqCbig)) )
!      elseif(Cbig.lt.0d0) then
!         sqCbig = sqrt(-Cbig)
!         Omega = 0.5d0*sqCbig*
!     1        log((a + b - 1 - 2*sqCbig)/
!     2        (a + b - 1 + 2*sqCbig))
!      else
!         Omega = 0
!      endif
!
!      return
!      end

*
***********************************************************************
*
      function phi(x,y,z)

c     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      implicit none
      double precision x,y,z,phi,pphi,myphi
      
      if(x.le.z.and.y.le.z) then
         pphi = myphi(x,y,z)
      elseif(z.le.x.and.y.le.x) then
         pphi = z/x*myphi(z,y,x)
      elseif(z.le.y.and.x.le.y) then
         pphi = z/y*myphi(z,x,y)
      endif

      phi = pphi
      
      end
      
      function myphi(x,y,z)
      
      implicit none

      double precision x,y,z,myphi
      double precision u,v
      double precision Pi,Li2
      double complex clam,cxp,cxm,CLI2,ccphi

      parameter (pi = 3.1415926535897932384626433832795029D0)

c     auxiliary variables

      u = x/z
      v = y/z
      
      if(u.le.1d-8) then
         
         if(v.ne.1d0) then
            myphi = (log(u)*log(v)+2d0*Li2(1d0-v))/(1d0-v)
         else
            myphi = 2d0-log(u)
         endif

      elseif(v.le.1d-8) then

         if(u.ne.1d0) then
            myphi = (log(v)*log(u)+2d0*Li2(1d0-u))/(1d0-u)
         else
            myphi = 2d0-log(v)
         endif

      else
         
         if((1d0-u-v)**2.ge.4d0*u*v) then         
            clam = DCMPLX(sqrt((1d0-u-v)**2 - 4d0*u*v),0d0)
         else
            clam = DCMPLX(0d0,sqrt(4d0*u*v - (1d0-u-v)**2))
         endif
         
         cxp = (1d0+(u-v)-clam)/2d0
         cxm = (1d0-(u-v)-clam)/2d0
         
c     phi function from eq. (A4)
            
         ccphi = (2d0*log(cxp)*log(cxm) - log(u)*log(v) - 
     &        2d0*(CLI2(cxp) + CLI2(cxm)) + Pi**2/3d0)/clam
         myphi = DBLE(ccphi)
                     
      endif
      
      return
      end

*
***********************************************************************
*

      function Li2(x)

      implicit none

      double complex CLI2,z
      double precision x,Li2

      z = DCMPLX(x,0d0)
      Li2 = DBLE(CLI2(z))

      return
      end

