
      subroutine DSgluglu(ih,mh,mz,mw,mt,mg,T1,T2,st,ct,q,mu,
     $     tanb,alpha,DRbar,SM,a1_full,a1,a2)

c     Two-loop SQCD contributions to the Wilson coefficient for gg->Higgs. 
c     Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
c     Based on G. Degrassi and P. Slavich, arXiv:0806.1495
c
c     06/04/2011  OS to DRbar to allow better POWHEG integration.
c     16/08/2010: collected TF as an overall factor (cosmetic change);
c                 allowed OS>1 (compatibility with bottom routine)
c     18/05/2010:  first release
c
c     I/O PARAMETERS:
c     ih = light (1) or heavy (2) CP-even Higgs
c     mh = Higgs mass (for the exact 1-loop result)
c     mz, mw = gauge boson masses 
c     mt = m_top, mg = m_gluino, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = ren. scale,
c     mu = Higgs mixing parameter, tanb = tan(beta), 
c     alpha = mixing angle in the CP-even Higgs sector
c     DRbar = renormalization scheme for 1-loop (1 = DRbar, 0 = On-Shell),
c     SM = include (1) or exclude (0) the purely-SM contributions
c     a1_full = exact result for the one-loop amplitude (top+stop)
c     a1 = one-loop part of the effective ggh vertex in the light-Higgs limit
c     a2 = two-loop part of the effective ggh vertex in the light-Higgs limit

      implicit none

      integer ih,DRbar,SM
      real*8 mh,mz,mw,mt,mg,T1,T2,st,ct,q,mu,tanb,alpha,da2
      real*8 TF,CF,CA,At,Xt,s2t,c2t,sina,cosa,sinb,cosb
      real*8 FF,GG,dFF,dFF_At,dGG,FFo,GGo,FFT,GGT,dFFT,dGGT,FFTo,GGTo
      real*8 DDo,DD,dDD,dL,dR,s2tw
      real*8 H1o,H2o,H1,H2,H1ct,H2ct
      complex*16 Gzero_T1,Gzero_T2,Ghalf_t,Gzero,Ghalf,a1_full,
     $     FFoc,GGoc,FFToc,GGToc,DDoc,H1oc,H2oc,a1,a2

      TF = 1d0/2d0
      CF = 4d0/3d0
      CA = 3d0

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      Xt = (T1-T2)*s2t/2d0/mt    
      At = Xt - mu/tanb           ! note our sign convention for mu

      sina = sin(alpha)
      cosa = cos(alpha)
      sinb = sin(atan(tanb))
      cosb = cos(atan(tanb))

      s2tw = 1-mw**2/mz**2
      dL = 1d0/2d0-2d0/3d0*s2tw
      dR = 2d0/3d0*s2tw

c      one-loop part of the form factor (full result)

      Gzero_T1 = Gzero(4*T1/mh**2)
      Gzero_T2 = Gzero(4*T2/mh**2)

      GGoc = 1/2d0*(Gzero_T1/T1 + GZero_T2/T2)
      FFoc = 1/2d0*(GZero_T1/T1 - Gzero_T2/T2)

      GGToc = GGoc
      FFToc = FFoc

      DDoc = (dL+dR)/2d0*GGToc + c2t*(dL-dR)/2d0*FFToc

      H1oc = (mt*mu*s2t*FFoc+2*mz**2*cosb*sinb*DDoc)/sinb
      H2oc = (mt*At*s2t*FFoc+2*mt**2*GGoc-2*mz**2*sinb**2*DDoc)/sinb

      if(SM.eq.1) then         ! include the SM contribution

         Ghalf_t = Ghalf(4*mt**2/mh**2)

         H2oc = H2oc + Ghalf_t/sinb

      endif

      if(ih.eq.1) then
         a1_full = TF*(-sina*H1oc+cosa*H2oc) ! light Higgs
      else
         a1_full = TF*(cosa*H1oc+sina*H2oc)  ! heavy Higgs
      endif    

c     one-loop part of the form factor (heavy-top limit)

      GGo = 1/2d0*(-1d0/3d0/T1 - 1d0/3d0/T2)
      FFo = 1/2d0*(-1d0/3d0/T1 + 1d0/3d0/T2)

      if(SM.eq.1) then 
         GGo = GGo + (1/2d0)*(-4d0/3d0/mt**2)
      endif

      GGTo = 1/2d0*(-1d0/3d0/T1 - 1d0/3d0/T2)
      FFTo = 1/2d0*(-1d0/3d0/T1 + 1d0/3d0/T2)

      DDo = (dL+dR)/2d0*GGTo + c2t*(dL-dR)/2d0*FFTo

      H1o = (mt*mu*s2t*FFo+2*mz**2*cosb*sinb*DDo)/sinb
      H2o = (mt*At*s2t*FFo+2*mt**2*GGo-2*mz**2*sinb**2*DDo)/sinb
      
      if(ih.eq.1) then
         a1 = TF*(-sina*H1o+cosa*H2o) ! light Higgs
      else
         a1 = TF*(cosa*H1o+sina*H2o)  ! heavy Higgs
      endif

c     two-loop part of the coefficient

      call functions(SM,mt,mg,T1,T2,s2t,c2t,q,CF,CA,GG,FF,GGT,FFT)

      call shifts(mt,mg,T1,T2,s2t,c2t,q,CF,Xt,
     $     dGG,dFF,dGGT,dFFT,dFF_At)

      if(s2t.ne.0d0) then
         
         DD = (dL+dR)/2d0*GGT + c2t*(dL-dR)/2d0*FFT
         
         H1 = (mt*mu*s2t*FF+2*mz**2*cosb*sinb*DD)/sinb
         H2 = (mt*At*s2t*FF+2*mt**2*GG-2*mz**2*sinb**2*DD)/sinb

         if(ih.eq.1) then
            a2 = TF*(-sina*H1+cosa*H2) ! light Higgs
         else
            a2 = TF*(cosa*H1+sina*H2) ! heavy Higgs
         endif

         if(DRbar.ne.1) then
            
            dDD = (dL+dR)/2d0*dGGT + c2t*(dL-dR)/2d0*dFFT 
            
            H1ct = (mt*mu*s2t*dFF+2*mz**2*cosb*sinb*dDD)/sinb
            H2ct = (mt*At*s2t*dFF+2*mt**2*dGG-2*mz**2*sinb**2*dDD)/sinb
            H2ct = H2ct + (mt*s2t*dFF_At)/sinb

            if(ih.eq.1) then
               da2 = TF*(-sina*H1ct+cosa*H2ct)
            else
               da2 = TF*(cosa*H1ct+sina*H2ct)
            endif

            a2 = a2 + da2

         endif

c     FF has poles in s2t=0, when necessary consider the residues:
         
      else
         
         DD = (dL+dR)/2d0*GGT + c2t*(dL-dR)/2d0*FFT
         
         H1 = (mt*mu*FF+2*mz**2*cosb*sinb*DD)/sinb
         H2 = (mt*At*FF+2*mt**2*GG-2*mz**2*sinb**2*DD)/sinb
         
         if(ih.eq.1) then
            a2 = TF*(-sina*H1+cosa*H2) ! light Higgs
         else
            a2 = TF*(cosa*H1+sina*H2) ! heavy Higgs
         endif
         
         if(DRbar.ne.1) then
            
            dDD = (dL+dR)/2d0*dGGT + c2t*(dL-dR)/2d0*dFFT 
            
            H1ct = (mt*mu*dFF+2*mz**2*cosb*sinb*dDD)/sinb
            H2ct = (mt*At*dFF+2*mt**2*dGG-2*mz**2*sinb**2*dDD)/sinb

            if(ih.eq.1) then
               da2 = TF*(-sina*H1ct+cosa*H2ct)
            else
               da2 = TF*(cosa*H1ct+sina*H2ct)
            endif

            a2 = a2 + da2

         endif
         
      endif

      return
      end

*
***********************************************************************
*
      
      subroutine functions(SM,mt,mg,T1,T2,s2t,c2t,q,CF,CA,
     $     GG,FF,GGDT,FFDT)

      implicit none

      integer SM
      real*8 mt,mg,T1,T2,s2t,c2t,q,CF,CA,GG,FF,GGDT,FFDT
      real*8 Gglu,Fglu,GgluDT,FgluDT,Gqua,Fqua,GquaDT,FquaDT,
     $     Gino,Fino,GinoDT,FinoDT
      real*8 DT1,DT2,Dt,Dc2t2,
     $     DT1CF,DtCF_1,Dc2t2CF_1,DT1CA,DtCA_1,Dc2t2CA_1,
     $     DT2CF,DtCF_2,Dc2t2CF_2,DT2CA,DtCA_2,Dc2t2CA_2

c     contributions from diagrams with gluons

      Dt = 1/2d0/mt**2*(CF - 5d0/3d0*CA)
      DT1 = 1/2d0/T1*(-3d0/4d0*CF-1d0/6d0*CA)
      DT2 = 1/2d0/T2*(-3d0/4d0*CF-1d0/6d0*CA)

      Gglu = DT1 + DT2 
      Fglu = DT1 - DT2

      GgluDT = DT1 + DT2
      FgluDT = DT1 - DT2

      if(SM.eq.1) Gglu = Gglu + Dt

c     contributions from diagrams with quartic stop coupling

      DT1 = -CF/24d0*((c2t**2*T1+s2t**2*T2)/T1**2
     $     +s2t**2/T1**2/T2*(T1**2*Log(T1/q**2)-T2**2*Log(T2/q**2)))

      DT2 = -CF/24d0*((c2t**2*T2+s2t**2*T1)/T2**2
     $     +s2t**2/T2**2/T1*(T2**2*Log(T2/q**2)-T1**2*Log(T1/q**2)))

      Dc2t2 = -CF/24d0*((T1-T2)**2/T1/T2
     $     - (T1-T2)/T2*Log(T1/q**2) - (T2-T1)/T1*Log(T2/q**2))
 
      Gqua = DT1 + DT2 
      Fqua = DT1 - DT2 - 4*c2t**2/(T1-T2)*Dc2t2

      GquaDT = DT1 + DT2
      FquaDT = DT1 - DT2 + 4*s2t**2/(T1-T2)*Dc2t2

c     contributions from the diagrams with gluinos

      call gluino(mt,mg,T1,s2t,q,
     $     DT1CF,DtCF_1,Dc2t2CF_1,DT1CA,DtCA_1,Dc2t2CA_1)

      call gluino(mt,mg,T2,-s2t,q,
     $     DT2CF,DtCF_2,Dc2t2CF_2,DT2CA,DtCA_2,Dc2t2CA_2)

      DT1 = CF*DT1CF+CA*DT1CA
      DT2 = CF*DT2CF+CA*DT2CA
      Dt = CF*(DtCF_1+DtCF_2)+CA*(DtCA_1+DtCA_2)
      Dc2t2 = CF*(Dc2t2CF_1+Dc2t2CF_2)+CA*(Dc2t2CA_1+Dc2t2CA_2)
      
      Gino = DT1 + DT2 + Dt
      Fino = DT1 - DT2 - 4*c2t**2/(T1-T2)*Dc2t2

      GinoDT = DT1 + DT2 
      FinoDT = DT1 - DT2 + 4*s2t**2/(T1-T2)*Dc2t2

c     sum all 

      GG = Gglu + Gqua + Gino 
      FF = Fglu + Fqua + Fino
           
      GGDT = GgluDT + GquaDT + GinoDT
      FFDT = FgluDT + FquaDT + FinoDT

      if(s2t.eq.0d0) then

         Dc2t2 = CF*(Dc2t2CF_1-Dc2t2CF_2)+CA*(Dc2t2CA_1-Dc2t2CA_2)
         FF = - 4*c2t**2/(T1-T2)*Dc2t2 ! only the residue

      endif

      return
      end

*     
***********************************************************************
*

      subroutine gluino(mt,mg,T1,s2t,q,
     $     DT1CF,DtCF,Dc2t2CF,DT1CA,DtCA,Dc2t2CA)
      
      implicit none

      real*8 mt,mg,T1,s2t,q,DT1CF,DtCF,Dc2t2CF,DT1CA,DtCA,Dc2t2CA
      real*8 q2,t,g,tq,gq,T1q,del,del2,del3,phigtT1,delt,phi

      q2 = q**2
      t = mt**2
      g = mg**2
      tq = t**2
      gq = g**2
      T1q = T1**2

      del = delt(g,t,T1)
      del2 = del**2
      del3 = del**3

      phigtT1 = phi(g,t,T1)

      DT1CF = 1d0/6d0/T1q/del2*((g+t)*del2+2*g*T1q*(del+10*g*t))
     $     -mg*s2t/6d0/mt/T1q/del2*(t*del2+2*T1q*(del+5*g*t)*(g+t-T1))
     $     -gq/6d0/T1q/del3*Log(g/t)*(
     $     (g-t-4*T1)*del2+2*T1*(-18*t*T1*del
     $     +((3*t-g)*del-30*g*t*T1)*(t-g+T1))
     $     -s2t/mt/mg*((2*T1*(T1+t)+t*(g-t))*del2-2*t*T1*(
     $     -9*g*T1*del+((2*g-9*T1)*del-30*g*t*T1)*(g-t+T1))))
     $     +g/6d0/del3*Log(T1/t)*(del2+12*g*t*del+
     $     (2*T1*del+60*g*t*T1)*(t+g-T1)
     $     -2*s2t/mg/mt*((g+t)*del2
     $     +g*t*(3*g+3*t+20*T1)*del+60*gq*tq*T1))
     $     -1/6d0/T1q*(g+t-s2t*mg*mt)*Log(t/q2)
     $     +gq*t/T1/del3*phigtT1*((g+t+3*T1)*del+20*g*t*T1
     $     -s2t/mg/mt*(del2+2*g*t*del+(3*T1*del+10*g*t*T1)*(g+t-T1)))

      DT1CA = 1/12d0/del2*(2*t*del-(del+20*g*t)*(g-t-T1))
     $     +mt*s2t/3d0/mg/del2*(2*g*del-(del+5*g*t)*(t-g-T1))
     $     +g/12d0/del3*Log(g/T1)*(
     $     del2+2*t*((13*g+9*t+15*T1)*del+120*g*t*T1)
     $     -2*s2t*mt/mg*(18*t*(t-T1)*del
     $     +((11*g-9*t+9*T1)*del+60*g*t*T1)*(g+t-T1)))
     $     -t/12d0/del3*Log(t/T1)*(12*t*(T1-t)*del
     $     +((15*g+13*t+3*T1)*del+120*g*t*T1)*(g+t-T1)
     $     -4*s2t/mg/mt*((3*g+t)*del2
     $     +g*t*((9*g+2*t+21*T1)*del+60*g*t*T1)))
     $     +g*t/2d0/T1/del3*phigtT1*(2*t*T1*del
     $     +((g+t+2*T1)*del+20*g*t*T1)*(t-g+T1)
     $     +s2t/mg/mt*((g-t+T1)*del2+2*t*(2*T1*(t-T1)*del
     $     +((g+5*T1)*del+10*g*t*T1)*(g-t+T1))))

      DtCF = -g/6d0/T1/del2*(4*T1*del+(del-10*g*T1)*(g-t-T1))
     $     -mg*s2t/12d0/mt/T1/del2*(del2+2*(5*g*T1*del
     $     +((2*T1-g)*del+10*g*t*T1)*(g-t+T1)))
     $     +gq/6d0/T1/del3*Log(g/t)*(del2
     $     -2*T1*((2*g+3*t+9*T1)*del+60*g*t*T1)
     $     -s2t/2d0/mg/mt**3*((tq-2*T1q+g*(t+2*T1))*del2 +t*T1*
     $     (20*g*t*del-((25*T1+9*g+11*t)*del+120*g*t*T1)*(g+t-T1))))
     $     +g*T1/3d0/del3*Log(T1/t)*(9*g*del+(del+30*g*t)*(g-t+T1)
     $     +s2t/2d0/mg/mt**3*((g+3*t-T1)*del2
     $     +g*t*(12*t*del+(17*del+60*g*t)*(t-g+T1))))
     $     +1/6d0/T1*(1-s2t/2d0*mg/mt)*Log(t/q2)
     $     -gq/del3*phigtT1*(2*t*del+(del+10*g*t)*(t-g+T1)
     $     +s2t/2d0/mg/mt*(del2+
     $     2*t*(3*g*del+(3*del+10*g*t)*(g-t+T1))))

      DtCA = -1/12d0/del2*((5*g+5*T1-t)*del+40*g*t*T1)
     $     -s2t/12d0/mg/mt/del2*(10*t*(T1-t)*del
     $     +((8*t-5*g-2*T1)*del-20*g*t*T1)*(t+g-T1))
     $     -g/12d0/del3*Log(g/t)*(del2+12*t*T1*del
     $     +((4*g+18*T1)*del+120*g*t*T1)*(t-g+T1)
     $     +s2t/mg/mt*((g+4*t+12*T1)*del2+4*t*(16*T1*(t-T1)*del
     $     +((t+24*T1)*del+30*g*t*T1)*(g-t+T1))))
     $     -T1/12d0/del3*Log(T1/t)*(8*t*(T1-t)*del
     $     +((19*g+9*t+3*T1)*del+120*g*t*T1)*(g+t-T1)
     $     -s2t/mg/mt*((11*g+2*t+2*T1)*del2
     $     +2*g*t*((21*g+3*t+43*T1)*del+120*g*t*T1)))
     $     -g/2d0/del3*phigtT1*(del**2+t*((7*g+7*T1+t)*del+40*g*t*T1)
     $     -s2t/2d0/mg/mt*((g+3*t-T1)*del2+2*t*(2*t*(T1-t)*del
     $     +((3*g+2*t+6*T1)*del+20*g*t*T1)*(g+t-T1))))

      Dc2t2CF = mt*mg/12d0/t1/del*(2*g*T1-del)
     $     +mg**3/12d0/mt/T1/del2*(4*t*T1*del
     $     +((2*T1-t)*del+6*g*t*T1)*(t-g+T1))*Log(g/t)
     $     +mg*T1/6d0/mt/del2*(del+3*g*t)*(g+t-T1)*Log(T1/t)
     $     +mg*mt/12d0/T1*Log(t/q2)
     $     +mg**3*mt/2d0/del2*(del+2*g*t)*phigtT1

      Dc2t2CA = mt*mg/12d0/del*(t+T1-g)
     $     +mg*mt/12d0/del2*((g+6*T1)*del+12*g*t*T1)*Log(g/t)
     $     -mt*T1/12d0/mg/del2*
     $     (3*g*del+(2*del+6*g*t)*(g-t+T1))*Log(T1/t)
     $     -mg*mt/4d0/del2*((del+2*g*t)*(g-t-T1))*phigtT1

      if(s2t.ne.0d0) then

         Dc2t2CF = Dc2t2CF/s2t
         Dc2t2CA = Dc2t2CA/s2t

      endif

      return
      end

*     
***********************************************************************
*

      subroutine shifts(mt,mg,T1,T2,s2t,c2t,qq,CF,Xt,
     $     dGG,dFF,dGGDT,dFFDT,dFF_A)
      
c     shift of the parameters from DRbar to On-Shell scheme
 
      implicit none      
      real*8 mt,mg,T1,T2,s2t,c2t,qq,CF,Xt
     $     ,dGG,dFF,dFF_A,dGGDT,dFFDT
      real*8 msdr,g,t,q,dT1,dT2,dmt,dAt,dth,ds2t,dc2t,myB0

      msdr = -5d0

      g = mg**2
      t = mt**2
      q = qq**2

      dmt =                     ! eq. (B2) of DSZ
     $     CF*mt*(3*Log(t/q) + msdr + .5d0*(2*g/t*(Log(g/q)-1)
     $     -T1/t*(Log(T1/q)-1) - T2/t*(Log(T2/q)-1)
     $     +(g+t-T1 - 2*s2t*mg*mt)/t*myB0(t,g,T1,q)
     $     +(g+t-T2 + 2*s2t*mg*mt)/t*myB0(t,g,T2,q)))
      
      dT1 =                     ! eq. (B3) of DSZ
     $     CF*T1*(3*Log(T1/q) - 7 - c2t**2*(Log(T1/q)-1)
     $     -s2t**2*T2/T1*(Log(T2/q)-1) + 2*(
     $     g/T1*(Log(g/q)-1) + t/T1*(Log(t/q)-1)
     $     +(T1-g-t + 2*s2t*mg*mt)/T1*myB0(T1,t,g,q)))
         
      dT2 =                     ! eq. (B4) of DSZ
     $     CF*T2*(3*Log(T2/q) - 7 - c2t**2*(Log(T2/q)-1)
     $     -s2t**2*T1/T2*(Log(T1/q)-1) + 2*(
     $     g/T2*(Log(g/q)-1) + t/T2*(Log(t/q)-1)
     $     +(T2-g-t - 2*s2t*mg*mt)/T2*myB0(T2,t,g,q)))
         
c$$$c     On-Shell theta-stop: asymmetric definition used in FeynHiggs
c$$$      dth = CF*(4d0*mg*mt*c2t*myB0(T1,t,g,q) +
c$$$     $     c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/(T1-T2)      

c     On-Shell theta-stop: eq. (B6)-(B7) of DSZ 
      dth = CF*(4d0*mg*mt*c2t*(myB0(T1,t,g,q)+myB0(T2,t,g,q)) +
     $     2d0*c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/
     $     2d0/(T1-T2)      

      ds2t = 2*c2t*dth
      dc2t = -2*s2t*dth

      dAt = ((dT1-dT2)/(T1-T2) + ds2t/s2t - dmt/mt)*Xt

c     now the shifts in the functions 
c     taking into account that the counterterms are multiplied by as/4/pi

      if(s2t.ne.0d0) then
         dFF = 1/6d0*(dT1/T1**2-dT2/T2**2
     $        - (dmt/mt+ds2t/s2t)*(1d0/T1-1d0/T2))/4d0
      else
         dFF = 1/6d0*(-ds2t)*(1d0/T1-1d0/T2)/4d0 ! just the residue
      endif

      dGG = 1/6d0*(dT1/T1**2+dT2/T2**2-2*dmt/mt*(1d0/T1+1d0/T2))/4d0

      dFFDT = 1/6d0*(dT1/T1**2-dT2/T2**2-dc2t/c2t*(1d0/T1-1d0/T2))/4d0

      dGGDT = 1/6d0*(dT1/T1**2+dT2/T2**2)/4d0

      dFF_A = -1/6d0*(1d0/T1-1d0/T2)*dAt/4d0

      return
      end

*
***********************************************************************
*

      double precision function myB0(q,m1,m2,mu2) 

c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
      
      double precision q,m1,m2,Omega,mu2

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
            
            myB0 = -(Log(q/mu2) - 2) ! cut the imaginary part (I Pi)
            
         else
            
            myB0 = -( log(q/mu2)-2.d0 + 
     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*log(m1/q) +
     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*log(m2/q) +
     3           2.d0*Omega(m1/q,m2/q))
            
         endif
         
      endif

      return
      end
      
c     function Omega(a,b) contained in myB0
c     optimized on 11/01/2010

      double precision function Omega(a,b)
      double precision a,b,cbig,sqCbig
      Cbig = 0.5d0*(a+b) - 0.25d0*(a-b)*(a-b) -0.25d0
      if(Cbig.gt.0d0) then
         sqCbig = sqrt(Cbig)
         Omega = sqCbig*
     1        (atan((1 + a - b)/(2*sqCbig)) +
     2        atan((1 - a + b)/(2*sqCbig)) )
      elseif(Cbig.lt.0d0) then
         sqCbig = sqrt(-Cbig)
         Omega = 0.5d0*sqCbig*
     1        log((a + b - 1 - 2*sqCbig)/
     2        (a + b - 1 + 2*sqCbig))
      else
         Omega = 0
      endif

      return
      end

*
***********************************************************************
*

      function delt(x,y,z)

      implicit none

      real*8 x,y,z,delt

      delt = x**2 + y**2 + z**2 - 2d0*(x*y + x*z + y*z)

      return
      end



      function Gzero(x)

      double precision x
      double complex Gzero,myff

      Gzero = x*(1 - x*myff(x))
      
      return
      end

      function Ghalf(x)

      double precision x
      double complex Ghalf,myff

      Ghalf = -2*x*(1 + (1- x)*myff(x))
      
      return
      end
      
      function myff(x)

      double precision x,pi,mylog
      double complex myff
      parameter (pi = 3.1415926535897932384626433832795029D0)
 
      if(x.ge.1d0) then
         
         myff = dcmplx(asin(1d0/sqrt(x))**2,0d0)

      else

         mylog = Log((1+sqrt(1-x))/(1-sqrt(1-x)))
         myff = dcmplx(-1d0/4d0*(mylog**2 - pi**2),1d0/2d0*pi*mylog)

      endif

      return
      end
         

      
