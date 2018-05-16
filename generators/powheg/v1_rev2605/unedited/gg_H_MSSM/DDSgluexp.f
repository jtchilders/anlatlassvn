
      subroutine DDSgluexp(ih,mh,mz,mw,mt,mg,T1,T2,st,ct,q,mu,
     $     tanb,alpha,DRBAR,SM,a1_full,a1,a2)

c     Two-loop SUSY contributions to the Wilson coefficient for gg->Higgs 
c     in an asymptotic expansion valid for mh,mt << Msusy
c     Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
c     Based on G. Degrassi, S. Di Vita and P. Slavich, arXiv:1204.1016
c
c     16/05/2012:  first release
c
c     I/O PARAMETERS:
c     ih = light (1) or heavy (2) CP-even Higgs
c     mh = Higgs mass
c     mz, mw = gauge boson masses 
c     mt = m_top, mg = m_gluino, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = ren. scale,
c     mu = Higgs mixing parameter, tanb = tan(beta), 
c     alpha = mixing angle in the CP-even Higgs sector
c     DRBAR = renormalization scheme for 1-loop (0 = On Shell, 1 = DR-bar),
c     SM = include (1) or exclude (0) the top contribution in one-loop part
c     a1_full = exact result for the one-loop amplitude (top+stop or stop only)
c     a1 = one-loop part of the ggh vertex (top+stop or stop only)
c     a2 = two-loop part of the ggh vertex (NOTE: SUSY contributions only!!!)

      implicit none

      integer ih,DRBAR,SM
      double precision mh,mz,mw,mt,mg,T1,T2,st,ct,q,mu,tanb,alpha
      double complex a1_full,a1,a2

      double precision TF,CF,CA,At,Xt,s2t,c2t,sina,cosa,sinb,cosb,
     $     dL,dR,s2tw

      double precision FFo,GGo,FFTo,GGTo,DDo

      double complex Gzero_T1,Gzero_T2,Ghalf_t,Gzero,Ghalf,
     $     FFoc,GGoc,FFToc,GGToc,DDoc,H1oc,H2oc,H1o,H2o,
     $     FF,GG,dFF,dFF_At,dGG,FFT,GGT,dFFT,dGGT,DD,dDD,
     $     H1,H2,H1ct,H2ct,da2

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

c     one-loop part of the form factor (heavy-squark limit)

      GGo = 1/2d0*(-1d0/3d0/T1 - 1d0/3d0/T2)
      FFo = 1/2d0*(-1d0/3d0/T1 + 1d0/3d0/T2)

      GGTo = GGo
      FFTo = FFo

      DDo = (dL+dR)/2d0*GGTo + c2t*(dL-dR)/2d0*FFTo

      H1o = (mt*mu*s2t*FFo+2*mz**2*cosb*sinb*DDo)/sinb
      H2o = (mt*At*s2t*FFo+2*mt**2*GGo-2*mz**2*sinb**2*DDo)/sinb
      
      if(SM.eq.1) then         ! include the SM contribution

         Ghalf_t = Ghalf(4*mt**2/mh**2)
         H2o = H2o + Ghalf_t/sinb

      endif

      if(ih.eq.1) then
         a1 = TF*(-sina*H1o+cosa*H2o) ! light Higgs
      else
         a1 = TF*(cosa*H1o+sina*H2o)  ! heavy Higgs
      endif

c     two-loop part of the coefficient

      call DDS_functions(mt,mh,mg,T1,T2,s2t,c2t,q,CF,CA,GG,FF,GGT,FFT)

      call DDS_shifts(mt,mh,mg,T1,T2,s2t,c2t,q,CF,Xt,
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

         if ((DRBAR.eq.0).or.(DRBAR.eq.2)) then
            
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
         
         if ((DRBAR.eq.0).or.(DRBAR.eq.2)) then
            
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
      
      subroutine DDS_functions(mt,mh,mg,T1,T2,s2t,c2t,q,CF,CA,
     $     GG,FF,GGDT,FFDT)

      implicit none

      double precision mt,mh,mg,T1,T2,s2t,c2t,q,CF,CA
      double complex GG,FF,GGDT,FFDT

      double precision q2

      double complex Gglu,Fglu,GgluDT,FgluDT,Gqua,Fqua,GquaDT,FquaDT,
     $     Gino,Fino,GinoDT,FinoDT

      double complex Dt,DT1,DT2,Dc2t2,
     $     Dt_1,DT1g(2),Dc2t2_1(2),Dt_2,DT2g(2),Dc2t2_2(2)

c     contributions from diagrams with squarks and gluons

      DT1 = -1/2d0/T1*(3d0/4d0*CF+CA/6d0)
      DT2 = -1/2d0/T2*(3d0/4d0*CF+CA/6d0)

      Gglu = DT1 + DT2 
      Fglu = DT1 - DT2

      GgluDT = DT1 + DT2
      FgluDT = DT1 - DT2

c     contributions from diagrams with quartic stop coupling

      q2 = q**2

      DT1 = -CF/24d0*((c2t**2*T1+s2t**2*T2)/T1**2
     $     +s2t**2/T1**2/T2*(T1**2*Log(T1/q2)-T2**2*Log(T2/q2)))

      DT2 = -CF/24d0*((c2t**2*T2+s2t**2*T1)/T2**2
     $     +s2t**2/T2**2/T1*(T2**2*Log(T2/q2)-T1**2*Log(T1/q2)))

      Dc2t2 = -CF/24d0*((T1-T2)**2/T1/T2
     $     - (T1-T2)/T2*Log(T1/q2) - (T2-T1)/T1*Log(T2/q2))
 
      Gqua = DT1 + DT2 
      Fqua = DT1 - DT2 - 4*c2t**2/(T1-T2)*Dc2t2

      GquaDT = DT1 + DT2
      FquaDT = DT1 - DT2 + 4*s2t**2/(T1-T2)*Dc2t2

c     contributions from the diagrams with gluinos

      call DDS_gluino(mt,mh,mg,T1,T2, s2t,q,CF,CA,DT1g,Dt_1,Dc2t2_1)
      call DDS_gluino(mt,mh,mg,T2,T1,-s2t,q,CF,CA,DT2g,Dt_2,Dc2t2_2)

      Dt =  Dt_1 + Dt_2 

      Dc2t2 = Dc2t2_1(1) + Dc2t2_2(1)

      Gino = DT1g(1) + DT2g(1) + Dt

      GinoDT = DT1g(1) + DT2g(1)

      FinoDT = DT1g(1) - DT2g(1) + 4*s2t**2/(T1-T2)*Dc2t2

c     add the O(mt/M) terms for the function F

      Dc2t2 = Dc2t2 + Dc2t2_1(2) + Dc2t2_2(2)

      Fino = (DT1g(1)+DT1g(2))-(DT2g(1)+DT2g(2))-4*c2t**2/(T1-T2)*Dc2t2

c     sum all 

      GG = Gglu + Gqua + Gino 
      FF = Fglu + Fqua + Fino
           
      GGDT = GgluDT + GquaDT + GinoDT
      FFDT = FgluDT + FquaDT + FinoDT

      if(s2t.eq.0d0) then

         Dc2t2 = (Dc2t2_1(1)+Dc2t2_1(2)) - (Dc2t2_2(1)+Dc2t2_2(2))
 
         FF = -4*c2t**2/(T1-T2)*Dc2t2 ! only the residue in 1/s2t

      endif

      return
      end

*     
***********************************************************************
*

      subroutine DDS_gluino(mt,mh,mg,T1,T2,s2t,q,CF,CA,
     $     YT1,Yt,Yc2t2)
      
      implicit none

      double precision mt,mh,mg,T1,T2,s2t,q,CF,CA
      double complex Yt,YT1(2),Yc2t2(2)

      double precision q2,h,t,g,p_Li2,tau,dmtovmt,dmtsusy,
     $     x,y,x2,x3,umx,umx2,umx3,umx4,logx,Lix,logtg,loggq
      double complex Ghalf,Khalf,Fhalf,myBB,G12,K12,BB,F12b,R1,R2

c     misc definitions

      q2 = q**2
      h = mh**2
      t = mt**2
      g = mg**2
      logtg = log(t/g)
      loggq = log(g/q2)

      x = T1/g
      y = T2/g
      x2 = x*x
      x3 = x2*x
      umx = 1-x
      umx2 = umx*umx
      umx3 = umx2*umx
      umx4 = umx2*umx2
      logx = log(x)
      Lix = p_Li2(1-x)
      
      tau = 4*t/h
      G12 = Ghalf(tau)
      K12 = Khalf(tau)
      F12b = Fhalf(tau)
      BB = myBB(tau)

      dmtovmt = dmtsusy(3,mt,mg,T1,s2t,q,CF)/mt

c     here are the various contributions

      R1 = CA/6/umx2*     
     $     (3*(1-x+x*logx)*(logtg - BB - K12/2d0 + 2)
     $     +6*x*Lix+2*x+2*x*(1+x)*logx-2)
     $     -CF/6d0/x/umx3*
     $     (3*(x-x3+2*x2*logx)*(logtg - BB - G12/4d0 - K12/2d0 + 2)
     $     +umx3*loggq+12*x2*Lix+5*x3-5*x2+x-1+2*(x3+2*x2)*logx)

      R2 = - CA/12d0/umx3*(
     $     3*(1-x2+2*x*logx)*(2*logtg - BB - K12/2d0 + 2)
     $     +24*x*Lix+1-x2+2*x*(3*x+10)*logx)
     $     + CF/18/x/umx4*
     $     (3*x*(umx*(5*x-x2+2)+6*x*logx)
     $     *(2*logtg - BB - G12/2d0 - K12/2d0 + 2)+6*umx4*loggq
     $     +72*x2*Lix-x*umx2*(11*x-26)-6*umx+6*x2*(2*x+9)*logx)

      Yt = (4d0/3d0*F12b*dmtovmt - CF/4d0*G12*mg/mt*s2t*x/umx*logx
     $     + s2t*mt/mg*R1 + t/g*R2)/2d0/t

      YT1(1) = (CF/4d0*s2t/mg/mt*G12-(2*CF+CA)/12/g)*(1d0/umx+logx/umx2)
     $     +CF/24d0/g/x2/umx3*
     $     (4*umx3*(1-loggq)-3*x2*G12*(umx*(3-x)+2*logx))
      
      YT1(2) = CF*s2t*mt/6d0/mg**3/x2/umx4*
     $     (3*x2*(umx*(x+5)+2*(2*x+1)*logx)*(G12/4d0-logtg)
     $     +umx4*loggq-12*x2*(2*x+1)*Lix
     $     -umx*(14*x2-3*x+1)-2*x2*(x2+18*x+5)*logx)
     $     + CA*s2t*mt/6d0/mg**3/umx3*
     $     (3*(2-2*x+(x+1)*logx)*(1+logtg)
     $     +6*(1+x)*Lix+2*x*umx+2*(6*x+1)*logx)
     $     + CF*s2t*h*G12/48d0/mt/mg**3/x/umx4*
     $     (umx*(x2-5*x-2)-6*x*logx)

      Yc2t2(1) = -CF*mg/8d0/mt*G12*x/umx*logx

      Yc2t2(2) = CF*mt/12d0/mg/x/umx3*
     $     (3*x*(x2-2*x*logx-1)*(G12/4-logtg)+umx3*loggq
     $     +12*x2*Lix-(1-2*x)*umx2+2*x2*(x+5)*logx)
     $     +CA*mt*x/12d0/mg/umx2*
     $     ((x-1-logx)*(3*logtg+1)-6*Lix-2*(x+2)*logx)
     $     +CF*h*G12/32d0/mg/mt/umx2*((umx*(x+y-2*x*y)/(1-y)/(x-y)
     $     +2*x*(x2+x*y-2*y)*logx/(x-y)**2))

      if(s2t.ne.0d0) then
         Yc2t2(1) = Yc2t2(1)/s2t
         Yc2t2(2) = Yc2t2(2)/s2t
      endif

      return
      end

*     
***********************************************************************
*

      double precision function dmtsusy(l,mt,mg,T1,s2t,q,CF)

      implicit none

      integer l
      double precision mt,mg,T1,s2t,q,CF

      integer i
      double precision dmtl(0:3),q2,t,g,x,x2,umx,umx2,umx3,logx
      
      q2 = q**2
      t = mt**2
      g = mg**2

      x = T1/g
      x2 = x*x
      umx = 1-x
      umx2 = umx*umx
      umx3 = umx2*umx
      logx = log(x)

      dmtl(0) = s2t*mg/mt*x/umx*logx
      
      dmtl(1) = log(g/q2)/2d0 + (x-3)/4d0/umx + x*(x-2)/2d0/umx2*logx

      dmtl(2) = s2t*mt/2d0/mg/umx3*(1-x2+2*x*logx)

      dmtl(3) = t/6d0/g/umx3*(x2-5*x-2-6*x/umx*logx)

      dmtsusy = 0d0

      do i = 0,l
         dmtsusy = dmtsusy - CF/4d0*mt*dmtl(i)
      enddo

      return
      end

*     
***********************************************************************
*

      subroutine DDS_shifts(mt,mh,mg,T1,T2,s2t,c2t,q,CF,Xt,
     $     dGG,dFF,dGGDT,dFFDT,dFF_A)
      
c     shift of the parameters from DRbar to On-Shell scheme
 
      implicit none      
      double precision mt,mh,mg,T1,T2,s2t,c2t,q,CF,Xt
      double complex dGG,dFF,dGGDT,dFFDT,dFF_A

      double precision msdr,g,t,h,q2,tau,
     $     dT10,dT11,dT20,dT21,dT1,dT2,
     $     dAt,dAtm1,dth,dth0,dth1,ds2t,dc2t,
     $     dmtsm,dmtsusy,dmt(0:3),
     $     logT1q,logT2q,logtq,loggq,logabsx1,logabsx2
      double complex Fhalf
      integer i
      
      msdr = -5d0

      g = mg**2
      t = mt**2
      h = mh**2
      q2 = q**2
      logtq = log(t/q2)
      loggq = log(g/q2)
      logT1q = log(T1/q2)
      logT2q = log(T2/q2)
      logabsx1 = log(abs(1-T1/g))
      logabsx2 = log(abs(1-T2/g))
      tau = 4*t/h

c     compute the shifts (expanded in mt, units of as/pi) 

      dmtsm = CF/4d0*mt*(3*logtq+msdr)

      do i = 0,3

         dmt(i) = dmtsusy(i,mt,mg,T1,s2t,q,CF) +
     $        dmtsusy(i,mt,mg,T2,-s2t,q,CF)

         if(i.gt.0) dmt(i) = dmt(i) + dmtsm

      enddo

      dT10 = CF/4d0*T1*(3*logT1q-3-c2t**2*(logT1q-1)
     $     -s2t**2*T2/T1*(logT2q-1)-6*g/T1-2*(1-2*g/T1)*loggq
     $     -2*(1-g/T1)**2*logabsx1)

      dT20 = CF/4d0*T2*(3*logT2q-3-c2t**2*(logT2q-1)
     $     -s2t**2*T1/T2*(logT1q-1)-6*g/T2-2*(1-2*g/T2)*loggq
     $     -2*(1-g/T2)**2*logabsx2)

      dT11 = -CF*s2t*mt*mg*(loggq+(1-g/T1)*logabsx1-2)
      dT21 =  CF*s2t*mt*mg*(loggq+(1-g/T2)*logabsx2-2)
      
      dth0 = CF/4d0*c2t*s2t/(T1-T2)*(T1*(logT1q-1)-T2*(logT2q-1))
      
      dth1 =-CF/2d0*c2t*mt*mg/(T1-T2)*(loggq+(1-g/T1)*logabsx1-2)
     $      +CF/2d0*c2t*mt*mg/(T2-T1)*(loggq+(1-g/T2)*logabsx2-2)

c     first compute dF and dAt (they need one more power of mt)

      dT1 = dT10 + dT11         
      dT2 = dT20 + dT21         
      dth = dth0 + dth1

      ds2t = 2*c2t*dth

      if(s2t.ne.0d0) then
         dFF = 1/6d0*(dT1/T1**2 - dT2/T2**2
     $        - (dmt(2)/mt + ds2t/s2t)*(1d0/T1 - 1d0/T2)
     $        -2d0/15d0*h*dmt(0)/mt*(1d0/T1**2 - 1d0/T2**2))
      else
         dFF = 1/6d0*(-ds2t)*(1d0/T1 - 1d0/T2) ! just the residue
      endif

      dAt = ((dT1-dT2)/(T1-T2) + ds2t/s2t - dmt(2)/mt)*Xt
      
      dAtm1 = -dmt(0)/mt*Xt
      
      dFF_A = -1/6d0*(dAt*(1d0/T1-1d0/T2)
     $     + 2d0/15d0*h*(1d0/T1**2-1d0/T2**2)*dAtm1)

c     now the shifts in the other functions (lower order in mt)

      dT1 = dT10
      dT2 = dT20         
      dth = dth0

      ds2t = 2*c2t*dth
      dc2t = -2*s2t*dth

      dGG = 1d0/6d0*(dT1/T1**2+dT2/T2**2-2*dmt(1)/mt*(1d0/T1+1d0/T2))
     $     -2d0/3d0*Fhalf(tau)*(dmt(3)-dmtsm)/mt**3

      dFFDT = 1d0/6d0*(dT1/T1**2-dT2/T2**2-dc2t/c2t*(1d0/T1-1d0/T2))

      dGGDT = 1d0/6d0*(dT1/T1**2+dT2/T2**2)

      return
      end

*
***********************************************************************
*

      double complex function myBB(tau)

      implicit none
      
      double precision tau,pi,mysqrt
      parameter (pi = 3.1415926535897932384626433832795029D0)

      if(tau.gt.1) then

         mysqrt = sqrt(tau-1)

         myBB = dcmplx(2-2*mysqrt*atan(1d0/mysqrt),0d0)
         
      else
         
         mysqrt = sqrt(1-tau)

         myBB = dcmplx(2-mysqrt*log((1+mysqrt)/(1-mysqrt)),mysqrt*pi)

      endif

      return
      end
      
*
***********************************************************************
*

      double complex function Fhalf(tau)

      implicit none
      
      double precision tau
      double complex Ghalf,Khalf,myBB

      Fhalf = -3d0/2d0*(2*Ghalf(tau) + tau*myBB(tau) - Khalf(tau))

      return
      end

c$$$*
c$$$***********************************************************************
c$$$*
c$$$
c$$$      function Gzero(x)
c$$$
c$$$      double precision x
c$$$      double complex Gzero,myff
c$$$
c$$$      Gzero = x*(1 - x*myff(x))
c$$$      
c$$$      return
c$$$      end
c$$$
c$$$      function Ghalf(x)
c$$$
c$$$      double precision x
c$$$      double complex Ghalf,myff
c$$$
c$$$      Ghalf = -2*x*(1 + (1- x)*myff(x))
c$$$      
c$$$      return
c$$$      end

      
c$$$      function myff(x)
c$$$
c$$$      double precision x,pi,mylog
c$$$      double complex myff
c$$$      parameter (pi = 3.1415926535897932384626433832795029D0)
c$$$ 
c$$$      if(x.ge.1d0) then
c$$$         
c$$$         myff = dcmplx(asin(1d0/sqrt(x))**2,0d0)
c$$$
c$$$      else
c$$$
c$$$         mylog = Log((1+sqrt(1-x))/(1-sqrt(1-x)))
c$$$         myff = dcmplx(-1d0/4d0*(mylog**2 - pi**2),1d0/2d0*pi*mylog)
c$$$
c$$$      endif
c$$$
c$$$      return
c$$$      end
c$$$         
c$$$*
c$$$***********************************************************************
c$$$*
c$$$
c$$$      function p_Li2(x)
c$$$
c$$$      implicit none
c$$$
c$$$      double complex CLI2,z
c$$$      double precision x,p_Li2
c$$$
c$$$      z = DCMPLX(x,0d0)
c$$$      p_Li2 = DBLE(CLI2(z))
c$$$
c$$$      return
c$$$      end
c$$$
c$$$*
c$$$***********************************************************************
c$$$*
c$$$
c$$$      DOUBLE COMPLEX FUNCTION CLI2(Z)
c$$$
c$$$c     just call the Dilog routine
c$$$      
c$$$      DOUBLE COMPLEX Z,Dilog
c$$$
c$$$      CLI2 = Dilog(Z)
c$$$
c$$$      return
c$$$      end
c$$$
c$$$*
c$$$**********************************************************************
c$$$*
c$$$* Dilog.F
c$$$* complex dilogarithm
c$$$* this file is part of FeynHiggs
c$$$* last modified 20 Oct 05 th
c$$$
c$$$
c$$$      double complex function Dilog(z)
c$$$      implicit none
c$$$      double complex z
c$$$      
c$$$      double complex Dilogsum
c$$$      external Dilogsum
c$$$      
c$$$      double precision absz, abs1z
c$$$      double complex t, mlogz
c$$$      
c$$$      double precision pi, zeta2
c$$$      parameter (pi = 3.1415926535897932384626433832795029D0)
c$$$      parameter (zeta2 = pi*pi/6D0)
c$$$      
c$$$      absz = abs(z)
c$$$      if( absz .lt. 1D-20 ) then
c$$$         Dilog = -log(1 - z)
c$$$         return
c$$$      endif
c$$$      
c$$$      abs1z = abs(1 - z)
c$$$	if( abs1z .lt. 1D-20 ) then
c$$$           Dilog = zeta2
c$$$           return
c$$$	endif
c$$$        
c$$$	if( DBLE(z) .gt. .5D0 ) then
c$$$           mlogz = -log(z)
c$$$           t = zeta2 + mlogz*log(1 - z)
c$$$           if( abs1z .gt. 1 ) then
c$$$              Dilog = Dilogsum(log(1 - 1/z)) + zeta2 +
c$$$     $             .5D0*log(z - 1)**2 + t
c$$$           else
c$$$	    Dilog = -Dilogsum(mlogz) + t
c$$$         endif
c$$$      else
c$$$         if( absz .gt. 1 ) then
c$$$	    Dilog = -Dilogsum(-log(1 - 1/z)) - zeta2 - .5D0*log(-z)**2
c$$$         else
c$$$	    Dilog = Dilogsum(-log(1 - z))
c$$$         endif
c$$$      endif
c$$$      end
c$$$      
c$$$
c$$$************************************************************************
c$$$
c$$$      double complex function Dilogsum(w)
c$$$      implicit none
c$$$      double complex w
c$$$      
c$$$      double complex u, t
c$$$      integer k
c$$$      
c$$$      double precision b2, b4, b6, b8, b10, b12, b14
c$$$      double precision b16, b18, b20, b22
c$$$      parameter (b2 = 1/6D0)
c$$$      parameter (b4 = -1/30D0)
c$$$      parameter (b6 = 1/42D0)
c$$$      parameter (b8 = -1/30D0)
c$$$      parameter (b10 = 5/66D0)
c$$$      parameter (b12 = -691/2730D0)
c$$$      parameter (b14 = 7/6D0)
c$$$      parameter (b16 = -3617/510D0)
c$$$      parameter (b18 = 43867/798D0)
c$$$      parameter (b20 = -174611/330D0)
c$$$      parameter (b22 = 854513/138D0)
c$$$      
c$$$      double precision bernoulliB(11)
c$$$      data bernoulliB /b2, b4, b6, b8, b10, b12, b14,
c$$$     &     b16, b18, b20, b22/
c$$$      
c$$$      Dilogsum = w*(1 - .25D0*w)
c$$$      if( abs(w) .lt. 1D-10 ) return
c$$$      
c$$$      u = w
c$$$      do k = 1, 11
c$$$         u = u*w**2/DBLE(2*k*(2*k + 1))
c$$$         t = u*bernoulliB(k)
c$$$         Dilogsum = Dilogsum + t
c$$$         if( abs(t) .lt. 1D-16*abs(Dilogsum) ) return
c$$$      enddo
c$$$
c$$$      end
