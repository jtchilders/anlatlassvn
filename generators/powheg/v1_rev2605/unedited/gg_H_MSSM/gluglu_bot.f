
      subroutine gluglu_bot(ih,mh,mz,mw,mbpole,mbsb,mbrun,mbmh,mg,B1,B2,
     $     sb,cb,q,mu,tanb,alpha,DRbar,flg_mb_mh,SM,a1_full,a1,a2)

c     Two-loop SQCD bottom/sbottom contributions to the Wilson coefficient 
c     for gg->Higgs. 
c     Routine written by P. Slavich (e-mail: slavich@lpthe.jussieu.fr).
c     Based on G. Degrassi and P. Slavich, arXiv:1007.3465
c
c     06/04/2011: Switched from OS to DRbar to allow better POWHEG integration
c     16/08/2010: collected TF as an overall factor (cosmetic change)
c     28/07/2010: first release
c
c     I/O PARAMETERS:
c     ih = light (1) or heavy (2) CP-even Higgs,
c     mh = Higgs mass, mz, mw = gauge boson masses, 
c     mbpole = pole bottom mass, 
c     mbsb = effective bottom "mass" used in the sbottom sector,
c     mbrun = MSSM DRbar bottom mass at q as given in input,
c     mbmh  = MSSM DRbar bottom mass at q=mh (for DRBAR=1 and flg_mb_m=1h)
c     mg = m_gluino, B1 = m_sbot1^2, B2 = m_sbot2^2, 
c     sb = sin(theta_sbot), cb = cos(theta_sbot), 
c     q = renormalization scale (relevant only for DR = 1 DR = 2),
c     mu = Higgs mixing parameter, tanb = tan(beta), 
c     alpha = mixing angle in the CP-even Higgs sector
c     DRbar = renormalization scheme (1 = DRbar, 0 = On-Shell, 2 = Mixed),
c     flg_mb_mh = 1 use mb(mh) for the bottom mass
c     SM = include (1) or exclude (0) the purely-SM contributions
c     a1_full = exact result for the one-loop amplitude (bot+sbot)
c     a1 = one-loop part of the effective ggh vertex in the light-Higgs limit
c     a2 = two-loop part of the effective ggh vertex in the light-Higgs limit

      implicit none

      integer ih,DRbar,SM,flg_mb_mh
      real*8 mh,mz,mw,mbpole,mbsb,mbrun,mbmh,mg,B1,B2,sb,cb,q,mu,
     $     tanb,alpha
      real*8 TF,CF,CA,Ab,Xb,s2b,c2b,sina,cosa,sinb,cosb,taub,mb,fac
     $     ,mbsm,qb

      real*8 FFo,GGo,FFTo,GGTo,DDo,dL,dR,s2tw,pi
      complex*16 FF,GG,FFT,GGT,DD,Gzero_1,Gzero_2,Ghalf_b,Gzero,Ghalf,
     $     FFoc,GGoc,FFToc,GGToc,DDoc,H1o,H2o,H1oc,H2oc,H1,H2,H1ct,H2ct,
     $     a1_full,a1,a2,da2,dH1ct

      parameter (pi = 3.1415926535897932384626433832795029D0)

      TF = 1d0/2d0
      CF = 4d0/3d0
      CA = 3d0

      s2b = 2d0*cb*sb
      c2b = cb**2 - sb**2

      Xb = (B1-B2)*s2b/2d0/mbsb    
      Ab = Xb - mu*tanb           ! note our sign convention for mu

      sina = sin(alpha)
      cosa = cos(alpha)
      sinb = sin(atan(tanb))
      cosb = cos(atan(tanb))

      s2tw = 1-mw**2/mz**2

      dL = -1d0/2d0+1d0/3d0*s2tw
      dR = -1d0/3d0*s2tw
c     DRBAR
      if(DRbar.eq.1) then
         mb = mbrun
      else
         mb = mbpole
      endif
c     Mixed
      if((DRbar.eq.2).or.(DRbar.eq.4).or.(DRbar.eq.5)) then          ! use running mass in the hbb vertex
         fac = mbrun/mbpole
      else
         fac = 1d0
      endif
c     flg_mb_mh  = 1
      if (flg_mb_mh.eq.1) then
         mbsm = mbmh
         qb = mh
      else
         mbsm = mb
         qb = q
      endif

      taub = 4*mbsm**2/mh**2

c      one-loop part of the form factor (full result)

      Gzero_1 = Gzero(4*B1/mh**2)
      Gzero_2 = Gzero(4*B2/mh**2)
         
      GGoc = 1/2d0*(Gzero_1/B1 + GZero_2/B2)
      FFoc = 1/2d0*(GZero_1/B1 - Gzero_2/B2)
         
      GGToc = GGoc
      FFToc = FFoc
         
      DDoc = (dL+dR)/2d0*GGToc + c2b*(dL-dR)/2d0*FFToc

      H1oc = (mbsb*Ab*s2b*FFoc+2*mbsb**2*GGoc+2*mz**2*cosb**2*DDoc)/cosb
      H2oc = (mbsb*mu*s2b*FFoc-2*mz**2*cosb*sinb*DDoc)/cosb

      if(SM.eq.1) then         ! include the SM contributions

         H1oc = H1oc + fac*Ghalf(taub)/cosb

      endif

      if(ih.eq.1) then
         a1_full = TF*(-sina*H1oc+cosa*H2oc) ! light Higgs
      else
         a1_full = TF*(cosa*H1oc+sina*H2oc)  ! heavy Higgs
      endif    

c     one-loop part of the form factor (light-Higgs, light-bottom limit)

      GGo = 1/2d0*(-1d0/3d0/B1 - 1d0/3d0/B2)
      FFo = 1/2d0*(-1d0/3d0/B1 + 1d0/3d0/B2)

      GGTo = GGo
      FFTo = FFo

      DDo = (dL+dR)/2d0*GGTo + c2b*(dL-dR)/2d0*FFTo

      H1o = (mbsb*Ab*s2b*FFo+2*mbsb**2*GGo+2*mz**2*cosb**2*DDo)/cosb
      H2o = (mbsb*mu*s2b*FFo-2*mz**2*cosb*sinb*DDo)/cosb

      Ghalf_b = -2*taub+taub/2d0*dcmplx(log(4/taub),-pi)**2

      if(SM.eq.1) then 

            H1o = H1o + fac*Ghalf_b/cosb

      endif
      
      if(ih.eq.1) then
         a1 = TF*(-sina*H1o+cosa*H2o) ! light Higgs
      else
         a1 = TF*(cosa*H1o+sina*H2o)  ! heavy Higgs
      endif

c     two-loop part of the coefficient

      call b_functions(SM,mb,mbsb,mbsm,fac,mh,mg,B1,B2,s2b,c2b,q,qb,
     $     CF,CA,GG,FF,GGT,FFT)

      DD = (dL+dR)/2d0*GGT + c2b*(dL-dR)/2d0*FFT
         
      if(s2b.ne.0d0) then
         
         H1 = (mbsb*Ab*s2b*FF+2*mb**2*GG+2*mz**2*cosb**2*DD)/cosb
         H2 = (mbsb*mu*s2b*FF-2*mz**2*cosb*sinb*DD)/cosb
         
      else                      ! consider the residues for s2b=0
         
         H1 = (mbsb*Ab*FF+2*mb**2*GG+2*mz**2*cosb**2*DD)/cosb
         H2 = (mbsb*mu*FF-2*mz**2*cosb*sinb*DD)/cosb

      endif

      if(ih.eq.1) then
         a2 = TF*(-sina*H1+cosa*H2)  ! light Higgs
      else
         a2 = TF*(cosa*H1+sina*H2)   ! heavy Higgs
      endif

c     shift to On-Shell parameters if we are onshell or in the mixed

      if ((DRbar.ne.1).and.(DRbar.ne.3)) then

         call shiftsbot(SM,mb,mbsb,fac,mh,mg,B1,B2,s2b,c2b,q,CF,
     $        Ab,mu,sinb,cosb,mz,dL,dR,H1ct,H2ct)

         if(DRbar.eq.2) then
            
            call extrashift(SM,q,fac,Ghalf_b,mb,mbsb,cosb,
     $           mg,B1,B2,s2b,CF,dH1ct)
            
            H1ct = H1ct + dH1ct
            
         endif
         
         if(ih.eq.1) then
            da2 = TF*(-sina*H1ct+cosa*H2ct)
         else
            da2 = TF*(cosa*H1ct+sina*H2ct)
         endif
         
         a2 = a2 + da2
         
      endif
      
      if(DRBAR.eq.4.or.DRBAR.eq.5) then
         call effshift(ih,DRBAR,fac,Ghalf_b,alpha,tanb,mu,mg,B1,B2,TF,
     $        CF,da2)

         a2 = a2 + da2
      endif
      
      return
      end

*
***********************************************************************
*
      
      subroutine b_functions(SM,mb,mbsb,mbsm,fac,mh,mg,B1,B2,s2b,c2b,
     $     q,qb,CF,CA,GG,FF,GGDT,FFDT)

      implicit none

      integer SM

      real*8 mb,mbsb,mbsm,fac,mh,mg,B1,B2,s2b,c2b,q,qb,CF,CA,pi

      complex*16 GG,FF,GGDT,FFDT

      complex*16 Gglu,Fglu,GgluDT,FgluDT,Gqua,Fqua,GquaDT,FquaDT,
     $     Gino,Fino,GinoDT,FinoDT,Fa,Fb,Ga,tau,logt,z2,z3

      complex*16 DB1,DB2,Db,Dc2b2,
     $     DB1CF,DbCF_1,Dc2b2CF_1,DB1CA,DbCA_1,Dc2b2CA_1,
     $     DB2CF,DbCF_2,Dc2b2CF_2,DB2CA,DbCA_2,Dc2b2CA_2

      parameter (pi = 3.1415926535897932384626433832795029D0)
      
c     contributions from diagrams with gluons

      z2 = 1.6449340668482264365d0
      z3 = 1.2020569031595942854d0

      tau = 4*mbsm**2/mh**2
      logt = dcmplx(log(mh**2/mbsm**2),-pi)

      Fa = -tau*(9+9/5d0*z2**2-z3-(1+z2+4*z3)*logt-(1-z2)*logt**2
     $     +logt**3/4d0+logt**4/48d0)

      Fb = 3*tau*(1+logt/2d0-logt**2/4d0)

      Ga = -tau*(3-8/5d0*z2**2-3*z3+3*z3*logt-(1+2*z2)*logt**2/4d0
     $     -logt**4/48d0)

      Db = 1/2d0/mb**2*(CF*(Fa+Fb*(Log(mbsm**2/qb**2)-1/3d0))+CA*Ga)

      DB1 = 1/2d0/B1*(-3d0/4d0*CF-1d0/6d0*CA)
      DB2 = 1/2d0/B2*(-3d0/4d0*CF-1d0/6d0*CA)

      Gglu = 0d0                ! the squarks enter at O(mb^2) 
      Fglu = DB1 - DB2

      GgluDT = DB1 + DB2
      FgluDT = DB1 - DB2

      if(SM.eq.1) Gglu = fac*Db

c     contributions from diagrams with quartic sbottom coupling

      DB1 = -CF/24d0*((c2b**2*B1+s2b**2*B2)/B1**2
     $     +s2b**2/B1**2/B2*(B1**2*Log(B1/q**2)-B2**2*Log(B2/q**2)))

      DB2 = -CF/24d0*((c2b**2*B2+s2b**2*B1)/B2**2
     $     +s2b**2/B2**2/B1*(B2**2*Log(B2/q**2)-B1**2*Log(B1/q**2)))

      Dc2b2 = -CF/24d0*((B1-B2)**2/B1/B2
     $     - (B1-B2)/B2*Log(B1/q**2) - (B2-B1)/B1*Log(B2/q**2))
 
      Gqua = 0d0                ! the squarks enter at O(mb^2)
      Fqua = DB1 - DB2 - 4*c2b**2/(B1-B2)*Dc2b2

      GquaDT = DB1 + DB2
      FquaDT = DB1 - DB2 + 4*s2b**2/(B1-B2)*Dc2b2

c     contributions from the diagrams with gluinos

      call b_gluino(mb,mbsb,mbsm,mh,mg,B1,s2b,q,qb,
     $     DB1CF,DbCF_1,Dc2b2CF_1,DB1CA,DbCA_1,Dc2b2CA_1)

      call b_gluino(mb,mbsb,mbsm,mh,mg,B2,-s2b,q,qb,
     $     DB2CF,DbCF_2,Dc2b2CF_2,DB2CA,DbCA_2,Dc2b2CA_2)

      DB1 = CF*DB1CF+CA*DB1CA
      DB2 = CF*DB2CF+CA*DB2CA
      Db = CF*(DbCF_1+DbCF_2)+CA*(DbCA_1+DbCA_2)
      Dc2b2 = CF*(Dc2b2CF_1+Dc2b2CF_2)+CA*(Dc2b2CA_1+Dc2b2CA_2)
      
      Gino = fac*Db 
      Fino = DB1 - DB2 - 4*c2b**2/(B1-B2)*Dc2b2

      GinoDT = DB1 + DB2 
      FinoDT = DB1 - DB2 + 4*s2b**2/(B1-B2)*Dc2b2

c     sum all 

      GG = Gglu + Gqua + Gino 
      FF = Fglu + Fqua + Fino
           
      GGDT = GgluDT + GquaDT + GinoDT
      FFDT = FgluDT + FquaDT + FinoDT

      if(s2b.eq.0d0) then

         Dc2b2 = CF*(Dc2b2CF_1-Dc2b2CF_2)+CA*(Dc2b2CA_1-Dc2b2CA_2)
         FF = - 4*c2b**2/(B1-B2)*Dc2b2 ! only the residue

      endif

      return
      end

*     
***********************************************************************
*

      subroutine b_gluino(mb,mbsb,mbsm,mh,mg,B1,s2b,q,qb,
     $     DB1CF,DbCF,Dc2b2CF,DB1CA,DbCA,Dc2b2CA)
      
      implicit none

      real*8 mb,mbsb,mbsm,mh,mg,B1,s2b,q,qb
      real*8 g,x,x2,x3,logx,loggq,taub,f1,p_Li2,LI2x,pi,dmbsusy

      complex*16 G12,F12,loghb,loghg,
     $     DB1CF,DbCF,Dc2b2CF,DB1CA,DbCA,Dc2b2CA

      parameter (pi = 3.1415926535897932384626433832795029D0)

      g = mg**2
      x = B1/g
      x2 = x**2
      x3 = x**3

      logx = log(x)
      Li2x = p_Li2(1-1/x)

      loggq = log(g/q**2)
      loghg = dcmplx(log(mh**2/g),-pi)
c      loghb = dcmplx(log(mh**2/mb**2),-pi)
      loghb = dcmplx(log(mh**2/mbsm**2),-pi)

      f1 = (x-3)/4d0/(1-x)+x*(x-2)/2d0/(1-x)**2*logx

c      taub = 4*mb**2/mh**2
      taub = 4*mbsm**2/mh**2

      G12 = -2*taub+taub/2d0*loghb**2
      F12 = 3*taub*(1+loghb/2d0-loghb**2/4d0)

c      dmbsusy = -(f1+mg/mbsb*s2b*x/(1-x)*logx+log(g/q**2)/2d0)/4d0
      dmbsusy = -(f1+mg/mbsb*s2b*x/(1-x)*logx+log(g/qb**2)/2d0)/4d0

      DB1CF = s2b/4d0/mbsb/mg*G12*(1/(1-x)+1/(1-x)**2*logx)
     $     -1/6d0/g*(1/(1-x)-1/x2+logx/(1-x)**2+loggq/x2)

      DB1CA = -1/12d0/g*(1/(1-x)+1/(1-x)**2*logx)

      DbCF = 4/3d0*F12*dmbsusy-G12/4d0*mg/mbsb*s2b*x/(1-x)*logx
     $     -mbsb/mg*s2b/6/x/(1-x)**3*
     $     ((1-x)**3*loggq+2*(x3+2*x2)*logx-3*(x3-x-2*x2*logx)*loghg
     $     +5*x3-5*x2+x-1-12*x2*Li2x-6*x2*logx**2)

      DbCA = mbsb/mg*s2b/6/(1-x)**2*(2*x*(1+x)*logx+2*x-2-6*x*Li2x
     $     -3*x*logx**2+3*(1-x+x*logx)*loghg)

      Dc2b2CF = -mg/mbsb/8d0*G12*x/(1-x)*logx

      Dc2b2CA = 0d0

      DbCF = DbCF/2d0/mb**2
      DbCA = DbCA/2d0/mb**2

      if(s2b.ne.0d0) then

         Dc2b2CF = Dc2b2CF/s2b

      endif

      return
      end

*     
***********************************************************************
*

      subroutine shiftsbot(SM,mb,mbsb,fac,mh,mg,B1,B2,s2b,c2b,qq,CF,
     $     Ab,mu,sinb,cosb,mz,dL,dR,H1ct,H2ct)
      
c     shift of the parameters from DRbar to On-Shell scheme
c     renormalization scheme as in BDSZ2

      implicit none
      integer SM
      real*8 mb,mbsb,fac,mh,mg,B1,B2,s2b,c2b,qq,CF,Ab,mu,
     $     sinb,cosb,mz,dL,dR
      complex*16 H1ct,H2ct

      real*8 msdr,g,q,dB1,dB2,dAb,dth,ds2b,dc2b,
     $     dmb_0,dmb_g,dhb,dmb,taub,FFol,pi,x1,x2,f1,f2
      complex*16 dGG,dFF,dGGDT,dFFDT,dDD,GGol_b,log4tau

      parameter (pi = 3.1415926535897932384626433832795029D0)

      msdr = -5d0

      g = mg**2
      q = qq**2

      taub = 4*mb**2/mh**2
      log4tau = dcmplx(Log(4/taub),-pi)

      dmb_g = CF*(3*Log(mb**2/q) + msdr)
       
      x1 = B1/g
      x2 = B2/g

      f1 = (x1-3)/4d0/(1-x1)+x1*(x1-2)/2d0/(1-x1)**2*log(x1)
      f2 = (x2-3)/4d0/(1-x2)+x2*(x2-2)/2d0/(1-x2)**2*log(x2)

      dmb_0 = -CF*(Log(g/q)+f1+f2
     $     +mg/mbsb*s2b*(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2)))
           
      dB1 = CF*B1*(3*Log(B1/q)-3-c2b**2*(Log(B1/q)-1)
     $     -s2b**2*B2/B1*(Log(B2/q)-1)-6*g/B1
     $     -2*(1-2*g/B1)*Log(g/q)-2*(1-g/B1)**2*Log(Abs(1-B1/g)))

      dB2 = CF*B2*(3*Log(B2/q)-3-c2b**2*(Log(B2/q)-1)
     $     -s2b**2*B1/B2*(Log(B1/q)-1)-6*g/B2
     $     -2*(1-2*g/B2)*Log(g/q)-2*(1-g/B2)**2*Log(Abs(1-B2/g)))

      dth = CF*c2b*s2b*(-1+(B1*Log(B1/q)-B2*Log(B2/q))/(B1-B2))      

      ds2b = 2*c2b*dth
      dc2b = -2*s2b*dth

      dhb = CF*(-4+2*Log(g/q)
     $     +2*B1/(B1-B2)*(2*Log(B1/g)-(1-g/B1)**2*Log(Abs(1-B1/g)))
     $     +2*B2/(B2-B1)*(2*Log(B2/g)-(1-g/B2)**2*Log(Abs(1-B2/g))))

      dAb = fac*2*CF*mg*(4-2*Log(g/q)
     $     -(1-g/B1)*Log(Abs(1-B1/g))-(1-g/B2)*Log(Abs(1-B2/g)))

c     shift the masses entering the functions 
c     taking into account that the counterterms are multiplied by as/4/pi

      if(SM.eq.1) then
         dmb = dmb_0 + dmb_g
      else         
         dmb = dmb_0   ! drop the SM contribution
      endif
      
      dGG = fac*(-taub/mb**2*log4tau*dmb)/4d0
      
      dFF = 1/6d0*(dB1/B1**2-dB2/B2**2)/4d0
      
      dFFDT = 1/6d0*(dB1/B1**2-dB2/B2**2)/4d0

      dGGDT = 1/6d0*(dB1/B1**2+dB2/B2**2)/4d0

      dDD = (dL+dR)/2d0*dGGDT + c2b*(dL-dR)/2d0*dFFDT

      H1ct = (mbsb*Ab*s2b*dFF+2*mb**2*dGG+2*mz**2*cosb**2*dDD)/cosb
      H2ct = (mbsb*mu*s2b*dFF-2*mz**2*cosb*sinb*dDD)/cosb

c     shift the coefficients of the functions 
c     taking into account that the counterterms are multiplied by as/4/pi

c     shift G

      GGol_b  = fac*(-2*taub+taub/2d0*log4tau**2)/2d0/mb**2

      H1ct = H1ct + 2*mb**2/cosb*GGol_b*2*dmb/4d0

c     shift F

      FFol = -1/6d0*(1/B1-1/B2)

      dmb = dhb

      if(s2b.ne.0d0) then

         H1ct = H1ct + mbsb*s2b*Ab/cosb*FFol*(dmb+ds2b/s2b+dAb/Ab)/4d0
     $        + mz**2*cosb*(dL-dR)*FFol*dc2b/4d0
         
         H2ct = H2ct + mbsb*s2b*mu/cosb*FFol*(dmb+ds2b/s2b)/4d0
     $        - mz**2*sinb*(dL-dR)*FFol*dc2b/4d0

      endif
         
      return
      end

*     
***********************************************************************
*

      subroutine extrashift(SM,q,fac,Ghalf_b,mb,mbsb,cosb,
     $     mg,B1,B2,s2b,CF,dH1)

c     additional shift for DR = 2, i.e. when the one-loop OS bottom
c     contribution is multiplied by mbrun/mbpole

      implicit none

      integer SM
      real*8 q,fac,mb,mbsb,cosb,mg,B1,B2,s2b,CF
      complex*16 Ghalf_b,dH1

      real*8 x1,x2,f1,f2,dmb

      x1 = B1/mg**2
      x2 = B2/mg**2

      f1 = (x1-3)/4d0/(1-x1)+x1*(x1-2)/2d0/(1-x1)**2*log(x1)
      f2 = (x2-3)/4d0/(1-x2)+x2*(x2-2)/2d0/(1-x2)**2*log(x2)

      dmb = -1/4d0*(2*log(mg/q)+f1+f2
     $     +mg/mbsb*s2b*(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2)))

      if(SM.eq.1) dmb = dmb + 3/2d0*log(mb/q) - 5/4d0

      dH1 = -fac*CF*Ghalf_b/cosb*dmb

      return
      end
*     
***********************************************************************
*

      subroutine effshift(ih,DRBAR,fac,Ghalf_b,alpha,tanb,mu,mg,
     $     B1,B2,TF,CF,dH)

c     additional shift for DRBAR = 4 or 5, i.e. when the one-loop OS bottom
c     contribution is rescaled as in the effective Lagrangian approach

      implicit none

      integer ih,DRBAR
      real*8 fac,alpha,tanb,mu,mg,B1,B2,TF,CF
      complex*16 Ghalf_b,H1l,dH

      real*8 x1,x2,Deltab,cosb

      x1 = B1/mg**2
      x2 = B2/mg**2

      cosb = 1d0/sqrt(1+tanb**2)

      Deltab = CF/2d0*mg*mu*tanb/(B1-B2)
     $     *(x1/(1-x1)*log(x1)-x2/(1-x2)*log(x2))

      if(ih.eq.1) then
         
         H1l = -TF*sin(alpha)/cosb*Ghalf_b

         if(DRBAR.eq.4) then
            dH = H1l*Deltab*fac
         else
            dH = H1l*Deltab*(fac+1d0/tan(alpha)/tanb)
         endif

      else

         H1l = TF*cos(alpha)/cosb*Ghalf_b

         if(DRBAR.eq.4) then
            dH = H1l*Deltab*fac
         else
            dH = H1l*Deltab*(fac-tan(alpha)/tanb)
         endif

      endif

      return
      end
C$$$*
c$$$***********************************************************************
c$$$*
c$$$
c$$$      double precision function myB0(q,m1,m2,mu2) 
c$$$
c$$$c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
c$$$      
c$$$      double precision q,m1,m2,Omega,mu2
c$$$
c$$$      if(q.eq.0d0) then
c$$$
c$$$         if(m1.eq.0d0.and.m2.ne.0d0) then
c$$$            myB0 = 1d0-Log(m2/mu2)
c$$$         elseif(m1.ne.0d0.and.m2.eq.0d0) then
c$$$            myB0 = 1d0-Log(m1/mu2)
c$$$         elseif(abs(m1-m2).le.1d-8) then
c$$$            myB0 = -Log(m1/mu2)
c$$$         else
c$$$            myB0 = 1d0 - Log(m2/mu2) + m1/(m1-m2)*Log(m2/m1)
c$$$         endif
c$$$         
c$$$      else
c$$$
c$$$         if(m1.eq.0d0.and.m2.ne.0d0) then
c$$$            
c$$$            if(m2.ne.q) then
c$$$               myB0 = -(Log(m2/mu2)-2-(m2/q-1d0)*Log(abs(1d0-q/m2))) 
c$$$            else 
c$$$               myB0 = -(Log(m2/mu2) - 2)
c$$$            endif
c$$$            
c$$$         elseif(m2.eq.0d0.and.m1.ne.0d0) then
c$$$            
c$$$            if(m1.ne.q) then
c$$$               myB0 = -(Log(m1/mu2)-2-(m1/q-1d0)*Log(abs(1d0-q/m1))) 
c$$$            else
c$$$               myB0 = -(Log(m1/mu2) - 2)
c$$$            endif
c$$$            
c$$$         elseif(m2.eq.0d0.and.m1.eq.0d0) then
c$$$            
c$$$            myB0 = -(Log(q/mu2) - 2) ! cut the imaginary part (I Pi)
c$$$            
c$$$         else
c$$$            
c$$$            myB0 = -( log(q/mu2)-2.d0 + 
c$$$     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*log(m1/q) +
c$$$     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*log(m2/q) +
c$$$     3           2.d0*Omega(m1/q,m2/q))
c$$$            
c$$$         endif
c$$$         
c$$$      endif
c$$$
c$$$      return
c$$$      end
c$$$      
c$$$c     function Omega(a,b) contained in myB0
c$$$c     optimized on 11/01/2010
c$$$
c$$$      double precision function Omega(a,b)
c$$$      double precision a,b,cbig,sqCbig
c$$$      Cbig = 0.5d0*(a+b) - 0.25d0*(a-b)*(a-b) -0.25d0
c$$$      if(Cbig.gt.0d0) then
c$$$         sqCbig = sqrt(Cbig)
c$$$         Omega = sqCbig*
c$$$     1        (atan((1 + a - b)/(2*sqCbig)) +
c$$$     2        atan((1 - a + b)/(2*sqCbig)) )
c$$$      elseif(Cbig.lt.0d0) then
c$$$         sqCbig = sqrt(-Cbig)
c$$$         Omega = 0.5d0*sqCbig*
c$$$     1        log((a + b - 1 - 2*sqCbig)/
c$$$     2        (a + b - 1 + 2*sqCbig))
c$$$      else
c$$$         Omega = 0
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
c$$$
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
c$$$      
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
         

      
