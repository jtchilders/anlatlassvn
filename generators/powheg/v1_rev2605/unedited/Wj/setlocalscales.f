      subroutine setlocalscales(iuborn,imode,rescfac)
c returns the rescaling factor including sudakov form factors and
c coupling rescaling, for Born (imode=1) and NLO corrections (imode=2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      integer iuborn,imode
      real * 8 rescfac,expsudakov,expsud,sudakov,pwhg_alphas
      real * 8 ptb2,mb2,mu2,alphas,b0,optb2,omb2,orescfac,omuf2
      real * 8 pb(0:3)
      integer oimode,i,flav
      save optb2,omb2,orescfac,oimode,omuf2
      data optb2/-1d0/
      logical ini
      data ini/.true./
      save ini
      real * 8 powheginput,factsc2min,frensc2min,as,y,b1,tmp
      save factsc2min,frensc2min,b0,b1
      integer imax
c      real * 8 rescfac1,rescfac2
c      common /crescfac/rescfac1,rescfac2
      
      if(ini) then
         factsc2min = powheginput("#factsc2min")
         frensc2min = powheginput("#frensc2min")
         if(factsc2min.lt.0) factsc2min = 0
         if(frensc2min.lt.0) frensc2min = 1
         call getq2min(1,tmp)
         write(*,*) ' ***** minimum Q of pdf:',sqrt(tmp)
         b0=(33d0-2d0*st_nlight)/(12*pi)
         b1=(153d0-19d0*st_nlight)/(24*pi**2)
         ini = .false.
      endif



      rescfac = 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   These lines are process dependent!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (flst_born(3,iuborn).eq.25 .and. abs(flst_born(4,iuborn)).le.6) 
     $     then
c     gg -> H production
c     Sudakov for a gluon 
         flav=0
      else  
c     Sudakov for a quark
         flav=1  ! any value different from zero
      endif

c     pb(0:3) is the colourless "boson" momentum
      pb(:)=0d0
c     sum over colourless particles (they must all come from a single boson decay)
      do i=3,nlegborn
c     the sequence of colourless particles is unchanged in the Born and in the real flavour list
         if (abs(flst_born(i,iuborn)).gt.6) then
            if(flg_minlo_real) then
               pb(:)=pb(:) + kn_cmpreal(:,i)
            else
               pb(:)=pb(:) + kn_cmpborn(:,i)
            endif
         endif
      enddo

      ptb2 = pb(1)**2 + pb(2)**2 ! transverse momentum squared
      mb2  = pb(0)**2 - pb(3)**2 - ptb2 ! invariant mass squared

      if(imode.eq.oimode.and.ptb2.eq.optb2.and.mb2.eq.omb2) then
         rescfac=orescfac
         st_mufact2=omuf2
         return
      else
         optb2=ptb2
         omb2=mb2
         oimode=imode
      endif
c      if(st_renfact**2*ptb2.lt.frensc2min*st_lambda5MSB**2) then
c         rescfac=0
c         return
c      endif

      as=pwhg_alphas(mb2,st_lambda5MSB,st_nlight)
      y = -as*b0*log(st_renfact**2*ptb2/mb2)

      if(y.ge.1) then
         rescfac = 0d0
         orescfac = 0d0
c         rescfac1 = rescfac
         return
      endif
      
      if(ptb2.gt.mb2) then
         rescfac = 1d0
         expsud = 0d0
      else
         rescfac = sudakov(ptb2,mb2,flav)**2
         expsud  = 2 * expsudakov(ptb2,mb2,flav)
      endif

c      rescfac1 = rescfac


c     alpha_s reweighting
      mu2=ptb2*st_renfact**2
      st_mufact2=max(ptb2*st_facfact**2,factsc2min)
      omuf2=st_mufact2
      alphas=as/(1-y)-as**2*b1*log(1-y)/(b0*(1-y)**2) 

c      write(*,*) "ratio ==>",alphas/pwhg_alphas(st_renfact**2*ptb2,
c     $      st_lambda5MSB,st_nlight)

      if(imode.eq.2) then
         rescfac = rescfac * (alphas/st_alpha)**2
      else
         rescfac = rescfac * (alphas/st_alpha)
         if(.not.flg_bornonly) then
            rescfac = rescfac * 
     1           (1+alphas*(expsud+b0*log(mu2/st_muren2)))
         endif
      endif
      orescfac=rescfac

c      rescfac2 = rescfac

      end


C     Inputs:                                      
C     *******                                      
C     q2h  - Upper node scale / bound on Sudakov   
C     q20  - Lower node scale / bound on Sudakov   
C     flav - flavour index for the evolving parton 
C                                                  
C     Outputs:                                     
C     ********                                     
C     sudakov - The Sudakov form factor         
      function sudakov(q20,q2h,flav)
      implicit none
      real * 8 sudakov,q2h,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
      real * 8 lam2,b0,y,as,pwhg_alphas
      logical isQuark
      real * 8 theExponentN
      real * 8 powheginput
      logical ini,sudscalevar
      data ini/.true./
      save ini,sudscalevar,b0

      if(ini) then
         if(powheginput("#sudscalevar").eq.0) then
            sudscalevar = .false.
         else
            sudscalevar = .true.
         endif
         b0 = (11*CA-2d0*st_nlight)/(12*pi)
         ini = .false.
      endif
      if(sudscalevar) then
         lam2=st_lambda5MSB**2 / st_renfact**2
      else
         lam2=st_lambda5MSB**2
      endif

c     running coupling coefficients 
      as=pwhg_alphas(q2h,st_lambda5MSB,st_nlight)
      y = -as*b0*log(q20*st_renfact**2/q2h)

      if (y .ge. 1) then 
         sudakov = 0d0 
         return 
      endif

      if(sudscalevar) then
         if(q2h.le.q20*st_renfact**2) then
            sudakov=1
            return
         endif
      else
         if(q2h.le.q20) then
            sudakov=1
            return
         endif
      endif
      if(flav.eq.0) then
         isQuark=.false.
      else
         isQuark=.true.
      endif

      if(flg_minlo_nnll) then
         call sudakov_exponent(q20,q2h,theExponentN,isQuark,1)
      else
         call sudakov_exponent(q20,q2h,theExponentN,isQuark,2)
      endif

      sudakov=exp(theExponentN)

      end

C     Inputs:                                      
C     *******                                      
C     q2h  : Upper node scale / bound on Sudakov   
C     q20  : Lower node scale / bound on Sudakov   
C     flav : flavour index for the evolving parton 
C                                               
C     Outputs:                                     
C     ********                                     
C     expsudakov : the expansion at order alpha_s of the exponent of the Sudakov form factor
c                  with a minus sign in front, without the alpha_s factor.
c     See Eq. (2.9) of arXiv:1212.4504
      function expsudakov(q20,q2h,flav)
      implicit none
      real * 8 expsudakov,q2h,q20
      integer flav
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'pwhg_flg.h'
      real * 8 b0,logf,m2,a1,b1,llog
      real * 8 powheginput
      logical ini,sudscalevar
      data ini/.true./
      save ini,sudscalevar,b0

      if(ini) then
         if(powheginput("#sudscalevar").eq.1) then
            sudscalevar = .true.
         else
            sudscalevar = .false.
         endif
         b0=(11*CA-2*st_nlight)/(12*pi)
         ini = .false.
      endif

      if(sudscalevar) then
         logf = log(st_renfact)
         m2 = q20*st_renfact**2
      else
         logf = 0
         m2 = q20
      endif

      if(q2h.le.m2.or.flg_bornonly) then
         expsudakov=0
         return
      endif

      if(flav.eq.0) then
         A1 =  CA
         B1 = -(11d0*CA-2d0*st_nlight)/6 
      else
         A1 = CF 
         B1 = -3d0/2*CF
      endif

c     we need the coefficients for A(as) and B(as) written as
c     A(as) = sum_n A_n (as)^n,     B(as) = sum_n B_n (as)^n,   i.e. no "pi" factors
      A1 = A1/pi
      B1 = B1/pi
c     introduce scale dependence
      B1 = B1 + 2*A1*logf 
c     Add an extra 1/2 in front of the expansion, since this is the sqrt of the final Sudakov
c     Add an overall minus sign, since this is to be subtracted to the Bbar function 
      llog = log(m2/q2h)
      expsudakov = -1d0/2 * (-1d0/2 * A1 * llog**2 + B1 * llog) ! Eq. (2.9) of arXiv:1212.4504
      end



C     The integral of the Sudakov exponent. See appendix A and B of arXiv:1212.4504      
C     isQuark = .true. for a quark propagator                        
C     theAccuracy = 0 for 1-loop alphaS and A2=B2=0,                 
C                 = 1 for 2-loop alphaS and Powheg A & B coeffs      
C                 = 2 for 2-loop alphaS and NLL A & B coeffs         
C                                                                     
      subroutine sudakov_exponent(q20,q2h,theExponent,
     $     isquark,theAccuracy)
      implicit none 
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      logical isQuark 
      integer theAccuracy
      real * 8 q2h,q20,theExponent
      real * 8 b0,be1,be2
      real * 8 A1,A2,A3
      real * 8 B1,B2
      real * 8 y, as
      real * 8 nf,K,logf,q20_lcl  
      real * 8 EulerGamma,zeta3  
      real * 8 f0,f1,f2,omy,lomy
      real * 8 pwhg_alphas,powheginput
      external pwhg_alphas,powheginput
      logical ini,sudscalevar
      data ini/.true./
      save ini,sudscalevar
      parameter (zeta3 = 1.2020569031595942854d0)
      parameter (EulerGamma = 0.57721566490153286061d0)
      
      if(ini) then
         if(powheginput("#sudscalevar").eq.1) then
            sudscalevar = .true.
         else
            sudscalevar = .false.
         endif
         ini = .false.
      endif

      if(sudscalevar) then
         logf = log(st_renfact)
         q20_lcl = q20*st_renfact**2
      else
         logf = 0d0
         q20_lcl = q20 
      endif


      nf = st_nlight         
c     running coupling coefficients 
c     d alpha/d log mu^2=-b0 alpha^2 - be1 alpha^3 -be2 alpha^4
      b0 = (11d0*CA-2d0*nf)/(12*pi)
      be1 = (153d0 - 19d0*nf)/(24*pi**2)
      be2 = 0d0                 ! not needed 
      K = (67d0/18d0-Pi**2/6d0)*CA-5d0/9d0*nf

      as=pwhg_alphas(q2h,st_lambda5MSB,st_nlight)
      y = -as*b0*log(q20_lcl/q2h)

      if (y .ge. 1) then 
         write(*,*) '-------> y',y
         theExponent = -1000d0 
         return 
      endif

c     Sudakov exponent coefficients for 
c     A(as) = sum_n A_n (as/pi)^n,     B(as) = sum_n B_n (as/pi)^n,    
c     A1, A2, B1, B2 taken from Catani, de Florian, Grazzini, hep-ph/0008184
      if (isQuark) then 
         A1 = CF  
         A2 = CF*K/2
         A3 = 0d0 ! not needed 
         B1 = -3d0/2*CF 
         B2 =  (CF**2 * (Pi**2-3d0/4-12*zeta3)
     $        + CF*CA * (11d0*Pi**2/9-193d0/12+6*zeta3)
     $        + CF*nf*0.5 * (17d0/3-4d0/9*Pi**2))/4
      else
         A1 = CA
         A2 = CA*K/2
         A3 = 0d0 ! not needed 
         B1 = -(11d0*CA-2d0*nf)/6
         B2 = ( CA**2 * (23d0/6+22*Pi**2/9-6*zeta3)
     $        + 2*CF*nf
     $        - CA*nf*(1d0/3+4d0*Pi**2/9) - 11d0/2*CA*CF)/4
      endif

c     we need the coefficients for A(as) and B(as) written as
c     A(as) = sum_n A_n (as)^n,     B(as) = sum_n B_n (as)^n,   i.e. no "pi" factors
      A1 = A1/pi
      B1 = B1/pi
      A2 = A2/pi**2
      B2 = B2/pi**2

c     qT space conversion
      B2 = B2 + 2*A1**2*zeta3

c     add logf dependence 
      B2 = B2 + 2*A2*logf + 2*b0*A1*logf**2 
      A2 = A2 + 2*b0*A1*logf
      B1 = B1 + 2*A1*logf 

c     scale variation induced by the scale variation of the H factor 
      B2 = B2 + 2*(st_bornorder-1)*b0**2*logf

c     switch off some coefficients according to accuracy wanted 
      if(theAccuracy.eq.0) then
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      elseif(theAccuracy.eq.2) then
         B2  = 0d0 ! NNLL coefficient
      elseif(theAccuracy.eq.3) then
         B1  = 0d0 ! NLL coefficient
         A2  = 0d0 ! NLL coefficient
         B2  = 0d0 ! NNLL coefficient
      endif

      omy = 1-y 
      lomy=log(omy)

      f0 = A1/b0**2*(y+lomy)

      f1=A1*be1/b0**3*(lomy**2/2+y/omy+lomy/omy)
     $     -A2/b0**2*(lomy+y/omy)+B1/b0*lomy  

      f2 = -B2/b0*y/omy+B1*be1/b0**2*(y+lomy)/omy
     $     -A3/(2*b0**2)*(y/omy)**2
     $     +A2*be1/(2*b0**3)*(3*y**2-2*y+(4*y-2)*lomy)/omy**2
     $     + A1/(2*b0**4)/omy**2*(
     $     be1**2*(1-2*y)*lomy**2
     $     +2*(b0*be2*omy**2+be1**2*y*omy)*lomy
     $     -3*b0*be2*y**2+be1**2*y**2+2*b0*be2*y)

      theExponent = 1d0/as*f0+f1+as*f2 

c     The Sudakov form factor computed up to now is the boson Sudakov form factor,
c     obtained in MiNLO by squaring the one returned by this subroutine.
c     For this reason, we have to divide the argument of the exponent by a factor of 2
      theExponent =  theExponent/2

      end




